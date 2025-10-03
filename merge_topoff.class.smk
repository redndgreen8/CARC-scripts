import os
from pathlib import Path
from collections import defaultdict
import glob
import logging
from typing import Dict, List, Tuple, Optional

# Configuration
configfile: "config.yaml"

# Convert paths to absolute
INPUT_DIRS = [os.path.abspath(d) for d in config["input_dirs"]]
TOPOFF_DIR = os.path.abspath(config["topoff_dir"])
OUTPUT_DIR = os.path.abspath(config["output_dir"])

# Create output directory and logs
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "logs"), exist_ok=True)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(os.path.join(OUTPUT_DIR, "logs", "workflow.log")),
        logging.StreamHandler()
    ]
)

class FileStats:
    def __init__(self):
        self.total_files = 0
        self.undetermined_files = 0
        self.processed_files = 0
        self.samples_with_topoff = 0
        self.samples_without_topoff = 0
        self.total_merges = 0
        self.total_links = 0

    def log_summary(self):
        logging.info("\nProcessing Statistics:")
        logging.info(f"  Total files found: {self.total_files}")
        logging.info(f"  Undetermined files skipped: {self.undetermined_files}")
        logging.info(f"  Files processed: {self.processed_files}")
        logging.info(f"  Samples with topoff: {self.samples_with_topoff}")
        logging.info(f"  Samples without topoff: {self.samples_without_topoff}")
        logging.info(f"  Total merge operations: {self.total_merges}")
        logging.info(f"  Total symlinks: {self.total_links}")

STATS = FileStats()

def parse_fastq_name(filepath: str) -> Optional[Tuple[str, str, str]]:
    """
    Parse FASTQ filename into components.
    Returns: (sample_id, lane, read_orientation) or None if invalid/undetermined
    Sample_id includes the S number to ensure proper pairing.
    """
    STATS.total_files += 1
    filename = os.path.basename(filepath)
    
    # Skip Undetermined
    if filename.startswith('Undetermined'):
        STATS.undetermined_files += 1
        return None
        
    try:
        parts = filename.split('_')
        
        # Find lane and R1/R2 indices
        lane_idx = next(i for i, part in enumerate(parts) if part.startswith('L00'))
        read_idx = next(i for i, part in enumerate(parts) if part.startswith('R'))
        
        # Include S number in sample_id to ensure proper grouping
        s_idx = next(i for i, part in enumerate(parts) if part.startswith('S'))
        sample_id = '_'.join(parts[:s_idx + 1])  # Include everything up to and including S number
        
        lane = parts[lane_idx]
        read_orientation = parts[read_idx]  # R1 or R2
        
        STATS.processed_files += 1
        return sample_id, lane, read_orientation
        
    except Exception as e:
        logging.warning(f"Failed to parse {filename}: {str(e)}")
        return None

def get_merge_groups() -> Dict[str, List[str]]:
    """
    Find all files that need to be merged across input dirs and topoff.
    Keeps R1 and R2 separate, only merging across lanes for each read orientation.
    Returns: dict mapping output files to input file lists
    """
    merge_dict = defaultdict(list)
    
    # Group files by sample ID and read orientation separately
    input_files = defaultdict(list)  # (sample_id, read_orientation) -> [files]
    for input_dir in INPUT_DIRS:
        for fastq in glob.glob(os.path.join(input_dir, "*.fastq.gz")):
            result = parse_fastq_name(fastq)
            if result:
                sample_id, lane, read_orientation = result
                # Keep separate entries for R1 and R2
                key = (sample_id, read_orientation)
                input_files[key].append(fastq)
    
    # Find matching topoff files, keeping R1 and R2 separate
    topoff_files = defaultdict(list)
    for fastq in glob.glob(os.path.join(TOPOFF_DIR, "*.fastq.gz")):
        result = parse_fastq_name(fastq)
        if result:
            sample_id, lane, read_orientation = result
            key = (sample_id, read_orientation)
            topoff_files[key].append(fastq)
    
    # Create merge groups, keeping R1 and R2 separate
    for key, files in input_files.items():
        sample_id, read_orientation = key
        matching_topoff = topoff_files[key]
        
        # Get lanes for logging
        input_lanes = sorted(set(parse_fastq_name(f)[1] for f in files))
        topoff_lanes = sorted(set(parse_fastq_name(f)[1] for f in matching_topoff))
        
        # Create output name with L001 regardless of input lanes
        base_file = sorted(files)[0]
        output_name = os.path.basename(base_file)
        # Explicitly force L001 in output name
        for lane in input_lanes:
            if lane in output_name:
                output_name = output_name.replace(lane, 'L001')
                break
        
        # Combine and sort all files for this read orientation
        all_files = sorted(files + matching_topoff)
        
        # Update stats
        if matching_topoff:
            STATS.samples_with_topoff += 1
            STATS.total_merges += 1
        else:
            STATS.samples_without_topoff += 1
            if len(files) > 1:
                STATS.total_merges += 1
            else:
                STATS.total_links += 1
            
        merge_dict[output_name] = all_files
        
        # Enhanced logging
        logging.info(f"\nSample: {sample_id} Read: {read_orientation}")
        logging.info(f"  Input lanes: {', '.join(input_lanes)}")
        logging.info(f"  Input files: {len(files)}")
        if matching_topoff:
            logging.info(f"  Topoff lanes: {', '.join(topoff_lanes)}")
            logging.info(f"  Topoff files: {len(matching_topoff)}")
        else:
            logging.warning(f"  No topoff files found")
        logging.info(f"  Output file: {output_name} (merged into L001)")
    
    return merge_dict

# Generate merge information
MERGE_INFO = get_merge_groups()
if not MERGE_INFO:
    raise ValueError("No valid files found to process!")

# Log final statistics
STATS.log_summary()

# Target rule
rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{output_file}"),
               output_file=MERGE_INFO.keys())

def merge_or_link_files(input_files: List[str], output_file: str, log_file: str):
    """Merge multiple files or create symlink for single file."""
    with open(log_file, 'w') as log:
        log.write(f"Processing: {os.path.basename(output_file)}\n")
        log.write(f"Input files ({len(input_files)}):\n")
        for idx, fname in enumerate(input_files, 1):
            log.write(f"{idx}. {fname}\n")
        
        try:
            if len(input_files) > 1:
                log.write("\nPerforming merge\n")
                with open(output_file, 'wb') as outfile:
                    for fname in input_files:
                        log.write(f"Reading: {fname}\n")
                        with open(fname, 'rb') as infile:
                            while chunk := infile.read(8192):
                                outfile.write(chunk)
            else:
                log.write("\nCreating symlink\n")
                if os.path.exists(output_file):
                    os.remove(output_file)
                os.symlink(os.path.abspath(input_files[0]), output_file)
                
        except Exception as e:
            log.write(f"Error: {str(e)}\n")
            raise

rule merge_fastq:
    input:
        lambda wildcards: MERGE_INFO[wildcards.output_file]
    output:
        os.path.join(OUTPUT_DIR, "{output_file}")
    log:
        os.path.join(OUTPUT_DIR, "logs", "{output_file}.log")
    resources:
        mem_mb=4000
    threads: 1
    run:
        merge_or_link_files(input, output[0], log[0])