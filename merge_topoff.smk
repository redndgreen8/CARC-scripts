import os
from pathlib import Path
from collections import defaultdict
import glob
import json
import yaml

current_dir = Path.cwd()  # or os.getcwd()

print(current_dir)

# Configuration
configfile: config["config_file"]



#print(f"Config file path: {configfile}")
print(config["input_dirs"])

# Convert input paths to absolute
INPUT_DIRS = [os.path.abspath(d) for d in config["input_dirs"]]
#INPUT_DIRS = os.path.abspath(config["input_dirs"])


TOPOFF_DIR = os.path.abspath(config["topoff_dir"])
OUTPUT_DIR = os.path.abspath(config["output_dir"])

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "logs"), exist_ok=True)

def find_fastq_files():
    """
    Find and group FASTQ files across input directories and topoff directory.
    Matches on first n fields (before lane) + read orientation.
    Merges across lanes (L001-L004) into L001 in output name.
    """
    merge_groups = defaultdict(lambda: {'input': [], 'topoff': []})
    total_input_files = 0
    total_topoff_files = 0
    
    # Process input directories
    for input_dir in INPUT_DIRS:
        if not os.path.exists(input_dir):
            print(f"Warning: Input directory not found: {input_dir}")
            continue
            
        fastq_files = [f for f in glob.glob(os.path.join(input_dir, "*.fastq.gz")) 
                      if "Undetermined" not in f]
        total_input_files += len(fastq_files)
            
        for fastq in fastq_files:
            abs_path = os.path.abspath(fastq)
            basename = os.path.basename(fastq)
            
            # Split into parts
            parts = basename.split('_')
            
            # Find read orientation and lane
            read_number = None
            lane_number = None
            for part in parts:
                if part.startswith('R'):
                    read_number = part
                elif part.startswith('L00'):
                    lane_number = part
            
            # Skip if missing required fields or not enough parts
            if not read_number or not lane_number or len(parts) < 8:
                print(f"Warning: Skipping {basename} - missing required fields or insufficient parts")
                continue
            
            # Create matching key using first 9 fields (before lane) + read orientation
            matching_key = f"{'_'.join(parts[:8])}_{read_number}"
            
            # Create output name (replace lane with L001)
            output_parts = parts.copy()
            for i, part in enumerate(output_parts):
                if part.startswith('L00'):
                    output_parts[i] = 'L001'
                    break
            output_name = '_'.join(output_parts)
            
            merge_groups[matching_key]['input'].append((abs_path, output_name))
    
    # Process topoff directory
    if os.path.exists(TOPOFF_DIR):
        fastq_files = [f for f in glob.glob(os.path.join(TOPOFF_DIR, "*.fastq.gz")) 
                      if "Undetermined" not in f]
        total_topoff_files = len(fastq_files)
        
        for fastq in fastq_files:
            abs_path = os.path.abspath(fastq)
            basename = os.path.basename(fastq)
            
            parts = basename.split('_')
            
            # Find read orientation and lane
            read_number = None
            lane_number = None
            for part in parts:
                if part.startswith('R'):
                    read_number = part
                elif part.startswith('L00'):
                    lane_number = part
            
            # Skip if missing required fields or not enough parts
            if not read_number or not lane_number or len(parts) < 8:
                print(f"Warning: Skipping {basename} - missing required fields or insufficient parts")
                continue
            
            # Create matching key using first 9 fields (before lane) + read orientation
            matching_key = f"{'_'.join(parts[:8])}_{read_number}"
            
            # Create output name (replace lane with L001)
            output_parts = parts.copy()
            for i, part in enumerate(output_parts):
                if part.startswith('L00'):
                    output_parts[i] = 'L001'
                    break
            output_name = '_'.join(output_parts)
            
            merge_groups[matching_key]['topoff'].append((abs_path, output_name))
    
    # Convert to final merge format and validate
    final_merge_groups = {}
    input_file_count = 0
    topoff_file_count = 0
    
    for matching_key, files in merge_groups.items():
        input_files = sorted(files['input'])  # List of (path, output_name) tuples
        topoff_files = sorted(files['topoff'])  # List of (path, output_name) tuples
        
        # Skip if missing either input or topoff files
        if not input_files or not topoff_files:
            continue
        
        # Get output name (should be same for all files in group except for lane number)
        output_name = input_files[0][1]
        
        # Add paths to final groups
        all_files = [f[0] for f in input_files + topoff_files]  # Just the paths
        final_merge_groups[output_name] = all_files
        
        # Update counts
        input_file_count += len(input_files)
        topoff_file_count += len(topoff_files)
    
    # Print summary
    print(f"\nSummary of files found:")
    print(f"Total input files found: {total_input_files}")
    print(f"Total topoff files found: {total_topoff_files}")
    print(f"\nSummary of files to be merged:")
    print(f"Number of sample groups with matching input and topoff files: {len(final_merge_groups)}")
    print(f"Input files to be merged: {input_file_count}")
    print(f"Topoff files to be merged: {topoff_file_count}")
    print(f"Total files to be merged: {input_file_count + topoff_file_count}")
    print(f"Input files without matching topoff: {total_input_files - input_file_count}")
    print(f"Topoff files without matching input: {total_topoff_files - topoff_file_count}\n")
    
    return final_merge_groups

def get_read_orientation(filepath):
    """Extract R1/R2 read orientation from filename."""
    basename = os.path.basename(filepath)
    for part in basename.split('_'):
        if part.startswith('R'):
            return part
    return None

# Generate merge information
MERGE_INFO = find_fastq_files()
if not MERGE_INFO:
    raise ValueError("No valid files found to process!")

# Define output directories
FASTQC_DIR = os.path.join(OUTPUT_DIR, "fastqc")
LOGS_DIR = os.path.join(OUTPUT_DIR, "logs")
MQC_DIR = os.path.join(OUTPUT_DIR, "multiqc")

# Ensure directories exist
os.makedirs(LOGS_DIR, exist_ok=True)
os.makedirs(FASTQC_DIR, exist_ok=True)
os.makedirs(MQC_DIR, exist_ok=True)



# Target rule
rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{output_file}"),
               output_file=MERGE_INFO.keys()),
        os.path.join(OUTPUT_DIR, ".fastqcCARCPass"),
        os.path.join(OUTPUT_DIR, ".multiqcCARCPass")

# Merge rule
rule merge_fastq:
    input:
        lambda wildcards: MERGE_INFO[wildcards.output_file]
    output:
        os.path.join(OUTPUT_DIR, "{output_file}")
    log:
        os.path.join(OUTPUT_DIR, "logs", "{output_file}.log")
    shell:
        """
        if [ "$(echo {input} | wc -w)" -gt 1 ]; then
            cat {input} > {output} 2> {log}
        else
            ln -sf {input} {output} 2> {log}
        fi

        """

# Flag rule to indicate all merges are complete
rule merge_complete:
    input:
        expand(os.path.join(OUTPUT_DIR, "{output_file}"),
               output_file=MERGE_INFO.keys())
    output:
        touch(os.path.join(OUTPUT_DIR, ".mergeFastqCARCPass"))
    log:
        os.path.join(LOGS_DIR, "merge_complete.log")
    shell:
        """
        echo "Completed merging of $(echo {input} | wc -w) files at $(date)" > {log}
        """



# FastQC rule with dependency on merge completion
rule fastqc:
    input:
        merged_files = expand(os.path.join(OUTPUT_DIR, "{output_file}"), 
                             output_file=MERGE_INFO.keys()),
        merge_flag = os.path.join(OUTPUT_DIR, ".mergeFastqCARCPass")
    output:
        outdir = directory(FASTQC_DIR),
        # Generate expected FastQC output files based on input patterns
        flag = touch(os.path.join(OUTPUT_DIR, ".fastqcCARCPass"))
    log:
        os.path.join(LOGS_DIR, "fastqc.log")
    shell:
        """
        mkdir -p {output.outdir}
        module load usc fastqc
        echo "Starting FastQC analysis at $(date)" > {log}
        fastqc --threads 16 --outdir {output.outdir} {input.merged_files} &>> {log}
        echo "Completed FastQC analysis at $(date)" >> {log}
       # touch {output.flag}
        """



rule multiQC:
    input:
        flag = os.path.join(OUTPUT_DIR, ".fastqcCARCPass"),
        fastqc_dir = FASTQC_DIR
    output:
        outdir = directory(MQC_DIR),
        flag = touch(os.path.join(OUTPUT_DIR, ".multiqcCARCPass"))
    log:
        os.path.join(LOGS_DIR, "multiqc.log")
    shell:
        """
        mkdir -p {output.outdir}
        export PYTHONPATH=/project/davidwcr_264/Packages/pythonModule/python-3.9.2
        export OPENBLAS_NUM_THREADS=4
        echo "Starting MultiQC analysis at $(date)" > {log}
        /project/davidwcr_264/Packages/pythonModule/python-3.9.2/bin/multiqc \
            --interactive \
            -n multiqc_merge \
            -o {output.outdir} \
            {input.fastqc_dir} \
            &>> {log}        
            echo "Completed MultiQC analysis at $(date)" >> {log}
        """