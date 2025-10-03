#!/usr/bin/env python3
"""
Simplified Bioinformatics File System Crawler

This script crawls a file system containing bioinformatics data, extracts metadata and 
content samples, and generates summaries suitable for analysis by an LLM.

Optimized for stability and simplicity with minimal dependencies.
"""

import os
import sys
import time
import json
import re
import gzip
import random
import logging
from datetime import datetime
from collections import Counter, defaultdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("bio_crawler.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define extensions to sample content from
BIOINFORMATICS_EXTENSIONS = {
    # Sequence data
    '.fastq', '.fq', '.fastq.gz', '.fq.gz',  # Raw sequencing reads
    '.fasta', '.fa', '.fasta.gz', '.fa.gz',  # Reference sequences
    '.sam', '.bam', '.cram',  # Aligned sequences
    
    # Variant data
    '.vcf', '.vcf.gz', '.bcf',  # Variant calls
    '.maf', '.maf.gz',  # Mutation annotation format
    
    # Analysis files
    '.bed', '.bed.gz', '.bedgraph', '.bigwig', '.bw',  # Genomic intervals
    
    # Workflow/Tool files
    '.json', '.yaml', '.yml', '.toml', '.config',  # Configuration
    '.cwl', '.wdl', '.nf', '.smk', '.snake',  # Workflow languages
    '.sh', '.py', '.r', '.pl',  # Scripts
    
    # Documentation
    '.md', '.txt', '.csv', '.tsv', '.log', '.out',  # Text files
    '.html', '.pdf'  # Documentation
}

# Common patterns for identifying bioinformatics tools
BIOINFORMATICS_TOOL_PATTERNS = [
    # Common variant callers
    r'mutect2?', r'varscan', r'strelka2?', r'manta', r'freebayes', r'gatk', r'lofreq',
    
    # Aligners
    r'bwa', r'bowtie2?', r'star', r'hisat2?', r'minimap2', r'tophat',
    
    # Other tools
    r'samtools', r'bcftools', r'bedtools', r'picard', r'trimmomatic', r'fastp',
    r'snpeff', r'vep', r'annovar', r'oncotator', r'funcotator',
    
    # Workflow tools
    r'nextflow', r'snakemake', r'cromwell', r'cwltool', r'wdltool', r'galaxy', r'wes',
    
    # Common analysis terms
    r'somatic', r'germline', r'cnv', r'rnaseq', r'chip-?seq', r'atac-?seq',
    r'tumor', r'normal', r'control', r'case', r'patient', r'sample'
]

# Compile patterns for efficiency
TOOL_PATTERN = re.compile('|'.join(BIOINFORMATICS_TOOL_PATTERNS), re.IGNORECASE)

class SimpleBioCrawler:
    def __init__(self, root_dir, output_dir="crawler_output", 
                 max_file_size=10*1024*1024, sample_size=4*1024,
                 max_files_per_dir=5, memory_dump_interval=10000):
        """
        Initialize the crawler.
        
        Args:
            root_dir: Root directory to start crawling from
            output_dir: Directory to save output files
            max_file_size: Maximum file size to sample content from (bytes)
            sample_size: Size of content sample to extract (bytes)
            max_files_per_dir: Maximum number of files to sample per directory
            memory_dump_interval: Number of directories between data dumps to disk
        """
        self.root_dir = os.path.abspath(root_dir)
        self.output_dir = output_dir
        self.max_file_size = max_file_size
        self.sample_size = sample_size
        self.max_files_per_dir = max_files_per_dir
        self.memory_dump_interval = memory_dump_interval
        
        # Statistics and aggregation
        self.stats = {
            "total_dirs": 0,
            "total_files": 0,
            "total_size": 0,
            "file_types": Counter(),
            "file_ages": defaultdict(int),
            "potential_tools": Counter(),
            "largest_files": [],
            "oldest_files": [],
            "newest_files": []
        }
        
        # Data stores
        self.directory_summaries = {}
        self.interesting_files = []
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, "data_dumps"), exist_ok=True)
        
        # Initialize counter for memory dumps
        self.dirs_since_dump = 0
        self.dump_counter = 0
        
    def crawl(self):
        """Main crawl method that orchestrates the filesystem traversal"""
        start_time = time.time()
        logger.info(f"Starting crawl of {self.root_dir}")
        
        try:
            # Simple recursive walk with standard os.walk
            for current_dir, dirs, files in os.walk(self.root_dir):
                # Skip hidden directories
                dirs[:] = [d for d in dirs if not d.startswith('.')]
                
                # Process this directory
                self._process_directory(current_dir, dirs, files)
                
                # Check if we should dump memory
                self.dirs_since_dump += 1
                if self.dirs_since_dump >= self.memory_dump_interval:
                    self._dump_memory()
                    self.dirs_since_dump = 0
                
        except KeyboardInterrupt:
            logger.warning("Crawl interrupted by user")
            self._dump_memory("interrupted")
        except Exception as e:
            logger.error(f"Error during crawl: {str(e)}", exc_info=True)
            self._dump_memory("error")
        
        # Final data processing
        self._finalize_crawl()
        
        # Report statistics
        elapsed = time.time() - start_time
        hrs, remainder = divmod(elapsed, 3600)
        mins, secs = divmod(remainder, 60)
        
        logger.info(f"Crawl completed in {int(hrs)}h {int(mins)}m {int(secs)}s")
        logger.info(f"Found {self.stats['total_files']:,} files in {self.stats['total_dirs']:,} directories")
        logger.info(f"Total size: {self._format_size(self.stats['total_size'])}")
        logger.info(f"Memory dumps performed: {self.dump_counter}")
        
        return self.get_summary()
    
    def _process_directory(self, current_dir, dirs, files):
        """Process a single directory and its files"""
        rel_path = os.path.relpath(current_dir, self.root_dir)
        self.stats["total_dirs"] += 1
        
        # Provide progress updates
        if self.stats["total_dirs"] % 100 == 0:
            elapsed = time.time() - (getattr(self, 'start_time', time.time()) or time.time())
            dirs_per_sec = self.stats["total_dirs"] / max(1, elapsed)
            logger.info(f"Processed {self.stats['total_dirs']:,} directories "
                        f"({dirs_per_sec:.1f}/s), "
                        f"{self.stats['total_files']:,} files")
        
        # Handle very large directories by sampling
        if len(files) > 1000:
            logger.info(f"Large directory with {len(files)} files: {current_dir}")
            # Sample a subset of files for very large directories
            files_sample = files[:50]  # First 50
            if len(files) > 50:
                try:
                    files_sample += random.sample(files[50:], min(50, len(files)-50))  # Random 50 more
                except:
                    # If random.sample fails, just take the first 100
                    files_sample = files[:100]
            files = files_sample
        
        # Process files sequentially
        file_entries = []
        dir_size = 0
        file_count = 0
        
        for filename in files:
            if filename.startswith('.'):
                continue
                
            file_entry = self._process_file(current_dir, filename)
            
            if file_entry:
                file_entries.append(file_entry)
                file_count += 1
                if 'size' in file_entry:
                    dir_size += file_entry['size']
        
        # Store directory summary if it has files
        if file_entries:
            self.directory_summaries[rel_path] = {
                "path": rel_path,
                "file_count": len(files),
                "dir_count": len(dirs),
                "dir_size": dir_size,
                "dir_size_human": self._format_size(dir_size),
                "samples": file_entries[:self.max_files_per_dir]
            }
        
        # Check directory name and path for tool mentions
        self._check_for_tools(current_dir)
    
    def _process_file(self, directory, filename):
        """Process a single file, extracting metadata and content samples if relevant"""
        filepath = os.path.join(directory, filename)
        
        try:
            # Skip symbolic links
            if os.path.islink(filepath):
                return None
            
            # Get file stats
            file_stat = os.stat(filepath)
            file_size = file_stat.st_size
            create_time = file_stat.st_ctime
            modify_time = file_stat.st_mtime
            
            # Update statistics
            self.stats["total_files"] += 1
            self.stats["total_size"] += file_size
            
            # Get file extension and update counters
            _, ext = os.path.splitext(filename.lower())
            if ext == '.gz' or ext == '.zip' or ext == '.bz2':
                # Handle double extensions like .vcf.gz
                base, first_ext = os.path.splitext(filename[:-len(ext)])
                if first_ext:
                    ext = first_ext + ext
            self.stats["file_types"][ext] += 1
            
            # Age category
            age_days = (time.time() - modify_time) / (60*60*24)
            if age_days < 7:
                age_category = "last_week"
            elif age_days < 30:
                age_category = "last_month"
            elif age_days < 90:
                age_category = "last_quarter"
            elif age_days < 365:
                age_category = "last_year"
            else:
                age_category = "older"
            self.stats["file_ages"][age_category] += 1
            
            # Format for output
            date_str = datetime.fromtimestamp(modify_time).strftime("%Y-%m-%d")
            file_entry = {
                "path": os.path.relpath(filepath, self.root_dir),
                "size": file_size,
                "size_human": self._format_size(file_size),
                "modified": date_str
            }
            
            # Keep track of extremes for summary
            self._update_extremes(file_entry)
            
            # Add file extension to entry for filtering
            file_entry["extension"] = ext
            
            # Check filename for tool mentions
            if TOOL_PATTERN.search(filename):
                matched_tools = TOOL_PATTERN.findall(filename)
                for tool in matched_tools:
                    self.stats["potential_tools"][tool.lower()] += 1
                file_entry["potential_tools"] = matched_tools
                
            # Sample content if it's a relevant file type and not too large
            if self._should_sample(filepath, file_size, ext):
                sample = self._get_file_sample(filepath, file_size, ext)
                if sample:
                    file_entry["sample"] = sample
                    self.interesting_files.append(file_entry)
            
            return file_entry
            
        except Exception as e:
            logger.error(f"Error processing file {filepath}: {str(e)}")
            return None
    
    def _should_sample(self, filepath, file_size, extension):
        """Determine if we should sample content from this file"""
        # Don't sample if the file is too large
        if file_size > self.max_file_size:
            return False
        
        # Sample if it's a bioinformatics file type
        for ext in BIOINFORMATICS_EXTENSIONS:
            if extension.endswith(ext):
                return True
        
        # Sample if the filename matches tool patterns
        if TOOL_PATTERN.search(os.path.basename(filepath)):
            return True
            
        return False
    
    def _get_file_sample(self, filepath, file_size, extension):
        """Extract a sample from the beginning of the file"""
        try:
            sample = ""
            sample_size = min(self.sample_size, file_size)
            
            # Handle different file types
            if extension.endswith('.gz'):
                try:
                    with gzip.open(filepath, 'rt', errors='replace') as f:
                        sample = f.read(sample_size)
                except:
                    sample = "Gzipped file, could not read content"
            elif extension in ('.sam', '.bam', '.cram'):
                # For binary files, just return a placeholder
                sample = f"Binary alignment file ({extension})"
            elif extension in ('.vcf', '.vcf.gz'):
                # Try to extract VCF header
                if extension == '.vcf':
                    try:
                        with open(filepath, 'r', errors='replace') as f:
                            lines = []
                            for i, line in enumerate(f):
                                if i > 50:  # Limit reading to 50 lines
                                    break
                                if line.startswith('#'):
                                    lines.append(line.strip())
                                else:
                                    # Add the first data line
                                    lines.append(line.strip())
                                    break
                            sample = '\n'.join(lines)
                    except:
                        sample = "VCF file, could not read content"
                else:  # .vcf.gz
                    try:
                        with gzip.open(filepath, 'rt', errors='replace') as f:
                            lines = []
                            for i, line in enumerate(f):
                                if i > 50:  # Limit reading to 50 lines
                                    break
                                if line.startswith('#'):
                                    lines.append(line.strip())
                                else:
                                    # Add the first data line
                                    lines.append(line.strip())
                                    break
                            sample = '\n'.join(lines)
                    except:
                        sample = "Gzipped VCF file, could not read content"
            else:
                # Regular text file
                try:
                    with open(filepath, 'r', errors='replace', encoding='utf-8') as f:
                        sample = f.read(sample_size)
                except:
                    # Fallback if UTF-8 fails
                    try:
                        with open(filepath, 'r', errors='replace', encoding='latin-1') as f:
                            sample = f.read(sample_size)
                    except:
                        sample = "Text file, could not read content"
            
            # Check if we found potential tools in the content
            if sample:
                self._check_for_tools(sample)
                
            return sample
            
        except Exception as e:
            logger.debug(f"Couldn't sample {filepath}: {str(e)}")
            return f"Could not read file: {str(e)}"
    
    def _check_for_tools(self, text):
        """Look for mentions of bioinformatics tools in text"""
        matches = TOOL_PATTERN.findall(str(text))
        for match in matches:
            self.stats["potential_tools"][match.lower()] += 1
    
    def _update_extremes(self, file_entry):
        """Update the lists of largest/oldest/newest files"""
        # Add to largest files
        self.stats["largest_files"].append(file_entry.copy())
        self.stats["largest_files"].sort(key=lambda x: x["size"], reverse=True)
        self.stats["largest_files"] = self.stats["largest_files"][:10]
        
        # Add to oldest files
        self.stats["oldest_files"].append(file_entry.copy())
        self.stats["oldest_files"].sort(key=lambda x: x["modified"])
        self.stats["oldest_files"] = self.stats["oldest_files"][:10]
        
        # Add to newest files
        self.stats["newest_files"].append(file_entry.copy())
        self.stats["newest_files"].sort(key=lambda x: x["modified"], reverse=True)
        self.stats["newest_files"] = self.stats["newest_files"][:10]
    
    def _dump_memory(self, dump_id=None):
        """Dump collected data to disk to free memory"""
        if dump_id is None:
            self.dump_counter += 1
            dump_id = f"dump_{self.dump_counter:04d}"
            
        dump_dir = os.path.join(self.output_dir, "data_dumps")
        dump_path = os.path.join(dump_dir, dump_id)
        
        # Save data to disk
        with open(f"{dump_path}_dirs.json", 'w') as f:
            json.dump(self.directory_summaries, f)
        
        with open(f"{dump_path}_files.json", 'w') as f:
            json.dump(self.interesting_files, f)
        
        # Update stats file with every dump
        with open(os.path.join(self.output_dir, "stats.json"), 'w') as f:
            json.dump(self.stats, f)
        
        logger.info(f"Memory dump completed: {dump_id}")
        logger.info(f"Dumped {len(self.directory_summaries)} directories and {len(self.interesting_files)} files")
        
        # Clear memory
        self.directory_summaries = {}
        self.interesting_files = []
    
    def _finalize_crawl(self):
        """Process collected data after crawl is complete"""
        # Save any remaining data
        if self.directory_summaries or self.interesting_files:
            self._dump_memory("final")
        
        # Load all the dumps to generate final outputs
        all_dirs = {}
        all_files = []
        
        # Load directory summaries
        dump_dir = os.path.join(self.output_dir, "data_dumps")
        for filename in os.listdir(dump_dir):
            if filename.endswith('_dirs.json'):
                with open(os.path.join(dump_dir, filename), 'r') as f:
                    dirs_data = json.load(f)
                    all_dirs.update(dirs_data)
        
        # Load interesting files (up to a reasonable limit)
        file_count = 0
        for filename in os.listdir(dump_dir):
            if filename.endswith('_files.json'):
                with open(os.path.join(dump_dir, filename), 'r') as f:
                    files_data = json.load(f)
                    # Limit the number of files to process
                    if file_count < 10000:
                        to_add = min(10000 - file_count, len(files_data))
                        all_files.extend(files_data[:to_add])
                        file_count += to_add
        
        # Generate summary document
        summary = self.get_summary()
        with open(os.path.join(self.output_dir, "summary.md"), 'w') as f:
            f.write(summary)
        
        # Generate a smaller LLM-ready version
        llm_input = self.get_llm_input(all_dirs, all_files)
        with open(os.path.join(self.output_dir, "llm_input.md"), 'w') as f:
            f.write(llm_input)
        
        logger.info(f"All outputs saved to {self.output_dir}")
    
    def get_summary(self):
        """Generate a detailed markdown summary of findings"""
        summary = [
            "# Bioinformatics Data Crawler Summary",
            "",
            f"**Root directory:** {self.root_dir}",
            f"**Crawl completed:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Overview Statistics",
            f"- Total directories: {self.stats['total_dirs']:,}",
            f"- Total files: {self.stats['total_files']:,}",
            f"- Total size: {self._format_size(self.stats['total_size'])}",
            "",
            "## File Types",
            ""
        ]
        
        # Add file types table
        summary.append("| Extension | Count | Percentage |")
        summary.append("|-----------|-------|------------|")
        for ext, count in sorted(self.stats["file_types"].items(), key=lambda x: x[1], reverse=True)[:20]:
            pct = (count / max(1, self.stats["total_files"])) * 100
            summary.append(f"| {ext or 'No extension'} | {count:,} | {pct:.1f}% |")
        
        # Add file ages
        summary.append("")
        summary.append("## File Age Distribution")
        summary.append("")
        summary.append("| Age | Count | Percentage |")
        summary.append("|-----|-------|------------|")
        for age, count in sorted(self.stats["file_ages"].items()):
            pct = (count / max(1, self.stats["total_files"])) * 100
            summary.append(f"| {age.replace('_', ' ').title()} | {count:,} | {pct:.1f}% |")
        
        # Add tools found
        summary.append("")
        summary.append("## Potential Bioinformatics Tools")
        summary.append("")
        summary.append("| Tool | Count |")
        summary.append("|------|-------|")
        for tool, count in sorted(self.stats["potential_tools"].items(), key=lambda x: x[1], reverse=True)[:20]:
            summary.append(f"| {tool} | {count:,} |")
        
        # Add largest files
        summary.append("")
        summary.append("## Largest Files")
        summary.append("")
        summary.append("| Path | Size | Modified |")
        summary.append("|------|------|----------|")
        for file in self.stats["largest_files"]:
            summary.append(f"| {file['path']} | {file['size_human']} | {file['modified']} |")
        
        # Add newest files
        summary.append("")
        summary.append("## Most Recent Files")
        summary.append("")
        summary.append("| Path | Size | Modified |")
        summary.append("|------|------|----------|")
        for file in self.stats["newest_files"]:
            summary.append(f"| {file['path']} | {file['size_human']} | {file['modified']} |")
        
        # Add oldest files
        summary.append("")
        summary.append("## Oldest Files")
        summary.append("")
        summary.append("| Path | Size | Modified |")
        summary.append("|------|------|----------|")
        for file in self.stats["oldest_files"]:
            summary.append(f"| {file['path']} | {file['size_human']} | {file['modified']} |")
        
        # Interesting bioinformatics patterns
        summary.append("")
        summary.append("## Interesting Findings")
        summary.append("")
        
        # Look for somatic variant calling
        if any("somatic" in tool.lower() for tool in self.stats["potential_tools"]):
            summary.append("- Found evidence of somatic variant calling")
        
        # Look for typical cancer/somatic analysis patterns
        if any(x in self.stats["potential_tools"] for x in ["tumor", "normal", "mutect", "strelka", "varscan"]):
            summary.append("- Found cancer genomics/somatic analysis tools")
        
        # Look for RNA-seq
        if any(x in self.stats["potential_tools"] for x in ["rnaseq", "star", "rsem", "salmon", "kallisto"]):
            summary.append("- Found RNA-seq analysis tools")
        
        # Look for workflow systems
        if any(x in self.stats["potential_tools"] for x in ["nextflow", "snakemake", "wdl", "cwl"]):
            workflow = next((x for x in ["nextflow", "snakemake", "wdl", "cwl"] 
                           if x in self.stats["potential_tools"]), None)
            if workflow:
                summary.append(f"- Workflow system detected: {workflow}")
        
        return "\n".join(summary)
    
    def get_llm_input(self, dirs_data=None, files_data=None):
        """Generate a compact version for LLM analysis"""
        # If no data is provided, use the current in-memory data
        dirs_data = dirs_data or self.directory_summaries
        files_data = files_data or self.interesting_files
        
        # Create a distilled version for the LLM that focuses on patterns
        llm_input = [
            "# Bioinformatics Dataset Analysis",
            "",
            "## Task",
            "You are analyzing a bioinformatics dataset. Below are key details about the filesystem:",
            "",
            f"- Total directories: {self.stats['total_dirs']:,}",
            f"- Total files: {self.stats['total_files']:,}",
            f"- Total size: {self._format_size(self.stats['total_size'])}",
            "",
            "## Most common file types:",
            ""
        ]
        
        # Add file types (limited)
        for ext, count in sorted(self.stats["file_types"].items(), key=lambda x: x[1], reverse=True)[:15]:
            pct = (count / max(1, self.stats["total_files"])) * 100
            llm_input.append(f"- {ext or 'No extension'}: {count:,} files ({pct:.1f}%)")
        
        # Add tools found
        llm_input.append("")
        llm_input.append("## Potential bioinformatics tools identified:")
        llm_input.append("")
        for tool, count in sorted(self.stats["potential_tools"].items(), key=lambda x: x[1], reverse=True)[:15]:
            llm_input.append(f"- {tool}: {count:,} mentions")
        
        # Add interesting file samples
        llm_input.append("")
        llm_input.append("## Samples from interesting files:")
        llm_input.append("")
        
        # Include samples from diverse file types
        added_extensions = set()
        max_samples = 10
        
        # Sort files by the frequency of their extension types
        sorted_files = sorted(
            files_data, 
            key=lambda x: self.stats["file_types"].get(x.get("extension", ""), 0),
            reverse=True
        )
        
        for file in sorted_files:
            ext = file.get("extension", "")
            
            # Ensure variety of file types
            if ext not in added_extensions and len(added_extensions) < max_samples:
                sample = file.get("sample", "")
                if sample:
                    # Truncate very long samples
                    if len(sample) > 500:
                        sample = sample[:500] + "...[truncated]"
                    
                    llm_input.append(f"### {file['path']}")
                    llm_input.append("```")
                    llm_input.append(sample)
                    llm_input.append("```")
                    llm_input.append("")
                    
                    added_extensions.add(ext)
        
        # Add analysis request
        llm_input.append("")
        llm_input.append("## Analysis Request:")
        llm_input.append("")
        llm_input.append("Based on the information above, please analyze this dataset and provide:")
        llm_input.append("")
        llm_input.append("1. What bioinformatics workflows and tools were likely used")
        llm_input.append("2. What kind of analysis was being performed (e.g., somatic variant calling, RNA-seq)")
        llm_input.append("3. The probable inputs and outputs of the analysis")
        llm_input.append("4. Any sequencing technologies or platforms that might have been used")
        llm_input.append("5. Recommendations for further exploration of this dataset")
        
        return "\n".join(llm_input)
    
    @staticmethod
    def _format_size(size_bytes):
        """Format file size in human-readable format"""
        if size_bytes < 1024:
            return f"{size_bytes} B"
        elif size_bytes < 1024 * 1024:
            return f"{size_bytes/1024:.1f} KB"
        elif size_bytes < 1024 * 1024 * 1024:
            return f"{size_bytes/(1024*1024):.1f} MB"
        elif size_bytes < 1024 * 1024 * 1024 * 1024:
            return f"{size_bytes/(1024*1024*1024):.1f} GB"
        else:
            return f"{size_bytes/(1024*1024*1024*1024):.1f} TB"


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Simplified Bioinformatics File System Crawler')
    parser.add_argument('root_dir', help='Root directory to crawl')
    parser.add_argument('--output-dir', '-o', default='crawler_output', 
                        help='Directory to save outputs')
    parser.add_argument('--max-file-size', '-m', type=int, default=10*1024*1024,
                        help='Maximum file size to sample (bytes)')
    parser.add_argument('--sample-size', '-s', type=int, default=4*1024,
                        help='Size of content sample (bytes)')
    parser.add_argument('--max-files', '-f', type=int, default=5,
                        help='Maximum files to sample per directory')
    parser.add_argument('--memory-dump-interval', '-d', type=int, default=5000,
                        help='Number of directories between memory dumps')
    
    args = parser.parse_args()
    
    # Validate root directory
    if not os.path.isdir(args.root_dir):
        print(f"Error: Root directory not found: {args.root_dir}")
        print("Make sure to provide an absolute path (e.g., '/home/user/data') or a valid relative path (e.g., './data')")
        sys.exit(1)
    
    try:
        start_time = time.time()
        
        # Start crawler
        crawler = SimpleBioCrawler(
            args.root_dir,
            output_dir=args.output_dir,
            max_file_size=args.max_file_size,
            sample_size=args.sample_size,
            max_files_per_dir=args.max_files,
            memory_dump_interval=args.memory_dump_interval
        )
        
        # Start the actual crawl
        crawler.start_time = start_time  # Store start time for progress reporting
        summary = crawler.crawl()
        
        print(f"Crawl complete. Outputs saved to {args.output_dir}")
        print("Summary file: summary.md")
        print("LLM input file: llm_input.md")
        
    except KeyboardInterrupt:
        print("\nCrawl interrupted. Partial results may have been saved.")
    except Exception as e:
        import traceback
        print(f"Error: {str(e)}")
        traceback.print_exc()