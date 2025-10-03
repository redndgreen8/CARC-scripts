#!/usr/bin/env python3
"""
Bioinformatics File System Crawler

This script crawls a file system containing bioinformatics data, extracts metadata and 
content samples, and generates summaries suitable for analysis by an LLM.

It's specifically designed for understanding a bioinformatics storage partition,
with special handling for common bioinformatics file formats.
"""
import sys
import os
import time
import json
import re
import mimetypes
import hashlib
import argparse
import gzip
import zipfile
import tarfile
from datetime import datetime
from collections import Counter, defaultdict
import concurrent.futures
import logging

HAVE_PSUTIL = False
HAVE_TQDM = False


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

# Additional patterns for identifying bioinformatics tools
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

class BioinformaticsCrawler:
    def __init__(self, root_dir, output_dir="crawler_output", 
                 max_file_size=10*1024*1024, sample_size=4*1024,
                 max_files_per_dir=5, max_workers=4, 
                 max_dir_workers=2, checkpoint_interval=5000,
                 memory_limit_gb=4, enable_progress=True, 
                 resume=False):
        """
        Initialize the crawler.
        
        Args:
            root_dir: Root directory to start crawling from
            output_dir: Directory to save output files
            max_file_size: Maximum file size to sample content from (bytes)
            sample_size: Size of content sample to extract (bytes)
            max_files_per_dir: Maximum number of files to sample per directory
            max_workers: Maximum number of worker threads for file processing
            max_dir_workers: Maximum number of worker threads for directory traversal
            checkpoint_interval: Number of directories between checkpoint saves
            memory_limit_gb: Memory limit in GB before flushing data to disk
            enable_progress: Enable progress tracking with tqdm
            resume: Whether to resume from previous checkpoint
        """
        # Validate root directory exists
        if not os.path.isdir(root_dir):
            raise ValueError(f"Root directory does not exist or is not accessible: {root_dir}")
            
        self.root_dir = os.path.abspath(root_dir)
        self.output_dir = output_dir
        self.max_file_size = max_file_size
        self.sample_size = sample_size
        self.max_files_per_dir = max_files_per_dir
        self.max_workers = max_workers
        self.max_dir_workers = max_dir_workers
        self.checkpoint_interval = checkpoint_interval
        self.memory_limit_bytes = memory_limit_gb * 1024 * 1024 * 1024  # Convert to bytes
        self.enable_progress = enable_progress and HAVE_TQDM
        self.resume = resume
        
        # Statistics and aggregation
        self.stats = {
            "total_dirs": 0,
            "total_files": 0,
            "total_size": 0,
            "processed_size": 0,
            "file_types": Counter(),
            "file_ages": defaultdict(int),
            "potential_tools": Counter(),
            "largest_files": [],
            "oldest_files": [],
            "newest_files": []
        }
        
        # Runtime tracking
        self.start_time = None
        self.checkpoint_time = None
        self.last_flush_time = None
        self.flush_count = 0
        
        # Initialize attributes that might be accessed later
        self.estimated_total_size = None
        self.estimated_dir_count = None
        self.pbar = None
        
        # Data stores with memory management
        self.directory_summaries = {}
        self.file_samples = {}
        self.interesting_files = []
        self.pending_dirs = queue.Queue()
        
        # Control flags
        self.running = False
        self.interrupted = False
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Setup checkpoint directory
        self.checkpoint_dir = os.path.join(output_dir, "checkpoints")
        os.makedirs(self.checkpoint_dir, exist_ok=True)
        
        # Register signal handlers for graceful shutdown
        try:
            signal.signal(signal.SIGINT, self._handle_interrupt)
            signal.signal(signal.SIGTERM, self._handle_interrupt)
        except (ValueError, AttributeError):
            # This can happen in environments where signals aren't supported
            logger.warning("Signal handlers could not be registered")
        
        # Try to load from checkpoint if resuming
        if resume:
            self._load_checkpoint()
        
    def _handle_interrupt(self, signum, frame):
        """Signal handler for graceful shutdown"""
        logger.warning(f"Received signal {signum}, initiating graceful shutdown...")
        self.interrupted = True
        if self.pbar:
            self.pbar.write("Interrupted! Saving checkpoint before exiting...")
        self._save_checkpoint("interrupt")
    
    def _monitor_memory(self):
        """Monitor memory usage and flush data if needed"""
        if not HAVE_PSUTIL:
            # If psutil isn't available, use a simpler heuristic
            if len(self.directory_summaries) > 10000 or len(self.interesting_files) > 1000:
                self._flush_data_to_disk()
                return True
            return False
            
        process = psutil.Process(os.getpid())
        mem_info = process.memory_info()
        
        # Check if memory usage exceeds limit
        if mem_info.rss > self.memory_limit_bytes:
            self._flush_data_to_disk()
            return True
        return False
    
    def _flush_data_to_disk(self):
        """Flush directory summaries to disk to free memory"""
        self.flush_count += 1
        flush_id = f"flush_{self.flush_count:04d}"
        flush_path = os.path.join(self.checkpoint_dir, flush_id)
        
        # Save data to disk
        with open(f"{flush_path}_dirs.pickle", 'wb') as f:
            pickle.dump(self.directory_summaries, f)
        
        with open(f"{flush_path}_files.pickle", 'wb') as f:
            pickle.dump(self.interesting_files, f)
        
        # Update metadata
        with open(f"{flush_path}_meta.json", 'w') as f:
            json.dump({
                'directories': len(self.directory_summaries),
                'files': len(self.interesting_files),
                'timestamp': datetime.now().isoformat()
            }, f)
        
        # Clear memory
        logger.info(f"Flushed data to disk (flush #{self.flush_count}): "
                    f"{len(self.directory_summaries)} directories, "
                    f"{len(self.interesting_files)} interesting files")
        
        # Keep the stats
        self.directory_summaries = {}
        self.file_samples = {}
        self.interesting_files = []
        
        # Force garbage collection
        import gc
        gc.collect()
        
        self.last_flush_time = time.time()
    
    def _save_checkpoint(self, checkpoint_id=None):
        """Save a checkpoint of the current state"""
        if checkpoint_id is None:
            checkpoint_id = f"checkpoint_{self.stats['total_dirs']:08d}"
        
        checkpoint_path = os.path.join(self.checkpoint_dir, checkpoint_id)
        
        # Flush any pending data first
        self._flush_data_to_disk()
        
        # Save overall stats
        with open(f"{checkpoint_path}_stats.json", 'w') as f:
            json.dump(self.stats, f, indent=2)
        
        logger.info(f"Checkpoint saved: {checkpoint_id}")
    
    def _load_checkpoint(self):
        """Load the latest checkpoint if available"""
        checkpoints = []
        for filename in os.listdir(self.checkpoint_dir):
            if filename.endswith('_stats.json'):
                checkpoint_id = filename.replace('_stats.json', '')
                checkpoints.append(checkpoint_id)
        
        if not checkpoints:
            logger.warning("No checkpoints found, starting fresh crawl")
            return False
        
        # Find latest checkpoint
        latest = sorted(checkpoints)[-1]
        checkpoint_path = os.path.join(self.checkpoint_dir, latest)
        
        # Load stats
        with open(f"{checkpoint_path}_stats.json", 'r') as f:
            self.stats = json.load(f)
        
        logger.info(f"Resumed from checkpoint: {latest}")
        logger.info(f"Loaded state: {self.stats['total_dirs']} directories, "
                    f"{self.stats['total_files']} files processed")
        
        return True
    
    def _estimate_total_size(self):
        """
        Make an initial estimate of total size and file count
        to provide progress information
        """
        # Use sampling to estimate total size if not resuming
        if self.estimated_total_size is None:
            # Sample a few directories to estimate
            sample_sizes = []
            sample_dirs = []
            
            logger.info("Estimating total size (sampling directories)...")
            
            # Quick sampling of directories
            for root, dirs, files in os.walk(self.root_dir, topdown=True):
                if len(sample_dirs) >= 20:  # Sample 20 directories
                    break
                
                try:
                    total = 0
                    for f in files:
                        try:
                            total += os.path.getsize(os.path.join(root, f))
                        except (OSError, IOError):
                            pass
                    
                    if total > 0:
                        sample_sizes.append((len(files), total))
                        sample_dirs.append(root)
                except Exception as e:
                    continue
            
            # Estimate based on average file size
            if sample_sizes:
                total_files = sum(f for f, _ in sample_sizes)
                total_size = sum(s for _, s in sample_sizes)
                avg_file_size = total_size / max(1, total_files)
                
                # Get total directory count estimate 
                # (this is fast even on large filesystems)
                est_dir_count = 0
                for _ in os.walk(self.root_dir, topdown=True):
                    est_dir_count += 1
                    if est_dir_count > 10000:  # Cap at 10000 for speed
                        est_dir_count = est_dir_count * 2  # Rough estimate
                        break
                
                # Very rough file count estimate based on directories
                est_file_count = est_dir_count * (total_files / max(1, len(sample_dirs)))
                self.estimated_total_size = est_file_count * avg_file_size
                
                logger.info(f"Size estimate: ~{self._format_size(self.estimated_total_size)} "
                            f"(~{int(est_file_count):,} files in ~{est_dir_count:,} directories)")
                
                return self.estimated_total_size
            
            # Fallback if estimation fails
            self.estimated_total_size = 1_000_000_000_000  # 1 TB default guess
            logger.warning("Could not estimate size, using default: 1 TB")
            return self.estimated_total_size
    
    def crawl(self):
        """Main crawl method that orchestrates the filesystem traversal"""
        self.start_time = time.time()
        self.checkpoint_time = self.start_time
        self.last_flush_time = self.start_time
        self.running = True
        
        logger.info(f"Starting crawl of {self.root_dir}")
        
        # Enable memory tracking if available
        if tracemalloc:
            tracemalloc.start()
        
        # Estimate total size for progress reporting
        if self.enable_progress:
            total_size = self._estimate_total_size()
            self.pbar = tqdm(
                total=total_size, 
                unit='B', 
                unit_scale=True,
                desc="Crawling filesystem"
            )
        
        try:
            # Setup directory traversal workers
            if self.max_dir_workers > 1:
                self._parallel_walk()
            else:
                self._sequential_walk()
                
        except KeyboardInterrupt:
            logger.warning("Crawl interrupted by user")
            if not self.interrupted:  # If not already handled by signal
                self._save_checkpoint("interrupt")
        except Exception as e:
            logger.error(f"Error during crawl: {str(e)}", exc_info=True)
            self._save_checkpoint("error")
        finally:
            # Clean up
            self.running = False
            if self.pbar:
                self.pbar.close()
            if tracemalloc:
                tracemalloc.stop()
        
        # Process collected data
        self._finalize_crawl()
        
        # Report statistics
        elapsed = time.time() - self.start_time
        hrs, remainder = divmod(elapsed, 3600)
        mins, secs = divmod(remainder, 60)
        
        logger.info(f"Crawl completed in {int(hrs)}h {int(mins)}m {int(secs)}s")
        logger.info(f"Found {self.stats['total_files']:,} files in {self.stats['total_dirs']:,} directories")
        logger.info(f"Total size: {self._format_size(self.stats['total_size'])}")
        logger.info(f"Memory flushes: {self.flush_count}")
        
        return self.get_summary()
    
    def _sequential_walk(self):
        """Sequential directory walk with memory management"""
        # Start with root directory
        self.pending_dirs.put((self.root_dir, [], []))
        
        # Process one directory at a time
        while not self.interrupted and not self.pending_dirs.empty():
            current_dir, _, _ = self.pending_dirs.get()
            
            try:
                # List directory contents
                entries = os.scandir(current_dir)
                subdirs = []
                files = []
                
                # Separate files and directories
                for entry in entries:
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            if not entry.name.startswith('.'):
                                subdirs.append(entry.name)
                        elif entry.is_file(follow_symlinks=False):
                            if not entry.name.startswith('.'):
                                files.append(entry.name)
                    except OSError:
                        continue
                
                # Process this directory
                self._process_directory(current_dir, subdirs, files)
                
                # Add subdirectories to the queue
                for d in subdirs:
                    dir_path = os.path.join(current_dir, d)
                    self.pending_dirs.put((dir_path, [], []))
                
                # Periodic checkpoint
                if self.stats['total_dirs'] % self.checkpoint_interval == 0:
                    self._save_checkpoint()
                
                # Check memory usage periodically
                if (time.time() - self.last_flush_time) > 300:  # Check every 5 minutes
                    self._monitor_memory()
                
            except OSError as e:
                logger.error(f"Error accessing {current_dir}: {str(e)}")
            except Exception as e:
                logger.error(f"Error processing {current_dir}: {str(e)}")
    
    def _parallel_walk(self):
        """Parallel directory traversal with worker threads"""
        # Start with root directory
        self.pending_dirs.put((self.root_dir, None))
        
        # Create stop event for signaling termination
        stop_event = threading.Event()
        
        # Create worker threads for directory traversal
        dir_workers = []
        for _ in range(self.max_dir_workers):
            worker = threading.Thread(
                target=self._dir_worker,
                args=(stop_event,)
            )
            worker.daemon = True
            worker.start()
            dir_workers.append(worker)
        
        # Monitor the queue and workers
        try:
            while any(w.is_alive() for w in dir_workers) and not self.interrupted:
                # Checkpoint periodically
                if (time.time() - self.checkpoint_time) > 600:  # Every 10 minutes
                    self._save_checkpoint()
                    self.checkpoint_time = time.time()
                
                # Check memory usage
                self._monitor_memory()
                
                # Print statistics periodically
                if self.pbar:
                    self.pbar.set_postfix({
                        'dirs': f"{self.stats['total_dirs']:,}",
                        'files': f"{self.stats['total_files']:,}"
                    })
                
                time.sleep(5)  # Check status every 5 seconds
                
                # Break if all queues are empty and workers are idle
                if self.pending_dirs.empty() and all(not getattr(w, 'is_processing', False) for w in dir_workers):
                    logger.info("No more directories to process, stopping workers")
                    break
                    
        except Exception as e:
            logger.error(f"Error in main thread: {str(e)}")
        finally:
            # Signal workers to stop
            stop_event.set()
            
            # Wait for workers to finish (with timeout)
            for worker in dir_workers:
                worker.join(timeout=5)
    
    def _dir_worker(self, stop_event):
        """Worker thread for directory traversal"""
        thread_id = threading.get_ident()
        
        while not stop_event.is_set() and not self.interrupted:
            try:
                # Get next directory with timeout
                try:
                    current_dir, parent_id = self.pending_dirs.get(timeout=2)
                    setattr(threading.current_thread(), 'is_processing', True)
                except queue.Empty:
                    setattr(threading.current_thread(), 'is_processing', False)
                    continue
                
                # List directory contents
                try:
                    entries = list(os.scandir(current_dir))
                    subdirs = []
                    files = []
                    
                    # Separate files and directories
                    for entry in entries:
                        try:
                            if entry.is_dir(follow_symlinks=False):
                                if not entry.name.startswith('.'):
                                    subdirs.append(entry.name)
                            elif entry.is_file(follow_symlinks=False):
                                if not entry.name.startswith('.'):
                                    files.append(entry.name)
                        except OSError:
                            continue
                    
                    # Process this directory
                    self._process_directory(current_dir, subdirs, files)
                    
                    # Add subdirectories to the queue
                    for d in subdirs:
                        dir_path = os.path.join(current_dir, d)
                        self.pending_dirs.put((dir_path, thread_id))
                    
                except OSError as e:
                    logger.error(f"Error accessing {current_dir}: {str(e)}")
                except Exception as e:
                    logger.error(f"Error in directory worker: {str(e)}")
                finally:
                    # Mark task as done
                    self.pending_dirs.task_done()
                    setattr(threading.current_thread(), 'is_processing', False)
                
            except Exception as e:
                logger.error(f"Unexpected error in directory worker: {str(e)}")
        
        logger.debug(f"Directory worker {thread_id} stopped")
    
    def _process_directory(self, current_dir, dirs, files):
        """Process a single directory and its files"""
        rel_path = os.path.relpath(current_dir, self.root_dir)
        self.stats["total_dirs"] += 1
        
        # Update progress reporting
        if self.stats["total_dirs"] % 100 == 0:
            dirs_per_sec = self.stats["total_dirs"] / max(1, time.time() - self.start_time)
            
            # Calculate ETA
            if dirs_per_sec > 0:
                # Make a very rough estimate of total directories based on what we've seen
                if not self.estimated_dir_count:
                    # Rough heuristic: directories often have similar depth patterns
                    depth = len(rel_path.split(os.sep))
                    if depth > 1:
                        # Estimate more directories at deeper levels
                        self.estimated_dir_count = self.stats["total_dirs"] * (2 ** min(depth, 5))
                    else:
                        self.estimated_dir_count = self.stats["total_dirs"] * 10
                
                remaining = max(0, self.estimated_dir_count - self.stats["total_dirs"])
                eta_seconds = remaining / dirs_per_sec
                eta = str(timedelta(seconds=int(eta_seconds)))
            else:
                eta = "unknown"
            
            logger.info(f"Processed {self.stats['total_dirs']:,} directories "
                        f"({dirs_per_sec:.1f}/s), "
                        f"{self.stats['total_files']:,} files, "
                        f"ETA: {eta}")
        
        # Update progress bar
        if self.pbar:
            self.pbar.update(self.stats["processed_size"] - self.pbar.n)
            
            # Update description with more information
            self.pbar.set_description(
                f"Crawling: {self.stats['total_dirs']:,} dirs, {self.stats['total_files']:,} files"
            )
        
        # Skip hidden directories
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        
        # Extract sample from files in this directory
        file_entries = []
        file_count = 0
        dir_size = 0
        
        # Special case for very large directories
        if len(files) > 1000:
            logger.info(f"Large directory with {len(files)} files: {current_dir}")
            # Sample a subset of files for very large directories
            files_sample = files[:100] + random.sample(files[100:], min(900, len(files)-100))
            files = files_sample
        
        # Process files in parallel for large directories
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Create a map of future to filename for tracking
            future_to_file = {
                executor.submit(self._process_file, current_dir, filename): filename
                for filename in files if not filename.startswith('.')
            }
            
            for future in concurrent.futures.as_completed(future_to_file):
                filename = future_to_file[future]
                try:
                    result = future.result()
                    if result:
                        file_entries.append(result)
                        file_count += 1
                        # Track directory size from file results
                        if 'size' in result:
                            dir_size += result['size']
                except Exception as e:
                    logger.error(f"Error processing {os.path.join(current_dir, filename)}: {str(e)}")
        
        # Update processed size for progress tracking
        self.stats["processed_size"] += dir_size
        
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
        
        # Check memory usage periodically
        if self.stats["total_dirs"] % 50 == 0:
            self._monitor_memory()
    
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
            
            # Keep track of largest and oldest/newest files
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
                with gzip.open(filepath, 'rt', errors='replace') as f:
                    sample = f.read(sample_size)
            elif extension.endswith(('.zip')):
                try:
                    with zipfile.ZipFile(filepath) as z:
                        with z.open(z.namelist()[0]) as f:
                            sample = f.read(sample_size).decode('utf-8', errors='replace')
                except:
                    sample = "ZIP file, but could not extract content"
            elif extension.endswith(('.tar', '.tar.gz', '.tgz')):
                try:
                    with tarfile.open(filepath) as t:
                        for member in t.getmembers():
                            if member.isfile():
                                with t.extractfile(member) as f:
                                    sample = f.read(sample_size).decode('utf-8', errors='replace')
                                    break
                except:
                    sample = "TAR file, but could not extract content"
            elif extension in ('.sam', '.bam', '.cram'):
                # For binary files, just return a placeholder
                sample = f"Binary alignment file ({extension})"
            elif extension in ('.vcf', '.vcf.gz'):
                # Try to extract VCF header
                if extension == '.vcf':
                    with open(filepath, 'r', errors='replace') as f:
                        lines = []
                        for line in f:
                            if line.startswith('#'):
                                lines.append(line.strip())
                            else:
                                # Add the first data line
                                lines.append(line.strip())
                                break
                        sample = '\n'.join(lines)
                else:  # .vcf.gz
                    with gzip.open(filepath, 'rt', errors='replace') as f:
                        lines = []
                        for line in f:
                            if line.startswith('#'):
                                lines.append(line.strip())
                            else:
                                # Add the first data line
                                lines.append(line.strip())
                                break
                        sample = '\n'.join(lines)
            else:
                # Regular text file
                with open(filepath, 'r', errors='replace') as f:
                    sample = f.read(sample_size)
            
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
    
    def _finalize_crawl(self):
        """Process collected data after crawl is complete"""
        # Save outputs
        self._save_json("file_stats.json", self.stats)
        self._save_json("directory_summaries.json", self.directory_summaries)
        self._save_json("interesting_files.json", self.interesting_files)
        
        # Generate summary document
        summary = self.get_summary()
        self._save_text("summary.md", summary)
        
        # Generate a smaller LLM-ready version
        llm_input = self.get_llm_input()
        self._save_text("llm_input.md", llm_input)
        
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
            pct = (count / self.stats["total_files"]) * 100
            summary.append(f"| {ext or 'No extension'} | {count:,} | {pct:.1f}% |")
        
        # Add file ages
        summary.append("")
        summary.append("## File Age Distribution")
        summary.append("")
        summary.append("| Age | Count | Percentage |")
        summary.append("|-----|-------|------------|")
        for age, count in sorted(self.stats["file_ages"].items()):
            pct = (count / self.stats["total_files"]) * 100
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
            workflow = next(x for x in ["nextflow", "snakemake", "wdl", "cwl"] 
                           if x in self.stats["potential_tools"])
            summary.append(f"- Workflow system detected: {workflow}")
        
        return "\n".join(summary)
    
    def get_llm_input(self):
        """Generate a compact version for LLM analysis"""
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
            pct = (count / self.stats["total_files"]) * 100
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
        samples_count = 0
        max_samples = 10
        
        for file in sorted(self.interesting_files, key=lambda x: self.stats["file_types"].get(x["extension"], 0), reverse=True):
            ext = file["extension"]
            
            # Ensure variety of file types
            if ext not in added_extensions and samples_count < max_samples:
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
                    samples_count += 1
        
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
    
    def _save_json(self, filename, data):
        """Save data to a JSON file"""
        filepath = os.path.join(self.output_dir, filename)
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
    
    def _save_text(self, filename, text):
        """Save text to a file"""
        filepath = os.path.join(self.output_dir, filename)
        with open(filepath, 'w') as f:
            f.write(text)
    
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
    parser = argparse.ArgumentParser(description='Bioinformatics File System Crawler')
    parser.add_argument('root_dir', help='Root directory to crawl')
    parser.add_argument('--output-dir', '-o', default='crawler_output', 
                        help='Directory to save outputs')
    parser.add_argument('--max-file-size', '-m', type=int, default=10*1024*1024,
                        help='Maximum file size to sample (bytes)')
    parser.add_argument('--sample-size', '-s', type=int, default=4*1024,
                        help='Size of content sample (bytes)')
    parser.add_argument('--max-files', '-f', type=int, default=5,
                        help='Maximum files to sample per directory')
    parser.add_argument('--workers', '-w', type=int, default=4,
                        help='Number of worker threads for file processing')
    parser.add_argument('--dir-workers', '-d', type=int, default=2,
                        help='Number of worker threads for directory traversal')
    parser.add_argument('--memory-limit', '-l', type=float, default=4.0,
                        help='Memory limit in GB before flushing to disk')
    parser.add_argument('--checkpoint-interval', '-c', type=int, default=5000,
                        help='Number of directories between checkpoints')
    parser.add_argument('--resume', '-r', action='store_true',
                        help='Resume from latest checkpoint if available')
    parser.add_argument('--no-progress', action='store_true',
                        help='Disable progress bar')
    
    args = parser.parse_args()
    
    # Validate root directory
    if not os.path.isdir(args.root_dir):
        logger.error(f"Root directory does not exist or is not accessible: {args.root_dir}")
        print(f"Error: Root directory not found: {args.root_dir}")
        print("Make sure to provide an absolute path (e.g., '/home/user/data') or a valid relative path (e.g., './data')")
        sys.exit(1)
    
    # Import missing module warnings
    if not HAVE_PSUTIL:
        print("Warning: For better memory management, install psutil: pip install psutil")
    if not HAVE_TQDM:
        print("Warning: For progress reporting, install tqdm: pip install tqdm")
    
    # Missing random import check
    try:
        import random
    except ImportError:
        print("Error: Unable to import 'random' module which should be part of standard library")
        sys.exit(1)
    
    # Recommended settings for very large filesystems
    if os.path.isdir(args.root_dir):
        total_size_estimate = 0
        try:
            # Get filesystem stats
            fs_stat = os.statvfs(args.root_dir)
            total_size_estimate = fs_stat.f_blocks * fs_stat.f_frsize
            used_size_estimate = (fs_stat.f_blocks - fs_stat.f_bfree) * fs_stat.f_frsize
            
            # Warn and adjust settings for large filesystems
            if total_size_estimate > 100 * 1024**4:  # > 100TB
                logger.warning(f"Large filesystem detected: {total_size_estimate / 1024**4:.1f} TB")
                logger.warning("Adjusting settings for optimal performance...")
                
                # Adjust workers based on cores
                cpu_count = os.cpu_count() or 4
                if args.workers == 4:  # Only if default wasn't changed
                    args.workers = min(max(cpu_count // 2, 2), 16)
                if args.dir_workers == 2:  # Only if default wasn't changed
                    args.dir_workers = min(max(cpu_count // 4, 1), 8)
                
                # Adjust memory limit based on system memory
                if HAVE_PSUTIL:
                    system_memory_gb = psutil.virtual_memory().total / 1024**3
                    if args.memory_limit == 4.0:  # Only if default wasn't changed
                        # Use 25% of system memory, but at least 4GB and at most 32GB
                        args.memory_limit = min(max(system_memory_gb * 0.25, 4), 32)
                
                logger.info(f"Using {args.workers} file workers, {args.dir_workers} directory workers")
                logger.info(f"Memory limit set to {args.memory_limit:.1f} GB")
        except Exception as e:
            logger.error(f"Error checking filesystem size: {str(e)}")
    
    try:
        crawler = BioinformaticsCrawler(
            args.root_dir,
            output_dir=args.output_dir,
            max_file_size=args.max_file_size,
            sample_size=args.sample_size,
            max_files_per_dir=args.max_files,
            max_workers=args.workers,
            max_dir_workers=args.dir_workers,
            checkpoint_interval=args.checkpoint_interval,
            memory_limit_gb=args.memory_limit,
            enable_progress=not args.no_progress,
            resume=args.resume
        )
        
        summary = crawler.crawl()
        print(f"Crawl complete. Outputs saved to {args.output_dir}")
        print("Summary file: summary.md")
        print("LLM input file: llm_input.md")
    except KeyboardInterrupt:
        print("\nCrawl interrupted. Partial results saved to checkpoints directory.")
        print("Run with --resume to continue from where you left off.")
    except Exception as e:
        import traceback
        print(f"Error: {str(e)}")
        traceback.print_exc()
        print("\nTry running with --resume to continue from last checkpoint.")
        sys.exit(1)