#!/usr/bin/env python3
"""
FastQC Parser and Visualization Tool

This script parses FastQC output files, extracts sequence information,
calculates sizes in GB, plots histograms, and identifies underperforming
samples based on user-defined thresholds.
"""

import os
import re
import glob
import zipfile
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


class FastQCParser:
    """Parse FastQC output files and extract relevant metrics."""
    
    def __init__(self, input_dir, output_dir):
        """
        Initialize the FastQC parser.
        
        Args:
            input_dir (str): Directory containing FastQC output files (.zip or unzipped folders)
            output_dir (str): Directory to save results and plots
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.samples = {}
        self.summary_df = None
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
    
    def find_fastqc_files(self):
        """Find all FastQC output files in the input directory."""
        # Look for .zip files
        zip_files = glob.glob(os.path.join(self.input_dir, "*_fastqc.zip"))
        
        # Look for unzipped directories
        directories = glob.glob(os.path.join(self.input_dir, "*_fastqc"))
        
        return zip_files, directories
    
    def parse_fastqc_data(self, filepath, is_zip=True):
        """
        Parse a single FastQC output file.
        
        Args:
            filepath (str): Path to FastQC output file or directory
            is_zip (bool): Whether the file is a zip file or directory
            
        Returns:
            dict: Dictionary containing extracted metrics
        """
        sample_name = os.path.basename(filepath).replace('_fastqc.zip', '').replace('_fastqc', '')
        metrics = {'sample_name': sample_name}
        
        # Get the path to the fastqc_data.txt file
        if is_zip:
            with zipfile.ZipFile(filepath, 'r') as zip_ref:
                data_file = [f for f in zip_ref.namelist() if f.endswith('fastqc_data.txt')][0]
                fastqc_data = zip_ref.read(data_file).decode('utf-8').splitlines()
        else:
            data_file = os.path.join(filepath, 'fastqc_data.txt')
            with open(data_file, 'r') as f:
                fastqc_data = f.readlines()
        
        # Extract basic statistics
        basic_stats_section = False
        for line in fastqc_data:
            line = line.strip()
            
            # Start of Basic Statistics section
            if line == ">>Basic Statistics	pass":
                basic_stats_section = True
                continue
            # End of Basic Statistics section
            elif line.startswith(">>END_MODULE") and basic_stats_section:
                basic_stats_section = False
                continue
            
            # Parse basic statistics
            if basic_stats_section and line and not line.startswith(">>"):
                key, value = line.split('\t')
                if key == "Filename":
                    metrics['filename'] = value
                elif key == "Total Sequences":
                    metrics['total_sequences'] = int(value)
                elif key == "Sequence length":
                    # Handle range format (e.g., "15-290")
                    if '-' in value:
                        min_len, max_len = value.split('-')
                        metrics['min_sequence_length'] = int(min_len)
                        metrics['max_sequence_length'] = int(max_len)
                        metrics['avg_sequence_length'] = (int(min_len) + int(max_len)) / 2
                    else:
                        metrics['min_sequence_length'] = int(value)
                        metrics['max_sequence_length'] = int(value)
                        metrics['avg_sequence_length'] = int(value)
                elif key == "%GC":
                    metrics['gc_content'] = float(value)
        
        # Calculate size in GB
        # Assuming 4 bytes per base (for FASTQ with quality scores)
        if 'total_sequences' in metrics and 'avg_sequence_length' in metrics:
            size_bytes = metrics['total_sequences'] * metrics['avg_sequence_length'] * 4
            metrics['size_gb'] = size_bytes / (1024**3)
        
        # Extract module status
        metrics['modules_status'] = {}
        for line in fastqc_data:
            line = line.strip()
            if line.startswith(">>") and not line.startswith(">>END_MODULE"):
                module_info = line.split('\t')
                if len(module_info) >= 2:
                    module_name = module_info[0][2:]  # Remove '>>' prefix
                    module_status = module_info[1]
                    metrics['modules_status'][module_name] = module_status
        
        return metrics
    
    def parse_all_files(self):
        """Parse all FastQC files in the input directory."""
        zip_files, directories = self.find_fastqc_files()
        
        # Parse zip files
        for zip_file in zip_files:
            sample_metrics = self.parse_fastqc_data(zip_file, is_zip=True)
            self.samples[sample_metrics['sample_name']] = sample_metrics
        
        # Parse directories
        for directory in directories:
            sample_metrics = self.parse_fastqc_data(directory, is_zip=False)
            self.samples[sample_metrics['sample_name']] = sample_metrics
        
        # Create summary DataFrame
        self.create_summary_dataframe()
        
        return self.samples
    
    def create_summary_dataframe(self):
        """Create a summary DataFrame from parsed samples."""
        if not self.samples:
            print("No samples were parsed. Run parse_all_files() first.")
            return None
        
        # Convert to DataFrame
        self.summary_df = pd.DataFrame.from_dict(self.samples, orient='index')
        
        # Extract module status columns
        if 'modules_status' in self.summary_df.columns:
            module_statuses = self.summary_df['modules_status'].apply(pd.Series)
            self.summary_df = pd.concat([self.summary_df.drop(['modules_status'], axis=1), module_statuses], axis=1)
        
        # Reset index to make sample_name a column
        self.summary_df.reset_index(drop=True, inplace=True)
        
        return self.summary_df
    
    def save_summary(self, filename='fastqc_summary.csv'):
        """Save the summary DataFrame to a CSV file."""
        if self.summary_df is None:
            print("No summary data available. Run parse_all_files() first.")
            return
        
        output_file = os.path.join(self.output_dir, filename)
        self.summary_df.to_csv(output_file, index=False)
        print(f"Summary saved to {output_file}")
    
    def plot_size_histogram(self, bin_width=0.5, figsize=(10, 6)):
        """
        Plot a histogram of sample sizes in GB.
        
        Args:
            bin_width (float): Width of histogram bins in GB
            figsize (tuple): Figure size (width, height) in inches
        """
        if self.summary_df is None:
            print("No summary data available. Run parse_all_files() first.")
            return
        
        plt.figure(figsize=figsize)
        sns.histplot(self.summary_df['size_gb'], bins=np.arange(
            self.summary_df['size_gb'].min(),
            self.summary_df['size_gb'].max() + bin_width,
            bin_width
        ))
        plt.xlabel('Size (GB)')
        plt.ylabel('Number of Samples')
        plt.title('Distribution of Sample Sizes')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Save the plot
        output_file = os.path.join(self.output_dir, 'size_histogram.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Size histogram saved to {output_file}")
    
    def plot_sequence_count_histogram(self, bin_width=1e6, figsize=(10, 6)):
        """
        Plot a histogram of sequence counts.
        
        Args:
            bin_width (float): Width of histogram bins (in number of sequences)
            figsize (tuple): Figure size (width, height) in inches
        """
        if self.summary_df is None:
            print("No summary data available. Run parse_all_files() first.")
            return
        
        plt.figure(figsize=figsize)
        sns.histplot(self.summary_df['total_sequences'], bins=np.arange(
            self.summary_df['total_sequences'].min(),
            self.summary_df['total_sequences'].max() + bin_width,
            bin_width
        ))
        plt.xlabel('Total Sequences (millions)')
        plt.xticklabels([f'{x/1e6:.1f}' for x in plt.xticks()[0]])
        plt.ylabel('Number of Samples')
        plt.title('Distribution of Sequence Counts')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Save the plot
        output_file = os.path.join(self.output_dir, 'sequence_count_histogram.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Sequence count histogram saved to {output_file}")
        
    def plot_underperforming_samples(self, underperforming_df, figsize=(12, 10)):
        """
        Create visualizations highlighting underperforming samples.
        
        Args:
            underperforming_df (DataFrame): DataFrame of underperforming samples
            figsize (tuple): Figure size (width, height) in inches
        """
        if underperforming_df.empty:
            print("No underperforming samples to plot.")
            return
            
        # Create a figure with multiple subplots
        fig, axs = plt.subplots(2, 1, figsize=figsize)
        
        # Plot 1: Size comparison (all samples vs underperforming)
        ax = axs[0]
        sns.histplot(data=self.summary_df, x='size_gb', bins=20, 
                     color='blue', alpha=0.5, label='All Samples', ax=ax)
        if len(underperforming_df) > 0:
            sns.histplot(data=underperforming_df, x='size_gb', bins=20, 
                         color='red', alpha=0.5, label='Underperforming', ax=ax)
        ax.set_xlabel('Size (GB)')
        ax.set_ylabel('Number of Samples')
        ax.set_title('Size Distribution: All vs. Underperforming Samples')
        ax.legend()
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Plot 2: Sequence count comparison
        ax = axs[1]
        sns.histplot(data=self.summary_df, x='total_sequences', bins=20, 
                     color='blue', alpha=0.5, label='All Samples', ax=ax)
        if len(underperforming_df) > 0:
            sns.histplot(data=underperforming_df, x='total_sequences', bins=20, 
                         color='red', alpha=0.5, label='Underperforming', ax=ax)
        ax.set_xlabel('Total Sequences')
        ax.set_ylabel('Number of Samples')
        ax.set_title('Sequence Count Distribution: All vs. Underperforming Samples')
        ax.legend()
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Adjust layout and save
        plt.tight_layout()
        output_file = os.path.join(self.output_dir, 'underperforming_analysis.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Underperforming samples analysis saved to {output_file}")
    
    def identify_underperforming_samples(self, size_threshold=None, sequence_threshold=None, module_failures=None):
        """
        Identify underperforming samples based on thresholds.
        
        Args:
            size_threshold (float): Minimum size threshold in GB (samples below this are underperforming)
            sequence_threshold (int): Minimum sequence count threshold (samples below this are underperforming)
            module_failures (list): List of modules that must pass (samples failing these are underperforming)
            
        Returns:
            DataFrame: Filtered DataFrame of underperforming samples
        """
        if self.summary_df is None:
            print("No summary data available. Run parse_all_files() first.")
            return None
        
        # Start with all samples
        filtered_df = self.summary_df.copy()
        underperforming = pd.DataFrame(columns=filtered_df.columns)
        reasons = []
        
        # Track which samples fail which criteria
        for idx, row in filtered_df.iterrows():
            sample_underperforming = False
            reason = []
            
            # Check size threshold
            if size_threshold is not None and row['size_gb'] < size_threshold:
                sample_underperforming = True
                reason.append(f"Size below threshold: {row['size_gb']:.2f}GB < {size_threshold}GB")
            
            # Check sequence threshold
            if sequence_threshold is not None and row['total_sequences'] < sequence_threshold:
                sample_underperforming = True
                reason.append(f"Sequence count below threshold: {row['total_sequences']} < {sequence_threshold}")
            
            # Check module failures
            if module_failures is not None:
                for module in module_failures:
                    if module in row and row[module] != 'pass':
                        sample_underperforming = True
                        reason.append(f"Module '{module}' status: {row[module]}")
            
            # Add underperforming sample to results
            if sample_underperforming:
                underperforming = pd.concat([underperforming, pd.DataFrame([row])], ignore_index=True)
                reasons.append("; ".join(reason))
        
        # Add failure reasons
        if not underperforming.empty:
            underperforming['failure_reason'] = reasons
        
        # Save the filtered results
        output_file = os.path.join(self.output_dir, 'underperforming_samples.csv')
        underperforming.to_csv(output_file, index=False)
        print(f"Identified {len(underperforming)} underperforming samples. Results saved to {output_file}")
        
        # Generate a report of underperforming samples
        if not underperforming.empty:
            report_file = os.path.join(self.output_dir, 'underperforming_report.txt')
            with open(report_file, 'w') as f:
                f.write("# Underperforming Samples Report\n\n")
                f.write(f"Total samples analyzed: {len(self.summary_df)}\n")
                f.write(f"Underperforming samples: {len(underperforming)} ({len(underperforming)/len(self.summary_df)*100:.1f}%)\n\n")
                
                if size_threshold is not None:
                    f.write(f"Size threshold: {size_threshold} GB\n")
                if sequence_threshold is not None:
                    f.write(f"Sequence count threshold: {sequence_threshold}\n")
                if module_failures is not None:
                    f.write(f"Module failure checks: {', '.join(module_failures)}\n")
                
                f.write("\n## Sample Breakdown\n\n")
                for i, row in underperforming.iterrows():
                    f.write(f"- {row['sample_name']}: {row['failure_reason']}\n")
            
            print(f"Detailed report saved to {report_file}")
            
            # Plot underperforming samples
            self.plot_underperforming_samples(underperforming)
        
        return underperforming


def main():
    """Main function to run the FastQC parser."""
    parser = argparse.ArgumentParser(description='Parse FastQC output files and visualize results')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing FastQC output files')
    parser.add_argument('-o', '--output_dir', default='fastqc_results', help='Directory to save results')
    parser.add_argument('--size_threshold', type=float, default=7.0, help='Minimum size threshold in GB (samples below this are underperforming)')
    parser.add_argument('--sequence_threshold', type=int, help='Minimum sequence count threshold (samples below this are underperforming)')
    parser.add_argument('--module_failures', nargs='+', help='List of modules that must pass (samples failing these are underperforming)')
    parser.add_argument('--bin_width_gb', type=float, default=0.5, help='Bin width for size histograms in GB')
    parser.add_argument('--bin_width_seq', type=float, default=1e6, help='Bin width for sequence count histograms')
    
    args = parser.parse_args()
    
    # Initialize and run parser
    fastqc_parser = FastQCParser(args.input_dir, args.output_dir)
    fastqc_parser.parse_all_files()
    
    # Save summary
    fastqc_parser.save_summary()
    
    # Generate plots
    fastqc_parser.plot_size_histogram(bin_width=args.bin_width_gb)
    fastqc_parser.plot_sequence_count_histogram(bin_width=args.bin_width_seq)
    
    # Always run the underperforming samples analysis
    # Default thresholds can be used for filtering
    underperforming = fastqc_parser.identify_underperforming_samples(
        size_threshold=args.size_threshold,
        sequence_threshold=args.sequence_threshold,
        module_failures=args.module_failures
    )
    
    # Print summary of results
    print("\nFastQC Analysis Summary:")
    print(f"Total samples analyzed: {len(fastqc_parser.summary_df)}")
    if underperforming is not None:
        print(f"Underperforming samples: {len(underperforming)}")
        if len(underperforming) > 0:
            print("Check 'underperforming_samples.csv' and 'underperforming_report.txt' for details.")


if __name__ == "__main__":
    main()