#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <primary_reads_dir> <topoff_reads_dir> <output_directory>"
    exit 1
fi

# Input directories
primary_dir="$1"
topoff_dir="$2"

# Output directory
output_dir="$3"
mkdir -p "$output_dir"

# Set up log file
log_file="$output_dir/merge_log.txt"
exec > >(tee -a "$log_file") 2>&1

# Function to extract fields from filename
get_fields() {
    echo "$1" | cut -d'_' -f1-10
}

# Process files
for primary_file in "$primary_dir"/*.fastq.gz; do
    # Extract the base name without the path
    base_name=$(basename "$primary_file")
    
    # Extract the first 10 fields
    name_part=$(get_fields "$base_name")
    
    # Extract R1 or R2 part
    r_part=$(echo "$base_name" | cut -d'_' -f13)

    # Extract terminal part
    terminal_part=$(echo "$base_name" | cut -d'_' -f14)
    
    # Extract S number
    s_part=$(echo "$base_name" | cut -d'_' -f11)
    
    # Find all matching files in primary directory
    primary_files=$(find "$primary_dir" -name "${name_part}*_L00*_${r_part}*.fastq.gz")
    
    # Find all matching files in topoff directory
    topoff_files=$(find "$topoff_dir" -name "${name_part}*_L00*_${r_part}*.fastq.gz")
    
    # Combine all matching files
    all_files="$primary_files $topoff_files"

    file_count=$(echo $all_files | wc -w)
    
    if [ "$file_count" -gt 1 ]; then
        echo "Concatenating files for $base_name"
        new_name=$(echo "$base_name" | sed 's/\(L00\)[0-9]/\11/')
        echo "Concatenating $all_files to $new_name"
        cat $all_files > "$output_dir/$new_name"
    else
        echo "No matching files found for $base_name. linking primary file only."
        ln -s "$primary_file" "$output_dir/$base_name"
    fi
done

echo "Concatenation complete. Check the $output_dir directory for results."
