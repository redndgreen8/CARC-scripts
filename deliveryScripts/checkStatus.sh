#!/bin/bash

# Check if two arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <root_directory> <assay_code>"
    exit 1
fi

root_directory=$1
assay_code=$2

# Define function to count directories and passed files for tools with similar structure
count_passed_files() {
    local dir=$1
    local pass_pattern=$2
    local num_dirs=$(ls -d "${dir}"/* | wc -l)
    local num_pass=$(ls "${dir}"/*/*${pass_pattern} | wc -l)
    echo "$num_dirs $num_pass $pass_pattern"
}

# Define function to count directories and passed files for Manta and Strelka
count_manta_strelka() {
    local dir=$1
    local pass_pattern=$2
    local num_dirs=$(ls -d "${dir}"/*/ | wc -l)
    local num_pass=$(ls -d "${dir}"/*${pass_pattern} | wc -l)
    echo "$num_dirs $num_pass $pass_pattern"
}

# Function to check status for a specific tool
check_tool_status() {
    local tool=$1
    local pass_pattern=$2
    local tool_dir="${root_directory}/${tool}"
    echo -n "$tool: "
    case $tool in
        "manta"|"strelka")
            count_manta_strelka "$tool_dir" "$pass_pattern"
            ;;
        *)
            count_passed_files "$tool_dir" "$pass_pattern"
            ;;
    esac
}

# Check status for each tool
check_tool_status "${assay_code}" "cpFqPass"
check_tool_status "${assay_code}" "fastqcPass"
check_tool_status "${assay_code}" "dnaAlignPass"
check_tool_status "${assay_code}" "recalibratePass"
check_tool_status "${assay_code}" "mergeBamPass"
check_tool_status "${assay_code}" "mdPass"



check_tool_status "cna" "cnaExomePass"
check_tool_status "trn" "trnPass"

check_tool_status "freebayes" "snpEffPass"
check_tool_status "hc" "snpEffPass"
check_tool_status "mutect" "snpEffPass"
check_tool_status "seurat" "snpEffPass"

check_tool_status "manta" "mantaPass"
check_tool_status "strelka" "strelkaPass"

