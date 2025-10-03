#!/usr/bin/env Rscript

# Check for command line arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Three arguments required: counts_file gtf_file output_file")
}

counts_file <- args[1]
gtf_file <- args[2]
output_file <- args[3]

library(rtracklayer)
library(data.table)

# Function to get gene lengths from GTF
get_gene_lengths <- function(gtf_file) {
    gtf <- import(gtf_file)
    gtf <- as.data.table(gtf)
    exon_lengths <- gtf[type == "exon", 
                       .(length = sum(end - start + 1)), 
                       by = gene_id]
    return(exon_lengths)
}

# Function to calculate TPM
calculate_tpm <- function(counts, lengths) {
    rpk <- counts / (lengths / 1000)
    scaling_factor <- sum(rpk) / 1e6
    tpm <- rpk / scaling_factor
    return(tpm)
}

# Read counts (assuming featureCounts format)
counts_data <- fread(counts_file, skip = 2)
gene_counts <- counts_data[, .(gene_id = V1, counts = V7)]

# Get gene lengths
gene_lengths <- get_gene_lengths(gtf_file)

# Merge counts with lengths
merged_data <- merge(gene_counts, gene_lengths, by = "gene_id")

# Calculate TPM
tpm_values <- calculate_tpm(merged_data$counts, merged_data$length)
result <- data.table(gene_id = merged_data$gene_id, TPM = tpm_values)

# Write output
fwrite(result, output_file, sep = "\t")