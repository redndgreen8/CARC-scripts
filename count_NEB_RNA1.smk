"""
Complete RNA-seq Pipeline for NEBNext Multiplex Oligos (Unique Dual Index UMI Adaptors RNA Set 1)
============================================================================================

Read Structure after demultiplexing:
- R1: Forward read (98bp)
- R2: UMI sequence (12bp)
- R3: Reverse read (101bp)

Note: i7 and i5 indices (8bp each) are used during demultiplexing 
but are not part of the output FASTQ files.

Pipeline Steps:
1. BCL to FASTQ conversion (with demultiplexing)
2. UMI extraction from R2
3. Quality control (FastQC)
4. Trimming (trim_galore)
5. Lane merging
6. Alignment (STAR)
7. BAM processing (samtools)
8. UMI deduplication
9. Gene counting
10. TPM calculation
"""

import os
from os import path
import glob





# Get input directory from config
INPUT_DIR = config["input_dir"]

# Function to find samples and their read files
FASTQS = glob.glob(os.path.join(INPUT_DIR, "**/*.fastq.gz"), recursive=True)

#print(f"Found fastq: {FASTQS}")

# Get sample names - just first 11 characters 
SAMPLES = list(set(['_'.join(os.path.basename(x).split('_')[:11]) for x in FASTQS]))

#print(f"Found samples: {SAMPLES}")

# Function to get actual file paths
def get_fastq_path(wildcards, read):
    pattern = os.path.join(INPUT_DIR, "**", "*", f"{wildcards.sample}_{wildcards.lane}_R{read}_001.fastq.gz")
    files = glob.glob(pattern, recursive=True)
    if not files:
        raise ValueError(f"No files found for {pattern}")
    return files[0]


LANES = ["L001", "L002"]  # Multiple lanes observed in FastQC
OUTDIR = config["out_dir"]

# Reference files
GENOME = config["genome"]
GTF = config["gtf"]


localrules : extract_umi, trim_galore, rename_and_merge, mark_duplicates_umi,calculate_tpm, multiqc, index_dedup

rule all:
    """
    Define pipeline outputs
    """
    input:
        # QC reports
        expand(OUTDIR + "/multiqc/multiqc_report.html"),
        # Final counts
        expand(OUTDIR + "/counts/{sample}_counts.txt", sample=SAMPLES),
        # Expression values
        expand(OUTDIR + "/expression/{sample}_tpm.txt", sample=SAMPLES)



# In the rule:
rule extract_umi:
    input:
        data = lambda wildcards: get_fastq_path(wildcards, wildcards.read),
        umi = lambda wildcards: get_fastq_path(wildcards, "2")  # UMI read
    output:
        temp(OUTDIR + "/umi_extracted/{sample}_{lane}_R{read}.fastq.gz")
    params:
        umi_pattern = lambda wildcards: "'(?P<umi_1>.{12})'"  # Added quotes around pattern
    log:
        "logs/umi_extract/{sample}_{lane}_R{read}.log"
    wildcard_constraints:
        read="[13]"  # Only allow R1 or R3
    shell:
        """
        umi_tools extract \
            -I {input.umi} \
            -S /dev/null \
            -p {params.umi_pattern} \
            --read2-in={input.data} \
            --read2-out={output} \
            --extract-method=regex \
            --log={log}
        """

rule fastqc:
    """
    Quality control of reads
    Run separate FastQC for each lane
    """
    input:
        r1 = OUTDIR + "/umi_extracted/{sample}_{lane}_R1.fastq.gz",
        r3 = OUTDIR + "/umi_extracted/{sample}_{lane}_R3.fastq.gz"
    output:
        html = expand(OUTDIR + "/fastqc/{{sample}}_{{lane}}_R{read}_fastqc.html", read=["1", "3"]),
        zip = expand(OUTDIR + "/fastqc/{{sample}}_{{lane}}_R{read}_fastqc.zip", read=["1", "3"])
    params:
        grid_opts = "grid_medium", 
        outdir = OUTDIR + "/fastqc"
    log:
        "logs/fastqc/{sample}_{lane}.log"
    shell:
        "fastqc -t 8 -o {params.outdir} {input.r1} {input.r3} 2> {log}"


rule trim_galore:
    """
    Trim reads using trim_galore
    Using complete NEBNext adapter sequences
    
    Read 1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    Read 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    """
    input:
        r1 = OUTDIR + "/umi_extracted/{sample}_{lane}_R1.fastq.gz",
        r3 = OUTDIR + "/umi_extracted/{sample}_{lane}_R3.fastq.gz"
    output:
        r1 = OUTDIR + "/trimmed/{sample}_{lane}_R1_val_1.fq.gz",
        r3 = OUTDIR + "/trimmed/{sample}_{lane}_R3_val_2.fq.gz"
    params:
        # Tool paths
        trim_galore = "/project/salhia_618/bin/trim_galore",
        cutadapt = "/project/salhia_618/bin/cutadapt",
        qual_flag = 30,
        outdir = directory(OUTDIR + "/trimmed"),
        # Full NEBNext adapter sequences
        adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    log:
        "logs/trim_galore/{sample}_{lane}.log"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run trim_galore with output directory specified
        {params.trim_galore} \
            --path_to_cutadapt {params.cutadapt} \
            --adapter {params.adapter1} \
            --adapter2 {params.adapter2} \
            -q {params.qual_flag} \
            --paired \
            {input.r1} {input.r3} \
            --output_dir {params.outdir} \
            2> {log}
        """

rule rename_and_merge:
    """
    Rename trimmed files and merge lanes in one step
    """
    input:
        r1 = expand(OUTDIR + "/trimmed/{{sample}}_{lane}_R1_val_1.fq.gz", lane=LANES),
        r3 = expand(OUTDIR + "/trimmed/{{sample}}_{lane}_R3_val_2.fq.gz", lane=LANES)
    output:
        r1 = OUTDIR + "/merged/{sample}_R1.fastq.gz",
        r3 = OUTDIR + "/merged/{sample}_R3.fastq.gz"
    params:
        outdir = directory(OUTDIR + "/merged")
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Directly concatenate the files to final output
        cat {input.r1} > {output.r1}
        cat {input.r3} > {output.r3}
        """

# STAR rule with multi-mapping
rule star_align:
    """
    Align merged, trimmed reads with STAR
    Allowing multi-mapped reads (max 20)
    """
    input:
        r1 = OUTDIR + "/merged/{sample}_R1.fastq.gz",
        r3 = OUTDIR + "/merged/{sample}_R3.fastq.gz",
        index = config["star_index"]  #read length 101
    output:
        bam = temp(OUTDIR + "/star/{sample}_Aligned.out.bam"),
        log = OUTDIR + "/star/{sample}_Log.final.out"
    params:
        grid_opts = "grid_large", 
        prefix = OUTDIR + "/star/{sample}_",
        # Parameters for multi-mapped RNA-seq
        outSAMattributes = "NH HI AS nM",    # Required tags for UMI-tools
        outFilterMultimapNmax = 20,          # Allow multi-mapping
        outFilterMismatchNoverReadLmax = 0.04,
        alignIntronMin = 20,
        alignIntronMax = 1000000
    threads: 16
    log:
        "logs/star/{sample}.log"
    shell:
        """
        STAR \
        --genomeDir {input.index} \
        --readFilesIn {input.r1} {input.r3} \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outSAMattributes {params.outSAMattributes} \
        --outFilterMultimapNmax {params.outFilterMultimapNmax} \
        --outFilterMismatchNoverReadLmax {params.outFilterMismatchNoverReadLmax} \
        --alignIntronMin {params.alignIntronMin} \
        --alignIntronMax {params.alignIntronMax} \
        --outFileNamePrefix {params.prefix} \
        --runThreadN {threads} \
        2> {log}
        """

rule sort_bam:
    """
    Sort BAM by coordinates and index
    Required for downstream UMI deduplication
    """
    input:
        OUTDIR + "/star/{sample}_Aligned.out.bam"
    output:
        bam = OUTDIR + "/sorted/{sample}.sorted.bam",
        bai = OUTDIR + "/sorted/{sample}.sorted.bam.bai"
    threads: 4
    params:
        grid_opts = "grid_largemem",
    log:
        "logs/samtools/{sample}.log"
    shell:
        """
        samtools sort -@ {threads} -m 3G -o {output.bam} {input} && \
        samtools index {output.bam} 2> {log}
        """

rule mark_duplicates_umi:
    """
    UMI-aware deduplication
    Uses the 12bp UMI from NEBNext kit
    Compatible with multi-mapped reads from STAR
    """
    input:
        OUTDIR + "/sorted/{sample}.sorted.bam"
    output:
        bam = OUTDIR + "/dedup/{sample}.dedup.bam",
        metrics = OUTDIR + "/dedup/{sample}.dedup.metrics"
    params:
        # Core parameters for UMI dedup
        method = "directional",                 # Best for RNA-seq
        edit_distance_threshold = 1,            # Allow 1 bp difference in UMIs
        multi_mapping_mode = "detect-all"       # Handle multi-mapped reads properly
    log:
        "logs/dedup/{sample}.log"
    shell:
        """
        umi_tools dedup \
        --stdin={input} \
        --stdout={output.bam} \
        --output-stats={output.metrics} \
        --method={params.method} \
        --edit-distance-threshold={params.edit_distance_threshold} \
        --paired \
        2> {log}
        """

rule index_dedup:
    """
    Index deduplicated BAM
    Required for efficient counting
    """
    input:
        OUTDIR + "/dedup/{sample}.dedup.bam"
    output:
        OUTDIR + "/dedup/{sample}.dedup.bam.bai"
    log:
        "logs/index/{sample}.log"
    shell:
        "samtools index {input} 2> {log}"

rule count_features:
    """
    Count reads per gene using featureCounts
    More flexible than htseq-count for paired-end data
    """
    input:
        bam = OUTDIR + "/dedup/{sample}.dedup.bam",
        bai = OUTDIR + "/dedup/{sample}.dedup.bam.bai",
        gtf = config["gtf"]
    output:
        counts = OUTDIR + "/counts/{sample}_counts.txt",
        summary = OUTDIR + "/counts/{sample}_counts.txt.summary"
    params:
        # Parameters for stranded RNA-seq
        grid_opts = "grid_medium", 
        stranded = config.get("stranded", "2"),  # 0=unstranded, 1=stranded, 2=reversely stranded
        min_quality = 30,
        feature_type = "exon",
        attribute_type = "gene_id"
    threads: 8
    log:
        "logs/counts/{sample}.log"
    shell:
        """
        featureCounts \
        -p -B -C \
        -T {threads} \
        -s {params.stranded} \
        -Q {params.min_quality} \
        -t {params.feature_type} \
        -g {params.attribute_type} \
        -a {input.gtf} \
        -o {output.counts} \
        {input.bam} \
        2> {log}
        """

rule calculate_tpm:
    input:
        counts = OUTDIR + "/counts/{sample}_counts.txt",
        gtf = config["gtf"]
    output:
        tpm = OUTDIR + "/expression/{sample}_tpm.txt"
    log:
        "logs/tpm/{sample}.log"
    shell:
        """
        Rscript scripts/calculate_tpm.R \
            {input.counts} \
            {input.gtf} \
            {output.tpm} \
            2> {log}
        """

rule multiqc:
    """
    Aggregate QC reports
    """
    input:
        fastqc = expand(OUTDIR + "/fastqc/{sample}_{lane}_R{read}_fastqc.zip", 
                       sample=SAMPLES, lane=LANES, read=["1", "3"]),
        star = expand(OUTDIR + "/star/{sample}_Log.final.out", 
                     sample=SAMPLES),
        dedup = expand(OUTDIR + "/dedup/{sample}.dedup.metrics", 
                      sample=SAMPLES),
        counts = expand(OUTDIR + "/counts/{sample}_counts.txt.summary", 
                       sample=SAMPLES),
        trims = expand(OUTDIR + "/trimmed/{sample}_{lane}_R{read}.fastq.gz_trimming_report.txt", 
            sample=SAMPLES, 
            lane=LANES,
            read=['1','3'])
    output:
        html = OUTDIR + "/multiqc/multiqc_report.html"
    params:
        outdir = OUTDIR + "/multiqc"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        """
        multiqc \
        --force \
        --outdir {params.outdir} \
        {OUTDIR} \
        2> {log}
        """
