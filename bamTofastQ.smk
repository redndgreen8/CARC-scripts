import os
import sys
import re

# Get samples and chunks from Snakemake config
samples_file = config.get('samples_file', None)
if samples_file:
    # Read samples from file
    with open(samples_file, 'r') as f:
        SAMPLES = [line.strip() for line in f if line.strip()]
else:
    # Use default samples if no file provided
    SAMPLES = config.get('samples', ["sample1", "sample2"])

NUM_CHUNKS = config.get('chunks', 4)

# Use SAMPLES as bam_files
bam_files = SAMPLES

def get_sample_number(bam):
    m = re.search(r'Samp(\d+)', bam)
    return m.group(1) if m else None

# Automatically assign A00xxx codes starting at 152
sample_numbers = []
for bam in bam_files:
    m = re.search(r'Samp(\d+)', bam)
    if m:
        sample_num = m.group(1)
        if sample_num not in sample_numbers:
            sample_numbers.append(sample_num)

# Assign A-codes based on sample order: first sample gets A00152, second gets A00153, etc.
sample_to_a_code = {}
for i, sample in enumerate(bam_files):
    sample_number = get_sample_number(sample)
    if sample_number:
        a_code = f"A00{152 + i}"  # A00152, A00153, A00154, etc.
        sample_to_a_code[sample_number] = a_code

def get_a_code(sample_number):
    return sample_to_a_code.get(sample_number, "UNKNOWN")

def make_prefix(sample_number, a_code):
    return f"KILMN5BASE_{int(sample_number):04d}_01_LibTube102_R00250S8A2M0000P0000_T1_KHWGSH_{a_code}_HV7THDSX5_TAAGGCGA-AGGCTTAG_S1_L001"

def generate_output_files(bam_files, num_chunks):
    chunks = list(range(1, num_chunks + 1))
    reads = [1, 2]
    output_files = []
    for bam in bam_files:
        sample_number = get_sample_number(bam)
        a_code = get_a_code(sample_number)
        prefix = make_prefix(sample_number, a_code)
        for chunk in chunks:
            for read in reads:
                output_files.append(f"{prefix}_R{read}_00{chunk:02d}.fastq.gz")
    return output_files

# Get output files
output_files = generate_output_files(bam_files, NUM_CHUNKS)

rule all:
    input:
        output_files

# Convert BAM to paired FASTQ files
rule bam_to_fastq:
    input:
        bam="data/{sample}.bam"
    output:
        r1=temp("output/temp/{sample}_R1.fastq.gz"),
        r2=temp("output/temp/{sample}_R2.fastq.gz")
    conda:
        "envs/samtools.yaml"
    log:
        "logs/{sample}_bam_to_fastq.log"
    threads: 4
    shell:
        """
        samtools fastq -@ {threads} \
            -1 {output.r1} -2 {output.r2} \
            -0 /dev/null -s /dev/null -n {input.bam} 2> {log}
        """

# Split FASTQ files using seqtk
rule split_fastq:
    input:
        r1="output/temp/{sample}_R1.fastq.gz",
        r2="output/temp/{sample}_R2.fastq.gz"
    output:
        r1_chunks=expand("output/fastq/{{sample}}_R1_chunk_{chunk:02d}.fastq.gz", chunk=CHUNKS),
        r2_chunks=expand("output/fastq/{{sample}}_R2_chunk_{chunk:02d}.fastq.gz", chunk=CHUNKS)
    conda:
        "envs/seqtk.yaml"
    log:
        "logs/{sample}_split_fastq.log"
    threads: 2
    shell:
        """
        # Count total reads
        TOTAL_READS=$(zcat {input.r1} | wc -l)
        READS_PER_PAIR=$((TOTAL_READS / 4))
        READS_PER_CHUNK=$((READS_PER_PAIR / """ + str(NUM_CHUNKS) + """))
        LINES_PER_CHUNK=$((READS_PER_CHUNK * 4))
        
        echo "Total reads: $READS_PER_PAIR, Reads per chunk: $READS_PER_CHUNK" > {log}
        
        # Create temp directory for splitting
        mkdir -p output/temp_split/{wildcards.sample}
        
        # Split R1
        zcat {input.r1} | split -l $LINES_PER_CHUNK \
            --numeric-suffixes=1 --suffix-length=2 \
            - output/temp_split/{wildcards.sample}/R1_chunk_
        
        # Split R2  
        zcat {input.r2} | split -l $LINES_PER_CHUNK \
            --numeric-suffixes=1 --suffix-length=2 \
            - output/temp_split/{wildcards.sample}/R2_chunk_
        
        # Compress and rename chunks
        for i in $(seq -f "%02g" 1 """ + str(NUM_CHUNKS) + """); do
            gzip -c output/temp_split/{wildcards.sample}/R1_chunk_$i > output/fastq/{wildcards.sample}_R1_chunk_$i.fastq.gz
            gzip -c output/temp_split/{wildcards.sample}/R2_chunk_$i > output/fastq/{wildcards.sample}_R2_chunk_$i.fastq.gz
        done
        
        # Cleanup temp files
        rm -rf output/temp_split/{wildcards.sample}/
        """

# Validate that chunks have paired reads
rule validate_chunks:
    input:
        r1_chunks=expand("output/fastq/{{sample}}_R1_chunk_{chunk:02d}.fastq.gz", chunk=CHUNKS),
        r2_chunks=expand("output/fastq/{{sample}}_R2_chunk_{chunk:02d}.fastq.gz", chunk=CHUNKS)
    output:
        validation="output/validation/{sample}_validation.txt"
    log:
        "logs/{sample}_validation.log"
    shell:
        """
        echo "Validating chunks for {wildcards.sample}" > {output.validation}
        
        for i in $(seq -f "%02g" 1 """ + str(NUM_CHUNKS) + """); do
            R1_COUNT=$(zcat output/fastq/{wildcards.sample}_R1_chunk_${{i}}.fastq.gz | wc -l)
            R2_COUNT=$(zcat output/fastq/{wildcards.sample}_R2_chunk_${{i}}.fastq.gz | wc -l)
            
            if [ $R1_COUNT -eq $R2_COUNT ]; then
                echo "Chunk ${{i}}: PASS (R1: $R1_COUNT, R2: $R2_COUNT lines)" >> {output.validation}
            else
                echo "Chunk ${{i}}: FAIL (R1: $R1_COUNT, R2: $R2_COUNT lines)" >> {output.validation}
            fi
        done
        """

