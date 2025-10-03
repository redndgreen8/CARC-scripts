# Snakefile

import os
from collections import defaultdict

# Configuration 
#FASTQ_DIR = "/project/davidwcr_264/Projects/KGP/KRRRUSCKG/240809_A00208_0364_AHJFNYDRX5/testFASTQs/"
FASTQ_DIR = "/project/davidwcr_264/Projects/KGP/KRRRUSCKG/240809_A00208_0364_AHJFNYDRX5/FASTQs/"
OUTPUT_DIR = "/scratch1/rdagnew/workingPipe/projectsHold/KRRRUSCKG"
STAR_PATH = "/home1/rdagnew/anaconda3/envs/rnaseq_umi/bin"
STAR_REF = "/project/salhia_618/genomes/Rnor_6/STARnew_genome_oh99"
STAR_GTF = "/project/salhia_618/genomes/Rnor_6/Rnor_6.gtf"
SAMTOOLS_PATH = "/home1/rdagnew/anaconda3/envs/rnaseq_umi/bin/samtools"
FEAT_CT_PATH = "/home1/rdagnew/anaconda3/envs/rnaseq_umi/bin/featureCounts"
R_PATH = "/home1/rdagnew/anaconda3/envs/rnaseq_umi/bin/Rscript"
MQC_PATH = "/home1/rdagnew/anaconda3/envs/rnaseq_umi/bin/multiqc"
SCRIPT_PATH = "/project/davidwcr_264/Projects/KGP/scripts"

# Function to get fastq groups
def get_fastq_groups():
   groups = defaultdict(lambda: {'R1': [], 'R2': []})
   for file in os.listdir(FASTQ_DIR):
       if file.endswith('.fastq.gz'):
           base = '_'.join(file.split('_')[:7])
           read_type = 'R1' if '_R1' in file else 'R2'
           groups[base][read_type].append(os.path.join(FASTQ_DIR, file))
   return groups

# Get sample information
GROUPS = get_fastq_groups()
SAMPLES = list(GROUPS.keys())

#print(SAMPLES)
rule all:
    input:
        expand("{outdir}/{sample}/alignedBAM/{sample}.Aligned.out.sorted.bam", outdir=OUTPUT_DIR, sample=SAMPLES),
        expand("{outdir}/{sample}/alignedBAM/{sample}.Aligned.out.sorted.bam.bai", outdir=OUTPUT_DIR, sample=SAMPLES),
#        expand("{outdir}/{sample}/chimericBAM/{sample}.Chimeric.out.sorted.bam", outdir=OUTPUT_DIR, sample=SAMPLES),
#        expand("{outdir}/{sample}/chimericBAM/{sample}.Chimeric.out.sorted.bai", outdir=OUTPUT_DIR, sample=SAMPLES),
#        expand("{outdir}/{sample}/junctions/{sample}.starChimeric.junctions", outdir=OUTPUT_DIR, sample=SAMPLES),
        expand("{outdir}/{sample}/junctions/{sample}.starAligned.junctions", outdir=OUTPUT_DIR, sample=SAMPLES),
        expand("{outdir}/{sample}/stats/{sample}.Log.final.out", outdir=OUTPUT_DIR, sample=SAMPLES),
        "multiqc/multiqc_report.html"















rule count_features:
    """
    Count reads per gene using featureCounts
    More flexible than htseq-count for paired-end data
    """
    input:
        bam = "{outdir}/{sample}/alignedBAM/{sample}.Aligned.out.sorted.bam",
        bai = "{outdir}/{sample}/alignedBAM/{sample}.Aligned.out.sorted.bam.bai",
        gtf = STAR_GTF,
    output:
        outd = "{outdir}/{sample}/counts",
        counts = "{outdir}/{sample}/counts/{sample}_counts.txt",
        summary = "{outdir}/{sample}/counts/{sample}_counts.txt.summary"
    params:
        fc_path = FEAT_CT_PATH,
        # Parameters for stranded RNA-seq
#        grid_opts = "grid_medium", 
        min_quality = 30,
    threads: 8
#    log:
  #      "logs/counts/{sample}.log"  2> {log}
    shell:
        """
        mkdir -p {output.outd}

        {params.fc_path} \
        -p -C \
        -T {threads} \
        -a {input.gtf} \
        -o {output.counts} \
        {input.bam} 
        

        touch {output.summary}
        """

rule calculate_tpm:
    input:
        counts = "{outdir}/{sample}/counts/{sample}_counts.txt",
        gtf = STAR_GTF,
    output:
        outd = "{outdir}/{sample}/expression",
        tpm = "{outdir}/{sample}/expression/{sample}_tpm.txt"
    params:
        script_path = SCRIPT_PATH,
        r_path = R_PATH
    shell:
        """
        mkdir -p {output.outd}

        {params.r_path} {params.script_path}/calculate_tpm.R \
            {input.counts} \
            {input.gtf} \
            {output.tpm} 

        """

rule multiqc:
    """
    Aggregate QC reports
    """
    input:
        star = expand("{outdir}/{sample}/stats/{sample}.Log.final.out", outdir=OUTPUT_DIR, sample=SAMPLES),
        counts = expand("{outdir}/{sample}/counts/{sample}_counts.txt.summary", outdir=OUTPUT_DIR, sample=SAMPLES),
    output:
        outd = "multiqc",
        html = "multiqc/multiqc_report.html"
   # log:
    #    "logs/multiqc/multiqc.log"
    shell:
        """
		mkdir -p {output.outd}
        
        multiqc \
        --force \
        --outdir {output.outd} \
        {input} 
        """
