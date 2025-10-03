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

rule merge_lanes:
    input:
        r1 = lambda wildcards: GROUPS[wildcards.sample]['R1'],
        r2 = lambda wildcards: GROUPS[wildcards.sample]['R2']
    output:
        r1 = "{outdir}/{sample}/temp/merged_R1.fastq.gz",
        r2 = "{outdir}/{sample}/temp/merged_R2.fastq.gz",
        #outlog = "{outdir}/logs", {output.outlog}
    params:
        star_path = STAR_PATH,
        star_ref = STAR_REF,
        temp_dir = lambda wildcards: f"{wildcards.outdir}/{wildcards.sample}/temp"
    resources:
        mem_mb = 8000,
        time = "01:00:00",
    shell:
        """
        mkdir -p {params.temp_dir} 

        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """

rule star_align:
    input:
        r1 = "{outdir}/{sample}/temp/merged_R1.fastq.gz",
        r2 = "{outdir}/{sample}/temp/merged_R2.fastq.gz"
    output:
        aligned_bam = "{outdir}/{sample}/temp/Aligned.out.bam",
#        chimeric_bam = temp("{outdir}/{sample}/temp/Chimeric.out.bam"),
#        chimeric_junction = temp("{outdir}/{sample}/temp/Chimeric.out.junction"),
        sj_tab = "{outdir}/{sample}/temp/SJ.out.tab",
        log_final = "{outdir}/{sample}/temp/Log.final.out",
      #  logs = "{outdir}/logs/star"
    params:
        star_path = STAR_PATH,
        star_ref = STAR_REF,
        temp_dir = lambda wildcards: f"{wildcards.outdir}/{wildcards.sample}/temp",
        rg_id = lambda wildcards: f"ID:{wildcards.sample}",
        lib_id = lambda wildcards: f"LB:{wildcards.sample}",
        sample_id = lambda wildcards: f"SM:{wildcards.sample}",
        pu_id = lambda wildcards: f"PU:{wildcards.sample}",
        cn_id = "USC-KGP"
    resources:
        mem_mb = 60000,
        time = "12:00:00",
    threads: 16
  #  log:
   #     "{outdir}/logs/star/{sample}.log"   {output.logs}
    shell:
        """
        mkdir -p {params.temp_dir} 
        cd {params.temp_dir}
        {params.star_path}/STAR --genomeDir {params.star_ref} \
            --runMode alignReads \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outBAMcompression -1 \
            --outSAMattrRGline {params.rg_id} {params.lib_id} PL:Illumina CN:{params.cn_id} {params.sample_id} {params.pu_id} \
            --runThreadN {threads} 
        """




rule post_process:
    input:
       aligned_bam = "{outdir}/{sample}/temp/Aligned.out.bam",
#       chimeric_bam = "{outdir}/{sample}/temp/Chimeric.out.bam",
#       chimeric_junction = "{outdir}/{sample}/temp/Chimeric.out.junction",
       sj_tab = "{outdir}/{sample}/temp/SJ.out.tab",
       log_final = "{outdir}/{sample}/temp/Log.final.out"
    output:
       aligned_sorted_bam = "{outdir}/{sample}/alignedBAM/{sample}.Aligned.out.sorted.bam",
       aligned_sorted_bai = "{outdir}/{sample}/alignedBAM/{sample}.Aligned.out.sorted.bam.bai",
#       chimeric_sorted_bam = "{outdir}/{sample}/chimericBAM/{sample}.Chimeric.out.sorted.bam",
#       chimeric_sorted_bai = "{outdir}/{sample}/chimericBAM/{sample}.Chimeric.out.sorted.bai",
#       chimeric_junctions = "{outdir}/{sample}/junctions/{sample}.starChimeric.junctions",
       aligned_junctions = "{outdir}/{sample}/junctions/{sample}.starAligned.junctions",
       log_final = "{outdir}/{sample}/stats/{sample}.Log.final.out"
    resources:
        mem_mb = 60000,
        time = "04:00:00",
    params:
        samtools_path = SAMTOOLS_PATH,
        aligned_dir = lambda wildcards: f"{wildcards.outdir}/{wildcards.sample}/alignedBAM",
#        chimeric_dir = lambda wildcards: f"{wildcards.outdir}/{wildcards.sample}/chimericBAM",
        junctions_dir = lambda wildcards: f"{wildcards.outdir}/{wildcards.sample}/junctions"
    threads: 4
    shell:
       """
       mkdir -p {params.aligned_dir}
       mkdir -p {params.junctions_dir}

       {params.samtools_path} sort -@{threads} -m8G {input.aligned_bam} -o {output.aligned_sorted_bam}

       {params.samtools_path} index {output.aligned_sorted_bam}

       mv {input.sj_tab} {output.aligned_junctions}

       mkdir -p $(dirname {output.log_final})
       cp {input.log_final} {output.log_final}
       
       #rm -rf $(dirname {input.aligned_bam})
       """



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
        outd = directory("{outdir}/{sample}/counts"),
        counts = "{outdir}/{sample}/counts/{sample}_counts.txt",
        summary = "{outdir}/{sample}/counts/{sample}_counts.txt.summary"
    params:
        fc_path = FEAT_CT_PATH,
        # Parameters for stranded RNA-seq
#        grid_opts = "grid_medium", 
        min_quality = 30,
    resources:
        mem_mb = 60000,
        time = "12:00:00",
    threads: 10
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
        outd = directory("{outdir}/{sample}/expression"),
        tpm = "{outdir}/{sample}/expression/{sample}_tpm.txt"
    params:
        script_path = SCRIPT_PATH,
        r_path = R_PATH
    resources:
        mem_mb = 8000,
        time = "02:00:00",
    threads: 1    
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
        outd = directory("multiqc"),
        html = "multiqc/multiqc_report.html"
   # log:
    #    "logs/multiqc/multiqc.log"
    resources:
        mem_mb = 8000,
        time = "02:00:00",
    threads: 1
    shell:
        """
        mkdir -p {output.outd}
        export PYTHONPATH=/project/davidwcr_264/Packages/pythonModule/python-3.9.2
        export OPENBLAS_NUM_THREADS=4

        /project/davidwcr_264/Packages/pythonModule/python-3.9.2/bin/multiqc \
        --force \
        --outdir {output.outd} \
        {input} 
        """


