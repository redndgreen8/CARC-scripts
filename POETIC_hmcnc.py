# Copy Number Variation Analysis Snakemake Workflow
# Best practices implementation

import os
from pathlib import Path



## in /project/salhia_618/POETIC_WGBS_Data/hmcnc

# Configuration
configfile: config["config_file"]

# Global variables
SAMPLES = config["samples"]
REFERENCE_BINS = config["reference_bins"]
GENE_ANNOTATIONS = config["gene_annotations"]
BIN_SIZE = config["bin_size"]
VITERBI_PARAMS = config["viterbi_params"]
TOOLS = config["tools"]

# Output directory structure
OUTDIR = config["output_dir"]
BAM_DIR = config["bam_dir"]

# Wildcards constraints
wildcard_constraints:
    sample="[A-Za-z0-9_-]+"

# Final target rule
rule all:
    input:
        expand("{outdir}/{sample}/{sample}.cn.bed", outdir=OUTDIR, sample=SAMPLES),
        expand("{outdir}/{sample}/{sample}.cn.DUP.bed", outdir=OUTDIR, sample=SAMPLES),
        expand("{outdir}/{sample}/{sample}.cn.DEL.bed", outdir=OUTDIR, sample=SAMPLES),
        expand("{outdir}/{sample}/{sample}.gene_count.bed", outdir=OUTDIR, sample=SAMPLES),
        expand("{outdir}/{sample}/{sample}.cn_gene.bed", outdir=OUTDIR, sample=SAMPLES)

# Rule 1: Convert BAM to BED using samToBed
rule bam_to_bed:
    input:
        bam="{bam_dir}/{sample}.bam"
    output:
        bed="{outdir}/{sample}/{sample}.samTo.bed"
    params:
        samtoBed=TOOLS["samToBed"],
        threads=config["threads"]
    threads: config["threads"]
    log:
        "{outdir}/{sample}/logs/{sample}.bam_to_bed.log"
    shell:
        """
        mkdir -p {wildcards.outdir}/{wildcards.sample}/logs
        samtools view -@ {params.threads} {input.bam} | \
        {params.samtoBed} /dev/stdin --useH --flag > {output.bed} 2> {log}
        """

# Rule 2: Count reads in genomic bins
rule count_reads_in_bins:
    input:
        bed="{outdir}/{sample}/{sample}.samTo.bed",
        ref_bins=REFERENCE_BINS
    output:
        counts="{outdir}/{sample}/{sample}.read_counts_{bin_size}bp.bed.gz".format(bin_size=BIN_SIZE)
    log:
        "{outdir}/{sample}/logs/{sample}.count_reads_in_bins.log"
    shell:
        """
        intersectBed -sorted -c -a {input.ref_bins} -b {input.bed} | \
        bgzip -c > {output.counts} 2> {log}
        """

# Rule 3: Calculate mean coverage
rule calculate_mean_coverage:
    input:
        counts="{outdir}/{sample}/{sample}.read_counts_{bin_size}bp.bed.gz".format(bin_size=BIN_SIZE)
    output:
        mean_cov="{outdir}/{sample}/{sample}.mean_cov.txt"
    log:
        "{outdir}/{sample}/logs/{sample}.calculate_mean_coverage.log"
    shell:
        """
        zcat {input.counts} | \
        awk 'BEGIN{{OFS="\\t";c=0;sum=0;}} {{sum=sum+$4;c=c+1;}} END{{print sum/c;}}' | \
        tail -1 > {output.mean_cov} 2> {log}
        """

# Rule 4: Run Viterbi algorithm for copy number calling
rule viterbi_analysis:
    input:
        counts="{outdir}/{sample}/{sample}.read_counts_{bin_size}bp.bed.gz".format(bin_size=BIN_SIZE),
        mean_cov="{outdir}/{sample}/{sample}.mean_cov.txt"
    output:
        viterbi_out="{outdir}/{sample}/{sample}.viterout.txt"
    params:
        viterbi_tool=TOOLS["viterbi"],
        bin_size=BIN_SIZE,
        viterbi_params=VITERBI_PARAMS,
        output_prefix="{outdir}/{sample}/{sample}"
    log:
        "{outdir}/{sample}/logs/{sample}.viterbi_analysis.log"
    shell:
        """
        mean_cov=$(cat {input.mean_cov})
        zcat {input.counts} | cut -f 4 | \
        {params.viterbi_tool} /dev/stdin $mean_cov {params.output_prefix} \
        {params.bin_size} {params.viterbi_params[0]} {params.viterbi_params[1]} > {log} 2>&1
        """

# Rule 5: Combine read counts with copy number calls
rule combine_counts_and_cn:
    input:
        counts="{outdir}/{sample}/{sample}.read_counts_{bin_size}bp.bed.gz".format(bin_size=BIN_SIZE),
        viterbi_out="{outdir}/{sample}/{sample}.viterout.txt"
    output:
        combined="{outdir}/{sample}/{sample}.cn_{bin_size}bp.bed".format(bin_size=BIN_SIZE)
    log:
        "{outdir}/{sample}/logs/{sample}.combine_counts_and_cn.log"
    shell:
        """
        paste <(zcat {input.counts}) <(cat {input.viterbi_out}) > {output.combined} 2> {log}
        """

# Rule 6: Merge adjacent bins with same copy number
rule merge_adjacent_bins:
    input:
        combined="{outdir}/{sample}/{sample}.cn_{bin_size}bp.bed".format(bin_size=BIN_SIZE)
    output:
        merged="{outdir}/{sample}/{sample}.cn.bed"
    log:
        "{outdir}/{sample}/logs/{sample}.merge_adjacent_bins.log"
    shell:
        """
        # Solution 1: bedtools cluster + groupby (cleanest)
        sort -k1,1 -k2,2n {input.combined} | \
        bedtools cluster -i stdin | \
        bedtools groupby -i stdin -g 1,5,6 -c 2,3 -o min,max | \
        awk 'BEGIN{{OFS="\\t"}} {{print $1,$4,$5,$2}}' | \
        sort -k1,1 -k2,2n > {output.merged} 2> {log}
        
        """

# Rule 7: Filter duplications (CN > 2)
rule filter_duplications:
    input:
        merged="{outdir}/{sample}/{sample}.cn.bed"
    output:
        duplications="{outdir}/{sample}/{sample}.cn.DUP.bed"
    log:
        "{outdir}/{sample}/logs/{sample}.filter_duplications.log"
    shell:
        """
        awk '$4>2' {input.merged} > {output.duplications} 2> {log}
        """

# Rule 8: Filter deletions (CN < 2)
rule filter_deletions:
    input:
        merged="{outdir}/{sample}/{sample}.cn.bed"
    output:
        deletions="{outdir}/{sample}/{sample}.cn.DEL.bed"
    log:
        "{outdir}/{sample}/logs/{sample}.filter_deletions.log"
    shell:
        """
        awk '$4<2' {input.merged} > {output.deletions} 2> {log}
        """

# Rule 9: Generate gene count annotation for CNVs
rule annotate_cnv_gene_counts:
    input:
        duplications="{outdir}/{sample}/{sample}.cn.DUP.bed",
        deletions="{outdir}/{sample}/{sample}.cn.DEL.bed",
        gene_bed=GENE_ANNOTATIONS
    output:
        gene_counts="{outdir}/{sample}/{sample}.gene_count.bed"
    log:
        "{outdir}/{sample}/logs/{sample}.annotate_cnv_gene_counts.log"
    shell:
        """
        cat {input.duplications} {input.deletions} | sort -k1,1 -k2,2n | \
        intersectBed -F 1 -wb -wa -a stdin -b {input.gene_bed} | \
        groupBy -g 1,2,3,6 -c 5,10,10 -o collapse,collapse,count > {output.gene_counts} 2> {log}
        """

# Rule 10: Generate copy number gene annotation
rule annotate_cnv_genes:
    input:
        duplications="{outdir}/{sample}/{sample}.cn.DUP.bed",
        deletions="{outdir}/{sample}/{sample}.cn.DEL.bed",
        gene_bed=GENE_ANNOTATIONS
    output:
        cn_genes="{outdir}/{sample}/{sample}.cn_gene.bed"
    log:
        "{outdir}/{sample}/logs/{sample}.annotate_cnv_genes.log"
    shell:
        """
        cat {input.duplications} {input.deletions} | sort -k1,1 -k2,2n | \
        intersectBed -F 1 -wb -wa -a stdin -b {input.gene_bed} | \
        groupBy -g 1,2,3,4 -c 8,8 -o collapse,count > {output.cn_genes} 2> {log}
        """

# Additional utility rules

# Rule to create output directories
rule create_directories:
    output:
        directory("{outdir}/{sample}/logs")
    shell:
        "mkdir -p {output}"

# Rule to validate input files
rule validate_inputs:
    input:
        bam="{bam_dir}/{sample}.bam",
        ref_bins=REFERENCE_BINS,
        gene_bed=GENE_ANNOTATIONS
    output:
        touch("{outdir}/{sample}/.validation_complete")
    shell:
        """
        if [ ! -f {input.bam} ]; then
            echo "Error: BAM file {input.bam} not found" >&2
            exit 1
        fi
        if [ ! -f {input.ref_bins} ]; then
            echo "Error: Reference bins file {input.ref_bins} not found" >&2
            exit 1
        fi
        if [ ! -f {input.gene_bed} ]; then
            echo "Error: Gene annotation file {input.gene_bed} not found" >&2
            exit 1
        fi
        touch {output}
        """

# Clean up rule
rule clean:
    shell:
        "rm -rf {OUTDIR}"

# Rule to generate summary statistics
rule generate_summary:
    input:
        expand("{outdir}/{sample}/{sample}.cn.bed", outdir=OUTDIR, sample=SAMPLES),
        expand("{outdir}/{sample}/{sample}.gene_count.bed", outdir=OUTDIR, sample=SAMPLES),
        expand("{outdir}/{sample}/{sample}.cn_gene.bed", outdir=OUTDIR, sample=SAMPLES)
    output:
        summary="{outdir}/summary_statistics.txt"
    params:
        samples=SAMPLES,
        outdir=OUTDIR
    shell:
        """
        echo "Copy Number Analysis Summary" > {output.summary}
        echo "============================" >> {output.summary}
        echo "" >> {output.summary}
        
        for sample in {params.samples}; do
            echo "Sample: $sample" >> {output.summary}
            dup_count=$(wc -l < {params.outdir}/$sample/$sample.cn.DUP.bed)
            del_count=$(wc -l < {params.outdir}/$sample/$sample.cn.DEL.bed)
            gene_count=$(wc -l < {params.outdir}/$sample/$sample.gene_count.bed)
            cn_gene_count=$(wc -l < {params.outdir}/$sample/$sample.cn_gene.bed)
            echo "  Duplications: $dup_count" >> {output.summary}
            echo "  Deletions: $del_count" >> {output.summary}
            echo "  Gene count annotations: $gene_count" >> {output.summary}
            echo "  Copy number gene annotations: $cn_gene_count" >> {output.summary}
            echo "" >> {output.summary}
        done
        """