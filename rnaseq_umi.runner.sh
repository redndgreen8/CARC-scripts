#!/usr/bin/bash

#conda init
#source ~/.bashrc

# Either use the conda shell hook
#CONDA_BASE=$(conda info --base)
#echo $CONDA_BASE
#source $CONDA_BASE/etc/profile.d/conda.sh

conda run -n rnaseq_umi snakemake -k -p -s $kgp/scripts/count_NEB_RNA1.smk -j 20 \
    --executor cluster     --cluster-config $kgp/grid.json     --cluster "{cluster.grid_opts}" \
    --rerun-incomplete --config input_dir=$1 out_dir=$2 \
    genome=/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa \
    gtf=/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/gtf/gencode.v29.primary_assembly.annotation.gtf \
    star_index=/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/star-index2
