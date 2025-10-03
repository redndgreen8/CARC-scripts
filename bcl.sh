#!/bin/sh

module purge
module load usc
module load intel/18.0.4 intel/19.0.4 bcl2fastq2/2.20.0.422

bcl2fastq --processing-threads 20 -R $1 --output-dir  bcl2fastq_out --sample-sheet SampleSheet.csv --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --use-bases-mask Y*,I8Y*,I8,Y*
