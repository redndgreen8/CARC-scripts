#!/usr/bin/env python3

import sys
import pandas as pd

def create_samplesheet():
    # Read filenames from stdin
    filenames = sys.stdin.read().splitlines()
    data = []
    
    for filename in filenames:
        if '_R1_' in filename:
            parts = filename.split('_')
            sample_id_pre = "_".join(parts[0:12])
            sample_id = sample_id_pre.split("/")[-1]
            r1_file = filename
            r2_file = filename.replace('_R1_', '_R2_')
            
            data.append({
                'sample': sample_id,
                'fastq_1': r1_file,
                'fastq_2': r2_file,
                'strandedness': 'auto'
            })
    
    df = pd.DataFrame(data)
    df.to_csv('samplesheet.csv', index=False)

if __name__ == "__main__":
    create_samplesheet()