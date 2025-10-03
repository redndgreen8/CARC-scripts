import sys
import os

def get_matching_fastqs(filename, read_num):
    # Get first 7 fields from: KRRRUSCKG_0033_01_RRRCONTR_R00154S8A2M0000P0000_C1_NEBKAP
    pattern = '_'.join(filename.split('_')[:7])
    
    base_dir = "/project/davidwcr_264/Projects/KGP/KRRRUSCKG/240809_A00208_0364_AHJFNYDRX5/FASTQs"
    # Single pattern - wildcard will match all lanes
    return f"{base_dir}/{pattern}*R{read_num}_001.fastq.gz"

def modify_line(line):
   # Split into pre-export and export parts
   pre_export = line[:line.find('--export=')]
   export_part = line[line.find('--export='):]
   
   # Split export part into variables at commas, keeping track of the script at end
   export_vars = export_part[9:].split(',')  # 9 to skip '--export='
   script_path = export_vars[-1].split()[-1]  # Get last word of last element
   export_vars[-1] = export_vars[-1].rsplit(None, 1)[0]  # Remove script path
   
   # Process each export variable
   for i, var in enumerate(export_vars):
       if var.startswith('RGTAGLIST='):
           parts = var[len('RGTAGLIST='):].split(';')
           new_parts = []
           for part in parts:
               new_parts.append(part.replace(' ', '_'))
           export_vars[i] = 'RGTAGLIST=' + ';'.join(new_parts)
       elif var.startswith('STARREF='):
           export_vars[i] = 'STARREF=/project/salhia_618/genomes/Rnor_6/STAR_genome_oh99'
       elif var.startswith('STARGTF='):
           export_vars[i] = 'STARGTF=/project/salhia_618/genomes/Rnor_6/Rnor_6.gtf'
       elif var.startswith('DIR='):
           export_vars[i] = 'DIR=/scratch1/rdagnew/workingPipe/projectsHold/KRRRUSCKG'
       elif var.startswith('FASTQL1='):
           parts = var.split('=', 1)
           filename = os.path.basename(parts[1].strip().split()[0])
           export_vars[i] = f"FASTQL1={get_matching_fastqs(filename, 1)}"
       elif var.startswith('FASTQL2='):
           parts = var.split('=', 1)
           filename = os.path.basename(parts[1].strip().split()[0])
           export_vars[i] = f"FASTQL2={get_matching_fastqs(filename, 2)}"
   
   # Reconstruct the line
   return f"{pre_export}--export={','.join(export_vars)} {script_path}"

# Read from stdin line by line
for line in sys.stdin:
   modified_line = modify_line(line.strip())
   print(modified_line)