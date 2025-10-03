import os
import re
import argparse

def rename_file(filename):
    # Split the filename into parts
    parts = filename.split('_')
    
    # Extract the number from the first part
    number = re.search(r'\d+', parts[0]).group()
    
    # Pad the number to 4 digits
    padded_number = number.zfill(4)
    
    # Construct the new filename
   # new_name = f"H23_{padded_number}_01_WBLCON1TD_R00270S8A2M0000P0000_C1_KHWGSH_A01310_HTTFKDMXY_AAGACGGA-CAATCTGA_S1_L{parts[3][1:].zfill(3)}_{parts[4][0]}1_001.fastq.gz"
        # Determine if it's R1 or R2
    #read_number = 'R1' if parts[4] == '1' else 'R2'
    read_number = 'R' + parts[-1].split('.')[0]
 
    # Construct the new filename
    new_name = f"H23-KRAS_{padded_number}_01_WBLCON1TD_R00270S8A2M0000P0000_C1_KHWGSH_A01310_HTTFKDMXY_AAGACGGA-CAATCTGA_S1_L{parts[3][1:].zfill(3)}_{read_number}_001.fastq.gz"
    
    return new_name

def main():
    parser = argparse.ArgumentParser(description="Rename files with option for dry run.")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without actually renaming files")
    args = parser.parse_args()    

    directory = '.'  # Current directory, change if needed
    
    for filename in os.listdir(directory):
        if filename.endswith('.fq.gz'):
            new_name = rename_file(filename)
            if args.dry_run:
                print(f"Dry Renaming: {filename} -> {new_name}")
            else:
                print(f"Renaming: {filename} -> {new_name}")            
            # Uncomment the following line to actually rename the files
                os.rename(os.path.join(directory, filename), os.path.join(directory, new_name))

if __name__ == "__main__":
    main()
