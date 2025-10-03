import os

# Configuration
if config.get("input_file"):
    with open(config["input_file"], "r") as f:
        config["run_folders"] = [line.strip() for line in f if line.strip()]
elif config.get("run_folders"):
    # If run_folders is provided as a string, convert it to a list
    if isinstance(config["run_folders"], str):
        config["run_folders"] = [folder.strip() for folder in config["run_folders"].split(",")]
else:
    raise ValueError("No input folders specified. Use either --config input_file=path/to/file or --config run_folders=folder1,folder2,...")

# Validate input folders
config["run_folders"] = [folder for folder in config["run_folders"] if os.path.isdir(folder)]

if not config["run_folders"]:
    raise ValueError("No valid input folders found.")

# Check for demux_snakefile in config
if "demux_snakefile" not in config:
    raise ValueError("demux_snakefile not specified. Use --config demux_snakefile=path/to/demultiplex_snakefile.smk")

if not os.path.isfile(config["demux_snakefile"]):
    raise ValueError(f"Specified demux_snakefile does not exist: {config['demux_snakefile']}")

print(f"Processing the following run folders: {config['run_folders']}")
print(f"Using demultiplex Snakefile: {config['demux_snakefile']}")

# Rule to process all specified run folders
rule all:
    input:
        expand("{run_folder}/demulOutCARC", run_folder=config["run_folders"]),
        expand("{run_folder}/statsCARC/fastqc.log", run_folder=config["run_folders"])

# Rule to process each run folder
rule process_run_folder:
    input:
        run_folder = "{run_folder}"
    output:
        demux = directory("{run_folder}/demulOutCARC"),
        fastqc = "{run_folder}/statsCARC/fastqc.log"
    params:
        demux_snakefile = config["demux_snakefile"]
    run:
        with open(os.path.join(input.run_folder, "config.yaml"), "w") as f:
            f.write(f"run_folder: {input.run_folder}\n")
        
        shell("snakemake -s {params.demux_snakefile} --cores 12 --config run_folder={input.run_folder}")

# Optional: Add a rule to generate a summary report
#rule generate_summary:
#    input:
#        expand("{run_folder}/statsCARC/fastqc.log", run_folder=config["run_folders"])
#    output:
#        "summary_report.txt"
#    run:
#        with open(output[0], "w") as out_f:
#            for log_file in input:
#                run_folder = os.path.dirname(os.path.dirname(log_file))
#                out_f.write(f"Run Folder: {run_folder}\n")
#                with open(log_file, "r") as in_f:
#                    out_f.write(in_f.read())
#                out_f.write("\n\n")