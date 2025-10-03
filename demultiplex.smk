import os
from datetime import datetime

#configfile: "config.yaml"

run_folder = config["run_folder"]
proj_name = os.path.basename(run_folder)
time = datetime.now().strftime("%d-%m-%Y-%H-%M")
ss = os.path.join(run_folder, "SampleSheet.csv")

rule all:
    input:
        "statsCARC/fastqc{time}.log".format(time=time),
        "demulLogsCARC/demultiplexCARC_{time}.log".format(time=time),
        "multiqcCARCPass"

rule check_samplesheet:
    input:
        samplesheet = ss
    output:
        ss = os.path.join(run_folder, "SampleSheetCARC.back"),
        checkpass = touch("check_samplesheet.done")
    params:
        run_folder = run_folder,
        proj_name = proj_name
    log:
        "logs/demulLogs/check.log"
    shell:
        """
        cd {params.run_folder}

        mkdir -p demulOutCARC demulLogsCARC

        cp -v {input.samplesheet} {output.ss}
        curl -X POST -H "Content-type: application/json" --data "{{\\"text\\": \\"{params.proj_name} Demultiplex InQueue\\"}}" 'https://hooks.slack.com/services/###'
        """

rule demultiplex:
    input:
        samplesheet = ss,
        checkpass = "check_samplesheet.done"
    output:
        log = "demulLogsCARC/demultiplexCARC_{time}.log".format(time=time),
        flag = touch("demultiplexCARCPass")
    params:
        proj_name = proj_name
    log:
        "logs/demulLogs/demultiplex.log"
    shell:
       """
        module purge
        #module load usc
        #module load intel/18.0.4 intel/19.0.4 bcl2fastq2/2.20.0.422

        module load gcc/13.3.0 bcl2fastq2/2.20.0.422

        mkdir -p demulOutCARC demulLogsCARC


        bcl2fastq --loading-threads 4 --processing-threads 10 --writing-threads 4 --barcode-mismatches 0 -o demulOutCARC &> {output.log}

        if [ $? -eq 0 ]; then
            curl -X POST -H "Content-type: application/json" --data "{{\\"text\\": \\"{params.proj_name} Demultiplex Pass\\"}}" 'https://hooks.slack.com/services/###'
        else
            curl -X POST -H "Content-type: application/json" --data "{{\\"text\\": \\\"!!!{params.proj_name} Demultiplex Fail\\"}}" 'https://hooks.slack.com/services/###'
            exit 1
        fi
        """

rule fastqc:
    input:
        flag = "demultiplexCARCPass",
    output:
        outdir = directory("statsCARC/fastqcOut"),
        log = "statsCARC/fastqc{time}.log".format(time=time),
        flag = touch("fastqcCARCPass")
    params:
        proj_name = proj_name
    log:
        "logs/demulLogs/fastqc.log"
    shell:
        """
        mkdir -p statsCARC
        mkdir -p statsCARC/fastqcOut

        module load usc fastqc
        fastqc --threads 16 --outdir {output.outdir} `find ./demulOutCARC -type f -name "*.fastq.gz" | tr "\\n" " "` &> {output.log}

        if [ $? -eq 0 ]; then
            curl -X POST -H "Content-type: application/json" --data "{{\\"text\\": \\"{params.proj_name} FastQC Pass\\"}}" 'https://hooks.slack.com/services/###'
        else
            curl -X POST -H "Content-type: application/json" --data "{{\\"text\\": \\\"!!!{params.proj_name} FastQC Fail\\"}}" 'https://hooks.slack.com/services/###'
            exit 1
        fi
        """

rule multiQC:
    input:
        flag = "fastqcCARCPass",
        fastqc_log = "statsCARC/fastqc{time}.log".format(time=time)
    output:
        outdir = directory("statsCARC/multiqc"),
        flag = touch("multiqcCARCPass")
    params:
        proj_name = proj_name
    log:
        "logs/demulLogs/multiqc.log"
    shell:
        """
        mkdir -p statsCARC/multiqc

        export PYTHONPATH=/project/davidwcr_264/Packages/pythonModule/python-3.9.2
        export OPENBLAS_NUM_THREADS=4

        /project/davidwcr_264/Packages/pythonModule/python-3.9.2/bin/multiqc \
            --interactive \
            -n multiqc_demul-fastqc_CARC \
            -o {output.outdir} \
            demulOutCARC/Stats/ \
            statsCARC/fastqcOut/

        if [ $? -eq 0 ]; then
            curl -X POST -H "Content-type: application/json" --data "{{\\"text\\": \\"{params.proj_name} multiQC Pass\\"}}" 'https://hooks.slack.com/services/###'
        else
            curl -X POST -H "Content-type: application/json" --data "{{\\"text\\": \\\"!!!{params.proj_name} multiQC Fail\\"}}" 'https://hooks.slack.com/services/###'
            exit 1
        fi
        """
