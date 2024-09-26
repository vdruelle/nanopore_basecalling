# Pipeline to basecall the raw data generated from our nanopore
import pathlib
import time
import os
import sys


# pass the path to the folder as and argument when calling snakemake: --congig input=PATH
DATA_DIR = os.path.normpath(config["run_dir"])
INPUT_DIR = os.path.join(DATA_DIR, "raw")
TMP_DIR = os.path.join(DATA_DIR, "tmp")
OUTPUT_DIR = os.path.join(DATA_DIR, "final")
STATISTICS_DIR = os.path.join(DATA_DIR, "statistics")
EXEC_TIME = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())
LOGFILE = os.path.join(DATA_DIR, "basecalling.log")

DORADO_BIN = "softwares/dorado-0.8.0-linux-x64/bin/dorado"
DORADO_MODEL = "softwares/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

# argument to define whether it was a 96 barcode run or not. Omitting this argument will default to 24 barcodes run
if "kit96" in config.keys():
    if config["kit96"] == True:  # if the argument is present and set to true
        NB_BARCODES = 96
    else:
        NB_BARCODES = 24  # if the argument is present and set to false
else:
    NB_BARCODES = 24  # if the argument is not present

NANOPORE_KIT = "SQK-RBK114-" + str(NB_BARCODES)
FLOW_CELL = "FLO-MIN114"
BARCODES = [str(ii).zfill(2) for ii in range(1, NB_BARCODES + 1)]

# create log directory if it does not exists
pathlib.Path("log").mkdir(exist_ok=True)


# Check that DATA_DIR exists
if not os.path.exists(DATA_DIR):
    sys.exit(
        f"Error: The data directory '{DATA_DIR}' does not exist. Make sure you gave the correct path for your folder."
    )

# Check that 'params.tsv' file exists in DATA_DIR
params_file = os.path.join(DATA_DIR, "params.tsv")
if not os.path.isfile(params_file):
    sys.exit(
        f"Error: The 'params.tsv' file is missing in '{DATA_DIR}'. Please create your params file and add it to your run folder."
    )

# Check that 'raw' directory exists in DATA_DIR
if not os.path.isdir(INPUT_DIR):
    sys.exit(
        f"Error: The 'raw' directory is missing in '{DATA_DIR}'. Make sure that your in the right folder and that your raw data folder is correctly named."
    )

# Optionally, check if the 'raw' directory is not empty
if not os.listdir(INPUT_DIR):
    sys.exit(
        f"Error: The 'raw' directory '{INPUT_DIR}' is empty. Please ensure that it contains the raw nanopore files."
    )


localrules:
    all,
    generate_log_file,
    clean,
    clean_all,


rule all:
    input:
        barcodes=expand(os.path.join(OUTPUT_DIR, "barcode_{barcode}.fastq.gz"), barcode=BARCODES),
        unclassified=os.path.join(OUTPUT_DIR, "unclassified.fastq.gz"),
        plot1=os.path.join(STATISTICS_DIR, "len_hist.png"),
        plot2=os.path.join(STATISTICS_DIR, "bp_per_barcode.png"),
        plot3=os.path.join(STATISTICS_DIR, "quality_mean.png"),
        plot4=os.path.join(STATISTICS_DIR, "quality_std.png"),
        clean=os.path.join(DATA_DIR, ".cleaned_dummy_file.txt"),  # comment for debugging


rule generate_log_file:
    output:
        LOGFILE,
    params:
        dorado=DORADO_BIN,
        model=DORADO_MODEL,
        flow_cell=FLOW_CELL,
        kit=NANOPORE_KIT,
        ex_time=EXEC_TIME,
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        python snakecommands.py generate-log-file {output} \
        {params.dorado} \
        {params.model} \
        {params.flow_cell} \
        {params.kit} \
        {params.ex_time}
        cat {DATA_DIR}/params.tsv >> {output}
        """


rule basecall:
    message:
        "Basecalling the reads using Dorado model {params.model} for the kit {params.kit}."
    input:
        input_dir=INPUT_DIR,
        logfile=LOGFILE,
    output:
        directory=directory(TMP_DIR + "/dorado_raw"),
        file=TMP_DIR + "/dorado_raw/basecalled.bam",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    params:
        kit=NANOPORE_KIT,
        model=DORADO_MODEL,
        dorado=DORADO_BIN,
    shell:
        """
        {params.dorado} basecaller {params.model} {input.input_dir} --kit-name {params.kit} > {output.file}
        """


rule demultiplex:
    message:
        "Splitting the reads based on detected barcodes, and removing the barcodes from the reads."
    input:
        rules.basecall.output.file,
    output:
        directory=directory(TMP_DIR + "/barcoded"),
        barcodes=expand(TMP_DIR + "/barcoded/barcode_{barcode}.fastq", barcode=BARCODES),
        unclassified=TMP_DIR + "/barcoded/unclassified.fastq",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    params:
        dorado=DORADO_BIN,
        kit=NANOPORE_KIT,
    threads: 4
    shell:
        """
        mkdir -p {output.directory}
        {params.dorado} demux --output-dir {output.directory} --no-classify {input} -t {threads} --emit-fastq
        cd {output.directory}
        for file in {params.kit}_barcode*.fastq; do mv "$file" "${{file/{params.kit}_barcode/barcode_}}"; done
        for bc in {BARCODES}; do
            if ! [[ -e barcode_$bc.fastq ]]; then
                touch barcode_$bc.fastq
            fi
        done
        """


rule compress:
    message:
        "Generating the final compressed file {output.output_file}."
    input:
        input_file=TMP_DIR + "/barcoded/{filename}.fastq",
    output:
        output_file=OUTPUT_DIR + "/{filename}.fastq.gz",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        gzip -c {input} > {output}
        """


rule stats:
    message:
        "Generating stats for {input.input_file}."
    input:
        input_file=TMP_DIR + "/barcoded/{filename}.fastq",
    output:
        output_file=TMP_DIR + "/stats/{filename}.tsv",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        python snakecommands.py generate-stats {input} {output}
        """


rule combine_stats:
    message:
        "Combining the stat files for all the barcodes."
    input:
        stat_files=expand(TMP_DIR + "/stats/barcode_{barcode}.tsv", barcode=BARCODES),
        stat_file_unclassified=TMP_DIR + "/stats/unclassified.tsv",
    output:
        output_lengths=STATISTICS_DIR + "/lengths.tsv",
        output_quality_mean=STATISTICS_DIR + "/quality.tsv",
        output_quality_std=STATISTICS_DIR + "/quality_std.tsv",
    params:
        nb_barcodes=NB_BARCODES,
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        python snakecommands.py combine-stats {TMP_DIR}/stats {output.output_lengths} {output.output_quality_mean} {output.output_quality_std} {params.nb_barcodes}
        """


rule make_plots_lengths:
    message:
        "Generating lenghts statistics plots."
    input:
        stats_file_lengths=rules.combine_stats.output.output_lengths,
    output:
        len_hist=STATISTICS_DIR + "/len_hist.png",
        bp_per_barcode=STATISTICS_DIR + "/bp_per_barcode.png",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        python snakecommands.py make-plots-lengths {input.stats_file_lengths}
        """


rule make_plots_quality:
    message:
        "Generating quality statistics plots."
    input:
        stats_file_quality_mean=rules.combine_stats.output.output_quality_mean,
        stats_file_quality_std=rules.combine_stats.output.output_quality_std,
    output:
        quality_mean_plot=STATISTICS_DIR + "/quality_mean.png",
        quality_std_plot=STATISTICS_DIR + "/quality_std.png",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        python snakecommands.py make-plots-quality {input.stats_file_quality_mean} {input.stats_file_quality_std}
        """


rule clean:
    message:
        "Cleaning up the output folder."
    input:
        rules.make_plots_lengths.output.len_hist,
        rules.make_plots_quality.output.quality_mean_plot,
        rules.make_plots_quality.output.quality_std_plot,
    output:
        DATA_DIR + "/.cleaned_dummy_file.txt",
    shell:
        """
        rm -rf {TMP_DIR}
        touch {output}
        """


rule clean_all:
    message:
        "Cleaning up the output folder."
    shell:
        """
        rm -rf {OUTPUT_DIR}
        rm -rf {STATISTICS_DIR}
        rm -rf {TMP_DIR}
        rm -rf log
        rm -f {OUTPUT_DIR}/basecalling.log
        rm -f {rules.clean.output}
        """
