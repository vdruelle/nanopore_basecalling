# Pipeline to basecall the raw data generated from our nanopore
import pathlib
import time
import os


# pass the path to the folder as and argument when calling snakemake: --congig input=PATH
DATA_DIR = config["run_dir"]
INPUT_DIR = DATA_DIR + "/raw"
TMP_DIR = DATA_DIR + "/tmp"
OUTPUT_DIR = DATA_DIR + "/final"
STATISTICS_DIR = DATA_DIR + "/statistics"
EXEC_TIME = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())
LOGFILE = DATA_DIR + "/basecalling.log"

DORADO_BIN = "softwares/dorado-0.5.0-linux-x64/bin/dorado"
DORADO_MODEL = "softwares/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0"

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


localrules:
    all,
    generate_log_file,
    clean,
    clean_all,


rule all:
    input:
        barcodes=expand(OUTPUT_DIR + "/barcode_{barcode}.fastq.gz", barcode=BARCODES),
        unclassified=OUTPUT_DIR + "/unclassified.fastq.gz",
        plot1=STATISTICS_DIR + "/len_hist.png",
        plot2=STATISTICS_DIR + "/bp_per_barcode.png",
        plot3=STATISTICS_DIR + "/quality_mean.png",
        plot4=STATISTICS_DIR + "/quality_std.png",
        clean=DATA_DIR + "/.cleaned_dummy_file.txt",  # comment for debugging


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


# download and extract dorado binaries. The URL is in the config file.
# rule download_dorado_bin:
#     output:
#         "softwares/dorado"
#     params:
#         url=config["dorado_url"]
#     shell:
#         """
#         mkdir -p softwares
#         wget {params.url} --directory softwares
#         BN=$(basename {params.url} .tar.gz)
#         FN=softwares/$BN.tar.gz
#         ln -s $BN/bin/dorado softwares/dorado
#         tar -xf $FN -C softwares --overwrite
#         rm $FN
#         """
# # download desired dorado model, as per config file specification
# rule download_dorado_model:
#     input:
#         drd=rules.download_dorado_bin.output
#     output:
#         f"softwares/dorado_models/{config['dorado_model']}"
#     params:
#         mdl=config["dorado_model"]
#     shell:
#         """
#         mkdir -p software/dorado_models
#         {input.drd} download --model {params.mdl} --directory {output}
#         """
