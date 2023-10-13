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

BARCODES = [str(ii).zfill(2) for ii in range(1, 25)]
DORADO_BIN = "softwares/dorado-0.4.1-linux-x64/bin/dorado"
DORADO_MODEL = "softwares/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
NANOPORE_KIT = "SQK-RBK114-24"
FLOW_CELL = "FLO-MIN114"

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
        cat test_data/params.tsv >> {output}
        """


rule basecall:
    message:
        "Basecalling the reads using Dorado model {params.model}."
    input:
        input_dir=INPUT_DIR,
        # logfile=LOGFILE,
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
        barcodes=expand(TMP_DIR + "/barcoded/barcode_{barcode}.bam", barcode=BARCODES),
        unclassified=TMP_DIR + "/barcoded/unclassified.bam",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    params:
        dorado=DORADO_BIN,
    threads: 4
    shell:
        """
        mkdir -p {output.directory}
        {params.dorado} demux --output-dir {output.directory} --no-classify {input} -t {threads}
        cd {output.directory}
        for file in SQK-RBK114-24_barcode*.bam; do mv "$file" "${{file/SQK-RBK114-24_barcode/barcode_}}"; done
        """


rule bam_to_fastq:
    message:
        "Converting the barcoded reads to fastq."
    input:
        barcodes=TMP_DIR + "/barcoded/{filename}.bam",
        unclassified=rules.demultiplex.output.unclassified,
    output:
        barcodes=TMP_DIR + "/fastq/{filename}.fastq",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        samtools fastq {input.barcodes} > {output.barcodes}
        """


rule compress:
    message:
        "Generating the final compressed file {output.output_file}."
    input:
        input_file=TMP_DIR + "/fastq/{filename}.fastq",
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
        input_file=TMP_DIR + "/fastq/{filename}.fastq",
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
        output_file=STATISTICS_DIR + "/statistics.tsv",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        python snakecommands.py combine-stats {TMP_DIR}/stats {output}
        """


rule make_plots:
    message:
        "Generating plots from {input.stats_file}."
    input:
        stats_file=rules.combine_stats.output.output_file,
    output:
        len_hist=STATISTICS_DIR + "/len_hist.png",
        bp_per_barcode=STATISTICS_DIR + "/bp_per_barcode.png",
    conda:
        "conda_envs/nanopore_basecalling.yml"
    shell:
        """
        python snakecommands.py make-plots {input.stats_file}
        """


rule clean:
    message:
        "Cleaning up the output folder."
    input:
        rules.make_plots.output.len_hist,
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
