# Pipeline to basecall the raw data generated from our nanopore

import time
import os

DATA_DIR = "test_data"
INPUT_DIR = DATA_DIR + "/raw"
TMP_DIR = DATA_DIR + "/tmp"
# TMP_DIR = os.environ.get("TMPDIR", DATA_DIR + "/tmp")  # Default to '/tmp' if $TMPDIR is not set
OUTPUT_DIR = DATA_DIR + "/final"
STATISTICS_DIR = DATA_DIR + "/statistics"
BARCODES = [str(ii) for ii in range(1, 25)]
TIME = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())
LOGFILE = DATA_DIR + "/params_" + TIME + ".log"

# Path of Dorado binary
DORADO_BIN = "softwares/dorado-0.4.0-linux-x64/bin/dorado"


rule all:
    input:
        barcodes=expand(OUTPUT_DIR + "/barcode_{barcode}.fastq.gz", barcode=BARCODES),
        unclassified=OUTPUT_DIR + "/unclassified.fastq.gz",
        plot1=STATISTICS_DIR + "/len_hist.png",
        plot2=STATISTICS_DIR + "/bp_per_barcode.png",


rule generate_log_file:
    output:
        LOGFILE,
    shell:
        f"""
        python snakecommands.py generate-log-file {DATA_DIR} {DORADO_BIN} {TIME}
        cat test_data/params.tsv >> {"test_data/params_" + TIME + ".log"}
        """


rule basecall:
    message:
        "Basecalling the reads using Dorado model {params.model}."
    input:
        input_dir=INPUT_DIR,
        logfile=LOGFILE,
    output:
        directory=directory(TMP_DIR + "/dorado_raw"),
        file=TMP_DIR + "/dorado_raw/basecalled.bam",
    params:
        model="softwares/dorado-0.4.0-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v4.2.0",
        kit="SQK-RBK114-24",
    shell:
        """
        {DORADO_BIN} basecaller {params.model} {input.input_dir} --kit-name {params.kit} > {output.file}
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
    params:
        kit="SQK-RBK114-24",
    threads: 4
    shell:  # TODO make this more robust for unclassified barcodes
        """
        mkdir -p {output.directory}
        samtools split {input} -f '{output.directory}/barcode_%#.bam' --threads {threads}
        mv {output.directory}/barcode_0.bam {output.directory}/unclassified.bam
        """


rule bam_to_fastq:
    message:
        "Converting the barcoded reads to fastq."
    input:
        barcodes=TMP_DIR + "/barcoded/{filename}.bam",
        unclassified=rules.demultiplex.output.unclassified,
    output:
        barcodes=TMP_DIR + "/fastq/{filename}.fastq",
    threads: 1
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
    threads: 1
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
    threads: 1
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
    threads: 1
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
    threads: 1
    shell:
        """
        python snakecommands.py make-plots {input.stats_file}
        """


rule clean:
    message:
        "Cleaning up the output folder."
    shell:
        """
        rm -rf {OUTPUT_DIR}
        rm -rf {STATISTICS_DIR}
        rm -rf {TMP_DIR}
        """
