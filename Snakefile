# Pipeline to basecall the raw data generated from our nanopore

import time
import os

DATA_DIR = "test_data"
INPUT_DIR = DATA_DIR + "/raw"
# TMP_DIR = DATA_DIR + "/tmp"
TMP_DIR = os.environ.get("TMPDIR", DATA_DIR + "/tmp")  # Default to '/tmp' if $TMPDIR is not set
OUTPUT_DIR = DATA_DIR + "/basecalled"
PARAMETER_FILE = DATA_DIR + "/params.tsv"  # Used as a log file for the run
BARCODES = [str(i).zfill(2) for i in range(1, 25)]

# Path of Dorado binary
DORADO_BIN = "softwares/dorado-0.3.2-linux-x64/bin/dorado"
GUPPY_BARCODER_BIN = "softwares/ont-guppy-cpu/bin/guppy_barcoder"


# --------- workflow ---------


# rule all:
#     input:
# barcodes=expand(
#     OUTPUT_DIR + "/barcode{barcode}.fastq.gz",
#     barcode=BARCODES,
# ),

TIME = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())
# TIME = "1h"


rule generate_log_file:
    output:
        DATA_DIR + "/params_" + TIME + ".log",
    shell:
        f"""
        python snakecommands.py generate-log-file {DATA_DIR} {DORADO_BIN} {TIME}
        cat test_data/params.tsv >> {"test_data/params_" + TIME + ".log"}
        """


rule basecall:
    message:
        "Basecalling the reads using Dorado model {params.model}."
    input:
        INPUT_DIR,
    output:
        directory=directory(TMP_DIR + "/dorado_raw"),
        file=TMP_DIR + "/dorado_raw/basecalled.fastq",
    params:
        model="softwares/dorado-0.3.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v4.2.0",
    shell:
        """
        {DORADO_BIN} basecaller {params.model} {input} --emit-fastq > {output.file}
        """


rule demultiplex:
    message:
        "Splitting the reads based on detected barcodes, and removing the barcodes from the reads."
    input:
        rules.basecall.output.directory,
    output:
        directory=directory(TMP_DIR + "/barcoded"),
    params:
        kit="SQK-RBK114-24",
    threads: 16
    shell:
        """
        {GUPPY_BARCODER_BIN} \
        -i {input} \
        -s {output} \
        --barcode_kits {params.kit} \
        --enable_trim_barcodes \
        --detect_mid_strand_barcodes \
        --num_barcoding_threads {threads}
        """


rule concatenate:
    message:
        "Concatenating the data for each barcode into a single file."
    input:
        input_dir=rules.demultiplex.output.directory,
    output:
        barcodes=expand(
            TMP_DIR + "/concatenated/barcode{barcode}.fastq",
            barcode=BARCODES,
        ),
        unclassified=TMP_DIR + "/concatenated/unclassified.fastq",
        output_dir=directory(TMP_DIR + "/concatenated"),
    shell:
        """
        python snakecommands.py concatenate {input.input_dir} {output.output_dir}
        """


rule compress:
    message:
        "Generating the final file {output.output_file}."
    input:
        input_file=TMP_DIR + "/concatenated/barcode{barcode}.fastq",
    output:
        output_file=OUTPUT_DIR + "/barcode{barcode}.fastq.gz",
    shell:
        """
        gzip -c {input} > {output}
        """


rule clean:
    message:
        "Cleaning up the output folder."
    shell:
        """
        rm -rf {DATA_DIR}/final
        """
