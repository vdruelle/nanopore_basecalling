# Pipeline to basecall the raw data generated from our nanopore

DATA_DIR = "test_data"
# TMP_DIR = DATA_DIR + "/tmp"
TMP_DIR = os.environ.get(
    "TMPDIR", DATA_DIR + "/tmp"
)  # Default to '/tmp' if $TMPDIR is not set
OUTPUT_DIR = DATA_DIR + "/basecalled"
PARAMETER_FILE = DATA_DIR + "/params.tsv"  # Used as a log file for the run
BARCODES = [str(i).zfill(2) for i in range(1, 25)]

# Path of Dorado binary
dorado_bin = "softwares/dorado-0.3.2-linux-x64/bin/dorado"
guppy_barcoder_bin = "softwares/ont-guppy-cpu/bin/guppy_barcoder"

# --------- parse parameters from file ---------


# # Function to read parameter file and create a dictionary
# def read_PARAMETER_FILE(file_path):
#     params = {}
#     with open(file_path) as f:
#         for line in f:
#             key, value = line.strip().split("\t")
#             params[key] = value
#     return params


# par_dict = read_PARAMETER_FILE(PARAMETER_FILE)

# # --------- produce log file ---------

# log_filename = f"{par_dict['parPrefix']}_{par_dict['timeNow']}.log"

# log_text = f"""
# Log-file for the basecalling executed by the Snakemake script.
# Execution time     : {par_dict['timeNow']}
# Snakemake run label: {snakemake.workflow.run_name}

# ------- CODE -------
# The code is stored in the repository: REPOREMOTE
# The current commit is: COMMITID

# ------- DORADO -------
# Dorado version: DORADOVER

# ------- BASECALLING PARAMETERS -------
# Parameter file: {PARAMETER_FILE}
# Barcodes: {par_dict['barcode_id']}
# Flowcell ID: {par_dict['flow_cell_id']}
# Flowcell type: {par_dict['flow_cell_type']}
# Kit: {par_dict['ligation_kit']}
# Dorado config: {par_dict['dorado_config_file']}
# Barcode kits: {par_dict['barcode_kits']}
# Nanopore data root dir: {par_dict['nanopore_data_root_dir']}

# ------- INPUT / OUTPUT DIRECTORIES -------
# Input dir: {INPUT_DIR}
# Output dir: {OUTPUT_DIR}
# """


# # Process that generates the log file
# rule generate_log_file:
#     output:
#         log_filename,
#     shell:
#         """
#         echo "{log_text}" > {log_filename}
#         REMOTE=$(git remote -v | head -n 1)
#         COMM=$(git rev-parse HEAD)
#         dorado=$({dorado_bin} -v | head -n 1)
#         sed -i "s|REPOREMOTE|${{REMOTE}}|g" {log_filename}
#         sed -i "s|COMMITID|${{COMM}}|g" {log_filename}
#         sed -i "s|DORADOVER|${{dorado}}|g" {log_filename}
#         """


# --------- workflow ---------


rule all:
    input:
        barcodes=expand(
            OUTPUT_DIR + "/barcode{barcode}.fastq.gz",
            barcode=BARCODES,
        ),


rule basecall:
    message:
        "Basecalling the reads using Dorado model {params.model}."
    input:
        DATA_DIR + "/raw",
    output:
        directory=directory(TMP_DIR + "/dorado_raw"),
        file=TMP_DIR + "/dorado_raw/basecalled.fastq",
    params:
        model="softwares/dorado-0.3.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v4.2.0",
    shell:
        """
        {dorado_bin} basecaller {params.model} {input} --emit-fastq > {output.file}
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
        {guppy_barcoder_bin} \
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
        echo lol
        """


rule clean:
    message:
        "Cleaning up the output folder."
    shell:
        """
        rm -rf {DATA_DIR}/final
        """
