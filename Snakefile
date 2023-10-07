# Pipeline to basecall the raw data generated from our nanopore

# --------- workflow parameters ---------

# Define input and output directories
DATA_DIR = "test_data"
INPUT_DIR = DATA_DIR + "/raw_big"
OUTPUT_DIR = DATA_DIR + "/basecalled"
PARAMETER_FILE = DATA_DIR + "/params.tsv"  # Used as a log file for the run

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


# # Function to generate input file paths
# def input_files(wildcards):
#     return expand(
#         "{INPUT_DIR}/{pod5_file}", INPUT_DIR=INPUT_DIR, pod5_file=wildcards.pod5_file
#     )


# # Function to generate output file paths
# def output_files(wildcards):
#     return expand(
#         "{OUTPUT_DIR}/{barcode}/fastq_pass/{{samples}}.fastq.gz",
#         OUTPUT_DIR=OUTPUT_DIR,
#         barcode=wildcards.barcode,
#         samples=par_dict["barcodes"].split(","),
#     )


rule all:
    input:
        OUTPUT_DIR,


rule basecall:
    message:
        "Basecalling the reads using Dorado model {params.model}."
    input:
        INPUT_DIR,
    output:
        directory=directory(DATA_DIR + "/dorado_raw"),
        file=DATA_DIR + "/dorado_raw/basecalled.fastq",
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
        directory=directory(OUTPUT_DIR),
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


rule clean:
    message:
        "Cleaning up the output folder."
    shell:
        """
        rm -rf {OUTPUT_DIR}
        rm -rf {DATA_DIR}/dorado_raw
        """


# # Process that concatenates and compresses fastq files by barcode
# rule concatenate_and_compress:
#     input:
#         barcode=lambda wildcards: wildcards.barcode,
#         reads=glob_wildcards(output_files).samples,
#     output:
#         "{barcode}.fastq.gz",
#     shell:
#         """
#         cat {{reads}} | gzip -c > {output}
#         chmod 444 {output}
#         """
# # Process for creating a CSV file with basecalling statistics
# rule basecalling_live_report:
#     input:
#         bc_files=expand("{par_dict['parDir']}/{bcstats_filename}", par_dict=par_dict),
#         reads=glob_wildcards(output_files).samples,
#     output:
#         "{par_dict['parDir']}/{bcstats_filename}",
#     run:
#         # Decompress and generate stats (modify as needed)
#         shell(
#             "gzip -dc {{reads}} > reads.fastq && "
#             "python3 {snakemake.wildcards.par_dict['parDir']}/basecall_stats.py reads.fastq && "
#             "tail -n +2 basecalling_stats.csv >> {output} && "
#             "rm reads.fastq"
#         )
# # Create a list of all barcodes
# barcodes = par_dict["barcode_id"].split(",")
# Workflow
# rule all:
#     input:
#         expand("concatenate_and_compress/{barcode}.fastq.gz", barcode=barcodes),
