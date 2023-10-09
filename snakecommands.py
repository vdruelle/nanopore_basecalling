import os
from glob import glob

import click


@click.group()
def cli():
    pass


@cli.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path())
def concatenate(input_dir, output_dir):
    # Get all the paths in directory
    paths = glob(os.path.join(input_dir, "*"))

    for path in paths:
        if os.path.isdir(path):
            # Get all the .fastq files in the folder
            fastq_files = glob(os.path.join(path, "*.fastq"))

            # Concatenate the .fastq files into one
            # if os.path.basename(path) != "unclassified":
            output_file = os.path.join(output_dir, os.path.basename(path) + ".fastq")
            concatenate_command = f"cat {' '.join(fastq_files)} > {output_file}"
            os.system(concatenate_command)

            # # Compress the concatenated file
            compressed_file = os.path.join(path, "concatenated.fastq.gz")
            compress_command = f"gzip {output_file} -c > {compressed_file}"
            os.system(compress_command)


@cli.command()
@click.argument("data_dir", type=click.Path(exists=True))
@click.argument("dorado_bin", type=click.Path(exists=True))
@click.argument("time", type=click.Path())
def generate_log_file(data_dir, dorado_bin, time):
    import subprocess
    logfile = data_dir + "/params_" + time + ".log"

    reporemote = subprocess.check_output("git remote -v | head -n 1", shell=True).decode()
    commitID = subprocess.check_output("git rev-parse HEAD", shell=True).decode()
    doradover = subprocess.check_output(dorado_bin + " -v 2>&1", shell=True).decode()

    log_text = f"""
    Log-file for the basecalling executed by the Snakemake script.
    Execution time: {time}

    The code is stored in the repository: {reporemote}
    The current commit is: {commitID}
    Dorado version: {doradover}
    Input dir: {os.path.join(data_dir, "raw")}
    Output dir: {os.path.join(data_dir, "basecalled")}

    Parameter file:

    """

    with open(logfile, "w") as file:
        file.write(log_text) # Writes the infos from the log_text


if __name__ == "__main__":
    cli()