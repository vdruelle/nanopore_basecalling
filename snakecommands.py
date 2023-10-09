import os
from glob import glob

import click


@click.group()
def cli():
    pass


@cli.command()
@click.argument("input_dir", type=click.Path(exists=True))
# @click.argument("output_dir", type=click.Path(exists=True))
@click.argument("output_dir")
def concatenate(input_dir, output_dir):
    # os.system(f"mkdir {output_dir}")

    # Get all the paths in directory
    paths = glob(os.path.join(input_dir, "*"))

    for path in paths:
        if os.path.isdir(path):
            # Get all the .fastq files in the folder
            fastq_files = glob(os.path.join(path, "*.fastq"))

            # Concatenate the .fastq files into one
            # if os.path.basename(path) != "unclassified":
            output_file = os.path.join(
                output_dir, os.path.basename(path) + ".fastq"
            )
            concatenate_command = f"cat {' '.join(fastq_files)} > {output_file}"
            os.system(concatenate_command)

            # # Compress the concatenated file
            compressed_file = os.path.join(path, "concatenated.fastq.gz")
            compress_command = f"gzip {output_file} -c > {compressed_file}"
            os.system(compress_command)


if __name__ == "__main__":
    cli()
