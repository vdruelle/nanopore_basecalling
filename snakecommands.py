import click


@click.group()
def cli():
    pass


@cli.command()
@click.argument("log_file", type=click.Path())
@click.argument("dorado_bin", type=click.Path(exists=True))
@click.argument("dorado_model", type=str)
@click.argument("flow_cell", type=str)
@click.argument("nanopore_kit", type=str)
@click.argument("time", type=str)
def generate_log_file(log_file, dorado_bin, dorado_model, flow_cell, nanopore_kit, time):
    import pathlib
    import subprocess

    logfile = pathlib.Path(log_file)
    data_dir = logfile.parent
    data_dir.mkdir(exist_ok=True)

    reporemote = subprocess.check_output("git remote -v | head -n 1", shell=True).decode()
    commitID = subprocess.check_output("git rev-parse HEAD", shell=True).decode()
    doradover = subprocess.check_output(dorado_bin + " -v 2>&1", shell=True).decode()

    log_text = f"""
Log-file for the basecalling executed by the Snakemake script.
Execution time: {time}

The code is stored in the repository: {reporemote}
The current commit is: {commitID}
Dorado version: {doradover}
Dorado model: {dorado_model}
Flow cell: {flow_cell}
Nanopore kit: {nanopore_kit}
Input dir: {data_dir / "raw"}
Output dir: {data_dir / "final"}

Parameter file:

"""

    with open(logfile, "w") as file:
        file.write(log_text) # Writes the infos from the log_text



@cli.command()
@click.argument("reads", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def generate_stats(reads, output):
    """
    Generate a stats file in .tsv format for a fastq file. 
    It contains the length, mean quality and standard deviation of the quality.
    """
    import numpy as np
    from Bio import SeqIO

    with open(reads, 'r') as f, open(output, 'w') as file:
        file.write("length\tmean_quality\tstd_quality\n")
        for record in SeqIO.parse(f, 'fastq'):
            length = len(record)
            mean_quality = np.mean(record.letter_annotations['phred_quality'])
            std_quality = np.std(record.letter_annotations['phred_quality'])
            tmp_string = str(length) + "\t" + str(mean_quality) + "\t" + str(std_quality) + "\n"
            file.write(tmp_string)


@cli.command()
@click.argument("stats_dir", type=click.Path(exists=True))
@click.argument("output_lengths", type=click.Path())
@click.argument("output_quality_mean", type=click.Path())
@click.argument("output_quality_std", type=click.Path())
@click.argument("nb_barcodes", type=int)
def combine_stats(stats_dir, output_lengths, output_quality_mean, output_quality_std, nb_barcodes):
    "Reads all csv files in the stats_dir and combine them into .tsv files (one for lengths, one for mean quality and one for std quality)."
    import os

    import pandas as pd
    
    file_list = [file for file in os.listdir(stats_dir) if file.endswith(".tsv")]
    df_lengths = pd.DataFrame()
    df_quality_mean = pd.DataFrame()
    df_quality_std = pd.DataFrame()

    # Iterate over each CSV file
    for file_name in file_list:
        file_path = os.path.join(stats_dir, file_name)
        header = os.path.splitext(file_name)[0]

        # Read the CSV file into a DataFrame
        df = pd.read_table(file_path)

        # Assign the DataFrame to the column with the header as the column label
        df_lengths = pd.concat([df_lengths,df["length"]], axis=1)
        df_lengths = df_lengths.rename(columns={"length": header})

        df_quality_mean = pd.concat([df_quality_mean,df["mean_quality"]], axis=1)
        df_quality_mean = df_quality_mean.rename(columns={"mean_quality": header})

        df_quality_std = pd.concat([df_quality_std,df["std_quality"]], axis=1)
        df_quality_std = df_quality_std.rename(columns={"std_quality": header})

    sorted_headers = [f"barcode_{str(ii).zfill(2)}" for ii in range(1, nb_barcodes+1)]
    sorted_headers += ["unclassified"]
    df_lengths = df_lengths[sorted_headers]
    df_quality_mean = df_quality_mean[sorted_headers]
    df_quality_std = df_quality_std[sorted_headers]

    df_lengths.to_csv(output_lengths, index=False, sep="\t")
    df_quality_mean.to_csv(output_quality_mean, index=False, sep="\t")
    df_quality_std.to_csv(output_quality_std, index=False, sep="\t")

@cli.command()
@click.argument("stats_file", type=click.Path(exists=True))
def make_plots_lengths(stats_file):
    "Generates the length statistics plots."
    import os

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    df = pd.read_csv(stats_file, sep="\t")
    output = os.path.split(stats_file)[:-1][0]

    # log-length distribution by barcode, normalized
    plt.figure(figsize=(12,8))
    sns.violinplot(data=df, orient="v", log_scale=True)
    plt.yscale("log")
    plt.ylabel("Length of reads")
    plt.xlabel("Barcode")  # Set the x-axis label
    plt.xticks(range(len(df.columns)), [f"{ii:02d}" for ii in range(len(df.columns))], rotation=80, fontsize="small")
    plt.tight_layout()
    plt.savefig(f"{output}/len_hist.png", facecolor="w", dpi=200)


    sum_values = df.sum() / 1e6

    # Bar plot of the total number of bp per barcode
    plt.figure(figsize=(12,8))
    sns.barplot(x=sum_values.index, y=sum_values.values)  # Create the bar plot
    plt.xlabel("Barcode")  # Set the x-axis label
    plt.ylabel("MBp")  # Set the y-axis label
    plt.title("Number of Basepairs for Each Barcode")  # Set the title of the plot
    plt.xticks(range(len(df.columns)), [f"{ii:02d}" for ii in range(len(df.columns))], rotation=80, fontsize="small")
    plt.tight_layout()
    plt.savefig(f"{output}/bp_per_barcode.png", facecolor="w", dpi=200)

@cli.command()
@click.argument("stats_file_mean", type=click.Path(exists=True))
@click.argument("stats_file_std", type=click.Path(exists=True))
def make_plots_quality(stats_file_mean, stats_file_std):
    "Generates the quality statistics plots."
    import os

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    df_mean = pd.read_csv(stats_file_mean, sep="\t")
    df_std = pd.read_csv(stats_file_std, sep="\t")
    output_dir = os.path.split(stats_file_mean)[:-1][0]

    plt.figure(figsize=(12,8))
    sns.violinplot(data=df_mean, orient="v")
    plt.ylabel("Mean quality")
    plt.xticks(range(len(df_mean.columns)), [f"{ii:02d}" for ii in range(len(df_mean.columns))], rotation=80, fontsize="small")
    plt.xlabel("Barcode")  # Set the x-axis label
    plt.tight_layout()
    plt.savefig(f"{output_dir}/quality_mean.png", facecolor="w", dpi=200)

    plt.figure(figsize=(12,8))
    sns.violinplot(data=df_std, orient="v")
    plt.ylabel("Std of quality")
    plt.xticks(range(len(df_mean.columns)), [f"{ii:02d}" for ii in range(len(df_mean.columns))], rotation=80, fontsize="small")
    plt.xlabel("Barcode")  # Set the x-axis label
    plt.tight_layout()
    plt.savefig(f"{output_dir}/quality_std.png", facecolor="w", dpi=200)

if __name__ == "__main__":
    cli()