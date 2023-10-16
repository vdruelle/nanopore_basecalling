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
Output dir: {data_dir / "basecalled"}

Parameter file:

"""

    with open(logfile, "w") as file:
        file.write(log_text) # Writes the infos from the log_text



@cli.command()
@click.argument("reads", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def generate_stats(reads, output):
    from Bio import SeqIO

    with open(reads, 'r') as f, open(output, 'w') as file:
        for record in SeqIO.parse(f, 'fastq'):
            length = str(len(record)) + "\n"
            file.write(length)


@cli.command()
@click.argument("stats_dir", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def combine_stats(stats_dir, output):
    "Load reads all csv files in the stats_dir and combine them into a single csv file"
    import os

    import pandas as pd
    
    file_list = [file for file in os.listdir(stats_dir) if file.endswith(".tsv")]
    combined_df = pd.DataFrame()

    # Iterate over each CSV file
    for file_name in file_list:
        file_path = os.path.join(stats_dir, file_name)
        header = os.path.splitext(file_name)[0]

        # Read the CSV file into a DataFrame
        df = pd.read_csv(file_path, names=[header])

        # Assign the DataFrame to the column with the header as the column label
        combined_df = pd.concat([combined_df,df], axis=1)

    sorted_headers = [f"barcode_{str(ii).zfill(2)}" for ii in range(1, 25)]
    sorted_headers += ["unclassified"]
    combined_df = combined_df[sorted_headers]
    combined_df.to_csv(output, index=False, sep="\t")

@cli.command()
@click.argument("stats_file", type=click.Path(exists=True))
def make_plots(stats_file):
    import os

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    df = pd.read_csv(stats_file, sep="\t")
    output = os.path.split(stats_file)[:-1][0]

    # log-length distribution by barcode, normalized
    plt.figure()
    sns.violinplot(data=df, orient="v", log_scale=True)
    plt.yscale("log")
    plt.ylabel("Length of reads")
    plt.xticks(rotation=80)
    plt.tight_layout()
    plt.savefig(f"{output}/len_hist.png", facecolor="w", dpi=200)

    sum_values = df.sum() / 1e6

    # Create the plot using seaborn
    plt.figure()
    sns.barplot(x=sum_values.index, y=sum_values.values)  # Create the bar plot
    plt.xlabel("Barcode")  # Set the x-axis label
    plt.ylabel("MBp")  # Set the y-axis label
    plt.title("Number of Basepairs for Each Barcode")  # Set the title of the plot
    plt.xticks(rotation=80)  # Rotate the x-axis labels if needed
    plt.tight_layout()
    plt.savefig(f"{output}/bp_per_barcode.png", facecolor="w", dpi=200)


if __name__ == "__main__":
    cli()