# nanopore_basecalling

## Introduction
This pipeline is used to basecall raw data generated from a nanopore sequencing run (pod5). It performs the following steps:

- Generate a log file with information about the basecalling process.
- Basecall the reads using Dorado.
- Split the reads based on detected barcodes and trim the barcodes from the reads (Dorado default behaviour).
- Convert the barcoded reads to FASTQ format.
- Generate statistics for the barcoded reads.
- Combine the statistics files for all the barcodes.
- Generate plots from the combined statistics file.
- Clean up the temporary files created along the way to avoid unecessary storage use.


## How to run

### Prerequisite
Start by cloning this repo (on the cluster, if aiming for cluster execution):
```
git clone https://github.com/vdruelle/nanopore_basecalling.git
```

Once this is done, you need to download dorado (https://github.com/nanoporetech/dorado, linux-x64 in our case), unpack it and move the folder called `dorado-0.4.1-linux-x64` to the directory `softwares/`.

You can then download the appropriate dorado model into the directory `softwares/doarado_models` by typing:
```
./dorado-0.4.1-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0 --directory sofwares/dorado_models
```

Last you need to create the conda environment for the pipeline:
```
conda env create -f conda_env/nanopore_basecalling.yml
```

### Running the pipeline
#### On the cluster (Scicore / SLURM) - recommended
We suggest running the pipeline from a tmux session. To do so start by launching a tmux session, then move to where you saved the repository and activate the conda environment with:

```
conda activate nanopore_basecalling
```

Your run folder must contain a subfolder named `raw` in which all your `.pod5` are located as well as the `params.tsv` file which you need to update. This is used for logging to keep track of the run. Start the pipeline, giving the path to your run folder as argument.
```
snakemake --profile cluster --config run_dir=<path to your run folder>
```
This command will launch the pipeline by submitting the appropriate jobs for cluster execution.
You can monitor the progress of the pipeline in the console output. For a good nanopore run (20Gbp), the pipeline should take around 1h30 to complete.

#### Local execution
If running locally (which necessitate a strong GPU, unless using faster models), also start by activating the conda environment. Then launch the pipeline with:

```
snakemake --config run_dir=<path to your run folder> --cores <number of cores you want to use>
```

### Output
Once the pipeline completes, should have more files and directories in your run folder. You should have:
- A directory `final`. It contains separate `.fastq` files for all the barcodes.
- A directory `statistics`. It contains a `.tsv` file for some statistics about the run, as well as two figures that can be used to get an idea of how the run went.
- A log file `basecalling.log`. This file describe when and with which parameter the basecalling happened.

The pipeline will also generate many intermediate outputs while it runs. They are automatically removed at the end.