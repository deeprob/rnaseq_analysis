# Running the rnaseq pipeline on SLURM using singularity

## Run data preparation using the provided shell script

1. Edit the *project_dir* variable in *prepare.sh* to your directory of choice. This directory will host all files required to run reads to counts part of the pipeline. __Make sure it has sufficient space__. 

2. Run the script *prepare.sh*
    ```bash
    $ sbatch prepare.sh
    ```

**Note: The script assumes that the following tools are previously installed/on path: ```wget, bwa, gunzip```.**


## Run reads to counts pipeline

1. In the project dir store all RNA-seq fastq files inside the *raw* subdir
```bash
$ cd /path/to/project_dir
$ mkdir raw
$ mv /path/to/read_files/* raw
```

2. Download docker image based on your host architecture
```bash
$ singularity pull docker://ghcr.io/deeprob/glrna-{arch}:latest
```

3. Edit sbatch variables and environmental variables marked as "TODO" in *reads_to_counts.sh*

4. Run reads to counts pipeline
```bash
$ sbatch reads_to_counts.sh
```

## Run counts to de pipeline

1. Download docker image based on your host architecture (if not previously downloaded)
```bash
$ singularity pull docker://ghcr.io/deeprob/glrna-{arch}:latest
```

2. Create design matrix and store it inside project dir
The design matrix contains information about each sample. The sample names provided in the design matrix must be exactly same the name in the sample counts file. Control samples info should be given first followed by the treated samples. Example design matrix file is provided in the *examples* folder.

3. Edit sbatch variables and environmental variables marked as "TODO" in *counts_to_de.sh*

4. Run counts to de pipeline
```bash
$ sbatch counts_to_de.sh
```
