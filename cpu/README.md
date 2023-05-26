# Running the rnaseq pipeline on CPU using docker

## Run data preparation using the provided shell script

1. Edit the *project_dir* variable in *prepare.sh* to your directory of choice. This directory will host all files required to run reads to counts part of the pipeline. __Make sure it has sufficient space__. 

2. Run the script *prepare.sh*
    ```bash
    $ bash prepare.sh
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
$ docker pull ghcr.io/deeprob/glrna-{arch}:latest
```

3. Run reads to counts 
```bash
$ docker run -v /path/to/project_dir:/data ghcr.io/deeprob/glrna-{arch}:latest python3 reads_to_counts.py read1_filename.fastq.gz --read_file2 read2_filename.fastq.gz --threads 4 --createstarindex
```

## Run counts to de pipeline

1. Download docker image based on your host architecture (if not previously downloaded)
```bash
$ docker pull ghcr.io/deeprob/glrna-{arch}:latest
```

2. Create design matrix and store it inside project dir
The design matrix contains information about each sample. The sample names provided in the design matrix must be exactly same the name in the sample counts file. Control samples info should be given first followed by the treated samples. Example design matrix file is provided in the *examples* folder.

3. Run counts to de 
```bash
$ docker run -v /path/to/project_dir:/data ghcr.io/deeprob/glrna-{arch}:latest python3 counts_to_de.py --treatment_names treated_rep1_colname treated_rep2_colname treated_rep3_colname --control_names control_rep1_colname control_rep2_colname control_rep3_colname --design_matrix /data/relative/path/inside/project_dir/to/design_matrix.csv
```
