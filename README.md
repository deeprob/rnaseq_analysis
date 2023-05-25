# RNASeq analysis steps
Steps: | Trim | Align | Count | DE | QC | Enrich |
|---|---|---|---|---|---|---|
Tools: | Trimmomatic | STAR | HTSeq | DESeq2 | sklearn & seaborn | clusterprofiler & enrichplot |

# Quickstart

## Running the pipeline
Step 1: Create a project dir and store all RNA-seq fastq files inside the *raw* subdir
```bash
$ cd /path/to/project_dir
$ mkdir raw
$ mv /path/to/read_files/* raw
```

Step 2: Download docker image
```bash
$ docker pull ghcr.io/deeprob/glrna:latest
```

Step 3: Run preparation code to get all required preliminary data
```bash
$ docker run -v /path/to/project_dir:/data ghcr.io/deeprob/glrna:latest prepare.sh
```
By default - reference genome version hg38, gencode annotation v42 and adapter sequence TruSeq3-PE.fa is downloaded. To change genome, gencode annotation version and adapter sequence look at Section: [Using alternate preliminary data](#using-alternate-preliminary-data)

Step 4: Run reads to counts 
```bash
$ docker run -v /path/to/project_dir:/data ghcr.io/deeprob/glrna:latest python3 reads_to_counts.py read1_filename.fastq.gz --read_file2 read2_filename.fastq.gz --threads 4 --createstarindex
```

Step 5: Create design matrix and store it inside project dir
The design matrix contains information about each sample. The sample names provided in the design matrix must be exactly same the name in the sample counts file. Control samples info should be given first followed by the treated samples. Example design matrix file is provided in the *examples* folder.


Step 6: Run counts to de 
```bash
$ docker run -v /path/to/project_dir:/data ghcr.io/deeprob/glrna:latest python3 counts_to_de.py --treatment_names treated_rep1_colname treated_rep2_colname treated_rep3_colname --control_names control_rep1_colname control_rep2_colname control_rep3_colname --design_matrix /data/relative/path/inside/project_dir/to/design_matrix.csv
```

# Test pipeline
Step 1: Download and untar test data
```bash
$ wget http://genomedata.org/rnaseq-tutorial/practical.tar
$ tar -xvf practical.tar
$ rm pratical.tar
```

Step 2: Follow [Quickstart](#quickstart)

# Using alternate preliminary data
By default our pipeline is set to use the human reference genome version hg38, gencode annotations v42 and illumina adapter sequence for trimming version TruSeq3-PE. However they can be easily changed by:

1. Downloading all preliminary data inside the project dir albeit maintaining the dir structure (as below) 
2. Specifying the changed file location while running the docker image

## Downloading alternate data while maintaining dir structure
1. Alternate genome file: Let's assume that you would like to use mm39 mouse genome.
    ```bash
    $ # inside project dir create genome folder with version number
    $ cd /path/to/project_dir
    $ mkdir genome/mm39
    $ cd genome/mm39
    $ # download genome fasta file from your link
    $ wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
    $ gunzip mm39.fa.gz
    $ # index genome fasta file using bwa (make sure it is on PATH)
    $ bwa index -a bwtsw mm39.fa
    ```

2. Alternate gtf file: Let's assume that you would like to use mm39 mouse genome.
    ```bash
    $ # inside project dir create genome folder with gene annotations subdir
    $ cd /path/to/project_dir
    $ mkdir genome/gene_annotations
    $ cd genome/gene_annotations
    $ # download gene annotations file from your link
    $ wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz
    $ gunzip refGene.gtf.gz
    ```

3. Alternate adapter sequence: Lets assume you would like to use TruSeq3 single end reads
    ```bash
    $ # inside project dir create adapters folder
    $ cd /path/to/project_dir
    $ mkdir adapters
    $ # download adapters file from your link
    $ wget https://raw.githubusercontent.com/deeprob/rnaseq_analysis/new-interface/examples/adapters/TruSeq3-SE.fa
    ```

# TODO:
1. Divide rnaseq counts and de completely (glrnacounts and glrnade)
2. No more conda env - use docker envs


**Most tools were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**
