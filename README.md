# RNA-seq analysis pipeline - Girirajan Lab
This repository holds the source code and instructions to run Girirajan Lab's RNA-seq analysis pipeline. To run the pipeline please refer to the README files present within the folder CPU or SLURM depending on your platform. 

The description of the files and folders present in this repository is as follows:

- cpu: This folder holds instructions to run the pipeline on local machine
- slurm: This folder holds instructions to run the pipeline on SLURM cluster
- docker: This folder holds the dockerfile used to create a docker image of the pipeline

# Current pipeline steps
Steps: | Trim | Align | Count | DE | QC | Enrich |
|---|---|---|---|---|---|---|
Tools: | Trimmomatic | STAR | HTSeq | DESeq2 | sklearn & seaborn | clusterprofiler & enrichplot |

# Quickstart
Step 0: Select platform (cpu or slurm cluster) and architecture (amd64 or arm64)

Step 1: Create a project dir and store all RNA-seq fastq files inside the *raw* subdir

Step 2: Download docker image based on your host architecture

Step 3: Run preparation code to get all required preliminary data inside project dir

Step 4: Run reads to counts 

Step 5: Create design matrix and store it inside project dir

Step 6: Run counts to de 

# Test pipeline
Step 1: Download and untar test data
```bash
$ wget http://genomedata.org/rnaseq-tutorial/practical.tar
$ tar -xvf practical.tar
$ rm pratical.tar
```

Step 2: Follow [Quickstart](#quickstart)

# Using alternate preliminary data
By default our pipeline downloads the human reference genome version hg38, gencode annotations v42 and illumina adapter sequence for trimming version TruSeq3-PE. However they can be easily changed by editing the *prepare.sh* script as given below:

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

# Running counts to de only with your own count matrix file



**Most tools were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**
