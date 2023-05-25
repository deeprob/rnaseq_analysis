# RNASeq analysis steps
Steps: | Trim | Align | Count | DE | QC | Enrich |
|---|---|---|---|---|---|---|
Tools: | Trimmomatic | STAR | HTSeq | DESeq2 | sklearn & seaborn | clusterprofiler & enrichplot |

# Quickstart
## Preparation steps
Step 1: Create raw dir and store fastq files here
```bash
$ cd /path/to/project/dir
$ mkdir raw
$ mv /path/to/read_files/* raw
```

Step 2: Create genome dir and store genome locations and gene annotations here
```bash
$ cd /path/to/project/dir
$ mkdir genome/hg38
$ cd genome/hg38
$ wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
$ gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
$ bwa index -a bwtsw GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
$ cd ../../
$ mkdir gene_annotations
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf.gz
$ gunzip gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf.gz
```

Step 3: Create adapters dir and store adapter sequence files here
```bash
$ cd /path/to/project/dir
$ mkdir adapters
$ 
```


## Running reads to counts
Step 1: Download docker image
```bash
$ docker pull ghcr.io/deeprob/glrnacounts:latest
```

Step 2: Run docker image 
```bash
$ docker run -v /path/to/data_dir:/data ghcr.io/deeprob/glrnacounts:latest --read_file1 hcc1395_normal_rep1_r1.fastq.gz --read_file2 hcc1395_normal_rep1_r2.fastq.gz --createstarindex --threads 4
```

## Running counts to de


Step 3: Download adapter sequence file for specific sequencers

We have provided the adapter sequences from Trimmomatic in the */examples/adapters/* folder. Either use those or you can add your own adapters in fasta format.

# Test pipeline

Step 1: Download and untar test data
```bash
$ wget http://genomedata.org/rnaseq-tutorial/practical.tar
$ tar -xvf practical.tar
$ rm pratical.tar
```

Step 2: 

# TODO:
1. Divide rnaseq counts and de completely (glrnacounts and glrnade)
2. No more conda env - use docker envs


# Steps

- Create metadata file
- Create conda environments
- Run 0_reads_to_counts.py
- Run 1_counts_to_de.py
- Run 2_analyze.py

**Most tools were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**
