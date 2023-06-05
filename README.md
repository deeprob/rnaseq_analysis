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

# Quickstart with example test data
0. Select platform (cpu or slurm cluster) and architecture (amd64 or arm64)

    Assuming cpu platform and arm64 architecture

1. Create a project dir and store all RNA-seq fastq files inside the *raw* subdir
    ```bash
    $ mkdir /path/to/project_dir
    $ cd  /path/to/project_dir
    $ mkdir raw
    $ cd raw
    $ wget http://genomedata.org/rnaseq-tutorial/practical.tar
    $ tar -xvf practical.tar
    $ rm practical.tar
    ```
    **Note**

    The following files will be inside raw dir after running the above commands: 
    ```
    hcc1395_normal_rep1_r1.fastq.gz hcc1395_normal_rep1_r2.fastq.gz
    hcc1395_normal_rep2_r1.fastq.gz hcc1395_normal_rep2_r2.fastq.gz
    hcc1395_normal_rep3_r1.fastq.gz hcc1395_normal_rep3_r2.fastq.gz
    hcc1395_tumor_rep1_r1.fastq.gz hcc1395_tumor_rep1_r2.fastq.gz
    hcc1395_tumor_rep2_r1.fastq.gz hcc1395_tumor_rep2_r2.fastq.gz
    hcc1395_tumor_rep3_r1.fastq.gz hcc1395_tumor_rep3_r2.fastq.gz
    ```

2. Run preparation code to get all required preliminary data inside project dir
    ```bash
    $ # change project dir path in prepare.sh
    $ bash prepare.sh
    ```

3. Download docker image based on your host architecture
    ```bash
    $ docker pull ghcr.io/deeprob/glrna-arm64:latest
    ```

4. Run reads to counts
    ```bash
    $ # mount project dir to the /data dir in container
    $ docker run -v /path/to/project_dir:/data ghcr.io/deeprob/glrna-arm64:latest python3 reads_to_counts.py hcc1395_normal_rep1_r1.fastq.gz --read_file2 hcc1395_normal_rep1_r2.fastq.gz --threads 4 --createstarindex
    ```
    **Note**

    Run the above command for all the above paired files without the createstarindex command since starindex needs to be created only once. It is both time and memory intensive. It requires at least 36 GB of RAM to run star indexing and alignment. 

5. Create design matrix and store it inside project dir

    The design matrix is a csv file that contains the read pair1 filename followed by "Aligned" as the first column and the type of condition of that sample as the second column. For the above example, it will look like this:
    ```bash
    ,condition
    hcc1395_normal_rep1_r1Aligned,control
    hcc1395_normal_rep1_r1Aligned,control
    hcc1395_normal_rep1_r1Aligned,control
    hcc1395_tumor_rep1_r1Aligned,treatment
    hcc1395_tumor_rep1_r1Aligned,treatment
    hcc1395_tumor_rep1_r1Aligned,treatment
    ```
    **Note**

    Store the design matrix inside the project dir. It is recommended to create a de dir inside project dir and store it with the name design_matrix.csv inside the de dir. Then the pipeline will detect it automatically. However, it can be stored anywhere inside the project dir but the --design_matrix argument should then be called to specify the path of the file within the container, **not** the host system.


6. Run counts to de
    ```bash
    $ docker run -v /path/to/project_dir:/data ghcr.io/deeprob/glrna-arm64:latest python3 counts_to_de.py --treatment_names hcc1395_tumor_rep1_r1Aligned hcc1395_tumor_rep1_r1Aligned hcc1395_tumor_rep1_r1Aligned --control_names hcc1395_normal_rep1_r1Aligned hcc1395_normal_rep1_r1Aligned hcc1395_normal_rep1_r1Aligned --design_matrix /data/relative/path/inside/project_dir/to/design_matrix.csv
    ```



# Using alternate preliminary data
By default our pipeline downloads the human reference genome version hg38, gencode annotations v42 and illumina adapter sequence for trimming version TruSeq3-PE. However they can be easily changed by editing the *prepare.sh* script as given below:

1. Alternate genome file: Let's assume that you would like to use mm39 mouse genome.
    ```bash
    $ # inside project dir create genome folder with version number
    $ cd /path/to/project_dir
    $ mkdir -p genome
    $ cd genome
    $ # download genome fasta file from your link
    $ wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
    $ mv mm39.fa.gz genome.fa.gz
    $ gunzip genome.fa.gz
    $ # index genome fasta file using bwa (make sure it is on PATH)
    $ bwa index -a bwtsw genome.fa
    ```

2. Alternate gtf file: Let's assume that you would like to use mm39 mouse genome.
    ```bash
    $ # inside project dir create genome folder with gene annotations subdir
    $ cd /path/to/project_dir
    $ mkdir -p gene_annotations
    $ cd gene_annotations
    $ # download gene annotations file from your link
    $ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf.gz
    $ mv gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf.gz annotation.gtf.gz
    $ gunzip annotation.gtf.gz
    ```

3. Alternate adapter sequence: Lets assume you would like to use TruSeq3 single end reads
    ```bash
    $ # inside project dir create adapters folder
    $ cd /path/to/project_dir
    $ mkdir -p adapters
    $ cd adapters
    $ # download adapters file from your link
    $ wget https://raw.githubusercontent.com/deeprob/rnaseq_analysis/new-interface/examples/adapters/TruSeq3-SE.fa
    $ mv TruSeq3-SE.fa adapter.fa
    ```

# Running counts to de only with your own count matrix file
*In progress*


**Most tools were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**
