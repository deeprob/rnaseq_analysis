#!/bin/bash

proj_dir=/data

# download and process genome
cd proj_dir
mkdir genome/hg38
cd genome/hg38
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
bwa index -a bwtsw GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

# download gene annotations
cd proj_dir
mkdir gene_annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf.gz
gunzip gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf.gz

# download adapter sequences
cd proj_dir
mkdir adapters
wget https://raw.githubusercontent.com/deeprob/rnaseq_analysis/new-interface/examples/adapters/TruSeq3-PE.fa
