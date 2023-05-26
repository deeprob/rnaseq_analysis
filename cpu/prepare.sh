#!/bin/bash

proj_dir=/data6/deepro/rna_cache # TODO: set dir to project dir

# download and process genome
cd $proj_dir
mkdir -p genome
cd genome
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
mv GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz genome.fasta.gz
gunzip genome.fasta.gz
bwa index -a bwtsw genome.fasta

# download gene annotations
cd $proj_dir
mkdir gene_annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf.gz
mv gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf.gz annotation.gtf.gz
gunzip annotation.gtf.gz

# download adapter sequences
cd $proj_dir
mkdir adapters
wget https://raw.githubusercontent.com/deeprob/rnaseq_analysis/new-interface/examples/adapters/TruSeq3-PE.fa
mv TruSeq3-PE.fa adapter.fa
