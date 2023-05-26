#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=glrna_prep 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data6/deepro/rna_cache # TODO: set dir to project dir
#SBATCH -o /data6/deepro/rna_cache/slurm/logs/out_prepare.log # TODO: set slurm output file
#SBATCH -e /data6/deepro/rna_cache/slurm/logs/err_prepare.log # TODO: set slurm input file
#SBATCH --exclude=durga,ramona,laila


echo `date` starting job on $HOSTNAME

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

echo `date` ending job on $HOSTNAME
