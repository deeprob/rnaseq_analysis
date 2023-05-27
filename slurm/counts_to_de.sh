#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=glrna_cd 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data6/deepro/rna_cache # TODO: set dir to project dir
#SBATCH -o /data6/deepro/rna_cache/slurm/logs/out_cd.log # TODO: set slurm output file
#SBATCH -e /data6/deepro/rna_cache/slurm/logs/err_cd.log # TODO: set slurm input file
#SBATCH --nodelist=sarah



echo `date` starting job on $HOSTNAME

cache_dir="/data6/deepro/rna_cache" # TODO: set project dir path
glrna_image="/data6/deepro/rna_cache/glrna-amd64_latest.sif" # TODO: set glrna pulled image path

# the workdir command does not work in singularity, hence absolute path to script is required
# TODO: change read file1 and read file 2 name.
# TODO: remove --createstarrindex flag if starrindex is already created and in the project dir
singularity exec --containall -H $cache_dir:/data -B $cache_dir:/data $glrna_image python3 /rnacounts/counts_to_de.py --treatment_names hcc1395_tumor_rep1_r1Aligned hcc1395_tumor_rep2_r1Aligned hcc1395_tumor_rep3_r1Aligned --control_names hcc1395_normal_rep1_r1Aligned hcc1395_normal_rep2_r1Aligned hcc1395_normal_rep3_r1Aligned 

echo `date` ending job on $HOSTNAME
