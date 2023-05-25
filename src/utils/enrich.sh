#!/bin/bash
# set -ue # don't set ue because it interrupts conda activate

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# activate conda environment
conda activate gseanew


gene_file=$1
gseaout_file=$2
keggout_file=$3
gseafigout_file=$4
keggfigout_file=$5
tmp_dir=$6

Rscript /data6/deepro/computational_pipelines/rnaseq_analysis/src/utils/enrich.R $gene_file $gseaout_file $keggout_file $gseafigout_file $keggfigout_file $tmp_dir
