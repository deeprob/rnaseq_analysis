







cache_dir="/data6/deepro/rna_cache" # TODO: set cache dir path
glrna_image="/data6/deepro/rna_cache/glrna-amd64_latest.sif" # TODO: set glrna pulled image path

singularity exec --containall -H $cache_dir:/data -B $cache_dir:/data $glrna_image bash /rnacounts/prepare.sh

