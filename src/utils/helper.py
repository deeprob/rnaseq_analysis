import re
import os
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def get_file_base(file_path):
    file_pattern = re.compile("(.+)\.(fastq|fq)\.gz")
    m = re.match(file_pattern, os.path.basename(file_path))
    return m.group(1)

def get_normalized_counts(counts_file):
    ncounts_file = os.path.splitext(counts_file)[0] + "_normalized.tsv"
    command = ["Rscript", os.path.join(CURRENT_DIR, "normalize.R"), counts_file, ncounts_file]
    subprocess.call(command)
    return ncounts_file

def get_log_cpm(ncounts_file):
    logcpm_file = os.path.splitext(ncounts_file)[0] + "_lognormalized.tsv"
    command = ["Rscript", os.path.join(CURRENT_DIR, "lognormalize.R"), ncounts_file, logcpm_file]
    subprocess.call(command)
    return logcpm_file

# TODO: visualization counts with PCA function
def get_pca_reps():
    return
