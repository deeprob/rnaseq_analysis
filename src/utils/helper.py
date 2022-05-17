import re
import os
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def get_file_base(file_path):
    file_pattern = re.compile("(.+)\.(fastq|fq)\.gz")
    m = re.match(file_pattern, os.path.basename(file_path))
    return m.group(1)

# def get_normalized_counts(counts_file):
#     ncounts_file = os.path.splitext(counts_file)[0] + "_normalized.tsv"
#     command = ["Rscript", os.path.join(CURRENT_DIR, "normalize.R"), counts_file, ncounts_file]
#     subprocess.call(command)
#     return ncounts_file

# def get_log_cpm(ncounts_file):
#     logcpm_file = os.path.splitext(ncounts_file)[0] + "_lognormalized.tsv"
#     command = ["Rscript", os.path.join(CURRENT_DIR, "lognormalize.R"), ncounts_file, logcpm_file]
#     subprocess.call(command)
#     return logcpm_file

def infer_counts_columns_helper(filename):
    filebasename = os.path.basename(filename)
    geno, sex, biorep, techrep, _ = filebasename.split("_")
    return geno, sex, biorep, techrep

def infer_counts_columns(filenames):
    column_list = []
    for fn in filenames:
        g, s, b, t = infer_counts_columns_helper(fn)
        column_list.append("_".join([g,s,b,t]))
    return column_list

def generate_design_matrix(colnames, savefile):
    # genotype, sex, biorep
    matrix = [list(map(str, colname.split("_")[:3])) for colname in colnames]
    with open(savefile, "w") as f:
        f.write(",genotype,sex,biorep\n")
        for col,line in zip(colnames, matrix):
            f.write(",".join([col]+line))
            f.write("\n")
    return
