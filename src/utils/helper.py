import re
import os
import json
import subprocess
from argparse import Namespace

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))


####################
# meta file parser #
####################

def load_json(json_file):
    with open(json_file, "r") as f:
        jdict = json.load(f)
    return jdict

def create_args(meta_file, libname):
    meta_dict = load_json(meta_file)
    args = Namespace(
        # from metadata file
        prefix = meta_dict[libname]["prefix"],
        reps = meta_dict[libname]["reps"],
        pairs= meta_dict[libname]["readpairs"],
        suffix = meta_dict[libname]["fastqsuffix"],
        shortform = meta_dict[libname]["shortform"],
        reference_genome = meta_dict["genome"]["fasta"],
        reference_genome_gtf = meta_dict["genome"]["gtf"],
    )
    return args


###################
# file processing #
###################

def create_dirs(root_dir, lib_short):
    trim_dir = os.path.join(root_dir, "trim", lib_short)
    align_dir = os.path.join(root_dir, "align", lib_short)
    count_dir = os.path.join(root_dir, "count", lib_short)
    tmp_dir = os.path.join(root_dir, "tmp", lib_short)

    os.makedirs(trim_dir, exist_ok=True)
    os.makedirs(align_dir, exist_ok=True)
    os.makedirs(count_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    return trim_dir, align_dir, count_dir, tmp_dir

def parse_rep_and_pairs(lib_reps, lib_pairs):
    lib_reps = lib_reps.split()
    lib_pairs = lib_pairs.split()
    return lib_reps, lib_pairs  



# def get_file_base(file_path):
#     file_pattern = re.compile("(.+)\.(fastq|fq)\.gz")
#     m = re.match(file_pattern, os.path.basename(file_path))
#     return m.group(1)

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

# def infer_counts_columns_helper(filename):
#     filebasename = os.path.basename(filename)
#     geno, sex, biorep, techrep, _ = filebasename.split("_")
#     return geno, sex, biorep, techrep

# def infer_counts_columns(filenames):
#     column_list = []
#     for fn in filenames:
#         g, s, b, t = infer_counts_columns_helper(fn)
#         column_list.append("_".join([g,s,b,t]))
#     return column_list

# def generate_design_matrix(colnames, savefile):
#     # genotype, sex, biorep
#     matrix = [list(map(str, colname.split("_")[:3])) for colname in colnames]
#     with open(savefile, "w") as f:
#         f.write(",genotype,sex,biorep\n")
#         for col,line in zip(colnames, matrix):
#             f.write(",".join([col]+line))
#             f.write("\n")
#     return
