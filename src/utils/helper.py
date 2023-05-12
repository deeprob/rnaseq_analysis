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

def get_factor_dict_from_meta_file(meta_file, libname, factornames):
    meta_dict = load_json(meta_file)
    factor_dict = {fn: meta_dict[libname][fn].split() for fn in factornames}
    return factor_dict

###################
# file processing #
###################

def create_dirs(store_dir):
    trim_dir = os.path.join(store_dir, "trim")
    align_dir = os.path.join(store_dir, "align")
    count_dir = os.path.join(store_dir, "count")
    tmp_dir = os.path.join(store_dir, "tmp")

    os.makedirs(trim_dir, exist_ok=True)
    os.makedirs(align_dir, exist_ok=True)
    os.makedirs(count_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    return trim_dir, align_dir, count_dir, tmp_dir


############
# trimming #
############

# TODO: add adapter sequence path option here
def trim_reads_helper(read_pair_1, read_pair_2, out_pair_1, out_pair_2, out_unpaired_1, out_unpaired_2, nthreads=64):
    # trim reads using trimmomatic
    # TODO: take adapter sequence file as input
    path_to_adapter_sequences_se = "/data5/deepro/miniconda3/envs/rnaseq/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa"
    path_to_adapter_sequences_pe = "/data5/deepro/miniconda3/envs/rnaseq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

    if read_pair_2:
        # call paired end trimmomatic
        subprocess.call(["trimmomatic", "PE", "-threads", f"{nthreads}", "-phred33", 
                        f"{read_pair_1}", f"{read_pair_2}", f"{out_pair_1}", f"{out_unpaired_1}", f"{out_pair_2}", f"{out_unpaired_2}",   
                        f"ILLUMINACLIP:{path_to_adapter_sequences_pe}:2:30:10:7:true", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:51"])
    else:
        # call single end trimmomatic
        subprocess.call(["trimmomatic", "SE", "-threads", f"{nthreads}", "-phred33", 
                        f"{read_pair_1}", f"{out_pair_1}",    
                        f"ILLUMINACLIP:{path_to_adapter_sequences_se}:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:51"])

    return

def trim_reads(in_dir, read_pair_1, read_pair_2, suff, trim_dir, threads):
    trim_paired_suff = "fastq.gz"
    trim_unpaired_suff = "unpaired.fastq.gz"
    for rep in reps:
        read_pair_1 = os.path.join(in_dir, read_pair_1)
        
        out_paired_1 = os.path.join(trim_dir, "_".join([pre, rep, pairs[0]]) + f".{trim_paired_suff}")
        out_unpaired_1 = os.path.join(trim_dir, "_".join([pre, rep, pairs[0]]) + f".{trim_unpaired_suff}")
        read_pair_2 = ""
        out_paired_2 = ""
        out_unpaired_2 = ""
        if read_pair_2:
            read_pair_2 = os.path.join(in_dir, "_".join([pre, rep, pairs[1]]) + f".{suff}")
            out_paired_2 = os.path.join(trim_dir, "_".join([pre, rep, pairs[1]]) + f".{trim_paired_suff}")
            out_unpaired_2 = os.path.join(trim_dir, "_".join([pre, rep, pairs[1]]) + f".{trim_unpaired_suff}")
        trim_reads_helper(read_pair_1, read_pair_2, out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2, threads)
    return trim_paired_suff


#############
# alignment #
#############

def create_star_index(genome_file, gtf_file, storage_dir, nthreads=64):
    # create index for star aligner
    subprocess.call(["STAR", "--runMode", "genomeGenerate",  "--genomeDir", f"{storage_dir}",
                    "--genomeFastaFiles", f"{genome_file}", "--sjdbGTFfile", f"{gtf_file}", "--sjdbOverhang", "100", 
                    "--runThreadN", f"{nthreads}"])
    return   


def align_helper(read_files, output_prefix, starindex_dir, nthreads=64):
    # align using star
    read_files_list = read_files.split()
    command = ["STAR", "--genomeDir", f"{starindex_dir}", "--runThreadN", f"{nthreads}", "--outFileNamePrefix", f"{output_prefix}", "--readFilesIn"]
    command += read_files_list
    command += ["--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within", "--outSAMattributes", "Standard", "--readFilesCommand", "zcat"]
    subprocess.call(command)
    return


def align(in_dir, pre, reps, pairs, suff, align_dir, index_dir, threads):
    align_suff = "Aligned.sortedByCoord.out.bam"
    for rep in reps:
        read_pair_1 = os.path.join(in_dir, "_".join([pre, rep, pairs[0]]) + f".{suff}")
        read_pair_2 = ""
        read_files = read_pair_1
        if len(pairs)==2:
            read_pair_2 = os.path.join(in_dir, "_".join([pre, rep, pairs[1]]) + f".{suff}")
            read_files = " ".join([read_pair_1, read_pair_2])
        align_out_prefix = os.path.join(align_dir, "_".join([pre, rep+"."]))
        align_helper(read_files, align_out_prefix, index_dir, threads)
    return align_suff
