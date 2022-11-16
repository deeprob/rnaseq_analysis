import os
import subprocess

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

def trim_reads(in_dir, pre, reps, pairs, suff, trim_dir, threads):
    trim_paired_suff = "fastq.gz"
    trim_unpaired_suff = "unpaired.fastq.gz"
    for rep in reps:
        read_pair_1 = os.path.join(in_dir, "_".join([pre, rep, pairs[0]]) + f".{suff}")
        out_paired_1 = os.path.join(trim_dir, "_".join([pre, rep, pairs[0]]) + f".{trim_paired_suff}")
        out_unpaired_1 = os.path.join(trim_dir, "_".join([pre, rep, pairs[0]]) + f".{trim_unpaired_suff}")
        read_pair_2 = ""
        out_paired_2 = ""
        out_unpaired_2 = ""
        if len(pairs)==2:
            read_pair_2 = os.path.join(in_dir, "_".join([pre, rep, pairs[1]]) + f".{suff}")
            out_paired_2 = os.path.join(trim_dir, "_".join([pre, rep, pairs[1]]) + f".{trim_paired_suff}")
            out_unpaired_2 = os.path.join(trim_dir, "_".join([pre, rep, pairs[1]]) + f".{trim_unpaired_suff}")
        trim_reads_helper(read_pair_1, read_pair_2, out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2, threads)
    return trim_paired_suff
