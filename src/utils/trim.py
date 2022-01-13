import subprocess

# TODO: add adapter sequence path option here
def trim_reads(read_pair_1, read_pair_2, out_pair_1, out_pair_2, out_unpaired_1, out_unpaired_2, nthreads=64):
    # trim reads using trimmomatic
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
