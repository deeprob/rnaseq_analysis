import re
import os
import json
import subprocess
import pandas as pd

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
deseq2_script_path = os.path.join(os.path.dirname(__file__), "deseq2.R")

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

def remove_ext(filename):
    find = re.compile(r"^(.*?)\..*")
    basename = re.match(find, filename).group(1)
    return basename

############
# trimming #
############

def trim_reads_helper(read_pair_1, read_pair_2, out_pair_1, out_pair_2, out_unpaired_1, out_unpaired_2, adapter_file, nthreads=64):
    # trim reads using trimmomatic
    if read_pair_2:
        # call paired end trimmomatic
        subprocess.call(["TrimmomaticPE", "-threads", f"{nthreads}", "-phred33", 
                        f"{read_pair_1}", f"{read_pair_2}", f"{out_pair_1}", f"{out_unpaired_1}", f"{out_pair_2}", f"{out_unpaired_2}",   
                        f"ILLUMINACLIP:{adapter_file}:2:30:10:7:true", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:51"])
    else:
        # call single end trimmomatic
        subprocess.call(["TrimmomaticSE", "-threads", f"{nthreads}", "-phred33", 
                        f"{read_pair_1}", f"{out_pair_1}",    
                        f"ILLUMINACLIP:{adapter_file}:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:51"])

    return


def trim_reads(in_dir, read_pair_1, read_pair_2, trim_dir, adapter_file, threads):
    trim_paired_suff = "fastq.gz"
    trim_unpaired_suff = "unpaired.fastq.gz"
    read_pair_1 = os.path.join(in_dir, read_pair_1)
    out_paired_1 = remove_ext(os.path.basename(read_pair_1))
    out_paired_1 = os.path.join(trim_dir, out_paired_1 + f".{trim_paired_suff}")
    out_unpaired_1 = os.path.join(trim_dir, out_paired_1 + f".{trim_unpaired_suff}")
    out_paired_2 = ""
    out_unpaired_2 = ""
    if read_pair_2:
        read_pair_2 = os.path.join(in_dir, read_pair_2)
        out_paired_2 = remove_ext(os.path.basename(read_pair_2))
        out_paired_2 = os.path.join(trim_dir, out_paired_2 + f".{trim_paired_suff}")
        out_unpaired_2 = os.path.join(trim_dir, out_paired_2 + f".{trim_unpaired_suff}")
    trim_reads_helper(read_pair_1, read_pair_2, out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2, adapter_file, threads)
    return out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2

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


def align(in_dir, read_pair_1, read_pair_2, align_dir, index_dir, threads):
    align_suff = "Aligned.sortedByCoord.out.bam"
    read_files = read_pair_1
    if read_pair_2:
        read_files = " ".join([read_pair_1, read_pair_2])
    align_out_prefix = remove_ext(os.path.basename(read_pair_1))
    align_out_prefix = os.path.join(align_dir, align_out_prefix)
    align_helper(read_files, align_out_prefix, index_dir, threads)
    aligned_file = align_out_prefix + align_suff
    return aligned_file

############
# counting #
############

def count_helper(alignment_files, counts_cols, gtf_file, output_file, paired=False, threads=64):
    # count using htseq
    with open(output_file, "w") as outfile:
        outfile.write("\t".join(["gene_id", "gene_name"] + counts_cols) + "\n")

    with open(output_file, "a") as outfile:
        if not paired:
            command = ["htseq-count", "-f", "bam", "-s", "no", "-r", "pos", "-i", "gene_id", "--additional-attr", "gene_name", "-n", f"{threads}"] + alignment_files + [gtf_file]
        else:
            command = ["htseq-count", "-f", "bam", "-s", "yes", "-r", "pos", "-i", "gene_id", "--additional-attr", "gene_name", "-n", f"{threads}"] + alignment_files + [gtf_file]
        subprocess.run(command, stdout=outfile)
    return

def count(in_dir, aligned_file, paired, gtf_file, count_dir, threads):
    aligned_files = [aligned_file]
    counts_cols = [remove_ext(os.path.basename(aligned_file))]
    count_out_file = os.path.join(count_dir, f"{counts_cols[0]}.tsv")
    count_helper(aligned_files, counts_cols, gtf_file, count_out_file, paired=paired, threads=threads)
    return count_out_file

####################
# meta counts file #
####################

def read_counts_matrix(counts_dir, counts_col):
    df = pd.read_csv(os.path.join(counts_dir, f"{counts_col}.tsv"), sep="\t", index_col=["gene_id", "gene_name"])
    return df

def make_meta_counts_mat_helper(counts_dir, counts_cols):
    counts_dfs = [read_counts_matrix(counts_dir, cc) for cc in counts_cols]
    meta_count_df = pd.concat(counts_dfs, axis=1)
    meta_count_df.columns = [f"X{c}" if c[0].isdigit() else c for c in meta_count_df.columns]
    return meta_count_df

def make_meta_counts(
        counts_dir, 
        counts_cols,
        meta_counts_outfile
        ):
    meta_count_df = make_meta_counts_mat_helper(counts_dir, counts_cols)
    # save meta counts file
    meta_count_df.iloc[:-5].to_csv(meta_counts_outfile, index=True, header=True)
    return

#################
# design matrix #
#################

def get_formula(design_matrix_file):
    df = pd.read_csv(design_matrix_file, nrows=1, index_col=0)
    return list(df.columns), f"~{'+'.join(list(df.columns))}"

#############
# deseq run #
#############

def run_deseq2(counts_file, counts_cols, designfile, design_formula, contrast, de_file, normcts_file):
	counts_cols_to_str = ",".join(counts_cols)
	cmd = ["Rscript", deseq2_script_path, counts_file, counts_cols_to_str, designfile, design_formula, contrast, de_file, normcts_file]
	subprocess.run(cmd)
	return
