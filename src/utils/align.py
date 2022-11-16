import os
import subprocess


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
