import os
import subprocess


def count_helper(alignment_files, counts_cols, gtf_file, output_file, paired=False, threads=64):
    # count using htseq
    with open(output_file, "w") as outfile:
        outfile.write("\t".join(["gene_name"] + counts_cols) + "\n")

    with open(output_file, "a") as outfile:
        if not paired:
            command = ["htseq-count", "-f", "bam", "-s", "no", "-r", "pos", "-i", "gene_name", "-n", f"{threads}"] + alignment_files + [gtf_file]
        else:
            command = ["htseq-count", "-f", "bam", "-s", "yes", "-r", "pos", "-i", "gene_name", "-n", f"{threads}"] + alignment_files + [gtf_file]
        subprocess.run(command, stdout=outfile)
    return

def count(in_dir, pre, reps, pairs, suff, gtf_file, count_dir, counts_matrix, threads):
    paired = False
    if len(pairs) == 2:
        paired = True
    aligned_files = [os.path.join(in_dir, "_".join([pre, rep]) + f".{suff}") for rep in reps]
    counts_cols = ["_".join([pre, rep]) for rep in reps]
    count_out_file = os.path.join(count_dir, counts_matrix)
    count_helper(aligned_files, counts_cols, gtf_file, count_out_file, paired=paired, threads=threads)
    return count_out_file
