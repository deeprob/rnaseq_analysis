import subprocess


def count(alignment_files, counts_cols, gtf_file, output_file, paired=False, threads=64):
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
