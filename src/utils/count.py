import subprocess


def count(alignment_files, gtf_file, output_file, paired=False):
    # count using htseq
    with open(output_file, "w") as outfile:
        if not paired:
            command = ["htseq-count", "-f", "bam", "-s", "no", "-r", "pos", "-i", "gene_name"] + alignment_files + [gtf_file]
        else:
            command = ["htseq-count", "-f", "bam", "-s", "yes", "-r", "pos", "-i", "gene_name"] + alignment_files + [gtf_file]
        subprocess.run(command, stdout=outfile)
    return
