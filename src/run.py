import os
import argparse
import utils.trim as utt
import utils.align as uta 
import utils.count as utc
import utils.helper as uth
from itertools import zip_longest


def create_dirs(root_dir):
    trim_dir = os.path.join(root_dir, "trim")
    align_dir = os.path.join(root_dir, "align")
    count_dir = os.path.join(root_dir, "count")
    de_dir = os.path.join(root_dir, "de")
    tmp_dir = os.path.join(root_dir, "tmp")

    os.makedirs(trim_dir, exist_ok=True)
    os.makedirs(align_dir, exist_ok=True)
    os.makedirs(count_dir, exist_ok=True)
    os.makedirs(de_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    return trim_dir, align_dir, count_dir, de_dir, tmp_dir


def main(root_dir, read_files_1, read_files_2, genome_file, gtf_file, debug=True):
    # create appropriate directories
    trim_dir, align_dir, count_dir, de_dir, tmp_dir = create_dirs(root_dir)

    if read_files_2:
        assert len(read_files_1) ==  len(read_files_2)
        paired=True

    # preprocess for align
    if not debug:
        uta.create_star_index(genome_file, gtf_file, tmp_dir)

    # get all read_columns
    counts_cols = uth.infer_counts_columns(read_files_1)

    align_out_files = []

    for read_file1, read_file2 in zip_longest(read_files_1, read_files_2):
        # base name of read files; all files including paired files will be named based on the prefix read file prefix
        file1_base = uth.get_file_base(read_file1)
        if read_file2:
            file2_base = uth.get_file_base(read_file2)
        # trim
        os.makedirs(trim_dir, exist_ok=True)
        read_file1_paired_out = os.path.join(trim_dir, f"{file1_base}.fastq.gz")
        read_file1_unpaired_out = os.path.join(trim_dir, f"{file1_base}.unpaired.fastq.gz")
        if read_file2:
            read_file2_paired_out = os.path.join(trim_dir, f"{file2_base}.fastq.gz")
            read_file2_unpaired_out = os.path.join(trim_dir, f"{file2_base}.unpaired.fastq.gz")
        else:
            read_file2_paired_out = ""
            read_file2_unpaired_out = ""
        if not debug:
            print("Trimming...")
            utt.trim_reads(read_file1, read_file2, read_file1_paired_out, read_file2_paired_out, read_file1_unpaired_out, read_file2_unpaired_out)
        # align
        if read_file2_paired_out:
            read_files = " ".join([read_file1_paired_out, read_file2_paired_out])
        else:
            read_files = read_file1_paired_out
        
        # aligned files will be stored under the name of read1
        align_out_prefix = os.path.join(align_dir, file1_base)
        if not debug:
            print("Aligning...")
            uta.align(read_files, align_out_prefix, tmp_dir)

        align_out_file = align_out_prefix + "Aligned.sortedByCoord.out.bam"
        align_out_files.append(align_out_file)
    
    assert len(counts_cols) == len(align_out_files)
    print(counts_cols, align_out_files)

    # count
    print("Counting...")
    os.makedirs(count_dir, exist_ok=True)
    count_out_file = os.path.join(count_dir, "counts.tsv")
    if not debug:
        utc.count(align_out_files, gtf_file, count_out_file, paired=paired)

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    # TODO: argument name with arguments
    parser.add_argument("root", type=str, help="Provide the root directory where RNASeq output will be stored")
    parser.add_argument("-f", "--fwd", type=str, help="RNASeq forward reads", nargs='+')
    parser.add_argument("-r", "--rev", type=str, help="RNASeq reverse read", nargs='+', default=[])
    parser.add_argument("-g", "--genome", type=str, help="The genome fasta file")
    parser.add_argument("-e", "--gtf", type=str, help="GTF file with ensemble reference genes")
    parser.add_argument("-d", "--debug", action="store_true")
    # TODO: adapter sequence file for trimming
    # TODO: number of threads to parallelize tasks

    args = parser.parse_args()
 
    main(args.root, args.fwd, args.rev, args.genome, args.gtf, args.debug)

    # example command: python run.py /data5/deepro/starrseq/rnaseq all /data5/deepro/starrseq/rnaseq/raw/all/16P12_1.fastq.gz -g /data5/deepro/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -e /data5/deepro/genomes/gencode.v38.annotation.gtf
    # example command: python run.py /data5/deepro/starrseq/rnaseq rep3 /data5/deepro/starrseq/rnaseq/raw/rep3/CC_R1_S22_R1_001.fastq.gz /data5/deepro/starrseq/rnaseq/raw/rep3/CC_R2_S23_R1_001.fastq.gz /data5/deepro/starrseq/rnaseq/raw/rep3/CC_R3_S24_R1_001.fastq.gz -d
    # example command: python run.py /data5/deepro/mao_lab/embryonic_brains_gyf2_ko/results