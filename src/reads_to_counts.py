import os
import argparse
import logging
import utils.helper as uth


def main(
    raw_dir,
    store_dir,
    index_dir,
    read1_file, 
    read2_file,
    genome_file,
    gtf_file,
    adapter_file,
    counts_matrix,
    index_flag,
    threads,
    ):

    # create appropriate directories
    trim_dir, align_dir, count_dir, tmp_dir = uth.create_dirs(store_dir)
    logfile = os.path.join(tmp_dir, "glrnacounts.log")
    logging.basicConfig(filename=logfile, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO, filemode="w")
    logging.info("Starting Girirajan Lab RNASeq reads to counts pipeline ...")

    # preprocess for align
    if index_flag:
        logging.info("Creating STAR index ...")
        uth.create_star_index(genome_file, gtf_file, index_dir)

    # trim
    logging.info("Trimming ...")
    out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2 = uth.trim_reads(raw_dir, read1_file, read2_file, trim_dir, adapter_file, threads)
    logging.info(f"Trimmed reads stored in {out_paired_1} & {out_paired_2}")

    # align
    logging.info("Aligning ...")
    aligned_file = uth.align(trim_dir, out_paired_1, out_paired_2, align_dir, index_dir, threads)
    logging.info(f"Aligned reads stored in {aligned_file}")

    # count
    logging.info("Counting ...")
    paired=False
    if read2_file:
        paired=True
    count_file = uth.count(align_dir, aligned_file, paired, gtf_file, count_dir, counts_matrix, threads)
    logging.info(f"Results stored in {count_file}.")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis read to counts pipeline.')
    parser.add_argument("--read_file1", type=str, help="Baseame of primary read file")
    parser.add_argument("--read_file2", type=str, help="Filename of paired end read file", default="")
    parser.add_argument("--rawdir", type=str, help="Provide the raw directory where RNASeq input fastq files are stored", default="/data/raw")
    parser.add_argument("--trim_adapter", type=str, help="Path of file with adapter sequences to be trimmed", default="/data/adapters/TruSeq3-PE.fa")
    parser.add_argument("--genome_file", type=str, help="Path to the reference genome file", default="/data/genome/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta")
    parser.add_argument("--gtf_file",  type=str, help="Path to the genome annotation file", default="/data/genome/gene_annotations/gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf")
    parser.add_argument("--storedir", type=str, help="Provide the project directory where RNASeq output will be stored", default="/data/")
    parser.add_argument("--indexdir", type=str, help="Provide the index directory where STARR index is stored or will be stored", default="/data/starrindex")
    parser.add_argument("--createstarindex", help="Create star index", action="store_true")
    parser.add_argument("--countsmat", type=str, help="count matrix filename",  default="counts.tsv")
    parser.add_argument("--threads", help="number of threads to use", default=64, type=int)

    cli_args = parser.parse_args()
 
    main(
        raw_dir=cli_args.rawdir,
        store_dir=cli_args.storedir, 
        index_dir=cli_args.indexdir,
        read1_file=cli_args.read_file1,
        read2_file=cli_args.read_file2,
        genome_file=cli_args.genome_file, 
        gtf_file=cli_args.gtf_file,
        adapter_file=cli_args.trim_adapter,
        counts_matrix=cli_args.countsmat,
        index_flag=cli_args.createstarindex,
        threads=cli_args.threads,
        )
