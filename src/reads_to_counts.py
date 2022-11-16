import os
import argparse
import utils.trim as utt
import utils.align as uta 
import utils.count as utc
import utils.helper as uth


def main(
    raw_dir,
    store_dir,
    index_dir,
    lib_prefix,
    lib_reps,
    lib_pairs,
    lib_suffix,
    lib_short,
    genome_file,
    gtf_file, 
    counts_matrix,
    index_flag,
    merge_flag,
    threads,
    ):

    # create appropriate directories
    trim_dir, align_dir, count_dir, tmp_dir = uth.create_dirs(store_dir, lib_short)

    # parse reps and pairs
    lib_reps, lib_pairs = uth.parse_rep_and_pairs(lib_reps, lib_pairs)

    # preprocess for align
    if index_flag:
        print("Creating STAR index ...")
        uta.create_star_index(genome_file, gtf_file, index_dir)

    # trim
    print("Trimming ...")
    trim_suff = utt.trim_reads(raw_dir, lib_prefix, lib_reps, lib_pairs, lib_suffix, trim_dir, threads)
        
    # align
    print("Aligning ...")
    align_suff = uta.align(trim_dir, lib_prefix, lib_reps, lib_pairs, trim_suff, align_dir, index_dir, threads)

    # merge technical replicates if merge flag is on 

    # count
    print("Counting ...")
    count_file = utc.count(align_dir, lib_prefix, lib_reps, lib_pairs, align_suff, gtf_file, count_dir, counts_matrix, threads)
    print(f"Results stored in {count_file}.")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    parser.add_argument("meta_file", type=str, help="the meta file path with library information")
    parser.add_argument("rawdir", type=str, help="Provide the raw directory where RNASeq input fastq files are stored")
    parser.add_argument("storedir", type=str, help="Provide the project directory where RNASeq output will be stored")
    parser.add_argument("indexdir", type=str, help="Provide the index directory where STARR index is stored or will be stored")
    parser.add_argument("library", type=str, help="the library name to extract from meta file")
    parser.add_argument("-c", "--countsmat", type=str, help="count matrix filename",  default="counts.tsv")
    parser.add_argument("-i", "--createstarindex", help="Create star index", action="store_true")
    parser.add_argument("-m", "--mergereps", help="Merge replicates", action="store_true")
    parser.add_argument("-p", "--threads", help="number of threads to use", default=64, type=int)

    # TODO: adapter sequence file for trimming
    # TODO: log file
    # TODO: merge technical replicates maybe?

    cli_args = parser.parse_args()
    lib_args = uth.create_args(cli_args.meta_file, cli_args.library)
 
    main(
        raw_dir=cli_args.rawdir,
        store_dir=cli_args.storedir, 
        index_dir=cli_args.indexdir,
        lib_prefix=lib_args.prefix, 
        lib_reps=lib_args.reps,
        lib_pairs=lib_args.pairs,
        lib_suffix=lib_args.suffix,
        lib_short=lib_args.shortform,
        genome_file=lib_args.reference_genome, 
        gtf_file=lib_args.reference_genome_gtf, 
        counts_matrix=cli_args.countsmat,
        index_flag=cli_args.createstarindex,
        merge_flag=cli_args.mergereps,
        threads=cli_args.threads,
        )
