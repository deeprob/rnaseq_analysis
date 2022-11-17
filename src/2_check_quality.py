import os
import argparse
import utils.qc as utq
import utils.helper as uth

def main(
    counts_dir,
    libraries,
    controls,
    treatments,
    countsmats,
    deoutfiles,
    store_dir
    ):
    # create appropriate directories
    qc_dir = utq.create_dirs(store_dir)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    parser.add_argument("meta_file", type=str, help="the meta file path with library information")
    parser.add_argument("counts_dir", type=str, help="Provide the raw directory where RNASeq counts files are stored")
    parser.add_argument("store_dir", type=str, help="Directory where RNASeq DE output will be stored")
    parser.add_argument("-l", "--libraries", type=str, help="All library names for pca plot", nargs="+")
    parser.add_argument("-c", "--controls", type=str, help="The control library names for volcano plots", nargs="+")
    parser.add_argument("-t", "--treatments", type=str, help="The treatment library names for volvano plots", nargs="+")
    parser.add_argument("--countsmats", type=str, help="count matrix filebasenames", nargs="+", default=["counts.tsv"])
    parser.add_argument("--deoutfiles", type=str, help="diff exp outfilenames", nargs="+", default=["de_results.csv"])

    cli_args = parser.parse_args()
    lib_args = [uth.create_args(cli_args.meta_file, lib) for lib in cli_args.libraries]
    control_args = [uth.create_args(cli_args.meta_file, lib) for lib in cli_args.controls]
    treatment_args = [uth.create_args(cli_args.meta_file, lib) for lib in cli_args.treatments]
