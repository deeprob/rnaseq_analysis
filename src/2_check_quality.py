import os
import argparse
import utils.qc as utq
import utils.helper as uth

def main(
    counts_dir,
    de_dir,
    libraries,
    controls,
    treatments,
    countsmats,
    deoutfiles,
    geneidtonamefiles,
    store_dir
    ):
    # create appropriate directories
    qc_dir = utq.create_dirs(store_dir)

    # match lengths for libraries and countsmats | controls and deoutfiles | controls and geneidtonamefiles
    libraries, countsmats = utq.match_lengths(libraries, countsmats)
    controls, deoutfiles = utq.match_lengths(controls, deoutfiles)
    controls, geneidtonamefiles = utq.match_lengths(controls, geneidtonamefiles)

    # create pca plot
    utq.save_pca_plot(counts_dir, libraries, countsmats, qc_dir)

    # create volcano plots
    utq.create_volcano_plots(de_dir, treatments, controls, deoutfiles, geneidtonamefiles, qc_dir)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    parser.add_argument("meta_file", type=str, help="the meta file path with library information")
    parser.add_argument("counts_dir", type=str, help="Provide the input directory where RNASeq counts files are stored")
    parser.add_argument("de_dir", type=str, help="Provide the input directory where RNASeq de files are stored")
    parser.add_argument("store_dir", type=str, help="Directory where RNASeq DE output will be stored")
    parser.add_argument("-l", "--libraries", type=str, help="All library names for pca plot", nargs="+")
    parser.add_argument("-c", "--controls", type=str, help="The control library names for volcano plots", nargs="+")
    parser.add_argument("-t", "--treatments", type=str, help="The treatment library names for volvano plots", nargs="+")
    parser.add_argument("--countsmats", type=str, help="count matrix filebasenames", nargs="+", default=["counts.tsv"])
    parser.add_argument("--deoutfiles", type=str, help="diff exp outfilenames", nargs="+", default=["de_results.csv"])
    parser.add_argument("--geneidtonamefiles", type=str, help="diff exp gene id to name mappings", nargs="+", default=["geneid2name.csv"])

    cli_args = parser.parse_args()
    lib_shorts = [uth.create_args(cli_args.meta_file, lib).shortform for lib in cli_args.libraries]
    control_shorts = [uth.create_args(cli_args.meta_file, lib).shortform for lib in cli_args.controls]
    treatment_shorts = [uth.create_args(cli_args.meta_file, lib).shortform for lib in cli_args.treatments]

    main(
        cli_args.counts_dir,
        cli_args.de_dir,
        lib_shorts,
        control_shorts,
        treatment_shorts,
        cli_args.countsmats,
        cli_args.deoutfiles,
        cli_args.geneidtonamefiles,
        cli_args.store_dir,
    )
