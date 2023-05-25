import os
import argparse
import utils.qc as utq
import utils.helper as uth
import utils.enrich as ute

def main(
    de_dir,
    control_short,
    treatment_short,
    normcountsmat,
    deoutfile,
    geneidtonamefile,
    store_dir,
    genenames,
    designmatrixfile,
    ):
    # create appropriate directories
    qc_dir, enrich_dir = utq.create_dirs(store_dir)

    # create volcano plot
    utq.create_volcano_plot(de_dir, treatment_short, control_short, deoutfile, geneidtonamefile, qc_dir, genenames)

    # create library specific pca plot with normalized read counts
    utq.create_lib_specific_pca_plot(de_dir, treatment_short, control_short, normcountsmat, designmatrixfile, qc_dir)

    # get enriched go terms
    degenes_dir = ute.get_de_genes(de_dir, treatment_short, control_short, deoutfile, geneidtonamefile, enrich_dir, thresh=0.01)
    ute.run_enrichment(degenes_dir, treatment_short, control_short, enrich_dir)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    parser.add_argument("meta_file", type=str, help="the meta file path with library information")
    parser.add_argument("de_dir", type=str, help="Provide the input directory where RNASeq de files are stored")
    parser.add_argument("store_dir", type=str, help="Directory where RNASeq DE output will be stored")
    parser.add_argument("control", type=str, help="The control library name to generate analysis results")
    parser.add_argument("treatment", type=str, help="The treatment library names to generate analysis results")
    parser.add_argument("--deoutfile", type=str, help="diff exp outfilename", default="de_results.csv")
    parser.add_argument("--normcountsmat", type=str, help="normalized count matrix filebasename", default="meta_norm_counts.csv")
    parser.add_argument("--geneidtonamefile", type=str, help="diff exp gene id to name mappings", default="geneid2name.csv")
    parser.add_argument("--designmatrixfile", help="design matrix filename required to plot factors in pca plot", type=str, default="design_mat.csv")
    parser.add_argument("--genenames", type=str, help="Gene names of interest that will be plotted in the volcano", nargs="+", default=[])

    cli_args = parser.parse_args()
    control_args = uth.create_args(cli_args.meta_file, cli_args.control)
    treatment_args = uth.create_args(cli_args.meta_file, cli_args.treatment)

    main(
        cli_args.de_dir,
        control_args.shortform,
        treatment_args.shortform,
        cli_args.normcountsmat,
        cli_args.deoutfile,
        cli_args.geneidtonamefile,
        cli_args.store_dir,
        cli_args.genenames,
        cli_args.designmatrixfile,
    )
