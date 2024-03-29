import os
import argparse
import logging
import utils.helper as uth


def main(
    treatment,
    control,
    treatment_cols,
    control_cols,
    counts_dir,
    store_dir,
    meta_counts_file,
    design_matrix_file, 
    de_results_file,
    meta_norm_counts_file,
    mingenecounts,
    ):
    # create appropriate directories
    de_dir = os.path.join(store_dir, f"{treatment}vs{control}")
    os.makedirs(de_dir, exist_ok=True)

    logfile = os.path.join(de_dir, "glrnade.log")
    logging.basicConfig(filename=logfile, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO, filemode="w")
    logging.info("Starting Girirajan Lab RNASeq counts to differential expression pipeline ...")
    logging.info(f"Created dir {de_dir} ...")

    # get count columns
    counts_columns = control_cols + treatment_cols

    # make meta counts matrix if not given
    if meta_counts_file=="":
        logging.info(f"Creating meta counts file ...")
        meta_counts_file = os.path.join(de_dir, "meta_counts.csv")
        count_col_dict = uth.make_meta_counts(counts_dir, counts_columns, meta_counts_file)
    else:
        logging.info(f"Using user provided meta counts file ...")
        count_col_dict = uth.rename_counts_columns_to_convention(counts_columns)

    counts_columns = [count_col_dict[c] for c in counts_columns]

    # get design formula components
    design_formula_components, design_formula, design_matrix_file = uth.get_formula(design_matrix_file, count_col_dict, de_dir)
    logging.info(f"The design formula is {design_formula}.")

    # create contrast
    contrast = ",".join([design_formula_components[0], treatment, control])
    logging.info(f"Genes with logFC > 0 are overexpressed in {treatment}.")
    
    # create out filenames
    results_file = os.path.join(de_dir, de_results_file)
    meta_norm_counts_file = os.path.join(de_dir, meta_norm_counts_file)

    # deseq2
    logging.info("Running Differential Expression using deseq2 ...")
    uth.run_deseq2(
        meta_counts_file, counts_columns, design_matrix_file, design_formula, contrast, results_file, meta_norm_counts_file, mingenecounts
        )
    
    # create pca and volcano plots
    uth.create_lib_specific_pca_plot(de_dir, meta_norm_counts_file, design_matrix_file)
    uth.create_volcano_plot(de_dir, results_file, meta_counts_file)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    parser.add_argument("--treatment_names", type=str, help="The treatment names as given in the count matrix", nargs="+")
    parser.add_argument("--control_names", type=str, help="The control names as given in the count matrix", nargs="+")
    parser.add_argument("--treatment", type=str, help="The treatment library name", default="treatment")
    parser.add_argument("--control", type=str, help="The control library name", default="control")
    parser.add_argument("--counts_dir", type=str, help="Provide the raw directory where RNASeq counts files are stored", default="/data/count")
    parser.add_argument("--store_dir", type=str, help="Directory where RNASeq DE output will be stored", default="/data/de")
    parser.add_argument("--design_matrix", type=str, help="design matrix file", default="/data/de/design_matrix.csv")
    parser.add_argument("--metacountsfile", type=str, help="meta counts matrix filename", default="")
    parser.add_argument("--deoutfile", type=str, help="diff exp outfilename", default="de_results.csv")
    parser.add_argument("--normctsoutfile", type=str, help="DESeq2 normalized counts out file", default="meta_norm_counts.csv")
    parser.add_argument("--mingenecounts", type=int, help="Minimum count of genes in all samples combined", default=10)


    cli_args = parser.parse_args()
 
    main(
        treatment=cli_args.treatment,
        control=cli_args.control,
        treatment_cols=cli_args.treatment_names,
        control_cols=cli_args.control_names,
        counts_dir=cli_args.counts_dir,
        store_dir=cli_args.store_dir,
        meta_counts_file=cli_args.metacountsfile,
        design_matrix_file=cli_args.design_matrix, 
        de_results_file=cli_args.deoutfile,
        meta_norm_counts_file=cli_args.normctsoutfile,
        mingenecounts=cli_args.mingenecounts,
        )
