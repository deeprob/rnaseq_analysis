import os
import argparse
import utils.de as utd
import utils.helper as uth


def main(
    counts_dir,
    store_dir,
    control_shortform,
    treatment_shortform,
    control_pre,
    treatment_pre,
    control_reps,
    treatment_reps,
    control_factor_dict,
    treatment_factor_dict,
    controlcountsmat,
    treatmentcountsmat,
    design_formula_components,
    meta_counts_file,
    design_matrix_file, 
    de_results_file,
    meta_norm_counts_file
    ):
    # create appropriate directories
    de_dir = utd.create_dirs(store_dir)

    # get formula
    design_formula = utd.get_formula(design_formula_components)

    # make meta counts matrix and get counts columns
    counts_columns, de_dir, meta_counts_outfile = utd.make_meta_counts(
        counts_dir, 
        control_shortform, treatment_shortform, 
        controlcountsmat, treatmentcountsmat,
        control_pre, treatment_pre,
        control_reps, treatment_reps,
        de_dir, meta_counts_file
        )

    # make design matrix
    design_matrix_file = os.path.join(de_dir, design_matrix_file)
    utd.generate_design_matrix(counts_columns, control_pre, treatment_pre, control_reps, treatment_reps, control_factor_dict, treatment_factor_dict, design_formula_components, design_matrix_file)

    # create contrast
    contrast = ",".join([design_formula_components[0], treatment_pre, control_pre])
    print(f"Genes with logFC > 0 are overexpressed in {treatment_pre}.")
    
    # create out filenames
    results_file = os.path.join(de_dir, de_results_file)
    meta_norm_counts_file = os.path.join(de_dir, meta_norm_counts_file)

    # deseq2
    print("Differential Expression ...")
    utd.run_deseq2(
        meta_counts_outfile, counts_columns, design_matrix_file, design_formula, contrast, results_file, meta_norm_counts_file
        )
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    parser.add_argument("meta_file", type=str, help="the meta file path with library information")
    parser.add_argument("counts_dir", type=str, help="Provide the raw directory where RNASeq counts files are stored")
    parser.add_argument("store_dir", type=str, help="Directory where RNASeq DE output will be stored")
    parser.add_argument("control", type=str, help="The control library name")
    parser.add_argument("treatment", type=str, help="The treatment library name")
    parser.add_argument("-c", "--controlcountsmat", type=str, help="control count matrix filename",  default="counts.tsv")
    parser.add_argument("-t", "--treatmentcountsmat", type=str, help="treatment count matrix filename",  default="counts.tsv")
    parser.add_argument("-x", "--formula", type=str, help="design formula components", nargs='+', default=["condition"])
    parser.add_argument("-n", "--metacountsoutfile", type=str, help="meta counts matrix filename", default="meta_counts.csv")
    parser.add_argument("-m", "--matrixoutfile", type=str, help="design matrix out filename", default="design_mat.csv")
    parser.add_argument("-o", "--deoutfile", type=str, help="diff exp outfilename", default="de_results.csv")
    parser.add_argument("-s", "--normctsoutfile", type=str, help="DESeq2 normalized counts out file", default="meta_norm_counts.csv")


    cli_args = parser.parse_args()
    control_args = uth.create_args(cli_args.meta_file, cli_args.control)
    treatment_args = uth.create_args(cli_args.meta_file, cli_args.treatment)

    control_factor_dict, treatment_factor_dict = dict(), dict()
    if len(cli_args.formula)>1:
        control_factor_dict = uth.get_factor_dict_from_meta_file(cli_args.meta_file, cli_args.control, cli_args.formula[1:])
        treatment_factor_dict = uth.get_factor_dict_from_meta_file(cli_args.meta_file, cli_args.treatment, cli_args.formula[1:])
 
    main(
        counts_dir=cli_args.counts_dir,
        store_dir=cli_args.store_dir,
        control_shortform=control_args.shortform,
        treatment_shortform=treatment_args.shortform,
        control_pre=control_args.prefix,
        treatment_pre=treatment_args.prefix,
        control_reps=control_args.reps,
        treatment_reps=treatment_args.reps,
        control_factor_dict=control_factor_dict,
        treatment_factor_dict=treatment_factor_dict,
        controlcountsmat=cli_args.controlcountsmat,
        treatmentcountsmat=cli_args.treatmentcountsmat,
        design_formula_components=cli_args.formula,
        meta_counts_file=cli_args.metacountsoutfile,
        design_matrix_file=cli_args.matrixoutfile, 
        de_results_file=cli_args.deoutfile,
        meta_norm_counts_file=cli_args.normctsoutfile
        )
