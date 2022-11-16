import os
import argparse
import utils.de as utd


def main(
    store_dir, 
    counts_matrices, 
    design_formula_components,
    contrast,
    meta_counts_file,
    design_matrix_file, 
    de_results_file
    ):
    # create appropriate directories
    de_dir = utd.create_dirs(store_dir)

    # get formula
    design_formula = utd.get_formula(design_formula_components)

    # make meta counts matrix and get counts columns
    counts_columns, de_dir, meta_counts_outfile, two_conditions = utd.make_meta_counts(counts_matrices, de_dir, meta_counts_file)

    # make design matrix
    design_matrix_file = os.path.join(de_dir, design_matrix_file)
    utd.generate_design_matrix(counts_columns, design_formula_components, design_matrix_file)
    results_file = os.path.join(de_dir, de_results_file)

    # create contrast
    if not contrast:
        contrast = ",".join([design_formula_components[0], two_conditions[1], two_conditions[0]])
        print(f"Genes with logFC > 0 are overexpressed in {two_conditions[1]}.")
    
    # deseq2
    print("Differential Expression ...")
    utd.run_deseq2(
        meta_counts_outfile, counts_columns, design_matrix_file, design_formula, contrast, results_file
        )
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    parser.add_argument("store_dir", type=str, help="Directory where RNASeq output will be stored")
    parser.add_argument("-c", "--countsmat", type=str, help="count matrix filepaths for each conditions in formula", nargs='+')
    parser.add_argument("-x", "--formula", type=str, help="design formula components", nargs='+', default=["condition"])
    parser.add_argument("-t", "--contrast", type=str, help="contrast for deseq2, by default it is calculated based on which count matrix comes first", default="")
    parser.add_argument("-n", "--metacountsoutfile", type=str, help="meta counts matrix filename", default="meta_counts.csv")
    parser.add_argument("-m", "--matrixoutfile", type=str, help="design matrix out filename", default="design_mat.csv")
    parser.add_argument("-o", "--deoutfile", type=str, help="diff exp outfilename", default="de_results.csv")


    cli_args = parser.parse_args()
 
    main(
        store_dir=cli_args.store_dir, 
        counts_matrices=cli_args.countsmat, 
        design_formula_components=cli_args.formula,
        contrast=cli_args.contrast,
        meta_counts_file=cli_args.metacountsoutfile,
        design_matrix_file=cli_args.matrixoutfile, 
        de_results_file=cli_args.deoutfile
        )
