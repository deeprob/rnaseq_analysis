import os
import argparse
import utils.de as utd
import utils.helper as uth


def create_dirs(root_dir):
    de_dir = os.path.join(root_dir, "de")

    os.makedirs(de_dir, exist_ok=True)
    return de_dir


def main(
    root_dir, counts_matrix, counts_mat_cols, 
    design_formula, design_matrix_file, contrast, 
    de_results_file
    ):
    # create appropriate directories
    de_dir = create_dirs(root_dir)

    # de 
    count_out_file = counts_matrix
    if counts_mat_cols:
        counts_cols = counts_mat_cols
    print("Differential Expression...")
    design_matrix_file = os.path.join(de_dir, design_matrix_file)
    uth.generate_design_matrix(counts_cols, design_matrix_file)
    results_file = os.path.join(de_dir, de_results_file)
    utd.run_deseq2(
        count_out_file, counts_cols, design_matrix_file, design_formula, contrast, results_file
        )
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='RNASeq Analysis pipeline.')
    # TODO: argument name with arguments
    parser.add_argument("root", type=str, help="Provide the root directory where RNASeq output will be stored")
    parser.add_argument("-c", "--countsmat", type=str, help="count matrix filename",  default="counts.tsv")
    parser.add_argument("-l", "--countsmatcols", type=str, help="count matrix columns to consider",  nargs='+', default=[])
    parser.add_argument("-x", "--formula", type=str, help="design formula", default="~genotype+sex")
    parser.add_argument("-m", "--matrix", type=str, help="design matrix out filename", default="design_mat.csv")
    parser.add_argument("-t", "--contrast", type=str, help="contrast for deseq2", default="genotype,Hetz,WT")
    parser.add_argument("-o", "--deoutfile", type=str, help="diff exp outfilename", default="de_results.csv")


    args = parser.parse_args()
 
    main(
        root_dir=args.root, counts_matrix=args.countsmat, counts_mat_cols=args.countsmatcols, 
        design_formula=args.formula, design_matrix_file=args.matrix, contrast=args.contrast, 
        de_results_file=args.deoutfile
        )
