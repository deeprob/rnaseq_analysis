import os
import subprocess
import pandas as pd

deseq2_script_path = os.path.join(os.path.dirname(__file__), "deseq2.R")


def create_dirs(store_dir):
    de_dir = os.path.join(store_dir, "de")
    os.makedirs(de_dir, exist_ok=True)
    return de_dir

def get_formula(components):
     return f"~{'+'.join(components)}"

def read_counts_matrix(mat_file):
     df = pd.read_csv(mat_file, sep="\t", index_col=0)
     return df

def make_meta_counts_mat_helper(counts_matrices):
    counts_dfs = [read_counts_matrix(cm) for cm in counts_matrices]
    meta_count_df = pd.concat(counts_dfs, axis=1)
    meta_count_df.columns = [f"X{c}" if c[0].isdigit() else c for c in meta_count_df.columns]
    return meta_count_df

def get_conditions_safe(columns):
    cols_conds = [c.split("_")[0] for c in columns]
    two_conditions = sorted(set(cols_conds), key=cols_conds.index)
    return two_conditions

def make_meta_counts(counts_matrices, store_dir, meta_counts_outfile):
    meta_count_df = make_meta_counts_mat_helper(counts_matrices)
    # get the two conditions being tested from counts df columns
    two_conditions = get_conditions_safe(meta_count_df.columns)
    # modify storedir to include the conditions
    store_dir = os.path.join(store_dir,"vs".join(two_conditions))
    # save meta counts file
    os.makedirs(store_dir, exist_ok=True)
    meta_counts_outfile = os.path.join(store_dir, meta_counts_outfile)
    meta_count_df.iloc[:-5].to_csv(meta_counts_outfile, index=True, header=True, sep="\t")
    return list(meta_count_df.columns), store_dir, meta_counts_outfile, two_conditions
    
def get_colname_conditions(colnames, design_formula_components):
     return [list(map(str, colname.split("_")[:len(design_formula_components)])) for colname in colnames]

def generate_design_matrix(colnames, design_formula_components, savefile):
    conditions = get_colname_conditions(colnames, design_formula_components)
    with open(savefile, "w") as f:
        f.write(f",{','.join(design_formula_components)}\n")
        for col,line in zip(colnames, conditions):
            f.write(",".join([col]+line))
            f.write("\n")
    return
	
def run_deseq2(counts_file, counts_cols, designfile, design_formula, contrast, de_file):
	counts_cols_to_str = ",".join(counts_cols)
	cmd = ["Rscript", deseq2_script_path, counts_file, counts_cols_to_str, designfile, design_formula, contrast, de_file]
	subprocess.run(cmd)
	return
