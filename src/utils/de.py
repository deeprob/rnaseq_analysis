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

def read_counts_matrix(dir, short, name, pre, rep):
    df = pd.read_csv(os.path.join(dir, short, name), sep="\t", index_col="gene_id")
    gene_id_to_name_map = df.loc[:, ["gene_name"]].reset_index()
    df = df.loc[:, ["_".join([pre, r]) for r in rep.split()]]
    return df, gene_id_to_name_map

def make_meta_counts_mat_helper(counts_dir, control_shortform, treatment_shortform, controlcountsmat, treatmentcountsmat, control_pre, treatment_pre, control_reps, treatment_reps):
    counts_dfs_info = [read_counts_matrix(counts_dir, s, n, p, r) for s,n,p,r in zip([control_shortform, treatment_shortform], [controlcountsmat, treatmentcountsmat], [control_pre, treatment_pre], [control_reps, treatment_reps])]
    meta_count_df = pd.concat([cdf for cdf,_ in counts_dfs_info], axis=1)
    gene_id2name_maps = [id2n for _,id2n in counts_dfs_info]
    assert gene_id2name_maps[0].equals(gene_id2name_maps[1])
    meta_count_df.columns = [f"X{c}" if c[0].isdigit() else c for c in meta_count_df.columns]
    return meta_count_df, gene_id2name_maps[0]

def make_meta_counts(
        counts_dir, 
        control_shortform, treatment_shortform, 
        controlcountsmat, treatmentcountsmat, 
        control_pre, treatment_pre,
        control_reps, treatment_reps,
        store_dir, meta_counts_outfile
        ):
    meta_count_df, gid2name_df = make_meta_counts_mat_helper(
            counts_dir, 
            control_shortform, treatment_shortform, 
            controlcountsmat, treatmentcountsmat,
            control_pre, treatment_pre,
            control_reps, treatment_reps,
            )
    # modify storedir to include the conditions
    store_dir = os.path.join(store_dir,"vs".join([treatment_shortform, control_shortform]))
    # save meta counts file
    os.makedirs(store_dir, exist_ok=True)
    meta_counts_outfile = os.path.join(store_dir, meta_counts_outfile)
    meta_count_df.iloc[:-5].to_csv(meta_counts_outfile, index=True, header=True, sep="\t")
    # save gene id to name mapping file
    gid2name_outfile = os.path.join(store_dir, "geneid2name.csv")
    gid2name_df.iloc[:-5].to_csv(gid2name_outfile, index=False, header=True)
    return list(meta_count_df.columns), store_dir, meta_counts_outfile
    
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
