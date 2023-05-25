import os
import subprocess
import pandas as pd

deseq2_script_path = os.path.join(os.path.dirname(__file__), "deseq2.R")


###############
# dir parsing #
###############

def create_dirs(store_dir):
    de_dir = store_dir
    os.makedirs(de_dir, exist_ok=True)
    return de_dir


####################
# meta counts file #
####################

def read_counts_matrix(counts_dir, counts_col):
    df = pd.read_csv(os.path.join(counts_dir, f"{counts_col}.tsv"), sep="\t", index_col="gene_id")
    gene_id_to_name_map = df.loc[:, ["gene_name"]].reset_index()
    df = df.drop(columns=["gene_name"])
    return df, gene_id_to_name_map

def make_meta_counts_mat_helper(counts_dir, counts_cols):
    counts_dfs_info = [read_counts_matrix(counts_dir, cc) for cc in counts_cols]
    meta_count_df = pd.concat([cdf for cdf,_ in counts_dfs_info], axis=1)
    gene_id2name_maps = [id2n for _,id2n in counts_dfs_info]
    assert gene_id2name_maps[0].equals(gene_id2name_maps[1])
    meta_count_df.columns = [f"X{c}" if c[0].isdigit() else c for c in meta_count_df.columns]
    return meta_count_df, gene_id2name_maps[0]

def make_meta_counts(
        counts_dir, 
        counts_cols,
        meta_counts_outfile
        ):
    meta_count_df, gid2name_df = make_meta_counts_mat_helper(
            counts_dir, 
            counts_cols
            )
    # save meta counts file
    meta_count_df.iloc[:-5].to_csv(meta_counts_outfile, index=True, header=True, sep="\t")
    # save gene id to name mapping file
    gid2name_outfile = os.path.join(os.path.dirname(meta_counts_outfile), "geneid2name.csv")
    gid2name_df.iloc[:-5].to_csv(gid2name_outfile, index=False, header=True)
    return


#################
# design matrix #
#################

def get_formula(design_matrix_file):
    df = pd.read_csv(design_matrix_file, nrows=1, index_col=0)
    return list(df.columns), f"~{'+'.join(list(df.columns))}"

#############
# deseq run #
#############

def run_deseq2(counts_file, counts_cols, designfile, design_formula, contrast, de_file, normcts_file):
	counts_cols_to_str = ",".join(counts_cols)
	cmd = ["Rscript", deseq2_script_path, counts_file, counts_cols_to_str, designfile, design_formula, contrast, de_file, normcts_file]
	subprocess.run(cmd)
	return
