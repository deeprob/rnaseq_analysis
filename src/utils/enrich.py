import os
import subprocess
import pandas as pd


CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))


################
# get de genes #
################

def read_deseq_results(de_dir, treatment, control, deseq_filebase, gid2n_filebase):
    deres_filename = os.path.join(de_dir, "vs".join([treatment, control]), deseq_filebase)
    gid2n_filename = os.path.join(de_dir, "vs".join([treatment, control]), gid2n_filebase)
    deres_df = pd.read_csv(deres_filename)
    gid2n_df = pd.read_csv(gid2n_filename, index_col=0)
    df = pd.concat((deres_df, gid2n_df), axis=1)
    # drop rows with na values
    df = df.dropna()
    return df

def get_de_genes_helper(de_dir, treatment, control, deseq_filebase, gid2n_filebase, degenes_dir, thresh=0.01):
    df = read_deseq_results(de_dir, treatment, control, deseq_filebase, gid2n_filebase)
    df = df.loc[df.padj<thresh]
    sig_degenes = [g.split(".")[0] for g in df.index]

    save_file = os.path.join(degenes_dir, f'{treatment}vs{control}_sig_genes.txt')
    with open(save_file, "w") as f:
        for g in sig_degenes:
            f.write(f"{g}\n")
    return

def get_de_genes(de_dir, treatments, controls, deseq_filebases, gid2n_filebases, store_dir, thresh=0.01):
    degenes_dir = os.path.join(store_dir, "tables")
    os.makedirs(degenes_dir, exist_ok=True)
    for t,c,d,g in zip(treatments, controls, deseq_filebases, gid2n_filebases):
        get_de_genes_helper(de_dir, t, c, d, g, degenes_dir, thresh)
    return degenes_dir


########################
# gsea kegg enrichment #
########################

def run_enrichment_helper(gene_in, gseaout_file, keggout_file, gseafigout_file, keggfigout_file, tmp_dir):
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/enrich.sh", 
        gene_in, gseaout_file, keggout_file, gseafigout_file, keggfigout_file, tmp_dir
        ]
    subprocess.run(cmd)
    return

def run_enrichment(de_genes_dir, treatments, controls, store_dir):
    table_dir = os.path.join(store_dir, "tables")
    os.makedirs(store_dir, exist_ok=True)
    fig_dir = os.path.join(store_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    tmp_dir = os.path.join(store_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    for t, c in zip(treatments, controls):
        gene_in = os.path.join(de_genes_dir, f'{t}vs{c}_sig_genes.txt')
        gseaout_file = os.path.join(table_dir, f'{t}vs{c}_gsea.csv')
        keggout_file = os.path.join(table_dir, f'{t}vs{c}_kegg.csv')
        gseafigout_file = os.path.join(fig_dir, f'{t}vs{c}_gsea.pdf')
        keggfigout_file = os.path.join(fig_dir, f'{t}vs{c}_kegg.pdf')
        run_enrichment_helper(gene_in, gseaout_file, keggout_file, gseafigout_file, keggfigout_file, tmp_dir)
    return

