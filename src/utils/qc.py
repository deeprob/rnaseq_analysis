import os
import pandas as pd
from bioinfokit import visuz
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.decomposition import PCA
import seaborn as sns


def save_pdf(fig, savefile):
    pdf_store = PdfPages(savefile)
    pdf_store.savefig(fig, bbox_inches='tight')
    pdf_store.close()
    return

def create_dirs(store_dir):
    qc_dir = os.path.join(store_dir, "qc")
    enrich_dir = os.path.join(store_dir, "enrich")
    os.makedirs(qc_dir, exist_ok=True)
    os.makedirs(enrich_dir, exist_ok=True)
    return qc_dir, enrich_dir

def match_lengths(iter1, iter2):
    if len(iter2) != len(iter1):
        iter2 = [iter2[0] for _ in range(len(iter1))]
    return iter1, iter2

############
# pca plot #
############

def read_counts_file(counts_dir, lib, filename):
    df = pd.read_csv(os.path.join(counts_dir, lib, filename), sep="\t", skipfooter=5, index_col=0, engine="python")
    return df.drop(columns="gene_name")

def get_meta_counts(counts_dir, libraries, count_filenames):
    df = pd.concat([read_counts_file(counts_dir, lib, cf) for lib, cf in zip(libraries, count_filenames)], axis=1)
    return df

def get_pca_components(meta_df):
    meta_df_norm = (meta_df-meta_df.mean())/meta_df.std()
    pca = PCA(n_components=2)
    components = pca.fit_transform(meta_df_norm.T)
    pca_df = pd.DataFrame({
        "lib_info": meta_df_norm.columns,
        "pca_comp1": components[:, 0],
        "pca_comp2": components[:, 1],
        "library": [c.split("_")[0] for c in meta_df_norm.columns],
    })
    return pca_df

def get_pca_plots_helper(meta_df):
    pca_df = get_pca_components(meta_df) 
    fig, ax = plt.subplots(1, 1, figsize=(8,4), sharex=True, sharey=True)
    sns.scatterplot(
        data=pca_df, 
        x="pca_comp1", 
        y="pca_comp2", 
        hue="library", 
        style="library", 
        legend=True, 
        ax=ax, 
        s=150,
        alpha=1.,
        linewidth=1.05,
        edgecolor="k",
        )
    ax.set_title("PCA plot")
    ax.legend(loc="lower center", markerscale=2,  prop={'size': 12})
    return fig

def save_pca_plot(counts_dir, libraries, count_filenames, store_dir):
    meta_df = get_meta_counts(counts_dir, libraries, count_filenames)
    f = get_pca_plots_helper(meta_df)
    store_dir = os.path.join(store_dir, "figures")
    os.makedirs(store_dir, exist_ok=True)
    save_file = os.path.join(store_dir, "library_pca.pdf")
    save_pdf(f, save_file)
    return


################
# volcano plot #
################

def read_deseq_results(de_dir, treatment, control, deseq_filebase, gid2n_filebase):
    deres_filename = os.path.join(de_dir, "vs".join([treatment, control]), deseq_filebase)
    gid2n_filename = os.path.join(de_dir, "vs".join([treatment, control]), gid2n_filebase)
    deres_df = pd.read_csv(deres_filename)
    gid2n_df = pd.read_csv(gid2n_filename, index_col=0)
    df = pd.concat((deres_df, gid2n_df), axis=1)
    # drop rows with na values
    df = df.dropna()
    # convert all 0 padj values to half of the min padj value thats greater than 0
    df.loc[df.padj==0, "padj"] = min(df.loc[df.padj>0].padj)/10
    return df


def create_volcano_plot_helper(de_dir, treatment, control, deseq_filebase, gid2n_filebase, store_dir):
    table_dir = os.path.join(store_dir, "tables")
    os.makedirs(table_dir, exist_ok=True)
    figure_dir = os.path.join(store_dir, "figures")
    os.makedirs(figure_dir, exist_ok=True)
    de_df = read_deseq_results(de_dir, treatment, control, deseq_filebase, gid2n_filebase)
    visuz.gene_exp.volcano(
        df=de_df,
        lfc='log2FoldChange', 
        pv='padj',
        geneid="gene_name",
        genenames=tuple(de_df.sort_values("padj", ascending=True).gene_name)[:5],
        gstyle=2,
        gfont=20,
        show=False,
        dotsize=4,
        dim=(8, 12),
        sign_line=True,
        plotlegend=True,
        axtickfontsize=12,
        axlabelfontsize=14,
        figtype="pdf",
        figname=os.path.join(figure_dir, f'{treatment}vs{control}_volcano'),
        )
    save_file = os.path.join(table_dir, f'{treatment}vs{control}_volcano_table.csv')
    de_df.to_csv(save_file, index=True, header=True)
    return

def create_volcano_plots(de_dir, treatments, controls, deseq_filebases, geneidtonamefiles, store_dir):
    for t,c,d,g in zip(treatments, controls, deseq_filebases, geneidtonamefiles):
        create_volcano_plot_helper(de_dir, t, c, d, g, store_dir)
    return

# ma plot