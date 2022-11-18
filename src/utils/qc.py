import os
import pandas as pd
import numpy as np
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
    return df


def parse_deseqres_for_volcano_plot(df):
    # drop rows with na values
    df = df.dropna()
    # convert all 0 padj values to half of the min padj value thats greater than 0
    df.loc[df.padj==0, "padj"] = min(df.loc[df.padj>0].padj)/10
    # create neglog10 val
    df["neglog10padj"] = -np.log10(df.padj)
    # create hue columns
    lfc_thresh = 1
    pv_thresh = 0.01
    df["hue"] = "Not Significant"
    df.loc[(df.log2FoldChange>lfc_thresh)&(df.padj<pv_thresh), "hue"] = "Significant Up"
    df.loc[(df.log2FoldChange<-lfc_thresh)&(df.padj<pv_thresh), "hue"] = "Significant Down"
    return df

def create_volcano_fig_raw(df_volcano, lfc_thresh=1, pv_thresh=0.01, gene_set=[]):
    fig, axes = plt.subplots(1, 1, figsize=(6, 8))
    sns.scatterplot(
        data=df_volcano, x="log2FoldChange", y="neglog10padj", 
        hue="hue", palette={"Significant Down": "green", "Not Significant": "grey", "Significant Up": "red"},
        ax=axes,
        legend=False,
        )
    axes.axvline(x=lfc_thresh, linestyle="--", color="k")
    axes.axvline(x=-lfc_thresh, linestyle="--", color="k")
    axes.axhline(y=-np.log10(pv_thresh), linestyle="--", color="k")

    if gene_set:
        xcoord_shift = 0.25
        ycoord_shift = 0.25
        y_coord = max(df_volcano.neglog10padj) + ycoord_shift
        x_coord = max(df_volcano.log2FoldChange) + xcoord_shift

        for info in df_volcano.loc[df_volcano.gene_name.isin(gene_set)].itertuples():
            axes.annotate(
                info.gene_name, 
                xy=(info.log2FoldChange, info.neglog10padj), xytext=(x_coord, y_coord), 
                arrowprops={"arrowstyle": "->", "lw": 1, "color": "black", "ls": "--", "relpos": (0,0.5)}, 
                fontsize=12)
            y_coord -= axes.get_ylim()[1]//20

        axes.set_xlim(axes.get_xlim()[0], axes.get_xlim()[1]+ xcoord_shift + axes.get_xlim()[1]/20)
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.collections[0].set_rasterized(True)
    return fig

def create_volcano_plot_helper(de_dir, treatment, control, deseq_filebase, gid2n_filebase, store_dir, gene_set=[]):
    # make appropriate dirs to store files
    table_dir = os.path.join(store_dir, "tables")
    os.makedirs(table_dir, exist_ok=True)
    figure_dir = os.path.join(store_dir, "figures")
    os.makedirs(figure_dir, exist_ok=True)
    # get deseq results
    de_df = read_deseq_results(de_dir, treatment, control, deseq_filebase, gid2n_filebase)
    # parse deseq results
    de_df = parse_deseqres_for_volcano_plot(de_df)
    # save volcano table
    save_file = os.path.join(table_dir, f'{treatment}vs{control}_volcano_table.csv')
    de_df.to_csv(save_file, index=True, header=True)
    # get volcano plot figure
    vfig = create_volcano_fig_raw(de_df, gene_set=gene_set)
    save_file = os.path.join(figure_dir, f'{treatment}vs{control}_volcano.pdf')
    save_pdf(vfig, save_file)
    return

def create_volcano_plots(de_dir, treatments, controls, deseq_filebases, geneidtonamefiles, store_dir, genenames):
    for t,c,d,g in zip(treatments, controls, deseq_filebases, geneidtonamefiles):
        create_volcano_plot_helper(de_dir, t, c, d, g, store_dir, genenames)
    return

# ma plot
