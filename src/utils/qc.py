import os
from itertools import zip_longest
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
    os.makedirs(qc_dir, exist_ok=True)
    return qc_dir


############
# pca plot #
############

def read_counts_file(counts_dir, lib, filename):
    df = pd.read_csv(os.path.join(counts_dir, lib, filename), sep="\t", skipfooter=5, index_col=0, engine="python")
    return df.drop(columns="gene_name")

def get_meta_counts(counts_dir, libraries, count_filenames):
    df = pd.concat([read_counts_file(counts_dir, lib, cf) for lib, cf in zip_longest(libraries, count_filenames, fillvalue=count_filenames[0])], axis=1)
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
    save_file = os.path.join(store_dir, "library_pca.pdf")
    save_pdf(f, save_file)
    return

# volcano plot
# ma plot