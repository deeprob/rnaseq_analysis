import re
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from sklearn.decomposition import PCA

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
deseq2_script_path = os.path.join(os.path.dirname(__file__), "deseq2.R")

###################
# file processing #
###################

def create_dirs(store_dir):
    trim_dir = os.path.join(store_dir, "trim")
    align_dir = os.path.join(store_dir, "align")
    count_dir = os.path.join(store_dir, "count")
    tmp_dir = os.path.join(store_dir, "tmp")

    os.makedirs(trim_dir, exist_ok=True)
    os.makedirs(align_dir, exist_ok=True)
    os.makedirs(count_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    return trim_dir, align_dir, count_dir, tmp_dir

def remove_ext(filename):
    find = re.compile(r"^(.*?)\..*")
    basename = re.match(find, filename).group(1)
    return basename

def save_pdf(fig, savefile):
    pdf_store = PdfPages(savefile)
    pdf_store.savefig(fig, bbox_inches='tight', dpi=300)
    pdf_store.close()
    return

############
# trimming #
############

def trim_reads_helper(read_pair_1, read_pair_2, out_pair_1, out_pair_2, out_unpaired_1, out_unpaired_2, adapter_file, nthreads=64):
    # trim reads using trimmomatic
    if read_pair_2:
        # call paired end trimmomatic
        subprocess.call(["TrimmomaticPE", "-threads", f"{nthreads}", "-phred33", 
                        f"{read_pair_1}", f"{read_pair_2}", f"{out_pair_1}", f"{out_unpaired_1}", f"{out_pair_2}", f"{out_unpaired_2}",   
                        f"ILLUMINACLIP:{adapter_file}:2:30:10:7:true", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:51"])
    else:
        # call single end trimmomatic
        subprocess.call(["TrimmomaticSE", "-threads", f"{nthreads}", "-phred33", 
                        f"{read_pair_1}", f"{out_pair_1}",    
                        f"ILLUMINACLIP:{adapter_file}:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:51"])

    return


def trim_reads(in_dir, read_pair_1, read_pair_2, trim_dir, adapter_file, threads):
    trim_paired_suff = "fastq.gz"
    trim_unpaired_suff = "unpaired.fastq.gz"
    read_pair_1 = os.path.join(in_dir, read_pair_1)
    out_paired_1 = remove_ext(os.path.basename(read_pair_1))
    out_paired_1 = os.path.join(trim_dir, out_paired_1 + f".{trim_paired_suff}")
    out_unpaired_1 = os.path.join(trim_dir, out_paired_1 + f".{trim_unpaired_suff}")
    out_paired_2 = ""
    out_unpaired_2 = ""
    if read_pair_2:
        read_pair_2 = os.path.join(in_dir, read_pair_2)
        out_paired_2 = remove_ext(os.path.basename(read_pair_2))
        out_paired_2 = os.path.join(trim_dir, out_paired_2 + f".{trim_paired_suff}")
        out_unpaired_2 = os.path.join(trim_dir, out_paired_2 + f".{trim_unpaired_suff}")
    trim_reads_helper(read_pair_1, read_pair_2, out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2, adapter_file, threads)
    return out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2

#############
# alignment #
#############

def create_star_index(genome_file, gtf_file, storage_dir, nthreads=64):
    # create index for star aligner
    subprocess.call(["STAR", "--runMode", "genomeGenerate",  "--genomeDir", f"{storage_dir}",
                    "--genomeFastaFiles", f"{genome_file}", "--sjdbGTFfile", f"{gtf_file}", "--sjdbOverhang", "100", 
                    "--runThreadN", f"{nthreads}"])
    return   


def align_helper(read_files, output_prefix, starindex_dir, nthreads=64):
    # align using star
    read_files_list = read_files.split()
    command = ["STAR", "--genomeDir", f"{starindex_dir}", "--runThreadN", f"{nthreads}", "--outFileNamePrefix", f"{output_prefix}", "--readFilesIn"]
    command += read_files_list
    command += ["--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within", "--outSAMattributes", "Standard", "--readFilesCommand", "zcat"]
    subprocess.call(command)
    return


def align(in_dir, read_pair_1, read_pair_2, align_dir, index_dir, threads):
    align_suff = "Aligned.sortedByCoord.out.bam"
    read_files = read_pair_1
    if read_pair_2:
        read_files = " ".join([read_pair_1, read_pair_2])
    align_out_prefix = remove_ext(os.path.basename(read_pair_1))
    align_out_prefix = os.path.join(align_dir, align_out_prefix)
    align_helper(read_files, align_out_prefix, index_dir, threads)
    aligned_file = align_out_prefix + align_suff
    return aligned_file

############
# counting #
############

def count_helper(alignment_files, counts_cols, gtf_file, output_file, paired=False, threads=64):
    # count using htseq
    with open(output_file, "w") as outfile:
        outfile.write("\t".join(["gene_id", "gene_name"] + counts_cols) + "\n")

    with open(output_file, "a") as outfile:
        if not paired:
            command = ["htseq-count", "-f", "bam", "-s", "no", "-r", "pos", "-i", "gene_id", "--additional-attr", "gene_name", "-n", f"{threads}"] + alignment_files + [gtf_file]
        else:
            command = ["htseq-count", "-f", "bam", "-s", "yes", "-r", "pos", "-i", "gene_id", "--additional-attr", "gene_name", "-n", f"{threads}"] + alignment_files + [gtf_file]
        subprocess.run(command, stdout=outfile)
    return

def count(in_dir, aligned_file, paired, gtf_file, count_dir, threads):
    aligned_files = [aligned_file]
    counts_cols = [remove_ext(os.path.basename(aligned_file))]
    count_out_file = os.path.join(count_dir, f"{counts_cols[0]}.tsv")
    count_helper(aligned_files, counts_cols, gtf_file, count_out_file, paired=paired, threads=threads)
    return count_out_file

###########################
# differential expression #
###########################

def read_counts_matrix(counts_dir, counts_col):
    df = pd.read_csv(os.path.join(counts_dir, f"{counts_col}.tsv"), sep="\t", index_col=["gene_id", "gene_name"])
    return df

def rename_counts_columns_to_convention(old_colnames):
    new_colnames = [f"X{c.replace('-', '_')}" if c[0].isdigit() else c for c in old_colnames]
    return dict(zip(old_colnames, new_colnames))

def make_meta_counts_mat_helper(counts_dir, counts_cols):
    counts_dfs = [read_counts_matrix(counts_dir, cc) for cc in counts_cols]
    meta_count_df = pd.concat(counts_dfs, axis=1)
    # rename all colnames starting with digits by adding X before (stupid R and deseq2)
    # replace all dashes with underscore (idiotic R and deseq2)
    meta_count_df.columns = [f"X{c.replace('-', '_')}" if c[0].isdigit() else c for c in meta_count_df.columns]
    return meta_count_df

def make_meta_counts(
        counts_dir, 
        counts_cols,
        meta_counts_outfile
        ):
    meta_count_df = make_meta_counts_mat_helper(counts_dir, counts_cols)
    # save meta counts file
    meta_count_df.iloc[:-5].to_csv(meta_counts_outfile, index=True, header=True)
    return dict(zip(counts_cols, list(meta_count_df.columns)))

def get_formula(design_matrix_file, count_col_dict, de_dir):
    df = pd.read_csv(design_matrix_file, index_col=0)
    df = df.rename(index=count_col_dict)
    dmf_savefile = os.path.join(de_dir, "design_matrix.csv")
    df.to_csv(dmf_savefile)
    return list(df.columns), f"~{'+'.join(list(df.columns))}", dmf_savefile

def run_deseq2(counts_file, counts_cols, designfile, design_formula, contrast, de_file, normcts_file, mingenecounts):
	counts_cols_to_str = ",".join(counts_cols)
	cmd = ["Rscript", deseq2_script_path, counts_file, counts_cols_to_str, designfile, design_formula, contrast, de_file, normcts_file, str(mingenecounts)]
	subprocess.run(cmd)
	return

############
# pca plot #
############

def read_design_matrix(filename):
    dm = pd.read_csv(filename, index_col=0)
    return dm

def get_pca_components(meta_df, designmatrixfile):
    meta_df_norm = (meta_df-meta_df.mean())/meta_df.std()
    dm_df = read_design_matrix(designmatrixfile)
    pca = PCA(n_components=2)
    components = pca.fit_transform(meta_df_norm.T)
    pca_df = pd.DataFrame({
        "lib_info": meta_df_norm.columns,
        "pca_comp1": components[:, 0],
        "pca_comp2": components[:, 1],
        "library": [dm_df.loc[c, dm_df.columns[0]] for c in meta_df_norm.columns],
    })
    if dm_df.shape[1]>1:
        pca_df["factor"] = [dm_df.loc[c, dm_df.columns[1]] for c in meta_df_norm.columns]
    return pca_df, dm_df.columns[1] if dm_df.shape[1]>1 else ""

def get_pca_plots_helper(meta_df, designmatrixfile):
    pca_df, factor = get_pca_components(meta_df, designmatrixfile) 
    fig, ax = plt.subplots(1, 1, figsize=(8,4), sharex=True, sharey=True)
    sns.scatterplot(
        data=pca_df, 
        x="pca_comp1", 
        y="pca_comp2", 
        hue="library", 
        style="factor" if factor else "library", 
        legend=True, 
        ax=ax, 
        s=100,
        alpha=1.,
        linewidth=1.05,
        edgecolor="k",
        )
    ax.set_title("PCA plot")
    ax.legend(loc="lower center", markerscale=2,  prop={'size': 12})
    return fig

def create_lib_specific_pca_plot(de_dir, norm_counts_file, designmatrixfile):
    norm_counts_df = pd.read_csv(norm_counts_file, index_col=0)
    designmatrixfile = os.path.join(designmatrixfile)
    f = get_pca_plots_helper(norm_counts_df, designmatrixfile)
    save_file = os.path.join(de_dir, f'pca.pdf')
    save_pdf(f, save_file)
    return

################
# volcano plot #
################

def read_deseq_results(deseq_outfile, meta_countsfile):
    deres_df = pd.read_csv(deseq_outfile)
    gid2n_df = pd.read_csv(meta_countsfile, index_col=0, usecols=["gene_id", "gene_name"])
    df = pd.concat((deres_df, gid2n_df), axis=1, join="inner")
    return df


def parse_deseqres_for_volcano_plot(df, lfc_thresh, pv_thresh):
    # drop rows with na values
    df = df.dropna()
    # convert all 0 padj values to half of the min padj value thats greater than 0
    df.loc[df.padj==0, "padj"] = min(df.loc[df.padj>0].padj)/10
    # create neglog10 val
    df["neglog10padj"] = -np.log10(df.padj)
    df["hue"] = "Not Significant"
    df.loc[(df.log2FoldChange>lfc_thresh)&(df.padj<pv_thresh), "hue"] = "Significant Up"
    df.loc[(df.log2FoldChange<-lfc_thresh)&(df.padj<pv_thresh), "hue"] = "Significant Down"
    return df

def create_volcano_fig_raw(df_volcano, lfc_thresh, pv_thresh, gene_set=[]):

    if not gene_set:
        gene_set = df_volcano.sort_values("padj").head().gene_name.values

    fig, axes = plt.subplots(1, 1, figsize=(6, 8))
    sns.scatterplot(
        data=df_volcano, x="log2FoldChange", y="neglog10padj", 
        hue="hue", palette={"Significant Down": "green", "Not Significant": "grey", "Significant Up": "red"},
        size="hue",
        sizes={"Significant Up": 20, "Not Significant":2, "Significant Down": 20}, 
        rasterized=True,
        ax=axes,
        legend=False,
        )
    axes.axvline(x=lfc_thresh, linestyle="--", color="k")
    axes.axvline(x=-lfc_thresh, linestyle="--", color="k")
    axes.axhline(y=-np.log10(pv_thresh), linestyle="--", color="k")

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
        y_coord -= 1 #axes.get_ylim()[1]//20 # division by 20 needs to be edited to fit the limits correctly

    axes.set_xlim(axes.get_xlim()[0], axes.get_xlim()[1]+ xcoord_shift + axes.get_xlim()[1]/20)
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.collections[0].set_rasterized(True)
    return fig

def create_volcano_plot(de_dir, deseq_outfile, meta_countsfile, gene_set=[], lfc_thresh=0.25, pv_thresh=0.05):
    # get deseq results
    de_df = read_deseq_results(deseq_outfile, meta_countsfile)
    # parse deseq results
    de_df = parse_deseqres_for_volcano_plot(de_df, lfc_thresh, pv_thresh)
    # save volcano table
    save_file = os.path.join(de_dir, f'volcano_table.csv')
    de_df.to_csv(save_file, index=True, header=True)
    # get volcano plot figure
    vfig = create_volcano_fig_raw(de_df, gene_set=gene_set, lfc_thresh=lfc_thresh, pv_thresh=pv_thresh)
    save_file = os.path.join(de_dir, f'volcano.pdf')
    save_pdf(vfig, save_file)
    return

