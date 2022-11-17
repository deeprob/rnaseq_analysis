library("biomaRt")
library("clusterProfiler")


args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=2) {
  stop("gene_list filename and output filename must be given", call.=FALSE)
}

gene_file = args[1]
out_file = args[2]

# get the genes and store as a vector
genes = read.table(gene_file, header=FALSE)

# convert ensemble id to entrez id :: link: https://support.bioconductor.org/p/114325/
mart <- useMart("ensembl","hsapiens_gene_ensembl")
entrez_genes <- getBM(c("ensembl_gene_id", "entrezgene_id"), "ensembl_gene_id", genes, mart)

# use kegg for go term analysis :: link: http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
kenrich = enrichKEGG(
    gene=entrez_genes[, 2],
    organism="hsa",
    pvalueCutoff=0.05,
    pAdjustMethod="BH"
)

# save to file .. 
write.table(kenrich, file=out_file, sep=",", row.names=TRUE, col.names=TRUE)

