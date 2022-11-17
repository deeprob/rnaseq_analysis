library("biomaRt")
library("clusterProfiler")
library("enrichplot")


args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=6) {
  stop("gene_list, gsea table, kegg table, gsea figure, kegg figure filenames and tmp_dir must be given", call.=FALSE)
}

gene_file = args[1]
gseaout_file = args[2]
keggout_file = args[3]
gseafigout_file = args[4]
keggfigout_file = args[5]
tmp_dir = args[6]

# get the genes and store as a vector
genes = read.table(gene_file, header=FALSE)

Sys.setenv(BIOMART_CACHE = tmp_dir)
print(biomartCacheInfo())

# convert ensemble id to entrez id :: link: https://support.bioconductor.org/p/114325/
mart <- useMart("ensembl","hsapiens_gene_ensembl")
entrez_genes <- getBM(c("ensembl_gene_id", "entrezgene_id"), "ensembl_gene_id", genes, mart)

# use enrichGO for go term analysis :: link: http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
goenrich = enrichGO(
    gene=entrez_genes[, 2],
    OrgDb='org.Hs.eg.db',
    pAdjustMethod="BH",
    pvalueCutoff=0.05,
    ont="BP"
)

pdf(gseafigout_file)
myplot <- dotplot(goenrich, showCategory=30) + ggtitle("dotplot for GSEA")
print(myplot)
dev.off()

# use kegg for go term analysis :: link: http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
kenrich = enrichKEGG(
    gene=entrez_genes[, 2],
    organism="hsa",
    pvalueCutoff=0.05,
    pAdjustMethod="BH"
)

pdf(keggfigout_file)
myplot <- dotplot(kenrich, showCategory=30) + ggtitle("dotplot for KEGG")
print(myplot)
dev.off()

# save to file .. 
write.table(goenrich, file=gseaout_file, sep=",", row.names=TRUE, col.names=TRUE)

# save to file .. 
write.table(kenrich, file=keggout_file, sep=",", row.names=TRUE, col.names=TRUE)

# make plots as well

