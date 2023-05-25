library(DESeq2)

## good example ##
# https://lashlock.github.io/compbio/R_presentation.html
# https://www.reneshbedre.com/blog/deseq2.html
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

# get the counts file
args = commandArgs(trailingOnly=TRUE)
# check to see that two argument is given
if (length(args)!=7) {
  stop("Input counts filename, in_columns, design matrix filename, design formula, contrast, output filename and normalized counts output filename must be given", call.=FALSE)
}

infilename = args[1]
incols = args[2]
designfilename = args[3]
designformula = args[4]
contraststr = args[5]
outfilename = args[6]
normcts_file = args[7]

# read count file
cnt_table <- read.table(infilename, sep=",", row.names=1, header=TRUE)
cnt_cols = unlist(strsplit(incols, ","))
cnt_table <- cnt_table[, cnt_cols]
# design matrix or treatment data construction
design_table <- read.table(designfilename, sep=",", row.names=1, header=TRUE)
design_table[] <- lapply(design_table, factor)
# construct deseq2 object
dds <- DESeqDataSetFromMatrix(countData = cnt_table,
                              colData = design_table,
                              design = as.formula(designformula))
# run DESeq pipeline, will normalize using median of ratios
resdds <- DESeq(dds)
# get the results
contrast_vec = unlist(strsplit(contraststr, ","))
res <- results(resdds, contrast=contrast_vec)
# save to file .. 
write.table(res, file=outfilename, sep=",", row.names=TRUE, col.names=TRUE)
# get the normalized counts from deseq2 to construct pca plots
norm <- estimateSizeFactors(dds)
normalized_counts <- counts(norm, normalized=TRUE)
write.table(normalized_counts, file=normcts_file, sep=",", row.names=TRUE, col.names=TRUE) #TODO: take norm filename save as input
