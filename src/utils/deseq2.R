library(DESeq2)

## good example ##
# https://lashlock.github.io/compbio/R_presentation.html

# get the counts file
args = commandArgs(trailingOnly=TRUE)
# check to see that two argument is given
if (length(args)!=6) {
  stop("Input counts filename, in_columns, design matrix filename, design formula contrast and output filename must be given", call.=FALSE)
}

infilename = args[1]
incols = args[2]
designfilename = args[3]
designformula = args[4]
contraststr = args[5]
outfilename = args[6]

# read count file
cnt_table <- read.table(infilename, sep="\t", row.names=1, header=TRUE)
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
dds <- DESeq(dds)
# get the results
contrast_vec = unlist(strsplit(contraststr, ","))
res <- results(dds, contrast=contrast_vec)
# save to file .. 
write.table(res, file=outfilename, sep=",", row.names=TRUE, col.names=TRUE)
