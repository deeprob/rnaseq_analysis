#### TMM normalization script ####
library(limma)
library(edgeR)

# get the counts file
args = commandArgs(trailingOnly=TRUE)
# check to see that only one argument is given
if (length(args)!=2) {
  stop("Two arguments must be supplied (input file).n", call.=FALSE)
}

infilename = args[1]
outfilename = args[2]

# read count file
tabla <- read.table(infilename,sep="\t",row.names=1)

# convert to edger dge object
dge <- DGEList(counts=tabla)
# calculate TMM normalized counts
dge <- calcNormFactors(dge, method="TMM")

# save to file .. 
write.table(dge$counts, file=outfilename, sep="\t", row.names=TRUE, col.names=F)
