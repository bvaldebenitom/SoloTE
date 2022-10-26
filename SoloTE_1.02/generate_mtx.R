args <- commandArgs(trailingOnly=TRUE)

barcodes_filename <- args[1]
genes_filename <- args[2]
counts_filename <- args[3]

outname <- args[4]

barcodes <- read.delim(barcodes_filename,header=FALSE)
genes <- read.delim(genes_filename,header=FALSE)
counts <- read.delim(counts_filename,header=FALSE,sep="\t")
counts$matchgene <- match(counts$V1,genes$V1)
counts$matchbarcode <- match(counts$V2,barcodes$V1)
write.table(file=outname,counts[,c(4,5,3)],sep=" ",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
