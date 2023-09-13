args <- commandArgs(trailingOnly=TRUE)

counts_filename <- args[1]
outname <- args[2]

counts <- read.delim(counts_filename,header=TRUE,sep="\t")

genes <- data.frame(c1=unique(counts$Gene_id),c2=unique(counts$Gene_id))
barcodes <- data.frame(c1=unique(counts$barcode))

counts$matchgene <- match(counts$Gene_id,genes$c1)
counts$matchbarcode <- match(counts$barcode,barcodes$c1)


dir.create(outname)

write.table(genes,file=paste0(outname,"/features.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
write.table(barcodes,file=paste0(outname,"/barcodes.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)

cat("%%MatrixMarket matrix coordinate integer general\n%\n",file=paste0(outname,"/matrix.mtx"))
cat(nrow(genes)," ",nrow(barcodes)," ",nrow(counts),"\n",file=paste0(outname,"/matrix.mtx"),append=TRUE,sep="")
write.table(counts[,c(4,5,3)],file=paste0(outname,"/matrix.mtx"),sep=" ",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)


