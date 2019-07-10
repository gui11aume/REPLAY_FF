genes = read.delim("EnsemblGenes_mm9.txt", as.is=TRUE)

start = tapply(X=genes$txStart, INDEX=genes$name2, min)
end = tapply(X=genes$txEnd, INDEX=genes$name2, max)
chr = tapply(X=genes$chrom, INDEX=genes$name2, unique)
strand = tapply(X=genes$strand, INDEX=genes$name2, unique)

stopifnot(length(unlist(chr)) == length(chr))
stopifnot(length(unlist(strand)) == length(strand))

ENSG = data.frame(ENSG = names(start), chrom = chr,
   start = start, end = end, strand = strand)

head(ENSG)
write.table(ENSG, file="ENSG_mm9.txt", sep="\t", quote=FALSE,
   row.names=FALSE)
