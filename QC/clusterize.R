compare = function(df1, df2) {
   # Compare two .co.gz files. The similarity score is the weighted
   # covariance among AT / (AT + GC), with weight AT + GC. The barcodes
   # are merged between experiments and these scores are computed. The
   # barcodes present in only one experiment are given weight 0.

   if (nrow(df1) < 10 || nrow(df2) < 10) {
      return (0)
   }

   # Columns are: barcode, FF, AT, GC.
   x1 = df1[,3] / (df1[,3] + df1[,4])
   x2 = df2[,3] / (df2[,3] + df2[,4])

   df1 = data.frame(bcd=df1[,1],
               x1 = df1[,3]/(df1[,3]+df1[,4]),
               w1 = df1[,3]+df1[,4])
   df2 = data.frame(bcd=df2[,1],
               x2 = df2[,3]/(df2[,3]+df2[,4]),
               w2 = df2[,3]+df2[,4])

   # Match barcodes across experiments.
   m = merge(df1, df2, by="bcd", all=TRUE)
   m[is.na(m)] = 0

   cov.wt(m[,c('x1','x2')], wt=m$w1+m$w2,
          center=FALSE, cor=TRUE)$cor[1,2]

}

sink("clusterize_report.txt")
fnames = commandArgs(trailingOnly = TRUE)

if (!file.exists("sim.rda")) {

   # Compute the similarity matrix.
   barcodes = list(
      CA = c(read.table("../mapping/CA1.ins", as.is=TRUE)$V1,
             read.table("../mapping/CA2.ins", as.is=TRUE)$V1),
      CT = c(read.table("../mapping/CT1.ins", as.is=TRUE)$V1,
             read.table("../mapping/CT2.ins", as.is=TRUE)$V1),
      GA = c(read.table("../mapping/GA1.ins", as.is=TRUE)$V1,
             read.table("../mapping/GA2.ins", as.is=TRUE)$V1),
      GT = c(read.table("../mapping/GT1.ins", as.is=TRUE)$V1,
             read.table("../mapping/GT2.ins", as.is=TRUE)$V1)
   )
   # Add other mismatch codes.
   barcodes[['24']] = barcodes[['48']] =
      barcodes[['LA']] = barcodes[['CT']]

   files = list()
   enames = sub(".*?([^/]+)\\.co\\.gz", "\\1", fnames)

   for (i in 1:length(fnames)) {
      # Filter contaminating barcodes.
      mmcode = toupper(substr(enames[i], 1,2))
      files[[i]] = subset(read.delim(fnames[i], as.is=TRUE),
                       barcode %in% barcodes[[mmcode]])
      # Check the overlap between the barcodes in the .co.gz files
      # and those in the .ins file. If it is too low, the sample does
      # not have the riht identity.
      if (nrow(files[[i]]) < 10) {
         cat(paste("sample swap:", enames[i], "\n"))
      }
   }

   # Create similarity matrix between .co.gz files.
   # The similarity score is the output of the 'compare()' function.
   sim = matrix(0, nrow=length(fnames), ncol=length(fnames))
   for (i in 1:(length(fnames)-1)) {
   for (j in (i+1):length(fnames)) {
      sim[i,j] = sim[j,i] = compare(files[[i]], files[[j]])
   }
   }

   colnames(sim) = enames
   save(sim, file="sim.rda")

} else {
   load("sim.rda")
}

h = hclust(as.dist(1-sim))

pdf("clusters.pdf", useDingbats=FALSE, width=18, height=8)
plot(h, xlab="Correlation between experiments")
#pdf("heatmap.pdf", useDingbats=FALSE)
#image(sim[h$order,h$order][1:35,1:35], col=colorRampPalette(c("white", "black"))(256))
dev.off()
sink()
