csize = read.table("csize.txt", row.names=1)

INDEX = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
   "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
   "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
   "chr19", "chrX", "chrY")

CA1 = subset(read.table("../mapping/CA1.ins", as.is=TRUE), V2 %in% INDEX)
CA2 = subset(read.table("../mapping/CA2.ins", as.is=TRUE), V2 %in% INDEX)
CT1 = subset(read.table("../mapping/CT1.ins", as.is=TRUE), V2 %in% INDEX)
CT2 = subset(read.table("../mapping/CT2.ins", as.is=TRUE), V2 %in% INDEX)
TC1 = subset(read.table("../mapping/TC1.ins", as.is=TRUE), V2 %in% INDEX)
TC2 = subset(read.table("../mapping/TC2.ins", as.is=TRUE), V2 %in% INDEX)
GA1 = subset(read.table("../mapping/GA1.ins", as.is=TRUE), V2 %in% INDEX)
GA2 = subset(read.table("../mapping/GA2.ins", as.is=TRUE), V2 %in% INDEX)
GT1 = subset(read.table("../mapping/GT1.ins", as.is=TRUE), V2 %in% INDEX)
GT2 = subset(read.table("../mapping/GT2.ins", as.is=TRUE), V2 %in% INDEX)

COL = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")

CA = table(c(CA1$V2, CA2$V2)) / (nrow(CA1) + nrow(CA2))
CT = table(c(CT1$V2, CT2$V2)) / (nrow(CT1) + nrow(CT2))
TC = table(c(TC1$V2, TC2$V2)) / (nrow(TC1) + nrow(TC2))
GA = table(c(GA1$V2, GA2$V2)) / (nrow(GA1) + nrow(GA2))
GT = table(c(GT1$V2, GT2$V2)) / (nrow(GT1) + nrow(GT2))

set.seed(123)

PCH=c(rep(20,19), 17, 18)

pdf("coverage.pdf", width=2)
par(mar=c(1,3.5,.5,.5))
plot(rnorm(n=21), (CA*sum(csize)/csize)$V2, col=COL[1], pch=PCH,
   bty="n", xaxt="n", yaxt="n", xlim=c(-3,3), col.lab="grey25",
   line=2.2, ylab="Relative insertion rate", xlab="")
points(rnorm(n=21), (CT*sum(csize)/csize)$V2, col=COL[2], pch=PCH)
points(rnorm(n=21), (GA*sum(csize)/csize)$V2, col=COL[3], pch=PCH)
points(rnorm(n=21), (GT*sum(csize)/csize)$V2, col=COL[4], pch=PCH)
points(rnorm(n=21), (TC*sum(csize)/csize)$V2, col=COL[5], pch=PCH)
text(x=c(0,-.5), y=c(.1, .47), labels=c("chrY", "chrX"),
   col="grey30", cex=.9)
axis(side=2, col="grey50", cex.axis=.8, col.axis="grey25",
   at=c(0,1,1.3))
dev.off()
