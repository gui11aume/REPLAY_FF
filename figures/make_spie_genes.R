library(GenomicRanges)
source("spie.R")

setwd("../mapping")
CA1 = subset(read.table("../mapping/CA1.ins"), V2 != "pT2")
CA2 = subset(read.table("../mapping/CA2.ins"), V2 != "pT2")
CT1 = subset(read.table("../mapping/CT1.ins"), V2 != "pT2")
CT2 = subset(read.table("../mapping/CT2.ins"), V2 != "pT2")
GA1 = subset(read.table("../mapping/GA1.ins"), V2 != "pT2")
GA2 = subset(read.table("../mapping/GA2.ins"), V2 != "pT2")
GT1 = subset(read.table("../mapping/GT1.ins"), V2 != "pT2")
GT2 = subset(read.table("../mapping/GT2.ins"), V2 != "pT2")
TC1 = subset(read.table("../mapping/TC1.ins"), V2 != "pT2")
TC2 = subset(read.table("../mapping/TC2.ins"), V2 != "pT2")

allins = rbind(CA1, CA2, CT1, CT2, GA1, GA2, GT1, GT2, TC1, TC2)
gallins = GRanges(Rle(allins$V2), IRanges(start=allins$V4, width=1))

# Get expression data in mouse ES cells.
setwd("../data/")
exprs = read.delim("GSE93238_gene.fpkm.txt.gz")
# Remove stupid dot from expression data set (ENSMUSG00000000125.5)
exprs$gene = sub("\\..*", "", exprs$gene)
score = exprs$ES1 + exprs$ES2
cutoff = 0
mean(score > cutoff)
higenes = exprs$gene[score > cutoff]
logenes = exprs$gene[score <= cutoff]
head(higenes)
head(logenes)

# Get the genes and their coordinates.
ENSG = read.delim("ENSG_mm9.txt")
ENSG = subset(ENSG, ENSG %in% exprs$gene)
head(ENSG)
gENSG = GRanges(Rle(ENSG$chrom),
            IRanges(start=ENSG$start, end=ENSG$end))

ENSG_hi = subset(ENSG, ENSG %in% higenes)
ENSG_lo = subset(ENSG, ENSG %in% logenes)
gENSG_hi = GRanges(Rle(ENSG_hi$chrom),
            IRanges(start=ENSG_hi$start, end=ENSG_hi$end))
gENSG_lo = GRanges(Rle(ENSG_lo$chrom),
            IRanges(start=ENSG_lo$start, end=ENSG_lo$end))

Olo = sum(countOverlaps(gallins, gENSG_lo))
Ohi = sum(countOverlaps(gallins, gENSG_hi))
Orest = nrow(allins) - Olo - Ohi

# Total size of mm9: 2780285931.
Elo = sum(sum(coverage(gENSG_lo)))
Ehi = sum(sum(coverage(gENSG_hi)))
Erest = 2780285931. - Elo - Ehi

obsv = c(Ohi, Olo, Orest)
expt = c(Ehi, Elo, Erest)

print(obsv)
print(expt / sum(expt) * sum(obsv))

COL=c("#0B2639", "#598078", "#E1DBC0")

setwd("../figures")
pdf("spie_genes.pdf", useDingbats=FALSE)
par(mar=c(0,0,0,0))
spiechart(expected=expt, observed=obsv, col=COL)
dev.off()
