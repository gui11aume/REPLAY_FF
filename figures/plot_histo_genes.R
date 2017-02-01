library(GenomicRanges)
library(ggplot2)

sem = function(x) { sd(x,na.rm=TRUE)/sqrt(length(x)) }

ensembl = subset(read.table("../data/GRCm38_genes.gtf", sep="\t"),
            grepl("ensembl", V2))
gid = sub(".*(ENSMUSG\\d{11}).*", "\\1", ensembl$V9)
gensembl = GRanges(Rle(paste("chr", ensembl$V1, sep="")),
            IRanges(start=ensembl$V4, end=ensembl$V5))

# Extract expressed genes.
exprs = read.delim("../data/GSE93238_gene.fpkm.txt.gz")
cutoff = median(exprs$ES2, na.rm=TRUE)
higenes = exprs$gene[exprs$ES2 > cutoff]

# Removed stupid dot from expression data set (ENSMUSG00000000125.5)
higenes = sub("\\..*", "", higenes)

ensembl_hi = subset(ensembl, gid %in% higenes)
gensembl_hi = GRanges(Rle(paste("chr", ensembl_hi$V1, sep="")),
            IRanges(start=ensembl_hi$V4, end=ensembl_hi$V5))

dat = read.table("moar_results.txt")
gdat = GRanges(Rle(dat$V10), IRanges(start=dat$V12, width=1))

is_in_a_gene = countOverlaps(gdat, gensembl) > 0
is_in_a_hi_gene = countOverlaps(gdat, gensembl_hi) > 0

dattest_in = subset(dat, V6 == "test" & is_in_a_hi_gene)
dattest_out = subset(dat, V6 == "test" & !is_in_a_hi_gene)

# Mean-aggregate scores per barcode.
test_in = aggregate(dattest_in$V2 == "AT", FUN=mean,
   by=list(brcd=dattest_in$V1, MM=dattest_in$V3,
           rep=dattest_in$V7, GC=dattest_in$V9))

test_out = aggregate(dattest_out$V2 == "AT", FUN=mean,
   by=list(brcd=dattest_out$V1, MM=dattest_out$V3,
           rep=dattest_out$V7, GC=dattest_out$V9))

testagg_in = aggregate(test_in$x, FUN=mean, na.rm=TRUE,
         by=list(MM=test_in$MM, rep=test_in$rep))
testagg_out = aggregate(test_out$x, FUN=mean, na.rm=TRUE,
         by=list(MM=test_out$MM, rep=test_out$rep))
testagg_out$rep = testagg_in$rep+2

testsd_in = aggregate(test_in$x, FUN=sem,
         by=list(MM=test_in$MM, rep=test_in$rep))
testsd_out = aggregate(test_out$x, FUN=sem,
         by=list(MM=test_out$MM, rep=test_out$rep))

testagg = rbind(testagg_out, testagg_in)
testagg$rep = as.factor(testagg$rep)


limits = aes(
   ymax = c(testagg_out$x + 1.96*testsd_out$x,
            testagg_in$x + 1.96*testsd_in$x),
   ymin = c(testagg_out$x - 1.96*testsd_out$x,
            testagg_in$x - 1.96*testsd_in$x),
)

COL = c("#E69F00", "#0072B2", "#CC79A7", "#009E73")

pdf("tmp.pdf", height=5, width=8, useDingbats=FALSE)
ggplot(testagg, aes(fill=MM, y=x, x=MM, group=rep)) +
   geom_bar(width=.8, position=position_dodge(.85), stat="identity") +
   scale_x_discrete(limits = c("CT", "CA", "GA")) +
   scale_fill_manual(values=COL) +
   geom_errorbar(limits, position=position_dodge(.85), width=0.2) +
   labs(x="Mismatch", y="Fraction A+T") +
   theme_bw() +
   theme(legend.position="none")
dev.off()
