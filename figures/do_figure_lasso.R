#agg = read.delim("aggfeat.txt")
agg = read.delim("aggfeat2.txt")
feat = read.delim("chromfeat.txt")

library(GenomicRanges)
gagg = GRanges(Rle(agg$chrom), IRanges(agg$pos, width=1))
gfeat = GRanges(Rle(feat$chrom),
               IRanges(start=feat$start, end=feat$end))

ov = findOverlaps(gagg, gfeat)
idxbrcd = queryHits(ov)
idxfeat = subjectHits(ov)

dat = data.frame(agg[idxbrcd,c(8,5,6,7)], feat[idxfeat,4:ncol(feat)])
dat = subset(dat, complete.cases(dat))

# 0-1 encoding of the factor (and remove intercept)
MM. = model.matrix(~ MM, data=dat)[,-1]

# Replace in the data.frame.
dat = data.frame(x=dat$x, MM., dat[,-c(1,2)])

# Remove profile P160, whic has correlation 1 with P159 (SUZ12).
dat = dat[,-match("P160", names(dat))]

write.table(dat, file="chromdata.txt", sep="\t",
      quote=FALSE, row.names=FALSE)

dat = read.delim("chromdata.txt")

COL = c("#0072B2", "#CC79A7", "#009E73")

# Discretize, center and scale variables.
dat[,7:ncol(dat)] = apply(MARGIN=2, dat[,7:ncol(dat)], function(x)
                          0+ (x > quantile(x, 0.9)))
dat = scale(dat, center=TRUE, scale=TRUE)

y = as.matrix(dat[,1])
x = as.matrix(dat[,-1])

# Do lasso (need library).
library(glmnet)

fit = glmnet(x,y)

pdf("figure_lasso.pdf", height=5.5, width=5)
plot(fit,
     col=c(COL, rep("#00000020", 41), 1, rep("#00000020",100)),
     lwd=c(2,2,2,rep(1,41), 2, rep(1,100)),
     bty="n", xlim=c(0,1))
dev.off()
