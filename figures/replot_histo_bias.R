library(ggplot2)
dat = read.table("features.txt")
# Aggregate the data.

dattest = subset(dat, V6 == "test")
dattest$GChi = dattest$V9 > quantile(dattest$V9, 0.25, na.rm=TRUE)
head(dattest)

# Mean-aggregate scores per barcode, also return the 
# length, which will be useful for computing the
# standard deviation.
aggmean = aggregate(dattest$V2, FUN=mean,
   by=list(brcd=dattest$V1, MM=dattest$V3, GC=dattest$GChi))
agglen = aggregate(dattest$V2, FUN=length,
   by=list(brcd=dattest$V1, MM=dattest$V3, GC=dattest$GChi))

testagg = aggregate(aggmean$x, FUN=mean,
         by=list(MM=aggmean$MM, GC=aggmean$GC))

testsd = aggregate(agglen$x, FUN=function(x)
                   sqrt(sum(1/x)/4) / length(x),
         by=list(MM=aggmean$MM, GC=aggmean$GC))

print(testagg)
print(testsd)

limits = aes(
   ymax = testagg$x + 1.96*testsd$x,
   ymin = testagg$x - 1.96*testsd$x
)


COL = c("#E69F00", "#0072B2", "#CC79A7", "#009E73")

pdf("tmp.pdf", height=5, width=8, useDingbats=FALSE)
ggplot(testagg, aes(fill=MM, y=x, x=MM, group=rep)) +
   geom_bar(width=.8, position=position_dodge(.85), stat="identity") +
   scale_x_discrete(limits = c("CT", "CA", "GA", "GT")) +
   scale_fill_manual(values=COL) +
   geom_errorbar(limits, position=position_dodge(.85), width=0.2) +
   labs(x="Mismatch", y="Fraction A+T") +
   theme_bw() +
   theme(legend.position="none")
dev.off()
