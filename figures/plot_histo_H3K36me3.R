library(ggplot2)
dat = read.delim("chromdata.txt")

# Reformat.
dat$H3K36me3 = dat[,45] > quantile(dat[,45], .9)
dat$MM = 1 + dat$MMCT + 2 * dat$MMGA + 3 * dat$MMGT
dat$MM = c("G:A", "G:T", "C:A", "C:T")[dat$MM]

# Aggregate the data.
aggmean = aggregate(dat$x, FUN=mean,
   by=list(MM=dat$MM, H3K36me3=dat$H3K36me3))

aggsd =  aggregate(dat$x, FUN=function(x) sd(x)/length(x),
   by=list(MM=dat$MM, H3K36me3=dat$H3K36me3))


print(aggmean)

limits = aes(
   ymax = aggmean$x + 1.96*aggsd$x,
   ymin = aggmean$x - 1.96*aggsd$x
)

COL = c("#CC79A7", "#009E73", "#E69F00", "#0072B2")

pdf("histo_H3K36me3.pdf", height=5, width=8, useDingbats=FALSE)
ggplot(aggmean, aes(fill=MM, y=x, x=MM, group=H3K36me3)) +
   geom_bar(width=.8, position=position_dodge(.85), stat="identity") +
   scale_x_discrete(limits = c("G:A", "G:T", "C:A", "C:T")) +
   scale_fill_manual(values=COL) +
   geom_errorbar(limits, position=position_dodge(.85), width=0.2) +
   labs(x="Mismatch", y="Fraction A+T") +
   theme_bw() +
   theme(legend.position="none")
dev.off()
