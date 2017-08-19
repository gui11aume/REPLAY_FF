library(ggplot2)
dat = read.table("moar_results.txt")
# Aggregate the data.

dattest = subset(dat, V6 == "test")
head(dattest)

# Mean-aggregate scores per barcode.
test = aggregate(dattest$V2, FUN=mean,
   by=list(brcd=dattest$V1, MM=dattest$V3,
           rep=as.factor(dattest$V7), GC=dattest$V9))

testagg = aggregate(test$x, FUN=mean, na.rm=TRUE,
         by=list(MM=test$MM, rep=test$rep))

testsd = aggregate(test$x, FUN=function(x)
                   sd(x,na.rm=TRUE)/sqrt(length(x)),
         by=list(MM=test$MM, rep=test$rep))

print(testagg)

limits = aes(
   ymax = testagg$x + 1.96*testsd$x,
   ymin = testagg$x - 1.96*testsd$x
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
