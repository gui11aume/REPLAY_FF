library(ggplot2)

sem = function(x) { sd(x,na.rm=TRUE)/sqrt(length(x)) }
dat = read.table("moar_results.txt")

dattestGChi = subset(dat, V6 == "test" & V9 >= .43)
dattestGClo = subset(dat, V6 == "test" & V9 < .43)

# Mean-aggregate scores per barcode.
testGChi = aggregate(dattestGChi$V2, FUN=mean,
   by=list(brcd=dattestGChi$V1, MM=dattestGChi$V3,
           rep=dattestGChi$V7, GC=dattestGChi$V9))

testGClo = aggregate(dattestGClo$V2, FUN=mean,
   by=list(brcd=dattestGClo$V1, MM=dattestGClo$V3,
           rep=dattestGClo$V7, GC=dattestGClo$V9))

testaggGChi = aggregate(testGChi$x, FUN=mean, na.rm=TRUE,
         by=list(MM=testGChi$MM, rep=testGChi$rep))
testaggGClo = aggregate(testGClo$x, FUN=mean, na.rm=TRUE,
         by=list(MM=testGClo$MM, rep=testGClo$rep))
testaggGClo$rep = testaggGChi$rep+2

testsdGChi = aggregate(testGChi$x, FUN=sem,
         by=list(MM=testGChi$MM, rep=testGChi$rep))
testsdGClo = aggregate(testGClo$x, FUN=sem,
         by=list(MM=testGClo$MM, rep=testGClo$rep))

testagg = rbind(testaggGClo, testaggGChi)
testagg$rep = as.factor(testagg$rep)

print(testagg)

limits = aes(
   ymax = c(testaggGClo$x + 1.96*testsdGClo$x,
            testaggGChi$x + 1.96*testsdGChi$x),
   ymin = c(testaggGClo$x - 1.96*testsdGClo$x,
            testaggGChi$x - 1.96*testsdGChi$x),
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
