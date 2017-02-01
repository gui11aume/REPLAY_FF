library(GenomicRanges)
library(ggplot2)
library(HMMt)

source("../helpers/domainify.R")

sem = function(x) { sd(x,na.rm=TRUE)/sqrt(length(x)) }

RT = read.delim("../data/RT_D3_EBM6_All.txt", skip=16)[,3:6]
# Remove outliers.
RT$Data_Value[RT$Data_Value < -2.5 | RT$Data_Value > 2.5] = NA

b = bridge(RT)
mod = BaumWelchT(b$x, b$series.length)

state = mod@ViterbiPath[b$nonvirtuals]
RT$Data_Value = state == 1
RTdom = domainify(RT)

gRT = GRanges(Rle(RTdom$block),
   IRanges(start=RTdom$start, end=RTdom$end))

dat = read.table("moar_results.txt")
gdat = GRanges(Rle(dat$V10), IRanges(start=dat$V12, width=1))

replicates_late = countOverlaps(gdat, gRT) > 0

dattest_late = subset(dat, V6 == "test" & replicates_late)
dattest_early = subset(dat, V6 == "test" & !replicates_late)

# Mean-aggregate scores per barcode.
test_late = aggregate(dattest_late$V2 == "AT", FUN=mean,
   by=list(brcd=dattest_late$V1, MM=dattest_late$V3,
           rep=dattest_late$V7, GC=dattest_late$V9))

test_early = aggregate(dattest_early$V2 == "AT", FUN=mean,
   by=list(brcd=dattest_early$V1, MM=dattest_early$V3,
           rep=dattest_early$V7, GC=dattest_early$V9))

testagg_late = aggregate(test_late$x, FUN=mean, na.rm=TRUE,
         by=list(MM=test_late$MM, rep=test_late$rep))
testagg_early = aggregate(test_early$x, FUN=mean, na.rm=TRUE,
         by=list(MM=test_early$MM, rep=test_early$rep))
testagg_early$rep = testagg_late$rep+2

testsd_late = aggregate(test_late$x, FUN=sem,
         by=list(MM=test_late$MM, rep=test_late$rep))
testsd_early = aggregate(test_early$x, FUN=sem,
         by=list(MM=test_early$MM, rep=test_early$rep))

testagg = rbind(testagg_early, testagg_late)
testagg$rep = as.factor(testagg$rep)


limits = aes(
   ymax = c(testagg_early$x + 1.96*testsd_early$x,
            testagg_late$x + 1.96*testsd_late$x),
   ymin = c(testagg_early$x - 1.96*testsd_early$x,
            testagg_late$x - 1.96*testsd_late$x),
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
