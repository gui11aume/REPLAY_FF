library(GenomicRanges)
library(RCircos)

# Overwrite RCircos default fuctions.
source("improvedRCircos.R")

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

CA = rbind(CA1, CA2)
CT = rbind(CT1, CT2)
GA = rbind(GA1, GA2)
GT = rbind(GT1, GT2)
TC = rbind(TC1, TC2)

gCA = GRanges(Rle(CA$V2), IRanges(start=CA$V4, width=1))
gCT = GRanges(Rle(CT$V2), IRanges(start=CT$V4, width=1))
gGA = GRanges(Rle(GA$V2), IRanges(start=GA$V4, width=1))
gGT = GRanges(Rle(GT$V2), IRanges(start=GT$V4, width=1))
gTC = GRanges(Rle(TC$V2), IRanges(start=TC$V4, width=1))

data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
chr.exclude = NULL
cyto.info = UCSC.Mouse.GRCm38.CytoBandIdeogram
tracks.inside = 5
tracks.outside = 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,
      tracks.inside, tracks.outside)

# Keep only insertions mapping to annotations.
gGRCm38 = makeGRangesFromDataFrame(UCSC.Mouse.GRCm38.CytoBandIdeogram,
   seqnames.field="Chromosome", start.field="ChromStart",
   end.field="ChromEnd")

CA = subset(CA, countOverlaps(gCA, gGRCm38) > 0)
CT = subset(CT, countOverlaps(gCT, gGRCm38) > 0)
GA = subset(GA, countOverlaps(gGA, gGRCm38) > 0)
GT = subset(GT, countOverlaps(gGT, gGRCm38) > 0)
TC = subset(TC, countOverlaps(gTC, gGRCm38) > 0)

print(nrow(CA))
print(nrow(CT))
print(nrow(GA))
print(nrow(GT))
print(nrow(TC))


# Take 1000 at random.
trackCA = CA[sample(nrow(CA),1000),c(2,4,4)]
trackCT = CT[sample(nrow(CT),1000),c(2,4,4)]
trackGA = GA[sample(nrow(GA),1000),c(2,4,4)]
trackGT = GT[sample(nrow(GT),1000),c(2,4,4)]
trackTC = TC[sample(nrow(TC),1000),c(2,4,4)]


rcircos.params = RCircos.Get.Plot.Parameters()
rcircos.params$track.height = 0.08
#rcircos.params$track.in.start = 1.05
rcircos.params$track.in.start = 1.2
rcircos.params$chrom.paddings = 500

RCircos.Reset.Plot.Parameters(rcircos.params)

# Previous palettes.
#COL = c("#81437410", "#51A39D10", "#B7695C10", "#CDBB7910")
# A color blind palette
# (taken from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
# cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73",
#  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#COL = c("#E69F0010", "#0072B210", "#CC79A710", "#009E7306")

# This one is taken from the web site below
# https://www.schemecolor.com/5-color-rainbow-color-scheme.php
COL = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")

setwd("../figures")
pdf("integ_circos.pdf", useDingbats=F)
par(mar=c(0,0,0,0))
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
Tile.Plot(trackCA, 1, "in", col=COL[1])
Tile.Plot(trackCT, 2, "in", col=COL[2])
Tile.Plot(trackGA, 3, "in", col=COL[3])
Tile.Plot(trackGT, 4, "in", col=COL[4])
Tile.Plot(trackTC, 5, "in", col=COL[5])
dev.off()
