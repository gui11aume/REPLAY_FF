REF1 = "GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC"
REF2 = "TGCAACGAATTCATTAGTGCGGATGATCTTGTCGGTGAAGATCACGCTGTCCTCGGGGAAGCCGGTGCCCACCACCTTGAAGTCGCCGATCA"
OK = total = 0
for (fname in dir(patt=".*_seq.txt.gz")) {
   #dat = subset(read.table(fname, as.is=TRUE), V3 > 10 & nchar(V2) == nchar(REF1))
   #dat = subset(read.table(fname, as.is=TRUE), V5 > 10 & nchar(V4) == nchar(REF2))
   #dat = subset(read.table(fname, as.is=TRUE), V7 > 10 & !grepl("N", V6) & V8 > .95 & nchar(V6) != nchar(REF2))
   #dat = subset(read.table(fname, as.is=TRUE), V4 > 10 & !grepl("N", V3) & V4 > .95)
   dat = subset(read.table(fname, as.is=TRUE), V7 > 10 & !grepl("N", V6) & V8 > .95)
   #dat = subset(read.table(fname, as.is=TRUE), V7 > 10 & !grepl("N", V6) & V8 > .95 & nchar(V6) == 93)
   #print(head((sample(dat$V6)), 10))
   #res = tapply(X=dat$V6, INDEX=dat$V2, FUN=function(x) mean(nchar(x) == nchar(REF2)))
   #res = tapply(X=dat$V6, INDEX=dat$V2, FUN=function(x) { y = table(nchar(x)); y/sum(y) })
   res = tapply(X=dat$V6, INDEX=dat$V2, FUN=function(x) mean(x == REF2))
   print(fname)
   #for (i in 1:25) {
   #   if (dat$V6[i] != REF2) {
   #      print(paste(dat$V2[i], dat$V6[i], sep="_"))
   #   }
   #}
   print(res)
   #OK = OK + (res["FF"] > res["GC"] && res["FF"] > res["AT"])
   #total = total + 1
}
print(c(OK, total))
