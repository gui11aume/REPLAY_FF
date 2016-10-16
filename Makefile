vpath %_read1.fastq.gz /data
vpath %_read2.fastq.gz /data

TARGETS= \
GA_24_15173_CGATGT.co.gz  \
GA_24_15173_GTTCAC.co.gz  \
GA_24_15173_TCTTAG.co.gz  \
GA_24_15173_TGCAGT.co.gz  \
GA_24c_15171_AGCTTG.co.gz \
GA_24c_15171_GCCAAG.co.gz \
GA_24c_15171_TGACCA.co.gz \
GA_24c_15171_TTAGGC.co.gz \
GA_48_15174_CTTGTA.co.gz  \
GA_48_15174_GCTGAA.co.gz  \
GA_48_15174_GGCTAC.co.gz  \
GA_48_15174_TAGCTT.co.gz  \
GA_48c_15172_ACAGTG.co.gz \
GA_48c_15172_ACTTGA.co.gz \
GA_48c_15172_ATGTCT.co.gz \
GA_48c_15172_GATCAG.co.gz \
GA_48c_15172_TCTACA.co.gz \
GT_24_15432_CGATGT.co.gz  \
GT_24_15432_GTTCAC.co.gz  \
GT_24_15432_TCTTAG.co.gz  \
GT_24_15432_TGCAGT.co.gz  \
GT_24c_15430_AGCTTG.co.gz \
GT_24c_15430_GCCAAG.co.gz \
GT_24c_15430_TGACCA.co.gz \
GT_24c_15430_TTAGGC.co.gz \
GT_48_15433_CTTGTA.co.gz  \
GT_48_15433_GCTGAA.co.gz  \
GT_48_15433_GGCTAC.co.gz  \
GT_48_15433_TAGCTT.co.gz  \
GT_48c_15431_ACAGTG.co.gz \
GT_48c_15431_ACTTGA.co.gz \
GT_48c_15431_ATGTCT.co.gz \
GT_48c_15431_GATCAG.co.gz \
CA_24_15639_CGATGT.co.gz  \
CA_24_15639_GTTCAC.co.gz  \
CA_24_15639_GTTCAC.co.gz  \
CA_24_15639_TCTTAG.co.gz  \
CA_24_15639_TGCAGT.co.gz  \
CA_24c_15637_AGCTTG.co.gz \
CA_24c_15637_GCCAAG.co.gz \
CA_24c_15637_TGACCA.co.gz \
CA_24c_15637_TTAGGC.co.gz \
CA_48_15640_CTTGTA.co.gz  \
CA_48_15640_GCTGAA.co.gz  \
CA_48_15640_GGCTAC.co.gz  \
CA_48_15640_TAGCTT.co.gz  \
CA_48c_15638_ACAGTG.co.gz \
CA_48c_15638_ACTTGA.co.gz \
CA_48c_15638_ATGTCT.co.gz \
CA_48c_15638_GATCAG.co.gz \
24_13073_GCCAAG.co.gz     \
24_13073_GTTCAC.co.gz     \
24_13073_TCTTAG.co.gz     \
24_13073_TGCAGT.co.gz     \
48_13074_ATGTCT.co.gz     \
48_13074_GCTGAA.co.gz     \
48_13074_TAGGAC.co.gz     \
48_13074_TCTACA.co.gz     \
con_13072_AGCTTG.co.gz    \
con_13072_CGATGT.co.gz    \
con_13072_CTGGAT.co.gz    \
con_13072_TTAGGC.co.gz    \
con_13072_TTAGGC.co.gz    \


all: $(TARGETS)

%.pps: %_read1.fastq.gz %_read2.fastq.gz
	python preprocess.py $^ > $@

%.stc: %.pps
	cut -f1 $< | starcode -d2 --print-clusters > $@

%.co.gz: %.pps %.stc
	python count.py $^ | gzip > $@ && rm $^

clean:
	rm -rf *.pps *.co.gz
