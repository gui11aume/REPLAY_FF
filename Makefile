vpath %_read1.fastq.gz /data
vpath %_read2.fastq.gz /data
#all: 24_1 24_2 24_3 24_4 48_1 48_2 48_3 48_4 con_1 con_2 con_3 con_424_1:
#	python preprocess.py /data/24_13073_GCCAAG_read1.fastq.gz /data/24_13073_GCCAAG_read2.fastq.gz  > 24_1
#
#24_2:
#	python preprocess.py /data/24_13073_GTTCAC_read1.fastq.gz /data/24_13073_GTTCAC_read2.fastq.gz > 24_2
#
#24_3:
#	python preprocess.py /data/24_13073_TCTTAG_read1.fastq.gz /data/24_13073_TCTTAG_read2.fastq.gz > 24_3
#
#24_4:
#	python preprocess.py /data/24_13073_TGCAGT_read1.fastq.gz /data/24_13073_TGCAGT_read2.fastq.gz > 24_4
#
#48_1:
#	python preprocess.py /data/48_13074_ATGTCT_read1.fastq.gz /data/48_13074_ATGTCT_read2.fastq.gz > 48_1
#
#48_2:
#	python preprocess.py /data/48_13074_GCTGAA_read1.fastq.gz /data/48_13074_GCTGAA_read2.fastq.gz > 48_2
#
#48_3:
#	python preprocess.py /data/48_13074_TAGGAC_read1.fastq.gz /data/48_13074_TAGGAC_read2.fastq.gz > 48_3
#
#48_4:
#	python preprocess.py /data/48_13074_TCTACA_read1.fastq.gz /data/48_13074_TCTACA_read2.fastq.gz > 48_4
#
#con_1:
#	python preprocess.py /data/con_13072_AGCTTG_read1.fastq.gz /data/con_13072_AGCTTG_read2.fastq.gz > con_1
#
#con_2:
#	python preprocess.py /data/con_13072_CGATGT_read1.fastq.gz /data/con_13072_CGATGT_read2.fastq.gz > con_2
#
#con_3:
#	python preprocess.py /data/con_13072_CTGGAT_read1.fastq.gz /data/con_13072_CTGGAT_read2.fastq.gz > con_3
#
#con_4:
#	python preprocess.py /data/con_13072_TTAGGC_read1.fastq.gz /data/con_13072_TTAGGC_read2.fastq.gz > con_4

TARGETS= \
GA_24_15173_CGATGT.co  \
GA_24_15173_GTTCAC.co  \
GA_24_15173_TCTTAG.co  \
GA_24_15173_TGCAGT.co  \
GA_24c_15171_AGCTTG.co \
GA_24c_15171_GCCAAG.co \
GA_24c_15171_TGACCA.co \
GA_24c_15171_TTAGGC.co \
GA_48_15174_CTTGTA.co  \
GA_48_15174_GCTGAA.co  \
GA_48_15174_GGCTAC.co  \
GA_48_15174_TAGCTT.co  \
GA_48c_15172_ACAGTG.co \
GA_48c_15172_ACTTGA.co \
GA_48c_15172_ATGTCT.co \
GA_48c_15172_GATCAG.co \
GA_48c_15172_TCTACA.co \
GT_24_15432_CGATGT.co  \
GT_24_15432_GTTCAC.co  \
GT_24_15432_TCTTAG.co  \
GT_24_15432_TGCAGT.co  \
GT_24c_15430_AGCTTG.co \
GT_24c_15430_GCCAAG.co \
GT_24c_15430_TGACCA.co \
GT_24c_15430_TTAGGC.co \
GT_48_15433_CTTGTA.co  \
GT_48_15433_GCTGAA.co  \
GT_48_15433_GGCTAC.co  \
GT_48_15433_TAGCTT.co  \
GT_48c_15431_ACAGTG.co \
GT_48c_15431_ACTTGA.co \
GT_48c_15431_ATGTCT.co \
GT_48c_15431_GATCAG.co 

all: $(TARGETS)

.INTERMEDIATE: %.pps
.INTERMEDIATE: %.stc

%.pps: %_read1.fastq.gz %_read2.fastq.gz
	python preprocess.py $^ > $@

%.stc: %.pps
	cut -f1 $< | starcode -d3 --print-clusters > $@

%.co: %.pps %.stc
	python count.py $^ > $@

clean:
	rm -rf *.pps *.co
