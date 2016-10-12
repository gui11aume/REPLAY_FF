all: 24_1 24_2 24_3 24_4 48_1 48_2 48_3 48_4 con_1 con_2 con_3 con_4

24_1:
	python preprocess.py /data/24_13073_GCCAAG_read1.fastq.gz /data/24_13073_GCCAAG_read2.fastq.gz  > 24_1

24_2:
	python preprocess.py /data/24_13073_GTTCAC_read1.fastq.gz /data/24_13073_GTTCAC_read2.fastq.gz > 24_2

24_3:
	python preprocess.py /data/24_13073_TCTTAG_read1.fastq.gz /data/24_13073_TCTTAG_read2.fastq.gz > 24_3

24_4:
	python preprocess.py /data/24_13073_TGCAGT_read1.fastq.gz /data/24_13073_TGCAGT_read2.fastq.gz > 24_4

48_1:
	python preprocess.py /data/48_13074_ATGTCT_read1.fastq.gz /data/48_13074_ATGTCT_read2.fastq.gz > 48_1

48_2:
	python preprocess.py /data/48_13074_GCTGAA_read1.fastq.gz /data/48_13074_GCTGAA_read2.fastq.gz > 48_2

48_3:
	python preprocess.py /data/48_13074_TAGGAC_read1.fastq.gz /data/48_13074_TAGGAC_read2.fastq.gz > 48_3

48_4:
	python preprocess.py /data/48_13074_TCTACA_read1.fastq.gz /data/48_13074_TCTACA_read2.fastq.gz > 48_4

con_1:
	python preprocess.py /data/con_13072_AGCTTG_read1.fastq.gz /data/con_13072_AGCTTG_read2.fastq.gz > con_1

con_2:
	python preprocess.py /data/con_13072_CGATGT_read1.fastq.gz /data/con_13072_CGATGT_read2.fastq.gz > con_2

con_3:
	python preprocess.py /data/con_13072_CTGGAT_read1.fastq.gz /data/con_13072_CTGGAT_read2.fastq.gz > con_3

con_4:
	python preprocess.py /data/con_13072_TTAGGC_read1.fastq.gz /data/con_13072_TTAGGC_read2.fastq.gz > con_4
