all:	AAC_CA.no \
	AGG_CA.no \
	TTC_CA.no \
	ACA_CA.no \
	ATC_GT.no \
	TGA_GT.no \
	ACT_GT.no \
	TGT_GT.no \
	ATT_GA.no \
	TAA_GA.no \
	TGC_GA.no \
	CCG_GA.no \

%.fasta:
	./tripeline /data/MM_M_14930_read1.fastq.gz \
		/data/MM_M_14930_read2.fastq.gz

%.map: %.fasta
	/GEM/bin/gem-mapper -I /data/mm10_pT2_unmasked.gem \
		-q ignore --unique-mapping -T4 -i $< > $@

%.stc: %.map
	cut -f1 $< | starcode -t4 -d2 --print-clusters -o $@

%.no: %.map
	./collect $*

cleantmp:
	rm -f *.fasta *.map
