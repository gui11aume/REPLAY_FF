all:	CA_AAC.no \
	CA_AGG.no \
	CA_TTC.no \
	CA_ACA.no \
	GT_ATC.no \
	GT_TGA.no \
	GT_ACT.no \
	GT_TGT.no \
	GA_ATT.no \
	GA_TAA.no \
	GA_TGC.no \
	GA_CCG.no \

%.fasta:
	./tripeline /data/MM_M_14930_read1.fastq.gz \
		/data/MM_M_14930_read2.fastq.gz

%.map: %.fasta
	gem-mapper -I /data/mm10_pT2_unmasked.gem \
		-q ignore --unique-mapping -T4 -i $< > $@

%.stc: %.map
	cut -f1 $< | starcode -t4 -d2 --print-clusters -o $@

%.no: %.map %.stc
	./collect $*

cleantmp:
	rm -f *.fasta *.map
