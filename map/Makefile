all: AAC_CA_.fasta ACT_GT_.fasta ATC_GT_.fasta CCG_GA_.fasta \
	TGA_GT_.fasta TGT_GT_.fasta ACA_CA_.fasta AGG_CA_.fasta \
	ATT_GA_.fasta TAA_GA_.fasta TGC_GA_.fasta TTC_CA_.fasta \
	AAC_CA_.map ACT_GT_.map ATC_GT_.map CCG_GA_.map TGA_GT_.map \
	TGT_GT_.map ACA_CA_.map AGG_CA_.map ATT_GA_.map TAA_GA_.map \
	TGC_GA_.map TTC_CA_.map \
	AAC_CA_.stc ACT_GT_.stc ATC_GT_.stc CCG_GA_.stc TGA_GT_.stc \
	TGT_GT_.stc ACA_CA_.stc AGG_CA_.stc ATT_GA_.stc TAA_GA_.stc \
	TGC_GA_.stc TTC_CA_.stc \
	AAC_CA_.no ACT_GT_.no ATC_GT_.no CCG_GA_.no TGA_GT_.no \
	TGT_GT_.no ACA_CA_.no AGG_CA_.no ATT_GA_.no TAA_GA_.no \
	TGC_GA_.no TTC_CA_.no

%.fasta:
	#./tripeline /data/MM_M_14930_read1.fastq.gz \
	#	/data/MM_M_14930_read2.fastq.gz
	echo "hello"

%.map: %.fasta
	#/GEM/bin/gem-mapper -I /data/mm10_pT2_unmasked.gem \
	#	-q ignore --unique-mapping -T4 -i $< > $@
	echo $<

%.stc: %.map
	cut -f1 $< | starcode -t4 -d2 --print-clusters -o $@

%.no: %.map
	./collect $*

clean:
	rm out *fasta
