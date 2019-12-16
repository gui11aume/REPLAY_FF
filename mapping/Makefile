INSERTIONS = \
	CT1.ins CA1.ins GT1.ins GA1.ins \
	CT2.ins CA2.ins GT2.ins GA2.ins \
	TC1.ins TC2.ins

TARGETS= $(INSERTIONS) black.lst

CT1= CT_ACG
CT2= CT_TCT
CA1= CA_AAC CA_TTC
CA2= CA_AGG CA_ACA 
GT1= GT_ATC GT_TGT
GT2= GT_TGA GT_ACT
GA1= GA_ATT GA_TAA
GA2= GA_TGC GA_CCG
TC1= TC_ACT TC_AAC
TC2= TC_TTC TC_CCG

# Note that we add a % sign here. This will generate
# a pattern rule with multiple targets, which is the
# only way to tell Make that all fasta files are generated
# by the same command.
FASTA= $(addsuffix %fasta, $(CA1) $(GT1) $(GA1) $(CT1) \
       $(CA2) $(GT2) $(GA2) $(CT2) $(TC1) $(TC2))

all: $(TARGETS)

# This creates all fasta files (see definition of 'FASTA').
$(FASTA):
	python makefasta_nextseq.py /data/iPCR1_read1.fastq.gz \
		/data/iPCR1_read2.fastq.gz
	python makefasta_hiseq.py /data/iPCR2_read1.fastq.gz \
		/data/iPCR2_read2.fastq.gz CA_GT_GA
	python makefasta_hiseq.py /data/iPCR3_read1.fastq.gz \
		/data/iPCR3_read2.fastq.gz TC
#	python makefasta_nextseq.py /data/GF2-iPCR-pT2FF_R1.fastq.gz \
#		/data/GF2-iPCR-pT2FF_R2.fastq.gz
#	python makefasta_hiseq.py /data/MM_M_14930_read1.fastq.gz \
#		/data/MM_M_14930_read2.fastq.gz CA_GT_GA
#	python makefasta_hiseq.py /data/lane6_lane6_read1.fastq.gz \
#		/data/lane6_lane6_read2.fastq.gz TC

#%.map: %.fasta
#	gem-mapper -I /data/mm10_pT2_unmasked.gem \
#		-q ignore --unique-mapping -T4 -i $< > $@
%.map.gz: %.fasta
	gem-mapper -I /data/mm9_pT2_unmasked.gem \
		-q ignore --unique-mapping -T4 -i $< | gzip > $@


.SECONDARY:
%.stc: %.map.gz
	zcat $< | cut -f1 | starcode -t4 -d2 --print-clusters -o $@

# The secondary expansions will allow to replace the matching
# stem of the rule by the variables declared above.
.SECONDEXPANSION:
%.ins: $$(addsuffix .map.gz, $$($$*)) $$(addsuffix .stc, $$($$*))
	python collect.py $^ > $@

black.lst: $(INSERTIONS)
	python make_blacklist.py $^ > $@

cleantmp:
	rm -f *.fasta *.map.gz *.stc

clean:
	rm -f *fasta *.map.gz *.stc *.ins *.bcd
