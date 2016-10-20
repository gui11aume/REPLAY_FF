TARGETS= CT.ins CA.ins GT.ins GA.ins

CT= CT
CA= CA_AAC CA_AGG CA_TTC CA_ACA 
GT= GT_ATC GT_TGA GT_ACT GT_TGT
GA= GA_ATT GA_TAA GA_TGC GA_CCG

# Note that we add a % sign here. This will generate
# a pattern rule with multiple targets, which is the
# only way to tell Make that all fasta files are generated
# by the same command.
FASTA= $(addsuffix %fasta, $(CA) $(GT) $(GA) $(CT))

all: $(TARGETS)

# This creates all fasta files (see definition of 'FASTA').
$(FASTA):
	python makefasta_nextseq.py /data/GF2-iPCR-pT2FF_R1.fastq.gz \
		/data/GF2-iPCR-pT2FF_R2.fastq.gz
	python makefasta_hiseq.py /data/MM_M_14930_read1.fastq.gz \
		/data/MM_M_14930_read2.fastq.gz

%.map: %.fasta
	gem-mapper -I /data/mm10_pT2_unmasked.gem \
		-q ignore --unique-mapping -T4 -i $< > $@

.SECONDARY:
%.stc: %.map
	cut -f1 $< | starcode -t4 -d2 --print-clusters -o $@

# The secondary expansions will allow to replace the matching
# stem of the rule by the variables declared above.
.SECONDEXPANSION:
%.ins: $$(addsuffix .map, $$($$*)) $$(addsuffix .stc, $$($$*))
	python collect.py $^ > $@

cleantmp:
	rm -f *.fasta *.map *.stc .fasta

clean:
	rm -f *fasta *.map *.stc *.ins
