#!/usr/bin/env python
# -*- coding:utf-8 -*-

# Standard library packages.
import sys

from itertools import izip

# Others.
import seeq

from gzopen import gzopen

def trimSuffix(matcher, txt):
   return matcher.matchPrefix(txt, False) or ''

def trimPrefix(matcher, txt):
   return matcher.matchSuffix(txt, False) or ''

########  Mapping Pipeline ###############################################

def extract_reads_from_PE_fastq(fname_iPCR_PE1, fname_iPCR_PE2):
   """This function takes the 2 pair-end sequencing files and extracts the
   barcode making sure that the other read contains the transposon."""

   # Those are the scarcodes that allow to identify which
   # experiment is sequenced (CA, GA or GT mismatch).
   matchers = {
      'CA': seeq.compile('GCTAGCAGTCAGGAATCATG', 3),
      'GA': seeq.compile('GCTAGCTCGTTGGAATCATG', 3),
      'GT': seeq.compile('GCTAGCTCCGCAGAATCATG', 3),
   }
   
   indexes = {
      'CA': frozenset(['AAC', 'ACA', 'AGG', 'TTC']),
      'GA': frozenset(['ATT', 'CCG', 'TAA', 'TGC']),
      'GT': frozenset(['ACT', 'ATC', 'TGA', 'TGT']),
   }

   # Assign all valid triplets to a single fasta file for
   # the CT mismatch. Other files can be properly demultiplexed.
   outfiles = {
      'AAC': open('CA_AAC.fasta', 'w'),
      'ACA': open('CA_ACA.fasta', 'w'),
      'AGG': open('CA_AGG.fasta', 'w'),
      'TTC': open('CA_TTC.fasta', 'w'),

      'ACT': open('GT_ACT.fasta', 'w'),
      'ATC': open('GT_ATC.fasta', 'w'),
      'TGA': open('GT_TGA.fasta', 'w'),
      'TGT': open('GT_TGT.fasta', 'w'),

      'ATT': open('GA_ATT.fasta', 'w'),
      'CCG': open('GA_CCG.fasta', 'w'),
      'TAA': open('GA_TAA.fasta', 'w'),
      'TGC': open('GA_TGC.fasta', 'w'),
   }

   # End of the pT2 transposon sequence.
   pT2 = seeq.compile('AAACTTCCGACTTCAACTGTA', 3)

   with gzopen(fname_iPCR_PE1) as f, gzopen(fname_iPCR_PE2) as g:
      # Aggregate iterator of f,g iterators -> izip(f,g).
      for lineno,(line1,line2) in enumerate(izip(f,g)):
         # Take sequence lines of the fastq file.
         if lineno % 4 != 1: continue

         # Use the scarcode to identify the experiment.
         for MM,matcher in matchers.items():
            brcd = trimSuffix(matcher, line1)
            # If we find a barcode between 13 and 25 nucleotides
            # then the scarcode must have been the right one.
            if len(brcd) < 13 or len(brcd) > 25: continue

            # Find pT2 on the reverse read. Abort if we cannot or
            # if the remaining sequence is unmappable (too short).
            genome = trimPrefix(pT2, line2.rstrip()).split('CATG')[0]
            if len(genome) < 18: continue

            # The first 3 nucleotides of the reverse read are the
            # index. Check that it belongs to the right group.
            idx = line2[:3]
            if idx in indexes[MM]:
               outf = outfiles[idx]
               outf.write('>%s\n%s\n' % (brcd,genome))

            # If the script reaches this point, the mismatch was
            # already identified (even though nothing may be
            # printed) so there is no need to check the others.
            break

if __name__ == '__main__':
   extract_reads_from_PE_fastq(sys.argv[1], sys.argv[2])
