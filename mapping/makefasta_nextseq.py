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

########  Mapping Pipeline ###############################################

def extract_reads_from_PE_fastq(fname_iPCR_PE1, fname_iPCR_PE2):
   """This function takes the 2 pair-end sequencing files and extracts the
   barcode making sure that the other read contains the transposon."""

   # Those are the scarcodes that allow to identify which
   # experiment is sequenced (CA, GA or GT mismatch).
   matcher = seeq.compile('CGCTAATTAATGGAATCATG', 3)

   outf = open('CT.fasta', 'w')

   with gzopen(fname_iPCR_PE1) as f, gzopen(fname_iPCR_PE2) as g:
      # Aggregate iterator of f,g iterators -> izip(f,g).
      for lineno,(line1,line2) in enumerate(izip(f,g)):
         # Take sequence lines of the fastq file.
         if lineno % 4 != 1: continue

         brcd = trimSuffix(matcher, line1)
         # If we find a barcode between 13 and 25 nucleotides
         # then the scarcode must have been the right one.
         if len(brcd) < 13 or len(brcd) > 25: continue

         # Remove first 25 nucleotides, split o "CATG" and take
         # the first fragment. If genome fragment is too short
         # for mapping then throw it away.
         genome = line2.rstrip()[25:].split('CATG')[0]
         if len(genome) < 18: continue

         outf.write('>%s\n%s\n' % (brcd,genome))

if __name__ == '__main__':
   extract_reads_from_PE_fastq(sys.argv[1], sys.argv[2])
