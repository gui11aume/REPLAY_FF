#!/usr/bin/env python
# -*- coding:utf-8 -*-

import seeq
import string
import sys

from itertools import izip

from gzopen import gzopen

comp = string.maketrans('gatcGATC', 'ctagCTAG')

# Forward stuff
anchorward = seeq.compile('CGCTAATTAATGGAATCATG', 3)
beforward = seeq.compile('CGCTACGAGGCCGGCCGC', 3)

# Reverse stuff
anchorev = seeq.compile('TGCAACGAATTCATTAG', 3)
beforev = seeq.compile('CACCTTGAAGTCGCCGATCA', 3)

def revcomp(seq):
   '''Reverse complement a DNA string.'''
   return seq[::-1].translate(comp)


def main(f, g):
   '''Top-level function to pre-process paired fastq files.'''

   linenumber = 0
   for (read1,read2) in izip(f,g):
      linenumber = linenumber + 1
      if linenumber % 4 == 2:
         # Reading sequence line.
         try:
            brcd = anchorward.matchPrefix(read1, False)
            umi = anchorev.matchPrefix(read2, False)
            SNP1 = beforward.matchSuffix(read1, False)[0]
            SNP2 = beforev.matchSuffix(read2, False)[0]
            print brcd, umi, SNP1, SNP2
         except TypeError:
            continue

if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f, gzopen(sys.argv[2]) as g:
      main(f, g)
