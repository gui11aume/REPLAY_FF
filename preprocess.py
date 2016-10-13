#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import seeq

from itertools import izip
from gzopen import gzopen

class AberrantReadException(Exception):
   pass

class Extractor:
   seq_after_tag = None
   seq_before_variant = None

   def extract_tag_and_variant(self, txt):
      '''Both reads have the same structure, with a tag (either a
      barcode or a UMI) immediately after the Illumina sequencing
      primer, and the variant towards the end of the read.'''

      # First extract the prefix and the suffix
      prefix = self.seq_after_tag.matchPrefix(txt, False)
      if not prefix:
         raise AberrantReadException

      # The first character of the suffix is the variant.
      suffix = self.seq_before_variant.matchSuffix(txt, False)
      if not suffix:
         raise AberrantReadException

      # The prefix is the tag, the first character
      # of the suffix is the variant.
      return (prefix, suffix[0])


class Read1Extractor(Extractor):
   def __init__(self):
      self.seq_after_tag = seeq.compile('CGCTAATTAATGGAATCATG', 3)
      self.seq_before_variant = seeq.compile('CGCTACGAGGCCGGCCGC', 3)


class Read2Extractor(Extractor):
   def __init__(self):
      self.seq_after_tag  = seeq.compile('TGCAACGAATTCATTAG', 3)
      self.seq_before_variant = seeq.compile('CACCTTGAAGTCGCCGATCA', 3)


def main(f, g):
   '''Top-level function to pre-process paired fastq files.'''

   Ex1 = Read1Extractor()
   Ex2 = Read2Extractor()

   linenumber = 0
   naberrant = 0
   for (read1,read2) in izip(f,g):
      linenumber = linenumber + 1
      if linenumber % 4 == 2:
         # Reading sequence line.
         try:
            BCD,SNP1 = Ex1.extract_tag_and_variant(read1)
            UMI,SNP2 = Ex2.extract_tag_and_variant(read2)
            print BCD, UMI, SNP1, SNP2
         except AberrantReadException:
            naberrant += 1
            continue

if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f, gzopen(sys.argv[2]) as g:
      main(f, g)
