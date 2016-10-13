#!/usr/bin/env python
# -*- coding:utf-8 -*-

'''
Relevant information.

Scarcodes:
   CT: CGCTAATTAATG
   CA: GCTAGCAGTCAG
   GA: GCTAGCTCGTTG
   GT: GCTAGCTCCGCA
'''

import sys
import datetime

import seeq


from itertools import izip
from gzopen import gzopen

class AberrantReadException(Exception):
   pass


class LaneInfo:
   def __init__(self, fname1, fname2):
      self.fname1    = fname1
      self.fname2    = fname2
      self.ntotal    = 0
      self.naberrant = 0

   def write_to_file(self, f):
      # Time stamp.
      dt = datetime.datetime
      f.write(dt.strftime(dt.now(),
         '%Y-%m-%d %H:%M:%S Preprocessing summary\n'))

      # Processed files.
      f.write('%s\n%s\n' % (self.fname1, self.fname2))

      # Total and percent reads lost.
      lost = float(self.naberrant)
      f.write('Reads lost:\t%d (%.2f %%)\n' % \
            (lost, 100 * lost / self.ntotal))

      # End of report.
      f.write('---\n')


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
      self.seq_after_tag = seeq.compile('GAATCATGAACACCCGCAT', 3)
      self.seq_before_variant = seeq.compile('CGCTACGAGGCCGGCCGC', 3)


class Read2Extractor(Extractor):
   def __init__(self):
      self.seq_after_tag  = seeq.compile('TGCAACGAATTCATTAG', 3)
      self.seq_before_variant = seeq.compile('CACCTTGAAGTCGCCGATCA', 3)


def main(f, g, info):
   '''Top-level function to pre-process paired fastq files.'''

   Ex1 = Read1Extractor()
   Ex2 = Read2Extractor()

   linenumber = 0
   for (read1,read2) in izip(f,g):
      linenumber = linenumber + 1
      if linenumber % 4 == 2:
         info.ntotal += 1
         try:
            BCD,SNP1 = Ex1.extract_tag_and_variant(read1)
            UMI,SNP2 = Ex2.extract_tag_and_variant(read2)
            # Concatenate the barcode and the UMI into a single
            # tag. To mark the distrinction, insertion 8 a's in
            # between. Starcode will keep the lower case.
            sys.stdout.write('%saaaaaaaa%s\t%s\t%s\n' % \
                  (BCD, UMI, SNP1, SNP2))
         except AberrantReadException:
            info.naberrant += 1
            continue


if __name__ == '__main__':
   info = LaneInfo(sys.argv[1], sys.argv[2])
   try:
      main(gzopen(sys.argv[1]), gzopen(sys.argv[2]), info)
   finally:
      info.write_to_file(open('pps_logs.txt', 'a'))
