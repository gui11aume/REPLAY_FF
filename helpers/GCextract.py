#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

FASTASEQ = "/data/mm10_pT2.fasta"

def read_genome(f):
   '''Read a fasta file and return a dictionary whose keys are the
   sequence names and values are the sequences in text format.
   Remove pT2 and weird chromosomes.'''

   genome = dict()
   txt = f.read()
   segments = txt.split('>')
   for segment in segments:
      if not segment: continue
      (header,seq) = segment.split('\n', 1)
      name = re.sub(r'\s.*', '', header)
      # Remove "chrUn_GL456385' etc. and pT2
      if '_' in name or 'pT2' in name: continue
      genome[name] = seq.replace('\n', '')
   return genome


def compute_GC(chrom, pos, genome):
   seq = genome[chrom]
   seq10kb = seq[(pos-5000):(pos+5000)].upper()
   seq1Mb = seq[(pos-500000):(pos+500000)].upper()

   GC10kb = (seq10kb.count("G") + seq10kb.count("C")) / \
      float(len(seq10kb) - seq10kb.count("N"))
   GC1Mb = (seq1Mb.count("G") + seq1Mb.count("C")) / \
      float(len(seq1Mb) - seq1Mb.count("N"))

   return (GC10kb,GC1Mb)


if __name__ == '__main__':
   with open(FASTASEQ) as f:
      genome = read_genome(f)
   for pos in range(0, 200000000, 1000000):
      try:
         (GC1,GC2) = compute_GC('chr1', pos, genome)
      except ZeroDivisionError:
         GC2 = 'NA'
      print pos, GC2
