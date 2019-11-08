#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

FASTASEQ = "/data/mm9_pT2.fasta"


def compute_GC(seqname, seq):
   for s in range(0,len(seq),5000):
      subseq = seq[s:(s+5000)]
      try:
         GC = (subseq.count("G") + subseq.count("C")) / \
            float(len(subseq) - subseq.count("N"))
      except ZeroDivisionError:
         GC = 'NA'
      print seqname, s+1, s+5000, GC


def read_genome_and_compute_GC(f):
   '''Read a fasta file and return a dictionary whose keys are the
   sequence names and values are the sequences in text format.
   Remove pT2 and weird chromosomes.'''

   txt = f.read()
   segments = txt.split('>')
   for segment in segments:
      if not segment: continue
      (header,seq) = segment.split('\n', 1)
      name = re.sub(r'\s.*', '', header)
      # Remove "chrUn_GL456385' etc. and pT2
      if '_' in name or 'pT2' in name: continue
      chromseq = seq.replace('\n', '')
      compute_GC(name, chromseq)

if __name__ == '__main__':
   with open(FASTASEQ) as f:
      read_genome_and_compute_GC(f)
