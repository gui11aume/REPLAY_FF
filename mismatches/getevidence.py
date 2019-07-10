#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen

COUNTER = defaultdict(lambda: defaultdict(lambda: defaultdict(int))) 
REF1 = 'GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC'

INDUCED = { 'indel':0, 'nndel':0 }
NNDUCED = { 'indel':0, 'nndel':0 }

def main(f):
   for line in f:
      items = line.split()
      if int(items[3]) < 10: continue
      brcd = items[0]
      mm   = items[1]
      seq  = items[2]
      COUNTER[brcd][mm][seq] += 1


def analyze():
   global INDUCED
   global NNDUCED
   for brcd,mm_counter in COUNTER.items():
      # COUNTER has been cleaned.
      seqset = set()
      for mm,seq_counter in mm_counter.items():
         # Get most frequent sequence for each MM code.
         seq = max(seq_counter, key=seq_counter.get)
         seqset.add(seq)
      if len(seqset) > 1:
         # Evidence for mutation induced by repair.
         seq1 = seqset.pop()
         seq2 = seqset.pop()
         if len(seq1) == len(seq2) == len(REF1): INDUCED['nndel'] += 1
         else:                                   INDUCED['indel'] += 1
      else:
         seq = seqset.pop()
         if seq != REF1:
            if len(seq) == len(REF1): NNDUCED['nndel'] += 1
            else:                     NNDUCED['indel'] += 1


def clean():
   for brcd in COUNTER.keys():
      mm_counter = COUNTER[brcd]
      # Remove all MM codes with less than 5 UMIs total.
      for mm in mm_counter.keys():
         if sum(mm_counter[mm].values()) < 5: mm_counter.pop(mm)
      # Only one MM code left (not interesting).
      if len(COUNTER[brcd].values()) < 2: COUNTER.pop(brcd)


if __name__ == '__main__':
   for fname in sys.argv[1:]:
      with gzopen(fname) as f: main(f)
   clean()
   analyze()
   print len(COUNTER)
   print INDUCED
   print NNDUCED
