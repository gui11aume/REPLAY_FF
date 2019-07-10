#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen

COUNTER = defaultdict(lambda: defaultdict(int))
MUT = defaultdict(int)

REF1 = 'GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC'

def main(f):
   for line in f:
      items = line.split()
      # Skip UMIs supported by less than 10 reads.
      if int(items[3]) < 10: continue
      brcd = items[0]
      seq  = items[2]
      COUNTER[brcd][seq] += 1


def getmut(seq):
   global MUT
   if len(seq) != len(REF1):
      MUT['indel'] += 1
      return
   for i in range(len(seq)):
      if seq[i] != REF1[i]: MUT[(REF1[i],seq[i])] += 1


def analyze():
   global COUNTER
   for brcd,seq_counter in COUNTER.items():
      for seq in seq_counter:
         if seq_counter[seq] == 1: getmut(seq)


if __name__ == '__main__':
   for fname in sys.argv[1:]:
      with gzopen(fname) as f: main(f)
   analyze()
   print MUT
