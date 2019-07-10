#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict

from gzopen import gzopen

REF1 = 'GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC'
MUT = defaultdict(int)

def compare(seq, MUT):
   for i in range(len(seq)):
      if seq[i] != REF1[i]:
         MUT[(REF1[i],seq[i])] += 1

def main(fname):
   with gzopen(fname) as f:
      for line in f:
         items = line.split()
         if int(items[3] < 10): continue
         if float(items[4] < .95): continue
         if len(items[2]) != len(REF1): continue
         compare(items[2], MUT)

if __name__ == '__main__':
   for fname in sys.argv[1:]: main(fname)
   print MUT
