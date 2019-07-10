#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from gzopen import gzopen

REF1 = 'GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC'
MUT = [0] * len(REF1)

def compare(seq1, seq2, MUT):
   for i in range(len(seq1)):
      if seq1[i] != seq2[i]: MUT[i] += 1

def main(fname):
   with gzopen(fname) as f:
      for line in f:
         items = line.split()
         if int(items[3] < 10): continue
         if float(items[4] < .95): continue
         if len(items[2]) != len(REF1): continue
         compare(items[2], REF1, MUT)

if __name__ == '__main__':
   for fname in sys.argv[1:]: main(fname)
   print MUT
