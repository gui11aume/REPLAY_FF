#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict

from gzopen import gzopen

REF1 = 'GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC'
MUT = defaultdict(int)

def compare(seq, MUT):
   for i in range(len(seq)):
      if seq[i] != REF1[i]: MUT[(REF1[i],seq[i])] += 1

def main():
   for line in sys.stdin:
      seq = line.rstrip()
      if len(seq) != len(REF1): continue
      compare(seq, MUT)

if __name__ == '__main__':
   main()
   print MUT
