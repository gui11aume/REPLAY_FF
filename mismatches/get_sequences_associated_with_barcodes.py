#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from itertools import combinations

from gzopen import gzopen

SEQ = defaultdict(lambda: defaultdict(set))

REF1 = "GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC"

def check(statusdict):
   S = set([e for L in statusdict.values() for e in L])
   if len(S) > 1:
      return "multiple"
   seq = S.pop()
   if seq == REF1:
      return "single WT"
   return "single MUT"


def find_intersection(statusdict):
   for S1,S2 in combinations(statusdict.values(),2):
      if S1 & S2: return (S1 & S2).pop()
   return False


def main(fname):
   with gzopen(fname) as f:
      for line in f:
         items = line.split()
         if int(items[3]) < 10 or float(items[4]) < .95: continue
         if items[2] == REF1: continue
         brd = items[0]
         status = items[1]
         # Add non WT sequences.
         SEQ[brd][fname].add((status, items[2]))

if __name__ == '__main__':
   for fname in sys.argv[1:]:
      main(fname)
   for brd,statusdict in SEQ.items():
      if len(statusdict) < 2: continue
      if check(statusdict) == "single MUT": continue
      res = find_intersection(statusdict)
      if res: print brd, res

