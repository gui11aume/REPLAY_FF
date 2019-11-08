#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen
from itertools import combinations

BRCD2REP = defaultdict(list)
TOREMOVE = set()

def remove_barcodes(f):
   _ = next(f)
   for line in f:
      bcd, ff, at, gc = line.split()
      if bcd in TOREMOVE: continue
      FF = float(ff)
      AT = float(at)
      GC = float(gc)
      if FF > AT and FF > GC: continue
      if AT + GC < 3: continue
      TOREMOVE.add(bcd)


def collect(f):
   try:
      _ = next(f)
   except StopIteration:
      return
   for line in f:
      bcd, ff, at, gc = line.split()
      if bcd in TOREMOVE: continue
      FF = float(ff)
      AT = float(at)
      GC = float(gc)
      if FF > AT and FF > GC: continue
      if AT + GC < 5: continue
      x = AT / (AT+GC)
      # Round here.
      BRCD2REP[bcd].append(int(round(x)))


if __name__ == '__main__':
   for fname in sys.argv[5:]:
      with gzopen(fname) as f:
         remove_barcodes(f)

   for fname in sys.argv[1:5]:
      with gzopen(fname) as f:
         collect(f)

   tot = defaultdict(int)
   for bcd, reps in BRCD2REP.items():
      if len(reps) < 2: continue
#      print '\t'.join([str(x) for x in reps])
      for a,b in combinations(reps, 2):
         tot[(a,b)] += 1
   print tot[(0,0)], tot[(0,1)]+tot[(1,0)], tot[(1,1)]
