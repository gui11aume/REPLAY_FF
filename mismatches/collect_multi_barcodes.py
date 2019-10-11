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
   _ = next(f)
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

   nsyn = 0
   npairs = 0
   for bcd, reps in BRCD2REP.items():
      if len(reps) < 2: continue
      for a,b in combinations(reps, 2):
         if a == 0 and b == 0: continue
         if a == 1 and b == 1:
            nsyn += 1
         npairs += 1
   print nsyn, npairs
