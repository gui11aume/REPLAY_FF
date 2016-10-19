#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen

def trim(counter, umicheck):
   for umi in counter.keys():
      if len(umicheck[umi]) > 1:
         counter.pop(umi)

def main(f,g):
   # Read canonical barcodes.
   canon = dict()
   for line in f:
      canonical,count,barcodes = line.split()
      for brcd in barcodes.split(','):
         canon[brcd] = canonical

   # Define counters.
   TA = defaultdict(int)
   CG = defaultdict(int)
   TG = defaultdict(int)
   counters = {
      ('T','A'): TA,
      ('C','G'): CG,
      ('T','G'): TG,
   }

   umicheck = defaultdict(set)

   for line in g:
      try:
         brcd,umi,A,B = line.split()
         if brcd == "None": continue
      except ValueError:
         continue
      try:
         counters[(A,B)][umi] += 1
         umicheck[umi].update([canon[brcd]])
      except KeyError:
         continue

   # Remove single read events.
   trim(TA, umicheck)
   trim(CG, umicheck)
   trim(TG, umicheck)

   print 'TA: %f' % (sum(TA.values()) / float(len(TA)))
   print 'CG: %f' % (sum(CG.values()) / float(len(CG)))
   print 'TG: %f' % (sum(TG.values()) / float(len(TG)))

if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f, gzopen(sys.argv[2]) as g:
      main(f,g)
