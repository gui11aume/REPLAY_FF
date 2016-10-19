#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen

def main(f,g):
   # Read canonical barcodes.
   canon = dict()
   for line in f:
      canonical,count,barcodes = line.split()
      for brcd in barcodes.split(','):
         canon[brcd] = canonical

   # Convenience.
   formatted = {
      ('T', 'A'): 'TA',
      ('C', 'G'): 'CG',
      ('T', 'G'): 'TG',
   }

   # Counter.
   dictofdict = defaultdict(dict)

   # UMI checks
   umicheck = defaultdict(set)
   sets = {('T','A'): set(), ('C','G'): set(), ('T','G'): set()}

   for line in g:
      try:
         brcd,umi,A,B = line.split()
      except ValueError:
         continue
      try:
         umicheck[umi].update([canon[brcd]])
         sets[(A,B)].update(umi)
         dictofdict[brcd][umi] = (A,B)
      except KeyError:
         continue

   ambiguous = (sets[('T','A')] ^ sets[('C','G')]) | \
               (sets[('T','A')] ^ sets[('T','G')]) | \
               (sets[('C','G')] ^ sets[('T','G')])

   for brcd in dictofdict:
      counter = defaultdict(int)
      for umi,var in dictofdict[brcd].items():
         # Skip UMIs associated with several barcodes.
         if len(umicheck[umi]) > 1: continue
         # Skip ambiguous UMIs.
         if umi in ambiguous: continue
         counter[var] += 1
      # Skip barcodes with too few UMIs.
      if 2 < sum(counter.values()):
         win = sorted(counter, key=counter.get, reverse=True)[0]
         # Check that the winner has > 90% UMIs.
         if counter[win] < 0.9 * sum(counter.values()): continue
         print brcd, formatted[win], counter[('T','A')], \
                  counter[('C','G')], counter[('T','G')]


if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f, gzopen(sys.argv[2]) as g:
      main(f,g)
