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

   # Define counters.
   counter = defaultdict(int)
   umicheck = defaultdict(set)

   for line in g:
      try:
         brcd,umi,A,B = line.split()
      except ValueError:
         continue
      try:
         counter[umi] += 1
         umicheck[umi].update([canon[brcd]])
      except KeyError:
         continue

   readslost_recomb = 0
   umilost_recomb = 0
   readslow = 0
   umilow = 0
   for umi in umicheck:
      if len(umicheck[umi]) > 1:
         readslost_recomb += counter[umi]
         umilost_recomb += 1
      if counter[umi] < 2:
         readslow += counter[umi]
         umilow += 1

   print 'total reads: %d' % sum(counter.values())
   print 'lost to recombination: %d' % readslost_recomb
   print 'low: %d' % readslow
   print 'total UMI: %d' % len(umicheck)
   print 'lost to recombination: %d' % umilost_recomb
   print 'low: %d' % umilow

if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f, gzopen(sys.argv[2]) as g:
      main(f,g)
