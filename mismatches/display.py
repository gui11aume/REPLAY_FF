#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
from gzopen import gzopen

TOREMOVE = set()

def remove_barcodes(f):
   for line in f:
      items = line.split()
      TOREMOVE.add(items[0])

   for line in f:
      brcd, ff, at, gc = line.split()
      if brcd in TOREMOVE: continue
      FF, AT, GC = float(ff), float(at), float(gc)
      if FF > 0 and FF > AT and FF > GC:
         FFtot += 1
         continue
      if AT + GC < 2: continue
      if AT > GC:
         ATwins += 1
      if GC > AT:
         GCwins += 1
      ATtot += AT
      GCtot += GC
      tot += 1

def zo(f):
   _ = next(f)
   for line in f:
      brcd, ff, at, gc = line.split()
      if brcd in TOREMOVE: continue
      FF, AT, GC = int(ff), int(at), int(gc)
      if FF > AT and FF > GC: continue
      if AT + GC < 5: continue
      print AT, GC
   return

def count(f):
   ATwins = 0
   GCwins = 0
   _ = next(f)
   for line in f:
      brcd, ff, at, gc = line.split()
      if brcd in TOREMOVE: continue
      FF, AT, GC = int(ff), int(at), int(gc)
      if FF > AT and FF > GC: continue
      if AT > GC and AT < 10:
         ATwins += 1
      if GC > AT and GC < 10:
         GCwins += 1
   return float(ATwins) / float(ATwins + GCwins)

def just_show(f):
   _ = next(f)
   for line in f:
      brcd, ff, at, gc = line.split()
      if brcd in TOREMOVE: continue
      sys.stdout.write(line)


def display(f):
   zo = 0
   ATwins = 0
   GCwins = 0
   FFtot = 0
   ATtot = 0
   GCtot = 0
   tot = 0
   _ = next(f)
   for line in f:
      brcd, ff, at, gc = line.split()
      if brcd in TOREMOVE: continue
      FF, AT, GC = float(ff), float(at), float(gc)
      if FF > 0 and FF > AT and FF > GC:
         FFtot += 1
#         continue
      if AT + GC < 2: continue
      if AT == 0. or GC == 0.:
         zo += 1
      if AT > GC:
         ATwins += 1
      if GC > AT:
         GCwins += 1
      ATtot += AT
      GCtot += GC
      tot += 1
   return float(zo) / float(tot), ATwins / float(ATwins + GCwins), ATtot / float(ATtot + GCtot), FFtot, tot

if __name__ == '__main__':
   for fname in sys.argv[2:]:
      with gzopen(fname) as f:
         remove_barcodes(f)
   with gzopen(sys.argv[1]) as f:
#      just_show(f)
      sys.stdout.write(sys.argv[1])
      sys.stdout.write('\t%.3f\t%.3f\t%.3f\t%d\t%d\n' % display(f))
#      sys.stdout.write('\t%.3f\n' % count(f))
