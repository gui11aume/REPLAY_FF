#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from gzopen import gzopen

def get_black_list(f):
   S = set()
   for line in f:
      S.add(line.rstrip())
   return S


def main(f, fname, BL):

   mmcode = { 'GA':'GA',
         'GT':'GT', 'CA':'CA' }.get(fname[:2].upper(), 'CT')
   tcode = 48 if '48' in fname else 24
   lacode = 'LA' if 'LA' in fname.upper() else '6xPCR'
   ctrl = 'ctrl' if '_c_' in fname or '24c' in fname or \
         '48c' in fname or 'con' in fname else 'test'

   mm = ['FF', 'AT', 'GC']

   rep = GC1 = GC2 = None
   # Discard header.
   next(f)
   for line in f: 
      #if line[0].isspace(): continue
      bcd,FF,AT,GC = line.split()
      # Skip barcodes in the black list.
      if bcd in BL: continue
      scores = [float(a) for a in (FF, AT, GC)]
      if scores[0] >= max(scores)-1: continue
      winner = max((0,1,2), key=lambda x: scores[x])
      if scores[winner] < 2: continue
      ratio = scores[1] / (scores[1] + scores[2])
      print "%s\t%.3f\t%s\t%d\t%s\t%s" % \
         (bcd, ratio, mmcode, tcode, lacode, ctrl)



if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      BL = get_black_list(f)

   with gzopen(sys.argv[2]) as f:
      main(f, sys.argv[2], BL)
