#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from gzopen import gzopen

def makeset(f):
   S = set()
   for line in f:
      S.add(line.rstrip())
   return S

def main(f, barcode_set_1, barcode_set_2):
   totals_1 = [0,0,0]
   totals_2 = [0,0,0]
   # Discard header.
   next(f)
   for line in f: 
      #if line[0].isspace(): continue
      bcd,FF,AT,GC = line.split()
      scores = [int(a) for a in (FF, AT, GC)]
      winner = max((0,1,2), key=lambda x: scores[x])
      if winner == 0: continue
      if bcd in barcode_set_1:
         totals_1[winner] += 1
      elif bcd in barcode_set_2:
         totals_2[winner] += 1
   return totals_1, totals_2



if __name__ == '__main__':
   with gzopen(sys.argv[2]) as f:
      barcode_set_1 = makeset(f)
   with gzopen(sys.argv[3]) as f:
      barcode_set_2 = makeset(f)

   with gzopen(sys.argv[1]) as f:
      totals_1, totals_2 = main(f, barcode_set_1, barcode_set_2)
   fname = sys.argv[1]
   mmcode = { 'GA':'GA',
         'GT':'GT', 'CA':'CA' }.get(fname[:2].upper(), 'CT')
   tcode = 24 if '24' in fname else 48
   lacode = 'LA' if 'LA' in fname.upper() else '6xPCR'
   ctrl = 'ctrl' if '_c_' in fname or '24c' in fname or \
         '48c' in fname or 'con' in fname else 'test'
   print "%.2f\t%s\t%d\t%s\t%s\t1" % (float(totals_1[1]) / sum(totals_1),
         mmcode, tcode, lacode, ctrl)
   print "%.2f\t%s\t%d\t%s\t%s\t2" % (float(totals_2[1]) / sum(totals_2),
         mmcode, tcode, lacode, ctrl)
