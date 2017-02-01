#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from gzopen import gzopen

def get_black_list(f):
   S = set()
   for line in f:
      S.add(line.rstrip())
   return S


def makedict(f):
   D = dict()
   for line in f:
      items = line.split()
      D[items[0]] = items
   return D


def main(f, fname, barcode_dict_1, barcode_dict_2, BL):

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
      scores = [int(a) for a in (FF, AT, GC)]
      winner = max((0,1,2), key=lambda x: scores[x])
      if winner == 0: continue
      if bcd in barcode_dict_1:
         the_dict = barcode_dict_1
         rep = 1
      elif bcd in barcode_dict_2:
         the_dict = barcode_dict_2
         rep = 2
      else:
         continue
      the_fields = the_dict[bcd]
      chrom  = the_fields[1]
      strand = the_fields[2]
      pos    = the_fields[3]
      GC1    = the_fields[5]
      GC2    = the_fields[6]
      print "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s" % \
         (bcd, mm[winner], mmcode, tcode, lacode,
               ctrl, rep, GC1, GC2, chrom, strand, pos)



if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      BL = get_black_list(f)
   with gzopen(sys.argv[3]) as f:
      barcode_dict_1 = makedict(f)
   with gzopen(sys.argv[4]) as f:
      barcode_dict_2 = makedict(f)

   with gzopen(sys.argv[2]) as f:
      main(f, sys.argv[2], barcode_dict_1, barcode_dict_2, BL)
