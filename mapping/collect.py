#!/usr/bin/env python
# -*- coding:utf-8 -*-

# Standard library packages.
import os
import re
import sys

from collections import defaultdict
from itertools import izip

# Others.
import seeq

from gzopen import gzopen

LOGFNAME = 'tripelog.txt'

class FormatException(Exception):
   pass

def collect_integrations(mapfnames, stcfnames):
   """This function reads the stacode outputs and changes all the
   barcodes mapped by their canonicals while it calculates the
   mapped distance rejecting multiple mapping integrations or
   unmmaped ones. It also counts the frequency that each barcode
   is found in the mapped data even for the non-mapping barcodes."""
   
   def radius(intlist):
      '''Convenience function to compute the radius from a
      list of insertion positions.'''
      intlist.sort()
      try:
         if intlist[0][0] != intlist[-1][0]: return float('inf')
         return intlist[-1][1] - intlist[0][1]
      except IndexError:
         return float('inf')

   def update(bcdict, f):
      '''Convenience function to update a dictionary
      from a starcode file.'''
      for line in f:
         items = line.split()
         for bcd in items[2].split(','):
            bcdict[bcd] = items[0]
   
   canonical = dict()
   # Open all provided starcode files and create a
   # replacement dictionary.
   for fname in stcfnames:
      with gzopen(fname) as f:
         update(canonical, f)


   # The keys of the 'counts' dictionary are the barcodes
   # and the values are couting dictionaries, where the keys
   # are positions (encoded as chromosome, location, strand)
   # and the values are counts.
   counts = defaultdict(lambda: defaultdict(int))

   for fname in mapfiles:
      with open(fname) as f:
         for line in f:
            items = line.split()
            try:
               barcode = canonical[items[0]]
            except KeyError:
               # Barcode has no canonical (can happen).
               continue
            if items[3] == '-':
               if items[2] == '!': position = ('!' , 0)
               else: continue
            else:
               pos = items[3].split(':')
               loc = int(pos[2]) if pos[1] == '+' else \
                     int(pos[2]) + len(items[1])
               position = (pos[0], loc, pos[1])
            counts[barcode][position] += 1
      
   # The keys of the 'integrations' dictionary are the barcodes
   # and the values are pairs of position and total reads.
   integrations = dict()

   for brcd,hist in counts.items():
       total = sum(hist.values())
       top = [pos for pos,count in hist.items() \
             if count > max(1, 0.1*total)]
       # Skip barcode in case of disagreement between top votes.
       if radius(top) > 30: continue
       ins = max(hist, key=hist.get)
       integrations[brcd] = (ins, total)

   # Print barcodes in sorted order in the genome.
   for brcd in sorted(integrations, key=lambda x: integrations.get(x)):
      try:
         (chrom,pos,strand),total = integrations[brcd]
      except ValueError:
         continue
      print '%s\t%s\t%s\t%d\t%d' % (brcd,chrom,strand,pos,total)
   

if __name__ == '__main__':
   # Expects an array of .map files followed by an array
   # of .stc files. Use the extensions to sort them out.
   mapfiles = [f for f in sys.argv[1:] if '.map' in f]
   stcfiles = [f for f in sys.argv[1:] if '.stc' in f]

   collect_integrations(mapfiles, stcfiles)
