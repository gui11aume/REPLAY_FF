#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict

from gzopen import gzopen

REF1 = 'GAATCATGAACACCCGCATCGAGAAGTACGAGGACGGCGGCGTGCTGCACGTGAGCTTCAGCTACCGCTACGAGGCCGGCCGC'
UMIs = defaultdict(lambda: defaultdict(int))

def main(fname):
   with gzopen(fname) as f:
      for line in f:
         items = line.split()
         if int(items[3] < 10): continue
         if float(items[4] < .95): continue
         brd = items[0]
         mm  = items[1]
         seq = items[2]
         UMIs[brd][(seq,mm)] += 1

def clean():
   '''Remove all the sequences that have only one UMI.'''
   for brd,counter in UMIs.items():
      for seq in counter.keys():
         if counter[seq] < 5: counter.pop(seq)

def compare(seq, mm):
   if len(seq) != len(REF1):
      print "%s indel" % mm
      return
   for i in range(len(REF1)):
      if REF1[i] != seq[i]:
         print '%s %s>%s [%d]' % (mm, REF1[i:(i+2)], seq[i:(i+2)], i)
         return

def analyze(counter, brd):
#   seqs = set([a for a,b in counter.keys()])
#   if len(seqs) != 2: return
#   a,b = seqs
#   if    a == REF1: pass
#   elif  b == REF1: (a,b) = (b,a)
#   else: return
#   compare(b, counter.values())
#   return
   cntREF1 = counter.pop((REF1,'FF'), 0) + \
         counter.pop((REF1,'AT'), 0) + counter.pop((REF1,'GC'), 0)
   cntrest = sum(counter.values())
   if cntREF1 < 0.05 * cntrest:
      seq, mm = max(counter, key=counter.get)
#      print cntREF1, cntrest, ','.join([a+'_'+b for a,b in counter.keys()])
      compare(seq, mm)
   return
   if sum(counter.values()) == 0:
      # There was only a WT sequence.
      return
   if cntrest > 2*cntREF1:
      print mm, cntREF1, cntrest, ','.join(counter.keys())

if __name__ == '__main__':
   for fname in sys.argv[1:]:
      main(fname)
      clean()
   total = defaultdict(int)
   for brd,counter in UMIs.items():
      if sum(counter.values()) < 10: continue
      seq, mm = max(counter, key=counter.get)
      total[mm] += 1
      # Count FF vs GC vs AT.
      analyze(counter, brd)
   print 'count %d %d %d' % (total['FF'], total['AT'], total['GC'])
