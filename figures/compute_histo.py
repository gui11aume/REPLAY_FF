#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen


def get_values(f):
   '''Parse counting files and extract the values of the
   AT / Total ratios.'''
   # Discard header.
   discard = next(f)
   values = list()
   for line in f:
      (FF,AT,GC) = [int(a) for a in line.split()[1:]]
      if FF > AT and FF > GC: continue
      if AT + GC < 20: continue
      values.append(float(AT) / float(AT + GC))
   return values

def compute_histo(values):
   # Divide the (0,1) interval in 20 categories.
   # For this just multiply by 20 and round down
   counter = defaultdict(int)
   for val in values:
      counter[int(20*val)] += 1
   return counter


if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      values = get_values(f)
   histo = compute_histo(values)
   for k in range(21):
      print histo[k]
