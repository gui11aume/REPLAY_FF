#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen

def main(f):
   dictofsets = defaultdict(set)
   for line in f:
      try:
         brcd,umi,A,B = line.split()
      except ValueError:
         continue
      dictofsets[umi].update([brcd])

   for umi in dictofsets:
      print umi, dictofsets[umi]


if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      main(f)
