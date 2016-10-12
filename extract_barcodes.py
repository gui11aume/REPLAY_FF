#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from gzopen import gzopen

def main(f):
   for line in f:
      try:
         brcd,umi,A,B = line.split()
      except ValueError:
         continue
      # This condition also eliminates 'None'.
      if len(brcd) > 15: print brcd


if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      main(f)
