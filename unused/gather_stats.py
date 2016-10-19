#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

def main(f):

   umiTA = list(); umiCG = list(); umiTG = list(); 

   for line in f:
      barcode,win,TA,CG,TG = line.split()
      if   win == 'TA': umiTA.append(int(TA))
      elif win == 'CG': umiCG.append(int(CG))
      else            : umiTG.append(int(TG))

   print 'TA %.2f' % (sum(umiTA) / float(len(umiTA)))
   print 'CG %.2f' % (sum(umiCG) / float(len(umiCG)))
   print 'TG %.2f' % (sum(umiTG) / float(len(umiTG)))
      

if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      main(f)
