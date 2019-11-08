#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import re
import sys

import torch
import torch.nn as nn
import torch.nn.functional as F

def readin(f):
   entries = list()
   _ = next(f)
   for line in f:
      # Remove first 5 fields.
      entries.append([float(x) for x in line.split()[5:]])
   return torch.tensor(entries)

if __name__ == '__main__':
   ifname = sys.argv[1]
   with open(ifname) as f:
      x = readin(f)
   ofname = re.sub(r'\.[^.]+$', '.tch', ifname)
   torch.save(x, ofname)
