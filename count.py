import sys

from collections import defaultdict

def display(counter):
   S = sorted(counter.items())
   return '\t'.join(['%s/%s:%d' % (a,b,c) for (a,b),c in S])


def main(fname1, fname2):

   # Read in starcode file.
   canonical = dict()

   with open(fname2) as f:
      for line in f:
         centroid,discard,others = line.split()
         for bcd in others.split(','):
            canonical[bcd] = centroid

   # Read in preprocessed file.
   counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

   with open(fname1) as f:
      for line in f:
         tag, SNP1, SNP2 = line.split()
         if tag not in canonical:
            continue
         centroid = canonical[tag]
         bcd,umi = centroid.split('aaaaaaaa')
         counts[bcd][umi][(SNP1,SNP2)] += 1

   for bcd,umicounts in counts.items():
      counter = defaultdict(int)
      for umi,SNPcounts in umicounts.items():
         if sum(SNPcounts.values()) < 2: continue
         SNP = max(SNPcounts, key=SNPcounts.get)
         counter[SNP] += 1
      if counter: print bcd, display(counter)



if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2])
