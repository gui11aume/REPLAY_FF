import sys

from collections import defaultdict


class SampleIDException(Exception):
   pass


class ScarcodeReader:
   '''The barcodes have a scar (the scarcode) that identifies the
   kind of mismatch that is generated during the DNA repair. This
   information can be extracted to specify the variants that are
   expected in a given set of reads.'''

   MM = {
      'CGCTAATTAATG': 'CT',
      'GCTAGCAGTCAG': 'CA',
      'GCTAGCTCGTTG': 'GA',
      'GCTAGCTCCGCA': 'GT',
   }

   @staticmethod
   def read(barcode_umi_pairs):

      # Reference scarcodes.
      refscars = {
         'CGCTAATTAATG': 0,
         'GCTAGCAGTCAG': 0,
         'GCTAGCTCGTTG': 0,
         'GCTAGCTCCGCA': 0,
      }

      for pair in barcode_umi_pairs:
         barcode,umi = pair.split('aaaaaaaa')
         scarcode = barcode[-12:]
         if scarcode in refscars:
            refscars[scarcode] += 1
         # Count maximum 1000 barcodes.
         if sum(refscars.values()) > 1000:
            break

      winner = max(refscars, key=refscars.get)
      # The winner must represent more than 90% of the barcodes.
      if refscars[winner] < 0.9 * float(sum(refscars.values())):
         raise SampleIDException

      return ScarcodeReader.MM[winner]


class Regularizer:
   '''In several cases the barcode/UMI pair is associated with multiple
   variants. Read errors and template switching during the PCR can
   cause this. In any event, a barcode/UMI pair corresponds to a unique
   molecule and therefore a single repair event that we can try to
   infer.'''

   def __init__(self):
      self.recors = {}

   def regularize(self, SNPcounts):
      '''Assign the barcode/UMI pair to the most frequent variant.'''
      if len(SNPcounts) > 1:
         variant = max(SNPcounts, key=SNPcounts.get)
      return variant


def display(counter):
#   S = sorted(counter.items())
#   return '\t'.join(['%s/%s:%d' % (a,b,c) for (a,b),c in S])
   S = [('A', 'C'), ('A', 'T'), ('G', 'C')]
   return '\t'.join(['%d' % counter[a] for a in S])


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
      if counter: print display(counter)

#   UMIS = defaultdict(set)
#   BCDS = defaultdict(set)
#   for bcd,umicounts in counts.items():
#      for umi,co in umicounts.items():
#         if sum(co.values()) > 1:
#            BCDS[umi].add(bcd)
#            UMIS[bcd].add(umi)

#   tot = 0
#   for umi in BCDS:
#      if len(BCDS[umi]) == 1: tot += 1
#   print '%d BCDs' % len(UMIS)
#   print '%d UMIs of which %d are unique' % \
#         (len(BCDS), tot)



if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2])
