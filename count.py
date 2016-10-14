import sys

from collections import defaultdict


class SampleIDException(Exception):
   pass

class AberrantTagException(Exception):
   pass


class ExperimentInfo:
   '''Basic information about the experiment.'''

   def __init__(self, code):

      assign = {
         'GA': (('C', 'A'), ('T', 'A'), ('C', 'G')),
      }

      self.FF, self.AT, self.GC = assign[code]
      self.varconflicts = []

   def get_variant(self, varcounts):
      '''In several cases the barcode/UMI pair is associated with
      multiple variants. Read errors and template switching during the
      PCR can cause this. In any event, a barcode/UMI pair corresponds
      to a unique molecule and therefore a single repair event that we
      can try to infer.'''

      # If variants are not unanimous for the barcode/UMI pair,
      # create an exception entry for the records.
      if len(varcounts) > 1:
         entry = [
            varcounts[self.FF],
            varcounts[self.AT],
            varcounts[self.GC],
         ]
         self.varconflicts.append(entry)

      # The most frequent variant is assigned to the pair.
      variant = max(varcounts, key=varcounts.get)
      return variant


class TagNormalizer:
   '''Correct reading errors.'''

   def __int__(self, fname):
      '''Construct a normalizer form Starcode file.'''

      self.canonical = dict()
      # Build a dictionary from Starcode file.
      with open(fname) as f:
         for line in f:
            centroid,discard,others = line.split()
            for bcd in others.split(','):
               self.canonical[bcd] = centroid

   def __iter__(self):
      '''Make this class iterable so that statements of the form
      'for tag in TagNormalizer' iterate over the dictionary
      of canonical tags. This allows to use a TagNormalizer
      as input of the ScarcodeReader.'''
      return self.canonical

   def normalize(self, tag):
      '''Replace the tag by its canonical sequence (where errors
      are reverted). Sperate the barcode from the UMI and returns
      both as a pair. In case of failure, an AberrantTag exception
      is raised.'''

      try:
         barcode, umi = self.canonical[tag].split('aaaaaaaa')
      except (KeyError, ValueError):
         raise AberrantTagException

      return barcode, umi



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


def display(counter):
#   S = sorted(counter.items())
#   return '\t'.join(['%s/%s:%d' % (a,b,c) for (a,b),c in S])
   S = [('A', 'C'), ('A', 'T'), ('G', 'C')]
   return '\t'.join(['%d' % counter[a] for a in S])


def main(fname1, fname2):

   # Instantiate basic analysis tools.
   normalizer = TagNormalizer(fname2)
   info = ExperimentInfo(ScarcodeReader.read(normalizer))

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
