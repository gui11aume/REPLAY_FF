#!/usr/bin/env python
# -*- coding:utf-8 -*-

import datetime
import seeq
import sys

from collections import defaultdict


class SampleIDException(Exception):
   pass

class AberrantTagException(Exception):
   pass


class EventCounter:
   '''Basic information about the experiment.'''

   # MMcodes and their scarcodes.
   MM = {
      'CGCTAATTAATG': 'CT',
      'GCTAGCAGTCAG': 'CA',
      'GCTAGCTCGTTG': 'GA',
      'GCTAGCTCCGCA': 'GT',
   }

   def __init__(self, normalizer, info):

      # Events are recorded in a nested counter.
      # barcode: umi: variant: count.
      self.events = defaultdict(lambda:
            defaultdict(lambda: defaultdict(int)))

      self.info = info

      self.normalizer = normalizer
      self.normalize_tag = normalizer.normalize

      assign = {
         #          FF          AT          GC
         'GA': (('A', 'C'), ('A', 'T'), ('G', 'C')),
      }

      clip = {
         'CT': seeq.compile('CGCTAATTAATG', 2),
         'CA': seeq.compile('GCTAGCAGTCAG', 2),
         'GA': seeq.compile('GCTAGCTCGTTG', 2),
         'GT': seeq.compile('GCTAGCTCCGCA', 2),
      }

      self.info.MMcode = self.get_MMcode(normalizer)

      self.MM = assign[self.info.MMcode]
      self.clip = clip[self.info.MMcode]


   def clip_barcode(self, bcd):
      bcd = self.clip.matchPrefix(bcd, False)
      if bcd is None:
         raise AberrantTagException
      return bcd


   @staticmethod
   def get_MMcode(tags):
      '''The barcodes have a scar (the scarcode) that identifies the
      kind of mismatch that is generated during the DNA repair. This
      information can be extracted to specify the variants that are
      expected in a given set of reads.'''

      # Scarcodes and the corresponding mismatches.

      refscars = {
         'CGCTAATTAATG': 0,
         'GCTAGCAGTCAG': 0,
         'GCTAGCTCGTTG': 0,
         'GCTAGCTCCGCA': 0,
      }

      for tag in tags:
         bcd,umi = tag.split('aaaaaaaa')
         scarcode = bcd[-12:]
         if scarcode in refscars:
            refscars[scarcode] += 1
         # Count maximum 1000 barcodes.
         if sum(refscars.values()) > 1000:
            break

      winner = max(refscars, key=refscars.get)
      # The winner must represent more than 90% of the barcodes.
      if not refscars[winner] > 0.9 * float(sum(refscars.values())):
         raise SampleIDException

      return EventCounter.MM[winner]


   def normalize_variant(self, dict_of_variants):
      '''In several cases the barcode/UMI pair is associated with
      multiple variants. Read errors and template switching during the
      PCR can cause this. In any event, a barcode/UMI pair corresponds
      to a unique molecule and therefore a single repair event that we
      can try to infer.'''

      # If variants are not unanimous for the barcode/UMI pair,
      # create an exception entry for the records.
      if len(dict_of_variants) > 1:
         entry = [dict_of_variants[a] for a in self.MM ]
         self.info.vart_conflicts.append(entry)

      # The most frequent variant is assigned to the pair.
      variant = max(dict_of_variants, key=dict_of_variants.get)
      return variant


   def count(self, f, outf=sys.stdout):
      '''Processing function to convert the reads to repair events.'''

      # Create a filter for UMIs used multiple times.
      reverse_lookup = defaultdict(set)

      for line in f:
         self.info.nreads += 1
         # Assume that preprocessed file is tab-separated.
         tag, V1, V2 = line.split()
         try:
            bcd, umi = self.normalize_tag(tag)
            bcd = self.clip_barcode(bcd)
            reverse_lookup[umi].add(bcd)
         except AberrantTagException:
            self.info.aberrant_tags += 1
            continue
         self.events[bcd][umi][(V1,V2)] += 1

      for bcd,dict_of_umis in self.events.items():
         counter = defaultdict(int)
         for umi,dict_of_variants in dict_of_umis.items():
            # Discard all unique reads.
            if sum(dict_of_variants.values()) < 2:
               self.info.thrown_reads += 1
               continue
            # Discard UMIs used multiple times.
            if len(reverse_lookup[umi]) > 1:
               self.info.thrown_reads += sum(dict_of_variants.values())
               continue
            variant = self.normalize_variant(dict_of_variants)
            counter[variant] += 1
         # If all UMIs were single read, the counter
         # is empty and there is nothing to show.
         if not counter:
            continue
         # Show counts.
         counts = '\t'.join(['%d' % counter[a] for a in self.MM])
         outf.write('%s\t%s\n' % (bcd, counts))



class TagNormalizer:
   '''Correct reading errors.'''

   def __init__(self, f):
      '''Construct a normalizer from open Starcode file.'''

      self.canonical = dict()
      # Build a dictionary from Starcode file.
      for line in f:
        centroid,discard,others = line.split()
        for bcd in others.split(','):
            self.canonical[bcd] = centroid

   def __iter__(self):
      '''Make this class iterable so that statements of the form
      'for tag in TagNormalizer' iterate over the canonical
      tags.'''
      return iter(set(self.canonical.values()))

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


class CountingInfo:

   def __init__(self, fname1, fname2):
      self.fname1 = fname1
      self.fname2 = fname2

      self.nreads = 0

      self.vart_conflicts = []
      self.aberrant_tags = 0
      self.thrown_reads = 0


   def write_to_file(self, f):
      '''Log the counting process.'''

      # Time stamp.
      dt = datetime.datetime
      f.write(dt.strftime(dt.now(),
         '%Y-%m-%d %H:%M:%S Counting summary\n'))

      # Processed files.
      f.write('%s\n' % self.fname1)
      f.write('%s\n' % self.fname2)

      # Mismatch type.
      f.write('MM type: %s\n' % self.MMcode)

      # Total and percent reads lost.
      f.write('Aberrant tags:\t%d\n' % self.aberrant_tags)
      percent = 100 * float(self.thrown_reads) / self.nreads
      f.write('Thrown reads:\t%d (%.2f%%)\n' % \
            (self.thrown_reads, percent))
      f.write('Recombined reads:\t%d\n' % len(self.vart_conflicts))

      # End of report.
      f.write('---\n')



def main(fname1, fname2, info):
   # Instantiate and run tools for the analysis.
   with open(fname1) as f, open(fname2) as g:
      normalizer = TagNormalizer(g)
      counter = EventCounter(normalizer, info)
      counter.count(f)


if __name__ == '__main__':
   info = CountingInfo(sys.argv[1], sys.argv[2])
   try:
      main(sys.argv[1], sys.argv[2], info)
   finally:
      info.write_to_file(open('counting_logs.txt', 'a'))
