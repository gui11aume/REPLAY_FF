#!/usr/bin/env python
# -*- coding:utf-8 -*-

# Standard library packages.
import os
import re
import subprocess
import sys
import tempfile

from collections import defaultdict
from itertools import izip

# Others.
import seeq

from gzopen import gzopen

LOGFNAME = 'tripelog.txt'

class FormatException(Exception):
   pass

def trimSuffix(matcher, txt):
   return matcher.matchPrefix(txt, False) or ''

def trimPrefix(matcher, txt):
   return matcher.matchSuffix(txt, False) or ''

########  Mapping Pipeline ###############################################

def extract_reads_from_PE_fastq(fname_iPCR_PE1, fname_iPCR_PE2):
   """This function takes the 2 pair-end sequencing files and extracts the
   barcode making sure that the other read contains the transposon."""

   # Those are the scarcodes that allow to identify which
   # experiment is sequenced (CA, GA or GT mismatch).
   matchers = {
      'CA': seeq.compile('GCTAGCAGTCAGGAATCATG', 3),
      'GA': seeq.compile('GCTAGCTCGTTGGAATCATG', 3),
      'GT': seeq.compile('GCTAGCTCCGCAGAATCATG', 3),
   }

   indexes = {
      'CA': frozenset(['AAC', 'ACA', 'AGG', 'TTC']),
      'GA': frozenset(['ATT', 'CCG', 'TAA', 'TGC']),
      'GT': frozenset(['ACT', 'ATC', 'TGA', 'TGT']),
   }

   outfiles = {
      'AAC': open('CA_AAC.fasta', 'w'),
      'ACA': open('CA_ACA.fasta', 'w'),
      'AGG': open('CA_AGG.fasta', 'w'),
      'TTC': open('CA_TTC.fasta', 'w'),

      'ACT': open('GT_ACT.fasta', 'w'),
      'ATC': open('GT_ATC.fasta', 'w'),
      'TGA': open('GT_TGA.fasta', 'w'),
      'TGT': open('GT_TGT.fasta', 'w'),

      'ATT': open('GA_ATT.fasta', 'w'),
      'CCG': open('GA_CCG.fasta', 'w'),
      'TAA': open('GA_TAA.fasta', 'w'),
      'TGC': open('GA_TGC.fasta', 'w'),
   }

   # End of the pT2 transposon sequence.
   pT2 = seeq.compile('AAACTTCCGACTTCAACTGTA', 3)

   with gzopen(fname_iPCR_PE1) as f, gzopen(fname_iPCR_PE2) as g:
      # Aggregate iterator of f,g iterators -> izip(f,g).
      for lineno,(line1,line2) in enumerate(izip(f,g)):
         # Take sequence lines of the fastq file.
         if lineno % 4 != 1: continue

         # Use the scarcode to identify the experiment.
         for MM,matcher in matchers.items():
            brcd = trimSuffix(matcher, line1)
            # If we find a barcode between 13 and 25 nucleotides
            # then the scarcode must have been the right one.
            if len(brcd) < 13 or len(brcd) > 25: continue

            # Find pT2 on the reverse read. Abort if we cannot or
            # if the remaining sequence is unmappable (too short).
            genome = trimPrefix(pT2, line2.rstrip()).split('CATG')[0]
            if len(genome) < 18: continue

            # The first 3 nucleotides of the reverse read are the
            # index. Check that it belongs to the right group.
            idx = line2[:3]
            if idx in indexes[MM]:
               outf = outfiles[idx]
               outf.write('>%s\n%s\n' % (brcd,genome))
            break



def call_starcode_on_fastq_file(fname_fastq):
   ''' Extracts the gDNA,cDNA reads and spikes and runs stracode on them.'''
   MIN_BRCD = 15
   MAX_BRCD = 25

   brcd_outfname = re.sub(r'\.fastq.*', '_starcode.txt', fname_fastq)
   spk_outfname = re.sub(r'\.fastq.*', '_spikes_starcode.txt', fname_fastq)
   if brcd_outfname == fname_fastq:
      brcd_outfname = fname_fastq + '_starcode.txt'
   if spk_outfname == fname_fastq:
      spk_outfname = fname_fastq + '_spikes_starcode.txt'

   if os.path.exists(brcd_outfname) and os.path.exists(spk_outfname):
      return (brcd_outfname, spk_outfname)

   GFP = seeq.compile('CATGCTAGTTGTGGTTTGTCCAAACT', 4)
   SPIKE = seeq.compile('CATGATTACCCTGTTATC', 2)
   barcode_tempf = tempfile.NamedTemporaryFile(delete=False)
   spike_tempf = tempfile.NamedTemporaryFile(delete=False)
   with gzopen(fname_fastq) as f:
      outf = None
      for lineno,line in enumerate(f):
         if lineno % 4 != 1: continue
         hit = GFP.match(line)
         if hit is not None:
            outf = barcode_tempf
         else:
            hit = SPIKE.match(line)
            if hit is not None:
               outf = spike_tempf
            else:
               continue
         pos = hit.matchlist[0][0]
         if MIN_BRCD <= pos <= MAX_BRCD:
            outf.write(line[:pos] + '\n')
   barcode_tempf.close()
   spike_tempf.close()

   # Skip if file exists.
   if not os.path.exists(brcd_outfname):
      # Call `starcode`.
      subprocess.call([
         'starcode',
         '-t4',
         '-i', barcode_tempf.name,
         '-o', brcd_outfname,
      ])

   if not os.path.exists(spk_outfname):
      subprocess.call([
         'starcode',
         '-t4',
         '-i', spike_tempf.name,
         '-o', spk_outfname,
      ])

   # Delete temporary files.
   os.unlink(barcode_tempf.name)
   os.unlink(spike_tempf.name)

   return (brcd_outfname, spk_outfname)


def collect_integrations(fname_starcode_out, fname_mapped, *args):
   """This function reads the stacode output and changes all the barcodes
   mapped by their canonicals while it calculates the mapped distance
   rejecting multiple mapping integrations or unmmaped ones. It also
   counts the frequency that each barcode is found in the mapped data
   even for the non-mapping barcodes."""
   
   fname_insertions_table = re.sub(r'\.map', '_insertions.txt',
          fname_mapped)
   # Substitution failed, append '_insertions.txt' to avoid name conflict.
   if fname_insertions_table == fname_mapped:
       fname_insertions_table = fname_mapped + '_insertions.txt'

   # Skip if file exists.
   if os.path.exists(fname_insertions_table): return

   def dist(intlist):
      intlist.sort()
      try:
         if intlist[0][0] != intlist[-1][0]: return float('inf')
         return intlist[-1][1] - intlist[0][1]
      except IndexError:
         return float('inf')
   
   canonical = dict()
   with open(fname_starcode_out) as f:
      for line in f:
         items = line.split()
         for brcd in items[2].split(','):
            canonical[brcd] = items[0]

   counts = defaultdict(lambda: defaultdict(int))
   with open(fname_mapped) as f:
      for line in f:
         items = line.split()
         try:
            barcode = canonical[items[0]]
         except KeyError:
            # Barcode has no canonical (can happen).
            continue
         if items[3] == '-':
            if items[2] == '!': position = ('!' , 0)
            else: continue
         else:
            pos = items[3].split(':')
            loc = int(pos[2]) if pos[1] == '+' else \
                  int(pos[2]) + len(items[1])
            position = (pos[0], loc, pos[1])
         counts[barcode][position] += 1
      
   integrations = dict()
   for brcd,hist in counts.items():
       total = sum(hist.values())
       top = [pos for pos,count in hist.items() \
             if count > max(1, 0.1*total)]
       # Skip barcode in case of disagreement between top votes.
       if dist(top) > 30: continue
       ins = max(hist, key=hist.get)
       integrations[brcd] = (ins, total)

   # Count reads from other files.
   reads = dict()
   # First item of tuple is barcode file, second is the spike's one
   for (fname,ignore) in args:
      reads[fname] = defaultdict(int)
      with open(fname) as f:
         for line in f:
            items = line.split('\t')
            try:
               reads[fname][items[0]] = int(items[1])
            except (IndexError, ValueError) as ex:
               raise FormatException("Input file with wrong format")
   with open(fname_insertions_table, 'w') as outf:
      unmapped = 0
      mapped = 0
      for brcd in sorted(integrations, key=lambda x: (integrations.get(x),x)):
         try:
            (chrom,pos,strand),total = integrations[brcd]
         except ValueError:
            continue
         mapped += 1
         outf.write('%s\t%s\t%s\t%d\t%d' % (brcd,chrom,strand,pos,total))
         for fname,ignore in args:
            outf.write('\t' + str(reads[fname][brcd]))
         outf.write('\n')

      # Now add the spikes if the experiment was spiked, otherwise continue.
      N = len(args)
      for i in range(N):
         (ignore,fname) = args[i]
         with open(fname) as f:
            for line in f:
               try:
                  items = line.rstrip().split('\t')
                  array = ['0'] * N
                  array[i] = items[1]
                  outf.write('%s\tspike\t*\t0\t0\t' % items[0])
                  outf.write('\t'.join(array) + '\n')
               except IndexError:
                  continue
   

def main(fname_fastq1, fname_fastq2, *args):
   extract_reads_from_PE_fastq(fname_fastq1, fname_fastq2)


if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2], *sys.argv[3:])