#!/usr/bin/env python
# -*- coding:utf-8 -*-

import unittest

from StringIO import StringIO
from textwrap import dedent

import count


class TestTagNormalizer(unittest.TestCase):

   def test_init(self):

      # Mini starcode file.
      f = StringIO(
         'AGATGCTACGCG\t3\tACATGCTACGGC,ACATGCTACGCC\n' \
         'CCATGCTACGAA\t9\tCCATGCTACGAC,CCATGCTACGCA'
      )

      # Make sure that the internal dictionary has been
      # updated upon construction.
      normalizer = count.TagNormalizer(f)
      self.assertNotEqual(normalizer.canonical, dict())

      # Make sure that an exception is raised when the input
      # file is not properly formatted.
      f = StringIO('AGATGCTACGCG')

      with self.assertRaises(ValueError):
         normalizer = count.TagNormalizer(f)


   def test_normalize(self):

      # Mini starcode file.
      f = StringIO(
         'AGATGCTACGCG\t3\tACATGCTACGGC,ACATGCTACGCC\n' \
         'CCATGCTACGAA\t9\tCCATGCTACGAC,CCATGCTACGCA'
      )

      normalizer = count.TagNormalizer(f)

      # Make sure that the normalizer can normalize tags.
      bcd,umi = normalizer.normalize('ACATGCTACGGC')
      self.assertEqual((bcd,umi), ('AG', 'CG'))

      bcd,umi = normalizer.normalize('ACATGCTACGCC')
      self.assertEqual((bcd,umi), ('AG', 'CG'))

      bcd,umi = normalizer.normalize('CCATGCTACGAC')
      self.assertEqual((bcd,umi), ('CC', 'AA'))

      bcd,umi = normalizer.normalize('CCATGCTACGCA')
      self.assertEqual((bcd,umi), ('CC', 'AA'))

      # Make sure that the normalizer raises the proper exception
      # when the tags are not recognized.
      with self.assertRaises(count.AberrantTagException):
         normalizer.normalize('GGATGCTACGGG')

   def test_iter(self):

      # Mini starcode file.
      f = StringIO(
         'AGATGCTACGCG\t3\tACATGCTACGGC,ACATGCTACGCC\n' \
         'CCATGCTACGAA\t9\tCCATGCTACGAC,CCATGCTACGCA'
      )

      normalizer = count.TagNormalizer(f)

      # Make sure that the normalizer can be called directly
      # in a for statement, and use list comprehension for the test.
      tags = sorted([tag for tag in normalizer])
      self.assertEqual(tags, ['AGATGCTACGCG', 'CCATGCTACGAA'])



class TestCountingInfo(unittest.TestCase):

   def test_init(self):
      info = count.CountingInfo('fname1', 'fname2')
      self.assertEqual(info.fname1, 'fname1')
      self.assertEqual(info.fname2, 'fname2')


   def test_get_MMcode(self):

      # Need a CountingInfo instance.
      info = count.CountingInfo('fname1', 'fname2')

      # Test case 1 (GA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCGTTGATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCGTTGATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCGTTGATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'GA')

      # Test case 2 (CT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGCGCTAATTAATGATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGCGCTAATTAATGATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCCGCTAATTAATGAATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'CT')

      # Test case 3 (CA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCAGTCAGATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCAGTCAGATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCAGTCAGATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'CA')

      # Test case 4 (GT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCCGCAATGCTACGTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCCGCAATGCTACGTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCCGCAATGCTACGATCCGTCGGGATACTAAC',
      ])
      mm = info.get_MMcode(tags)
      self.assertEqual(mm, 'GT')


   def test_write_to_file(self):
      info = count.CountingInfo('fname1', 'fname2')
      info.MMcode = 'GA'
      info.nreads = 100
      info.aberrant_tags = 2
      info.thrown_reads = 3
      info.vart_conflicts = [(1,2,3)]
      info.prop_rightMM = 0.9

      buffer = StringIO()
      info.write_to_file(buffer)

      txt = '''fname1
         fname2
         MM type: GA
         Right scarcodes: 90.00%
         Aberrant tags:\t2
         Thrown reads:\t3 (3.00%)
         Recombined reads:\t1
         ---'''

      check = '\n'.join(buffer.getvalue().splitlines()[1:])
      self.assertEqual(check, dedent(' '*9 + txt))


class TestEventCounter(unittest.TestCase):

   def test_init(self):

      # Mini starcode file (GA).
      f = StringIO(
         'GCTAGCTCGTTGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCGTTGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCGTTGATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      info = count.CountingInfo('fname1', 'fname2')
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'GA') 
      
      # Mini starcode file (CT).
      f = StringIO(
         'CGCTAATTAATGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'CGCTAATTAATGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'CGCTAATTAATGATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'CT') 

      # Mini starcode file (CA).
      f = StringIO(
         'GCTAGCAGTCAGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCAGTCAGATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCAGTCAGATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'CA') 

      # Mini starcode file (GT).
      f = StringIO(
         'GCTAGCTCCGCAATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCCGCAATGCTACGT\t1\tGCTAGCAGTCAGATGCTACGT\n' \
         'GCTAGCTCCGCAATGCTACGA\t1\tGCTAGCAGTCAGATGCTACGA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'GT') 


   def test_clip_barcode(self):

      # Mini starcode file (GA).
      f = StringIO(
         'GATGCTAGCTCGTTGATGCTACGTAC\t1\tGATGCTAGCTCGTTGATGCTACGTAC\n' \
         'GATGCTAGCTCGTTGATGCTACGAAA\t1\tGATGCTAGCTCGTTGATGCTACGAAA'
      )

      normalizer = count.TagNormalizer(f)
      info = count.CountingInfo('fname1', 'fname2')
      counter = count.EventCounter(normalizer, info)

      bcd = counter.clip_barcode('GATGCTAGCTCGTTG')
      self.assertEqual(bcd, 'GAT')

      bcd = counter.clip_barcode('GATGCTAGCTCgTTG')
      self.assertEqual(bcd, 'GAT')

      bcd = counter.clip_barcode('GATGCTAGCTCTTG')
      self.assertEqual(bcd, 'GAT')

      with self.assertRaises(count.AberrantTagException):
         counter.clip_barcode('AAAAAAAAAAAAA')


   def test_count(self):
      
      # Mini starcode file (GA).
      f = StringIO(
         'GATGCTAGCTCGTTGATGCTACGTAC\t1\tGATGCTAGCTCGTTGATGCTACGTAC\n' \
         'GATGCTAGCTCGTTGATGCTACGGGG\t1\tGATGCTAGCTCGTTGATGCTACGGGG\n' \
         'GATGCTAGCTCGTTGATGCTACGAAA\t1\tGATGCTAGCTCGTTGATGCTACGAAA'
      )

      normalizer = count.TagNormalizer(f)
      info = count.CountingInfo('fname1', 'fname2')
      counter = count.EventCounter(normalizer, info)

      # Mini pps file.
      f = StringIO(
         'ATGCTACGATGCTACGATGCTACGaa\tA\tC\n' \
         'GATGCTAGCTCGTTGATGCTACGTAC\tA\tC\n' \
         'GATGCTAGCTCGTTGATGCTACGTAC\tA\tC\n' \
         'GATGCTAGCTCGTTGATGCTACGTAC\tA\tT\n' \
         'GATGCTAGCTCGTTGATGCTACGAAA\tA\tC\n' \
         'GATGCTAGCTCGTTGATGCTACGAAA\tA\tC\n' \
         'GATGCTAGCTCGTTGATGCTACGGGG\tA\tC\n'
      )

      out = StringIO()
      counter.count(f, out)

      self.assertEqual(out.getvalue(), 'GAT\t2\t0\t0\n')

      # Test info gathering.
      self.assertEqual(info.MMcode, 'GA')
      self.assertEqual(info.aberrant_tags, 1)
      self.assertEqual(info.thrown_reads, 2)
      self.assertEqual(len(info.vart_conflicts), 1)


if __name__ == '__main__':
   unittest.main()
