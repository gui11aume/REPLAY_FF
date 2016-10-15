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
         'AGaaaaaaaaCG\t3\tACaaaaaaaaGC,ACaaaaaaaaCC\n' \
         'CCaaaaaaaaAA\t9\tCCaaaaaaaaAC,CCaaaaaaaaCA'
      )

      # Make sure that the internal dictionary has been
      # updated upon construction.
      normalizer = count.TagNormalizer(f)
      self.assertNotEqual(normalizer.canonical, dict())

      # Make sure that an exception is raised when the input
      # file is not properly formatted.
      f = StringIO('AGaaaaaaaaCG')

      with self.assertRaises(ValueError):
         normalizer = count.TagNormalizer(f)


   def test_normalize(self):

      # Mini starcode file.
      f = StringIO(
         'AGaaaaaaaaCG\t3\tACaaaaaaaaGC,ACaaaaaaaaCC\n' \
         'CCaaaaaaaaAA\t9\tCCaaaaaaaaAC,CCaaaaaaaaCA'
      )

      normalizer = count.TagNormalizer(f)

      # Make sure that the normalizer can normalize tags.
      bcd,umi = normalizer.normalize('ACaaaaaaaaGC')
      self.assertEqual((bcd,umi), ('AG', 'CG'))

      bcd,umi = normalizer.normalize('ACaaaaaaaaCC')
      self.assertEqual((bcd,umi), ('AG', 'CG'))

      bcd,umi = normalizer.normalize('CCaaaaaaaaAC')
      self.assertEqual((bcd,umi), ('CC', 'AA'))

      bcd,umi = normalizer.normalize('CCaaaaaaaaCA')
      self.assertEqual((bcd,umi), ('CC', 'AA'))

      # Make sure that the normalizer raises the proper exception
      # when the tags are not recognized.
      with self.assertRaises(count.AberrantTagException):
         normalizer.normalize('GGaaaaaaaaGG')

   def test_iter(self):

      # Mini starcode file.
      f = StringIO(
         'AGaaaaaaaaCG\t3\tACaaaaaaaaGC,ACaaaaaaaaCC\n' \
         'CCaaaaaaaaAA\t9\tCCaaaaaaaaAC,CCaaaaaaaaCA'
      )

      normalizer = count.TagNormalizer(f)

      # Make sure that the normalizer can be called directly
      # in a for statement, and use list comprehension for the test.
      tags = sorted([tag for tag in normalizer])
      self.assertEqual(tags, ['AGaaaaaaaaCG', 'CCaaaaaaaaAA'])


   def test_read_on_TagNormalizer(self):
      
      # Mini starcode file (CA).
      f = StringIO(
         'GCTAGCAGTCAGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCAGTCAGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCAGTCAGaaaaaaaaA\t1\tGCTAGCAGTCAGaaaaaaaaA'
      )

      normalizer = count.TagNormalizer(f)
      mm = count.EventCounter.get_MMcode(normalizer)
      self.assertEqual(mm, 'CA')


class TestCountingInfo(unittest.TestCase):

   def test_init(self):
      info = count.CountingInfo('fname1', 'fname2')
      self.assertEqual(info.fname1, 'fname1')
      self.assertEqual(info.fname2, 'fname2')

   def test_write_to_file(self):
      info = count.CountingInfo('fname1', 'fname2')
      info.MMcode = 'GA'
      info.nreads = 100
      info.aberrant_tags = 2
      info.thrown_reads = 3
      info.vart_conflicts = [(1,2,3)]

      buffer = StringIO()
      info.write_to_file(buffer)

      txt = '''fname1
         fname2
         MM type: GA
         Aberrant tags:\t2
         Thrown reads:\t3 (3.00%)
         Recombined reads:\t1
         ---'''

      check = '\n'.join(buffer.getvalue().splitlines()[1:])
      self.assertEqual(check, dedent(' '*9 + txt))


class TestEventCounter(unittest.TestCase):

   def test_init(self):

      # Mini starcode file.
      f = StringIO(
         'AGaaaaaaaaCG\t3\tACaaaaaaaaGC,ACaaaaaaaaCC\n' \
         'CCaaaaaaaaAA\t9\tCCaaaaaaaaAC,CCaaaaaaaaCA'
      )

      normalizer = count.TagNormalizer(f)
      info = count.CountingInfo('fname1', 'fname2')

      # Make sure that the proper exception is raised when
      # reading the wrong file format.
      with self.assertRaises(count.SampleIDException):
         count.EventCounter(normalizer, info)

      # Mini starcode file (GA).
      f = StringIO(
         'GCTAGCTCGTTGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCTCGTTGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCTCGTTGaaaaaaaaA\t1\tGCTAGCAGTCAGaaaaaaaaA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'GA') 
      
      # Mini starcode file (CT).
      f = StringIO(
         'CGCTAATTAATGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'CGCTAATTAATGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'CGCTAATTAATGaaaaaaaaA\t1\tGCTAGCAGTCAGaaaaaaaaA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'CT') 

      # Mini starcode file (CA).
      f = StringIO(
         'GCTAGCAGTCAGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCAGTCAGaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCAGTCAGaaaaaaaaA\t1\tGCTAGCAGTCAGaaaaaaaaA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'CA') 

      # Mini starcode file (GT).
      f = StringIO(
         'GCTAGCTCCGCAaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCTCCGCAaaaaaaaaT\t1\tGCTAGCAGTCAGaaaaaaaaT\n' \
         'GCTAGCTCCGCAaaaaaaaaA\t1\tGCTAGCAGTCAGaaaaaaaaA'
      )

      normalizer = count.TagNormalizer(f)
      counter = count.EventCounter(normalizer, info)

      self.assertEqual(counter.info.MMcode, 'GT') 


   def test_get_MMcode(self):

      # Test case 1 (GA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCGTTGaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCGTTGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCGTTGaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.EventCounter.get_MMcode(tags)
      self.assertEqual(mm, 'GA')

      # Test case 2 (CT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGCGCTAATTAATGaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGCGCTAATTAATGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCCGCTAATTAATGAaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.EventCounter.get_MMcode(tags)
      self.assertEqual(mm, 'CT')

      # Test case 3 (CA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCAGTCAGaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCAGTCAGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCAGTCAGaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.EventCounter.get_MMcode(tags)
      self.assertEqual(mm, 'CA')

      # Test case 4 (GT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCCGCAaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCCGCAaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCCGCAaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.EventCounter.get_MMcode(tags)
      self.assertEqual(mm, 'GT')

      # Make sure an exception is raised if more than 10% of the
      # scarcodes are from different samples.
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCCGCAaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCAGTCAGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCGTTGaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      with self.assertRaises(count.SampleIDException):
         count.EventCounter.get_MMcode(tags)


   def test_clip_barcode(self):

      # Mini starcode file (GA).
      f = StringIO(
         'GATGCTAGCTCGTTGaaaaaaaaTAC\t1\tGATGCTAGCTCGTTGaaaaaaaaTAC\n' \
         'GATGCTAGCTCGTTGaaaaaaaaAAA\t1\tGATGCTAGCTCGTTGaaaaaaaaAAA'
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
         'GATGCTAGCTCGTTGaaaaaaaaTAC\t1\tGATGCTAGCTCGTTGaaaaaaaaTAC\n' \
         'GATGCTAGCTCGTTGaaaaaaaaGGG\t1\tGATGCTAGCTCGTTGaaaaaaaaGGG\n' \
         'GATGCTAGCTCGTTGaaaaaaaaAAA\t1\tGATGCTAGCTCGTTGaaaaaaaaAAA'
      )

      normalizer = count.TagNormalizer(f)
      info = count.CountingInfo('fname1', 'fname2')
      counter = count.EventCounter(normalizer, info)

      # Mini pps file.
      f = StringIO(
         'aaaaaaaaaaaaaaaaaaaaaaaaaa\tA\tC\n' \
         'GATGCTAGCTCGTTGaaaaaaaaTAC\tA\tC\n' \
         'GATGCTAGCTCGTTGaaaaaaaaTAC\tA\tC\n' \
         'GATGCTAGCTCGTTGaaaaaaaaTAC\tA\tT\n' \
         'GATGCTAGCTCGTTGaaaaaaaaAAA\tA\tC\n' \
         'GATGCTAGCTCGTTGaaaaaaaaAAA\tA\tC\n' \
         'GATGCTAGCTCGTTGaaaaaaaaGGG\tA\tC\n'
      )

      out = StringIO()
      counter.count(f, out)

      self.assertEqual(out.getvalue(), 'GAT\t2\t0\t0\n')

      # Test info gathering.
      self.assertEqual(info.MMcode, 'GA')
      self.assertEqual(info.aberrant_tags, 1)
      self.assertEqual(info.thrown_reads, 1)
      self.assertEqual(len(info.vart_conflicts), 1)


if __name__ == '__main__':
   unittest.main()
