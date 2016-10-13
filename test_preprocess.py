#!/usr/bin/env python
# -*- coding:utf-8 -*-

import unittest

from StringIO import StringIO
from textwrap import dedent

import preprocess


class TestLaneInfo(unittest.TestCase):

   def test_constructor(self):
      info = preprocess.LaneInfo('fname1', 'fname2')
      self.assertEqual(info.fname1, 'fname1')
      self.assertEqual(info.fname2, 'fname2')
      self.assertEqual(info.ntotal, 0)
      self.assertEqual(info.naberrant, 0)

   def test_write_to_file(self):
      info = preprocess.LaneInfo('fname1', 'fname2')
      info.ntotal = 1000000
      info.naberrant = 100000

      buffer = StringIO()
      info.write_to_file(buffer)

      txt = '''fname1
         fname2
         Reads lost:\t100000 (10.00 %)
         ---'''

      check = '\n'.join(buffer.getvalue().splitlines()[1:])
      self.assertEqual(check, dedent(' '*9 + txt))

      buffer.close()


class TestExtractor(unittest.TestCase):

   def test_init(self):
      ex = preprocess.Extractor()
      self.assertIsNone(ex.seq_after_tag)
      self.assertIsNone(ex.seq_before_variant)


   def test_Read1Extractor(self):
      ex = preprocess.Read1Extractor()

      # Test case 1.
      seq = 'aaaaaaaaGAATCATGAACACCCGCATCGCTACGAGGCCGGCCGCgc'

      BCD,SNP1 = ex.extract_tag_and_variant(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')

      # Test case 2.
      seq = 'aaaaaaaaGAATCATGAACACCCGCATttCGCTACGAGGCCGGCCGCgc'

      BCD,SNP1 = ex.extract_tag_and_variant(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')

      # Test case 3.
      seq = 'aaaaaaaaGAAaCATtAACAgCCGCATttCGCTACGAcGCgcGCCGCgc'

      BCD,SNP1 = ex.extract_tag_and_variant(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')

      # Test case 4.
      # Structure of a read 2.
      seq = 'aaaaaaaaTGCAACGAATTCATTAGCACCTTGAAGTCGCCGATCAgc'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_tag_and_variant(seq)

      # Test case 5.
      # Extra mutation compared to test case 3.
      seq = 'aaaaaaaaGAAaCATtgACAgCCGCATttCGCTACGAcGCgcGCCGCgc'
      #                      ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_tag_and_variant(seq)

      # Test case 6.
      # Extra mutation compared to test case 3.
      seq = 'aaaaaaaaGAAaCATtAACAgCCGCATttCGaTACGAcGCgcGCCGCgc'
      #                                     ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_tag_and_variant(seq)


   def test_Read2Extractor(self):
      ex = preprocess.Read2Extractor()

      # Test case 1.
      seq = 'aaaaaaaaTGCAACGAATTCATTAGCACCTTGAAGTCGCCGATCAgc'

      UMI,SNP2 = ex.extract_tag_and_variant(seq)
      self.assertEqual(UMI, 'aaaaaaaa')
      self.assertEqual(SNP2, 'g')

      # Test case 2.
      seq = 'aaaaaaaaTGCAACGAATTCATTAGttCACCTTGAAGTCGCCGATCAgc'

      UMI,SNP2 = ex.extract_tag_and_variant(seq)
      self.assertEqual(UMI, 'aaaaaaaa')
      self.assertEqual(SNP2, 'g')

      # Test case 3.
      seq = 'aaaaaaaaTGaAACGtATaCATTAGttCACCTaGAAcaCGCCGATCAgc'

      UMI,SNP2 = ex.extract_tag_and_variant(seq)
      self.assertEqual(UMI, 'aaaaaaaa')
      self.assertEqual(SNP2, 'g')

      # Test case 4.
      # Structure of a read 1.
      seq = 'aaaaaaaaGAATCATGAACACCCGCATCGCTACGAGGCCGGCCGCgc'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_tag_and_variant(seq)

      # Test case 5.
      # Extra mutation compared to test case 3.
      seq = 'aaaaaaaaTGaAACttATaCATTAGttCACCTaGAAcaCGCCGATCAgc'
      #                    ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_tag_and_variant(seq)

      # Test case 6.
      # Extra mutation compared to test case 3.
      seq = 'aaaaaaaaTGaAACttATaCATTAGttCACCTaGAtcaCGCCGATCAgc'
      #                                         ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_tag_and_variant(seq)


if __name__ == '__main__':
   unittest.main()
