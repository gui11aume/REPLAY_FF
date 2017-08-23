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


class TestExtractor(unittest.TestCase):

   def test_init(self):
      ex = preprocess.Extractor()
      self.assertIsNone(ex.seq_after_tag)
      self.assertIsNone(ex.seq_before_variant)


   def test_Read1Extractor(self):
      ex = preprocess.Read1Extractor()

      # Test case 1.
      seq = 'aaaaaaaaGAATCATGAACACCCGCATCGCTACGAGGCCGGCCGCgc'

      BCD,SNP1,SEQ1 = ex.extract_all(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')
      self.assertEqual(SEQ1, 'GAATCATGAACACCCGCATCGCTACGAGGCCGGCCGC')

      # Test case 2.
      seq = 'aaaaaaaaGAATCATGAACACCCGCATttCGCTACGAGGCCGGCCGCgc'

      BCD,SNP1,SEQ1 = ex.extract_all(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')
      self.assertEqual(SEQ1, 'GAATCATGAACACCCGCATttCGCTACGAGGCCGGCCGC')

      # Test case 3.
      seq = 'aaaaaaaaGAAaCATtAACAgCCGCATttCGCTACGAcGCgcGCCGCgc'

      BCD,SNP1,SEQ1 = ex.extract_all(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')
      self.assertEqual(SEQ1, 'GAAaCATtAACAgCCGCATttCGCTACGAcGCgcGCCGC')

      # Test case 4.
      # Structure of a read 2.
      seq = 'aaaaaaaaTGCAACGAATTCATTAGCACCTTGAAGTCGCCGATCAgc'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 5.
      seq = 'GAATCATGAACACCCGCATCGCTACGAGGCCGGCCGCgc'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 6.
      seq = 'aaaaaaaaGAATCATGAACACCCGCATCGCTACGAGGCCGGCCGC'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 7.
      # Extra mutations compared to test case 3.
      seq = 'aaaaaaaaGAAaCATtgACAgCCaCATttCGCTACGAcGCgcGCCGCgc'
      #                      ^      ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 8.
      # Extra mutations compared to test case 3.
      seq = 'aaaaaaaaGAAaCATtAAtAgCCGCATttCGCTACGAcGtgcGCtGCgc'
      #                        ^                    ^    ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)


   def test_Read2Extractor(self):
      ex = preprocess.Read2Extractor()

      # Test case 1.
      seq = 'aaaaaaaaTGCAACGAATTCATTAGCACCTTGAAGTCGCCGATCAgc'

      UMI,SNP2,SEQ2 = ex.extract_all(seq)
      self.assertEqual(UMI, 'aaaaaaaa')
      self.assertEqual(SNP2, 'g')
      self.assertEqual(SEQ2, 'TGCAACGAATTCATTAGCACCTTGAAGTCGCCGATCA')

      # Test case 2.
      seq = 'aaaaaaaaTGCAACGAATTCATTAGttCACCTTGAAGTCGCCGATCAgc'

      UMI,SNP2,SEQ2 = ex.extract_all(seq)
      self.assertEqual(UMI, 'aaaaaaaa')
      self.assertEqual(SNP2, 'g')
      self.assertEqual(SEQ2, 'TGCAACGAATTCATTAGttCACCTTGAAGTCGCCGATCA')

      # Test case 3.
      seq = 'aaaaaaaaTGaAACGtATaCATTAGttCACCTaGAAcaCGCCGATCAgc'

      UMI,SNP2,SEQ2 = ex.extract_all(seq)
      self.assertEqual(UMI, 'aaaaaaaa')
      self.assertEqual(SNP2, 'g')
      self.assertEqual(SEQ2, 'TGaAACGtATaCATTAGttCACCTaGAAcaCGCCGATCA')

      # Test case 4.
      # Structure of a read 1.
      seq = 'aaaaaaaaGAATCATGAACACCCGCATCGCTACGAGGCCGGCCGCgc'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 5.
      seq = 'TGCAACGAATTCATTAGCACCTTGAAGTCGCCGATCAgc'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)
      
      # Test case 6.
      seq = 'aaaaaaaaTGCAACGAATTCATTAGCACCTTGAAGTCGCCGATCA'

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 7.
      # Extra mutations compared to test case 3.
      seq = 'aaaaaaaaTGaAACttAcaCATTAGttCACCTaGAAcaCGCCGATCAgc'
      #                    ^  ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)

      # Test case 8.
      # Extra mutations compared to test case 3.
      seq = 'aaaaaaaaTGaAACttATaCATTAGttCACCTaGAtcaCGCCcATCAgc'
      #                                         ^      ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_all(seq)


if __name__ == '__main__':
   unittest.main(verbosity=2)
