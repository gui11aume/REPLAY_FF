#!/usr/bin/env python
# -*- coding:utf-8 -*-
import unittest

import preprocess

class TestExtractor(unittest.TestCase):

   def test_init(self):
      ex = preprocess.Extractor()
      self.assertIsNone(ex.seq_after_tag)
      self.assertIsNone(ex.seq_before_variant)


   def test_Read1Extractor(self):
      ex = preprocess.Read1Extractor()

      # Test case 1.
      seq = 'aaaaaaaaCGCTAATTAATGGAATCATGCGCTACGAGGCCGGCCGCgc'

      BCD,SNP1 = ex.extract_tag_and_variant(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')

      # Test case 2.
      seq = 'aaaaaaaaCGCTAATTAATGGAATCATGttCGCTACGAGGCCGGCCGCgc'

      BCD,SNP1 = ex.extract_tag_and_variant(seq)
      self.assertEqual(BCD, 'aaaaaaaa')
      self.assertEqual(SNP1, 'g')

      # Test case 3.
      seq = 'aaaaaaaaCGCTtATTAAgGGAAaCATGttCGCTACGAcGCgcGCCGCgc'

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
      seq = 'aaaaaaaaCGCTtATTtAgGGAAaCATGttCGCTACGAcGCgcGCCGCgc'
      #                      ^

      with self.assertRaises(preprocess.AberrantReadException):
         ex.extract_tag_and_variant(seq)

      # Test case 6.
      # Extra mutation compared to test case 3.
      seq = 'aaaaaaaaCGCTtATTAAgGGAAaCATGttCcCTACGAcGCgcGCCGCgc'
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
      seq = 'aaaaaaaaCGCTAATTAATGGAATCATGCGCTACGAGGCCGGCCGCgc'

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
