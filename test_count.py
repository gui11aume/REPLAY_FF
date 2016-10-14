#!/usr/bin/env python
# -*- coding:utf-8 -*-

import unittest

import count


class TestScarcodeReader(unittest.TestCase):

   def test_read(self):

      # Test case 1 (GA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCGTTGaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCGTTGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCGTTGaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.ScarcodeReader.read(tags)
      self.assertEqual(mm, 'GA')

      # Test case 2 (CT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGCGCTAATTAATGaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGCGCTAATTAATGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCCGCTAATTAATGAaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.ScarcodeReader.read(tags)
      self.assertEqual(mm, 'CT')

      # Test case 3 (CA).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCAGTCAGaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCAGTCAGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCAGTCAGaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.ScarcodeReader.read(tags)
      self.assertEqual(mm, 'CA')

      # Test case 4 (GT).
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCCGCAaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCTCCGCAaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCCGCAaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      mm = count.ScarcodeReader.read(tags)
      self.assertEqual(mm, 'GT')

      # Make sure an exception is raised if more than 10% of the
      # scarcodes are from different samples.
      tags = set([
         'TTCGTGAGATAAATCAGTTGGCTAGCTCCGCAaaaaaaaaTGCGTCGGACAGCGACGC',
         'AAACTCATCTAAACGTTTTGGCTAGCAGTCAGaaaaaaaaTATCTGGCTTCCCGGCCA',
         'CACGCTCTGCATGTTTCCCAGCTAGCTCGTTGaaaaaaaaATCCGTCGGGATACTAAC',
      ])
      with self.assertRaises(count.SampleIDException):
         count.ScarcodeReader.read(tags)

if __name__ == '__main__':
   unittest.main()
