import unittest
from yaps2.utils import empty_gzipped_vcf
import tempfile
import os

class TestRegion(unittest.TestCase):

    def setUp(self):
        self.test_data_dir = os.path.dirname(os.path.abspath(__file__))

    def test_empty_file(self):
        self.assertTrue(empty_gzipped_vcf(os.path.join(self.test_data_dir, 'empty_file.vcf.gz')))

    def test_empty_gzipped_file(self):
        self.assertTrue(empty_gzipped_vcf(os.path.join(self.test_data_dir, 'empty_gzip.vcf.gz')))

    def test_empty_vcf_file(self):
        self.assertTrue(empty_gzipped_vcf(os.path.join(self.test_data_dir, 'empty.vcf.gz')))

    def test_non_empty_vcf_file(self):
        self.assertFalse(empty_gzipped_vcf(os.path.join(self.test_data_dir, 'non_empty.vcf.gz')))

