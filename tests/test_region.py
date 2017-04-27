import unittest
from yaps2.utils import Region
import tempfile
import os

class TestRegion(unittest.TestCase):

    def setUp(self):
        temp_fai_descriptor, temp_fai_path = tempfile.mkstemp(suffix='.fai')
        temp_file = os.fdopen(temp_fai_descriptor, 'w')
        temp_file.write('\n'.join(['chr22\t22', 'chrHLA:1:2:3:4\t222']))
        temp_file.close()
        self.fai_path = temp_fai_path

    def tearDown(self):
        os.remove(self.fai_path)

    def test_chromosome_only(self):
        regions = [ 'chr22', 'chrHLA:1:2:3:4' ]
        for region in regions:
            r = Region(self.fai_path, region)
            self.assertEqual(r.chrom, region)
            self.assertIsNone(r.start)
            self.assertIsNone(r.end)

    def test_invalid_chromosome(self):
        with self.assertRaises(RuntimeError):
            r = Region(self.fai_path, 'gobbeldygook')

    def test_short_region(self):
        regions = { 'chr22': '100', 'chrHLA:1:2:3:4': '1000' }
        for chrom, start in regions.iteritems():
            r = Region(self.fai_path, ':'.join((chrom, str(start))))
            self.assertEqual(r.chrom, chrom)
            self.assertEqual(r.start, start)
            self.assertIsNone(r.end)

    def test_long_region(self):
        regions = { 'chr22': ('100', '999'), 'chrHLA:1:2:3:4': ('1000', '1100') }
        for chrom, positions in regions.iteritems():
            r = Region(self.fai_path, ':'.join((chrom, '-'.join(map(str, positions)))))
            self.assertEqual(r.chrom, chrom)
            self.assertEqual(r.start, positions[0])
            self.assertEqual(r.end, positions[1])





