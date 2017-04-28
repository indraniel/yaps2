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

    def test_start_lt(self):
        self.assertFalse(Region._start_lt('22', '10'))
        self.assertFalse(Region._start_lt('22', None))
        self.assertFalse(Region._start_lt(None, None))
        self.assertTrue(Region._start_lt(None, '1'))
        self.assertTrue(Region._start_lt('1', '10'))

    def test_end_lt(self):
        self.assertFalse(Region._end_lt('22', '10'))
        self.assertTrue(Region._end_lt('22', None))
        self.assertFalse(Region._end_lt(None, None))
        self.assertFalse(Region._end_lt(None, '1'))
        self.assertTrue(Region._end_lt('1', '10'))

    def test_hash(self):
        r1 = Region(self.fai_path, 'chr22:1-10000')
        r2 = Region(self.fai_path, 'chr22:1-10000')
        self.assertEqual(r1.__hash__(), r2.__hash__())

    def test_eq(self):
        r1 = Region(self.fai_path, 'chr22:1-10000')
        r2 = Region(self.fai_path, 'chr22:1-10000')
        self.assertTrue(r1 == r2)

        r3 = Region(self.fai_path, 'chr22')
        self.assertFalse(r2 == r3)

    def test_lt(self):
        r1 = Region(self.fai_path, 'chr22:1-10000')
        r2 = Region(self.fai_path, 'chrHLA:1:2:3:4:1-10000')
        self.assertTrue(r1 < r2)
        
        r3 = Region(self.fai_path, 'chr22')
        self.assertTrue(r3 < r1)

        r4 = Region(self.fai_path, 'chr22:10000')
        self.assertTrue(r3 < r4)
        self.assertTrue(r1 < r4)

    def test_sort(self):
        r1 = Region(self.fai_path, 'chr22:1-10000')
        r2 = Region(self.fai_path, 'chrHLA:1:2:3:4:1-10000')
        r3 = Region(self.fai_path, 'chr22:10001-20000')
        r4 = Region(self.fai_path, 'chr22:10001-20001')
        self.assertEqual([ r1, r3, r4, r2 ], sorted([r4, r3, r2, r1]))

