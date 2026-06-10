import unittest
import sys
import pysam
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
#from Scripts.load_tabix import tb_resource_manager,TabixResource
from Scripts.group_annotation import generate_chrom_ranges
from Scripts.data_access.db import Variant
from Scripts.autoreporting_utils import Region

class TestAnnotations(unittest.TestCase):

    def test_generate_ranges(self):
        test_data = [Variant(str(i//15+1),i,"A","C") for i in range(1,31)]
        out = generate_chrom_ranges(set(test_data))
        validation = {'1': [Region(chrom='1', start=0, end=14)],
            '2': [Region(chrom='2', start=14, end=29)],
            '3': [Region(chrom='3', start=29, end=30)]}
        self.assertEqual(out,validation)
        out2 = generate_chrom_ranges(set(test_data),maximum_range_length=4)
        # each wrapped subrange starts at pos-1 (not pos): tabix fetch is half-open [start,end),
        # so the subrange-starting variant must be inside it. (The old expected values started
        # subranges at pos and silently dropped every subrange-starting variant from annotation.)
        validation2 = {'2': [
                Region(chrom='2', start=14, end=17),
                Region(chrom='2', start=17, end=20),
                Region(chrom='2', start=20, end=23),
                Region(chrom='2', start=23, end=26),
                Region(chrom='2', start=26, end=29)],
            '1': [
                Region(chrom='1', start=0, end=3),
                Region(chrom='1', start=3, end=6),
                Region(chrom='1', start=6, end=9),
                Region(chrom='1', start=9, end=12),
                Region(chrom='1', start=12, end=14)],
            '3': [
                Region(chrom='3', start=29, end=30)]}
        self.assertEqual(out2,validation2)

    def test_ranges_cover_all_variants(self):
        """Every variant must be fetchable from some subrange regardless of chunk size.
        tabix fetch(chrom,start,end) is half-open [start,end) in 0-based coords, so a variant
        at 1-based pos (0-based pos-1) is returned iff start < pos <= end. Regression guard for
        the off-by-one that dropped variants starting a wrapped subrange."""
        # positions spanning many chunk boundaries, including pairs that straddle a boundary
        positions = [1, 2, 4, 5, 9, 10, 11, 16, 17, 30, 31, 100]
        variants = {Variant("1", p, "A", "C") for p in positions}
        for maxlen in (1, 2, 3, 4, 5, 7, 100):
            regions = generate_chrom_ranges(variants, maximum_range_length=maxlen)["1"]
            for p in positions:
                self.assertTrue(any(r.start < p <= r.end for r in regions),
                                f"pos {p} not covered by any subrange at maxlen={maxlen}: {regions}")




if __name__=="__main__":
    unittest.main()