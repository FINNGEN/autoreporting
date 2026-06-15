import unittest
import sys
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts.load_tabix import tb_resource_manager
from Scripts.data_access.db import Variant
import Scripts.grouping as grouping

RESOURCE = "testing/db_resources/test_vcf.vcf.gz"


class TestLoadVariantsByChrom(unittest.TestCase):
    def _all_variants(self, r):
        out = []
        for cols in r.fetch_all_tuples():
            out.append(Variant(cols[r.hdi["#CHROM"]], int(cols[r.hdi["POS"]]),
                               cols[r.hdi["REF"]], cols[r.hdi["ALT"]]))
        return out

    def _per_variant(self, r, variants, cols):
        ref = {}
        for v in variants:
            d = r.load_region(v.chrom, v.pos, v.pos + 1, cols)
            if v in d:
                ref[v] = d[v]
        return ref

    def test_batched_matches_per_variant(self):
        cols = ["POS", "REF"]
        with tb_resource_manager(RESOURCE, "#CHROM", "POS", "REF", "ALT") as r:
            variants = self._all_variants(r)
            ref = self._per_variant(r, variants, cols)
            # small gap -> many clusters; large gap -> a single cluster per chrom
            for max_gap in (5, 1000, 10 ** 9):
                batched = grouping._load_variants_by_chrom(r, variants, cols, max_gap=max_gap)
                self.assertEqual(ref, batched, f"mismatch at max_gap={max_gap}")

    def test_subset_only_returns_requested(self):
        cols = ["POS"]
        with tb_resource_manager(RESOURCE, "#CHROM", "POS", "REF", "ALT") as r:
            variants = self._all_variants(r)
            subset = variants[::3]
            batched = grouping._load_variants_by_chrom(r, subset, cols, max_gap=5)
            self.assertEqual(set(batched.keys()), set(subset))


if __name__ == "__main__":
    unittest.main()
