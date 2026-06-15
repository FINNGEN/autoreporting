import unittest
import sys
import os
import shutil
import tempfile
import pysam

sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import data_access.linkage as linkage
from data_access.db import Variant, LDData

# LD rows per chromosome: (#chrom, pos, variant1, variant2, r2). pos == variant1 position,
# i.e. the fixture is indexed by variant1's position so the narrow-fetch path can be tested.
DATA_ROWS = {
    "1": [
        ("1", "1000", "chr1_1000_A_G", "chr1_1200_C_T", "0.8"),
        ("1", "1000", "chr1_1000_A_G", "chr1_1500_G_A", "0.3"),
        ("1", "1000", "chr1_1000_A_G", "chr1_9000_T_C", "0.9"),  # out of a 2kb window
        ("1", "1200", "chr1_1200_C_T", "chr1_1000_A_G", "0.8"),
        ("1", "1200", "chr1_1200_C_T", "chr1_1400_T_A", "0.6"),
    ],
    "X": [
        ("X", "5000", "chrX_5000_A_T", "chrX_5100_G_C", "0.7"),
    ],
}

ALL_CHROMS = [str(i) for i in range(1, 23)] + ["X"]


def build_fixture(tmpdir):
    """Write a bgzipped+tabix-indexed LD file per chromosome, return the path template."""
    for chrom in ALL_CHROMS:
        tsv = os.path.join(tmpdir, f"test_ld_chr{chrom}.tsv")
        rows = DATA_ROWS.get(chrom)
        if not rows:
            # filler row so empty chromosomes still index; never fetched by the tests
            rows = [(chrom, "1", f"chr{chrom}_1_A_G", f"chr{chrom}_1_A_G", "1.0")]
        with open(tsv, "w") as f:
            f.write("#chrom\tpos\tvariant1\tvariant2\tr2\n")
            for r in rows:
                f.write("\t".join(r) + "\n")
        # seq/start/end columns are 0-based; produces tsv+".gz" and the .tbi index
        pysam.tabix_index(tsv, seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
    return os.path.join(tmpdir, "test_ld_chr{CHROM}.tsv.gz")


class TestTabixLD(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="ld_fixture_")
        cls.template = build_fixture(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_get_range_basic(self):
        ld = linkage.TabixLD(self.template)
        res = ld.get_range(Variant("1", 1000, "A", "G"), 2000)
        ld.close()
        expected = [
            LDData(Variant("1", 1000, "A", "G"), Variant("1", 1200, "C", "T"), 0.8),
            LDData(Variant("1", 1000, "A", "G"), Variant("1", 1500, "G", "A"), 0.3),
            LDData(Variant("1", 1000, "A", "G"), Variant("1", 1000, "A", "G"), 1.0),
        ]
        self.assertEqual(res, expected)

    def test_get_range_threshold(self):
        ld = linkage.TabixLD(self.template)
        res = ld.get_range(Variant("1", 1000, "A", "G"), 2000, 0.5)
        ld.close()
        expected = [
            LDData(Variant("1", 1000, "A", "G"), Variant("1", 1200, "C", "T"), 0.8),
            LDData(Variant("1", 1000, "A", "G"), Variant("1", 1000, "A", "G"), 1.0),
        ]
        self.assertEqual(res, expected)

    def test_chrx_normalization(self):
        ld = linkage.TabixLD(self.template)
        res = ld.get_range(Variant("23", 5000, "A", "T"), 2000)
        ld.close()
        expected = [
            LDData(Variant("23", 5000, "A", "T"), Variant("23", 5100, "G", "C"), 0.7),
            LDData(Variant("23", 5000, "A", "T"), Variant("23", 5000, "A", "T"), 1.0),
        ]
        self.assertEqual(res, expected)

    def test_narrow_vs_wide_equivalence(self):
        """The narrow-fetch optimization must return exactly the wide-fetch result."""
        wide = linkage.TabixLD(self.template, assume_variant1_indexed=False)
        narrow = linkage.TabixLD(self.template, assume_variant1_indexed=True)
        leads = [
            Variant("1", 1000, "A", "G"),
            Variant("1", 1200, "C", "T"),
            Variant("23", 5000, "A", "T"),
        ]
        for lead in leads:
            self.assertEqual(wide.get_range(lead, 2000), narrow.get_range(lead, 2000))
        wide.close()
        narrow.close()

    def test_get_ranges_serial_matches_loop(self):
        ld = linkage.TabixLD(self.template)
        leads = [Variant("1", 1000, "A", "G"), Variant("1", 1200, "C", "T")]
        batch = ld.get_ranges(leads, 2000, workers=1)
        loop = {v: ld.get_range(v, 2000) for v in leads}
        ld.close()
        self.assertEqual(batch, loop)

    def test_get_ranges_parallel_matches_serial(self):
        ld = linkage.TabixLD(self.template)
        leads = [Variant("1", 1000, "A", "G"), Variant("1", 1200, "C", "T")]
        serial = ld.get_ranges(leads, 2000, workers=1)
        parallel = ld.get_ranges(leads, 2000, workers=2)
        ld.close()
        self.assertEqual(serial, parallel)

    def test_get_ranges_per_lead_ranges(self):
        ld = linkage.TabixLD(self.template)
        leads = [Variant("1", 1000, "A", "G")]
        # a 100bp range drops the 1200/1500 partners; only the self-pair remains
        res = ld.get_ranges(leads, 2000, bp_ranges=[100], workers=1)
        ld.close()
        self.assertEqual(res[leads[0]], [LDData(leads[0], leads[0], 1.0)])


if __name__ == "__main__":
    unittest.main()
