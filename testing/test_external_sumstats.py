import unittest
import sys
import os
import io
import math
from contextlib import redirect_stdout

sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')

from Scripts.grouping import ld_threshold, _passes_threshold_eq, _sort_most_significant_last
from Scripts.grouping_model import LDMode, Var, PeakLocus, Grouping
from Scripts.group_annotation import ExtraColAnnotation, FGAnnotation
from Scripts.load_tabix import TabixOptions
from Scripts.data_access.db import Variant

RESOURCE_DIR = os.path.join(os.path.dirname(__file__), "annotate_resources")


class TestMlog10pThresholds(unittest.TestCase):
    def test_ld_threshold_equivalent_raw_vs_mlog10p(self):
        # a dynamic r2 threshold must be identical whether the pval is given raw or as -log10(p)
        raw_pval = 1e-8
        mlog10p = -math.log10(raw_pval)
        raw = ld_threshold(5.0, LDMode.DYNAMIC, raw_pval, pval_is_mlog10p=False)
        mlog = ld_threshold(5.0, LDMode.DYNAMIC, mlog10p, pval_is_mlog10p=True)
        self.assertAlmostEqual(raw, mlog, places=10)

    def test_ld_threshold_static_unaffected(self):
        # constant mode ignores the pval entirely
        self.assertEqual(ld_threshold(0.4, LDMode.CONSTANT, 3.0, pval_is_mlog10p=True), 0.4)

    def test_ld_threshold_no_underflow_past_1e324(self):
        # -log10(p) beyond ~324 underflows 10**(-pval) to 0; the log-space impl must
        # keep resolving thresholds instead of flatlining at a floor value.
        t324 = ld_threshold(5.0, LDMode.DYNAMIC, 324.0, pval_is_mlog10p=True)
        t400 = ld_threshold(5.0, LDMode.DYNAMIC, 400.0, pval_is_mlog10p=True)
        t1000 = ld_threshold(5.0, LDMode.DYNAMIC, 1000.0, pval_is_mlog10p=True)
        # strictly decreasing: a more significant hit gets a stricter (smaller) r2 threshold
        self.assertGreater(t324, t400)
        self.assertGreater(t400, t1000)

    def test_threshold_conversion_roundtrip(self):
        for thresh in (5e-8, 1e-6, 5e-2):
            self.assertAlmostEqual(10 ** (-(-math.log10(thresh))), thresh, places=15)

    def test_passes_threshold_eq_raw(self):
        # raw pval: smaller is more significant, inclusive at the boundary
        self.assertTrue(_passes_threshold_eq(1e-9, 5e-8, False))
        self.assertTrue(_passes_threshold_eq(5e-8, 5e-8, False))
        self.assertFalse(_passes_threshold_eq(1e-7, 5e-8, False))

    def test_passes_threshold_eq_mlog10p(self):
        # mlog10p: larger is more significant, inclusive at the boundary
        thr = -math.log10(5e-8)
        self.assertTrue(_passes_threshold_eq(9.0, thr, True))
        self.assertTrue(_passes_threshold_eq(thr, thr, True))
        self.assertFalse(_passes_threshold_eq(6.0, thr, True))

    def test_sort_most_significant_last(self):
        items = [3.0, 1.0, 2.0]
        # raw pval -> smallest last
        self.assertEqual(_sort_most_significant_last(items, key=lambda x: x, pval_is_mlog10p=False)[-1], 1.0)
        # mlog10p -> largest last
        self.assertEqual(_sort_most_significant_last(items, key=lambda x: x, pval_is_mlog10p=True)[-1], 3.0)


class TestTolerantExtraColumns(unittest.TestCase):
    def _opts(self):
        return TabixOptions(
            os.path.join(RESOURCE_DIR, "finngen_anno.tsv.gz"), "chr", "pos", "ref", "alt")

    def test_missing_column_is_skipped_with_warning(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            anno = ExtraColAnnotation(self._opts(), ["AF", "column_that_does_not_exist"])
        self.assertEqual(anno.extra_columns, ["AF"])
        self.assertIn("skipping", buf.getvalue())
        self.assertIn("column_that_does_not_exist", buf.getvalue())

    def test_all_columns_present_kept(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            anno = ExtraColAnnotation(self._opts(), ["AF", "gene"])
        self.assertEqual(anno.extra_columns, ["AF", "gene"])
        self.assertNotIn("skipping", buf.getvalue())


class TestFinngenVariantsOnly(unittest.TestCase):
    def _peak(self, chrom, pos):
        v = Var(Variant(chrom, pos, "A", "C"), 1e-9, 0.1, None)
        return PeakLocus(v, None, Grouping.NONE)

    def test_lead_variant_membership_filter(self):
        # replicate the --finngen-variants-only lead-variant drop used in main.py
        fg_variants = {Variant("1", 100, "A", "C")}
        loci = [self._peak("1", 100), self._peak("2", 200)]
        kept = [l for l in loci if l.get_vars().lead.id in fg_variants]
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0].get_vars().lead.id.chrom, "1")


if __name__ == "__main__":
    unittest.main()
