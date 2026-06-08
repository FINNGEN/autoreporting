import unittest
import sys
import random

sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')

import grouping
from grouping import _greedy_ld_group
from grouping_model import (Var, GroupingOptions, SummstatColumns, Grouping, LDMode,
                            PeakLocus)
from data_access.db import Variant, LDData


def reference_ld_group(p1_piles, p2_piles, ld_cache, options):
    """Faithful copy of the original per-iteration list-rebuild greedy loop, kept here as
    the equivalence oracle for the optimized _greedy_ld_group."""
    output = []
    for seq in p1_piles.keys():
        p1_pile = grouping._sort_most_significant_last(p1_piles[seq], key=lambda x: x.pval,
                                                       pval_is_mlog10p=options.pval_is_mlog10p)
        p2_pile = list(p2_piles[seq])
        while p1_pile:
            lead_var = p1_pile.pop()
            lead_variant = lead_var.id
            ld_thresh = grouping.ld_threshold(options.r2_threshold, options.ld_mode,
                                              lead_var.pval, options.pval_is_mlog10p)
            ld_data = [a for a in ld_cache[lead_variant] if a.r2 > ld_thresh]
            ld_data = [a for a in ld_data if (a.variant1 == lead_variant) and (a.variant2 != lead_variant)]
            ld_data_varset = set([a.variant2 for a in ld_data])
            ld_partners = [a for a in p2_pile if a.id in ld_data_varset]
            ld_index_dict = {a.variant2: i for i, a in enumerate(ld_data)}
            ld_partners = [Var(a.id, a.pval, a.beta, ld_data[ld_index_dict[a.id]].r2) for a in ld_partners]
            p1_pile = [a for a in p1_pile if a.id not in ld_data_varset]
            if not options.overlap:
                p2_pile = [a for a in p2_pile if a.id not in ld_data_varset]
            output.append(PeakLocus(lead_var, ld_partners, Grouping.LD))
    return output


def make_options(overlap, r2_threshold=0.2):
    cols = SummstatColumns("chrom", "pos", "ref", "alt", "pval", "beta")
    # LDMode.CONSTANT makes ld_threshold deterministic (== r2_threshold)
    return GroupingOptions("x", Grouping.LD, cols, 1_000_000, LDMode.CONSTANT,
                           r2_threshold, 5e-8, 1e-5, overlap)


def loci_key(loci):
    # full equality including partner list order, ids, pval/beta and r2_to_lead
    return [(l.lead, tuple(l.ld_partners)) for l in loci]


def random_scenario(rng, n_chrom=2, n_var=14):
    """Build (p1_piles, p2_piles, ld_cache) for a random scenario.

    p2 is the superset of variants; p1 is a random subset of leads (all also in p2). Each
    lead's cache entry holds random LD partners (drawn from all variants) plus the self-pair.
    """
    p1_piles, p2_piles, ld_cache = {}, {}, {}
    for c in range(1, n_chrom + 1):
        seq = str(c)
        variants = [Variant(seq, 1000 + 10 * i, "A", "G") for i in range(n_var)]
        vars_pb = {v: (rng.uniform(1e-12, 1e-4), rng.uniform(-1, 1)) for v in variants}
        p2 = [Var(v, vars_pb[v][0], vars_pb[v][1], 1.0) for v in variants]
        rng.shuffle(p2)  # exercise file-order independence
        leads = [v for v in p2 if rng.random() < 0.5]
        p1_piles[seq] = list(leads)
        p2_piles[seq] = p2
        for lv in leads:
            partners = [v for v in variants if v != lv.id and rng.random() < 0.4]
            entry = [LDData(lv.id, v2, round(rng.uniform(0.0, 1.0), 3)) for v2 in partners]
            entry.append(LDData(lv.id, lv.id, 1.0))
            ld_cache[lv.id] = entry
    return p1_piles, p2_piles, ld_cache


class TestGreedyEquivalence(unittest.TestCase):
    def _check(self, p1, p2, ld, options):
        # pass independent copies so neither implementation can perturb the other
        ref = reference_ld_group({k: list(v) for k, v in p1.items()},
                                 {k: list(v) for k, v in p2.items()}, ld, options)
        new = _greedy_ld_group({k: list(v) for k, v in p1.items()},
                               {k: list(v) for k, v in p2.items()}, ld, options)
        self.assertEqual(loci_key(ref), loci_key(new))

    def test_random_overlap(self):
        rng = random.Random(1234)
        for _ in range(200):
            p1, p2, ld = random_scenario(rng)
            self._check(p1, p2, ld, make_options(overlap=True))

    def test_random_no_overlap(self):
        rng = random.Random(5678)
        for _ in range(200):
            p1, p2, ld = random_scenario(rng)
            self._check(p1, p2, ld, make_options(overlap=False))

    def test_random_thresholds(self):
        rng = random.Random(99)
        for _ in range(200):
            p1, p2, ld = random_scenario(rng)
            thr = rng.choice([0.0, 0.1, 0.5, 0.9])
            ov = rng.random() < 0.5
            self._check(p1, p2, ld, make_options(overlap=ov, r2_threshold=thr))


if __name__ == "__main__":
    unittest.main()
