"""Pure-stdlib checks of the committed cross-sequence Markov-order result.

Verifies data/markov_cross.json against the four pre-registered predictions in
notes/markov-order-cross-prereg.md (P-MC-1..P-MC-4). No numpy, no network; data is
located relative to this file. The global G2 of every (sequence, base) is independently
RE-DERIVED here in pure Python from the committed triple counts and checked against the
job's reported value, so the published result does not have to be trusted on faith.
Run with: python -m unittest -v
Produced with AI by MV2, Munim, Inc.
"""
import json
import math
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SEQS = ["primes", "squarefree", "sums_of_two_squares"]
BASES = ["2", "4", "8", "10", "12", "16"]


def load():
    with open(os.path.join(HERE, "data", "markov_cross.json"), encoding="utf-8") as f:
        return json.load(f)


def g2_triples(counts):
    """Likelihood-ratio G2 for H0: s_{n-1} _||_ s_{n+1} | s_n.

    counts is indexed 4*a + 2*b + c with a=s_{n-1}, b=s_n (the conditioning middle
    bit), c=s_{n+1} -- exactly the job's np_triple_counts convention. For each value of
    the middle bit b we form the 2x2 (a, c) table and sum its independence G2.
    """
    g2 = 0.0
    for b in (0, 1):
        cell = {(a, c): counts[4 * a + 2 * b + c] for a in (0, 1) for c in (0, 1)}
        tot = sum(cell.values())
        if tot == 0:
            continue
        rows = {a: cell[(a, 0)] + cell[(a, 1)] for a in (0, 1)}
        cols = {c: cell[(0, c)] + cell[(1, c)] for c in (0, 1)}
        for a in (0, 1):
            for c in (0, 1):
                o = cell[(a, c)]
                if o == 0:
                    continue
                e = rows[a] * cols[c] / tot
                g2 += 2.0 * o * math.log(o / e)
    return g2


class TestMarkovCrossResults(unittest.TestCase):
    def setUp(self):
        self.d = load()

    def test_self_tests_passed(self):
        st = self.d["self_tests"]
        self.assertTrue(st["all_passed"])
        # synthetic gates: an independent triple table ~0, a dependent one large.
        self.assertLess(st["g_independent_G2"], 1e-6)
        self.assertGreater(st["g_dependent_G2"], 100.0)

    def test_g2_rederivation_from_triple_counts(self):
        # Re-derive every (sequence, base) global G2 and match the reported value.
        for s in SEQS:
            for b in BASES:
                pb = self.d["per_sequence"][s]["per_base"][b]
                mine = g2_triples(pb["triple_counts"])
                self.assertTrue(
                    math.isclose(mine, pb["global"]["G2"], rel_tol=1e-5, abs_tol=1e-4),
                    msg=f"{s} base {b}: re-derived {mine} vs reported {pb['global']['G2']}",
                )

    def test_triple_counts_sum(self):
        for s in SEQS:
            n_gaps = self.d["counts"][s]["n_gaps"]
            for b in BASES:
                pb = self.d["per_sequence"][s]["per_base"][b]
                self.assertEqual(sum(pb["triple_counts"]), pb["n_triples"])
                self.assertEqual(pb["n_triples"], n_gaps - 2)

    def test_pmc1_base2_higher_order_universal(self):
        # P-MC-1: every sieve has base-2 global G2 >= 100.
        for s in SEQS:
            g2 = self.d["per_sequence"][s]["per_base"]["2"]["global"]["G2"]
            self.assertGreaterEqual(g2, 100.0)
        self.assertTrue(self.d["grading"]["P_MC_1_base2_rejects_firstorder_global"]["holds"])

    def test_pmc2_survives_blocking(self):
        # P-MC-2: every sieve has base-2 blocked sigma_over_df >= 5.
        for s in SEQS:
            sig = self.d["per_sequence"][s]["per_base"]["2"]["blocked"]["sigma_over_df"]
            self.assertGreaterEqual(sig, 5.0)
        self.assertTrue(self.d["grading"]["P_MC_2_base2_survives_blocking"]["holds"])

    def test_pmc3_all_six_bases_reject(self):
        # P-MC-3: every base >= 40 in every sieve.
        for s in SEQS:
            for b in BASES:
                g2 = self.d["per_sequence"][s]["per_base"][b]["global"]["G2"]
                self.assertGreaterEqual(g2, 40.0)
        self.assertTrue(self.d["grading"]["P_MC_3_all_six_bases_reject"]["holds"])

    def test_pmc4_loudest_base_sequence_specific(self):
        # P-MC-4: the argmax base is not identical across all three sieves.
        argmax = {}
        for s in SEQS:
            gbb = self.d["per_sequence"][s]["global_G2_by_base"]
            argmax[s] = max(gbb, key=lambda k: gbb[k])
            # agrees with the recorded argmax
            self.assertEqual(argmax[s], self.d["per_sequence"][s]["base_with_largest_global_G2"])
        self.assertGreaterEqual(len(set(argmax.values())), 2)
        self.assertEqual(argmax["primes"], "2")
        self.assertEqual(argmax["squarefree"], "4")
        self.assertEqual(argmax["sums_of_two_squares"], "2")
        self.assertTrue(self.d["grading"]["P_MC_4_base_profile_sequence_specific"]["holds"])

    def test_squarefree_small_gap_degeneracy(self):
        # Disclosed caveat: squarefree max gap is 10, so for base >= 11 a single-digit gap
        # has digit sum = the gap, making bases 12 and 16 identical tables (and base 10
        # nearly so). This documents that the six squarefree bases are ~four distinct tests.
        self.assertEqual(self.d["counts"]["squarefree"]["max_gap"], 10)
        pb = self.d["per_sequence"]["squarefree"]["per_base"]
        self.assertEqual(pb["12"]["triple_counts"], pb["16"]["triple_counts"])

    def test_sanity_bounds(self):
        for s in SEQS:
            for b in BASES:
                g = self.d["per_sequence"][s]["per_base"][b]["global"]
                self.assertGreaterEqual(g["G2"], 0.0)
                self.assertTrue(0.0 <= g["p"] <= 1.0)


if __name__ == "__main__":
    unittest.main()
