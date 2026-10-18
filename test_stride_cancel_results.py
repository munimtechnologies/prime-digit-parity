# SPDX-License-Identifier: Apache-2.0
"""Pure-standard-library verification of the committed stride-cancellation result.

Runs offline in a few milliseconds. It (1) re-derives the exact gap identity
P_2(g) = P_4(g) XOR P_4(floor(g/2)) used to define the three streams, (2) checks the
stored z-scores are internally consistent with their stored null moments, (3) checks the
grading block's abs_z matches the per-sieve z magnitudes, and (4) re-applies the
pre-registered P1/P2/P3 thresholds to the stored numbers and asserts the published
verdicts. No network, no numpy. Produced with AI by MV2, Munim, Inc.
"""
import json
import math
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "stride_cancel_results.json")

STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]
SIEVES = ["primes", "squarefree", "sums_of_two_squares"]


def digit_sum_parity(n, b):
    s = 0
    while n:
        s += n % b
        n //= b
    return s & 1


class TestStrideIdentity(unittest.TestCase):
    def test_xor_identity_and_thue_morse(self):
        # P_2(g) = P_4(g) XOR P_4(floor(g/2)), and P_2(g) == popcount parity.
        for g in range(0, 8192):
            a = digit_sum_parity(g, 4)
            b = digit_sum_parity(g >> 1, 4)
            self.assertEqual(a ^ b, digit_sum_parity(g, 2), g)
            self.assertEqual(digit_sum_parity(g, 2), bin(g).count("1") % 2, g)


class TestCommittedResult(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(DATA, encoding="utf-8") as f:
            cls.d = json.load(f)

    def test_self_tests_passed(self):
        self.assertTrue(self.d["self_tests"]["all_passed"])

    def test_shape(self):
        self.assertEqual(self.d["n_max"], 100000000)
        for sv in SIEVES:
            for st in STREAMS:
                self.assertIn(st, self.d["per_sieve"][sv])

    def test_z_consistent_with_null_moments(self):
        # z == (S1 - E_null) / sqrt(var_null) for every stored stream.
        for sv in SIEVES:
            for st in STREAMS:
                r = self.d["per_sieve"][sv][st]
                z = (r["S1"] - r["E_null"]) / math.sqrt(r["var_null"])
                self.assertAlmostEqual(z, r["z"], places=6, msg=(sv, st))

    def test_s_matches_marginal(self):
        # s = (#0) - (#1) implies marginal_frac_odd = (n_eff - s) / 2 / n_eff.
        for sv in SIEVES:
            for st in STREAMS:
                r = self.d["per_sieve"][sv][st]
                frac = (r["n_eff"] - r["s"]) / 2 / r["n_eff"]
                self.assertAlmostEqual(frac, r["marginal_frac_odd"], places=9)

    def test_absz_block_matches_per_sieve(self):
        az = self.d["grading"]["abs_z"]
        keymap = {"a": "a_base4_even", "b": "b_base4_odd", "x_base2": "x_base2"}
        for sv in SIEVES:
            for short, full in keymap.items():
                self.assertAlmostEqual(
                    az[sv][short],
                    abs(self.d["per_sieve"][sv][full]["z"]),
                    places=6, msg=(sv, short),
                )

    def test_regrade_matches_published_verdicts(self):
        # Re-apply the pre-registered thresholds (30 / 100 / 10) independently.
        az = self.d["grading"]["abs_z"]
        s2a, s2b, s2x = az["sums_of_two_squares"]["a"], az["sums_of_two_squares"]["b"], az["sums_of_two_squares"]["x_base2"]
        pa, pb, px = az["primes"]["a"], az["primes"]["b"], az["primes"]["x_base2"]
        sa, sb, sx = az["squarefree"]["a"], az["squarefree"]["b"], az["squarefree"]["x_base2"]

        p1_cancel = max(s2a, s2b) >= 30 and s2x <= 30
        p1_falsified = max(s2a, s2b) < 30
        p2 = px >= 100 and min(pa, pb) >= 10
        p3 = sx >= 100 and min(sa, sb) >= 10

        g = self.d["grading"]
        self.assertEqual(p1_cancel, g["P1_sumsq_cancellation"])
        self.assertEqual(p1_falsified, g["P1_sumsq_falsified_no_signal"])
        self.assertEqual(p2, g["P2_primes_both_streams_contribute"])
        self.assertEqual(p3, g["P3_squarefree_both_streams_contribute"])

        # The published verdicts: cancellation confirmed, not falsified, P2 and P3 hold.
        self.assertTrue(p1_cancel)
        self.assertFalse(p1_falsified)
        self.assertTrue(p2)
        self.assertTrue(p3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
