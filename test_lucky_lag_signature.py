# SPDX-License-Identifier: Apache-2.0
"""Pure-standard-library verification of the committed lucky-lag-signature result.

Runs offline in well under a second (no numpy, no network). It independently:
  1. re-derives the exact gap identity P_2(g) = P_4(g) XOR P_4(floor(g/2)) that defines the
     three streams, and that P_2 == popcount parity (Thue-Morse);
  2. brute-forces the exact marginal-preserving permutation-null moments on tiny sequences
     and checks them against the closed-form formula used below (so the formula itself is
     validated, not assumed);
  3. for EVERY sieve/stream/lag, re-derives E_null and var_null from only the stored
     (n_eff, s, S) using that validated formula, and checks they match the committed
     numbers and that z = (S - E_null)/sqrt(var_null) -- a full re-derivation, not a
     re-division;
  4. re-runs the scale-free lag-2 unique-quiet signature detector on the stored x_base2
     z-scores and checks it reproduces the committed signature_present for all four sieves;
  5. re-applies the pre-registered P-L-1/2/3 thresholds to the stored numbers and asserts
     the published verdicts (all three hold).

Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.
"""
import json
import math
import os
import statistics
import unittest
from fractions import Fraction
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "lucky_lag_signature.json")

STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]
SIEVES = ["primes", "squarefree", "sums_of_two_squares", "lucky"]
LAGS = [1, 2, 3, 4, 5, 6, 7, 8]


def digit_sum_parity(n, b):
    s = 0
    while n:
        s += n % b
        n //= b
    return s & 1


def null_moments(n_eff, s, k):
    """Exact marginal-preserving permutation-null E[S^(k)] and Var[S^(k)] (floats)."""
    N_ = n_eff
    M = N_ - k
    n_shared = max(0, N_ - 2 * k)
    n_disjoint = M * (M - 1) // 2 - n_shared
    rho2 = Fraction(s * s - N_, N_ * (N_ - 1))
    rho4 = Fraction(
        s ** 4 - 6 * N_ * s * s + 8 * s * s + 3 * N_ * N_ - 6 * N_,
        N_ * (N_ - 1) * (N_ - 2) * (N_ - 3),
    )
    eS = M * rho2
    eS2 = M + 2 * n_shared * rho2 + 2 * n_disjoint * rho4
    var = eS2 - eS * eS
    return float(eS), float(var)


def lag2_signature(x_lags):
    """Scale-free 'lag-2 unique-quiet' detector on the x_base2 z-profile.
    Present iff the smallest |z| among lags 1..8 is at lag 2 AND |z(2)| <= 0.1 x the
    median |z| over the other seven lags (a ratio test, independent of sample size)."""
    az = {k: abs(x_lags[str(k)]["z"]) for k in LAGS}
    quiet_lag = min(az, key=az.get)
    others = [az[k] for k in LAGS if k != 2]
    med_others = statistics.median(others)
    present = bool(quiet_lag == 2 and az[2] <= 0.1 * med_others)
    return present, quiet_lag, az[2], med_others


class TestStrideIdentity(unittest.TestCase):
    def test_xor_identity_and_thue_morse(self):
        for g in range(0, 8192):
            a = digit_sum_parity(g, 4)
            b = digit_sum_parity(g >> 1, 4)
            self.assertEqual(a ^ b, digit_sum_parity(g, 2), g)
            self.assertEqual(digit_sum_parity(g, 2), bin(g).count("1") % 2, g)


class TestNullFormula(unittest.TestCase):
    def test_bruteforce_matches_formula(self):
        # Exhaustive permutation average/variance vs the closed form, k = 1..4.
        for seq in ([1, 1, 1, -1, -1, -1, 1, -1], [1, 1, 1, 1, -1, -1, 1]):
            N_ = len(seq)
            s = sum(seq)
            for k in (1, 2, 3, 4):
                tot = tot2 = cnt = 0
                for perm in permutations(seq):
                    Sk = sum(perm[i] * perm[i + k] for i in range(N_ - k))
                    tot += Sk
                    tot2 += Sk * Sk
                    cnt += 1
                emp_mean = tot / cnt
                emp_var = tot2 / cnt - emp_mean * emp_mean
                eS, var = null_moments(N_, s, k)
                self.assertAlmostEqual(emp_mean, eS, places=9, msg=(seq, k))
                self.assertAlmostEqual(emp_var, var, places=9, msg=(seq, k))


class TestCommittedResult(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(DATA, encoding="utf-8") as f:
            cls.d = json.load(f)

    def test_self_tests_passed(self):
        self.assertTrue(self.d["self_tests"]["all_passed"])

    def test_shape(self):
        self.assertEqual(self.d["n_max"], 20000000)
        for sv in SIEVES:
            for st in STREAMS:
                self.assertIn(st, self.d["per_sieve"][sv])
                for k in LAGS:
                    self.assertIn(str(k), self.d["per_sieve"][sv][st]["lags"])

    def test_null_moments_and_z_rederived(self):
        # Re-derive E_null, var_null and z from only (n_eff, s, S) for every cell.
        for sv in SIEVES:
            for st in STREAMS:
                r = self.d["per_sieve"][sv][st]
                n_eff, s = r["n_eff"], r["s"]
                for k in LAGS:
                    cell = r["lags"][str(k)]
                    eS, var = null_moments(n_eff, s, k)
                    self.assertAlmostEqual(eS, cell["E_null"], places=4, msg=(sv, st, k))
                    self.assertAlmostEqual(var, cell["var_null"], delta=abs(var) * 1e-9 + 1e-3,
                                           msg=(sv, st, k))
                    z = (cell["S"] - eS) / math.sqrt(var)
                    self.assertAlmostEqual(z, cell["z"], places=4, msg=(sv, st, k))

    def test_s_matches_marginal(self):
        for sv in SIEVES:
            for st in STREAMS:
                r = self.d["per_sieve"][sv][st]
                frac = (r["n_eff"] - r["s"]) / 2 / r["n_eff"]
                self.assertAlmostEqual(frac, r["marginal_frac_odd"], places=9)

    def test_signature_detector_matches_grading(self):
        g = self.d["grading"]["signature_all_sieves"]
        for sv in SIEVES:
            present, quiet, azl2, med = lag2_signature(self.d["per_sieve"][sv]["x_base2"]["lags"])
            self.assertEqual(present, g[sv]["signature_present"], sv)
            self.assertEqual(quiet, g[sv]["quiet_lag"], sv)
            self.assertAlmostEqual(azl2, g[sv]["abs_z_lag2"], places=6, msg=sv)

    def test_decisive_split(self):
        # The point of the study: signature present for the two mu^2-type multiplicative
        # sieves, absent for primes AND the positional lucky sieve of matched density.
        g = self.d["grading"]["signature_all_sieves"]
        self.assertTrue(g["squarefree"]["signature_present"])
        self.assertTrue(g["sums_of_two_squares"]["signature_present"])
        self.assertFalse(g["primes"]["signature_present"])
        self.assertFalse(g["lucky"]["signature_present"])

    def test_regrade_matches_published_verdicts(self):
        gr = self.d["grading"]
        lucky = self.d["per_sieve"]["lucky"]

        # P-L-1: lucky base-2 lag-1 anti-persistence (z <= -8).
        z1 = lucky["x_base2"]["lags"]["1"]["z"]
        pl1 = z1 <= -8.0
        self.assertEqual(pl1, gr["P_L_1_lucky_base2_antipersistence"]["holds"])

        # P-L-2 (PRIMARY): lucky does NOT show the lag-2 unique-quiet signature.
        present, _, _, _ = lag2_signature(lucky["x_base2"]["lags"])
        pl2 = not present
        self.assertEqual(pl2, gr["P_L_2_lucky_no_lag2_signature"]["holds"])

        # P-L-3: both lucky base-4 strides loud at lag 1 (|z| >= 5).
        za = abs(lucky["a_base4_even"]["lags"]["1"]["z"])
        zb = abs(lucky["b_base4_odd"]["lags"]["1"]["z"])
        pl3 = za >= 5.0 and zb >= 5.0
        self.assertEqual(pl3, gr["P_L_3_lucky_base4_strides_loud"]["holds"])

        # All three pre-registered predictions held.
        self.assertTrue(pl1)
        self.assertTrue(pl2)
        self.assertTrue(pl3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
