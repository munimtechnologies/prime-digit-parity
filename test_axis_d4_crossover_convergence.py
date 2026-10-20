# SPDX-License-Identifier: Apache-2.0
"""Pure-standard-library verification of the committed axis-d4-crossover-convergence result.

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
     z-scores and checks it reproduces the committed signature_present for all six sieves;
  5. re-applies the pre-registered P-D4-sf/ex17/ex19/ex23/odd/primes verdicts to the stored
     numbers and asserts the published outcome exactly as graded by the frozen detector:

       * P-D4-sf     (control, full squarefree A005117)       : predicted present -> HELD (present)
       * P-D4-ex17   (PRIMARY, squarefree-except-17)           : predicted present -> HELD (present)
       * P-D4-ex19   (replication, squarefree-except-19)       : predicted present -> HELD (present)
       * P-D4-ex23   (replication, squarefree-except-23)       : predicted present -> HELD (present)
       * P-D4-odd    (anchor, odd squarefree A056911)          : predicted absent  -> HELD (absent)
       * P-D4-primes (negative control, A000040)               : predicted absent  -> HELD (absent)

     This is the CONVERGENCE result that follows the Axis D3 crossover. Axis D (v0.24.0) and
     Axis D2 (v0.26.0) found that relaxing the strict squarefree (mu^2) exponent cap at a single
     prime 2, 3, 5 or 7 destroys the base-2 gap-parity lag-2 unique-quiet signature, with the
     frozen ratio r(p) = |z(2)|/median|z(k!=2)| falling monotonically (1.142, 0.773, 0.342,
     0.169) but staying above the 0.1 detector threshold. Axis D3 (v0.27.0) found a CROSSOVER:
     at the larger primes 11 and 13 the ratio crossed BELOW 0.1 (r=0.0777, 0.0603), so the
     signature RETURNED, falsifying the "destroyed at every single prime" reading. Axis D4
     pre-registered the next question: past the crossover, does r(p) keep falling monotonically
     toward the control (convergence), or is the signature destroyed again at a larger prime (a
     second crossover)? Here the cap is relaxed at 17 (allow 289|n not 4913|n), 19 (allow 361|n
     not 6859|n) and 23 (allow 529|n not 12167|n), all keeping every even integer. The ratio
     continues straight down -- r(17)=0.0420, r(19)=0.0366, r(23)=0.0291 -- all below 0.1, so
     the detector scores all three PRESENT. The pre-registered convergence prediction HOLDS:
     `convergence_to_control_confirmed` is True and `second_crossover_or_redestruction` is
     False. The honest headline is that past the crossover at 11 the single-prime relaxation
     monotonically converges toward the squarefree control r=0.0149 with no second crossover
     through prime 23.

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
DATA = os.path.join(HERE, "data", "axis_d4_crossover_convergence.json")

STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]
SIEVES = ["squarefree", "squarefree_except_17", "squarefree_except_19",
          "squarefree_except_23", "odd_squarefree", "primes"]
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

    def test_convergence_decisive(self):
        # The point of the study: past the Axis D3 crossover at prime 11, does the frozen ratio
        # r(p) = |z(2)|/median|z(k!=2)| keep falling monotonically toward the squarefree control
        # (convergence-to-control), or does a second crossover / re-destruction occur? At the
        # LARGER primes 17, 19, 23 the ratio continues straight down and stays BELOW the 0.1
        # detector threshold, so the lag-2 unique-quiet signature stays PRESENT for all three.
        # Control stays present; anchor (odd squarefree) and the negative control (primes) stay
        # absent. The pre-registered convergence prediction HOLDS; no second crossover appears.
        g = self.d["grading"]["signature_all_sieves"]
        self.assertTrue(g["squarefree"]["signature_present"])
        self.assertTrue(g["squarefree_except_17"]["signature_present"])
        self.assertTrue(g["squarefree_except_19"]["signature_present"])
        self.assertTrue(g["squarefree_except_23"]["signature_present"])
        self.assertFalse(g["odd_squarefree"]["signature_present"])
        self.assertFalse(g["primes"]["signature_present"])
        gr = self.d["grading"]
        self.assertTrue(gr["convergence_to_control_confirmed"])
        self.assertFalse(gr["second_crossover_or_redestruction"])
        self.assertTrue(gr["monotone_trend_continues"])
        self.assertTrue(gr["all_three_below_threshold"])
        # All three relaxed-prime ratios are below the 0.1 threshold (that is what "present" means).
        self.assertLessEqual(g["squarefree_except_17"]["ratio_lag2_over_median_others"], 0.1)
        self.assertLessEqual(g["squarefree_except_19"]["ratio_lag2_over_median_others"], 0.1)
        self.assertLessEqual(g["squarefree_except_23"]["ratio_lag2_over_median_others"], 0.1)

    def test_monotone_trend(self):
        # r(p) falls monotonically across every relaxed prime 2 < 3 < 5 < 7 < 11 < 13 < 17 < 19 <
        # 23, continuing the Axis D2/D3 trend down through the crossover and on toward the control.
        rbp = self.d["grading"]["ratio_by_relaxed_prime"]
        ordered = [rbp[k] for k in ["2", "3", "5", "7", "11", "13", "17", "19", "23"]]
        for a, b in zip(ordered, ordered[1:]):
            self.assertGreater(a, b)
        # The committed Axis D2/D3 anchor values are preserved, and 17/19/23 land below r(13)<0.1.
        self.assertAlmostEqual(rbp["5"], 0.342, places=3)
        self.assertAlmostEqual(rbp["7"], 0.169, places=3)
        self.assertLess(rbp["11"], 0.1)
        self.assertLess(rbp["13"], rbp["11"])
        self.assertLess(rbp["17"], rbp["13"])
        self.assertLess(rbp["19"], rbp["17"])
        self.assertLess(rbp["23"], rbp["19"])

    def test_regrade_matches_published_verdicts(self):
        gr = self.d["grading"]
        ps = self.d["per_sieve"]

        def present(sv):
            p, _, _, _ = lag2_signature(ps[sv]["x_base2"]["lags"])
            return p

        # Control: full squarefree predicted PRESENT; it IS present, so the prediction HELD.
        self.assertEqual(present("squarefree"),
                         gr["P_D4_sf_control"]["squarefree_signature_present"])
        self.assertTrue(gr["P_D4_sf_control"]["holds"])
        self.assertTrue(present("squarefree"))

        # PRIMARY: squarefree-except-17 predicted PRESENT (convergence reading); it IS present,
        # so the primary prediction HELD. This is the headline: no second crossover at 17.
        self.assertEqual(present("squarefree_except_17"),
                         gr["P_D4_ex17"]["squarefree_except_17_signature_present"])
        self.assertTrue(gr["P_D4_ex17"]["holds"])
        self.assertTrue(present("squarefree_except_17"))

        # Replication: squarefree-except-19 predicted PRESENT; it IS present -> HELD.
        self.assertEqual(present("squarefree_except_19"),
                         gr["P_D4_ex19"]["squarefree_except_19_signature_present"])
        self.assertTrue(gr["P_D4_ex19"]["holds"])
        self.assertTrue(present("squarefree_except_19"))

        # Replication: squarefree-except-23 predicted PRESENT; it IS present -> HELD.
        self.assertEqual(present("squarefree_except_23"),
                         gr["P_D4_ex23"]["squarefree_except_23_signature_present"])
        self.assertTrue(gr["P_D4_ex23"]["holds"])
        self.assertTrue(present("squarefree_except_23"))

        # Anchor: odd squarefree predicted ABSENT; it IS absent, so the anchor HELD.
        self.assertEqual(present("odd_squarefree"),
                         gr["P_D4_odd_anchor"]["odd_squarefree_signature_present"])
        self.assertTrue(gr["P_D4_odd_anchor"]["holds"])
        self.assertFalse(present("odd_squarefree"))

        # Negative control: primes predicted ABSENT; it IS absent, so the control HELD.
        self.assertEqual(present("primes"),
                         gr["P_D4_primes_negative_control"]["primes_signature_present"])
        self.assertTrue(gr["P_D4_primes_negative_control"]["holds"])
        self.assertFalse(present("primes"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
