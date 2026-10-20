# SPDX-License-Identifier: Apache-2.0
"""Pure-standard-library verification of the committed axis-d2-generality-lag-signature result.

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
     z-scores and checks it reproduces the committed signature_present for all five sieves;
  5. re-applies the pre-registered P-D2-sf/ex5/ex7/odd/primes verdicts to the stored numbers
     and asserts the published outcome exactly as graded by the frozen detector:

       * P-D2-sf    (control, full squarefree A005117)       : predicted present -> HELD (present)
       * P-D2-ex5   (PRIMARY, squarefree-except-5)            : predicted absent  -> HELD (absent)
       * P-D2-ex7   (replication, squarefree-except-7)        : predicted absent  -> HELD (absent)
       * P-D2-odd   (anchor, odd squarefree A056911)          : predicted absent  -> HELD (absent)
       * P-D2-primes (negative control, A000040)              : predicted absent  -> HELD (absent)

     This is the GENERALITY replication of Axis D (v0.24.0). Axis D relaxed the squarefree
     exponent cap at the two smallest relevant primes 2 and 3 and found the lag-2 unique-quiet
     signature destroyed; the open question was whether that destruction is a general property
     of relaxing the cap at ANY single prime or an artefact of 2/3. Here the cap is relaxed at
     the LARGER primes 5 (allow 25|n but not 125|n) and 7 (allow 49|n but not 343|n), both
     keeping every even integer. Following the Axis D result the pre-registration predicted
     both would again be ABSENT; both are, so the predictions HELD and
     `strict_mu2_requirement_generalizes_to_larger_primes` is True: the signature needs EXACTLY
     the strict squarefree (mu^2) congruence structure, not merely the small primes 2/3.

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
DATA = os.path.join(HERE, "data", "axis_d2_generality_lag_signature.json")

STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]
SIEVES = ["squarefree", "squarefree_except_5", "squarefree_except_7",
          "odd_squarefree", "primes"]
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
        # The point of the study: the lag-2 unique-quiet signature is present for the full strict
        # squarefree sieve, but ABSENT for every single-prime exponent-cap relaxation at the LARGER
        # primes 5 (allow 25|n) and 7 (allow 49|n) -- exactly as it was for the small primes 2/3 in
        # Axis D -- even though both keep every even integer. So relaxing the mu^2 cap at any single
        # prime destroys it; the destruction is NOT an artefact of 2/3. The anchor (odd squarefree,
        # no evens) and the negative control (primes) are absent too. The strict-mu^2 requirement
        # therefore GENERALIZES to larger primes.
        g = self.d["grading"]["signature_all_sieves"]
        self.assertTrue(g["squarefree"]["signature_present"])
        self.assertFalse(g["squarefree_except_5"]["signature_present"])
        self.assertFalse(g["squarefree_except_7"]["signature_present"])
        self.assertFalse(g["odd_squarefree"]["signature_present"])
        self.assertFalse(g["primes"]["signature_present"])
        self.assertTrue(self.d["grading"]["strict_mu2_requirement_generalizes_to_larger_primes"])

    def test_regrade_matches_published_verdicts(self):
        gr = self.d["grading"]
        ps = self.d["per_sieve"]

        def present(sv):
            p, _, _, _ = lag2_signature(ps[sv]["x_base2"]["lags"])
            return p

        # Control: full squarefree predicted PRESENT; it IS present, so the prediction HELD.
        self.assertEqual(present("squarefree"),
                         gr["P_D2_sf_control"]["squarefree_signature_present"])
        self.assertTrue(gr["P_D2_sf_control"]["holds"])
        self.assertTrue(present("squarefree"))

        # PRIMARY: squarefree-except-5 (relax the cap at the larger prime 5, keep evens) predicted
        # ABSENT following Axis D; it IS absent, so the primary prediction HELD.
        self.assertEqual(present("squarefree_except_5"),
                         gr["P_D2_ex5_primary"]["squarefree_except_5_signature_present"])
        self.assertTrue(gr["P_D2_ex5_primary"]["holds"])
        self.assertFalse(present("squarefree_except_5"))

        # Replication: squarefree-except-7 predicted ABSENT; it IS absent -> HELD.
        self.assertEqual(present("squarefree_except_7"),
                         gr["P_D2_ex7_replication"]["squarefree_except_7_signature_present"])
        self.assertTrue(gr["P_D2_ex7_replication"]["holds"])
        self.assertFalse(present("squarefree_except_7"))

        # Anchor: odd squarefree predicted ABSENT; it IS absent, so the anchor HELD.
        self.assertEqual(present("odd_squarefree"),
                         gr["P_D2_odd_anchor"]["odd_squarefree_signature_present"])
        self.assertTrue(gr["P_D2_odd_anchor"]["holds"])
        self.assertFalse(present("odd_squarefree"))

        # Negative control: primes predicted ABSENT; it IS absent, so the control HELD.
        self.assertEqual(present("primes"),
                         gr["P_D2_primes_negative_control"]["primes_signature_present"])
        self.assertTrue(gr["P_D2_primes_negative_control"]["holds"])
        self.assertFalse(present("primes"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
