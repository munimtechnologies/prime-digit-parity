# SPDX-License-Identifier: Apache-2.0
"""Pure-standard-library verification of the committed axis-d-exponent-cap-lag-signature result.

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
  5. re-applies the pre-registered P-D-sf/ex2/ex3/odd/primes verdicts to the stored numbers
     and asserts the published outcome exactly as graded by the frozen detector:

       * P-D-sf   (control, full squarefree A005117)   : predicted present -> HELD   (present)
       * P-D-ex2  (PRIMARY, squarefree-except-2)        : predicted present -> FALSIFIED (absent)
       * P-D-ex3  (replication, squarefree-except-3)    : predicted present -> FALSIFIED (absent)
       * P-D-odd  (anchor, odd squarefree A056911)      : predicted absent  -> HELD   (absent)
       * P-D-primes (negative control, A000040)         : predicted absent  -> HELD   (absent)

     The PRIMARY prediction was FALSIFIED: relaxing the squarefree exponent cap at a SINGLE
     small prime -- squarefree-except-2 (allow 4|n but not 8|n) and squarefree-except-3 (allow
     9|n but not 27|n) -- keeps every even integer and every odd gap, yet the lag-2 unique-quiet
     signature is GONE (quiet lag 5 / 8, ratio 1.14 / 0.77). Both sieves still admit even
     integers, so even-admission (shown NECESSARY in the odd-squarefree anchor) is confirmed
     NOT SUFFICIENT from the enlarging side too. Together with the Axis-C thinning result, the
     signature tolerates NEITHER thinning a congruence class NOR enlarging the exponent cap: it
     needs EXACTLY the strict squarefree (mu^2) congruence structure.
     `signature_tolerates_exponent_cap_relaxation` is accordingly False.

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
DATA = os.path.join(HERE, "data", "axis_d_exponent_cap_lag_signature.json")

STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]
SIEVES = ["squarefree", "squarefree_except_2", "squarefree_except_3",
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
        # The point of the study and its FALSIFICATION: the lag-2 unique-quiet signature is
        # present for the full strict squarefree sieve, but ABSENT for every single-prime
        # exponent-cap relaxation (squarefree-except-2 allows 4|n; squarefree-except-3 allows
        # 9|n) even though both KEEP every even integer and every odd gap. So relaxing the mu^2
        # cap at even one prime destroys it, just as thinning a congruence class did in Axis C.
        # Even-admission is confirmed necessary but not sufficient from the enlarging side too;
        # the signature needs EXACTLY the strict squarefree congruence structure. The "tolerates
        # exponent-cap relaxation" flag is therefore False.
        g = self.d["grading"]["signature_all_sieves"]
        self.assertTrue(g["squarefree"]["signature_present"])
        self.assertFalse(g["squarefree_except_2"]["signature_present"])
        self.assertFalse(g["squarefree_except_3"]["signature_present"])
        self.assertFalse(g["odd_squarefree"]["signature_present"])
        self.assertFalse(g["primes"]["signature_present"])
        self.assertFalse(self.d["grading"]["signature_tolerates_exponent_cap_relaxation"])

    def test_regrade_matches_published_verdicts(self):
        gr = self.d["grading"]
        ps = self.d["per_sieve"]

        def present(sv):
            p, _, _, _ = lag2_signature(ps[sv]["x_base2"]["lags"])
            return p

        # Control: full squarefree predicted PRESENT; it IS present, so the prediction HELD.
        self.assertEqual(present("squarefree"),
                         gr["P_D_sf_control"]["squarefree_signature_present"])
        self.assertTrue(gr["P_D_sf_control"]["holds"])
        self.assertTrue(present("squarefree"))

        # PRIMARY: squarefree-except-2 predicted PRESENT (relax the cap at 2, keep evens); it is
        # ABSENT, so the primary prediction was FALSIFIED (holds == False).
        self.assertEqual(present("squarefree_except_2"),
                         gr["P_D_ex2_primary"]["squarefree_except_2_signature_present"])
        self.assertFalse(gr["P_D_ex2_primary"]["holds"])
        self.assertFalse(present("squarefree_except_2"))

        # Replication: squarefree-except-3 (relax the cap at an ODD prime) predicted PRESENT;
        # it is ABSENT -> FALSIFIED.
        self.assertEqual(present("squarefree_except_3"),
                         gr["P_D_ex3_replication"]["squarefree_except_3_signature_present"])
        self.assertFalse(gr["P_D_ex3_replication"]["holds"])
        self.assertFalse(present("squarefree_except_3"))

        # Anchor: odd squarefree predicted ABSENT; it IS absent, so the anchor HELD.
        self.assertEqual(present("odd_squarefree"),
                         gr["P_D_odd_anchor"]["odd_squarefree_signature_present"])
        self.assertTrue(gr["P_D_odd_anchor"]["holds"])
        self.assertFalse(present("odd_squarefree"))

        # Negative control: primes predicted ABSENT; it IS absent, so the control HELD.
        self.assertEqual(present("primes"),
                         gr["P_D_primes_negative_control"]["primes_signature_present"])
        self.assertTrue(gr["P_D_primes_negative_control"]["holds"])
        self.assertFalse(present("primes"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
