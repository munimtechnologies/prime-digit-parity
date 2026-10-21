# SPDX-License-Identifier: Apache-2.0
"""Pure-standard-library verification of the committed axis-d6-density-additivity result.

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
  5. re-applies the pre-registered P-D6-sf/pair-17_19/pair-23_29/pair-31_37/odd/primes
     verdicts to the stored numbers and asserts the published outcome exactly as graded by
     the frozen detector:

       * P-D6-sf          (control, full squarefree A005117)    : predicted present -> HELD
       * P-D6-pair-17_19  (squarefree-except-{17,19})           : predicted present -> HELD
       * P-D6-pair-23_29  (squarefree-except-{23,29})           : predicted present -> HELD
       * P-D6-pair-31_37  (squarefree-except-{31,37})           : predicted present -> HELD
       * P-D6-odd         (anchor, odd squarefree A056911)      : predicted absent  -> HELD
       * P-D6-primes      (negative control, A000040)           : predicted absent  -> HELD

  6. and -- the new, out-of-sample test of this axis -- re-derives the density-ADDITIVITY
     coefficient k({p,q}) = (r({p,q}) - r_control) / (d(p) + d(q)) with d(p) = 1/p^2 - 1/p^3
     from the stored ratios, checks it reproduces the committed coefficient_k_by_pair, and
     asserts that k({17,19}), k({23,29}), k({31,37}) each land inside the pre-registered
     acceptance band [5, 11] (additivity_law_holds). A pass is positive evidence that two
     independent single-prime relaxations perturb the signature ADDITIVELY in the admitted
     density -- the first evidence that the Axis-D5 density scaling law composes linearly.

     This CONTINUES the Axis D3 crossover (v0.27.0) / Axis D4 convergence (v0.28.0) /
     Axis D5 single-prime density scaling law (v0.29.0) story. Axis D5 found that relaxing
     the strict squarefree (mu^2) exponent cap at a SINGLE prime p perturbs the frozen ratio
     r(p) = |z(2)|/median|z(k!=2)| above the squarefree control by C*(1/p^2 - 1/p^3) with
     C ~ 8. Axis D6 relaxes the cap at TWO primes at once (pairs {17,19}, {23,29}, {31,37})
     and tests whether the two-prime excess is the SUM of the single-prime excesses:
     r({p,q}) - r_control ~ C*(d(p) + d(q)). All three pairs stay present (r below 0.1) and
     k({p,q}) = 8.33, 7.91, 8.06 -- all in [5, 11] and clustered on the same C ~ 8 as the
     single-prime law. convergence_holds and additivity_law_holds are both True.

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
DATA = os.path.join(HERE, "data", "axis_d6_density_additivity.json")

STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]
SIEVES = ["squarefree", "squarefree_except_17_19", "squarefree_except_23_29",
          "squarefree_except_31_37", "odd_squarefree", "primes"]
LAGS = [1, 2, 3, 4, 5, 6, 7, 8]
# (sieve key, (p, q)) for the three two-prime relaxations.
PAIRS = {"17_19": (17, 19), "23_29": (23, 29), "31_37": (31, 37)}


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


def density(p):
    """Density of the admitted extra class {p^2 * m : p does not divide m}: 1/p^2 - 1/p^3."""
    return 1.0 / (p * p) - 1.0 / (p * p * p)


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
        # Relaxing the strict squarefree exponent cap at TWO large primes at once keeps the
        # base-2 gap-parity lag-2 unique-quiet signature PRESENT for all three pairs: every
        # pair's ratio stays BELOW the 0.1 threshold. Control stays present; the anchor (odd
        # squarefree) and the negative control (primes) stay absent. No pair is destroyed.
        g = self.d["grading"]["signature_all_sieves"]
        self.assertTrue(g["squarefree"]["signature_present"])
        for key in PAIRS:
            sv = "squarefree_except_" + key
            self.assertTrue(g[sv]["signature_present"], sv)
            self.assertLessEqual(g[sv]["ratio_lag2_over_median_others"], 0.1, sv)
        self.assertFalse(g["odd_squarefree"]["signature_present"])
        self.assertFalse(g["primes"]["signature_present"])
        gr = self.d["grading"]
        self.assertTrue(gr["convergence_holds"])
        self.assertTrue(gr["all_pairs_present"])
        self.assertTrue(gr["all_pairs_below_threshold"])
        self.assertFalse(gr["some_pair_destroyed"])
        self.assertEqual(gr["threshold"], 0.1)

    def test_density_additivity_law(self):
        # The new, predictive test. If two independent single-prime relaxations perturb the
        # signature additively in the admitted density, the two-prime super-control excess
        # r({p,q}) - r_control is proportional to the SUM of the single-prime densities,
        # d(p) + d(q) with d(p) = 1/p^2 - 1/p^3. Re-derive the dimensionless coefficient
        # k({p,q}) = (r({p,q}) - r_control)/(d(p)+d(q)) from the stored ratios, check it
        # matches the committed coefficient_k_by_pair, and assert all three pairs land in the
        # pre-registered band [5, 11] -- the same band validated for the single-prime law.
        dal = self.d["grading"]["density_additivity_law"]
        sas = self.d["grading"]["signature_all_sieves"]
        r_control = dal["r_control_measured"]
        # r_control must equal the measured control ratio (not a hard-coded constant).
        self.assertAlmostEqual(
            r_control, sas["squarefree"]["ratio_lag2_over_median_others"], places=12)
        band_lo, band_hi = dal["preregistered_band"]
        self.assertEqual([band_lo, band_hi], [5.0, 11.0])
        for key, (p, q) in PAIRS.items():
            d_sum = density(p) + density(q)
            self.assertAlmostEqual(d_sum, dal["summed_density_by_pair"][key], places=12, msg=key)
            r_pair = sas["squarefree_except_" + key]["ratio_lag2_over_median_others"]
            self.assertAlmostEqual(r_pair, dal["ratio_by_pair"][key], places=12, msg=key)
            k_pair = (r_pair - r_control) / d_sum
            self.assertAlmostEqual(k_pair, dal["coefficient_k_by_pair"][key], places=6, msg=key)
            # Out-of-sample acceptance: each pair's coefficient inside the pre-registered band.
            self.assertGreaterEqual(k_pair, band_lo, key)
            self.assertLessEqual(k_pair, band_hi, key)
        self.assertTrue(dal["additivity_law_holds"])

    def test_regrade_matches_published_verdicts(self):
        gr = self.d["grading"]
        ps = self.d["per_sieve"]

        def present(sv):
            p, _, _, _ = lag2_signature(ps[sv]["x_base2"]["lags"])
            return p

        # Control: full squarefree predicted PRESENT; it IS present, so the prediction HELD.
        self.assertEqual(present("squarefree"),
                         gr["P_D6_sf_control"]["squarefree_signature_present"])
        self.assertTrue(gr["P_D6_sf_control"]["holds"])
        self.assertTrue(present("squarefree"))

        # Each two-prime relaxation predicted PRESENT (additive convergence); each IS present.
        for key in PAIRS:
            sv = "squarefree_except_" + key
            gkey = "P_D6_pair_" + key
            self.assertEqual(present(sv), gr[gkey][sv + "_signature_present"])
            self.assertTrue(gr[gkey]["holds"], gkey)
            self.assertTrue(present(sv), sv)

        # Anchor: odd squarefree predicted ABSENT; it IS absent, so the anchor HELD.
        self.assertEqual(present("odd_squarefree"),
                         gr["P_D6_odd_anchor"]["odd_squarefree_signature_present"])
        self.assertTrue(gr["P_D6_odd_anchor"]["holds"])
        self.assertFalse(present("odd_squarefree"))

        # Negative control: primes predicted ABSENT; it IS absent, so the control HELD.
        self.assertEqual(present("primes"),
                         gr["P_D6_primes_negative_control"]["primes_signature_present"])
        self.assertTrue(gr["P_D6_primes_negative_control"]["holds"])
        self.assertFalse(present("primes"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
