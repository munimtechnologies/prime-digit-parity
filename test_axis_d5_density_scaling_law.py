# SPDX-License-Identifier: Apache-2.0
"""Pure-standard-library verification of the committed axis-d5-density-scaling-law result.

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
  5. re-applies the pre-registered P-D5-sf/ex29/ex31/ex37/odd/primes verdicts to the stored
     numbers and asserts the published outcome exactly as graded by the frozen detector:

       * P-D5-sf     (control, full squarefree A005117)       : predicted present -> HELD (present)
       * P-D5-ex29   (squarefree-except-29)                    : predicted present -> HELD (present)
       * P-D5-ex31   (squarefree-except-31)                    : predicted present -> HELD (present)
       * P-D5-ex37   (squarefree-except-37)                    : predicted present -> HELD (present)
       * P-D5-odd    (anchor, odd squarefree A056911)          : predicted absent  -> HELD (absent)
       * P-D5-primes (negative control, A000040)               : predicted absent  -> HELD (absent)

  6. and -- the new, out-of-sample test of this axis -- re-derives the density-scaling
     coefficient k(p) = (r(p) - r_control) / d(p) with d(p) = 1/p^2 - 1/p^3 from the stored
     ratios, checks it reproduces the committed coefficient_k_by_prime, and asserts that
     k(29), k(31), k(37) each land inside the pre-registered acceptance band [5, 11]
     (density_law_holds).

     This CONTINUES the Axis D3 crossover (v0.27.0) / Axis D4 convergence (v0.28.0) story.
     Axis D (v0.24.0) and Axis D2 (v0.26.0) found that relaxing the strict squarefree (mu^2)
     exponent cap at a single prime 2, 3, 5 or 7 destroys the base-2 gap-parity lag-2
     unique-quiet signature, the frozen ratio r(p) = |z(2)|/median|z(k!=2)| falling
     monotonically (1.142, 0.773, 0.342, 0.169) but staying above 0.1. Axis D3 found a
     CROSSOVER at 11/13 (r=0.0777, 0.0603, signature RETURNS); Axis D4 confirmed monotone
     CONVERGENCE at 17/19/23 (r=0.0420, 0.0366, 0.0291) with no second crossover. Axis D5
     extends the convergence to 29/31/37 (r=0.0239, 0.0228, 0.0208, all present) AND tests a
     quantitative DENSITY SCALING LAW r(p) - r_control ~ C*(1/p^2 - 1/p^3), C ~ 8, estimated
     from p <= 23 and confirmed out-of-sample: k(29)=7.86, k(31)=7.84, k(37)=8.42, all in
     [5, 11]. Both convergence_to_control_confirmed and density_law_holds are True.

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
DATA = os.path.join(HERE, "data", "axis_d5_density_scaling_law.json")

STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]
SIEVES = ["squarefree", "squarefree_except_29", "squarefree_except_31",
          "squarefree_except_37", "odd_squarefree", "primes"]
LAGS = [1, 2, 3, 4, 5, 6, 7, 8]
RELAXED_PRIMES = ["2", "3", "5", "7", "11", "13", "17", "19", "23", "29", "31", "37"]


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
        # Past the Axis D3 crossover (11) and the Axis D4 convergence (17/19/23), does the frozen
        # ratio r(p) keep falling toward the squarefree control at 29/31/37? It does: all three stay
        # BELOW the 0.1 threshold, so the lag-2 unique-quiet signature stays PRESENT for all three.
        # Control stays present; anchor (odd squarefree) and the negative control (primes) stay
        # absent. The pre-registered convergence prediction HOLDS; no second crossover appears.
        g = self.d["grading"]["signature_all_sieves"]
        self.assertTrue(g["squarefree"]["signature_present"])
        self.assertTrue(g["squarefree_except_29"]["signature_present"])
        self.assertTrue(g["squarefree_except_31"]["signature_present"])
        self.assertTrue(g["squarefree_except_37"]["signature_present"])
        self.assertFalse(g["odd_squarefree"]["signature_present"])
        self.assertFalse(g["primes"]["signature_present"])
        gr = self.d["grading"]
        self.assertTrue(gr["convergence_to_control_confirmed"])
        self.assertFalse(gr["second_crossover_or_redestruction"])
        self.assertTrue(gr["monotone_trend_continues"])
        self.assertTrue(gr["all_three_below_threshold"])
        self.assertLessEqual(g["squarefree_except_29"]["ratio_lag2_over_median_others"], 0.1)
        self.assertLessEqual(g["squarefree_except_31"]["ratio_lag2_over_median_others"], 0.1)
        self.assertLessEqual(g["squarefree_except_37"]["ratio_lag2_over_median_others"], 0.1)

    def test_monotone_trend(self):
        # r(p) falls monotonically across every relaxed prime 2 < 3 < 5 < 7 < 11 < ... < 37,
        # continuing the Axis D2/D3/D4 trend down through the crossover and on toward the control.
        rbp = self.d["grading"]["ratio_by_relaxed_prime"]
        ordered = [rbp[k] for k in RELAXED_PRIMES]
        for a, b in zip(ordered, ordered[1:]):
            self.assertGreater(a, b)
        # The committed Axis D2/D3/D4 anchor values are preserved, and 29/31/37 land below r(23).
        self.assertAlmostEqual(rbp["5"], 0.342, places=3)
        self.assertAlmostEqual(rbp["7"], 0.169, places=3)
        self.assertLess(rbp["11"], 0.1)
        self.assertLess(rbp["23"], rbp["19"])
        self.assertLess(rbp["29"], rbp["23"])
        self.assertLess(rbp["31"], rbp["29"])
        self.assertLess(rbp["37"], rbp["31"])

    def test_density_scaling_law(self):
        # The new, predictive test. The super-control excess r(p) - r_control is modelled as
        # proportional to the admitted-class density d(p) = 1/p^2 - 1/p^3, i.e. r(p) - r_control
        # ~ C*(1/p^2 - 1/p^3). Re-derive the dimensionless coefficient k(p) = (r(p) - r_control)/d(p)
        # from the stored ratios, check it matches the committed coefficient_k_by_prime, and assert
        # the three OUT-OF-SAMPLE primes 29, 31, 37 each land in the pre-registered band [5, 11].
        dsl = self.d["grading"]["density_scaling_law"]
        rbp = self.d["grading"]["ratio_by_relaxed_prime"]
        r_control = dsl["r_control_measured"]
        # r_control must equal the measured control ratio (not a hard-coded constant).
        self.assertAlmostEqual(
            r_control,
            self.d["grading"]["signature_all_sieves"]["squarefree"]["ratio_lag2_over_median_others"],
            places=12)
        band_lo, band_hi = dsl["preregistered_band"]
        self.assertEqual([band_lo, band_hi], [5.0, 11.0])
        for p_str in RELAXED_PRIMES:
            p = int(p_str)
            d_p = density(p)
            self.assertAlmostEqual(d_p, dsl["density_by_prime"][p_str], places=12, msg=p_str)
            k_p = (rbp[p_str] - r_control) / d_p
            self.assertAlmostEqual(k_p, dsl["coefficient_k_by_prime"][p_str], places=6, msg=p_str)
        # Out-of-sample acceptance: k(29), k(31), k(37) all inside the pre-registered band.
        for p_str in ["29", "31", "37"]:
            k_p = (rbp[p_str] - r_control) / density(int(p_str))
            self.assertGreaterEqual(k_p, band_lo, p_str)
            self.assertLessEqual(k_p, band_hi, p_str)
        self.assertEqual(dsl["k_29_31_37"],
                         [dsl["coefficient_k_by_prime"][p] for p in ["29", "31", "37"]])
        self.assertTrue(dsl["density_law_holds"])

    def test_regrade_matches_published_verdicts(self):
        gr = self.d["grading"]
        ps = self.d["per_sieve"]

        def present(sv):
            p, _, _, _ = lag2_signature(ps[sv]["x_base2"]["lags"])
            return p

        # Control: full squarefree predicted PRESENT; it IS present, so the prediction HELD.
        self.assertEqual(present("squarefree"),
                         gr["P_D5_sf_control"]["squarefree_signature_present"])
        self.assertTrue(gr["P_D5_sf_control"]["holds"])
        self.assertTrue(present("squarefree"))

        # squarefree-except-29 predicted PRESENT (convergence continues); it IS present -> HELD.
        self.assertEqual(present("squarefree_except_29"),
                         gr["P_D5_ex29"]["squarefree_except_29_signature_present"])
        self.assertTrue(gr["P_D5_ex29"]["holds"])
        self.assertTrue(present("squarefree_except_29"))

        # squarefree-except-31 predicted PRESENT; it IS present -> HELD.
        self.assertEqual(present("squarefree_except_31"),
                         gr["P_D5_ex31"]["squarefree_except_31_signature_present"])
        self.assertTrue(gr["P_D5_ex31"]["holds"])
        self.assertTrue(present("squarefree_except_31"))

        # squarefree-except-37 predicted PRESENT; it IS present -> HELD.
        self.assertEqual(present("squarefree_except_37"),
                         gr["P_D5_ex37"]["squarefree_except_37_signature_present"])
        self.assertTrue(gr["P_D5_ex37"]["holds"])
        self.assertTrue(present("squarefree_except_37"))

        # Anchor: odd squarefree predicted ABSENT; it IS absent, so the anchor HELD.
        self.assertEqual(present("odd_squarefree"),
                         gr["P_D5_odd_anchor"]["odd_squarefree_signature_present"])
        self.assertTrue(gr["P_D5_odd_anchor"]["holds"])
        self.assertFalse(present("odd_squarefree"))

        # Negative control: primes predicted ABSENT; it IS absent, so the control HELD.
        self.assertEqual(present("primes"),
                         gr["P_D5_primes_negative_control"]["primes_signature_present"])
        self.assertTrue(gr["P_D5_primes_negative_control"]["holds"])
        self.assertFalse(present("primes"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
