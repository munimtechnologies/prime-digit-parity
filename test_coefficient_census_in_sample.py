#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Produced with AI by MV2 of Munim, Inc.
"""In-sample coefficient census + internal-consistency audit of the density scaling law.

This is NOT a new sealed experiment. It is a self-contained, pure-standard-library re-analysis of
values ALREADY PUBLISHED in this repository, giving the single-prime density-law coefficient
    k(p) = (r(p) - r_control) / d(p),   d(p) = 1/p^2 - 1/p^3
an in-sample baseline over the twelve primes 2..37 that were measured one-at-a-time across Axes
D(exponent-cap)..D5, before Axis D8's out-of-sample dense sweep (primes 41..97) is graded.

It verifies three things, each re-derivable offline with no network and no large data:

  1. INTERNAL CONSISTENCY (exact). For p in {29, 31, 37} the exact full-precision r(p) and r_control
     stored in data/axis_d5_density_scaling_law.json, run back through the SAME density() formula the
     job used, must reproduce the k(p) published in that file to float precision. This catches any
     transcription drift between the stored ratios and the stored coefficients.

  2. ROUNDED-TABLE CONSISTENCY. The 3-significant-figure r(p) table quoted in SYNTHESIS.md must, for
     those same three primes, give a k(p) within the table's own rounding of the exact value.

  3. IN-SAMPLE CENSUS (honest, bounded). Across the eight primes 11..37 (where 1/p is small enough
     that the asymptotic density law applies), every k(p) lands in [7.5, 8.7], the spread is small,
     the OLS slope of k vs p is consistent with zero, and the median sits at ~8.1 -- i.e. the
     coefficient is ALREADY near-constant in-sample. For the small primes 3 and 5, k(p) > 9,
     recording honestly that the law is asymptotic with O(1/p) corrections at small p. This is an
     in-sample baseline only; it does NOT settle universality-vs-drift, which needs Axis D8's
     out-of-sample high-N sweep at larger primes.

Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.
"""
import math
import unittest


def density(p):
    """d(p) = 1/p^2 - 1/p^3, identical to the formula in the Axis D5/D8 jobs."""
    return 1.0 / (p * p) - 1.0 / (p * p * p)


def coefficient(r, r_control, p):
    return (r - r_control) / density(p)


def ols_slope(xs, ys):
    """Ordinary-least-squares slope of y = a + b*x (same estimator as the Axis D8 job)."""
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    return sxy / sxx


# ---- Exact values stored in data/axis_d5_density_scaling_law.json (full float precision) ----
R_CONTROL_EXACT = 0.014860313424210063
R_EXACT = {29: 0.02388594394194274, 31: 0.022755927382039735, 37: 0.020843655330602732}
K_PUBLISHED = {29: 7.8616465248922225, 31: 7.840607847256792, 37: 8.418728266236329}

# ---- 3-significant-figure r(p) table as quoted in SYNTHESIS.md, with r_control = 0.01486 ----
R_CONTROL_SYN = 0.014860
R_SYN = {
    2: 1.142, 3: 0.773, 5: 0.342, 7: 0.169,
    11: 0.0777, 13: 0.0603, 17: 0.0420, 19: 0.0366, 23: 0.0291,
    29: 0.0239, 31: 0.0228, 37: 0.0208,
}
LARGE_PRIMES = [11, 13, 17, 19, 23, 29, 31, 37]  # 1/p small: asymptotic regime


class CoefficientCensus(unittest.TestCase):
    def test_1_internal_consistency_exact(self):
        """The stored exact r(p)/r_control reproduce the stored k(p) to float precision."""
        for p in (29, 31, 37):
            k = coefficient(R_EXACT[p], R_CONTROL_EXACT, p)
            self.assertAlmostEqual(k, K_PUBLISHED[p], places=9,
                                   msg=f"k({p}) recomputed {k} != published {K_PUBLISHED[p]}")

    def test_2_rounded_table_consistency(self):
        """The 3-sig-fig SYNTHESIS r(p) table agrees with the exact k(p) within its own rounding.

        The tolerance is 0.08 (not tighter) on purpose: because r(p)-r_control is smallest for the
        largest primes, a half-ulp of 3-sig-fig rounding in r(37) already moves k(37) by ~0.06. That
        is itself the point -- rounding alone perturbs k by an amount comparable to a plausible drift,
        which is exactly why settling universality-vs-drift needs Axis D8's full-precision high-N sweep.
        """
        for p in (29, 31, 37):
            k_syn = coefficient(R_SYN[p], R_CONTROL_SYN, p)
            self.assertLess(abs(k_syn - K_PUBLISHED[p]), 0.08,
                            msg=f"SYNTHESIS k({p})={k_syn} vs exact {K_PUBLISHED[p]}")

    def test_3_large_prime_band(self):
        """In-sample: every k(p) for 11<=p<=37 lands inside [7.5, 8.7]."""
        for p in LARGE_PRIMES:
            k = coefficient(R_SYN[p], R_CONTROL_SYN, p)
            self.assertTrue(7.5 <= k <= 8.7, msg=f"k({p})={k} outside [7.5, 8.7]")

    def test_4_cluster_tight_and_centred(self):
        """The large-prime coefficients cluster tightly and are centred near 8.1."""
        ks = [coefficient(R_SYN[p], R_CONTROL_SYN, p) for p in LARGE_PRIMES]
        spread = max(ks) - min(ks)
        self.assertLess(spread, 0.6, msg=f"spread {spread} too large")
        median = sorted(ks)[len(ks) // 2 - 1:len(ks) // 2 + 1]
        med = sum(median) / 2.0
        self.assertTrue(7.8 <= med <= 8.4, msg=f"median k {med} outside [7.8, 8.4]")

    def test_5_slope_consistent_with_zero(self):
        """OLS slope of k vs p over 11..37 is small -- no monotone drift in-sample."""
        ks = [coefficient(R_SYN[p], R_CONTROL_SYN, p) for p in LARGE_PRIMES]
        b = ols_slope([float(p) for p in LARGE_PRIMES], ks)
        self.assertLess(abs(b), 0.05, msg=f"in-sample slope {b} not near zero")

    def test_6_small_prime_asymptotic_correction(self):
        """Honest limit: at small p (3, 5) the coefficient is clearly above 9 -- the law is asymptotic."""
        for p in (3, 5):
            k = coefficient(R_SYN[p], R_CONTROL_SYN, p)
            self.assertGreater(k, 9.0, msg=f"k({p})={k} not above 9")

    def test_7_density_formula_matches_job(self):
        """density() matches the closed form d(p) = (p-1)/p^3 used throughout the arc."""
        for p in (2, 3, 5, 7, 11, 13, 41, 97):
            self.assertAlmostEqual(density(p), (p - 1) / float(p ** 3), places=15)
        self.assertTrue(density(41) > density(43) > density(97))


if __name__ == "__main__":
    unittest.main(verbosity=2)
