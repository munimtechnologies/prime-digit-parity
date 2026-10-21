"""Proof-by-enumeration that the D9 density-only prediction rho = (p+1)/p is the
admitted-class density ratio fixed by the sieve, and nothing fitted.

Self-contained, pure-stdlib verifier. No network, no sealed compute job. It
grounds, in exact arithmetic, what Axis D9 (cap-2 vs cap-3, queued) WOULD measure,
so the density-only null has a proven baseline even while D9's job has not run.

Mechanism (read from this repo's code, not assumed):
  * A member n of the strict squarefree (mu^2) sieve satisfies
    v_p(n) <= caps.get(p, 1) for every prime p, with the default cap 1.
  * "Relaxing the exponent cap at p to c" sets caps[p] = c and forbids p^(c+1).
    - cap 2 forbids p^3, so relative to the cap-1 control it newly admits
      EXACTLY {n : v_p(n) == 2}, density d2(p) = 1/p^2 - 1/p^3.
    - cap 3 forbids p^4, so relative to the same control it newly admits
      EXACTLY {n : v_p(n) in {2, 3}}, density d3(p) = 1/p^2 - 1/p^4.
  * Axis D9 measures r_cap2 - r_control and r_cap3 - r_control and forms
    rho(p) = (r_cap3 - r_control) / (r_cap2 - r_control). Under the density-only
    null (response proportional to newly-admitted density) rho(p) = d3/d2.

This test proves, by exact counting and exact rational arithmetic, five things:

  (1) Counting identity (cap 2): #{n<=N : v_p(n)==2} == floor(N/p^2)-floor(N/p^3).
  (2) Counting identity (cap 3): #{n<=N : v_p(n) in {2,3}}
                                   == floor(N/p^2) - floor(N/p^4).
  (3) The v_p==3 class alone has density exactly 1/p^3 - 1/p^4 = (p-1)/p^4, and
      d3(p) = d2(p) + (1/p^3 - 1/p^4) -- the cap-3 admission is the cap-2 class
      PLUS the v_p==3 class, with no overlap.
  (4) EXACT ratio identity (Fraction, no float): d3(p)/d2(p) == (p+1)/p for every
      prime p. Hence the density-only null for D9 is rho(p) = (p+1)/p exactly, and
      the repo's structure statistic C3/C2 = p*(rho-1) equals exactly 1 under it.
  (5) Density limits: count/N -> d2, d3 with |count/N - density| < (#terms)/N,
      an a-priori floor bound (proven arithmetic, not a fitted tolerance).

Honest limit: these are elementary identities. The point is that the D9 null is
the admitted-class density ratio fixed by the sieve definition, so if D9 ever
runs, any rho(p) that deviates from (p+1)/p beyond resolution is structure, not
density -- the decision rule is pinned to proven arithmetic, not a guessed curve.
This grades no live data; it only fixes the baseline. Pre-registered constants and
tolerances are fixed below, before any count is taken.

Produced by AI agents working as MV2, Munim, Inc. Apache-2.0.
"""

import unittest
from fractions import Fraction

# --- pre-registered constants (fixed before measurement) ---------------------
D5_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
D9_PRIMES = [13, 17, 19, 23, 29, 31, 37]   # Axis D9's cap-2 vs cap-3 sweep range
OUT_OF_SAMPLE_PRIMES = [41, 43, 47]
ALL_PRIMES = D5_PRIMES + OUT_OF_SAMPLE_PRIMES

N_ENUM = 200_000        # window for the exact enumeration <-> formula cross-check
N_DENSITY = 5_000_000   # window for the density-limit check (formula only)

DENSITY_LIMIT_SLACK = 1e-15   # float slack on top of the exact floor bound
RATIO_ID_PLACES = 12          # float cross-check of the (p+1)/p ratio


def v_p(n, p):
    """p-adic valuation of n (n >= 1)."""
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def d2(p):
    """Density newly admitted at cap 2 (vs cap 1): {v_p == 2}."""
    return 1.0 / p ** 2 - 1.0 / p ** 3


def d3(p):
    """Density newly admitted at cap 3 (vs cap 1): {v_p in {2, 3}}."""
    return 1.0 / p ** 2 - 1.0 / p ** 4


def d2_exact(p):
    return Fraction(1, p ** 2) - Fraction(1, p ** 3)


def d3_exact(p):
    return Fraction(1, p ** 2) - Fraction(1, p ** 4)


class TestCountingIdentities(unittest.TestCase):
    def test_cap2_enumeration_matches_floor_formula(self):
        # (1) #{n<=N : v_p(n)==2} == floor(N/p^2) - floor(N/p^3)
        for p in ALL_PRIMES:
            enum = sum(1 for n in range(1, N_ENUM + 1) if v_p(n, p) == 2)
            formula = N_ENUM // p ** 2 - N_ENUM // p ** 3
            self.assertEqual(
                enum, formula,
                msg="cap-2 counting identity failed at p={}: {} vs {}".format(
                    p, enum, formula))

    def test_cap3_enumeration_matches_floor_formula(self):
        # (2) #{n<=N : v_p(n) in {2,3}} == floor(N/p^2) - floor(N/p^4)
        for p in ALL_PRIMES:
            enum = sum(1 for n in range(1, N_ENUM + 1) if v_p(n, p) in (2, 3))
            formula = N_ENUM // p ** 2 - N_ENUM // p ** 4
            self.assertEqual(
                enum, formula,
                msg="cap-3 counting identity failed at p={}: {} vs {}".format(
                    p, enum, formula))


class TestClassDecomposition(unittest.TestCase):
    def test_vp3_class_density_and_disjoint_sum(self):
        # (3) {v_p==3} has density (p-1)/p^4 and d3 == d2 + that, exactly.
        for p in ALL_PRIMES:
            vp3 = Fraction(1, p ** 3) - Fraction(1, p ** 4)
            self.assertEqual(
                vp3, Fraction(p - 1, p ** 4),
                msg="v_p==3 density wrong at p={}".format(p))
            self.assertEqual(
                d3_exact(p), d2_exact(p) + vp3,
                msg="cap-3 class != cap-2 class + v_p==3 class at p={}".format(p))

    def test_vp3_enumeration_matches_formula(self):
        # The v_p==3 class count == floor(N/p^3) - floor(N/p^4) exactly.
        for p in ALL_PRIMES:
            enum = sum(1 for n in range(1, N_ENUM + 1) if v_p(n, p) == 3)
            formula = N_ENUM // p ** 3 - N_ENUM // p ** 4
            self.assertEqual(
                enum, formula,
                msg="v_p==3 counting identity failed at p={}: {} vs {}".format(
                    p, enum, formula))


class TestDensityRatioIsExact(unittest.TestCase):
    def test_ratio_equals_p_plus_one_over_p_exactly(self):
        # (4) EXACT: d3(p)/d2(p) == (p+1)/p for every prime p (Fraction, no float).
        for p in ALL_PRIMES:
            ratio = d3_exact(p) / d2_exact(p)
            self.assertEqual(
                ratio, Fraction(p + 1, p),
                msg="exact ratio != (p+1)/p at p={}: got {}".format(p, ratio))

    def test_structure_statistic_is_one_under_density_null(self):
        # (4 cont.) D9's C3/C2 = p*(rho-1); with rho=(p+1)/p this is exactly 1.
        for p in D9_PRIMES:
            rho = Fraction(p + 1, p)
            self.assertEqual(
                p * (rho - 1), Fraction(1, 1),
                msg="density-null C3/C2 != 1 at p={}".format(p))

    def test_ratio_float_crosscheck(self):
        # Float sanity check of the same ratio to many places.
        for p in ALL_PRIMES:
            self.assertAlmostEqual(
                d3(p) / d2(p), (p + 1) / p, places=RATIO_ID_PLACES,
                msg="float ratio mismatch at p={}".format(p))


class TestDensityLimits(unittest.TestCase):
    def test_cap2_and_cap3_counts_converge(self):
        # (5) |count/N - density| < (#floor terms)/N, an a-priori bound.
        for p in ALL_PRIMES:
            c2 = N_DENSITY // p ** 2 - N_DENSITY // p ** 3
            c3 = N_DENSITY // p ** 2 - N_DENSITY // p ** 4
            self.assertLess(
                abs(c2 / N_DENSITY - d2(p)), 2.0 / N_DENSITY + DENSITY_LIMIT_SLACK,
                msg="cap-2 density limit failed at p={}".format(p))
            self.assertLess(
                abs(c3 / N_DENSITY - d3(p)), 2.0 / N_DENSITY + DENSITY_LIMIT_SLACK,
                msg="cap-3 density limit failed at p={}".format(p))


if __name__ == "__main__":
    for p in D9_PRIMES:
        print("p={:>2}  d2={:.10f}  d3={:.10f}  rho=d3/d2={}  (p+1)/p={:.10f}".format(
            p, d2(p), d3(p), d3_exact(p) / d2_exact(p), (p + 1) / p))
    unittest.main(verbosity=2)
