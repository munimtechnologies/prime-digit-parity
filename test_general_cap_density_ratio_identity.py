"""Proof-by-enumeration that the Axis-D9 density-only ratio generalises to every
relaxation depth as a truncated geometric series: rho_c(p) = sum_{j=0}^{c-2} p^-j.

Self-contained, pure-stdlib verifier. No network, no sealed compute job. It
extends the cap-3 result (rho_3 = (p+1)/p, shipped the previous run) to an
arbitrary exponent cap c, so the density-only null baseline for ANY cap-c vs
cap-2 relaxation is pinned in exact arithmetic, not a guessed curve.

Mechanism (read from this repo's code, not assumed):
  * A member n of the strict squarefree (mu^2) sieve satisfies
    v_p(n) <= caps.get(p, 1) for every prime p, with the default cap 1.
  * "Relaxing the exponent cap at p to c" sets caps[p] = c and forbids p^(c+1).
    Relative to the cap-1 control it newly admits EXACTLY {n : 2 <= v_p(n) <= c},
    the union of the disjoint classes {v_p == 2}, ..., {v_p == c}. Its natural
    density is
        d_c(p) = density{v_p >= 2} - density{v_p >= c+1} = 1/p^2 - 1/p^(c+1).
    (cap 2 -> 1/p^2 - 1/p^3; cap 3 -> 1/p^2 - 1/p^4, recovering last run.)
  * Axis D9 compares a cap-c relaxation with the cap-2 relaxation and, under the
    density-only null (response proportional to newly-admitted density), predicts
        rho_c(p) = d_c(p) / d_2(p).

This test proves, by exact counting and exact rational arithmetic, that:

  (1) Counting identity: #{n<=N : 2 <= v_p(n) <= c} == floor(N/p^2)-floor(N/p^(c+1)).
  (2) Incremental class {v_p == c} has density exactly (p-1)/p^(c+1), count
      floor(N/p^c) - floor(N/p^(c+1)), and d_c == d_{c-1} + that (telescoping).
  (3) EXACT closed form (Fraction, no float): the two algebraically-distinct forms
          rho_c(p) = (p^(c-1) - 1) / (p^(c-2) * (p-1))      [ratio of densities ]
          rho_c(p) = sum_{j=0}^{c-2} p^-j                    [truncated geometric]
      are equal, and for c == 3 both equal (p+1)/p (the previous run's result).
  (4) Generalised structure statistic: the cap-(c-1) -> cap-c increment satisfies
          p^(c-2) * (rho_c - rho_{c-1}) == 1   exactly, for every c >= 3,
      recovering the repo's C3/C2 = p*(rho_3 - 1) = 1 at c == 3.
  (5) Monotone bounded limit: rho_c is strictly increasing in c, bounded above by
      p/(p-1), with EXACT tail  p/(p-1) - rho_c == 1 / (p^(c-2) * (p-1)).
  (6) Density limits: count/N -> d_c with |count/N - d_c| < 2/N, an a-priori floor
      bound (proven arithmetic, not a fitted tolerance), plus a 12-place float
      cross-check of each rho_c.

Honest limit: these are elementary p-adic valuation identities. The point is not
the arithmetic but that the ENTIRE family of D9-style density nulls -- one per cap
depth -- is now fixed by the sieve definition as a geometric series in 1/p, so if a
cap-c experiment ever runs, any rho_c that deviates from sum_{j=0}^{c-2} p^-j beyond
resolution is structure, not density. This grades no live data; it only fixes the
baselines. Pre-registered constants and tolerances are fixed below, before any
count is taken.

Produced by AI agents working as MV2, Munim, Inc. Apache-2.0.
"""

import unittest
from fractions import Fraction

# --- pre-registered constants (fixed before measurement) ---------------------
D5_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
OUT_OF_SAMPLE_PRIMES = [41, 43, 47]
ALL_PRIMES = D5_PRIMES + OUT_OF_SAMPLE_PRIMES

CAPS = [2, 3, 4, 5, 6]        # relaxation depths to pin; cap c forbids p^(c+1)
CMAX = max(CAPS)

N_ENUM = 200_000              # window for exact enumeration <-> formula cross-check
N_DENSITY = 5_000_000         # window for the density-limit check (formula only)

DENSITY_LIMIT_SLACK = 1e-15   # float slack on top of the exact floor bound
RATIO_FLOAT_PLACES = 12       # float cross-check of each rho_c


def v_p(n, p):
    """p-adic valuation of n (n >= 1)."""
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def d_exact(p, c):
    """Density newly admitted at cap c (vs cap 1): {n : 2 <= v_p(n) <= c}."""
    return Fraction(1, p ** 2) - Fraction(1, p ** (c + 1))


def rho_ratio_form(p, c):
    """rho_c as the ratio of densities: (p^(c-1) - 1) / (p^(c-2) * (p-1))."""
    return Fraction(p ** (c - 1) - 1, p ** (c - 2) * (p - 1))


def rho_geometric_form(p, c):
    """rho_c as the truncated geometric series: sum_{j=0}^{c-2} p^-j."""
    return sum((Fraction(1, p ** j) for j in range(0, c - 1)), Fraction(0))


def valuation_histogram(p, N, top):
    """counts[j] = #{1 <= n <= N : v_p(n) == j} for j in 0..top (exact)."""
    counts = [0] * (top + 1)
    for n in range(1, N + 1):
        e = v_p(n, p)
        if e > top:
            e = top
        counts[e] += 1
    return counts


class TestCountingIdentity(unittest.TestCase):
    def test_cap_c_enumeration_matches_floor_formula(self):
        # (1) #{n<=N : 2 <= v_p(n) <= c} == floor(N/p^2) - floor(N/p^(c+1))
        for p in ALL_PRIMES:
            counts = valuation_histogram(p, N_ENUM, CMAX + 1)
            for c in CAPS:
                enum = sum(counts[j] for j in range(2, c + 1))
                formula = N_ENUM // p ** 2 - N_ENUM // p ** (c + 1)
                self.assertEqual(
                    enum, formula,
                    msg="cap-{} counting identity failed at p={}: {} vs {}".format(
                        c, p, enum, formula))


class TestIncrementalClass(unittest.TestCase):
    def test_vp_equals_c_density_count_and_telescoping(self):
        # (2) {v_p == c} has density (p-1)/p^(c+1), count floor(N/p^c)-floor(N/p^(c+1)),
        #     and d_c == d_{c-1} + that, exactly.
        for p in ALL_PRIMES:
            counts = valuation_histogram(p, N_ENUM, CMAX + 1)
            for c in CAPS:
                self.assertEqual(
                    counts[c], N_ENUM // p ** c - N_ENUM // p ** (c + 1),
                    msg="{{v_p==c}} count identity failed at p={}, c={}".format(p, c))
                inc = Fraction(1, p ** c) - Fraction(1, p ** (c + 1))
                self.assertEqual(
                    inc, Fraction(p - 1, p ** (c + 1)),
                    msg="{{v_p==c}} density wrong at p={}, c={}".format(p, c))
                self.assertEqual(
                    d_exact(p, c), d_exact(p, c - 1) + inc,
                    msg="telescoping d_c != d_(c-1)+inc at p={}, c={}".format(p, c))


class TestExactClosedForm(unittest.TestCase):
    def test_two_forms_agree_and_equal_density_ratio(self):
        # (3) rho_c = d_c/d_2 = ratio form = geometric form, all exact.
        for p in ALL_PRIMES:
            d2 = d_exact(p, 2)
            for c in CAPS:
                rho_from_density = d_exact(p, c) / d2
                self.assertEqual(
                    rho_from_density, rho_ratio_form(p, c),
                    msg="d_c/d_2 != ratio form at p={}, c={}".format(p, c))
                self.assertEqual(
                    rho_ratio_form(p, c), rho_geometric_form(p, c),
                    msg="ratio form != geometric form at p={}, c={}".format(p, c))

    def test_c3_recovers_previous_run(self):
        # (3 cont.) rho_3 == (p+1)/p, the cap-3 result shipped last run.
        for p in ALL_PRIMES:
            self.assertEqual(
                rho_geometric_form(p, 3), Fraction(p + 1, p),
                msg="rho_3 != (p+1)/p at p={}".format(p))


class TestStructureStatistic(unittest.TestCase):
    def test_increment_structure_statistic_is_one(self):
        # (4) p^(c-2) * (rho_c - rho_{c-1}) == 1 exactly, for every c >= 3.
        for p in ALL_PRIMES:
            for c in CAPS:
                if c < 3:
                    continue
                rho_c = rho_geometric_form(p, c)
                rho_cm1 = rho_geometric_form(p, c - 1)
                self.assertEqual(
                    p ** (c - 2) * (rho_c - rho_cm1), Fraction(1, 1),
                    msg="structure statistic != 1 at p={}, c={}".format(p, c))

    def test_c3_recovers_repo_C3_over_C2(self):
        # (4 cont.) at c == 3 this is the repo's C3/C2 = p*(rho_3 - 1) = 1.
        for p in ALL_PRIMES:
            rho_3 = rho_geometric_form(p, 3)
            self.assertEqual(
                p * (rho_3 - 1), Fraction(1, 1),
                msg="p*(rho_3 - 1) != 1 at p={}".format(p))


class TestMonotoneBoundedLimit(unittest.TestCase):
    def test_strictly_increasing_bounded_with_exact_tail(self):
        # (5) rho_c strictly increasing, rho_c < p/(p-1), and the tail is exact.
        for p in ALL_PRIMES:
            limit = Fraction(p, p - 1)
            prev = None
            for c in CAPS:
                rho_c = rho_geometric_form(p, c)
                self.assertLess(
                    rho_c, limit,
                    msg="rho_c not below p/(p-1) at p={}, c={}".format(p, c))
                if prev is not None:
                    self.assertGreater(
                        rho_c, prev,
                        msg="rho_c not strictly increasing at p={}, c={}".format(p, c))
                self.assertEqual(
                    limit - rho_c, Fraction(1, p ** (c - 2) * (p - 1)),
                    msg="exact tail wrong at p={}, c={}".format(p, c))
                prev = rho_c


class TestDensityLimitsAndFloat(unittest.TestCase):
    def test_density_limit_and_float_crosscheck(self):
        # (6) |count/N - d_c| < 2/N (a-priori floor bound) + 12-place float check.
        for p in ALL_PRIMES:
            d2f = 1.0 / p ** 2 - 1.0 / p ** 3
            for c in CAPS:
                count = N_DENSITY // p ** 2 - N_DENSITY // p ** (c + 1)
                d_c = float(d_exact(p, c))
                self.assertLess(
                    abs(count / N_DENSITY - d_c),
                    2.0 / N_DENSITY + DENSITY_LIMIT_SLACK,
                    msg="density limit failed at p={}, c={}".format(p, c))
                self.assertAlmostEqual(
                    d_c / d2f, float(rho_geometric_form(p, c)),
                    places=RATIO_FLOAT_PLACES,
                    msg="float rho_c mismatch at p={}, c={}".format(p, c))


if __name__ == "__main__":
    for p in [13, 17, 19, 23]:
        row = "  ".join(
            "c={}:{}".format(c, rho_geometric_form(p, c)) for c in CAPS)
        print("p={:>2}  limit p/(p-1)={}  {}".format(
            p, Fraction(p, p - 1), row))
    unittest.main(verbosity=2)
