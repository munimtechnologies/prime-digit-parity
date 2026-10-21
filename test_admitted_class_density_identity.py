"""Proof-by-enumeration that the D5 density regressor is the admitted-class density.

Self-contained, pure-stdlib verifier. It establishes, with no network and no
sealed compute job, WHY the single-prime density scaling law of this repository
(Axis D5, v0.29.0) uses the regressor d(p) = 1/p^2 - 1/p^3 and no other form.

Mechanism (read directly from this repo's code, not assumed):
  * classes.py / axis_d5_density_scaling_law_job.py: a member n of the strict
    squarefree (mu^2) sieve satisfies v_p(n) <= caps.get(p, 1) for every prime p,
    with the default cap 1 (squarefree). "Relaxing the exponent cap at a single
    prime p" sets caps[p] = 2 and marks the forbidden power q = p^(cap+1) = p^3.
  * So relaxing at p newly admits EXACTLY the integers with v_p(n) == 2, i.e. the
    class A_p = {p^2 * m : p does not divide m} -- p^2 | n but p^3 does not.
  * The job's own density(p) returns 1.0/p**2 - 1.0/p**3 and feeds it into
    k(p) = (r(p) - r_control) / density(p).

This test proves, by exact counting over finite windows, four things:

  (1) Counting identity:  #{n <= N : v_p(n) == 2} == floor(N/p^2) - floor(N/p^3).
  (2) Density limit:      that count / N -> 1/p^2 - 1/p^3 = (p-1)/p^3, with
                          |count/N - d(p)| < 1/N (an O(1/N) approach).
  (3) Mechanism discriminator: the alternative "remove the cap entirely" (admit
                          any v_p >= 2) has density exactly 1/p^2, which exceeds
                          d(p) by exactly 1/p^3. So the regressor's FORM
                          (1/p^2 - 1/p^3, not 1/p^2) encodes that the relaxation
                          caps the exponent at 2 rather than removing it.
  (4) Repo consistency:   the density column published in
                          data/axis_d5_density_scaling_law.json equals (p-1)/p^3
                          to 1e-12 for all twelve measured primes.

Honest limit: items (1), (2), (4) are elementary identities; the point is that
the regressor is the admitted-class density fixed by the sieve definition, not a
fitted ansatz, and that its form distinguishes "cap at 2" from "no cap". The
2.7%-9% gap in (3) is below the resolution of the current single-prime ratio data
(see the companion note), so the FORM is pinned by the code, and the data is only
consistent with it. Pre-registered constants and tolerances are fixed below,
before any count is taken.

Produced by AI agents working as MV2, Munim, Inc. Apache-2.0.
"""

import unittest

# --- pre-registered constants (fixed before measurement) ---------------------
D5_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
OUT_OF_SAMPLE_PRIMES = [41, 43, 47]
ALL_PRIMES = D5_PRIMES + OUT_OF_SAMPLE_PRIMES

N_ENUM = 200_000        # window for the exact enumeration <-> formula cross-check
N_DENSITY = 5_000_000   # window for the density-limit check (formula only)

# Density column of data/axis_d5_density_scaling_law.json (exact published values).
PUBLISHED_DENSITY = {
    2: 0.125,
    3: 0.07407407407407407,
    5: 0.032,
    7: 0.017492711370262388,
    11: 0.007513148009015778,
    13: 0.005461993627674101,
    17: 0.003256665988194586,
    19: 0.0026242892549934393,
    23: 0.0018081696391879674,
    29: 0.0011480585509861003,
    31: 0.0010070155416065255,
    37: 0.0007107180226245237,
}

DENSITY_LIMIT_SLACK = 1e-15   # float slack on the exact 1/N approach bound
DISCRIMINATOR_PLACES = 12     # places for the exact 1/p^3 gap identity
PUBLISHED_ABS_TOL = 1e-12     # published density column vs closed form


def v_p(n, p):
    """p-adic valuation of n (n >= 1)."""
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def d(p):
    """Natural density of the newly admitted class A_p = {n : v_p(n) == 2}."""
    return 1.0 / (p * p) - 1.0 / (p * p * p)


class TestCountingIdentity(unittest.TestCase):
    def test_enumeration_matches_floor_formula(self):
        # (1) exact: #{n<=N : v_p(n)==2} == floor(N/p^2) - floor(N/p^3)
        for p in ALL_PRIMES:
            enum = sum(1 for n in range(1, N_ENUM + 1) if v_p(n, p) == 2)
            formula = N_ENUM // (p * p) - N_ENUM // (p * p * p)
            self.assertEqual(
                enum, formula,
                msg="counting identity failed at p={}: enum {} vs formula {}".format(
                    p, enum, formula))


class TestDensityLimit(unittest.TestCase):
    def test_count_over_N_converges_to_d(self):
        # (2) |count/N - d(p)| < 1/N (each floor truncation costs < 1/N)
        for p in ALL_PRIMES:
            count = N_DENSITY // (p * p) - N_DENSITY // (p * p * p)
            emp = count / N_DENSITY
            self.assertLess(
                abs(emp - d(p)), 1.0 / N_DENSITY + DENSITY_LIMIT_SLACK,
                msg="density limit failed at p={}: emp {} vs d {}".format(p, emp, d(p)))


class TestMechanismDiscriminator(unittest.TestCase):
    def test_unbounded_relaxation_exceeds_by_one_over_p_cubed(self):
        # (3) "admit any v_p>=2" has density 1/p^2; gap over d(p) is exactly 1/p^3.
        for p in ALL_PRIMES:
            any_power = 1.0 / (p * p)
            self.assertAlmostEqual(
                any_power - d(p), 1.0 / (p ** 3), places=DISCRIMINATOR_PLACES,
                msg="discriminator gap wrong at p={}".format(p))

    def test_classes_really_differ_in_window(self):
        # Where p^3 <= N, the "exactly 2" and "at least 2" counts must differ,
        # so the two mechanism hypotheses are genuinely distinct, not a tautology.
        for p in ALL_PRIMES:
            exactly2 = N_ENUM // (p * p) - N_ENUM // (p * p * p)
            atleast2 = N_ENUM // (p * p)
            if p ** 3 <= N_ENUM:
                self.assertGreater(
                    atleast2, exactly2,
                    msg="classes coincide at p={} inside the window".format(p))


class TestPublishedConsistency(unittest.TestCase):
    def test_published_density_equals_closed_form(self):
        # (4) the repo's regressor column equals (p-1)/p^3 == 1/p^2-1/p^3.
        for p, val in PUBLISHED_DENSITY.items():
            self.assertAlmostEqual(
                val, (p - 1) / (p ** 3), delta=PUBLISHED_ABS_TOL,
                msg="published density != (p-1)/p^3 at p={}".format(p))
            self.assertAlmostEqual(
                val, d(p), delta=PUBLISHED_ABS_TOL,
                msg="published density != 1/p^2-1/p^3 at p={}".format(p))


if __name__ == "__main__":
    # Print the realized counts for the record, then run the assertions.
    for p in ALL_PRIMES:
        count = N_DENSITY // (p * p) - N_DENSITY // (p * p * p)
        print("p={:>2}  count/N={:.12f}  d(p)={:.12f}  gap_vs_nocap=1/p^3={:.3e}".format(
            p, count / N_DENSITY, d(p), 1.0 / p ** 3))
    unittest.main(verbosity=2)
