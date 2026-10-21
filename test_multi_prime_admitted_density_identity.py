"""Proof-by-enumeration that the multi-prime admitted-class density is
sub-additive by exactly the product of the single-prime densities.

Self-contained, pure-stdlib verifier. No network, no sealed compute job, no
repository data file is read. It grounds, at the level of exact integer
counting, WHY the multi-prime density terms of this repository's density
scaling law add to first order (Axis D6 pairs v0.30.0, Axis D7 triples
v0.31.0, and the run-150 additivity audit) -- and it quantifies the exact
deviation from pure additivity.

Mechanism (established in the companion note admitted-class-density-identity.md,
read from classes.py / axis_d5_density_scaling_law_job.py): relaxing the
exponent cap at a single prime p (caps[p]=2, forbidding p^3) newly admits
exactly the class A_p = {n : v_p(n) = 2} = {p^2*m : p does not divide m}, of
natural density d(p) = 1/p^2 - 1/p^3. Relaxing at several primes admits the
UNION of their per-prime classes.

The single-prime valuations v_p and v_q at distinct primes p, q are independent
by the Chinese Remainder Theorem, so:

  * density(A_p INTERSECT A_q) = d(p)*d(q)          (independence)
  * density(A_p UNION A_q)     = d(p) + d(q) - d(p)*d(q)
  * so the additive approximation d(p)+d(q) OVERCOUNTS the true admitted
    density by EXACTLY d(p)*d(q) > 0 (strict sub-additivity). That cross term
    is O(1/(p^2 q^2)); for distinct primes both >= 11 it is < 0.33% of the
    sum, which is why Axes D6/D7 and the run-150 audit measured excess = sum
    of single excesses "to <1%" with no detectable cross term.

The general m-prime statement, by inclusion-exclusion over CRT-independent
classes, is density(UNION A_{p_i}) = 1 - prod_i (1 - d(p_i)).

This test proves, by exact counting over finite windows:

  (A) Single prime (continuity with the companion note):
      #{n<=N : v_p(n)==2} == floor(N/p^2) - floor(N/p^3), by enumeration.
  (B) Pair exact count: #{n<=N : v_p==2 OR v_q==2} equals the inclusion-
      exclusion formula, by direct enumeration (exact integer equality).
  (C) Pair decomposition: that union count == A_p + A_q - A_pq, where
      A_pq = floor(N/p^2q^2) - floor(N/p^3q^2) - floor(N/p^2q^3) + floor(N/p^3q^3).
  (D) Triple exact count: #{v_p==2 OR v_q==2 OR v_s==2} equals the full
      inclusion-exclusion formula, by direct enumeration.
  (E) Larger-window confirmation of (B)/(D) at N = 1,000,000.
  (F) Density limits: union count / N -> 1 - prod(1 - d(p_i)); intersection
      count / N -> prod d(p_i) (CRT independence), each within an a-priori
      floor bound (num_floor_terms + 1)/N.
  (G) Cross-term identity: (d(p)+d(q)) - density(union) == d(p)*d(q) to 1e-12,
      density(union) < d(p)+d(q) strictly, and the relative cross term
      d(p)d(q)/(d(p)+d(q)) < 0.01 for distinct pairs with both primes >= 11.
  (H) Three-prime grounding of D7: additive overcount Sum d_i - density(union)
      == Sum_{i<j} d_i d_j - d_p d_q d_s (dominated by the pairwise term,
      O(1/p^2 q^2)); so "additive to third order, no cross term" is the exact
      statement that the first correction is a product of two densities.

Honest limit: (A)-(F) are elementary consequences of CRT independence of
p-adic valuations; the contribution is to PIN the multi-prime density terms of
this program to that independence and to make the empirical "additive to <1%"
(D6/D7/run-150) exact -- the cross correction is -C*d(p)*d(q) in the excess,
negligible at the tested primes, not zero. Pre-registered cases and tolerances
are fixed below, before any count is taken.

Produced by AI agents working as MV2, Munim, Inc. Apache-2.0.
"""

import unittest
from math import comb
from itertools import combinations

# --- pre-registered cases and tolerances (fixed before measurement) ----------
PRIMES_IN_SAMPLE = [11, 13, 17, 19, 23, 29, 31, 37]
PRIMES_OUT = [41, 43, 47]
ALL_SINGLE_PRIMES = PRIMES_IN_SAMPLE + PRIMES_OUT

PAIR_CASES = [(11, 13), (11, 37), (13, 17), (29, 31), (17, 41), (13, 43)]
TRIPLE_CASES = [(11, 13, 17), (17, 19, 23), (11, 13, 41)]

N_ENUM = 300_000          # window for exact enumeration <-> formula cross-checks
N_BIG = 1_000_000         # larger-window confirmation
BIG_PAIR = (11, 13)
BIG_TRIPLE = (11, 13, 17)
N_DENSITY = 10_000_000    # window for the density-limit checks (formula only)

REL_CROSS_BOUND = 0.01    # relative cross term for distinct pairs, both primes >= 11
ALG_PLACES = 12           # places for the exact algebraic cross-term identity


def v_p(n, p):
    """p-adic valuation of n (n >= 1)."""
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def d(p):
    """Natural density of A_p = {n : v_p(n) == 2}."""
    return 1.0 / (p * p) - 1.0 / (p * p * p)


def count_intersection_exactly2(primes, N):
    """#{n<=N : v_p(n)==2 for EVERY p in primes}, via inclusion-exclusion over
    the cube conditions. base = prod p^2 (p^2 | n for all); subtract windows
    where some p^3 | n."""
    base = 1
    for p in primes:
        base *= p * p
    total = 0
    k = len(primes)
    for r in range(k + 1):
        for S in combinations(primes, r):
            m = base
            for p in S:
                m *= p
            total += ((-1) ** r) * (N // m)
    return total


def count_union_exactly2(primes, N):
    """#{n<=N : v_p(n)==2 for AT LEAST ONE p in primes}, by inclusion-exclusion
    over non-empty subsets of the per-prime classes."""
    total = 0
    k = len(primes)
    for r in range(1, k + 1):
        for T in combinations(primes, r):
            total += ((-1) ** (r + 1)) * count_intersection_exactly2(list(T), N)
    return total


def enum_union_exactly2(primes, N):
    """Brute-force count of n<=N having v_p(n)==2 for at least one listed prime."""
    c = 0
    for n in range(1, N + 1):
        for p in primes:
            if v_p(n, p) == 2:
                c += 1
                break
    return c


def density_union(primes):
    """Limit density of the union = 1 - prod(1 - d(p))."""
    prod = 1.0
    for p in primes:
        prod *= (1.0 - d(p))
    return 1.0 - prod


def num_floor_terms(primes):
    """Count of floor() terms in the full subset inclusion-exclusion for the
    union, used to bound the finite-N truncation error a priori."""
    k = len(primes)
    return sum(comb(k, r) * (2 ** r) for r in range(1, k + 1))


class TestSinglePrime(unittest.TestCase):
    def test_single_prime_counting_identity(self):
        # (A) continuity with the companion note, by enumeration.
        for p in ALL_SINGLE_PRIMES:
            enum = sum(1 for n in range(1, N_ENUM + 1) if v_p(n, p) == 2)
            formula = count_intersection_exactly2([p], N_ENUM)
            self.assertEqual(formula, N_ENUM // (p * p) - N_ENUM // (p * p * p))
            self.assertEqual(
                enum, formula,
                msg="single-prime identity failed at p={}: enum {} vs formula {}".format(
                    p, enum, formula))


class TestPairExactCount(unittest.TestCase):
    def test_pair_union_enumeration_matches_formula(self):
        # (B) exact integer equality by direct enumeration.
        for (p, q) in PAIR_CASES:
            enum = enum_union_exactly2([p, q], N_ENUM)
            formula = count_union_exactly2([p, q], N_ENUM)
            self.assertEqual(
                enum, formula,
                msg="pair union count mismatch at ({},{}): enum {} vs formula {}".format(
                    p, q, enum, formula))

    def test_pair_inclusion_exclusion_decomposition(self):
        # (C) union == A_p + A_q - A_pq, exactly.
        for (p, q) in PAIR_CASES:
            A_p = count_intersection_exactly2([p], N_ENUM)
            A_q = count_intersection_exactly2([q], N_ENUM)
            A_pq = count_intersection_exactly2([p, q], N_ENUM)
            self.assertEqual(
                count_union_exactly2([p, q], N_ENUM), A_p + A_q - A_pq,
                msg="pair decomposition failed at ({},{})".format(p, q))
            # A_pq is the explicit 4-term floor expression.
            explicit = (N_ENUM // (p * p * q * q)
                        - N_ENUM // (p * p * p * q * q)
                        - N_ENUM // (p * p * q * q * q)
                        + N_ENUM // (p * p * p * q * q * q))
            self.assertEqual(A_pq, explicit,
                             msg="A_pq floor expression wrong at ({},{})".format(p, q))


class TestTripleExactCount(unittest.TestCase):
    def test_triple_union_enumeration_matches_formula(self):
        # (D) exact integer equality by direct enumeration.
        for (p, q, s) in TRIPLE_CASES:
            enum = enum_union_exactly2([p, q, s], N_ENUM)
            formula = count_union_exactly2([p, q, s], N_ENUM)
            self.assertEqual(
                enum, formula,
                msg="triple union count mismatch at ({},{},{}): enum {} vs formula {}".format(
                    p, q, s, enum, formula))


class TestLargerWindow(unittest.TestCase):
    def test_pair_and_triple_at_one_million(self):
        # (E) confirm (B)/(D) survive a larger window.
        self.assertEqual(
            enum_union_exactly2(list(BIG_PAIR), N_BIG),
            count_union_exactly2(list(BIG_PAIR), N_BIG),
            msg="big-pair union mismatch at N={}".format(N_BIG))
        self.assertEqual(
            enum_union_exactly2(list(BIG_TRIPLE), N_BIG),
            count_union_exactly2(list(BIG_TRIPLE), N_BIG),
            msg="big-triple union mismatch at N={}".format(N_BIG))


class TestDensityLimits(unittest.TestCase):
    def test_intersection_density_is_product(self):
        # (F) A_pq / N -> d(p)*d(q): CRT independence of the valuations.
        for (p, q) in PAIR_CASES:
            count = count_intersection_exactly2([p, q], N_DENSITY)
            emp = count / N_DENSITY
            bound = (num_floor_terms([p, q]) + 1) / N_DENSITY
            self.assertLess(
                abs(emp - d(p) * d(q)), bound,
                msg="intersection density != d(p)d(q) at ({},{})".format(p, q))

    def test_union_density_is_one_minus_product(self):
        # (F) union count / N -> 1 - prod(1 - d(p_i)).
        for case in PAIR_CASES + TRIPLE_CASES:
            primes = list(case)
            count = count_union_exactly2(primes, N_DENSITY)
            emp = count / N_DENSITY
            bound = (num_floor_terms(primes) + 1) / N_DENSITY
            self.assertLess(
                abs(emp - density_union(primes)), bound,
                msg="union density != 1-prod(1-d) at {}".format(case))


class TestCrossTerm(unittest.TestCase):
    def test_pair_cross_term_identity_and_sign(self):
        # (G) (d_p + d_q) - density(union) == d_p*d_q, strictly positive.
        max_rel = 0.0
        for (p, q) in PAIR_CASES:
            du = density_union([p, q])
            self.assertAlmostEqual(
                (d(p) + d(q)) - du, d(p) * d(q), places=ALG_PLACES,
                msg="cross-term identity wrong at ({},{})".format(p, q))
            self.assertLess(du, d(p) + d(q),
                            msg="union not sub-additive at ({},{})".format(p, q))
            if p >= 11 and q >= 11:
                rel = d(p) * d(q) / (d(p) + d(q))
                max_rel = max(max_rel, rel)
        self.assertLess(
            max_rel, REL_CROSS_BOUND,
            msg="relative cross term {} exceeds pre-registered bound".format(max_rel))


class TestTripleGrounding(unittest.TestCase):
    def test_additive_overcount_is_pairwise_minus_triple(self):
        # (H) Sum d_i - density(union) == Sum_{i<j} d_i d_j - d_p d_q d_s.
        for (p, q, s) in TRIPLE_CASES:
            ds = [d(p), d(q), d(s)]
            sum_singles = sum(ds)
            overcount = sum_singles - density_union([p, q, s])
            pairwise = d(p) * d(q) + d(p) * d(s) + d(q) * d(s)
            triple = d(p) * d(q) * d(s)
            self.assertAlmostEqual(
                overcount, pairwise - triple, places=ALG_PLACES,
                msg="triple overcount identity wrong at ({},{},{})".format(p, q, s))
            self.assertGreater(overcount, 0.0)


if __name__ == "__main__":
    for case in PAIR_CASES + TRIPLE_CASES:
        primes = list(case)
        emp = count_union_exactly2(primes, N_DENSITY) / N_DENSITY
        lim = density_union(primes)
        add = sum(d(p) for p in primes)
        print("{}: union/N={:.10f}  limit={:.10f}  additive={:.10f}  "
              "overcount={:.3e}".format(case, emp, lim, add, add - lim))
    unittest.main(verbosity=2)
