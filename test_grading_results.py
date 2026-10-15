# test_grading_results.py -- Self-contained checks for GRADING-RESULTS.md.
# Stdlib only; redefines everything inline (no repo imports, no file reads).
# It (1) sieves the primes below 1e6, rebuilds the three base-8 parity columns
# from the bit positions, and verifies they reproduce the published anchors
# E_2(1e6)=4764 and E_8(1e6)=-510 and the exact moment-cumulant reconstruction;
# then (2) gates the hardcoded 1e8 octal-triple counts on their anchors and
# checks the singleton-graded (k=0/1/2/3) split of E_2(1e8)=70509.
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.

import unittest
from fractions import Fraction


def primes_upto(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(2, limit + 1) if sieve[i]]


def col_parity(n, r):
    """Parity (0 even, 1 odd) of the bits of n at positions r, r+3, r+6, ...

    This is the parity of the r-th base-8 digit-sum column: s_8(n>>r) mod 2
    equals the XOR of the lowest bit of every octal digit of n>>r, i.e. the
    bits of n at positions congruent to r modulo 3.
    """
    nn = n >> r
    p = 0
    while nn:
        p ^= (nn & 1)
        nn >>= 3
    return p


def counts_from_primes(primes):
    """Eight octal-triple counts, index 4*a+2*b+c with a=col0,b=col1,c=col2."""
    counts = [0] * 8
    for p in primes:
        i = (col_parity(p, 0) << 2) | (col_parity(p, 1) << 1) | col_parity(p, 2)
        counts[i] += 1
    return counts


def sign(i, which):
    return 1 if ((i >> (2 - which)) & 1) == 0 else -1


def moment(counts, n, subset):
    tot = Fraction(0)
    for i, c in enumerate(counts):
        s = 1
        for w in subset:
            s *= sign(i, w)
        tot += Fraction(s * c, n)
    return tot


def cumulant(counts, n, subset):
    s = tuple(subset)
    m = lambda ss: moment(counts, n, ss)
    if len(s) == 1:
        return m(s)
    if len(s) == 2:
        x, y = s
        return m((x, y)) - m((x,)) * m((y,))
    x, y, z = s
    return (m((x, y, z)) - m((x,)) * m((y, z)) - m((y,)) * m((x, z))
            - m((z,)) * m((x, y)) + 2 * m((x,)) * m((y,)) * m((z,)))


def set_partitions(items):
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for smaller in set_partitions(rest):
        for k in range(len(smaller)):
            yield smaller[:k] + [(first,) + smaller[k]] + smaller[k + 1:]
        yield [(first,)] + smaller


def graded(counts, n):
    buckets = {0: Fraction(0), 1: Fraction(0), 2: Fraction(0), 3: Fraction(0)}
    for part in set_partitions((0, 1, 2)):
        k = sum(1 for b in part if len(b) == 1)
        term = Fraction(1)
        for b in part:
            term *= cumulant(counts, n, b)
        buckets[k] += term
    return buckets


COUNTS_1E8 = [693845, 736247, 734665, 716635, 736921, 719439, 715554, 708149]


class TestGradingResults(unittest.TestCase):

    def test_independent_sieve_reproduces_1e6_anchors(self):
        primes = primes_upto(1_000_000)
        self.assertEqual(len(primes), 78498)
        counts = counts_from_primes(primes)
        self.assertEqual(sum(counts), 78498)
        n = sum(counts)
        sum_abc = int(moment(counts, n, (0, 1, 2)) * n)
        sum_a = int(moment(counts, n, (0,)) * n)
        # E_2 = -sum(ABC) (tower reduction chi_2 = A*B*C); E_8 = -sum(A).
        self.assertEqual(-sum_abc, 4764)   # E_2(1e6) anchor
        self.assertEqual(-sum_a, -510)     # E_8(1e6) anchor
        # moment-cumulant reconstruction is exact at 1e6 too.
        buckets = graded(counts, n)
        self.assertEqual(sum(buckets.values()), moment(counts, n, (0, 1, 2)))
        self.assertEqual(buckets[2], Fraction(0))

    def test_1e8_counts_gate_on_anchors(self):
        self.assertEqual(sum(COUNTS_1E8), 5761455)  # pi(1e8)
        n = sum(COUNTS_1E8)
        sum_a = int(moment(COUNTS_1E8, n, (0,)) * n)
        sum_abc = int(moment(COUNTS_1E8, n, (0, 1, 2)) * n)
        self.assertEqual(sum_a, 1329)        # -> E_8(1e8) = -1329
        self.assertEqual(sum_abc, -70509)    # -> E_2(1e8) = 70509

    def test_singleton_graded_split_is_exact_and_matches_prose(self):
        n = sum(COUNTS_1E8)
        buckets = graded(COUNTS_1E8, n)
        # Orders sum exactly to the whole base-2 mean E[ABC] = sum(ABC)/N.
        self.assertEqual(sum(buckets.values()),
                         moment(COUNTS_1E8, n, (0, 1, 2)))
        self.assertEqual(sum(buckets.values()), Fraction(-70509, n))
        # P(3,2) = 0: there is no order-2 (two-singleton) term, exactly.
        self.assertEqual(buckets[2], Fraction(0))
        # Contributions to E_2 = -N * (mean contribution).
        e2 = {k: -buckets[k] * n for k in buckets}
        self.assertEqual(round(float(e2[0])), 70395)  # k=0 singleton-free 3pt
        self.assertEqual(round(float(e2[1])), 114)     # k=1 three marginal terms
        self.assertEqual(round(float(e2[3])), 0)       # k=3 fully disconnected
        self.assertEqual(sum(e2.values()), Fraction(70509))
        # k=0 carries 99.8% of the excess; the whole correction is < 0.2%.
        share0 = float(e2[0]) / 70509
        self.assertGreater(share0, 0.998)
        self.assertLess(share0, 0.999)

    def test_k1_breakdown_three_terms(self):
        n = sum(COUNTS_1E8)
        m = lambda ss: moment(COUNTS_1E8, n, ss)
        # The three one-singleton terms: <single> * kappa(pair).
        termA = m((0,)) * (m((1, 2)) - m((1,)) * m((2,)))  # A singleton, BC pair
        termB = m((1,)) * (m((0, 2)) - m((0,)) * m((2,)))  # B singleton, AC pair
        termC = m((2,)) * (m((0, 1)) - m((0,)) * m((1,)))  # C singleton, AB pair
        buckets = graded(COUNTS_1E8, n)
        self.assertEqual(termA + termB + termC, buckets[1])
        # Measured E_2 shares of each (B is the largest marginal).
        self.assertEqual(round(float(-termA * n)), 12)
        self.assertEqual(round(float(-termB * n)), 98)
        self.assertEqual(round(float(-termC * n)), 5)


if __name__ == "__main__":
    unittest.main(verbosity=2)
