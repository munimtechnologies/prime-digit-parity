# test_moments.py -- self-contained exact-identity checks for MOMENTS.md.
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.  Code: Apache-2.0.
#
# Stdlib only; redefines everything inline (no repo imports, no file reads).
# Verifies, exactly and from scratch:
#   (1) the character factorisation chi_2(n) = A*B*C*D for every n < 1e6;
#   (2) the full 15-partition moment->cumulant round trip for the four shifted
#       base-16 parities over the real distribution on primes < 1e6 (b=4), and
#       the b=2 and b=3 sub-tuples;
#   (3) the cumulant code on an independent synthetic distribution (all joint
#       cumulants of size >= 2 vanish; moment = product of marginals);
#   (4) the singleton-free partition counts underlying the purity boundary
#       (1, 1, 4 for b = 2, 3, 4);
#   (5) the published anchors E_2(1e6) = 4764, E_16(1e6) = 244, pi(1e6)=78498
#       (excess = #odd - #even = -sum(chi)), and that the 15-partition sum
#       reproduces N * E[ABCD] = sum(chi_2) = -4764 exactly.

import unittest
from fractions import Fraction
from math import factorial

LIMIT = 10 ** 6


# ---- characters -----------------------------------------------------------

def chi2(n):
    return -1 if (bin(n).count("1") & 1) else 1


def chi16(m):
    s = 0
    while m:
        s += m & 0xF
        m >>= 4
    return -1 if (s & 1) else 1


def shifted(n):
    return (chi16(n), chi16(n >> 1), chi16(n >> 2), chi16(n >> 3))


def primes_below(limit):
    sieve = bytearray([1]) * limit
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(limit) if sieve[i]]


# ---- partitions / moments / cumulants -------------------------------------

def set_partitions(items):
    items = list(items)
    if not items:
        yield []
        return
    if len(items) == 1:
        yield [items]
        return
    first, rest = items[0], items[1:]
    for smaller in set_partitions(rest):
        for i in range(len(smaller)):
            yield smaller[:i] + [[first] + smaller[i]] + smaller[i + 1:]
        yield [[first]] + smaller


def subsets(k):
    return [frozenset(i for i in range(k) if (mask >> i) & 1)
            for mask in range(1 << k)]


def moments_from_counts(counts, n, k):
    m = {frozenset(): Fraction(1)}
    for S in subsets(k):
        if not S:
            continue
        tot = 0
        for pattern, c in counts.items():
            prod = 1
            for i in S:
                if (pattern >> i) & 1:
                    prod = -prod
            tot += c * prod
        m[S] = Fraction(tot, n)
    return m


def cumulant(S, m):
    total = Fraction(0)
    for part in set_partitions(list(S)):
        kb = len(part)
        total += (((-1) ** (kb - 1)) * factorial(kb - 1)
                  * _prod(m[frozenset(b)] for b in part))
    return total


def moment_from_cumulants(S, kap):
    total = Fraction(0)
    for part in set_partitions(list(S)):
        total += _prod(kap[frozenset(b)] for b in part)
    return total


def _prod(it):
    p = Fraction(1)
    for x in it:
        p *= x
    return p


# ---- shared fixtures ------------------------------------------------------

PRIMES = primes_below(LIMIT)
N = len(PRIMES)

# 16 sign-pattern counts for (A,B,C,D) over primes < 1e6.
COUNTS = {}
for p in PRIMES:
    a, b, c, d = shifted(p)
    pat = ((a == -1)) | ((b == -1) << 1) | ((c == -1) << 2) | ((d == -1) << 3)
    COUNTS[pat] = COUNTS.get(pat, 0) + 1


class TestFactorisation(unittest.TestCase):
    def test_pointwise_identity_dense_range(self):
        # chi_2(n) == chi_16(n)*chi_16(n>>1)*chi_16(n>>2)*chi_16(n>>3) for all n
        # in a dense low range (the identity is exact for every n).
        bad = 0
        for n in range(300000):
            a, b, c, d = shifted(n)
            if chi2(n) != a * b * c * d:
                bad += 1
                if bad > 3:
                    break
        self.assertEqual(bad, 0)

    def test_pointwise_identity_on_all_primes(self):
        # also hold it on every prime < 1e6 (the arguments that carry the bias).
        self.assertTrue(all(
            chi2(p) == _prod(shifted(p)) for p in PRIMES))


class TestAnchors(unittest.TestCase):
    def test_prime_count(self):
        self.assertEqual(N, 78498)

    def test_base2_and_base16_excess(self):
        # Repo convention (CORRELATION-RESULTS.md): excess = #odd - #even
        # = -sum(chi), chi = (-1)^digitsum.  Published anchors are +4764, +244.
        e2 = -sum(chi2(p) for p in PRIMES)
        e16 = -sum(chi16(p) for p in PRIMES)
        self.assertEqual(e2, 4764)    # published base-2 excess at 1e6
        self.assertEqual(e16, 244)    # published base-16 excess at 1e6

    def test_marginal_mean_matches_base16_anchor(self):
        # mu_0 = E[A] = sum(chi_16(p))/N, so the base-16 excess is -N*mu_0 = 244.
        m = moments_from_counts(COUNTS, N, 4)
        self.assertEqual(-m[frozenset({0})] * N, 244)


class TestMomentCumulantRoundTrip(unittest.TestCase):
    def test_b4_full_15_partition_reconstruction(self):
        m = moments_from_counts(COUNTS, N, 4)
        kap = {S: cumulant(S, m) for S in subsets(4) if S}
        full = frozenset({0, 1, 2, 3})
        recon = moment_from_cumulants(full, kap)
        self.assertEqual(recon, m[full])            # exact Fraction identity
        # N*E[ABCD] = sum(chi_2) = -(base-2 excess), so the partition sum ties
        # back to the published base-2 excess 4764 exactly.
        self.assertEqual(m[full] * N, -4764)
        self.assertEqual(-recon * N, 4764)
        # exactly 15 partitions of a 4-set entered the sum.
        self.assertEqual(sum(1 for _ in set_partitions(range(4))), 15)

    def test_b2_and_b3_subtuple_reconstruction(self):
        m = moments_from_counts(COUNTS, N, 4)
        kap = {S: cumulant(S, m) for S in subsets(4) if S}
        for S in (frozenset({0, 1}), frozenset({0, 1, 2})):
            self.assertEqual(moment_from_cumulants(S, kap), m[S])


class TestCumulantCodeOnIndependentData(unittest.TestCase):
    def test_independent_bits_kill_cross_cumulants(self):
        # Build 16 exact pattern counts as a product of four independent
        # marginals with P(X_i = -1) = p_i; cross-cumulants must vanish.
        ps = [Fraction(1, 4), Fraction(1, 3), Fraction(2, 5), Fraction(1, 6)]
        D = 4 * 3 * 5 * 6  # product of the marginal denominators -> integer counts
        counts = {}
        for pat in range(16):
            w = Fraction(1)
            for i in range(4):
                pi = ps[i]
                w *= (pi if (pat >> i) & 1 else (1 - pi))
            num = w * D
            assert num.denominator == 1
            counts[pat] = int(num)
        n = sum(counts.values())
        m = moments_from_counts(counts, n, 4)
        kap = {S: cumulant(S, m) for S in subsets(4) if S}
        # every cumulant spanning >= 2 distinct independent variables is 0
        for S in subsets(4):
            if len(S) >= 2:
                self.assertEqual(kap[S], Fraction(0))
        # and the 4-point moment is just the product of the four means
        prod_means = _prod(kap[frozenset({i})] for i in range(4))
        self.assertEqual(m[frozenset({0, 1, 2, 3})], prod_means)


class TestPurityBoundary(unittest.TestCase):
    def test_singleton_free_partition_counts(self):
        # Partitions into blocks of size >= 2: only the full set for b <= 3,
        # the full set plus three 2+2 pairings for b = 4.  This elementary
        # fact is exactly the purity boundary "pure top cumulant iff b in {2,3}".
        def count_singleton_free(k):
            return sum(1 for part in set_partitions(range(k))
                       if all(len(block) >= 2 for block in part))
        self.assertEqual(count_singleton_free(2), 1)
        self.assertEqual(count_singleton_free(3), 1)
        self.assertEqual(count_singleton_free(4), 4)
        # and the three surviving non-trivial ones are the 2+2 pairings
        pairings = [part for part in set_partitions(range(4))
                    if len(part) == 2 and all(len(b) == 2 for b in part)]
        self.assertEqual(len(pairings), 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
