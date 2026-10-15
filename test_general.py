# test_general.py -- self-contained, stdlib-only, exact tests for GENERAL.md.
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.
#
# Everything here is redefined inline: no repo imports, no file reads, stdlib
# only.  The heavy pointwise checks sieve to 1e6 (78498 primes) once and read
# every number back through independent digit-sum definitions.  All arithmetic
# that matters is exact (Fraction / integer); nothing is floating point.

import unittest
from collections import Counter
from fractions import Fraction
from math import comb, factorial


# ----------------------------------------------------------------------------
# primes and digit-sum parity characters (independent definitions)
# ----------------------------------------------------------------------------
def primes_below(n):
    sieve = bytearray([1]) * n
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    return [i for i in range(n) if sieve[i]]


def digit_sum_parity(m, base):
    """Parity of the sum of the base-`base` digits of m (independent of bits)."""
    s = 0
    while m:
        s += m % base
        m //= base
    return s & 1


def chi(m, base):
    """(-1)^{digit-sum of m in `base`}: +1 even digit-sum, -1 odd."""
    return 1 - 2 * digit_sum_parity(m, base)


# ----------------------------------------------------------------------------
# set partitions (recursive, exact enumeration)
# ----------------------------------------------------------------------------
def set_partitions(elements):
    elements = list(elements)
    if not elements:
        yield []
        return
    first = elements[0]
    for rest in set_partitions(elements[1:]):
        for i in range(len(rest)):
            yield rest[:i] + [[first] + rest[i]] + rest[i + 1:]
        yield [[first]] + rest


def bell_numbers(upto):
    """B_0 .. B_upto via the Bell triangle."""
    row = [1]
    bells = [1]
    for _ in range(upto):
        nxt = [row[-1]]
        for x in row:
            nxt.append(nxt[-1] + x)
        row = nxt
        bells.append(row[0])
    return bells


# ----------------------------------------------------------------------------
# 1. The general character factorisation (tower reduction) holds for all b
# ----------------------------------------------------------------------------
class TestTowerFactorisation(unittest.TestCase):
    def test_chi2_is_product_of_b_shifted_base2b_parities(self):
        # chi_2(n) = prod_{r=0}^{b-1} chi_{2^b}(n >> r), checked pointwise with
        # the INDEPENDENT base-2^b digit-sum definition (not via bit masks).
        for b in (2, 3, 4, 5):
            base = 1 << b
            for n in range(2, 20000):
                prod = 1
                for r in range(b):
                    prod *= chi(n >> r, base)
                self.assertEqual(prod, chi(n, 2))


# ----------------------------------------------------------------------------
# 2. Anchor gate: the published excesses E_b = #odd - #even = -sum chi_b @1e6
# ----------------------------------------------------------------------------
class TestAnchors(unittest.TestCase):
    def test_excess_anchors_at_1e6(self):
        ps = primes_below(1_000_000)
        self.assertEqual(len(ps), 78498)
        # powers-of-two bases (the tower bases) at 1e6, from the repo's data
        for base, anchor in ((2, 4764), (4, 2902), (8, -510), (16, 244)):
            s = sum(chi(p, base) for p in ps)
            self.assertEqual(-s, anchor)   # sign convention: E = -sum(chi)


# ----------------------------------------------------------------------------
# 3. General-b singleton-free partition counts = OEIS A000296, three ways
# ----------------------------------------------------------------------------
class TestSingletonFreeCounts(unittest.TestCase):
    # A000296(n) for n = 2..8 (partitions of an n-set into blocks of size >= 2).
    EXPECTED = {2: 1, 3: 1, 4: 4, 5: 11, 6: 41, 7: 162, 8: 715}

    def test_three_independent_derivations_agree(self):
        bells = bell_numbers(8)
        self.assertEqual(bells, [1, 1, 2, 5, 15, 52, 203, 877, 4140])
        for n, exp in self.EXPECTED.items():
            # (a) brute enumeration of all set partitions, keep singleton-free
            cnt = sum(
                1 for pi in set_partitions(list(range(n)))
                if all(len(block) >= 2 for block in pi)
            )
            # (b) Bell-number inclusion-exclusion:
            #     a(n) = sum_{k=0}^{n} (-1)^k C(n,k) B_{n-k}
            incl = sum((-1) ** k * comb(n, k) * bells[n - k] for k in range(n + 1))
            self.assertEqual(cnt, exp)
            self.assertEqual(incl, exp)

    def test_purity_boundary_is_a_corollary(self):
        # "pure top cumulant" <=> exactly ONE singleton-free partition (the
        # whole set) <=> a(n) == 1 <=> n in {2, 3}.  First break at n = 4.
        a = {
            n: sum(
                1 for pi in set_partitions(list(range(n)))
                if all(len(block) >= 2 for block in pi)
            )
            for n in range(2, 9)
        }
        self.assertEqual([n for n in a if a[n] == 1], [2, 3])
        self.assertTrue(all(a[n] > 1 for n in range(4, 9)))
        self.assertEqual(a[4], 4)   # kappa_4 + three kappa_2.kappa_2 pairings


# ----------------------------------------------------------------------------
# 4. Explicit cumulant-structure tables for b = 5 and b = 6 (by block sizes)
# ----------------------------------------------------------------------------
class TestBlockTables(unittest.TestCase):
    @staticmethod
    def size_table(n):
        table = Counter()
        for pi in set_partitions(list(range(n))):
            if all(len(block) >= 2 for block in pi):
                table[tuple(sorted(len(block) for block in pi))] += 1
        return dict(table)

    def test_b5(self):
        # kappa_5 (1) + ten kappa_3.kappa_2 (3+2 splits) = 11
        self.assertEqual(self.size_table(5), {(5,): 1, (2, 3): 10})

    def test_b6(self):
        # kappa_6 + 15 kappa_4.kappa_2 + 10 kappa_3.kappa_3 + 15 kappa_2^3 = 41
        self.assertEqual(
            self.size_table(6),
            {(6,): 1, (2, 4): 15, (3, 3): 10, (2, 2, 2): 15},
        )
        self.assertEqual(sum(self.size_table(6).values()), 41)


# ----------------------------------------------------------------------------
# 5. Full b = 5 moment->cumulant round trip on real prime data (exact)
# ----------------------------------------------------------------------------
class TestB5RoundTrip(unittest.TestCase):
    def test_roundtrip_reproduces_the_anchor(self):
        B = 5
        base = 1 << B  # 32
        ps = primes_below(1_000_000)
        N = len(ps)
        self.assertEqual(N, 78498)

        # Histogram over the 5-bit vector of shifted base-32 digit-sum parities.
        # Also check pointwise that the product of the five columns is chi_2.
        f = [0] * (1 << B)
        structural_ok = True
        for p in ps:
            v = 0
            for r in range(B):
                v |= digit_sum_parity(p >> r, base) << r
            f[v] += 1
            if (bin(v).count("1") & 1) != (bin(p).count("1") & 1):
                structural_ok = False
        self.assertTrue(structural_ok)

        def moment_num(mask):
            tot = 0
            for v in range(1 << B):
                fv = f[v]
                if fv:
                    tot += -fv if (bin(v & mask).count("1") & 1) else fv
            return tot

        full = (1 << B) - 1
        # Gate: the full product reproduces sum chi_2 = -E_2(1e6) = -4764.
        self.assertEqual(moment_num(full), -4764)

        # Joint moments for every subset, exact.
        mom = {}
        for mask in range(1 << B):
            key = frozenset(i for i in range(B) if (mask >> i) & 1)
            mom[key] = Fraction(moment_num(mask), N)

        # Joint cumulants via Mobius inversion over partitions.
        def kappa(block):
            tot = Fraction(0)
            for pi in set_partitions(list(block)):
                coef = (-1) ** (len(pi) - 1) * factorial(len(pi) - 1)
                prod = Fraction(1)
                for blk in pi:
                    prod *= mom[frozenset(blk)]
                tot += coef * prod
            return tot

        kap = {}
        for mask in range(1, 1 << B):
            key = frozenset(i for i in range(B) if (mask >> i) & 1)
            kap[key] = kappa(key)

        # Reconstruct E[prod of all five] from cumulants: the Bell(5)=52
        # partition sum.  This is an exact algebraic round trip.
        recon = Fraction(0)
        npart = 0
        for pi in set_partitions(list(range(B))):
            npart += 1
            prod = Fraction(1)
            for blk in pi:
                prod *= kap[frozenset(blk)]
            recon += prod
        self.assertEqual(npart, 52)
        self.assertEqual(recon, mom[frozenset(range(B))])
        self.assertEqual(recon, Fraction(-4764, N))

        # Of the 52 partitions exactly 11 are singleton-free: one kappa_5 and
        # ten kappa_3.kappa_2.  The other 41 all carry >=1 single-parity
        # marginal mean and are therefore leading-order suppressed.
        sf = [
            pi for pi in set_partitions(list(range(B)))
            if all(len(blk) >= 2 for blk in pi)
        ]
        self.assertEqual(len(sf), 11)
        self.assertEqual(sum(1 for pi in sf if len(pi) == 1), 1)
        self.assertEqual(sum(1 for pi in sf if len(pi) == 2), 10)
        self.assertTrue(
            all(sorted(len(b) for b in pi) == [2, 3] for pi in sf if len(pi) == 2)
        )

        # The singleton marginals here are genuinely small, which is WHY the 41
        # singleton-carrying partitions are suppressed.  A loose, certainly-true
        # bound (every single-character excess at 1e6 is in the hundreds, far
        # below 0.1 * N) keeps the test robust; the measured values are ~1e-3.
        for i in range(B):
            self.assertLess(abs(mom[frozenset([i])]), Fraction(1, 10))


if __name__ == "__main__":
    unittest.main(verbosity=2)
