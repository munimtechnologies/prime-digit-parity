"""
Self-contained exact-algebra tests for the cumulant tower (HIERARCHY.md).

No imports from the repo, no file reads, stdlib only. Every assertion is an
exact identity in integers or Fractions, so each test passes by mathematics,
not by any sieve output. Work by Munim, Inc., produced with AI. Apache-2.0.
"""
import unittest
from fractions import Fraction


def digit_sum(n, base):
    s = 0
    while n:
        s += n % base
        n //= base
    return s


def chi(n, q):
    """chi_q(n) = (-1)^{s_q(n)} in {+1,-1}."""
    return 1 - 2 * (digit_sum(n, q) & 1)


def mean(col):
    return Fraction(sum(col), len(col))


def jmom(*cols):
    n = len(cols[0])
    tot = 0
    for i in range(n):
        p = 1
        for c in cols:
            p *= c[i]
        tot += p
    return Fraction(tot, n)


def set_partitions(items):
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def _fact(k):
    r = 1
    for i in range(2, k + 1):
        r *= i
    return r


def kappa(cols):
    """Joint cumulant by Leonov-Shiryaev, exact Fraction."""
    k = len(cols)
    tot = Fraction(0)
    for part in set_partitions(range(k)):
        m = len(part)
        term = Fraction((-1) ** (m - 1)) * _fact(m - 1)
        for block in part:
            term *= jmom(*[cols[j] for j in block])
        tot += term
    return tot


class TestReductionIdentity(unittest.TestCase):
    def test_general_identity_pointwise(self):
        # chi_{2^a}(n) == prod_r chi_{2^{ab}}(n >> r*a), identity (star).
        for a in (1, 2, 3):
            for b in (1, 2, 3, 4):
                Q = 1 << (a * b)
                q = 1 << a
                for n in range(1 << 13):
                    prod = 1
                    for r in range(b):
                        prod *= chi(n >> (r * a), Q)
                    self.assertEqual(prod, chi(n, q))

    def test_base2_tower(self):
        # The SAME chi_2(n) factors three ways: base-4 (b=2), base-8 (b=3),
        # base-16 (b=4). All three equal chi_2(n) pointwise.
        for n in range(1 << 16):
            t2 = chi(n, 4) * chi(n >> 1, 4)
            t3 = chi(n, 8) * chi(n >> 1, 8) * chi(n >> 2, 8)
            t4 = (chi(n, 16) * chi(n >> 1, 16)
                  * chi(n >> 2, 16) * chi(n >> 3, 16))
            self.assertEqual(t2, chi(n, 2))
            self.assertEqual(t3, chi(n, 2))
            self.assertEqual(t4, chi(n, 2))

    def test_base4_from_base16(self):
        # a=2, b=2: chi_4(n) == chi_16(n) * chi_16(n >> 2).
        for n in range(1 << 16):
            self.assertEqual(chi(n, 16) * chi(n >> 2, 16), chi(n, 4))


class TestBlockVanishing(unittest.TestCase):
    def test_complete_block_sums_are_zero(self):
        # Over a complete base-Q block the digit-sum parity is balanced, so the
        # (unnormalised) joint correlation is EXACTLY zero over the integers.
        n4 = 16 * 16                      # [0, 16^2) = complete base-16 block
        self.assertEqual(sum(chi(n, 4) for n in range(n4)), 0)
        self.assertEqual(
            sum(chi(n, 16) * chi(n >> 2, 16) for n in range(n4)), 0)
        n2 = 8 * 8                        # [0, 8^2) = complete base-8 block
        self.assertEqual(sum(chi(n, 2) for n in range(n2)), 0)
        self.assertEqual(
            sum(chi(n, 8) * chi(n >> 1, 8) * chi(n >> 2, 8)
                for n in range(n2)), 0)


class TestMomentCumulant(unittest.TestCase):
    """Exact moment<->cumulant algebra on deterministic +/-1 columns."""

    def setUp(self):
        # Real, deterministic columns from shifted base-16 parities on [0,256).
        self.cols4 = [
            [chi(n, 16) for n in range(256)],
            [chi(n >> 1, 16) for n in range(256)],
            [chi(n >> 2, 16) for n in range(256)],
            [chi(n >> 3, 16) for n in range(256)],
        ]

    def test_moment_rebuilds_from_cumulants(self):
        # E[prod] == sum over partitions of prod of block cumulants, for k=2,3,4.
        for k in (2, 3, 4):
            cols = self.cols4[:k]
            rebuilt = Fraction(0)
            for part in set_partitions(range(k)):
                term = Fraction(1)
                for block in part:
                    term *= kappa([cols[j] for j in block])
                rebuilt += term
            self.assertEqual(jmom(*cols), rebuilt)

    def test_pairwise_decomposition(self):
        # b=2: E[XY] = kappa_2 + E[X]E[Y], exactly.
        X, Y = self.cols4[0], self.cols4[1]
        k2 = jmom(X, Y) - mean(X) * mean(Y)
        self.assertEqual(jmom(X, Y), k2 + mean(X) * mean(Y))
        self.assertEqual(k2, kappa([X, Y]))

    def test_triple_decomposition(self):
        # b=3: E[XYZ] = kappa_3 + sum(marginal * pair-cov) + prod(marginals).
        X, Y, Z = self.cols4[0], self.cols4[1], self.cols4[2]
        mX, mY, mZ = mean(X), mean(Y), mean(Z)
        k2XY = jmom(X, Y) - mX * mY
        k2XZ = jmom(X, Z) - mX * mZ
        k2YZ = jmom(Y, Z) - mY * mZ
        k3 = kappa([X, Y, Z])
        rhs = k3 + mX * k2YZ + mY * k2XZ + mZ * k2XY + mX * mY * mZ
        self.assertEqual(jmom(X, Y, Z), rhs)

    def test_purity_boundary_b2_b3(self):
        # No-singleton partitions: for b=2 it is exactly kappa_2, for b=3 it is
        # exactly kappa_3 -- the bias is a PURE top cumulant at leading order.
        for k in (2, 3):
            cols = self.cols4[:k]
            no_single = Fraction(0)
            for part in set_partitions(range(k)):
                if any(len(block) == 1 for block in part):
                    continue
                term = Fraction(1)
                for block in part:
                    term *= kappa([cols[j] for j in block])
                no_single += term
            self.assertEqual(no_single, kappa(cols))

    def test_purity_breaks_at_b4(self):
        # b=4: no-singleton sum = kappa_4 + three pairing products kappa_2*kappa_2.
        # Use columns with guaranteed non-zero pairings so the extra term is real.
        W = [1, 1, -1, -1] * 4
        X = [1, 1, -1, -1] * 4          # X == W  -> kappa_2(W,X) != 0
        Y = [1, -1, 1, -1] * 4
        Z = [1, -1, 1, -1] * 4          # Z == Y  -> kappa_2(Y,Z) != 0
        cols = [W, X, Y, Z]
        no_single = Fraction(0)
        for part in set_partitions(range(4)):
            if any(len(block) == 1 for block in part):
                continue
            term = Fraction(1)
            for block in part:
                term *= kappa([cols[j] for j in block])
            no_single += term
        pairings = (kappa([W, X]) * kappa([Y, Z])
                    + kappa([W, Y]) * kappa([X, Z])
                    + kappa([W, Z]) * kappa([X, Y]))
        self.assertEqual(no_single, kappa(cols) + pairings)
        # The pairing part is non-zero, so b=4 is NOT a pure top cumulant.
        self.assertNotEqual(pairings, Fraction(0))


if __name__ == "__main__":
    unittest.main()
