# test_hierarchy_results.py -- self-contained check of HIERARCHY-RESULTS.md.
#
# Measures the b=2 rung of the cumulant tower over every prime below 1e6:
# the two character factorisations, the published 1e6 excess anchors under the
# exact sign convention, the split kappa2 + m0*m1 = E[joint], and the dominance
# of the connected pairwise cumulant over the marginal product. Also re-checks
# arithmetically that the published 1e8 splits sum to the anchor excesses.
#
# Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.
#
# Stdlib only. No repo imports, no file reads. python -I safe.

import unittest


def sieve(limit):
    """Primes up to `limit` inclusive, as a list."""
    s = bytearray([1]) * (limit + 1)
    s[0] = s[1] = 0
    i = 2
    while i * i <= limit:
        if s[i]:
            s[i * i::i] = b"\x00" * len(s[i * i::i])
        i += 1
    return [i for i in range(2, limit + 1) if s[i]]


def build_mask(a, top_bit=40):
    """Bits at positions 0, a, 2a, ... up to top_bit."""
    m = 0
    k = 0
    while k * a <= top_bit:
        m |= 1 << (k * a)
        k += 1
    return m


MASK = {a: build_mask(a) for a in (1, 2, 3, 4)}


def chi(n, a):
    """chi_{2^a}(n) = (-1)^{s_{2^a}(n)} in {+1,-1}.

    By the positional lemma, s_{2^a}(n) mod 2 is the parity of the binary
    digits of n at positions == 0 (mod a).  popcount via bin().count for
    portability across Python versions.
    """
    return 1 - 2 * (bin(n & MASK[a]).count("1") & 1)


class HierarchyResults(unittest.TestCase):
    LIMIT = 10 ** 6
    # Published excess anchors at 1e6 (E_b = -sum_p chi_b(p)).
    ANCHORS = {2: 4764, 4: 2902, 8: -510, 16: 244}

    @classmethod
    def setUpClass(cls):
        cls.primes = sieve(cls.LIMIT)

    def test_prime_count_1e6(self):
        self.assertEqual(len(self.primes), 78498)  # pi(1e6)

    def test_character_factorisations_pointwise(self):
        # (star) at the two b=2 rungs, exact for every prime < 1e6.
        for p in self.primes:
            # chi_2(p) = chi_4(p) * chi_4(p>>1)
            self.assertEqual(chi(p, 1), chi(p, 2) * chi(p >> 1, 2))
            # chi_4(p) = chi_16(p) * chi_16(p>>2)
            self.assertEqual(chi(p, 2), chi(p, 4) * chi(p >> 2, 4))

    def test_excess_anchors_and_sign(self):
        # E_b = -sum_p chi_b(p) must match the published anchors exactly.
        for a, b in ((1, 2), (2, 4), (3, 8), (4, 16)):
            total = sum(chi(p, a) for p in self.primes)
            self.assertEqual(-total, self.ANCHORS[b],
                             f"base {b} excess mismatch")

    def _rung(self, a_joint, a_src, shift):
        """Measure a b=2 rung X0=chi_{2^a_src}(p), X1=chi(p>>shift) whose
        product is chi_{2^a_joint}(p).  Returns (m0, m1, Ejoint, k2, mp)."""
        n = len(self.primes)
        s0 = s1 = sj = 0
        for p in self.primes:
            x0 = chi(p, a_src)
            x1 = chi(p >> shift, a_src)
            s0 += x0
            s1 += x1
            sj += x0 * x1
        m0, m1, ej = s0 / n, s1 / n, sj / n
        k2 = ej - m0 * m1
        return m0, m1, ej, k2, m0 * m1

    def test_b2_split_identity_and_dominance(self):
        n = len(self.primes)
        # base 2 from base 4: a_src=2 (chi_4), shift=1.
        m0, m1, ej, k2, mp = self._rung(1, 2, 1)
        # the product X0*X1 must equal chi_2(p), so E[joint] == -E2/n.
        self.assertAlmostEqual(ej, -self.ANCHORS[2] / n, places=12)
        # exact moment-cumulant split for b=2.
        self.assertAlmostEqual(k2 + mp, ej, places=12)
        # connected cumulant dominates the marginal product.
        self.assertLess(abs(mp), abs(k2))
        self.assertLess(abs(mp / k2), 0.5)
        print(f"[base2_from_base4 @1e6] k2={k2:.8f} m0m1={mp:.8f} "
              f"ratio={mp / k2:.6f}")

        # base 4 from base 16: a_src=4 (chi_16), shift=2.
        m0, m1, ej, k2, mp = self._rung(2, 4, 2)
        self.assertAlmostEqual(ej, -self.ANCHORS[4] / n, places=12)
        self.assertAlmostEqual(k2 + mp, ej, places=12)
        self.assertLess(abs(mp), abs(k2))
        self.assertLess(abs(mp / k2), 0.5)
        print(f"[base4_from_base16 @1e6] k2={k2:.8f} m0m1={mp:.8f} "
              f"ratio={mp / k2:.6f}")

    def test_published_1e8_splits_sum_to_anchors(self):
        # The headline 1e8 numbers (from the sealed sandbox job) must satisfy
        # excess = connected_excess + marginal_excess exactly.
        rows = [
            # (excess, connected, marginal)
            (70509, 71506, -997),   # base2_from_base4
            (65633, 66230, -597),   # base4_from_base16
        ]
        for excess, connected, marginal in rows:
            self.assertEqual(excess, connected + marginal)
        # And the connected share is within [0.99, 1.02] of the excess.
        for excess, connected, _ in rows:
            share = connected / excess
            self.assertGreater(share, 0.99)
            self.assertLess(share, 1.02)


if __name__ == "__main__":
    unittest.main(verbosity=2)
