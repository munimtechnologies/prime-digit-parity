"""Self-contained tests for the even/odd decade law (DECADES.md, decades.py).

Work by Munim, Inc., produced with AI. License: Apache-2.0.

No repo imports and no file reads: every constant is a gated sandbox output
embedded here, and every assertion is an exact algebraic identity or an
exhaustive check. If these pass, the numbers quoted in DECADES.md are internally
consistent with one another and with the stated mathematical relations.
"""

import unittest
from math import sqrt

# Gated value-cutoff data (same literals as decades.py / the CSV).
EXCESS = {1: 2, 2: -1, 3: 6, 4: 75, 5: 112,
          6: 1916, 7: 2967, 8: 62627, 9: 62582}
NPRIMES = {6: 78498, 7: 664579, 8: 5761455, 9: 50847534,
           10: 455052511, 11: 4118054813}
Z_PUB = {6: 6.84, 7: 3.64, 8: 26.09, 9: 8.78, 10: 98.93}


def digit_sum(n):
    s = 0
    while n:
        s += n % 10
        n //= 10
    return s


def chi10(n):
    return -1 if digit_sum(n) & 1 else 1


class TestZConsistency(unittest.TestCase):
    """z = excess / sqrt(n) must reproduce every published z from the excess."""

    def test_z_matches_excess_over_sqrt_n(self):
        for m in (6, 7, 8, 9):
            z = EXCESS[m] / sqrt(NPRIMES[m])
            self.assertAlmostEqual(z, Z_PUB[m], places=1)

    def test_1e9_is_a_plateau_not_a_fall(self):
        # The raw excess barely moves across 10^8 -> 10^9 ...
        self.assertLess(abs(EXCESS[9] - EXCESS[8]), 100)
        # ... yet z drops, purely because sqrt(n) grows ~8.8x.
        self.assertGreater(NPRIMES[9] / NPRIMES[8], 8.0)
        self.assertLess(Z_PUB[9], Z_PUB[8] / 2)


class TestEvenOddDecadeLaw(unittest.TestCase):
    """The discovery: ratio > 10 iff m even, ratio < 2 iff m odd, for m=4..9."""

    def ratio(self, m):
        return EXCESS[m] / EXCESS[m - 1]

    def test_even_decades_jump(self):
        for m in (4, 6, 8):
            self.assertGreater(self.ratio(m), 10.0,
                               f"decade ending 10^{m} should be a jump")

    def test_odd_decades_flat(self):
        for m in (5, 7, 9):
            self.assertLess(self.ratio(m), 2.0,
                            f"decade ending 10^{m} should be flat")

    def test_law_is_exactly_even_vs_odd(self):
        for m in range(4, 10):
            r = self.ratio(m)
            is_jump = r > 10.0
            self.assertEqual(is_jump, (m % 2 == 0),
                             f"m={m} ratio={r} broke the even/odd law")

    def test_even_increment_dwarfs_preceding_odd_increment(self):
        # Each even decade adds far more new excess than the odd decade before.
        for m in (6, 8):
            even_incr = EXCESS[m] - EXCESS[m - 1]
            odd_incr = EXCESS[m - 1] - EXCESS[m - 2]
            self.assertGreater(abs(even_incr), 20 * abs(odd_incr))


class TestPredictionMethod(unittest.TestCase):
    """Flat-decade z prediction, validated on the known 10^8 -> 10^9 decade."""

    def test_flat_prediction_reproduces_1e9(self):
        pred = Z_PUB[8] * sqrt(NPRIMES[8] / NPRIMES[9])
        self.assertAlmostEqual(pred, Z_PUB[9], places=1)

    def test_preregistered_1e11_prediction_is_a_big_drop(self):
        # [10^10, 10^11) is an ODD decade; if the law holds the excess stays put
        # and z should fall well below its 10^10 value of 98.93.
        pred = Z_PUB[10] * sqrt(NPRIMES[10] / NPRIMES[11])
        self.assertTrue(30.0 < pred < 36.0)
        self.assertLess(pred, Z_PUB[10] / 2)


class TestLeadingDigitAlgebra(unittest.TestCase):
    """Exact identities behind the mechanism discussion in DECADES.md."""

    def test_digit_sum_is_additive_over_leading_digit(self):
        # chi_10(p) = (-1)^d * chi_10(p - d*10^(k)) where d is the leading digit.
        # Checked exhaustively for all n in [1, 20000].
        for n in range(1, 20001):
            k = len(str(n)) - 1
            d = n // (10 ** k)
            low = n - d * (10 ** k)
            self.assertEqual(chi10(n), ((-1) ** d) * chi10(low))

    def test_leading_digit_sign_sum(self):
        # sum_{d=1..9} (-1)^d = -1: four even d give +1, five odd d give -1.
        self.assertEqual(sum((-1) ** d for d in range(1, 10)), -1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
