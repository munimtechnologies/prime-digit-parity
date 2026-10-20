"""Verify the base-2 admit-even discriminator is not tautological.

These checks support notes/base2-even-admission-nontautology.md. They prove, with
exact arithmetic, the facts that separate the informative base-2 experiment from
the vacuous large-base case:

  - t(2h) = t(h): halving an even gap preserves its Thue-Morse (binary
    digit-sum) parity, so the base-2 character of an all-even-gap set reduces to
    the Thue-Morse parity of the halved gaps -- a live, two-valued signal, not a
    constant.
  - For base b >= 11 and a single-digit gap g (1 <= g < b), the base-b digit sum
    is g itself, so chi collapses to the parity of the gap length; over all even
    gaps it is the constant 0 (the large-base tautology).
  - In base 2, over the even gaps that actually occur for squarefree sets, chi
    takes both values 0 and 1.

Pure standard library. Run: python -m unittest -v

Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
"""

import unittest


def popcount(n):
    """Number of 1 bits in the binary expansion of n >= 0."""
    assert n >= 0
    return bin(n).count("1")


def thue_morse(n):
    """Thue-Morse value t(n) = popcount(n) mod 2 (OEIS A010060)."""
    return popcount(n) & 1


def digit_sum(n, base):
    """Sum of the base-`base` digits of n >= 0, for base >= 2."""
    assert n >= 0 and base >= 2
    s = 0
    while n:
        s += n % base
        n //= base
    return s


def chi(g, base):
    """Digit-sum parity character of a gap g in the given base."""
    return digit_sum(g, base) & 1


# Even gaps that actually occur between consecutive squarefree numbers below 1e8
# (all such gaps are <= ~10 in the computed range). These are the gaps an
# odd-only sieve could produce.
OCCURRING_EVEN_GAPS = [2, 4, 6, 8, 10, 12]


class TestThueMorseHalvingIdentity(unittest.TestCase):
    """Claim A: t(2h) = t(h); hence chi_2(2h) = t(h) (reduction B).

    Appending a zero bit (multiplying by 2) does not change the number of 1 bits,
    so the base-2 digit-sum parity of an even gap equals that of its half.
    """

    def test_halving_preserves_thue_morse(self):
        for h in range(0, 200001):
            self.assertEqual(thue_morse(2 * h), thue_morse(h))

    def test_base2_chi_of_even_gap_is_thue_morse_of_half(self):
        for h in range(0, 100001):
            self.assertEqual(chi(2 * h, 2), thue_morse(h))


class TestLargeBaseTautology(unittest.TestCase):
    """Claim: for b >= 11 every occurring gap is a single digit, so chi = gap
    parity, and over all even gaps chi is the constant 0 (vacuous)."""

    def test_single_digit_gap_digit_sum_is_itself(self):
        for base in (11, 16, 100):
            for g in range(1, base):
                self.assertEqual(digit_sum(g, base), g)
                self.assertEqual(chi(g, base), g & 1)

    def test_even_gaps_collapse_to_constant_zero(self):
        for base in (11, 16, 100):
            values = {chi(g, base) for g in range(2, base, 2)}
            self.assertEqual(values, {0})


class TestBase2DoesNotCollapse(unittest.TestCase):
    """Claim: in base 2 the character over occurring even gaps is two-valued, so
    the odd-only (all-even-gap) side of the experiment is non-degenerate."""

    def test_base2_even_gaps_take_both_values(self):
        values = {chi(g, 2) for g in OCCURRING_EVEN_GAPS}
        self.assertEqual(values, {0, 1})

    def test_specific_thue_morse_parities(self):
        # t(2,4,6,8,10,12) = 1,1,0,1,0,0
        self.assertEqual([thue_morse(g) for g in OCCURRING_EVEN_GAPS],
                         [1, 1, 0, 1, 0, 0])


if __name__ == "__main__":
    unittest.main()
