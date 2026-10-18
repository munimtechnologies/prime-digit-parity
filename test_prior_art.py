"""Self-contained checks for the technical identities asserted in PRIOR-ART.md.

No third-party dependencies, no network. Run: python -m unittest -v

The prior-art note positions the "serial gap-parity" statistic against the
literature and makes three small exact claims about the digit-sum-parity
character chi_b. These tests verify each of those claims by direct computation
so the note is grounded rather than merely asserted.

Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
"""
import unittest


def digit_sum(n, base):
    """Sum of the base-`base` digits of the non-negative integer n."""
    assert n >= 0 and base >= 2
    s = 0
    while n:
        s += n % base
        n //= base
    return s


def thue_morse(n):
    """The Thue-Morse value t(n): parity of the number of 1 bits in n."""
    return bin(n).count("1") & 1


def chi(g, base):
    """The +/-1 character used throughout the arc: (-1)^(base-b digit sum of g)."""
    return 1 if digit_sum(g, base) % 2 == 0 else -1


def strided_bit_parity(n, k):
    """Parity of the binary digits of n at positions 0, k, 2k, 3k, ..."""
    s = 0
    pos = 0
    while n:
        if pos % k == 0:
            s += n & 1
        n >>= 1
        pos += 1
    return s % 2


# First 32 terms of the Thue-Morse sequence, OEIS A010060.
A010060 = [0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
           1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1]


class TestThueMorseIdentity(unittest.TestCase):
    """Claim 1: in base 2, chi_2(g) = (-1)^t(g), the Thue-Morse sequence."""

    def test_base2_parity_equals_thue_morse(self):
        for g in range(0, 500):
            self.assertEqual(digit_sum(g, 2) % 2, thue_morse(g))

    def test_matches_oeis_a010060(self):
        for n, expected in enumerate(A010060):
            self.assertEqual(thue_morse(n), expected)

    def test_chi_base2_is_pm_one_thue_morse(self):
        for g in range(0, 500):
            self.assertEqual(chi(g, 2), (-1) ** thue_morse(g))


class TestOddBaseCollapse(unittest.TestCase):
    """Claim 2: for any ODD base b, s_b(g) = g (mod 2), so chi_b carries no
    base-specific information -- it is just the parity of the gap length.
    (This is why only even bases are studied; see THEORY.md.)"""

    def test_odd_base_digit_sum_parity_equals_value_parity(self):
        for base in (3, 5, 7, 9, 11, 15, 21):
            for n in range(0, 600):
                self.assertEqual(digit_sum(n, base) % 2, n % 2)


class TestPowerOfTwoSubsampling(unittest.TestCase):
    """Claim 3: in base 2^k, s_{2^k}(g) = parity of the binary digits of g at
    positions 0, k, 2k, ... -- a strided subsample of the binary expansion.
    Hence comparing bases 2, 4, 8, 16 compares serial correlations of strided
    Thue-Morse subsamples, making the 'subsequence' framing precise."""

    def test_power_of_two_is_strided_binary_parity(self):
        for k in (1, 2, 3, 4):
            base = 2 ** k
            for n in range(0, 600):
                self.assertEqual(digit_sum(n, base) % 2, strided_bit_parity(n, k))

    def test_base4_concrete_values(self):
        # 6 = '12' base 4 -> digit sum 3 (odd); bits 0,2 of 110b are 0,1 -> 1.
        self.assertEqual(digit_sum(6, 4) % 2, 1)
        self.assertEqual(strided_bit_parity(6, 2), 1)
        # 5 = '11' base 4 -> digit sum 2 (even); bits 0,2 of 101b are 1,1 -> 0.
        self.assertEqual(digit_sum(5, 4) % 2, 0)
        self.assertEqual(strided_bit_parity(5, 2), 0)


if __name__ == "__main__":
    unittest.main()
