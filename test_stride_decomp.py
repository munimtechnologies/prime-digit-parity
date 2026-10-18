# SPDX-License-Identifier: Apache-2.0
"""Exhaustive verification of the stride-decomposition identities.

Pure Python, no third-party dependencies. Checks every integer in [0, N) with
N = 2^14, cross-checking an independent digit-sum implementation against the
bit-level formulas.
"""
import unittest

from stride_decomp import (
    digit_sum,
    P,
    chi,
    P_pow2_stride,
    P2_from_streams,
    base4_streams,
)

N = 1 << 14  # 16384: exhaustive over all 14-bit integers


def popcount_parity(g):
    """XOR of all bits of g, computed independently of stride_decomp."""
    x = 0
    while g:
        x ^= g & 1
        g >>= 1
    return x


def digit_sum_ref(g, b):
    """Independent reference base-b digit sum."""
    digits = []
    while g:
        digits.append(g % b)
        g //= b
    return sum(digits)


class TestStrideDecomposition(unittest.TestCase):
    def test_digit_sum_matches_reference(self):
        for b in (2, 3, 4, 5, 7, 8, 9, 16):
            for g in range(N):
                self.assertEqual(digit_sum(g, b), digit_sum_ref(g, b))

    def test_base2_is_thue_morse(self):
        # Identity 1: P_2(g) = popcount(g) mod 2 (Thue-Morse, OEIS A010060).
        for g in range(N):
            self.assertEqual(P(g, 2), popcount_parity(g))

    def test_pow2_stride_identity(self):
        # Identity (BIT-STRIDE): P_{2^k}(g) = XOR of bits 0, k, 2k, ...
        for k in (1, 2, 3, 4):
            b = 1 << k
            for g in range(N):
                self.assertEqual(P(g, b), P_pow2_stride(g, k))

    def test_decomposition(self):
        # NEW: P_2(g) = XOR_{r=0}^{k-1} P_{2^k}(g >> r), for k = 1..4.
        for k in (1, 2, 3, 4):
            for g in range(N):
                self.assertEqual(P(g, 2), P2_from_streams(g, k))

    def test_base4_corollary(self):
        # NEW: P_2(g) = P_4(g) XOR P_4(g >> 1).
        for g in range(N):
            a, b = base4_streams(g)
            self.assertEqual(P(g, 2), a ^ b)

    def test_base4_corollary_character_form(self):
        # chi_2(g) = chi_4(g) * chi_4(floor(g / 2)) in +/-1 form.
        for g in range(N):
            self.assertEqual(chi(g, 2), chi(g, 4) * chi(g >> 1, 4))

    def test_odd_base_collapse(self):
        # Identity 2: for any odd base b, P_b(g) = g mod 2 (gap-length parity),
        # because s_b(g) == g (mod b-1) and 2 | (b-1) when b is odd.
        for b in (3, 5, 7, 9, 11, 15):
            for g in range(N):
                self.assertEqual(P(g, b), g & 1)

    def test_chi_consistent_with_P(self):
        for b in (2, 4, 8, 16):
            for g in range(0, N, 7):  # sparse sample; relation is definitional
                self.assertEqual(chi(g, b), 1 if P(g, b) == 0 else -1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
