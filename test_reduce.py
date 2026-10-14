# Copyright 2026 Munim, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Work by Munim, Inc., produced with AI.
"""Exhaustive, self-contained test of the inter-base reduction identity.

Self-contained on purpose: every function it checks is redefined inline, so the
test depends on nothing but the Python standard library and cannot fail on an
import path or working directory. It re-derives s_b(n) parity from scratch (plain
base conversion) and checks the reduction against it for all n < 2**16.
"""

import unittest


def digit_sum_parity(n, base):
    """Parity (0/1) of the sum of base-`base` digits of n, from first principles."""
    s = 0
    while n > 0:
        s += n % base
        n //= base
    return s & 1


def reduction_rhs(n, a, b):
    """sum_{r=0}^{b-1} s_{2**(a*b)}( n >> (r*a) )  (mod 2)."""
    big = 1 << (a * b)
    total = 0
    for r in range(b):
        total += digit_sum_parity(n >> (r * a), big)
    return total & 1


def subsample_parity(n, k):
    """Parity of the binary 1-bits of n at positions 0, k, 2k, ... ."""
    bits = 0
    pos = 0
    while n > 0:
        if pos % k == 0:
            bits += n & 1
        n >>= 1
        pos += 1
    return bits & 1


class TestInterBaseReduction(unittest.TestCase):
    LIMIT = 1 << 16  # 65536

    def test_reduction_exhaustive(self):
        # For every pair (a, b), s_{2**a}(n) must equal the reduction RHS for all n.
        pairs = [(1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (3, 2), (1, 5), (2, 4)]
        for a, b in pairs:
            small = 1 << a
            for n in range(self.LIMIT):
                lhs = digit_sum_parity(n, small)
                rhs = reduction_rhs(n, a, b)
                self.assertEqual(
                    lhs, rhs, msg=f"a={a} b={b} n={n}: lhs={lhs} rhs={rhs}"
                )

    def test_subsample_characterization(self):
        # Cross-check the underlying subsample identity (SUBSAMPLE.md) directly.
        for k in range(1, 7):
            base = 1 << k
            for n in range(1 << 14):
                self.assertEqual(digit_sum_parity(n, base), subsample_parity(n, k))

    def test_named_special_cases(self):
        for n in range(self.LIMIT):
            # base 4 = s_16(n) XOR s_16(n // 4)
            self.assertEqual(
                digit_sum_parity(n, 4),
                (digit_sum_parity(n, 16) + digit_sum_parity(n >> 2, 16)) & 1,
            )
            # Thue-Morse (base 2) = XOR of three shifted base-8 parities
            self.assertEqual(
                digit_sum_parity(n, 2),
                (
                    digit_sum_parity(n, 8)
                    + digit_sum_parity(n >> 1, 8)
                    + digit_sum_parity(n >> 2, 8)
                )
                & 1,
            )

    def test_degenerate_b_equals_1(self):
        # b = 1 is the identity: the RHS is just s_{2**a}(n) itself.
        for a in (1, 2, 3):
            for n in range(1 << 12):
                self.assertEqual(digit_sum_parity(n, 1 << a), reduction_rhs(n, a, 1))


if __name__ == "__main__":
    unittest.main()
