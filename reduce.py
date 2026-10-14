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
"""Inter-base reduction for power-of-two digit-sum parity.

Companion to SUBSAMPLE.md. The subsample identity states that for b = 2**k the
parity of the base-b digit sum of n equals the parity of the binary 1-bits of n
at positions 0, k, 2k, ... (a density-1/k subsample of the binary expansion).

This module adds the *inter-base* reduction, which expresses the digit-sum parity
in a smaller power-of-two base as an exact XOR-sum of the parity in a larger
power-of-two base, evaluated at right-shifted copies of n:

    s_{2**a}(n) == sum_{r=0}^{b-1} s_{2**(a*b)}( n >> (r*a) )   (mod 2)

Everything here is stdlib-only and exact; see test_reduce.py for an exhaustive
check over all n < 2**16 and several (a, b) pairs.
"""

from __future__ import annotations


def digit_sum_parity(n: int, base: int) -> int:
    """Return the parity (0 or 1) of the sum of the base-`base` digits of n.

    n must be a non-negative integer and base >= 2.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if base < 2:
        raise ValueError("base must be at least 2")
    s = 0
    while n > 0:
        s += n % base
        n //= base
    return s & 1


def reduction_rhs(n: int, a: int, b: int) -> int:
    """Return the right-hand side of the inter-base reduction (mod 2).

    sum_{r=0}^{b-1} s_{2**(a*b)}( floor(n / 2**(r*a)) )  (mod 2)

    By the theorem this equals digit_sum_parity(n, 2**a) for every n >= 0.
    """
    if a < 1 or b < 1:
        raise ValueError("a and b must be positive integers")
    big = 1 << (a * b)
    total = 0
    for r in range(b):
        total += digit_sum_parity(n >> (r * a), big)
    return total & 1


def base4_from_base16(n: int) -> int:
    """s_4(n) parity expressed via base 16: s_16(n) XOR s_16(floor(n/4))."""
    return (digit_sum_parity(n, 16) + digit_sum_parity(n >> 2, 16)) & 1


def thue_morse_from_base8(n: int) -> int:
    """s_2(n) parity (the Thue-Morse value) as an XOR of three shifted base-8
    parities: s_8(n) XOR s_8(floor(n/2)) XOR s_8(floor(n/4))."""
    return (
        digit_sum_parity(n, 8)
        + digit_sum_parity(n >> 1, 8)
        + digit_sum_parity(n >> 2, 8)
    ) & 1


if __name__ == "__main__":
    # Tiny self-check for every n < 2**14 and the two named special cases.
    for m in range(1 << 14):
        assert digit_sum_parity(m, 4) == base4_from_base16(m)
        assert digit_sum_parity(m, 2) == thue_morse_from_base8(m)
    print("inter-base reduction self-check passed for n < 2**14")
