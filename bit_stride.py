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
"""Bit-stride explanation of the base-structured gap-parity lag profile.

This module makes precise *why* the lag-k serial correlation of prime-gap
digit-sum parity depends so strongly on the base (measured in LAGK-RESULTS.md /
release v0.9.0). The key is an elementary identity for digit-sum parity in a
power-of-two base, plus one number-theoretic fact about the size of prime gaps
below 1e8.

Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
"""

from __future__ import annotations


def digit_sum_parity(g: int, base: int) -> int:
    """Return (sum of base-`base` digits of g) mod 2, as 0 or 1."""
    if base < 2:
        raise ValueError("base must be >= 2")
    if g < 0:
        raise ValueError("g must be >= 0")
    s = 0
    while g:
        s ^= (g % base) & 1
        g //= base
    return s


def chi(g: int, base: int) -> int:
    """(-1) ** digit_sum_parity(g, base), i.e. +1 or -1."""
    return 1 - 2 * digit_sum_parity(g, base)


def bit_stride_parity(g: int, j: int) -> int:
    """XOR of the binary bits of g at positions 0, j, 2j, 3j, ... (mod 2).

    This equals digit_sum_parity(g, 2**j) by the bit-stride identity:
    the m-th base-2**j digit is floor(g / 2**(m*j)) mod 2**j, whose own
    parity is precisely bit (m*j) of g; summing digit parities mod 2 keeps
    only one bit per digit.
    """
    if j < 1:
        raise ValueError("j must be >= 1")
    s = 0
    pos = 0
    while (g >> pos) != 0:
        s ^= (g >> pos) & 1
        pos += j
    return s


def chi_power_of_two(g: int, j: int) -> int:
    """(-1) ** bit_stride_parity(g, j); equals chi(g, 2**j)."""
    return 1 - 2 * bit_stride_parity(g, j)


# --- Collapse corollaries, valid only for 0 <= g < 256 (all prime gaps
# --- below 1e8 are <= 220 < 256). For even g, bit 0 = 0, so only the high
# --- sampled bits remain. b_i(g) = (g >> i) & 1.

def collapse_parity(g: int, base: int) -> int:
    """Digit-sum parity of g in a power-of-two base, specialized to g < 256.

    Returns the same value as digit_sum_parity(g, base) for every g in
    [0, 256), written as an explicit sum of a few high binary bits so the
    dependence on gap *magnitude* is visible. Raises for g >= 256.
    """
    if not (0 <= g < 256):
        raise ValueError("collapse corollaries only hold for 0 <= g < 256")

    def b(i: int) -> int:
        return (g >> i) & 1

    if base == 2:          # j=1: every bit -> popcount parity (Thue-Morse)
        return bin(g).count("1") & 1
    if base == 4:          # j=2: bits 0,2,4,6
        return (b(0) ^ b(2) ^ b(4) ^ b(6))
    if base == 8:          # j=3: bits 0,3,6
        return (b(0) ^ b(3) ^ b(6))
    if base == 16:         # j=4: bits 0,4  -> for even g, just bit 4
        return (b(0) ^ b(4))
    raise ValueError("base must be one of 2, 4, 8, 16")


def base16_band_parity(g: int) -> int:
    """Parity of floor(g / 16). For every prime gap below 1e8 (even, <= 220),
    this equals digit_sum_parity(g, 16): base-16 'digit parity' is exactly a
    16-wide banding of the gap magnitude over the studied range."""
    return (g // 16) & 1


def sieve_primes(limit: int) -> list[int]:
    """Primes p with 2 <= p <= limit (simple Eratosthenes; sandbox use)."""
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0] = sieve[1] = 0
    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = b"\x00" * len(sieve[p * p::p])
    return [i for i in range(2, limit + 1) if sieve[i]]
