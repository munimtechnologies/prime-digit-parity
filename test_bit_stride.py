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
"""Self-contained proof-by-enumeration for BIT-STRIDE.md (release v0.10.0).

Stdlib only; no repo imports; no file reads. Verifies the bit-stride identity
for power-of-two bases, the g<256 collapse corollaries, and the base-16 =
floor(g/16)-parity reduction over real prime gaps. Run: python -m unittest -v
"""

import unittest


def digit_sum_parity(g, base):
    s = 0
    while g:
        s ^= (g % base) & 1
        g //= base
    return s


def bit_stride_parity(g, j):
    s = 0
    pos = 0
    while (g >> pos) != 0:
        s ^= (g >> pos) & 1
        pos += j
    return s


def sieve_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0] = sieve[1] = 0
    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = b"\x00" * len(sieve[p * p::p])
    return [i for i in range(2, limit + 1) if sieve[i]]


class TestBitStrideIdentity(unittest.TestCase):
    """chi_{2^j}(g) digit-sum parity == XOR of bits at stride j, for ALL g."""

    def test_identity_all_g_below_65536(self):
        for j in (1, 2, 3, 4, 5):
            base = 2 ** j
            for g in range(0, 1 << 16):
                self.assertEqual(
                    digit_sum_parity(g, base),
                    bit_stride_parity(g, j),
                    msg="mismatch at g=%d, base=%d" % (g, base),
                )


class TestCollapseCorollaries(unittest.TestCase):
    """For 0 <= g < 256 the digit-sum parity reduces to a few high bits."""

    def b(self, g, i):
        return (g >> i) & 1

    def test_base2_is_popcount(self):
        for g in range(256):
            self.assertEqual(digit_sum_parity(g, 2), bin(g).count("1") & 1)

    def test_base4_bits_0_2_4_6(self):
        for g in range(256):
            expect = self.b(g, 0) ^ self.b(g, 2) ^ self.b(g, 4) ^ self.b(g, 6)
            self.assertEqual(digit_sum_parity(g, 4), expect)

    def test_base8_bits_0_3_6(self):
        for g in range(256):
            expect = self.b(g, 0) ^ self.b(g, 3) ^ self.b(g, 6)
            self.assertEqual(digit_sum_parity(g, 8), expect)

    def test_base16_bits_0_4(self):
        for g in range(256):
            expect = self.b(g, 0) ^ self.b(g, 4)
            self.assertEqual(digit_sum_parity(g, 16), expect)

    def test_base16_even_g_is_band_parity(self):
        # even g < 256: bit 0 = 0, so chi_16 == parity of floor(g/16)
        for g in range(0, 256, 2):
            self.assertEqual(digit_sum_parity(g, 16), (g // 16) & 1)


class TestOnRealPrimeGaps(unittest.TestCase):
    """Tie the reduction to the actual object: prime gaps."""

    def setUp(self):
        self.primes = sieve_primes(1_000_000)
        self.gaps = [self.primes[i + 1] - self.primes[i]
                     for i in range(len(self.primes) - 1)]

    def test_all_gaps_below_256(self):
        # Max prime gap below 1e6 is 114; well under 256 (= 2**8).
        self.assertLess(max(self.gaps), 256)

    def test_base16_parity_equals_band_parity_on_all_gaps(self):
        # The headline reduction: base-16 gap digit-sum parity IS floor(g/16)
        # parity for every (even) gap in the range. The first gap (3-2=1) is
        # odd, so bit 0 is set and it is excluded -- every later gap is even.
        for g in self.gaps[1:]:
            self.assertEqual(digit_sum_parity(g, 16), (g // 16) & 1)

    def test_even_gaps_except_first(self):
        # Every gap after the first (3-2=1) is even, so bit 0 = 0.
        self.assertEqual(self.gaps[0], 1)
        for g in self.gaps[1:]:
            self.assertEqual(g & 1, 0)


if __name__ == "__main__":
    unittest.main()
