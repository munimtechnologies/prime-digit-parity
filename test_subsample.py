# Self-contained verification of the bit-subsampling identity (SUBSAMPLE.md).
# Stdlib only; no repo imports and no file reads, so cwd/path cannot affect it.
#
# Identity: for base b = 2**k, the parity of the base-b digit sum of n equals
# the parity of the number of 1-bits of n at positions 0, k, 2k, ...
#
# Work by Munim, Inc., produced with AI. Apache-2.0.

import unittest


def _digit_sum_parity(n, b):
    s = 0
    while n:
        s += n % b
        n //= b
    return s & 1


def _subsample_bit_parity(n, k):
    c = 0
    pos = 0
    while (n >> pos) != 0:
        c += (n >> pos) & 1
        pos += k
    return c & 1


class TestSubsampleIdentity(unittest.TestCase):

    def test_explicit_small_cases(self):
        # Hand-checked cases from SUBSAMPLE.md.
        # n=5 = 101_2 = "5" base 8: digit sum 5 (odd); bits at pos {0,3} = {1,0}.
        self.assertEqual(_digit_sum_parity(5, 8), 1)
        self.assertEqual(_subsample_bit_parity(5, 3), 1)
        # n=9 = 1001_2 = "11" base 8: digit sum 2 (even); bits at pos {0,3}={1,1}.
        self.assertEqual(_digit_sum_parity(9, 8), 0)
        self.assertEqual(_subsample_bit_parity(9, 3), 0)
        # n=6 = 110_2 = "12" base 4: digit sum 3 (odd); even bits {0,2} = {0,1}.
        self.assertEqual(_digit_sum_parity(6, 4), 1)
        self.assertEqual(_subsample_bit_parity(6, 2), 1)
        # n=64 = 2^6 = "100" base 8: digit sum 1 (odd); bit at pos 6 (==0 mod 3).
        self.assertEqual(_digit_sum_parity(64, 8), 1)
        self.assertEqual(_subsample_bit_parity(64, 3), 1)

    def test_identity_exhaustive(self):
        # Brute-force every n < 2^16 for the power-of-two bases 2,4,8,16,32,64.
        for k in range(1, 7):
            b = 1 << k
            for n in range(1 << 16):
                self.assertEqual(
                    _digit_sum_parity(n, b),
                    _subsample_bit_parity(n, k),
                    msg=f"identity failed at n={n}, k={k} (base {b})",
                )

    def test_base2_is_thue_morse(self):
        # k=1 selects every bit, so the parity is exactly the Thue-Morse word.
        # Thue-Morse t(n) for n=0..15.
        tm = [0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0]
        for n in range(16):
            self.assertEqual(_digit_sum_parity(n, 2), tm[n])
            self.assertEqual(_subsample_bit_parity(n, 1), tm[n])

    def test_density_is_one_over_k(self):
        # Sanity: among n < 2^18, the subsampled set size tracks ceil(bits/k).
        # Here we just confirm k=2 uses roughly half the bits that k=1 does by
        # comparing how many positions are inspected for a worst-case n.
        n = (1 << 18) - 1  # all 18 bits set
        inspected = {}
        for k in (1, 2, 3):
            pos, count = 0, 0
            while (n >> pos) != 0:
                count += 1
                pos += k
            inspected[k] = count
        self.assertEqual(inspected[1], 18)
        self.assertEqual(inspected[2], 9)
        self.assertEqual(inspected[3], 6)


if __name__ == "__main__":
    unittest.main()
