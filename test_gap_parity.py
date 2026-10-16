"""Self-contained anchor checks for the prime-gap digit-sum parity result.

Re-derives, with an independent sieve to 1e6, the exact anchors that gate the
gap-parity measurement -- none of which depend on the (new, unknown) numeric
excesses, so this test validates the METHOD and passes standalone:

  * pi(1e6) = 78498;
  * telescoping: sum of gaps between consecutive primes <= 1e6 == p_M - 2;
  * exactly one odd gap (g_1 = 1);
  * half-gap Thue-Morse identity chi_2(g) == chi_2(g/2) for every even gap, so
    the half-gap E_2 tally (over A028334) equals the direct E_2;
  * two independent base-b digit-sum-parity routines agree for every gap and base.

Run: python -m unittest -v test_gap_parity
Work by Munim, Inc., produced with AI. MV2 of Munim, Inc. Apache-2.0.
"""
import unittest

BASES = (2, 4, 6, 8, 10, 16)
LIMIT = 1_000_000
PI_1E6 = 78498


def sieve(limit):
    s = bytearray(b"\x01") * (limit + 1)
    s[0:2] = b"\x00\x00"
    i = 2
    while i * i <= limit:
        if s[i]:
            s[i * i::i] = b"\x00" * len(range(i * i, limit + 1, i))
        i += 1
    return s


def parity_div(m, b):
    s = 0
    while m:
        s ^= (m % b) & 1
        m //= b
    return s


def parity_str(m, b):
    if m == 0:
        return 0
    out = 0
    while m:
        out += m % b
        m //= b
    return out & 1


class TestGapParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        s = sieve(LIMIT)
        cls.primes = [n for n in range(2, LIMIT + 1) if s[n]]

    def test_pi(self):
        self.assertEqual(len(self.primes), PI_1E6)

    def test_anchors_and_identities(self):
        primes = self.primes
        gaps = [primes[i + 1] - primes[i] for i in range(len(primes) - 1)]
        # telescoping
        self.assertEqual(sum(gaps), primes[-1] - 2)
        # exactly one odd gap, which is g_1 = 1
        odd = [g for g in gaps if g & 1]
        self.assertEqual(odd, [1])
        # half-gap identity + two-routine guard + direct E_b tally
        excess = {b: 0 for b in BASES}
        half_e2 = 0
        for g in gaps:
            for b in BASES:
                pd = parity_div(g, b)
                self.assertEqual(pd, parity_str(g, b),
                                 msg=f"routine mismatch g={g} b={b}")
                excess[b] += 1 if pd else -1
            h = g >> 1 if (g & 1) == 0 else g
            # Thue-Morse: s_2(2h) == s_2(h)
            if g & 1 == 0:
                self.assertEqual(parity_div(g, 2), parity_div(h, 2))
            half_e2 += 1 if parity_div(h, 2) else -1
        self.assertEqual(half_e2, excess[2])

    def test_chi2_halfgap_small(self):
        # chi_2(2h) = chi_2(h) directly on small even numbers.
        for h in range(0, 50000):
            self.assertEqual(parity_div(2 * h, 2), parity_div(h, 2))


if __name__ == "__main__":
    unittest.main(verbosity=2)
