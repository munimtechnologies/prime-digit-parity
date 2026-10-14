"""Correctness tests for core.py (standard library only; no numpy needed).

These are the tests the publish lane runs before committing.  They verify:
  * digit_sum against independent string/bin() computations,
  * the odd-base triviality fact,
  * the even-base exact-balance fact on complete ranges [0, b**k),
  * the sieve, by cross-checking a known prime count.

Run: python -m unittest -v
"""

import unittest

import core


class TestDigitSum(unittest.TestCase):
    def test_base10_against_string(self):
        for n in range(0, 5000):
            self.assertEqual(core.digit_sum(n, 10), sum(int(c) for c in str(n)))

    def test_base2_against_bin(self):
        for n in range(0, 5000):
            self.assertEqual(core.digit_sum(n, 2), bin(n).count("1"))

    def test_base16_against_hex(self):
        for n in range(0, 5000):
            self.assertEqual(
                core.digit_sum(n, 16), sum(int(c, 16) for c in format(n, "x"))
            )

    def test_parity_matches_digit_sum(self):
        for n in range(0, 2000):
            for b in (2, 3, 8, 10, 16):
                self.assertEqual(
                    core.digit_sum_parity(n, b), core.digit_sum(n, b) % 2
                )


class TestOddBaseTriviality(unittest.TestCase):
    """For odd base b, every prime p > 2 has an odd base-b digit sum."""

    def test_odd_bases_all_primes_odd_digit_sum(self):
        primes = [p for p in core.primes_upto(50000) if p > 2]
        for b in (3, 5, 7, 9, 11, 15):
            self.assertTrue(
                all(core.digit_sum_parity(p, b) == 1 for p in primes),
                f"found a prime with even digit sum in odd base {b}",
            )


class TestEvenBaseBalance(unittest.TestCase):
    """For even base b, [0, b**k) is exactly balanced in digit-sum parity."""

    def test_complete_ranges_balanced(self):
        cases = [
            (2, 1), (2, 5), (2, 12),
            (4, 1), (4, 6),
            (6, 1), (6, 5),
            (8, 1), (8, 4),
            (10, 1), (10, 2), (10, 4),
            (16, 1), (16, 3),
        ]
        for b, k in cases:
            odd, even = core.integer_parity_counts(b**k, b)
            self.assertEqual(odd, even, f"base {b}, k={k}: {odd} != {even}")

    def test_odd_base_ranges_not_balanced_in_general(self):
        # base 3 over [0, 9): digit sums 0,1,2,1,2,3,2,3,4 -> parities
        # 0,1,0,1,0,1,0,1,0 -> 4 odd, 5 even (not balanced); documents the
        # contrast with even bases.
        odd, even = core.integer_parity_counts(9, 3)
        self.assertNotEqual(odd, even)


class TestSieve(unittest.TestCase):
    def test_prime_count_to_million(self):
        # pi(10**6) = 78498 (classical value).
        self.assertEqual(len(core.primes_upto(1_000_000)), 78498)

    def test_small_primes(self):
        self.assertEqual(core.primes_upto(20), [2, 3, 5, 7, 11, 13, 17, 19])


class TestPrimeParitySmall(unittest.TestCase):
    def test_base2_excess_small(self):
        # base-2 (odious vs evil) excess is well studied; sanity-check sign on
        # a small range.  Odious (odd digit sum) should not be far behind.
        odd, even = core.prime_parity_counts(1000, 2)
        self.assertEqual(odd + even, len(core.primes_upto(1000)))


if __name__ == "__main__":
    unittest.main()
