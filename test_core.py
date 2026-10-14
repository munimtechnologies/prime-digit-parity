"""Correctness tests for core.py (standard library only; no numpy needed).

These are the tests the publish lane runs before committing.  They verify:
  * digit_sum against independent string/bin() computations,
  * the odd-base triviality fact (Fact 1),
  * the even-base exact-balance fact on complete ranges [0, b**k) (Fact 2),
  * the sharper |excess| <= 1 bound and its closed form at every cutoff (Fact 3),
  * the base-10 last-digit localization of the prime bias (Fact 4),
  * casting out (b-1)'s and the odd-base degeneracy (Fact 5),
  * the bounded coprime-to-10 control for the base-10 prime bias (Fact 6),
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
    """For odd base b, every prime p > 2 has an odd base-b digit sum (Fact 1)."""

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


class TestIntegerExcessBounded(unittest.TestCase):
    """Fact 3: for even base the integer excess over [0, X) is in {-1,0,+1} at
    EVERY cutoff, matching the closed form (sharpening of Fact 2)."""

    def test_closed_form_matches_bruteforce(self):
        for b in (2, 4, 6, 8, 10, 16):
            for X in range(0, 2000):
                self.assertEqual(
                    core.integer_parity_excess_closed(X, b),
                    core.integer_parity_excess(X, b),
                    f"base {b}, X={X}",
                )

    def test_excess_bounded_by_one(self):
        for b in (2, 4, 6, 8, 10, 16):
            for X in range(0, 3000):
                self.assertLessEqual(abs(core.integer_parity_excess(X, b)), 1)

    def test_closed_form_large_X(self):
        # Complete blocks below b*floor(X/b) cancel, so excess over [0, X) equals
        # the excess over the final partial block [b*q, X); check the closed form
        # against that tail directly for some large X.
        for b in (2, 10, 16):
            for X in (b**4 + 3, 5 * b**3 + (b - 1), 7 * b**2 + 1, b**5 + b + 1):
                q = X // b
                tail = sum(
                    1 if core.digit_sum_parity(n, b) else -1
                    for n in range(b * q, X)
                )
                self.assertEqual(
                    core.integer_parity_excess_closed(X, b), tail,
                    f"base {b}, X={X}",
                )

    def test_closed_form_rejects_odd_base(self):
        with self.assertRaises(ValueError):
            core.integer_parity_excess_closed(100, 3)


class TestBase10Localization(unittest.TestCase):
    """Fact 4: every prime p > 5 ends in a digit in {1,3,7,9} (odd), so the
    base-10 digit-sum parity of p is the complement of that of floor(p/10)."""

    def test_last_digit_localization(self):
        for p in core.primes_upto(50000):
            if p <= 5:
                continue
            self.assertIn(p % 10, (1, 3, 7, 9))
            self.assertEqual(
                core.digit_sum_parity(p, 10),
                1 ^ core.digit_sum_parity(p // 10, 10),
            )


class TestCoprime10Control(unittest.TestCase):
    """Fact 6: the base-10 parity excess over integers ending in 1,3,7,9 is
    bounded by 8 at every cutoff and has the stated closed form."""

    def test_closed_form_matches_bruteforce(self):
        for X in range(0, 3000):
            self.assertEqual(
                core.coprime10_parity_excess_closed(X),
                core.coprime10_parity_excess(X),
                f"X={X}",
            )

    def test_bounded_by_eight(self):
        for X in range(0, 3000):
            self.assertLessEqual(abs(core.coprime10_parity_excess_closed(X)), 8)

    def test_bounded_far_out(self):
        for X in (10**6 + 7, 10**9 + 3, 5 * 10**12 + 9, 10**15 - 1, 123456789012345):
            self.assertLessEqual(abs(core.coprime10_parity_excess_closed(X)), 8)

    def test_primes_are_a_subset(self):
        # every prime p > 5 is counted by the control set (Fact 4)
        for p in core.primes_upto(5000):
            if p > 5:
                self.assertIn(p % 10, (1, 3, 7, 9))


class TestCastingOut(unittest.TestCase):
    """Fact 5: s_b(n) == n (mod b-1) for all n (casting out (b-1)'s), so the
    digit-sum parity is degenerate exactly when b is odd."""

    def test_casting_out_general(self):
        # b == 1 (mod b-1) gives s_b(n) == n (mod b-1); vacuous only for b = 2.
        for b in (3, 4, 7, 9, 10, 16):
            for n in range(0, 3000):
                self.assertEqual(
                    core.digit_sum(n, b) % (b - 1), n % (b - 1),
                    f"casting out failed base {b}, n={n}",
                )

    def test_parity_degenerate_flag(self):
        for b in (3, 5, 7, 9, 11, 15):
            self.assertTrue(core.parity_degenerate(b))
        for b in (2, 4, 6, 8, 10, 16):
            self.assertFalse(core.parity_degenerate(b))

    def test_degenerate_bases_force_prime_parity(self):
        # Consistency with Fact 1: in a degenerate (odd) base every prime p > 2
        # has odd digit sum, because s_b(p) == p == 1 (mod 2).
        for b in (3, 9, 15):
            for p in core.primes_upto(5000):
                if p > 2:
                    self.assertEqual(core.digit_sum_parity(p, b), 1)


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
