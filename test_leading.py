# test_leading.py -- self-contained checks for the positional (leading-digit)
# decomposition of the base-10 prime digit-sum excess.
#
# Work by Munim, Inc., produced with AI. Apache-2.0.
#
# Everything is redefined inline (no repo imports, no file reads). The tests
# assert exact algebraic identities and gate a stdlib sieve against the
# published anchors E(10^m) for m <= 6 (E(10^3)=6, 10^4=75, 10^5=112,
# 10^6=1916), which independently reproduces Gupta & Sridhar's phenomenon.

import unittest


def s10(n):
    t, acc = n, 0
    while t:
        acc += t % 10
        t //= 10
    return acc


def chi(n):
    return 1 if s10(n) % 2 == 0 else -1


def chi_via_truncations(n):
    prod, t = 1, n
    while True:
        if t % 2:
            prod = -prod
        if t == 0:
            break
        t //= 10
    return prod


def chi_leading(n):
    if n == 0:
        return 1
    m = len(str(n))
    p = 10 ** (m - 1)
    return (1 if (n // p) % 2 == 0 else -1) * chi(n % p)


def primes_upto(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    i = 2
    while i * i <= limit:
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        i += 1
    return [i for i in range(2, limit + 1) if sieve[i]]


class TestCharacterIdentities(unittest.TestCase):
    def test_truncation_product_matches_digit_sum(self):
        # chi(n) = prod_{i>=0} (-1)^{floor(n/10^i)} for all n < 100000.
        for n in range(0, 100000):
            self.assertEqual(chi(n), chi_via_truncations(n))

    def test_leading_digit_recursion(self):
        # chi(n) = (-1)^d chi(n mod 10^(m-1)) for all n < 100000.
        for n in range(0, 100000):
            self.assertEqual(chi(n), chi_leading(n))

    def test_units_factor_is_constant_on_odd_numbers(self):
        # For odd n, chi(n) = - prod_{i>=1} (-1)^{floor(n/10^i)}: the units
        # digit of an odd number is odd, contributing a fixed -1.
        for n in range(1, 100000, 2):
            higher = 1
            t = n // 10
            while True:
                if t % 2:
                    higher = -higher
                if t == 0:
                    break
                t //= 10
            self.assertEqual(chi(n), -higher)


class TestIntegerBaseline(unittest.TestCase):
    def test_full_decade_character_sum_is_zero(self):
        # sum_{n in [10^(m-1), 10^m)} chi(n) = 0 for m = 2..6.
        for m in range(2, 7):
            lo, hi = 10 ** (m - 1), 10 ** m
            self.assertEqual(sum(chi(n) for n in range(lo, hi)), 0)


class TestPrimeExcessDecomposition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.primes = primes_upto(10 ** 6)

    def _excess_upto(self, x):
        acc = 0
        for p in self.primes:
            if p > x:
                break
            acc -= chi(p)
        return acc

    def test_published_anchors(self):
        # Gate the sieve against the published base-10 excess anchors.
        self.assertEqual(self._excess_upto(10 ** 3), 6)
        self.assertEqual(self._excess_upto(10 ** 4), 75)
        self.assertEqual(self._excess_upto(10 ** 5), 112)
        self.assertEqual(self._excess_upto(10 ** 6), 1916)

    def test_leading_digit_decomposition_reconstructs_each_decade(self):
        # E_decade(m) = E(10^m) - E(10^(m-1)) = sum_{d=1}^9 (-1)^(d+1) T_d(m),
        # with T_d(m) = sum over primes p in [10^(m-1),10^m) of leading digit
        # d of chi(p - d*10^(m-1)).
        for m in range(2, 7):
            lo, hi = 10 ** (m - 1), 10 ** m
            base = 10 ** (m - 1)
            T = {d: 0 for d in range(1, 10)}
            for p in self.primes:
                if p < lo:
                    continue
                if p >= hi:
                    break
                d = p // base
                T[d] += chi(p - d * base)
            from_profile = sum((1 if d % 2 else -1) * T[d] for d in range(1, 10))
            direct = self._excess_upto(hi - 1) - self._excess_upto(lo - 1)
            self.assertEqual(from_profile, direct)

    def test_decade_parity_law_holds_in_measured_range(self):
        # The v0.4.6 law: |E_decade(m)| is large when m is even and small
        # when m is odd. Check the ratio E(10^m)/E(10^(m-1)) for the
        # verified decades m = 4, 5, 6 (even m -> big jump; odd m -> flat).
        e = {m: self._excess_upto(10 ** m) for m in range(3, 7)}
        # even m = 4, 6: at least a 5x jump
        self.assertGreater(e[4] / e[3], 5.0)
        self.assertGreater(e[6] / e[5], 5.0)
        # odd m = 5: nearly flat (ratio below 2)
        self.assertLess(e[5] / e[4], 2.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
