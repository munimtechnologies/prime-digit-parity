"""Self-contained exact-identity tests for the uniform even-base reduction.

These redefine every function inline (no repo import, no file reads, no sandbox
numbers): each test passes purely by an algebraic identity verified against
brute force. Work by Munim, Inc., produced with AI. Licence: Apache-2.0.
"""
import unittest

EVEN_BASES = [2, 4, 6, 8, 10, 12, 16]


def _digit_sum_parity(n, b):
    s = 0
    while n:
        s ^= (n % b) & 1
        n //= b
    return s


def _truncation_parity(n, b):
    s = 0
    q = n
    while q:
        s ^= q & 1
        q //= b
    return s


class TestEvenBaseReduction(unittest.TestCase):
    def test_identity_E_for_all_even_bases(self):
        # (E): s_b(n) == sum_{i>=0} floor(n/b^i) (mod 2), for even b.
        for b in EVEN_BASES:
            for n in range(0, 20000):
                self.assertEqual(
                    _digit_sum_parity(n, b), _truncation_parity(n, b),
                    msg=f"(E) failed at n={n}, b={b}")

    def test_identity_Eprime_product_form(self):
        # (E'): chi_b(n) = prod_{i>=0} eps(floor(n/b^i)), eps(x)=(-1)^x.
        for b in EVEN_BASES:
            for n in range(0, 15000):
                prod = 1
                q = n
                while q:
                    prod *= 1 - 2 * (q & 1)
                    q //= b
                chi = 1 - 2 * _digit_sum_parity(n, b)
                # for n == 0 the empty product is 1 and chi_b(0) = +1
                self.assertEqual(prod, chi, msg=f"(E') failed at n={n}, b={b}")

    def test_legendre_consistency(self):
        # (E) is the mod-2 shadow of Legendre:
        #   S_b(n) = n - (b-1) * sum_{i>=1} floor(n/b^i)   (the actual digit sum)
        for b in EVEN_BASES:
            for n in range(0, 6000):
                tail = 0
                q = n // b
                while q:
                    tail += q
                    q //= b
                legendre_sum = n - (b - 1) * tail
                self.assertEqual(legendre_sum % 2, _digit_sum_parity(n, b),
                                 msg=f"Legendre parity mismatch n={n}, b={b}")

    def test_prime_corollary_sign(self):
        # For odd n (every prime > 2), the i=0 factor eps(n) = -1, so
        #   chi_b(n) = - prod_{i>=1} eps(floor(n/b^i)).
        for b in EVEN_BASES:
            for n in range(1, 12000, 2):  # odd n only
                prod = 1
                q = n // b
                while q:
                    prod *= 1 - 2 * (q & 1)
                    q //= b
                chi = 1 - 2 * _digit_sum_parity(n, b)
                self.assertEqual(chi, -prod,
                                 msg=f"prime corollary failed n={n}, b={b}")

    def test_identity_is_even_specific(self):
        # (E) is specific to even b: exhibit a counterexample at an odd base,
        # so the identity is not vacuously true for every base.
        found = any(_digit_sum_parity(n, 3) != _truncation_parity(n, 3)
                    for n in range(0, 3000))
        self.assertTrue(found, "expected (E) to fail somewhere at base 3")

    def test_dyadic_factors_are_single_bits(self):
        # For b = 2^a the i-th truncation lsb equals bit a*i of n; for a
        # non-dyadic even base the factors are not single bits of n.
        for a in (1, 2, 3, 4):
            b = 1 << a
            for n in range(0, 4000):
                q = n
                i = 0
                while q:
                    self.assertEqual(q & 1, (n >> (a * i)) & 1)
                    q >>= a
                    i += 1


if __name__ == "__main__":
    unittest.main()
