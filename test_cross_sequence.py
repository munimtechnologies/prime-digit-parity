"""Self-contained gating test for CROSS-SEQUENCE.md.

Work by Munim, Inc., produced with AI (MV2 of Munim, Inc.). Apache-2.0.

This test is INTENTIONALLY self-contained: stdlib only, no repo imports, no file
reads. It re-sieves to 1e6 and recomputes the base-8 digit-sum parity excess E8 for
four sequences -- primes, squarefree numbers, semiprimes (Omega==2), and all integers
-- so any reader can verify the 1e6 slice of the cross-sequence comparison on their own
machine in about a second, independently of the compute-sandbox job that measures the
larger scales.

Four assertions are hard facts that must hold on publish:

1. PRIME ANCHOR. E8(primes <= 1e6) == -510. Measured by the prime-parity-trajectory-v1
   sandbox job and recorded as a gate anchor; reproducing it here cross-validates this
   file's sieve + octal-digit-parity logic end to end.

2. INTEGER NULL IDENTITY (hand-proved, exact). E8(all integers in [0, N)) == 0 for every
   N divisible by 8. Proof: write n = 8a + d with 0 <= d < 8; then s_8(n) = s_8(a) + d
   (d is the last octal digit, no carry), so chi_8(n) = chi_8(a) * (-1)^d. Summing
   d = 0..7 gives sum (-1)^d = 0, hence sum_{n=0}^{8k-1} chi_8(n) = 0, and E8 = -sum = 0.
   N = 10**6 is divisible by 8, so the integer null excess is EXACTLY zero over [0, N).
   (The sandbox counts [1, N] and so reports the boundary term 1 - chi_8(1e6) = 2; this
   test fixes the [0, N) convention in which the null is exactly 0.)

3. SQUAREFREE 1e6 SLICE. E8(squarefree <= 1e6) == 13272, matching cross-seq-parity-v1.

4. SEMIPRIME 1e6 SLICE. E8(Omega==2, n <= 1e6) == 2449, matching cross-seq-parity-v1.

Sign convention (repo-wide): E_b = #odd - #even = -sum_p chi_b, chi_b(n) = (-1)^{s_b(n)}.
"""
import unittest

N = 10 ** 6  # divisible by 8


def s8_parity(n):
    """Parity (0/1) of the base-8 digit sum of n."""
    s = 0
    while n:
        s += n & 7
        n >>= 3
    return s & 1


def chi8(n):
    return 1 if s8_parity(n) == 0 else -1


def sieves(limit):
    """Return (is_prime, omega, is_squarefree) arrays over [0, limit).

    omega[n] = number of prime factors of n counted WITH multiplicity (Omega).
    is_squarefree[n] = 1 iff no prime square divides n (1 counts as squarefree).
    """
    is_prime = bytearray([1]) * limit
    is_prime[0] = 0
    if limit > 1:
        is_prime[1] = 0
    omega = [0] * limit
    has_square = bytearray(limit)  # 1 iff divisible by some p^2
    for p in range(2, limit):
        if is_prime[p]:
            for m in range(p + p, limit, p):
                is_prime[m] = 0
            pk = p
            while pk < limit:
                for m in range(pk, limit, pk):
                    omega[m] += 1
                    if pk != p:  # pk = p^e, e >= 2  ->  m divisible by a square
                        has_square[m] = 1
                pk *= p
    is_squarefree = bytearray([1]) * limit
    is_squarefree[0] = 0
    for n in range(1, limit):
        if has_square[n]:
            is_squarefree[n] = 0
    return is_prime, omega, is_squarefree


def e8_over(predicate, limit):
    """E8 = #odd - #even = -sum chi8, over n in [0, limit) with predicate(n) true."""
    s = 0
    for n in range(limit):
        if predicate(n):
            s += chi8(n)
    return -s


class CrossSequenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.is_prime, cls.omega, cls.is_squarefree = sieves(N)

    def test_prime_anchor(self):
        e8_primes = e8_over(lambda n: self.is_prime[n], N)
        self.assertEqual(e8_primes, -510,
                         "base-8 prime parity excess at 1e6 must match the gate anchor")

    def test_integer_null_identity(self):
        # exact zero at any multiple of 8, by the proof in the module docstring
        self.assertEqual(N % 8, 0)
        total = sum(chi8(n) for n in range(N))
        self.assertEqual(total, 0)
        self.assertEqual(-total, 0)  # E8(integers in [0, N)) == 0 exactly

    def test_squarefree_1e6_slice(self):
        e8_sqfree = e8_over(lambda n: self.is_squarefree[n], N)
        self.assertEqual(e8_sqfree, 13272,
                         "base-8 squarefree excess at 1e6 must match cross-seq-parity-v1")

    def test_semiprime_1e6_slice(self):
        # Semiprime = Omega(n) == 2 (p*q with p!=q AND p^2 both count), the job's definition.
        e8_semi_all = e8_over(lambda n: self.omega[n] == 2, N)
        self.assertEqual(e8_semi_all, 2449,
                         "base-8 semiprime (Omega==2) excess at 1e6 must match cross-seq-parity-v1")


if __name__ == "__main__":
    unittest.main(verbosity=2)
