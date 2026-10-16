"""Self-contained verification of the COMPLETED null baseline for ALL bases.

Let chi_b(n) = (-1)^{base-b digit sum of n} and S_b(N) = sum_{n=0}^{N-1} chi_b(n), the
all-integer digit-sum-parity summatory function (repo excess E_b = -S_b).

ODD base b (Theorem B -- exact identity):
    S_b(N) = (1 - chi_b(N)) / 2 = [ s_b(N) is odd ]  in {0, 1}
    i.e. the running parity-excess of all integers equals the base-b digit-sum parity of N
    itself -- it toggles between 0 and 1 and never grows.  (So |S_b(N)| <= 1.)

EVEN base b (Theorem A):
    S_b(N) = chi_b(N // b) if (N % b) is odd else 0,  in {-1, 0, 1}.  (So |S_b(N)| <= 1.)

Unified: |S_b(N)| <= 1 for EVERY base b >= 2.

NOTE (honesty): a naive unrolling bound suggests odd-base growth ~ log_b N; that bound is
not tight -- the telescoped terms cancel to {0,1}, as Theorem B proves and this test checks.

Part of the CROSS-SEQUENCE release for munimtechnologies/prime-digit-parity.
Work by Munim, Inc., produced with AI (MV2 of Munim, Inc.). Apache-2.0.

Stdlib only; no repo imports, no file reads. Sign convention: chi_b(n) = (-1)^{s_b(n)}.
"""
import unittest


def chi(n, b):
    """(-1) ** (base-b digit sum of n)."""
    s = 0
    while n:
        s += n % b
        n //= b
    return 1 if s % 2 == 0 else -1


def digit_sum(n, b):
    s = 0
    while n:
        s += n % b
        n //= b
    return s


class TestNullBaselineAllBases(unittest.TestCase):
    M = 200000
    ODD_BASES = (3, 5, 7, 9, 15)
    EVEN_BASES = (2, 4, 6, 8, 10, 12, 16)

    def _running(self, b):
        """S_b(0..M) as a list, built incrementally straight from the definition."""
        S = [0] * (self.M + 1)  # S[0] = S_b(0) = 0 (empty sum)
        run = 0
        for n in range(self.M):
            run += chi(n, b)
            S[n + 1] = run
        return S

    def test_odd_exact_identity(self):
        # Theorem B: S_b(N) = (1 - chi_b(N))/2 = [s_b(N) odd], in {0,1}, for odd b.
        for b in self.ODD_BASES:
            S = self._running(b)
            for N in range(self.M + 1):
                identity = (1 - chi(N, b)) // 2
                self.assertEqual(S[N], identity,
                                 f"odd base {b}: S_b({N})={S[N]} != (1-chi)/2={identity}")
                self.assertEqual(S[N], digit_sum(N, b) % 2,
                                 f"odd base {b}: S_b({N})={S[N]} != s_b(N) mod 2")
                self.assertIn(S[N], (0, 1),
                              f"odd base {b}: S_b({N})={S[N]} left {{0,1}}")

    def test_even_exact_form_and_bound(self):
        # Theorem A: even b -> chi_b(N//b) if (N%b) odd else 0, in {-1,0,1}.
        for b in self.EVEN_BASES:
            S = self._running(b)
            for N in range(self.M + 1):
                predicted = chi(N // b, b) if (N % b) % 2 == 1 else 0
                self.assertEqual(S[N], predicted,
                                 f"even base {b}: S_b({N})={S[N]} != closed form {predicted}")
                self.assertIn(S[N], (-1, 0, 1),
                              f"even base {b}: S_b({N})={S[N]} left {{-1,0,1}}")

    def test_unified_bound_one(self):
        # The headline unified statement: |S_b(N)| <= 1 for every base tested.
        for b in self.ODD_BASES + self.EVEN_BASES:
            run = 0
            self.assertLessEqual(abs(run), 1)
            for n in range(self.M):
                run += chi(n, b)
                self.assertLessEqual(abs(run), 1,
                                     f"base {b}: |S_b({n + 1})|={abs(run)} exceeds 1")


if __name__ == "__main__":
    unittest.main(verbosity=2)
