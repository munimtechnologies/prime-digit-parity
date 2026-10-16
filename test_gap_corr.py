"""Self-contained checks for the lag-1 gap-parity correlation result (stdlib only).

Validates, independently of the (new) measured numbers:
  * the v0.6.2 method anchors on primes < 1e6 (pi, telescoping, one odd gap,
    half-gap Thue-Morse identity, two-routine digit guard);
  * the S_b adjacent-pair sum on a hand-checked tiny sequence;
  * the closed-form permutation-null MEAN E_perm[S] = (E^2 - N)/N, by brute force
    over EVERY permutation of a small +/-1 multiset;
  * the cancellation-free permutation-null VARIANCE closed form
    Var[S] = (N-1) - 2 s^2 / N + s^4 / (N^2 (N-1)) - 1/(N-1)  (s = sum chi),
    against exact enumeration with fractions.Fraction (the form cited in
    GAP-CORR.md and matching null_moments.perm_var shipped in v0.6.3).

Run: python -m unittest -v test_gap_corr
Work by Munim, Inc., produced with AI. MV2 of Munim, Inc. Apache-2.0.
"""
import itertools
import unittest
from fractions import Fraction

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


def S_of(chi):
    return sum(chi[i] * chi[i + 1] for i in range(len(chi) - 1))


def perm_var_closed(s, N):
    """Cancellation-free exact permutation variance of S for a +/-1 multiset.

    Var[S] = (N-1) - 2 s^2 / N + s^4 / (N^2 (N-1)) - 1/(N-1),  s = sum chi.
    Returns an exact Fraction. (Identical value to null_moments.perm_var.)
    """
    s2 = Fraction(s * s)
    return (Fraction(N - 1)
            - 2 * s2 / N
            + s2 * s2 / (Fraction(N) ** 2 * (N - 1))
            - Fraction(1, N - 1))


class TestGapCorr(unittest.TestCase):
    def test_S_small(self):
        # +1 +1 -1 +1 -> pairs (+1)(+1)=+1, (+1)(-1)=-1, (-1)(+1)=-1 => -1
        self.assertEqual(S_of([1, 1, -1, 1]), -1)
        self.assertEqual(S_of([1, 1, 1, 1]), 3)     # all same
        self.assertEqual(S_of([1, -1, 1, -1]), -3)  # all alternate

    def test_perm_null_mean_bruteforce(self):
        # For every permutation of a fixed +/-1 multiset, the average S equals
        # (E^2 - N)/N with E = sum(chi), N = len.  (E = -excess, but E^2 only.)
        for multiset in ([1, 1, 1, -1, -1], [1, -1, -1, -1], [1, 1, -1, -1, -1, -1]):
            N = len(multiset)
            E = sum(multiset)
            closed = (E * E - N) / N
            perms = list(set(itertools.permutations(multiset)))
            avg = sum(S_of(list(p)) for p in perms) / len(perms)
            self.assertAlmostEqual(avg, closed, places=9,
                                   msg=f"multiset {multiset}")

    def test_perm_null_variance_bruteforce(self):
        # The cancellation-free variance closed form (GAP-CORR.md) matches the
        # exact permutation variance enumerated with Fraction.  Includes the
        # hand-checked anchors Var = 8/9, 1, 0 and an imbalanced N=7 case.
        cases = {
            (1, 1, -1): Fraction(8, 9),       # N=3, s=1
            (1, -1, -1, -1): Fraction(1),     # N=4, s=-2
            (1, 1, 1, 1): Fraction(0),        # N=4, s=4 (one sign -> var 0)
            (1, -1): Fraction(0),             # N=2, s=0 single edge, one layout
        }
        extra = [[1, 1, 1, -1, -1], [1, 1, -1, -1, -1, -1], [1, 1, 1, 1, 1, -1, -1]]
        for multiset in [list(k) for k in cases] + extra:
            N = len(multiset)
            s = sum(multiset)
            perms = list(set(itertools.permutations(multiset)))
            mean = Fraction(sum(S_of(list(p)) for p in perms), len(perms))
            ex2 = Fraction(sum(S_of(list(p)) ** 2 for p in perms), len(perms))
            var_exact = ex2 - mean * mean
            var_closed = perm_var_closed(s, N)
            self.assertEqual(var_closed, var_exact, msg=f"multiset {multiset}")
            if tuple(multiset) in cases:
                self.assertEqual(var_closed, cases[tuple(multiset)],
                                 msg=f"anchor {multiset}")

    def test_method_anchors(self):
        s = sieve(LIMIT)
        primes = [n for n in range(2, LIMIT + 1) if s[n]]
        self.assertEqual(len(primes), PI_1E6)
        gaps = [primes[i + 1] - primes[i] for i in range(len(primes) - 1)]
        self.assertEqual(sum(gaps), primes[-1] - 2)          # telescoping
        self.assertEqual([g for g in gaps if g & 1], [1])    # single odd gap
        excess2 = 0
        half_e2 = 0
        for g in gaps:
            for b in BASES:
                self.assertEqual(parity_div(g, b), parity_str(g, b))
            excess2 += 1 if parity_div(g, 2) else -1
            h = g >> 1 if (g & 1) == 0 else g
            if g & 1 == 0:
                self.assertEqual(parity_div(g, 2), parity_div(h, 2))
            half_e2 += 1 if parity_div(h, 2) else -1
        self.assertEqual(half_e2, excess2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
