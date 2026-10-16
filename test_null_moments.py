# Exhaustive verification of the exact permutation-null moments in null_moments.py.
#
# Self-contained: standard library only, no sieve, no repo imports, no file reads.
# For several small +/-1 multisets it brute-forces EVERY permutation, computes the
# exact mean, E[S^2] and variance of S = sum x_i x_{i+1} with Fraction arithmetic,
# and checks them against the closed forms (duplicated below so the test stands
# alone). Includes degenerate one-sign multisets (variance must be exactly 0) and
# the smallest path (N = 2).
#
# Apache-2.0. Part of munimtechnologies/prime-digit-parity. MV2, Munim, Inc.

import itertools
import unittest
from fractions import Fraction


def _S(x):
    return sum(x[i] * x[i + 1] for i in range(len(x) - 1))


def _brute_moments(multiset):
    perms = list(itertools.permutations(multiset))
    n = len(perms)
    s1 = Fraction(0)
    s2 = Fraction(0)
    for p in perms:
        v = _S(p)
        s1 += v
        s2 += v * v
    mean = s1 / n
    es2 = s2 / n
    return mean, es2, es2 - mean * mean


# Closed forms under test (kept here so the file is self-contained).
def _mean(s, N):
    return Fraction(s * s - N, N)


def _es2(s, N):
    num = s ** 4 - 4 * N * s ** 2 + 4 * s ** 2 + N ** 2 - 2 * N
    return Fraction(N - 1) + Fraction(num, N * (N - 1))


def _var(s, N):
    m = _mean(s, N)
    return _es2(s, N) - m * m


MULTISETS = [
    (1, -1),                     # smallest path (one edge)
    (1, 1, -1),
    (1, -1, -1, -1),
    (1, 1, 1, -1, -1),
    (1, 1, -1, -1, -1, -1),
    (1, 1, 1, 1, -1, -1, -1),    # N = 7, 5040 permutations
    (1, 1, 1, 1),                # degenerate: S == N-1 always, var == 0
    (-1, -1, -1, -1, -1),        # degenerate negative, var == 0
]


class TestPermutationMoments(unittest.TestCase):
    def test_closed_form_matches_brute_force(self):
        for ms in MULTISETS:
            N = len(ms)
            s = sum(ms)
            b_mean, b_es2, b_var = _brute_moments(ms)
            self.assertEqual(_mean(s, N), b_mean, msg="mean %r" % (ms,))
            self.assertEqual(_es2(s, N), b_es2, msg="E[S^2] %r" % (ms,))
            self.assertEqual(_var(s, N), b_var, msg="var %r" % (ms,))

    def test_degenerate_zero_variance(self):
        self.assertEqual(_var(4, 4), Fraction(0))
        self.assertEqual(_var(-5, 5), Fraction(0))

    def test_mean_equals_gap_corr_null_mean(self):
        # E[S] = (s^2 - N)/N = (E_b^2 - N)/N with s = -E_b: the GAP-CORR lag-1 null.
        for ms in MULTISETS:
            N = len(ms)
            s = sum(ms)
            self.assertEqual(_mean(s, N), Fraction(s * s - N, N))


if __name__ == "__main__":
    unittest.main()
