# Exhaustive verification of the exact lag-k permutation-null moments in
# null_moments_lagk.py.
#
# Self-contained: standard library only, no sieve, no repo imports, no file reads.
# For several small +/-1 multisets and every lag k from 1 to N-1 it brute-forces
# EVERY permutation, computes the exact mean, E[S^2] and variance of
# S^(k) = sum_i x_i x_{i+k} with Fraction arithmetic, and checks them against the
# closed forms (duplicated below so the test stands alone). It also checks that the
# lag-k closed form reduces at k=1 to the independently derived lag-1 form shipped
# in v0.6.3, and that one-sign multisets have variance exactly 0 at every lag.
#
# Apache-2.0. Part of munimtechnologies/prime-digit-parity. MV2, Munim, Inc.

import itertools
import unittest
from fractions import Fraction


# ---- brute force (independent of the closed forms) -------------------------
def _S(x, k):
    return sum(x[i] * x[i + k] for i in range(len(x) - k))


def _brute_moments(multiset, k):
    perms = list(itertools.permutations(multiset))
    n = len(perms)
    s1 = Fraction(0)
    s2 = Fraction(0)
    for p in perms:
        v = _S(p, k)
        s1 += v
        s2 += v * v
    mean = s1 / n
    es2 = s2 / n
    return mean, es2, es2 - mean * mean


# ---- closed forms under test (duplicated so the file is self-contained) ----
def _rho2(s, N):
    return Fraction(s * s - N, N * (N - 1))


def _rho4(s, N):
    num = s ** 4 - 6 * N * s ** 2 + 8 * s ** 2 + 3 * N ** 2 - 6 * N
    return Fraction(num, N * (N - 1) * (N - 2) * (N - 3))


def _counts(N, k):
    M = max(0, N - k)
    if M < 2:
        return M, 0, 0
    n_shared = max(0, N - 2 * k)
    return M, n_shared, M * (M - 1) // 2 - n_shared


def _mean(s, N, k):
    M, _, _ = _counts(N, k)
    return Fraction(0) if M <= 0 else M * _rho2(s, N)


def _es2(s, N, k):
    M, ns, nd = _counts(N, k)
    if M <= 0:
        return Fraction(0)
    total = Fraction(M) + 2 * ns * _rho2(s, N)
    if nd > 0:
        total += 2 * nd * _rho4(s, N)
    return total


def _var(s, N, k):
    m = _mean(s, N, k)
    return _es2(s, N, k) - m * m


# lag-1 closed form shipped in v0.6.3 (null_moments.py), reproduced for the
# k=1 cross-check.
def _es2_lag1_v063(s, N):
    num = s ** 4 - 4 * N * s ** 2 + 4 * s ** 2 + N ** 2 - 2 * N
    return Fraction(N - 1) + Fraction(num, N * (N - 1))


MULTISETS = [
    (1, -1),                       # smallest path (one edge)
    (1, 1, -1),
    (1, -1, -1, -1),
    (1, 1, 1, -1, -1),
    (1, 1, -1, -1, -1, -1),
    (1, 1, 1, 1, -1, -1, -1),      # N = 7, 5040 permutations, lags 1..6
    (1, 1, 1, 1),                  # degenerate: var == 0 at every lag
    (-1, -1, -1, -1, -1),          # degenerate negative, var == 0
]


class TestLagKPermutationMoments(unittest.TestCase):
    def test_closed_form_matches_brute_force_all_lags(self):
        for ms in MULTISETS:
            N = len(ms)
            s = sum(ms)
            for k in range(1, N):
                b_mean, b_es2, b_var = _brute_moments(ms, k)
                self.assertEqual(_mean(s, N, k), b_mean, msg="mean %r k=%d" % (ms, k))
                self.assertEqual(_es2(s, N, k), b_es2, msg="E[S^2] %r k=%d" % (ms, k))
                self.assertEqual(_var(s, N, k), b_var, msg="var %r k=%d" % (ms, k))

    def test_lagk_reduces_to_v063_lag1(self):
        for ms in MULTISETS:
            N = len(ms)
            s = sum(ms)
            self.assertEqual(_es2(s, N, 1), _es2_lag1_v063(s, N),
                             msg="lag-1 mismatch %r" % (ms,))
            self.assertEqual(_mean(s, N, 1), Fraction(s * s - N, N),
                             msg="lag-1 mean mismatch %r" % (ms,))

    def test_degenerate_zero_variance_all_lags(self):
        for ms in [(1, 1, 1, 1), (-1, -1, -1, -1, -1)]:
            N = len(ms)
            s = sum(ms)
            for k in range(1, N):
                self.assertEqual(_var(s, N, k), Fraction(0),
                                 msg="one-sign var nonzero %r k=%d" % (ms, k))

    def test_empty_statistic_is_zero(self):
        # k >= N: no terms, both moments identically zero.
        self.assertEqual(_mean(0, 4, 4), Fraction(0))
        self.assertEqual(_es2(0, 4, 5), Fraction(0))
        self.assertEqual(_var(0, 4, 4), Fraction(0))


if __name__ == "__main__":
    unittest.main()
