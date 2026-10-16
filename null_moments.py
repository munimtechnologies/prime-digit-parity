# Exact permutation-null moments of the linear lag-1 sign-correlation statistic.
#
# S = sum_{i=1..N-1} x_i * x_{i+1}, where (x_1,...,x_N) is a uniform random
# permutation of a fixed multiset of values in {+1, -1}. Writing s = sum x_i
# (so sum x_i^2 = N, since each x_i^2 = 1), the first two permutation moments
# have the closed forms below. This is the linear (path, non-wrap) specialization
# to a +/-1 multiset of the classical serial-correlation permutation moments of
# Wald & Wolfowitz (1943) and Pitman (1937); see NULL-MOMENTS.md for the
# derivation and exhaustive verification (test_null_moments.py).
#
# Apache-2.0. Part of munimtechnologies/prime-digit-parity. MV2, Munim, Inc.

from fractions import Fraction


def perm_mean(s, N):
    """E[S] under a uniform permutation of a +/-1 multiset on a path of N nodes.

    s = sum x_i, N = number of values (>= 2). Exact (Fraction).
    Equals (E_b^2 - N) / N with s = -E_b, matching the GAP-CORR lag-1 null mean.
    """
    return Fraction(s * s - N, N)


def perm_e_s2(s, N):
    """E[S^2] under the same null. Exact (Fraction). Requires N >= 2."""
    num = s ** 4 - 4 * N * s ** 2 + 4 * s ** 2 + N ** 2 - 2 * N
    return Fraction(N - 1) + Fraction(num, N * (N - 1))


def perm_var(s, N):
    """Var[S] = E[S^2] - E[S]^2. Exact (Fraction). Zero for a one-sign multiset."""
    m = perm_mean(s, N)
    return perm_e_s2(s, N) - m * m
