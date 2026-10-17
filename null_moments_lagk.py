# Exact permutation-null moments of the linear lag-k sign-correlation statistic.
#
# S^(k) = sum_{i=1..N-k} x_i * x_{i+k}, where (x_1,...,x_N) is a uniform random
# permutation of a fixed multiset of values in {+1, -1}. Writing s = sum x_i
# (so sum x_i^2 = N, since each x_i^2 = 1), the first two permutation moments
# have the closed forms below, for every lag k >= 1.
#
# This generalizes null_moments.py (the lag-1 case, v0.6.3) to arbitrary lag.
# The lag-k "edge set" {(i, i+k) : 1 <= i <= N-k} is a disjoint union of k paths
# (one per residue class of position mod k), so the only structural quantity that
# changes with k is how many edge pairs share a vertex. It is the linear (path,
# non-wrap), +/-1-multiset specialization of the classical serial-correlation
# permutation moments of Wald & Wolfowitz (1943) and Pitman (1937); see
# NULL-MOMENTS-LAGK.md for the derivation and exhaustive verification
# (test_null_moments_lagk.py). No novelty is claimed for the method itself.
#
# Apache-2.0. Part of munimtechnologies/prime-digit-parity. MV2, Munim, Inc.

from fractions import Fraction


def _rho2(s, N):
    """E[x_i x_j] over distinct ordered position pairs. Requires N >= 2."""
    return Fraction(s * s - N, N * (N - 1))


def _rho4(s, N):
    """E[x_i x_j x_k x_l] over distinct ordered 4-tuples. Requires N >= 4.

    From Newton's identities with power sums p1=s, p2=N, p3=s, p4=N the sum over
    distinct ordered 4-tuples is M4 = s^4 - 6 N s^2 + 8 s^2 + 3 N^2 - 6 N.
    """
    num = s ** 4 - 6 * N * s ** 2 + 8 * s ** 2 + 3 * N ** 2 - 6 * N
    return Fraction(num, N * (N - 1) * (N - 2) * (N - 3))


def edge_counts(N, k):
    """(M, n_shared, n_disjoint) for the lag-k path-union on N nodes.

    M        = number of edges (i, i+k)            = max(0, N - k)
    n_shared = unordered edge pairs sharing a vertex = max(0, N - 2k)   (if M >= 2)
    n_disjoint = C(M, 2) - n_shared.
    """
    M = max(0, N - k)
    if M < 2:
        return M, 0, 0
    n_shared = max(0, N - 2 * k)
    n_disjoint = M * (M - 1) // 2 - n_shared
    return M, n_shared, n_disjoint


def perm_mean(s, N, k=1):
    """E[S^(k)] under a uniform permutation of a +/-1 multiset. Exact (Fraction).

    s = sum x_i, N = number of values (>= 2), lag k >= 1. Zero when N - k <= 0.
    Equals (N-k)(E_b^2 - N)/(N(N-1)) with s = -E_b, matching the GAP-CORR / lag-k
    null mean. For k = 1 this is (s^2 - N)/N.
    """
    M, _, _ = edge_counts(N, k)
    if M <= 0:
        return Fraction(0)
    return M * _rho2(s, N)


def perm_e_s2(s, N, k=1):
    """E[(S^(k))^2] under the same null. Exact (Fraction). Requires N >= 2."""
    M, n_shared, n_disjoint = edge_counts(N, k)
    if M <= 0:
        return Fraction(0)
    total = Fraction(M) + 2 * n_shared * _rho2(s, N)
    if n_disjoint > 0:  # needs N >= 4; guarded so rho4 is never formed otherwise
        total += 2 * n_disjoint * _rho4(s, N)
    return total


def perm_var(s, N, k=1):
    """Var[S^(k)] = E[(S^(k))^2] - E[S^(k)]^2. Exact (Fraction).

    Zero for a one-sign multiset (|s| == N) and when N - k <= 0.
    """
    m = perm_mean(s, N, k)
    return perm_e_s2(s, N, k) - m * m
