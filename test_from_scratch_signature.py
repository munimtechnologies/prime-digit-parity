# SPDX-License-Identifier: Apache-2.0
#
# From-scratch, zero-external-data reproduction of the base-2 gap digit-sum-parity
# "lag-2 unique-quiet" signature that is specific to the squarefree numbers.
#
# Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.
#
# Every other statistical test in this repository re-derives its numbers from a
# committed JSON file that a sealed compute job produced at N = 2e7. This test is
# different: it trusts NO committed data. It sieves the squarefree numbers (and the
# primes, as a negative control) from scratch at a modest N, forms consecutive gaps,
# maps each gap to its base-2 digit-sum parity (the Thue-Morse bit), and recomputes
# the per-lag agreement z-scores under an EXACT marginal-preserving permutation null
# whose closed-form first two moments are verified here against exhaustive brute force.
# It then runs the published "lag-2 unique-quiet" detector and checks that the
# squarefree stream shows the signature while the prime control does not.
#
# Because it recomputes everything, a reader can run `python -m unittest
# test_from_scratch_signature` and watch the signature emerge in a few seconds with
# no downloads and no reliance on any previously committed result.
#
# Licence: code Apache-2.0; any text/data CC BY 4.0.

import math
import statistics
import unittest
from itertools import permutations


# ---------------------------------------------------------------------------
# Thue-Morse / base-b digit-sum parity
# ---------------------------------------------------------------------------
def digit_sum_parity(n, b):
    """Parity (0/1) of the sum of the base-b digits of n."""
    s = 0
    while n:
        s += n % b
        n //= b
    return s & 1


def thue_morse(n):
    """Thue-Morse bit of n via the standard recurrence."""
    t = 0
    while n:
        t ^= n & 1
        n >>= 1
    return t


# ---------------------------------------------------------------------------
# Sieves -> consecutive gaps (pure standard library, slice-assignment sieves)
# ---------------------------------------------------------------------------
def _small_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytes(len(sieve[i * i :: i]))
    return [i for i in range(2, limit + 1) if sieve[i]]


def squarefree_gaps(N):
    """Gaps between consecutive squarefree integers in [1, N]."""
    is_sf = bytearray(b"\x01") * (N + 1)
    is_sf[0] = 0
    for p in _small_primes(int(N ** 0.5) + 1):
        p2 = p * p
        if p2 > N:
            break
        is_sf[p2 :: p2] = bytes(len(is_sf[p2 :: p2]))
    gaps = []
    ap = gaps.append
    prev = -1
    for n in range(1, N + 1):
        if is_sf[n]:
            if prev >= 0:
                ap(n - prev)
            prev = n
    return gaps


def prime_gaps(N):
    """Gaps between consecutive primes in [2, N]."""
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytes(len(sieve[i * i :: i]))
    gaps = []
    ap = gaps.append
    prev = -1
    for n in range(2, N + 1):
        if sieve[n]:
            if prev >= 0:
                ap(n - prev)
            prev = n
    return gaps


# ---------------------------------------------------------------------------
# Lag-k agreement statistic + exact permutation-null z-score
# ---------------------------------------------------------------------------
def parity_bits(gaps):
    return [bin(g).count("1") & 1 for g in gaps]


def _bits_to_int(bits):
    return int("".join("1" if b else "0" for b in bits), 2)


def equal_pairs(X, L, k):
    """Number of lag-k position pairs whose bits agree, and the pair count."""
    m = L - k
    if m <= 0:
        return 0, 0
    mask = (1 << m) - 1
    diff = (X ^ (X >> k)) & mask
    unequal = bin(diff).count("1")
    return m - unequal, m


def _ff(x, j):
    """Falling factorial x*(x-1)*...*(x-j+1)."""
    r = 1
    for t in range(j):
        r *= (x - t)
    return r


def null_moments(L, n1, k):
    """
    Exact first two moments of T_k = (number of agreeing lag-k pairs) under a
    uniformly random permutation of the fixed bit multiset (n1 ones among L bits).

    The pairs {(i, i+k)} split into k chains (by residue mod k). Consecutive
    edges of a chain share a vertex (adjacent); all other edge pairs are
    vertex-disjoint. The covariance only depends on which of these two cases holds.
    """
    n = L
    n0 = n - n1
    m = L - k  # number of lag-k edges
    if m <= 1:
        p_eq = (n1 * (n1 - 1) + n0 * (n0 - 1)) / (n * (n - 1))
        return m * p_eq, (m * p_eq * (1 - p_eq) if m == 1 else 0.0)

    # adjacency count A (edge pairs sharing a vertex), robust for any L
    A = 0
    for r in range(k):
        c = len(range(r, L, k))  # positions in residue class r
        edges = c - 1
        if edges >= 1:
            A += edges - 1
    total_pairs = m * (m - 1) // 2
    D = total_pairs - A  # vertex-disjoint edge pairs

    p_eq = (n1 * (n1 - 1) + n0 * (n0 - 1)) / (n * (n - 1))
    p_adj = (_ff(n1, 3) + _ff(n0, 3)) / _ff(n, 3)
    p_disj = (_ff(n1, 4) + 2 * _ff(n1, 2) * _ff(n0, 2) + _ff(n0, 4)) / _ff(n, 4)

    e_t = m * p_eq
    var = (
        m * p_eq * (1 - p_eq)
        + 2 * A * (p_adj - p_eq * p_eq)
        + 2 * D * (p_disj - p_eq * p_eq)
    )
    return e_t, var


def lag_z(X, L, n1, k):
    t_obs, _ = equal_pairs(X, L, k)
    e_t, var = null_moments(L, n1, k)
    if var <= 0:
        return 0.0
    return (t_obs - e_t) / math.sqrt(var)


def lag2_signature(bits, max_lag=12):
    """Published detector: lag 2 is the uniquely quiet lag and is <= 10% of the
    median |z| of the other lags."""
    X = _bits_to_int(bits)
    L = len(bits)
    n1 = sum(bits)
    az = {k: abs(lag_z(X, L, n1, k)) for k in range(1, max_lag + 1)}
    quiet = min(az, key=az.get)
    med_others = statistics.median([az[k] for k in az if k != quiet])
    present = (quiet == 2) and (az[2] <= 0.1 * med_others)
    return quiet, az, present


# Brute-force permutation null for tiny cases (ground truth for the formula)
def brute_moments(seq, k):
    ts = []
    for p in permutations(seq):
        ts.append(sum(1 for i in range(len(p) - k) if p[i] == p[i + k]))
    return statistics.mean(ts), statistics.pvariance(ts)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestThueMorseIdentity(unittest.TestCase):
    def test_base2_digit_sum_parity_is_thue_morse(self):
        for n in range(0, 20001):
            self.assertEqual(digit_sum_parity(n, 2), bin(n).count("1") & 1)
            self.assertEqual(digit_sum_parity(n, 2), thue_morse(n))

    def test_thue_morse_recurrence(self):
        for n in range(0, 10001):
            self.assertEqual(thue_morse(2 * n), thue_morse(n))
            self.assertEqual(thue_morse(2 * n + 1), 1 - thue_morse(n))


class TestNullMomentFormula(unittest.TestCase):
    """The closed-form permutation-null moments must match exhaustive brute force."""

    CASES = [
        ([1, 1, 0, 0, 0], 2),
        ([1, 1, 1, 0, 0, 0], 2),
        ([1, 1, 1, 0, 0, 0], 3),
        ([1, 1, 1, 1, 0, 0, 0], 2),
        ([1, 1, 1, 1, 0, 0, 0], 3),
        ([1, 0, 0, 0, 0, 0], 2),
    ]

    def test_moments_match_brute_force(self):
        for seq, k in self.CASES:
            L = len(seq)
            n1 = sum(seq)
            e_formula, var_formula = null_moments(L, n1, k)
            e_brute, var_brute = brute_moments(seq, k)
            self.assertAlmostEqual(
                e_formula, e_brute, places=9, msg=f"E mismatch seq={seq} k={k}"
            )
            self.assertAlmostEqual(
                var_formula, var_brute, places=9, msg=f"Var mismatch seq={seq} k={k}"
            )


class TestEqualPairsCounter(unittest.TestCase):
    def test_counts_match_direct(self):
        for bits in ([1, 0, 1, 1, 0, 0, 1], [0, 0, 0, 1, 1, 1], [1, 1, 1, 1]):
            X = _bits_to_int(bits)
            L = len(bits)
            for k in range(1, L):
                t, m = equal_pairs(X, L, k)
                direct = sum(1 for i in range(L - k) if bits[i] == bits[i + k])
                self.assertEqual(t, direct, msg=f"bits={bits} k={k}")
                self.assertEqual(m, L - k)


class TestFromScratchSignature(unittest.TestCase):
    """The headline reproduction: no committed data, everything recomputed."""

    N = 2_000_000

    @classmethod
    def setUpClass(cls):
        cls.sf_bits = parity_bits(squarefree_gaps(cls.N))
        cls.pr_bits = parity_bits(prime_gaps(cls.N))

    def test_enough_data(self):
        # squarefree density 6/pi^2 ~ 0.6079; sanity-check the sieve produced data
        self.assertGreater(len(self.sf_bits), 1_000_000)
        self.assertGreater(len(self.pr_bits), 100_000)

    def test_squarefree_shows_lag2_signature(self):
        quiet, az, present = lag2_signature(self.sf_bits)
        self.assertEqual(
            quiet, 2, msg=f"expected lag 2 quietest; got {quiet}; az={az}"
        )
        self.assertTrue(
            present,
            msg=f"squarefree lag-2 signature absent at N={self.N}; az={az}",
        )

    def test_primes_do_not_show_lag2_signature(self):
        quiet, az, present = lag2_signature(self.pr_bits)
        self.assertFalse(
            present,
            msg=f"prime control unexpectedly shows the signature; az={az}",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
