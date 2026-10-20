# SPDX-License-Identifier: Apache-2.0
#
# From-scratch, zero-external-data reproduction of the AXIS-D (exponent-cap) result of the
# base-2 gap digit-sum-parity "lag-2 unique-quiet" signature
# (AXIS-D-EXPONENT-CAP-LAG-SIGNATURE-RESULTS.md, release v0.24.0):
#
#   The lag-2 unique-quiet base-2 gap-parity signature is PRESENT for the strict squarefree
#   numbers A005117, but ABSENT once the squarefree exponent cap is relaxed at a SINGLE small
#   prime:
#
#     * squarefree-except-2  (allow 4 | n, i.e. v_2(n) <= 2, but not 8 | n; every odd prime
#                             still capped at exponent 1)  -> signature ABSENT
#     * squarefree-except-3  (allow 9 | n, i.e. v_3(n) <= 2, but not 27 | n; every other prime
#                             incl. 2 still capped at exponent 1) -> signature ABSENT
#
#   Both relaxed sieves still ADMIT even integers (unlike the odd-squarefree Axis-C arm, which
#   lifts the prime 2 and is loud-and-uniform), yet the lag-2 cancellation is gone. Together
#   with the Axis-C thinning negative, this shows the signature tolerates NEITHER thinning a
#   congruence class NOR enlarging the exponent cap: it needs EXACTLY the strict squarefree
#   (mu^2) congruence structure.
#
# Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.
#
# Those verdicts were established by a single sealed compute job at N = 2e7 and then audited by
# a test that re-derives the numbers from that job's committed JSON. This file trusts NO
# committed data. With the Python standard library only it sieves the strict squarefree numbers
# and the two single-prime exponent-cap relaxations from scratch at a modest N, maps each
# consecutive gap to its base-2 digit-sum parity (Thue-Morse bit), recomputes the per-lag
# agreement z-scores under an EXACT marginal-preserving permutation null whose closed-form
# moments are checked here against exhaustive brute force, and verifies the three fingerprints
# that underlie the published Axis-D conclusion. Each sieve is itself checked for exactness
# against independent trial-division factorisation.
#
# ROBUSTNESS NOTE (why this test asserts fingerprints, not absolute |z| thresholds):
# the autocorrelation z-scores scale like sqrt(N), so their absolute magnitudes at the test's
# N = 2e6 are smaller than at the sealed job's N = 2e7. The published conclusions are
# QUALITATIVE fingerprints (fires / does-not-fire; lag-2 quiet / not-quiet) and their RELATIVE
# ordering, which are N-robust; hard magnitude constants calibrated at 2e7 are not. So every
# Axis-D assertion below is either a detector verdict (present / absent) or a scale-free ratio
# comparison (az[2] vs the median of the other lags) -- never a bare |z| threshold. The
# data-volume floors and the exactness checks (each sieve vs trial division; null moments vs
# brute force) remain exact.
#
# Run: python -m unittest test_from_scratch_axis_d
#
# Licence: code Apache-2.0; any text/data CC BY 4.0.

import math
import statistics
import unittest
from itertools import permutations


# ---------------------------------------------------------------------------
# Base-2 digit-sum parity (Thue-Morse bit)
# ---------------------------------------------------------------------------
def parity_bits(gaps):
    return [bin(g).count("1") & 1 for g in gaps]


def _bits_to_int(bits):
    return int("".join("1" if b else "0" for b in bits), 2)


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


def _squarefree_flags(N):
    """is_sf[n] == 1 iff n in [1, N] is strictly squarefree (v_p(n) <= 1 for every prime p)."""
    is_sf = bytearray(b"\x01") * (N + 1)
    is_sf[0] = 0
    for p in _small_primes(int(N ** 0.5) + 1):
        p2 = p * p
        if p2 > N:
            break
        is_sf[p2 :: p2] = bytes(len(is_sf[p2 :: p2]))
    return is_sf


def _squarefree_except_flags(N, relaxed):
    """is_ok[n] == 1 iff v_p(n) <= 1 for every prime p != relaxed, and v_relaxed(n) <= 2.

    I.e. allow relaxed^2 | n but not relaxed^3 | n, with every other prime still capped at
    exponent 1. Built by striking out multiples of p^2 for every prime p != relaxed and
    multiples of relaxed^3 for the relaxed prime. `relaxed` must be prime."""
    flags = bytearray(b"\x01") * (N + 1)
    flags[0] = 0
    for p in _small_primes(int(N ** 0.5) + 1):
        if p == relaxed:
            p3 = p ** 3
            if p3 <= N:
                flags[p3 :: p3] = bytes(len(flags[p3 :: p3]))
        else:
            p2 = p * p
            if p2 <= N:
                flags[p2 :: p2] = bytes(len(flags[p2 :: p2]))
    return flags


def _gaps_from_flags(flags, lo, hi, step=1):
    gaps = []
    ap = gaps.append
    prev = -1
    for n in range(lo, hi + 1, step):
        if flags[n]:
            if prev >= 0:
                ap(n - prev)
            prev = n
    return gaps


def squarefree_gaps(N):
    """Gaps between consecutive strictly squarefree integers in [1, N]."""
    return _gaps_from_flags(_squarefree_flags(N), 1, N)


def squarefree_except_gaps(N, relaxed):
    """Gaps between consecutive squarefree-except-`relaxed` integers in [1, N]."""
    return _gaps_from_flags(_squarefree_except_flags(N, relaxed), 1, N)


# ---------------------------------------------------------------------------
# Independent ground truth: classify n by trial-division factorisation.
# ---------------------------------------------------------------------------
def _is_squarefree_trial(n):
    m = n
    d = 2
    while d * d <= m:
        if m % d == 0:
            v = 0
            while m % d == 0:
                m //= d
                v += 1
            if v > 1:
                return 0
        d += 1
    return 1


def _admit_except_trial(n, relaxed):
    """Ground truth: v_p(n) <= 1 for p != relaxed, and v_relaxed(n) <= 2."""
    m = n
    d = 2
    while d * d <= m:
        if m % d == 0:
            v = 0
            while m % d == 0:
                m //= d
                v += 1
            cap = 2 if d == relaxed else 1
            if v > cap:
                return 0
        d += 1
    # any leftover m > 1 is a single prime factor at exponent 1 -> within every cap
    return 1


# ---------------------------------------------------------------------------
# Lag-k agreement statistic + exact permutation-null z-score
# ---------------------------------------------------------------------------
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
    """Exact first two moments of T_k = (number of agreeing lag-k pairs) under a uniformly
    random permutation of the fixed bit multiset (n1 ones among L bits)."""
    n = L
    n0 = n - n1
    m = L - k
    if m <= 1:
        p_eq = (n1 * (n1 - 1) + n0 * (n0 - 1)) / (n * (n - 1))
        return m * p_eq, (m * p_eq * (1 - p_eq) if m == 1 else 0.0)

    A = 0
    for r in range(k):
        c = len(range(r, L, k))
        edges = c - 1
        if edges >= 1:
            A += edges - 1
    total_pairs = m * (m - 1) // 2
    D = total_pairs - A

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


def z_profile(bits, max_lag=12):
    X = _bits_to_int(bits)
    L = len(bits)
    n1 = sum(bits)
    return {k: abs(lag_z(X, L, n1, k)) for k in range(1, max_lag + 1)}


def lag2_detail(bits, max_lag=12):
    """Published detector plus the scale-free lag-2 ratio.
    present iff lag 2 is the uniquely quiet lag AND az[2] <= 0.1 x the median |z| of the
    other lags. ratio2 = az[2] / median(|z| over the other lags) is sample-size-independent."""
    az = z_profile(bits, max_lag)
    quiet = min(az, key=az.get)
    med_excl2 = statistics.median([az[k] for k in az if k != 2])
    ratio2 = az[2] / med_excl2 if med_excl2 > 0 else float("inf")
    present = (quiet == 2) and (ratio2 <= 0.1)
    return quiet, az, present, ratio2


def brute_moments(seq, k):
    ts = []
    for p in permutations(seq):
        ts.append(sum(1 for i in range(len(p) - k) if p[i] == p[i + k]))
    return statistics.mean(ts), statistics.pvariance(ts)


# ===========================================================================
# Tests
# ===========================================================================
class TestSieveExactness(unittest.TestCase):
    """Each slice-assignment sieve must match independent trial-division factorisation,
    so the Axis-D streams are provably the intended congruence classes."""

    LIMIT = 5000

    def test_squarefree_matches_trial(self):
        flags = _squarefree_flags(self.LIMIT)
        for n in range(1, self.LIMIT + 1):
            self.assertEqual(flags[n], _is_squarefree_trial(n), msg=f"squarefree at n={n}")

    def test_except_2_matches_trial(self):
        flags = _squarefree_except_flags(self.LIMIT, 2)
        for n in range(1, self.LIMIT + 1):
            self.assertEqual(flags[n], _admit_except_trial(n, 2),
                             msg=f"squarefree-except-2 at n={n}")

    def test_except_3_matches_trial(self):
        flags = _squarefree_except_flags(self.LIMIT, 3)
        for n in range(1, self.LIMIT + 1):
            self.assertEqual(flags[n], _admit_except_trial(n, 3),
                             msg=f"squarefree-except-3 at n={n}")

    def test_relaxation_is_strictly_wider(self):
        """Sanity: relaxing the cap at a prime admits a strict superset of the squarefree
        numbers, and genuinely admits some non-squarefree ones (4, 9, ...)."""
        sf = _squarefree_flags(self.LIMIT)
        ex2 = _squarefree_except_flags(self.LIMIT, 2)
        ex3 = _squarefree_except_flags(self.LIMIT, 3)
        for n in range(1, self.LIMIT + 1):
            if sf[n]:
                self.assertTrue(ex2[n] and ex3[n], msg=f"superset broken at n={n}")
        self.assertEqual(ex2[4], 1)   # 4 = 2^2 admitted by except-2
        self.assertEqual(ex2[8], 0)   # 8 = 2^3 excluded
        self.assertEqual(ex2[9], 0)   # 9 = 3^2 still excluded by except-2
        self.assertEqual(ex3[9], 1)   # 9 = 3^2 admitted by except-3
        self.assertEqual(ex3[27], 0)  # 27 = 3^3 excluded
        self.assertEqual(ex3[4], 0)   # 4 = 2^2 still excluded by except-3


class TestNullMomentFormula(unittest.TestCase):
    """The closed-form permutation-null moments must match exhaustive brute force."""

    CASES = [
        ([1, 1, 0, 0, 0], 2),
        ([1, 1, 1, 0, 0, 0], 2),
        ([1, 1, 1, 0, 0, 0], 3),
        ([1, 1, 1, 1, 0, 0, 0], 2),
        ([1, 0, 0, 0, 0, 0], 2),
    ]

    def test_moments_match_brute_force(self):
        for seq, k in self.CASES:
            L = len(seq)
            n1 = sum(seq)
            e_formula, var_formula = null_moments(L, n1, k)
            e_brute, var_brute = brute_moments(seq, k)
            self.assertAlmostEqual(e_formula, e_brute, places=9, msg=f"E {seq} k={k}")
            self.assertAlmostEqual(var_formula, var_brute, places=9, msg=f"V {seq} k={k}")

    def test_thue_morse_identity(self):
        """Base-2 digit-sum parity == popcount parity (Thue-Morse) for every gap value."""
        for g in range(1, 4096):
            self.assertEqual(parity_bits([g])[0], bin(g).count("1") % 2, msg=f"g={g}")


class TestAxisDFingerprints(unittest.TestCase):
    """The headline Axis-D reproduction: the strict squarefree control fires, but relaxing
    the exponent cap at a single prime (except-2 or except-3) does not -- recomputed from
    scratch with no committed data.

    Every assertion is a detector verdict (present / absent) or the scale-free lag-2 ratio;
    never a bare |z| threshold calibrated at a different N."""

    N = 2_000_000

    @classmethod
    def setUpClass(cls):
        cls.sf = parity_bits(squarefree_gaps(cls.N))
        cls.ex2 = parity_bits(squarefree_except_gaps(cls.N, 2))
        cls.ex3 = parity_bits(squarefree_except_gaps(cls.N, 3))

    def test_enough_data(self):
        self.assertGreater(len(self.sf), 1_000_000)   # density 6/pi^2 ~ 0.6079
        self.assertGreater(len(self.ex2), 1_000_000)  # strictly wider than squarefree
        self.assertGreater(len(self.ex3), 1_000_000)

    def test_squarefree_fires_control(self):
        """Positive control: the strict squarefree signature fires -- lag 2 is the uniquely
        quiet lag and anomalously quiet (ratio2 <= 0.1)."""
        quiet, az, present, ratio2 = lag2_detail(self.sf)
        self.assertEqual(quiet, 2, msg=f"expected lag 2 quietest; az={az}")
        self.assertTrue(present, msg=f"squarefree signature absent; az={az}")
        self.assertLessEqual(ratio2, 0.1, msg=f"lag-2 should be anomalously quiet; ratio2={ratio2}")

    def test_except_2_does_not_fire(self):
        """Relaxing the cap at the prime 2 (allow 4|n, keep every even integer) destroys the
        signature: the detector does not fire and lag 2 is NOT anomalously quiet."""
        quiet, az, present, ratio2 = lag2_detail(self.ex2)
        self.assertFalse(present, msg=f"squarefree-except-2 must not fire; az={az}")
        self.assertGreater(ratio2, 0.1,
                           msg=f"except-2 lag-2 should not be anomalously quiet; ratio2={ratio2}")

    def test_except_3_does_not_fire(self):
        """Relaxing the cap at the odd prime 3 (allow 9|n) likewise destroys the signature."""
        quiet, az, present, ratio2 = lag2_detail(self.ex3)
        self.assertFalse(present, msg=f"squarefree-except-3 must not fire; az={az}")
        self.assertGreater(ratio2, 0.1,
                           msg=f"except-3 lag-2 should not be anomalously quiet; ratio2={ratio2}")

    def test_decisive_split(self):
        """The whole Axis-D point in one place: among three sieves that ALL admit even
        integers, only the strict squarefree one carries the lag-2 unique-quiet signature;
        enlarging the exponent cap at even a single prime removes it. So the signature needs
        EXACTLY the strict squarefree (mu^2) congruence structure -- it tolerates neither
        thinning (Axis C) nor enlarging (Axis D)."""
        _, _, sf_present, sf_ratio = lag2_detail(self.sf)
        _, _, ex2_present, ex2_ratio = lag2_detail(self.ex2)
        _, _, ex3_present, ex3_ratio = lag2_detail(self.ex3)
        self.assertTrue(sf_present)
        self.assertFalse(ex2_present)
        self.assertFalse(ex3_present)
        # Scale-free: the strict control's lag-2 is far quieter (relative to its own profile)
        # than either relaxed sieve's lag-2.
        self.assertLess(sf_ratio, ex2_ratio, msg=f"sf={sf_ratio}, ex2={ex2_ratio}")
        self.assertLess(sf_ratio, ex3_ratio, msg=f"sf={sf_ratio}, ex3={ex3_ratio}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
