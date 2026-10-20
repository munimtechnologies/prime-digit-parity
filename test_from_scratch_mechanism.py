# SPDX-License-Identifier: Apache-2.0
#
# From-scratch, zero-external-data reproduction of the two MECHANISM results of the
# base-2 gap digit-sum-parity "lag-2 unique-quiet" signature:
#
#   Axis B  (LIOUVILLE-LAG-SIGNATURE-RESULTS.md, v0.20.0): the Liouville +1 level set
#           A028260 (integers with Omega(n) even) does NOT carry the signature -- its
#           whole base-2 gap-parity z-profile is FLAT. Finding: the signature needs a
#           multiplicative INDICATOR, not merely a multiplicative origin (lambda is
#           multiplicative, but its +1 level-set indicator is not).
#
#   Axis C  (ODD-SQUAREFREE-LAG-SIGNATURE-RESULTS.md, v0.21.0): the odd squarefree
#           numbers A056911 (prime 2 lifted, every gap even) do NOT carry the signature
#           either -- they are LOUD and UNIFORM at every lag. Finding: the driver is the
#           sieve ADMITTING even integers (odd gaps + half-shifted stride phase), not
#           multiplicativity.
#
# Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.
#
# Both of those results were established by a single sealed compute job at N = 2e7 and
# then audited by a test that re-derives the numbers from that job's committed JSON.
# This file trusts NO committed data. With the Python standard library only it sieves
# the squarefree numbers, the odd squarefree numbers, the Liouville +1 level set and the
# primes from scratch at a modest N, maps each consecutive gap to its base-2 digit-sum
# parity (Thue-Morse bit), recomputes the per-lag agreement z-scores under an EXACT
# marginal-preserving permutation null whose closed-form moments are checked here against
# exhaustive brute force, and verifies the four distinct fingerprints that underlie the
# published mechanism conclusions. The Omega sieve is itself checked against independent
# trial-division factorisation.
#
# ROBUSTNESS NOTE (why this test asserts fingerprints, not absolute |z| thresholds):
# the autocorrelation z-scores scale like sqrt(N), so their absolute magnitudes at the
# test's N = 2e6 are smaller than at the sealed job's N = 2e7. The published conclusions
# are QUALITATIVE fingerprints (fires / loud-uniform / flat / not-quietest) and their
# RELATIVE ordering, which are N-robust; hard magnitude constants calibrated at 2e7 are
# not. So every mechanism assertion below is either a detector verdict (present / which
# lag is quietest) or a relative comparison between two sieves -- never a bare |z|
# threshold. The data-volume floors and the exactness checks (Omega sieve vs trial
# division; null moments vs brute force) remain exact.
#
# Run: python -m unittest test_from_scratch_mechanism
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
    """is_sf[n] == 1 iff n in [1, N] is squarefree."""
    is_sf = bytearray(b"\x01") * (N + 1)
    is_sf[0] = 0
    for p in _small_primes(int(N ** 0.5) + 1):
        p2 = p * p
        if p2 > N:
            break
        is_sf[p2 :: p2] = bytes(len(is_sf[p2 :: p2]))
    return is_sf


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
    """Gaps between consecutive squarefree integers in [1, N]."""
    return _gaps_from_flags(_squarefree_flags(N), 1, N)


def odd_squarefree_gaps(N):
    """Gaps between consecutive ODD squarefree integers in [1, N] (A056911).
    The prime 2 is lifted, so every gap is even."""
    return _gaps_from_flags(_squarefree_flags(N), 1, N, step=2)


def prime_gaps(N):
    """Gaps between consecutive primes in [2, N]."""
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytes(len(sieve[i * i :: i]))
    return _gaps_from_flags(sieve, 2, N)


# ---------------------------------------------------------------------------
# Big Omega (prime factors with multiplicity) via an additive prime-power sieve,
# and the Liouville +1 level set A028260 = {n : Omega(n) even} = {n : lambda(n)=+1}.
# ---------------------------------------------------------------------------
def omega_counts(N):
    """Omega(n) for n in [0, N]. Omega(1)=0; for N up to a few 1e6 the maximum value
    is about 20 (<< 256), so a bytearray is exact."""
    omega = bytearray(N + 1)
    for p in _small_primes(N):
        pk = p
        while pk <= N:
            for m in range(pk, N + 1, pk):
                omega[m] += 1
            pk *= p
    return omega


def liouville_plus_gaps(N):
    """Gaps between consecutive integers n in [1, N] with lambda(n) = (-1)^Omega(n) = +1
    (Omega even), i.e. the A028260 level set. Includes n=1 (Omega(1)=0, even)."""
    omega = omega_counts(N)
    even_flags = bytearray(N + 1)
    for n in range(1, N + 1):
        if omega[n] % 2 == 0:
            even_flags[n] = 1
    return _gaps_from_flags(even_flags, 1, N)


def _omega_trial(n):
    """Independent Omega(n) by trial division (ground truth for the sieve)."""
    c = 0
    d = 2
    while d * d <= n:
        while n % d == 0:
            n //= d
            c += 1
        d += 1
    if n > 1:
        c += 1
    return c


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
    """
    Exact first two moments of T_k = (number of agreeing lag-k pairs) under a uniformly
    random permutation of the fixed bit multiset (n1 ones among L bits). The lag-k edges
    split into k residue chains; consecutive edges of a chain share a vertex (adjacent),
    every other edge pair is vertex-disjoint, and the covariance depends only on which
    case holds, giving a closed form in falling factorials of n1 and n0.
    """
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


def lag2_signature(bits, max_lag=12):
    """Published detector: lag 2 is the uniquely quiet lag and is <= 10% of the median
    |z| of the other lags."""
    az = z_profile(bits, max_lag)
    quiet = min(az, key=az.get)
    med_others = statistics.median([az[k] for k in az if k != quiet])
    present = (quiet == 2) and (az[2] <= 0.1 * med_others)
    return quiet, az, present


def brute_moments(seq, k):
    ts = []
    for p in permutations(seq):
        ts.append(sum(1 for i in range(len(p) - k) if p[i] == p[i + k]))
    return statistics.mean(ts), statistics.pvariance(ts)


# ===========================================================================
# Tests
# ===========================================================================
class TestOmegaSieve(unittest.TestCase):
    """The additive Omega sieve must match independent trial-division factorisation,
    and must reproduce the Liouville lambda sign on textbook values."""

    def test_omega_matches_trial_division(self):
        omega = omega_counts(5000)
        for n in range(1, 5001):
            self.assertEqual(omega[n], _omega_trial(n), msg=f"Omega mismatch at n={n}")

    def test_liouville_sign_on_known_values(self):
        omega = omega_counts(64)
        # lambda(n) = (-1)^Omega(n); +1 <=> Omega even
        self.assertEqual(omega[1] % 2, 0)   # lambda(1)=+1
        self.assertEqual(omega[2] % 2, 1)   # 2         -> lambda=-1
        self.assertEqual(omega[4] % 2, 0)   # 2^2       -> lambda=+1
        self.assertEqual(omega[6] % 2, 0)   # 2*3       -> lambda=+1
        self.assertEqual(omega[8] % 2, 1)   # 2^3       -> lambda=-1
        self.assertEqual(omega[12] % 2, 1)  # 2^2*3     -> lambda=-1
        self.assertEqual(omega[30] % 2, 1)  # 2*3*5     -> lambda=-1
        self.assertEqual(omega[60] % 2, 0)  # 2^2*3*5   -> lambda=+1


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


class TestMechanismFingerprints(unittest.TestCase):
    """The headline reproduction: four sieves, four distinct base-2 gap-parity
    fingerprints, recomputed from scratch with no committed data.

    Every mechanism assertion here is a detector VERDICT (present / which lag is
    quietest) or a RELATIVE comparison between sieves -- never a bare |z| threshold
    calibrated at a different N. The |z| magnitudes shrink like sqrt(N); the
    fingerprints and their ordering do not."""

    N = 2_000_000

    @classmethod
    def setUpClass(cls):
        cls.sf = parity_bits(squarefree_gaps(cls.N))
        cls.osf = parity_bits(odd_squarefree_gaps(cls.N))
        cls.lio = parity_bits(liouville_plus_gaps(cls.N))
        cls.pr = parity_bits(prime_gaps(cls.N))
        cls.sf_az = z_profile(cls.sf)
        cls.osf_az = z_profile(cls.osf)
        cls.lio_az = z_profile(cls.lio)
        cls.pr_az = z_profile(cls.pr)

    def test_enough_data(self):
        self.assertGreater(len(self.sf), 1_000_000)   # density 6/pi^2
        self.assertGreater(len(self.osf), 300_000)    # odd squarefree (~4/pi^2 of odds)
        self.assertGreater(len(self.lio), 400_000)    # density 1/2
        self.assertGreater(len(self.pr), 100_000)

    def test_squarefree_fires_control(self):
        """Positive control: the squarefree signature fires -- lag 2 is the uniquely
        quiet lag (detector verdict). This is the same fingerprint the v0.19.0
        from-scratch test confirms at this N."""
        quiet, az, present = lag2_signature(self.sf)
        self.assertEqual(quiet, 2, msg=f"expected lag 2 quietest; az={az}")
        self.assertTrue(present, msg=f"squarefree signature absent; az={az}")
        # Relative, not absolute: lag 1 is loud vs the quiet lag 2 (same stream).
        self.assertGreater(az[1], az[2], msg=f"lag-1 should exceed lag-2; az={az}")

    def test_odd_squarefree_does_not_fire(self):
        """Axis C: lifting the prime 2 (every gap even) destroys the lag-2 quiet.
        The detector no longer fires -- there is no unique quiet lag at 2."""
        quiet, az, present = lag2_signature(self.osf)
        self.assertFalse(present, msg=f"odd-squarefree must not fire; az={az}")
        # Relative: odd-squarefree lag 2 is loud compared with the squarefree lag 2
        # (same sieve, prime 2 lifted) -- the loud-everywhere fingerprint.
        self.assertGreater(
            az[2], self.sf_az[2],
            msg=f"odd-sf lag-2 should exceed sf lag-2; osf={az}, sf2={self.sf_az[2]}",
        )

    def test_liouville_plus_is_flat(self):
        """Axis B: the Liouville +1 level set has a FLAT profile. Expressed robustly:
        its loudest lag is far below the squarefree control's loudest lag, and the
        detector does not fire."""
        quiet, az, present = lag2_signature(self.lio)
        self.assertFalse(present, msg=f"Liouville +1 must not fire; az={az}")
        self.assertLess(
            max(self.lio_az.values()), max(self.sf_az.values()),
            msg=f"Liouville profile must be far below squarefree; lio={self.lio_az}",
        )

    def test_primes_do_not_fire_control(self):
        """Negative control: for primes, lag 2 is not the uniquely quiet lag."""
        quiet, az, present = lag2_signature(self.pr)
        self.assertFalse(present, msg=f"prime control should not fire; az={az}")

    def test_contrast_is_sharp(self):
        """The decisive contrast in one place: the SAME squarefree sieve has a quiet
        lag-2 when even integers are present (detector fires) and a loud lag-2 once the
        prime 2 is lifted (odd squarefree; detector does not fire) -- isolating that
        'admits even integers' is what the lag-2 cancellation needs. Stated as detector
        verdicts plus a relative lag-2 comparison, both N-robust."""
        _, _, sf_present = lag2_signature(self.sf)
        _, _, osf_present = lag2_signature(self.osf)
        self.assertTrue(sf_present, msg="squarefree (even admitted) should fire")
        self.assertFalse(osf_present, msg="odd squarefree (2 lifted) should not fire")
        self.assertGreater(
            self.osf_az[2], self.sf_az[2],
            msg=f"lag-2 flips quiet->loud; sf2={self.sf_az[2]}, osf2={self.osf_az[2]}",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
