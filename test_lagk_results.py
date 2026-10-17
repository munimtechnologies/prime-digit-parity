"""Self-contained checks for the lag-k gap-parity correlation result (stdlib only).

Validates, independently of the (new) measured numbers:
  * the v0.6.2/v0.7.0 method anchors on primes < 1e6 (pi, telescoping, one odd
    gap, half-gap Thue-Morse identity, two-routine digit guard);
  * the lag-k adjacent-pair sum S^(k) on hand-checked tiny sequences;
  * the generalized closed-form permutation-null mean
        E_perm[S^(k)] = (N - k) * (E^2 - N) / (N * (N - 1)),
    by brute force over EVERY permutation of small +/-1 multisets, for k = 1,2,3;
  * that this general form reduces to the v0.7.0 lag-1 form (E^2 - N) / N at k = 1.

No file reads, no repo imports, no network: everything here is recomputed from
scratch, so a pass certifies the method the sealed sandbox job used.

Run: python -m unittest -v test_lagk_results
Work by Munim, Inc., produced with AI. MV2 of Munim, Inc. Apache-2.0.
"""
import itertools
import unittest

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


def S_lag(chi, k):
    return sum(chi[i] * chi[i + k] for i in range(len(chi) - k))


class TestGapLagK(unittest.TestCase):
    def test_S_lag_small(self):
        # +1 +1 -1 +1
        # lag 1 pairs: (+1)(+1), (+1)(-1), (-1)(+1) => 1 -1 -1 = -1
        self.assertEqual(S_lag([1, 1, -1, 1], 1), -1)
        # lag 2 pairs: chi0*chi2, chi1*chi3 => (+1)(-1)+(+1)(+1) = -1+1 = 0
        self.assertEqual(S_lag([1, 1, -1, 1], 2), 0)
        # lag 3: chi0*chi3 => (+1)(+1) = 1
        self.assertEqual(S_lag([1, 1, -1, 1], 3), 1)
        self.assertEqual(S_lag([1, 1, 1, 1], 1), 3)      # all same, lag1
        self.assertEqual(S_lag([1, 1, 1, 1], 2), 2)      # all same, lag2

    def test_perm_null_mean_bruteforce(self):
        # For every permutation of a fixed +/-1 multiset, the average S^(k) equals
        # (N-k)(E^2-N)/(N(N-1)) with E = sum(chi), N = len.
        multisets = (
            [1, 1, -1],
            [1, 1, 1, -1, -1],
            [1, -1, -1, -1],
            [1, 1, -1, -1, -1, -1],
        )
        for multiset in multisets:
            N = len(multiset)
            E = sum(multiset)
            perms = list(set(itertools.permutations(multiset)))
            for k in (1, 2, 3):
                if N - k <= 0:
                    continue
                closed = (N - k) * (E * E - N) / (N * (N - 1))
                avg = sum(S_lag(list(p), k) for p in perms) / len(perms)
                self.assertAlmostEqual(avg, closed, places=9,
                                       msg=f"multiset {multiset} lag {k}")

    def test_lag1_matches_v070_form(self):
        # Sanity: the general form at k=1 must equal the v0.7.0 lag-1 form (E^2-N)/N.
        for multiset in ([1, 1, 1, -1, -1], [1, -1, -1, -1, -1]):
            N = len(multiset)
            E = sum(multiset)
            general = (N - 1) * (E * E - N) / (N * (N - 1))
            v070 = (E * E - N) / N
            self.assertAlmostEqual(general, v070, places=12)

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
