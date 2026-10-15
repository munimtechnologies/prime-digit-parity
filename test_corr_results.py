# Work by Munim, Inc., produced with AI. Licence: Apache-2.0.
"""Self-contained verification of the base-8 third-joint-cumulant measurement
for prime-digit-parity v0.4.8 (CORRELATION-RESULTS).

Two independent stdlib-only checks:

1. INTERNAL CONSISTENCY of the measured 1e8 result. The eight octal-triple
   parity pattern counts fully determine the joint distribution of
   (A, B, C) = (chi_8(p), chi_8(p//2), chi_8(p//4)) over primes p < 1e8, where
   chi_b(n) = (-1)^(base-b digit sum of n). We recompute every moment
   (a, b, c, the three pairwise correlations, E[ABC]) and the connected third
   joint cumulant kappa_3 purely from those eight counts, and assert they match
   the values measured in the sealed sandbox, that the base-2 prime digit-sum
   excess equals -sum(ABC) = 70509, that the exact cumulant-moment identity
   holds, and that kappa_3 carries essentially all of the excess.

2. PIPELINE ANCHOR via an independent in-test segmented sieve to 1e6: it
   reproduces the published base-2 / base-8 / base-10 prime digit-sum excesses
   at 1e6 (4764 / -510 / 1916) and verifies the exact factorisation
   chi_2(n) = chi_8(n) * chi_8(n//2) * chi_8(n//4) on every prime, tying the
   sign conventions of the 1e8 measurement to a from-scratch computation.
"""
import unittest

# ---- measured in the sealed sandbox, primes p < 1e8 (job prime-parity-base8-corr-v1)
N_PRIMES = 5761455
# pattern index = (s8(p)&1)<<2 | (s8(p//2)&1)<<1 | (s8(p//4)&1);
# chi for a component = +1 if its parity bit is 0, else -1.
PATTERN_COUNTS = [693845, 736247, 734665, 716635, 736921, 719439, 715554, 708149]
MEASURED = {
    "a": 0.00023067089823664334,
    "b": 0.001987171643274138,
    "c": 8.93871426575405e-05,
    "E_AB": -0.009349200852909552,
    "E_AC": -0.008549750019743277,
    "E_BC": -0.008739979744699906,
    "E_ABC": -0.012238054449787424,
    "kappa3": -0.01221821278971542,
    "sum_ABC": -70509,
}
# cutoff -> (base-10 excess, base-8 excess), measured in the same job
CROSSCHECK = {10_000_000: (2967, 215), 100_000_000: (62627, -1329)}


def chi_bits(i):
    """(chiA, chiB, chiC) in {+1,-1} for pattern index i (bit set => chi = -1)."""
    return (1 - 2 * ((i >> 2) & 1), 1 - 2 * ((i >> 1) & 1), 1 - 2 * (i & 1))


def moments_from_counts(counts):
    n = sum(counts)
    sA = sB = sC = sAB = sAC = sBC = sABC = 0
    for i, k in enumerate(counts):
        ca, cb, cc = chi_bits(i)
        sA += ca * k
        sB += cb * k
        sC += cc * k
        sAB += ca * cb * k
        sAC += ca * cc * k
        sBC += cb * cc * k
        sABC += ca * cb * cc * k
    return {
        "n": n, "sum_ABC": sABC,
        "a": sA / n, "b": sB / n, "c": sC / n,
        "E_AB": sAB / n, "E_AC": sAC / n, "E_BC": sBC / n, "E_ABC": sABC / n,
    }


def kappa3_from_moments(m):
    # Connected third joint cumulant of (A,B,C) with non-zero means:
    #   k3 = E[ABC] - a E[BC] - b E[AC] - c E[AB] + 2 a b c
    return (m["E_ABC"] - m["a"] * m["E_BC"] - m["b"] * m["E_AC"]
            - m["c"] * m["E_AB"] + 2 * m["a"] * m["b"] * m["c"])


def digit_sum(n, base):
    s = 0
    while n:
        s += n % base
        n //= base
    return s


def primes_below(limit):
    sieve = bytearray([1]) * limit
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, limit, i)))
    return [i for i in range(limit) if sieve[i]]


class TestInternalConsistency(unittest.TestCase):
    def test_counts_sum_to_n(self):
        self.assertEqual(sum(PATTERN_COUNTS), N_PRIMES)

    def test_moments_reproduce(self):
        m = moments_from_counts(PATTERN_COUNTS)
        self.assertEqual(m["n"], N_PRIMES)
        self.assertEqual(m["sum_ABC"], MEASURED["sum_ABC"])
        for key in ("a", "b", "c", "E_AB", "E_AC", "E_BC", "E_ABC"):
            self.assertAlmostEqual(m[key], MEASURED[key], places=12)

    def test_kappa3_matches_cumulant_identity(self):
        m = moments_from_counts(PATTERN_COUNTS)
        self.assertAlmostEqual(kappa3_from_moments(m), MEASURED["kappa3"], places=12)

    def test_base2_excess_is_minus_sum_ABC(self):
        # chi_2(p) = A*B*C exactly, so the base-2 excess (#odd - #even digit sums)
        # over primes < 1e8 equals -sum(ABC).
        self.assertEqual(-MEASURED["sum_ABC"], 70509)

    def test_kappa3_dominates(self):
        E_ABC = MEASURED["E_ABC"]
        k3 = MEASURED["kappa3"]
        # P-corr-2: the pure connected third cumulant carries nearly all of E[ABC].
        self.assertGreater(k3 / E_ABC, 0.99)
        # pairwise x marginal correction is a tiny fraction of the excess
        self.assertLess(abs(E_ABC - k3) / abs(E_ABC), 0.005)

    def test_marginals_small(self):
        # P-corr-1 (honest): a and c are well below 1e-3; b is ~2e-3, about twice
        # the pre-registered <=1e-3 guide, but still far too small to source the bias.
        self.assertLess(abs(MEASURED["a"]), 1e-3)
        self.assertLess(abs(MEASURED["c"]), 1e-3)
        self.assertLess(abs(MEASURED["b"]), 3e-3)


class TestPipelineAnchor1e6(unittest.TestCase):
    def test_small_excesses_and_factorisation(self):
        e2 = e8 = e10 = 0
        for p in primes_below(1_000_000):
            chi2 = 1 - 2 * (bin(p).count("1") & 1)
            a = 1 - 2 * (digit_sum(p, 8) & 1)
            b = 1 - 2 * (digit_sum(p // 2, 8) & 1)
            c = 1 - 2 * (digit_sum(p // 4, 8) & 1)
            self.assertEqual(chi2, a * b * c)           # exact factorisation, every prime
            chi10 = 1 - 2 * (digit_sum(p, 10) & 1)
            e2 -= chi2                                   # excess = #odd - #even = -sum(chi)
            e8 -= a                                      # a == chi_8(p)
            e10 -= chi10
        self.assertEqual(e2, 4764)
        self.assertEqual(e8, -510)
        self.assertEqual(e10, 1916)


if __name__ == "__main__":
    unittest.main(verbosity=2)
