# test_grading.py -- Self-contained verification of the singleton grading of
# the prime digit-sum parity correlation.  Redefines everything inline; stdlib
# only; no repo imports or file reads.
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.
#
# Two independent pillars:
#   (A) Combinatorics -- the graded counting law P(b,k) = C(b,k)*A000296(b-k),
#       its row sums (Bell), its leading column (A000296), the vanishing
#       (b-1)-order P(b,b-1)=0, the EGF recurrence, and the Bernhart identity.
#   (B) Data -- the exact moment-cumulant reconstruction of the base-8 (b=3)
#       and base-4 (b=2) block-parity correlations over the primes below 10^6,
#       graded by singletons, gated on the published anchors E2=4764,
#       E4=2902, E8=-510 and on the exact tower factorisation of chi_2.

import unittest
from fractions import Fraction
from math import comb, isqrt


# ----------------------------- combinatorics -------------------------------

def set_partitions(items):
    items = list(items)
    if len(items) == 0:
        yield []
        return
    if len(items) == 1:
        yield [items]
        return
    first = items[0]
    for rest in set_partitions(items[1:]):
        for i in range(len(rest)):
            yield rest[:i] + [[first] + rest[i]] + rest[i + 1:]
        yield [[first]] + rest


def num_singletons(partition):
    return sum(1 for block in partition if len(block) == 1)


def bell_numbers(upto):
    row = [1]
    bells = [1]
    for _ in range(upto):
        nxt = [row[-1]]
        for x in row:
            nxt.append(nxt[-1] + x)
        row = nxt
        bells.append(row[0])
    return bells


def a000296_enum(n):
    return sum(1 for p in set_partitions(range(n)) if num_singletons(p) == 0)


def a000296_closed(n, bells):
    return sum((-1) ** (n - j) * comb(n, j) * bells[j] for j in range(n + 1))


def a000296_recurrence(upto):
    a = [1]
    for n in range(0, upto):
        a.append(sum(comb(n, m) * a[m] for m in range(0, n)))
    return a


def singleton_triangle(bmax):
    tri = []
    for b in range(bmax + 1):
        counts = [0] * (b + 1)
        for p in set_partitions(range(b)):
            counts[num_singletons(p)] += 1
        tri.append(counts)
    return tri


# ------------------------------ digit parities -----------------------------

def popcount(m):
    return bin(m).count("1")


def chi2(m):
    return 1 if popcount(m) % 2 == 0 else -1


def chi_pow2(m, k):
    s = 0
    i = 0
    while m:
        if i % k == 0:
            s ^= (m & 1)
        m >>= 1
        i += 1
    return 1 if s == 0 else -1


def primes_upto(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, limit + 1, i)))
    return [i for i in range(2, limit + 1) if sieve[i]]


# Published anchors at 10^6 (excess E_b = #odd - #even = -sum_p chi_b(p)).
N_PRIMES_1E6 = 78498
E2_1E6 = 4764
E4_1E6 = 2902
E8_1E6 = -510

# A000296, b = 0..9
A000296 = [1, 0, 1, 1, 4, 11, 41, 162, 715, 3425]


class TestSingletonGrading(unittest.TestCase):
    BMAX = 8

    def test_a000296_three_ways(self):
        bells = bell_numbers(9)
        rec = a000296_recurrence(9)
        for n in range(10):
            enum = a000296_enum(n)
            self.assertEqual(enum, A000296[n], f"enum A000296({n})")
            self.assertEqual(a000296_closed(n, bells), A000296[n], f"closed({n})")
            self.assertEqual(rec[n], A000296[n], f"recurrence({n})")

    def test_egf_recurrence(self):
        # a(n+1) = sum_{m=0}^{n-1} C(n,m) a(m), from A'(x) = (e^x - 1) A(x),
        # A(x) = exp(e^x - 1 - x).  Checked against brute-force enumeration.
        for n in range(0, 9):
            rhs = sum(comb(n, m) * a000296_enum(m) for m in range(0, n))
            self.assertEqual(a000296_enum(n + 1), rhs, f"recurrence at n={n}")

    def test_singleton_counting_law(self):
        # P(b,k) = C(b,k) * A000296(b-k) -- the u^k coefficient of the bivariate
        # EGF exp(u x + e^x - 1 - x) at x^b/b!.
        tri = singleton_triangle(self.BMAX)
        bells = bell_numbers(self.BMAX)
        for b in range(self.BMAX + 1):
            row = tri[b]
            for k in range(b + 1):
                self.assertEqual(
                    row[k], comb(b, k) * a000296_enum(b - k),
                    f"P({b},{k})")
            # row sum = Bell(b)  (u = 1, full moment expansion -- MOMENTS.md)
            self.assertEqual(sum(row), bells[b], f"row sum b={b}")
            # leading column = A000296(b)  (u = 0, singleton-free -- GENERAL.md)
            self.assertEqual(row[0], a000296_enum(b), f"leading col b={b}")

    def test_order_b_minus_one_vanishes(self):
        # Structural corollary: no set partition of [b] has exactly b-1
        # singletons, so the delta^{b-1} order of every block-parity
        # correlation vanishes identically.
        tri = singleton_triangle(self.BMAX)
        for b in range(1, self.BMAX + 1):
            self.assertEqual(tri[b][b - 1], 0, f"P({b},{b-1}) must be 0")
            self.assertEqual(comb(b, b - 1) * a000296_enum(1), 0)

    def test_bernhart_identity(self):
        # B_n = A000296(n) + A000296(n+1)  and  B_n = sum_j C(n,j) A000296(j).
        bells = bell_numbers(9)
        for n in range(9):
            self.assertEqual(bells[n], A000296[n] + A000296[n + 1])
            self.assertEqual(
                bells[n], sum(comb(n, j) * A000296[j] for j in range(n + 1)))

    # --------------------------- data pillar -------------------------------

    def test_primes_anchors_and_tower(self):
        primes = primes_upto(10 ** 6)
        self.assertEqual(len(primes), N_PRIMES_1E6)

        s_chi2 = s_chi4 = s_chi8 = 0
        s_prod2 = s_prod3 = 0  # chi_2 reconstructed via base-4 / base-8 towers
        for p in primes:
            c2 = chi2(p)
            a0 = chi_pow2(p, 2); a1 = chi_pow2(p >> 1, 2)
            x0 = chi_pow2(p, 3); x1 = chi_pow2(p >> 1, 3); x2 = chi_pow2(p >> 2, 3)
            # tower factorisation of chi_2 (EVENBASE.md), checked pointwise
            self.assertEqual(a0 * a1, c2)
            self.assertEqual(x0 * x1 * x2, c2)
            s_chi2 += c2
            s_chi4 += a0
            s_chi8 += x0
            s_prod2 += a0 * a1
            s_prod3 += x0 * x1 * x2

        # anchors: E_b = -sum chi_b
        self.assertEqual(-s_chi2, E2_1E6)
        self.assertEqual(-s_chi4, E4_1E6)
        self.assertEqual(-s_chi8, E8_1E6)
        # the two towers reproduce the base-2 excess exactly
        self.assertEqual(-s_prod2, E2_1E6)
        self.assertEqual(-s_prod3, E2_1E6)

    def test_grading_reconstruction_b2_b3(self):
        primes = primes_upto(10 ** 6)
        N = len(primes)

        # integer moment sums for the base-4 (b=2) and base-8 (b=3) towers
        sA = sB = sAB = 0
        sX = sY = sZ = 0
        sXY = sXZ = sYZ = sXYZ = 0
        for p in primes:
            A = chi_pow2(p, 2); B = chi_pow2(p >> 1, 2)
            X = chi_pow2(p, 3); Y = chi_pow2(p >> 1, 3); Z = chi_pow2(p >> 2, 3)
            sA += A; sB += B; sAB += A * B
            sX += X; sY += Y; sZ += Z
            sXY += X * Y; sXZ += X * Z; sYZ += Y * Z; sXYZ += X * Y * Z

        # exact rational moments
        EA, EB, EAB = Fraction(sA, N), Fraction(sB, N), Fraction(sAB, N)
        EX, EY, EZ = Fraction(sX, N), Fraction(sY, N), Fraction(sZ, N)
        EXY, EXZ, EYZ = Fraction(sXY, N), Fraction(sXZ, N), Fraction(sYZ, N)
        EXYZ = Fraction(sXYZ, N)

        # --- b = 2 : E[AB] = kappa2(A,B) + E[A]E[B] -----------------------
        # k=0 (singleton-free, A000296(2)=1): the connected pair cumulant
        kappa2 = EAB - EA * EB
        # k=2 (all singletons, C(2,2)A000296(0)=1): product of marginals
        allsing2 = EA * EB
        # k=1 is empty: P(2,1) = C(2,1)A000296(1) = 0
        self.assertEqual(kappa2 + allsing2, EAB)              # exact identity
        self.assertEqual(EXYZ, EAB)                           # both equal mean chi_2
        self.assertLess(abs(allsing2), abs(EAB))             # top order suppressed

        # --- b = 3 : E[XYZ] graded by singletons --------------------------
        # k=0 (A000296(3)=1): the connected triple cumulant
        kappa3 = (EXYZ - EX * EYZ - EY * EXZ - EZ * EXY
                  + 2 * EX * EY * EZ)
        # k=1 (C(3,1)A000296(2)=3): singleton x connected pair cumulant
        kxy = EXY - EX * EY
        kxz = EXZ - EX * EZ
        kyz = EYZ - EY * EZ
        order1 = EX * kyz + EY * kxz + EZ * kxy
        # k=2 is empty: P(3,2) = C(3,2)A000296(1) = 0
        # k=3 (C(3,3)A000296(0)=1): all three marginals
        allsing3 = EX * EY * EZ

        # exact graded reconstruction of the 3-point moment
        self.assertEqual(kappa3 + order1 + allsing3, EXYZ)
        # the highest-order (all-singleton) term is far below the signal
        self.assertLess(abs(allsing3), abs(EXYZ))
        # and below the leading connected cumulant
        self.assertLess(abs(allsing3), abs(kappa3))
        # ties to the anchor: mean chi_2 = -E2/N
        self.assertEqual(EXYZ, Fraction(-E2_1E6, N))


if __name__ == "__main__":
    unittest.main(verbosity=2)
