"""Proof-by-enumeration that the Axis-D6/D7 (multi-prime) and D9 (cap-c) density
nulls combine into ONE identity: distinct primes relaxed to INDEPENDENT exponent
caps.

Self-contained, pure-stdlib verifier. No network, no sealed compute job, no repo
data files. It crosses two results already shipped in this repo:

  * the multi-prime density identity (cap 2 only): for distinct primes relaxed to
    cap 2, density(union of newly-admitted classes) = 1 - prod(1 - d_2(p)), and the
    additive approximation sum d_2(p) overcounts the union by exactly the
    inclusion-exclusion remainder (leading pairwise term d_2(p)*d_2(q));
  * the general single-prime cap-c ratio: d_c(p) = 1/p^2 - 1/p^(c+1), with
    d_c(p)/d_2(p) = sum_{j=0}^{c-2} p^-j.

The new statement proved here is the FULL cross-product: several distinct primes
relaxed to SEPARATE caps at once.

Mechanism (read from this repo's code, not assumed):
  * A member n of the strict squarefree (mu^2) sieve satisfies
    v_p(n) <= caps.get(p, 1) for every prime p, with default cap 1.
  * Relaxing the cap at prime p to c (c >= 2) newly admits EXACTLY
    A_p^(c) = {n : 2 <= v_p(n) <= c}, the disjoint union {v_p == 2} u ... u {v_p == c}.
    Its natural density telescopes to
        d_c(p) = density{v_p >= 2} - density{v_p >= c+1} = 1/p^2 - 1/p^(c+1).

Claim (new, this run):
  Relaxing distinct primes p_1,...,p_k to caps c_1,...,c_k are INDEPENDENT events
  under natural density -- A_{p_i}^(c_i) is a condition on residues modulo
  p_i^(c_i+1), and distinct prime-power moduli are coprime (CRT). Hence
    (U)  density(union_i A_{p_i}^(c_i)) = 1 - prod_i (1 - d_{c_i}(p_i))
    (I)  density(intersection_i A_{p_i}^(c_i)) = prod_i d_{c_i}(p_i)      [exact]
    (O)  the additive approximation sum_i d_{c_i}(p_i) OVERCOUNTS the union by
         exactly the inclusion-exclusion remainder
         sum_{|S|>=2} (-1)^(|S|+1) prod_{i in S} d_{c_i}(p_i);
         for two primes this remainder is exactly d_c(p)*d_{c'}(q).
  Setting c_i == 2 for all i recovers this repo's multi-prime (cap-2) identity.

Everything is checked two independent ways:
  (A) EXACT rational arithmetic with fractions.Fraction (N-independent);
  (B) DIRECT ENUMERATION over n = 1..N, cross-checked against the closed forms by
      (i) a tolerance-free integer identity -- the enumerated count of every
          valuation class equals an inclusion-exclusion of floor(N / modulus) terms
          exactly -- and (ii) an a-priori floor-rounding bound on the density (NOT a
          fitted tolerance): a count of a class built from m floor(N/.) terms differs
          from N * density by fewer than m such terms.

License: Apache-2.0. Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
"""

import unittest
from fractions import Fraction
from functools import reduce
from itertools import combinations, product

N = 300000  # enumeration budget (fixed before any residual was inspected)


# ---------------------------------------------------------------------------
# exact building blocks
# ---------------------------------------------------------------------------

def d_cap(p, c):
    """Exact natural density of {n : 2 <= v_p(n) <= c}: d_c(p) = 1/p^2 - 1/p^(c+1)."""
    assert c >= 2
    return Fraction(1, p * p) - Fraction(1, p ** (c + 1))


def union_product_form(caps):
    """1 - prod(1 - d) over a list of (p, c) pairs."""
    prod_ = reduce(lambda a, b: a * b,
                   (1 - d_cap(p, c) for (p, c) in caps), Fraction(1))
    return 1 - prod_


def union_inclusion_exclusion(caps):
    """sum_{nonempty S} (-1)^(|S|+1) prod_{i in S} d_i."""
    ds = [d_cap(p, c) for (p, c) in caps]
    total = Fraction(0)
    k = len(ds)
    for r in range(1, k + 1):
        sign = 1 if (r % 2 == 1) else -1
        for combo in combinations(range(k), r):
            term = reduce(lambda a, b: a * b, (ds[i] for i in combo), Fraction(1))
            total += sign * term
    return total


# ---------------------------------------------------------------------------
# tolerance-free floor identities for the enumerated counts
# ---------------------------------------------------------------------------

def floor_class_count(caps, Nn):
    """Exact count of {n <= Nn : 2 <= v_{p_i}(n) <= c_i for every i} as an
    inclusion-exclusion over floor(Nn / modulus) terms.  Each prime contributes a
    factor (+ p^2  - p^(c+1)); expanding the product and using that distinct prime
    powers combine multiplicatively in the modulus gives the count."""
    total = 0
    # per prime: choices are (exponent 2, sign +1) and (exponent c+1, sign -1)
    choice_sets = [((2, +1), (c + 1, -1)) for (p, c) in caps]
    for choice in product(*choice_sets):
        modulus = 1
        sign = 1
        for (p, c), (e, s) in zip(caps, choice):
            modulus *= p ** e
            sign *= s
        total += sign * (Nn // modulus)
    return total


# ---------------------------------------------------------------------------
# direct enumeration (the independent verification)
# ---------------------------------------------------------------------------

def vp(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def enumerate_counts(caps, Nn):
    """Direct enumeration over n = 1..Nn.  Returns singles, all pairwise
    intersections, the full intersection and the full union of the classes
    A_i = {2 <= v_{p_i} <= c_i}."""
    k = len(caps)
    single = [0] * k
    pair = {(i, j): 0 for i in range(k) for j in range(i + 1, k)}
    inter = 0
    union = 0
    for n in range(1, Nn + 1):
        mem = [2 <= vp(n, p) <= c for (p, c) in caps]
        for i in range(k):
            if mem[i]:
                single[i] += 1
        for i in range(k):
            for j in range(i + 1, k):
                if mem[i] and mem[j]:
                    pair[(i, j)] += 1
        if all(mem):
            inter += 1
        if any(mem):
            union += 1
    return {'single': single, 'pair': pair, 'inter': inter, 'union': union}


# ===========================================================================
# EXACT (rational, N-independent) tests
# ===========================================================================

PRIMES = [2, 3, 5, 7, 11, 13, 17, 37]
CAPS = [2, 3, 4, 5]


class ExactIdentities(unittest.TestCase):

    def test_single_prime_density_telescopes(self):
        for p in PRIMES:
            for c in CAPS:
                disjoint = sum((Fraction(1, p ** k) - Fraction(1, p ** (k + 1))
                                for k in range(2, c + 1)), Fraction(0))
                self.assertEqual(disjoint, d_cap(p, c))
                self.assertEqual(d_cap(p, c),
                                 Fraction(1, p * p) - Fraction(1, p ** (c + 1)))

    def test_bridge_to_cap_c_ratio(self):
        # d_c(p)/d_2(p) must equal the run-154 truncated geometric series.
        for p in PRIMES:
            for c in CAPS:
                ratio = d_cap(p, c) / d_cap(p, 2)
                geom = sum((Fraction(1, p ** j) for j in range(0, c - 1)),
                           Fraction(0))
                self.assertEqual(ratio, geom)

    def test_two_prime_union_and_overcount(self):
        for p in PRIMES:
            for q in PRIMES:
                if q <= p:
                    continue
                for c in CAPS:
                    for cc in CAPS:
                        dp, dq = d_cap(p, c), d_cap(q, cc)
                        union = union_product_form([(p, c), (q, cc)])
                        self.assertEqual(union, dp + dq - dp * dq)
                        additive = dp + dq
                        overcount = additive - union
                        self.assertEqual(overcount, dp * dq)

    def test_product_form_equals_inclusion_exclusion(self):
        families = [
            [(2, 3), (3, 2)],
            [(2, 4), (5, 3)],
            [(2, 3), (3, 2), (5, 2)],
            [(3, 2), (5, 4), (7, 3)],
            [(2, 2), (3, 3), (5, 4), (7, 5)],
        ]
        for caps in families:
            self.assertEqual(union_product_form(caps),
                             union_inclusion_exclusion(caps))

    def test_recovers_cap2_multiprime(self):
        # all caps == 2 must reproduce 1 - prod(1 - (1/p^2 - 1/p^3)).
        caps = [(p, 2) for p in [11, 13, 17, 19, 23]]
        ref_prod = reduce(lambda a, b: a * b,
                          (1 - (Fraction(1, p * p) - Fraction(1, p ** 3))
                           for (p, _) in caps), Fraction(1))
        self.assertEqual(union_product_form(caps), 1 - ref_prod)


# ===========================================================================
# ENUMERATION (direct, with tolerance-free integer identities + a-priori floor)
# ===========================================================================

class EnumerationTwoPrimes(unittest.TestCase):

    CASES = [(2, 3, 3, 2), (2, 4, 5, 3), (3, 3, 7, 2),
             (2, 2, 3, 2), (5, 2, 7, 3)]

    def test_two_prime_counts(self):
        for (p, c, q, cc) in self.CASES:
            counts = enumerate_counts([(p, c), (q, cc)], N)
            sp, sq = counts['single']
            inter, union = counts['inter'], counts['union']

            # (i) tolerance-free integer identities
            self.assertEqual(sp, floor_class_count([(p, c)], N))
            self.assertEqual(sq, floor_class_count([(q, cc)], N))
            self.assertEqual(inter, floor_class_count([(p, c), (q, cc)], N))
            self.assertEqual(union, sp + sq - inter)

            # (ii) a-priori floor-rounding bounds on the density
            dp, dq = d_cap(p, c), d_cap(q, cc)
            self.assertLess(abs(sp / N - float(dp)), 2 / N + 1e-12)   # 2 floor terms
            self.assertLess(abs(sq / N - float(dq)), 2 / N + 1e-12)
            # INDEPENDENCE: intersection density -> product of the two densities
            self.assertLess(abs(inter / N - float(dp * dq)), 4 / N + 1e-12)
            u2 = dp + dq - dp * dq
            self.assertLess(abs(union / N - float(u2)), 8 / N + 1e-12)

            # overcount is exactly the product (restate, exact)
            self.assertEqual((dp + dq) - u2, dp * dq)


class EnumerationThreePrimes(unittest.TestCase):

    def test_three_prime_union_and_independence(self):
        caps = [(2, 3), (3, 2), (5, 2)]
        counts = enumerate_counts(caps, N)
        s = counts['single']
        inter, union, pair = counts['inter'], counts['union'], counts['pair']

        # tolerance-free integer identities for every class
        for i, pc in enumerate(caps):
            self.assertEqual(s[i], floor_class_count([pc], N))
        for (i, j), cnt in pair.items():
            self.assertEqual(cnt, floor_class_count([caps[i], caps[j]], N))
        self.assertEqual(inter, floor_class_count(caps, N))

        # union density -> 1 - prod(1 - d); 3 singles + 3 pairs + 1 triple
        # => fewer than 3*2 + 3*4 + 8 = 26 floor terms of slack
        u3 = union_product_form(caps)
        self.assertLess(abs(union / N - float(u3)), 26 / N + 1e-12)

        # pairwise and triple independence (density -> product)
        ds = [d_cap(p, c) for (p, c) in caps]
        for (i, j), cnt in pair.items():
            self.assertLess(abs(cnt / N - float(ds[i] * ds[j])), 4 / N + 1e-12)
        self.assertLess(abs(inter / N - float(ds[0] * ds[1] * ds[2])),
                        8 / N + 1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
