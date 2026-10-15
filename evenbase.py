"""evenbase.py -- Uniform even-base digit-sum parity reduction.

Work by Munim, Inc., produced with AI. Author: MV2 of Munim, Inc.
Licence: Apache-2.0 (code). Prose and data in EVENBASE.md are CC BY 4.0.

Main identity. For EVERY even base b and every integer n >= 0,

    s_b(n)  ==  sum_{i>=0} floor(n / b^i)   (mod 2),                     (E)

where s_b(n) is the parity of the base-b digit sum. Equivalently, with the
least-significant-bit character eps(x) = (-1)^x,

    chi_b(n) = (-1)^{s_b(n)} = prod_{i>=0} eps(floor(n / b^i)).           (E')

Why it is true. Legendre's classical formula for the base-b digit sum is
    S_b(n) = n - (b-1) * sum_{i>=1} floor(n / b^i).
Reducing mod 2 and using that b is even (so b-1 is odd and drops out) gives
    S_b(n) == n + sum_{i>=1} floor(n/b^i) == sum_{i>=0} floor(n/b^i)  (mod 2),
which is (E). Equivalently and directly: for even b, (x mod b) mod 2 = x mod 2,
so the parity of the i-th base-b digit equals the parity of floor(n/b^i).

What it is for. (E')/(E) express the digit-sum parity of ANY even base -- not
only the powers of two of HIERARCHY.md -- as a product of lsb characters of the
successive b-ary truncations floor(n/b^i). For b = 2^a the i-th truncation's lsb
is exactly bit a*i of n, recovering the dyadic positional lemma. For a
non-dyadic even base b = 2^a * m with m > 1 odd (bases 6 = 2*3, 10 = 2*5,
12 = 4*3), the truncations floor(n/b^i) genuinely mix the odd part m, so these
factors are NOT single bits of n -- the structural difference studied in
EVENBASE.md. The elementary identity (E) is a one-line consequence of Legendre;
its novelty here is the uniform even-base framing and its use to compare bases.
"""


def digit_sum_parity(n: int, b: int) -> int:
    """s_b(n) mod 2, computed directly from the base-b digits (LHS of (E))."""
    if b < 2:
        raise ValueError("base must be >= 2")
    if n < 0:
        raise ValueError("n must be >= 0")
    s = 0
    while n:
        s ^= (n % b) & 1
        n //= b
    return s


def truncation_parity(n: int, b: int) -> int:
    """RHS of (E) for even b: parity of #{ i>=0 : floor(n/b^i) is odd }."""
    if b % 2 != 0:
        raise ValueError("identity (E) holds only for even b")
    if n < 0:
        raise ValueError("n must be >= 0")
    s = 0
    q = n
    while q:
        s ^= q & 1
        q //= b
    return s


def chi(n: int, b: int) -> int:
    """chi_b(n) = (-1)^{s_b(n)} in {+1, -1}."""
    return 1 - 2 * digit_sum_parity(n, b)


def truncation_factors(n: int, b: int):
    """The lsb characters eps(floor(n/b^i)) for i >= 0 as a list in {+1,-1}.

    Their product is chi_b(n) when b is even (identity (E')).
    """
    if b % 2 != 0:
        raise ValueError("identity (E') holds only for even b")
    out = []
    q = n
    while q:
        out.append(1 - 2 * (q & 1))
        q //= b
    return out


def num_factors(n: int, b: int) -> int:
    """Number of i >= 0 with b^i <= n, i.e. floor(log_b n) + 1 for n >= 1.

    For an even base this is the number of non-vanishing lsb factors in (E').
    The i = 0 factor is eps(n), which is -1 for every odd n (so for every
    prime > 2 in every even base).
    """
    if n <= 0:
        return 0
    k = 0
    q = n
    while q:
        k += 1
        q //= b
    return k


if __name__ == "__main__":
    # Quick self-check of (E) and (E') over a small range and several even bases.
    for base in (2, 4, 6, 8, 10, 12, 16):
        for m in range(0, 5000):
            assert digit_sum_parity(m, base) == truncation_parity(m, base)
            prod = 1
            for f in truncation_factors(m, base):
                prod *= f
            if m:
                assert prod == chi(m, base)
    print("evenbase self-check passed for bases 2,4,6,8,10,12,16 up to n=5000")
