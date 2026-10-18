# SPDX-License-Identifier: Apache-2.0
"""Exact stride-decomposition of base-2 digit-sum parity into base-2^k streams.

Part of the prime-digit-parity project by MV2 of Munim, Inc. Produced with AI.

The gap-parity study observes, per sieve, which even base carries the strongest
serial correlation of the gap digit-sum parity: primes and squarefree numbers
are base-2 dominant, while sums of two squares are base-2 *quiet* and base-4
*loud*. The notes BIT-STRIDE.md and PRIOR-ART.md establish one direction of the
base-2^k structure: for b = 2^k, the digit-sum parity is the XOR of a single
strided subsample of the binary expansion (bits 0, k, 2k, ...).

This module adds the complementary direction, which those notes do not state:
base-2 parity is the exact XOR (equivalently, product in +/-1 form) of the k
strided base-2^k streams obtained by shifting g by r = 0..k-1 before taking the
base-2^k parity. The k shifts partition the bit positions by residue mod k, so
their XOR recovers every bit, i.e. base-2 parity. In particular base-2 parity is
the XOR of exactly two base-4 streams, which is the algebraic handle behind the
"base-4 loud but base-2 quiet" signature of sums of two squares (see
STRIDE-DECOMP.md).

All identities here are elementary and classical; the contribution is stating the
decomposition explicitly and tying it to the per-sieve base signatures.

Notation (g a non-negative integer, bit_i(g) = (g >> i) & 1):
  s_b(g)   = sum of the base-b digits of g
  P_b(g)   = s_b(g) mod 2          (GF(2) / {0,1} form)
  chi_b(g) = (-1) ** s_b(g)        (character / {+1,-1} form)
"""


def digit_sum(g, b):
    """Sum of the base-b digits of a non-negative integer g."""
    if b < 2:
        raise ValueError("base must be >= 2")
    if g < 0:
        raise ValueError("g must be >= 0")
    s = 0
    while g:
        s += g % b
        g //= b
    return s


def P(g, b):
    """Base-b digit-sum parity in {0, 1}."""
    return digit_sum(g, b) & 1


def chi(g, b):
    """Base-b digit-sum character in {+1, -1}."""
    return 1 if (digit_sum(g, b) & 1) == 0 else -1


def P_pow2_stride(g, k):
    """P_{2^k}(g) via the strided-bit identity: XOR of bits 0, k, 2k, ...

    For b = 2^k each base-b digit occupies bits [jk, jk+k-1]; the digit's parity
    is its lowest bit, bit_{jk}(g). Summing digit parities mod 2 therefore keeps
    one bit per digit: s_{2^k}(g) mod 2 = XOR_{j>=0} bit_{jk}(g).
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    x = 0
    pos = 0
    while (g >> pos) > 0:
        x ^= (g >> pos) & 1
        pos += k
    return x


def P2_from_streams(g, k):
    """Reconstruct P_2(g) from the k complementary base-2^k strided streams.

    Identity: P_2(g) = XOR_{r=0}^{k-1} P_{2^k}( floor(g / 2^r) ).

    P_{2^k}(g >> r) is the XOR of the original bits at positions congruent to r
    (mod k). The k residues r = 0..k-1 partition all bit positions, so the XOR of
    the k streams is the XOR of every bit, which is exactly P_2(g).
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    x = 0
    for r in range(k):
        x ^= P(g >> r, 1 << k)
    return x


def base4_streams(g):
    """The two base-4 strided streams whose XOR is base-2 parity.

    Returns (a, b) with a = P_4(g) (even-position bits) and b = P_4(g >> 1)
    (odd-position bits). Then P_2(g) = a XOR b, and chi_2(g) = chi_4(g) *
    chi_4(floor(g / 2)) in the +/-1 character form.
    """
    return P(g, 4), P(g >> 1, 4)
