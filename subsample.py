# Bit-subsampling identity for digit-sum parity in power-of-two bases.
#
# Work by Munim, Inc., produced with AI. Apache-2.0.
#
# For a base b = 2**k, the parity of the base-b digit sum of n depends ONLY on a
# density-1/k subsample of the binary digits of n: the bits at positions
# 0, k, 2k, 3k, ...  This recasts the cross-base comparison among powers of two
# (studied in RESULTS.md) as a comparison of Thue-Morse-type subsequences.
#
# See SUBSAMPLE.md for the statement, the one-line proof, and why it does NOT
# explain the empirical base-8 near-cancellation on its own.

def digit_sum(n, b):
    """Sum of the base-b digits of the non-negative integer n."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if b < 2:
        raise ValueError("b must be >= 2")
    s = 0
    while n:
        s += n % b
        n //= b
    return s


def digit_sum_parity(n, b):
    """Parity (0/1) of the base-b digit sum of n."""
    return digit_sum(n, b) & 1


def subsample_bit_parity(n, k):
    """Parity (0/1) of the number of 1-bits of n at positions 0, k, 2k, ...

    This is the right-hand side of the identity proved in SUBSAMPLE.md:
    for b = 2**k, digit_sum_parity(n, b) == subsample_bit_parity(n, k).
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    c = 0
    pos = 0
    while (n >> pos) != 0:
        c += (n >> pos) & 1
        pos += k
    return c & 1


def identity_holds(n, k):
    """True iff the identity holds for this n and k (b = 2**k)."""
    return digit_sum_parity(n, 1 << k) == subsample_bit_parity(n, k)


if __name__ == "__main__":
    # Small demonstration across the power-of-two bases studied in RESULTS.md.
    for k, b in [(1, 2), (2, 4), (3, 8), (4, 16)]:
        ok = all(identity_holds(n, k) for n in range(1 << 16))
        frac = 1.0 / k
        print(f"base {b:2d} (k={k}): identity holds for all n < 2^16 = {ok}; "
              f"samples fraction {frac:.3f} of the binary digits")
