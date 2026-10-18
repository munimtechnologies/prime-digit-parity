# Stride decomposition: base-2 parity is the XOR of the base-2^k streams

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. This note states
and proves an exact identity and sets up a falsifiable test; it reports no new
empirical result. As of this writing the project has earned $0 and has no stars,
forks or comments. It is a reputation and curiosity bet, not a money signal.*

## Why this note exists

The gap-parity study measures, for each sieve, the serial correlation of the
base-`b` digit-sum parity of consecutive gaps. The central unexplained finding is
that **different sieves favour different bases**: primes and squarefree numbers
are base-2 dominant, while sums of two squares are base-2 *quiet* and base-4
*loud* (see `../squarefree-gap-parity`, `CROSS-SEQUENCE.md`).

`BIT-STRIDE.md` and `PRIOR-ART.md` already establish one direction of the
base-`2^k` structure: for `b = 2^k`, the digit-sum parity is the XOR of a single
strided subsample of the binary expansion,

```
P_{2^k}(g) = XOR_{j>=0} bit_{jk}(g)          (bits 0, k, 2k, 3k, ...)
```

What those notes do **not** state, and what `THEORY.md` is explicitly missing, is
the complementary direction: how base-2 parity is *built back up* from the
base-`2^k` streams. That is what this note adds.

## The decomposition identity

Write `bit_i(g) = (g >> i) & 1`, `P_b(g) = (sum of base-b digits of g) mod 2`,
and `chi_b(g) = (-1)^(sum of base-b digits of g)`.

**Claim.** For every integer `k >= 1` and every non-negative integer `g`,

```
P_2(g) = XOR_{r=0}^{k-1} P_{2^k}( floor(g / 2^r) ).
```

Equivalently, in the +/-1 character form,

```
chi_2(g) = product_{r=0}^{k-1} chi_{2^k}( floor(g / 2^r) ).
```

**Proof.** Shifting right by `r` maps bit position `p` of `floor(g/2^r)` to bit
`p+r` of `g`. Applying the strided identity to `floor(g/2^r)`,

```
P_{2^k}(floor(g/2^r)) = XOR_{j>=0} bit_{jk}(floor(g/2^r)) = XOR_{j>=0} bit_{jk+r}(g),
```

which is the XOR of the bits of `g` at positions congruent to `r (mod k)`. The `k`
residues `r = 0, 1, ..., k-1` partition every bit position exactly once, so the
XOR over `r` is the XOR of all bits of `g`, i.e. `P_2(g)` (Thue-Morse). ∎

### The base-4 corollary (the case that matters)

Taking `k = 2`:

```
P_2(g) = P_4(g) XOR P_4(floor(g/2)),      chi_2(g) = chi_4(g) * chi_4(floor(g/2)).
```

So the base-2 observable is the XOR of exactly two base-4 streams:
`a(g) = P_4(g)` (the even-position bits) and `b(g) = P_4(floor(g/2))` (the
odd-position bits).

## Why this is the right handle for the per-sieve signatures

Base-2 and base-4 are **not independent observables**. Base-2 parity is the XOR
of the base-4 even-bit stream and the base-4 odd-bit stream. Consequently a sieve
can be **base-4 loud yet base-2 quiet** precisely when the serial signals of the
two base-4 streams `a(g_n)` and `b(g_n)` *cancel* once XORed together, even though
one of them (base-4) is individually strong. This is a concrete, mechanistic
reason the sums-of-two-squares gaps could be base-4 loud and base-2 quiet, and it
is testable: decompose each sieve's base-2 lag-1 signal into the two base-4
streams and measure them separately. The pre-registration for that test is in
`notes/stride-cancellation-prereg.md`. **No such measurement is claimed here.**

## Novelty, honestly

The algebra is elementary and classical (Thue-Morse is a 2-automatic sequence;
strided subsequences of automatic sequences are automatic). The strided identity
`P_{2^k}(g) = XOR bit_{jk}(g)` is folklore and already in this repo. The pieces
that are new *to this project* are (1) the explicit complementary-stride
decomposition `P_2 = XOR_r P_{2^k}(g >> r)` as the exact inverse of the strided
view, (2) the base-4 corollary cast as "base-2 is the XOR of two base-4 streams",
and (3) the falsifiable cancellation mechanism it implies for the per-sieve base
signatures. I claim no new mathematics, only a precise framing and a test plan.

## Verification

`stride_decomp.py` implements every function from the definitions.
`test_stride_decomp.py` checks, **exhaustively for all integers `g` in
`[0, 2^14)`**, with an independent digit-sum reference:

- `P_2(g)` equals popcount parity (Thue-Morse);
- `P_{2^k}(g)` equals the strided-bit formula for `k = 1..4`;
- the decomposition `P_2(g) = XOR_{r} P_{2^k}(g >> r)` for `k = 1..4`;
- the base-4 corollary in both `{0,1}` and `{+1,-1}` forms;
- the odd-base collapse `P_b(g) = g mod 2` for `b = 3,5,7,9,11,15`.

Pure Python, no third-party dependencies, runs in seconds.

## Limits

This is an exact identity plus a test plan, not an empirical finding. Whether the
two base-4 streams actually cancel for sums of two squares (and how they combine
for primes and squarefree numbers) is an open, pre-registered question awaiting a
sandbox run; see the pre-registration. Nothing here has been refereed.
