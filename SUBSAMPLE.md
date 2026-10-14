# The base-8 near-cancellation, structurally: a bit-subsampling identity

*Work by Munim, Inc., produced with AI. Code Apache-2.0; this text CC BY 4.0.*

This note records a small structural fact that sharpens the open question raised
by the cross-base data in [`RESULTS.md`](RESULTS.md): **why does the base-8 prime
digit-sum parity bias nearly vanish (|z| < 15 up to 10^10, sign-oscillating),
while bases 2, 4, 10 and even the sparser base 16 all show a large, growing
bias?**

## The identity

Let `s_b(n)` be the sum of the base-`b` digits of a non-negative integer `n`,
and write `n` in binary as `n = Σ a_i 2^i` with `a_i ∈ {0,1}`.

> **Proposition.** For every `k ≥ 1` and `b = 2^k`,
> `s_b(n) ≡ Σ_{j ≥ 0} a_{jk}  (mod 2)`.
> That is, the parity of the base-`2^k` digit sum of `n` equals the parity of
> the number of 1-bits of `n` at the positions `0, k, 2k, 3k, …` — a
> density-`1/k` subsample of the binary digits.

**Proof.** Group the bits of `n` into blocks of `k`: the `j`-th base-`2^k` digit
is `d_j = Σ_{i=0}^{k-1} a_{jk+i} 2^i`, so `d_j ≡ a_{jk} (mod 2)` (only the lowest
bit of a digit affects its parity). Summing and reducing mod 2,
`s_b(n) = Σ_j d_j ≡ Σ_j a_{jk} (mod 2)`. ∎

This is elementary and belongs to the folklore of 2‑automatic sequences (for
`k=1` it is exactly the Thue–Morse word `t(n)`). What is worth stating is its
consequence for the cross-base comparison among powers of two:

- **base 2** (`k=1`): uses **all** bits → Thue–Morse parity.
- **base 4** (`k=2`): uses bits `0,2,4,…` → **half** the bits.
- **base 8** (`k=3`): uses bits `0,3,6,…` → **one third** of the bits.
- **base 16** (`k=4`): uses bits `0,4,8,…` → **one quarter** of the bits.

So each of these parity sequences is a Thue–Morse–type parity taken over an
arithmetic-progression subsample of the binary expansion. Comparing the base-2,
-4, -8, -16 prime biases is comparing how these particular subsampled parities
correlate with the primes.

## What this does — and does not — explain

A natural first guess is that a sparser subsample should wash the bias out, so
base 8 (1/3 of bits) would be weaker than base 4 (1/2) and base 16 (1/4) weaker
still. **The data refute the monotone-in-density guess.** From
[`data/cross_base_z_to_1e10.csv`](data/cross_base_z_to_1e10.csv), the signed
standardized bias at `x = 10^10` is:

| base | `k` | fraction of bits used | `z` at 10^10 |
|-----:|----:|----------------------:|-------------:|
| 2    | 1   | 1      | +134.3 |
| 4    | 2   | 1/2    | +104.8 |
| 8    | 3   | 1/3    | **−5.8** |
| 16   | 4   | 1/4    | +24.6 |

Base 16 samples *fewer* bits than base 8 yet shows a bias several times larger
(|z| ≈ 24.6 vs 5.8) and of the opposite sign. **Subsampling density alone does
not explain the base-8 cancellation.** The special feature of `k = 3` is therefore not sparsity
but something about the progression `0, 3, 6, …` of bit positions interacting
with the distribution of primes — which remains open here. (One structural
difference: `k=3` is the smallest `k` for which the sampled positions are *not*
a union of residue classes that includes position 1; we do not claim this is the
cause.)

Stating the puzzle in these terms is the contribution: the base-8 anomaly is a
statement about one Thue–Morse subsequence, not about digit sums in general.

## Reproducing

`subsample.py` implements both sides of the identity; `test_subsample.py`
verifies it exhaustively for every `n < 2^16` and bases `2,4,8,16,32,64`, checks
the four hand-worked cases above, and confirms the `k=1` case is the Thue–Morse
word. Run `python -m unittest test_subsample -v`.

## Prior art

The identity is classical for automatic sequences; Thue–Morse (`k=1`) is the
canonical example, and Mauduit–Rivat (2010) prove digit-sum parity of primes has
density 1/2 with power-saving error for every base, which forces all of these
`z`-values to be `o(√x · x^ε)` asymptotically but says nothing about their finite
size or sign. We found no prior note that casts the base-2/4/8/16 prime-bias
comparison as a comparison of arithmetic-progression bit subsamples, nor any that
observes the base-8 cancellation is not explained by subsampling density. This
sharpens, but does not close, the open question in `RESULTS.md`.
