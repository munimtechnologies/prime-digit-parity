# The base-2 prime digit-sum bias is a third cumulant of base-8 parities

*Work by Munim, Inc., produced with AI. Signed: MV2 of Munim, Inc.*
*Text and data: CC BY 4.0. Code: Apache-2.0.*

This note sharpens the open question left by `REDUCE.md` into a **quantitative,
measurable** statement. It adds no new empirical number by itself; the numbers
it frames are produced by an independent stdlib compute job (see "Measurement"
below) and will be recorded in `data/base8_correlation.csv` when that job
returns. What is proved here is an exact algebraic decomposition, verified for
all `n < 2^16` by `test_corr.py`.

## Setup (recalled from REDUCE.md)

Let `s_b(n)` be the base-`b` digit sum and write the ±1 parity character
`chi_b(n) = (-1)^{s_b(n)}`. The folklore subsample identity, specialised in
`REDUCE.md`, gives the exact pointwise factorisation

```
chi_2(n) = chi_8(n) * chi_8(n // 2) * chi_8(n // 4)      for all n >= 0.
```

Summing over primes `p < x` and using `#{odd} - #{even} = -sum chi`, the
standardised base-2 excess that Gupta & Sridhar (arXiv:1812.11841) found
"peculiar" (our verified `z ~ +134` to 1e10) is **exactly**

```
E_2(x) = - sum_{p < x} chi_8(p) * chi_8(p//2) * chi_8(p//4).
```

Each base-8 factor is near-balanced along primes (our verified base-8
`|z| < 6` to 1e10, sign-oscillating), yet their product is strongly biased.
`REDUCE.md` observed this is *logically compatible* but did not say which order
of correlation carries the bias. That is what this note pins down.

## The decomposition (new)

Treat the primes `p < x` as a finite sample and write expectations as sample
means over that set. Put the three ±1 random variables

```
A = chi_8(p),   B = chi_8(p//2),   C = chi_8(p//4).
```

Their raw moments are the marginals `a = E[A], b = E[B], c = E[C]`, the three
pairwise means `E[AB], E[AC], E[BC]`, and the triple mean `E[ABC]`. Note
`E[ABC] = -E_2(x)/pi(x)` by the identity above. The standard third **joint
cumulant** (connected 3-point correlation) is

```
kappa_3(A,B,C) = E[ABC] - a*E[BC] - b*E[AC] - c*E[AB] + 2*a*b*c.
```

An elementary expansion (proved and machine-checked in `test_corr.py`) shows
this equals the centred third moment exactly:

```
kappa_3 = E[(A-a)(B-b)(C-c)].
```

Rearranging gives the decomposition of the triple correlation that drives the
base-2 bias:

```
E[ABC] = kappa_3  +  a*E[BC] + b*E[AC] + c*E[AB]  -  2*a*b*c.
```

**Consequence.** Along primes the marginals `a, b, c` are each near zero
(base-8 near-balance), so every term on the right except `kappa_3` is a
*product that contains at least one small marginal factor*. Therefore, to
leading order in the small marginals,

```
E[ABC]  ≈  kappa_3(A,B,C),       i.e.   E_2(x) ≈ - pi(x) * kappa_3.
```

The large, famous base-2 excess is, to leading order, the **genuine connected
third-order correlation** of the three shifted base-8 digit-sum parities along
the primes — not an artefact of any single biased digit statistic, and not
reducible to the pairwise correlations alone. This is a strictly sharper
statement than REDUCE's "compatibility": it names *which* cumulant carries the
anomaly and makes the claim falsifiable by measurement.

## The full-integer baseline (new, exact)

Over a complete octal block `n in [0, 8^m) = [0, 2^{3m})` the identity gives,
**exactly**,

```
sum chi_8(n) chi_8(n//2) chi_8(n//4) = sum chi_2(n) = 0,
sum chi_8(n) = 0,
```

because the Thue-Morse sequence is balanced on every dyadic block `[0, 2^k)`,
`k >= 1`, and `[0, 2^{3m})` is also a complete octal block. (The shifted
marginals `sum chi_8(n//2)`, `sum chi_8(n//4)` and the pairwise sums are equal
to sums over *incomplete* octal blocks and so vanish only in density, not
exactly.) The key fact is that the triple correlation `E[ABC]` — the very
quantity that equals the base-2 bias — is **exactly zero over the full
integers**. Its large nonzero value over the primes is therefore entirely a
*primes* effect, not an artefact of digit combinatorics. (Both exact zeros are
verified in `test_corr.py` on `[0, 2^15) = [0, 8^5)`.)

## Measurement (pre-registered)

`prime-parity-base8-corr-v1` (stdlib segmented sieve, 2 cpu / 2 GB, network
off) computes over primes `p < 1e8`:

1. the eight pattern counts of `(s_8(p), s_8(p//2), s_8(p//4))`;
2. from them `a, b, c`, the three pairwise means, `E[ABC]` and `kappa_3`;
3. as an **independent cross-check**, the base-10 and base-8 odd−even excess at
   value cutoffs `1e7` and `1e8`, with a *different* sieve than
   `prime-digit-parity-v3`.

Pre-registered predictions (to be scored honestly, negative results kept):

- **P-corr-1:** `a, b, c` are each small (`|mean| ≲ 1e-3`), consistent with the
  base-8 near-balance; `kappa_3` is of the same sign and order as
  `-E_2(1e8)/pi(1e8)` (i.e. the connected term dominates).
- **P-corr-2:** the pairwise×marginal terms contribute only a small fraction of
  `E[ABC]` (`kappa_3 / E[ABC]` close to 1). If instead the pairwise terms carry
  most of the signal, P-corr-1 is false and the bias is a 2-point, not 3-point,
  effect — a clean negative result that would still refine the picture.
- **P-corr-3 (cross-check):** base-10 excess at `1e8` ≈ 62627 and base-8 excess
  ≈ −1329 (matching `prime-parity-trajectory-v1`). Disagreement means a pipeline
  bug to fix before trusting the correlation numbers.

## Honest limits and prior art

- The cumulant decomposition is elementary algebra; the novelty is its
  *application* — identifying the Gupta–Sridhar / Thue–Morse prime bias with a
  third joint cumulant of coarser (base-8) digit parities, and the leading-order
  claim that this connected term, not the marginals or pairwise terms, carries
  the anomaly. A prior-art pass (Mauduit–Rivat 2010 Annals; Gupta–Sridhar
  arXiv:1812.11841; the 2026 papers arXiv:2606.04677 and arXiv:2607.10654; OEIS
  A000069/A130911/A380192/A111309) found no statement of this decomposition.
- The leading-order claim `E[ABC] ≈ kappa_3` is quantitative only once the
  marginals and pairwise terms are measured; the job above provides those
  numbers. Until it returns, the "≈" is a prediction, not a result.
- This note does **not** explain *why* `kappa_3` is large — i.e. why primality
  induces a strong 3-way dependence among `s_8(n), s_8(n//2), s_8(n//4)` while
  leaving the marginals balanced. That mechanism remains the open problem; this
  note only localises it to the third cumulant and makes it measurable.
