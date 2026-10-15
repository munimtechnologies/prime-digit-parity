# CORRELATION-RESULTS: the base-2 prime digit-sum bias is measured to be 99.8% a connected 3-point correlation

*Work by Munim, Inc., produced with AI. Code Apache-2.0; this text and the
accompanying dataset CC BY 4.0. — MV2 of Munim, Inc.*

## What this adds

`CORRELATION.md` (v0.4.2) proved an exact algebraic identity and then
**pre-registered** three numbered predictions for a queued measurement. This
file reports that measurement, run over every prime below 1e8 in a sealed
offline sandbox, and scores the predictions honestly. The headline:

> At 1e8 the base-2 (Thue–Morse) prime digit-sum excess is **70509**. Of that,
> **70395 (99.84%)** is the pure connected third joint cumulant `kappa_3` of the
> three shifted base-8 digit-sum parities; the marginal- and pairwise-driven
> corrections together account for only **114 (0.16%)**.

This turns the earlier "to leading order" claim into a measured number with a
reproducible, self-contained check.

## The identity being tested

For `chi_b(n) = (-1)^(sum of base-b digits of n)`, the dyadic reduction gives,
for every `n`, the exact factorisation

```
chi_2(n) = chi_8(n) * chi_8(n//2) * chi_8(n//4)
```

(base-8 digit-sum parity picks out bit positions `= 0 mod 3`; the two shifts
pick out `1 mod 3` and `2 mod 3`; together they are the full popcount parity).
Writing `A = chi_8(p)`, `B = chi_8(p//2)`, `C = chi_8(p//4)` over primes `p < x`,
the standardised base-2 excess is `E_2(x) = -sum_p A B C = -pi(x) * E[ABC]`, and
the moment splits into its connected third cumulant plus lower-order pieces:

```
E[ABC] = kappa_3 + a*E[BC] + b*E[AC] + c*E[AB] - 2abc
```

with marginals `a = E[A]`, `b = E[B]`, `c = E[C]`. If `a, b, c` are near zero
(base-8 is nearly balanced along primes) and the pairwise terms stay small, the
bias is carried by `kappa_3` — a genuinely 3-point, connected effect.

## Measured at 1e8 (primes < 1e8, pi = 5,761,455)

| quantity | value |
|---|---|
| base-8 marginals `a, b, c` | 2.307e-4, 1.987e-3, 8.939e-5 |
| pairwise `E[AB], E[AC], E[BC]` | -9.349e-3, -8.550e-3, -8.740e-3 |
| `E[ABC]` | -1.223805e-2 |
| `kappa_3` | -1.221821e-2 |
| `sum_p ABC` | -70509 |
| base-2 excess `-sum ABC` | **70509** (matches the independent anchor E(b=2)@1e8) |

Derived:

- `kappa_3 / E[ABC] = 0.99838` — the connected third cumulant is essentially all
  of the moment.
- Excess attributed to `kappa_3`: `pi * |kappa_3| = 70395` (**99.84%** of 70509).
- Excess from the marginal/pairwise corrections: `70509 - 70395 = 114` (**0.16%**).

The full joint distribution of `(A,B,C)` is published as the eight octal-triple
pattern counts in `data/base8_cumulant.csv`; every number above is recomputed
from those counts in `test_corr_results.py`, so the result is checkable without
rerunning the sieve.

## Scoring the pre-registered predictions

**P-corr-1 (marginals small; `kappa_3` dominates) — CONFIRMED, with one honest
caveat.** `a = 2.3e-4` and `c = 8.9e-5` are comfortably below the pre-registered
`|mean| <= 1e-3`; `b = 1.99e-3` is about twice that guide and is the one number
that missed its stated bound. It does not change the conclusion: the largest
single correction, `b*E[AC]`, is `1.7e-5` against `E[ABC] = -1.22e-2`, i.e. 0.14%.
`kappa_3` carries the same sign and the same magnitude as `-E_2/pi`, as predicted.

**P-corr-2 (`kappa_3 / E[ABC]` near 1) — CONFIRMED.** Measured 0.9984. Had the
pairwise terms carried the signal, this ratio would have been far from 1 and the
bias would have been a 2-point, not a 3-point, effect — the pre-registered
negative outcome. That did not happen.

**P-corr-3 (pipeline cross-check) — CONFIRMED exactly.** The base-10 odd−even
prime excess at 1e8 came out **62627** and the base-8 excess **−1329**, matching
the independent trajectory anchors to the integer. The same cross-check at 1e7
gives **2967** and **215**. The sign conventions are therefore sound.

## Verification (no trust required)

`test_corr_results.py` is stdlib-only and self-contained:

1. **Internal consistency.** From the eight pattern counts alone it recomputes
   `a, b, c`, all three pairwise correlations, `E[ABC]` and `kappa_3`, and
   asserts they match the sandbox values to 12 places; it checks the exact
   cumulant–moment identity, `sum ABC = -70509`, and `kappa_3/E[ABC] > 0.99`
   with the correction below 0.5%.
2. **From-scratch pipeline anchor.** An in-test segmented sieve to 1e6
   reproduces the published base-2 / base-8 / base-10 prime excesses
   (4764 / −510 / 1916) and verifies `chi_2 = chi_8(n)·chi_8(n//2)·chi_8(n//4)`
   on every prime, independently of the 1e8 run.

## Honest limits

- Single range (1e8). The 99.84% share is reported at one cutoff; the 1e7
  cross-check confirms the pipeline, not the share's stability. Whether the
  `kappa_3` fraction drifts with `x` is open (and tied to the even/odd decade
  law in `DECADES.md`).
- `b = 1.99e-3` exceeded its pre-registered `<=1e-3` guide. Stated plainly above;
  it leaves the conclusion intact but the bound itself was not met.
- This localises the anomaly to the connected 3-point correlation; it does not
  explain *why* `kappa_3` is non-zero along primes. `HIERARCHY.md` shows the
  same object is a pure top-order cumulant only for block widths 2 and 3.

## Prior art

Gupta & Sridhar (arXiv:1812.11841) report a single base-10 cutoff; Mauduit–Rivat
(Annals 2010) and arXiv:2606.04677 give asymptotic balance with power-saving
error only. None decomposes the Thue–Morse prime bias into a connected cumulant
of coarser-base digit parities, and none reports this measured 99.8% third-cumulant
share. A fresh search (2026-10-15) found no prior statement of this result.
