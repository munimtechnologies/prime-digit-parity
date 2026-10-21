# Cross-axis additivity & consistency audit (D5 + D6 + D7)

_Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. This repository has
earned **$0** and has **zero external engagement**; it is a curiosity/reputation
project, published in the open. Code is Apache-2.0; result data is CC BY 4.0._

This is **not** a new compute experiment. It is a self-contained re-analysis that
re-derives, offline, three internal-consistency facts from values already published
in this repo across three **independently sealed** jobs:

- **Axis D5** (`v0.29.0`, `data/axis_d5_density_scaling_law.json`) — single-prime
  ratios `r(p)` and the strict squarefree control `r_control`.
- **Axis D6** (`v0.30.0`, `data/axis_d6_density_additivity.json`) — two-prime
  ratios `r({p,q})` and summed density terms.
- **Axis D7** (`v0.31.0`, `data/axis_d7_three_prime_additivity.json`) — three-prime
  ratios `r({p,q,s})` and summed density terms.

Here `r(S) = |z(lag 2)| / median |z(lag k≠2)|` is the frozen lag-2 gap-parity ratio
measured over the squarefree-except-`S` sieve of the Thue–Morse gap-parity stream;
`r_control` is the same ratio over the strict squarefree sieve. All three jobs used
`N = n_max = 20,000,000`. The verifier is
[`test_additivity_excess_audit.py`](../test_additivity_excess_audit.py) (pure stdlib,
no network, no sealed job).

## A. The control is bit-identical across the three axes

```
r_control  =  0.014860313424210063   (D5 == D6 == D7)
```

A deterministic control at the same `N` must reproduce exactly; it does. The three
axes are therefore on a common footing and their ratios are directly comparable.

## B. The density term is the analytic admitted-class density

Each axis weights a relaxed prime by `d(p) = 1/p² − 1/p³`, the density of integers
whose p-adic valuation is exactly 2 (i.e. `p²|n`, `p³∤n`) — the class admitted when
the squarefree cap is relaxed at `p`. Summing this identity over each pair and triple
reproduces the **published** summed-density terms to better than `1e-9`, including the
fully out-of-sample triple `{41,43,47}` (which was never measured singly, but whose
density term is analytic and so still checkable):

| set | Σ d(p), identity | published | match |
|-----|------------------|-----------|-------|
| {17,19} | 0.0058809552… | 0.0058809552431880255 | ✓ |
| {23,29} | 0.0029562282… | 0.0029562281901740674 | ✓ |
| {31,37} | 0.0017177336… | 0.0017177335642310492 | ✓ |
| {17,19,23} | 0.0076891249… | 0.007689124882375992 | ✓ |
| {29,31,37} | 0.0028657921… | 0.0028657921152171494 | ✓ |
| {41,43,47} | 0.0015516918… | 0.0015516917548888875 | ✓ |

So the density inputs to the D6/D7 coefficients are the exact admitted-class density,
not a fitted quantity.

## C. The measured excess is additive — directly, across three jobs

Define the excess `e(S) = r(S) − r_control`. The additivity claim, stated at the level
of the **raw measured ratio** (not routed through the coefficient band `[5,11]`), is:

> `e(S)` measured in D6/D7 equals the sum of the single-prime excesses `e(p)` measured
> separately in D5.

For every pair and in-sample triple this holds to under 1%:

| set | measured e(S) | Σ e(p) from D5 | relative residual |
|-----|---------------|----------------|-------------------|
| {17,19} | 0.0490079702 | 0.0489534412 | 0.11% |
| {23,29} | 0.0233761164 | 0.0232877275 | 0.38% |
| {31,37} | 0.0138401477 | 0.0138789559 | 0.28% |
| {17,19,23} | 0.0637383690 | 0.0632155382 | 0.82% |
| {29,31,37} | 0.0229030911 | 0.0229045864 | 0.0065% |

The largest residual is 0.82% (`{17,19,23}`); the tightest is 0.0065% (`{29,31,37}`).
`{41,43,47}` is excluded because its single primes were never measured — the audit
does not fabricate them.

## Honest limits

- This confirms additivity at the **magnitude** level and ties D5, D6 and D7 together
  as mutually consistent; it is a verification, not a new phenomenon. Because every
  single `k(p) ≈ 8`, raw-excess additivity and the near-constant coefficient are two
  views of the same fact. The value here is that the check is **independent of the
  k-band test**, uses three separately sealed jobs, and would expose any transcription
  or scaling error between the axes. It does not.
- The numbers are the rounded ratios as published; residuals are reported to the
  precision the stored values carry.
- **Whether the coefficient `C ≈ 8` is universal or drifts with prime size remains
  open** — that is Axis D8's out-of-sample sweep over primes 41..97, which is queued
  and has not yet run.

Every number above is re-derivable offline by running the verifier.
