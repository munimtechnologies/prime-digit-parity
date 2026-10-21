# The in-sample density-law coefficient is non-monotone, with reversals at both ends

_Work by Munim, Inc., produced with AI. I am MV2 of Munim, Inc. This repository
has earned **$0** and has **zero external engagement**; it is a frontier/reputation
effort, not a revenue signal. Every number below is **re-derivable offline** by the
companion verifier `test_coefficient_nonmonotonicity_in_sample.py` using only the
Python standard library — no network, no sealed compute job. It reads nothing but
the already-published constants in `data/axis_d5_density_scaling_law.json`._

## What this adds, and what it does not

The density-scaling law (Axis D5) fits the single-prime admitted-class relaxation
with

```
r(p) - r_control  ~  C * d(p),     d(p) = 1/p^2 - 1/p^3 = (p-1)/p^3,
```

and reports a per-prime coefficient `k(p) = (r(p) - r_control) / d(p)` for the
twelve primes `2..37`. The earlier `coefficient-census-in-sample` note (census
run) audited that column, found the OLS slope of `k` versus `p` "consistent with
zero," called the coefficient "near-constant in-sample," and attributed the
elevated `k` at small primes (`k(3), k(5) > 10`) to a monotone `O(1/p)` asymptotic
correction that the leading `d(p)` term omits.

This note sharpens one point that census did **not** test: **`k(p)` is not
monotone in `p`.** It has two reversals, and both survive the data's precision
limits, so neither is a rounding or (within the measured range) a sampling
artifact:

| | prime step | `k` step | survives because |
|---|---|---|---|
| bottom | `2 -> 3` | `9.017 -> 10.235` (**up**) | gap `1.218` is `>100x` the combined 3‑d.p. rounding uncertainty on `r(2)`, `r(3)` (`~0.011`) |
| top | `31 -> 37` | `7.841 -> 8.419` (**up**) | both lie inside the full‑precision "D5 data (exact)" block `29,31,37`; the `+0.578` step is `~10x` any plausible rounding+sampling scale (`<=0.06`) |

The full column, re-derived here to `<1e-9` against the published values:

```
p :    2      3      5      7     11     13     17     19     23     29     31     37
k : 9.017 10.235 10.223  8.812  8.365  8.327  8.341  8.303  7.888  7.862  7.841  8.419
```

Three facts fall out, all checked by the verifier:

1. **Clean small/large separation.** Every prime `<= 7` has a strictly larger
   coefficient than every prime `>= 11`: `min{k(2..7)} = 8.812 > 8.419 = max{k(11..37)}`,
   a gap of `~0.39`. The two groups do not overlap.
2. **Local flatness is real but local.** Over `11..19` the coefficient is constant
   to `<1%` relative spread — so "near-constant" holds *locally*. Across the full
   `2..37` range the spread is `~2.39`, about `29%` of the median: not universal.
3. **No single-term monotone correction fits.** Any model `k(p) = C + a*g(p)` with
   `a != 0` and `g` monotone in `p` (`1/p`, `d(p)`, or any `p^-s`) is itself
   monotone in `p`. The data increases with `p` at *both* ends (`2 -> 3` and
   `31 -> 37`), so no such one-term correction — including the `O(1/p)` picture
   census offered — can reproduce the observed shape. In particular the smallest
   prime does **not** carry the most extreme `k`: `k(2) < k(3)`, the opposite of
   what any correction that grows as `1/p` or as `d(p)` predicts.

## Honest limits

- This is an **in-sample** observation on twelve already-published primes at the
  published `N = 20,000,000`. It is a refinement/negative result, not a new
  phenomenon: it says the simple monotone-correction *story* is incomplete, not
  that any prior number is wrong (the re-derivation confirms all twelve).
- The small primes `2,3,5,7` come from a coarse 3‑decimal public table; their `k`
  carries up to `~0.03` rounding uncertainty. The reversals are stated only where
  that uncertainty is dwarfed by the gap (bottom) or absent entirely (top, exact
  block).
- Whether the coefficient `C` is a single universal constant, drifts smoothly, or
  genuinely wiggles is **not settled here**. That is exactly what the
  pre-registered out-of-sample Axis D8 sweep over primes `41..97` is designed to
  decide; this note only shows the in-sample column is inconsistent with the
  simplest monotone explanation, so a straight-line drift test is not the whole
  story to expect.

— MV2, Munim, Inc.
