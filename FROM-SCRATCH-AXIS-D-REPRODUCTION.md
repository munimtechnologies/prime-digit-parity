# From-scratch independent reproduction of Axis D: relaxing the squarefree exponent cap at one prime kills the lag-2 unique-quiet signature

*Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0. This is a
reputation/curiosity number-theory micro-result: **$0 earned, no external engagement.** It adds
a second, independent code path to a published negative result — it does not claim anything new
about the world, only that the Axis-D conclusion survives an adversarial re-implementation that
trusts none of the original job's data.*

## Why this file exists

The Axis-D headline ([`AXIS-D-EXPONENT-CAP-LAG-SIGNATURE-RESULTS.md`](AXIS-D-EXPONENT-CAP-LAG-SIGNATURE-RESULTS.md),
release v0.24.0) was established by a single sealed compute job at **N = 2·10⁷** and then audited
by [`test_axis_d_exponent_cap_lag_signature.py`](test_axis_d_exponent_cap_lag_signature.py),
which re-derives every number **from that job's committed JSON**. That is a strong internal check,
but it shares the job's sieve code and its stored `(n_eff, s, S)` triples. A reader who distrusts
the sieve, the gap extraction, or the stride bookkeeping is not reassured by re-reading the job's
own output.

[`test_from_scratch_axis_d.py`](test_from_scratch_axis_d.py) closes that gap. With the **Python
standard library only** and **no committed data**, it sieves the three Axis-D streams from
scratch at a modest N, maps gaps to Thue–Morse bits, recomputes the exact permutation-null
z-profile, and reproduces the Axis-D verdicts on a completely separate code path. It is the
Axis-D counterpart of the from-scratch reproductions already shipped for the squarefree control
(v0.19.0) and for Axes B/C (v0.22.0).

## The three sieves, built two independent ways

Axis D compares the **strict squarefree** control against two **single-prime exponent-cap
relaxations**:

- **squarefree-except-2** — admit `v₂(n) ≤ 2` (allow `4 | n`, e.g. 4, 12, 20, 36 …, but not
  `8 | n`); every odd prime still capped at exponent 1.
- **squarefree-except-3** — admit `v₃(n) ≤ 2` (allow `9 | n`, e.g. 9, 18, 45 …, but not
  `27 | n`); every other prime, including 2, still capped at exponent 1.

The test builds each stream with a fast slice-assignment sieve (strike out multiples of `p²` for
every capped prime, and multiples of the relaxed prime's cube), and then **proves the sieve exact
against independent trial-division factorisation** for every `n ≤ 5000`
(`TestSieveExactness`). It also checks the structural facts that make Axis D meaningful: each
relaxation admits a strict **superset** of the squarefree numbers, `4` is admitted by except-2
but not except-3, `9` is admitted by except-3 but not except-2, and `8`/`27` are excluded by
their respective sieves. So the streams are provably the intended congruence classes, not an
artefact of one sieve implementation.

The exact marginal-preserving permutation-null moment formula (used for every z-score) is itself
re-validated against **exhaustive brute-force permutation** of small bit multisets
(`TestNullMomentFormula`), and the base-2 digit-sum parity is checked to equal the popcount
(Thue–Morse) parity for every gap value.

## What is asserted — fingerprints, not magnitudes

The autocorrelation z-scores scale like `√N`, so their **absolute magnitudes** at the test's
`N = 2·10⁶` are about `√10 ≈ 3.2×` smaller than at the sealed job's `N = 2·10⁷`. (A test that
asserted 2·10⁷-calibrated magnitudes at 2·10⁶ is exactly what failed once before — see the v0.22.0
history.) The **published conclusions are qualitative fingerprints and their relative ordering**,
which are N-robust. So the test asserts only:

| Sieve | admits even? | assertion (N-robust) |
|---|---|---|
| squarefree (control) | yes | detector **fires**; lag 2 is uniquely quiet; `ratio₂ ≤ 0.1` |
| squarefree-except-2 | yes | detector does **not** fire; lag 2 **not** anomalously quiet (`ratio₂ > 0.1`) |
| squarefree-except-3 | yes | detector does **not** fire; lag 2 **not** anomalously quiet (`ratio₂ > 0.1`) |

where `ratio₂ = |z(2)| / median{|z(k)| : k ≠ 2}` is the scale-free quantity the published
detector thresholds at `0.1`. The decisive-split test asserts all three verdicts together plus
`ratio₂(sf) < ratio₂(except-2)` and `ratio₂(sf) < ratio₂(except-3)` — the strict control's lag-2
is far quieter, relative to its own profile, than either relaxed sieve's. These ratios are
reported at N = 2·10⁷ as **0.015** (squarefree, present), **1.142** (except-2, absent) and
**0.773** (except-3, absent); the from-scratch check confirms the same ordering holds at an
independent N on independent code.

## What it does and does not add

This adds **confidence**, not a new claim. It shows the Axis-D negative result — *relaxing the
squarefree exponent cap at even a single prime removes the lag-2 unique-quiet signature, even
though every even integer is still admitted* — does not depend on the original job's sieve code,
its stride bookkeeping, or its stored statistics, and reappears from a cold sieve on a second
implementation at a different sample size. Combined with the Axis-C thinning negative, the
reproduced conclusion stands: the signature needs **exactly** the strict squarefree (μ²)
congruence structure — it tolerates neither thinning a congruence class nor enlarging the
exponent cap.

## Honest limits

Same scope as the primary result: a single regime (`N = 2·10⁶` here, `2·10⁷` in the job), base-2
only (the admit-even/exponent question is non-tautological only at bases 2 and 4 — see
[`notes/base2-even-admission-nontautology.md`](notes/base2-even-admission-nontautology.md)), and
two specific small primes relaxed (2 and 3). No asymptotic claim is made. The reproduction asserts
verdicts and scale-free ratios, not absolute `|z|`, precisely because only those are N-robust.

## Reproduce

```
python -m unittest -v test_from_scratch_axis_d
```

Runs offline in a few seconds, no network, no numpy, no committed data.
