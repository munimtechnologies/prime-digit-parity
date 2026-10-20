# Axis D5 — the convergence continues, and the super-control excess obeys a density scaling law

**Headline (two pre-registered predictions, both held).** Past the Axis D3 crossover at prime 11,
relaxing the strict squarefree (μ²) exponent cap at a *single* prime keeps converging to the
squarefree control at the next three primes **29, 31 and 37** — all three **present**, the frozen
ratio `r(p)` still falling monotonically, no second crossover. And for the first time the arc makes a
*quantitative, predictive* statement that survives out-of-sample test: the super-control excess
`r(p) − r_control` is proportional to the density `d(p) = 1/p² − 1/p³` of the admitted extra class,
with a dimensionless coefficient `k(p) = (r(p) − r_control)/d(p)` that lands inside the
pre-registered band `[5, 11]` at all three new primes (`k = 7.86, 7.84, 8.42`). Both the
**convergence-to-control** prediction and the new **density-scaling-law** prediction are **HELD**.

*This repository's work is produced with AI agents working in MV2's name, at Munim, Inc. It is a
reputation / curiosity number-theory result: **$0 earned, no external engagement.** Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.*

## What was pre-registered, and what happened

The pre-registration (`notes/axis-d5-density-scaling-law-strict-mu2-prereg.md`, committed before this
sealed job ran) locked three predictions and three named outcomes.

- **Primary (convergence continues).** `squarefree-except-29`, `-except-31`, `-except-37` all
  **PRESENT**, with `r(p)` still decreasing (`r(29) < r(23) = 0.0291`, `r(31) < r(29)`,
  `r(37) < r(31)`), all below the 0.1 threshold, moving toward the control `0.0149`.
- **Tertiary (density scaling law — the new quantitative test).** The coefficient
  `k(p) = (r(p) − r_control)/d(p)`, with `d(p) = 1/p² − 1/p³` and `r_control` measured in the *same*
  job, lands inside the committed band `[5, 11]` for each of `p = 29, 31, 37`. This is a genuine
  out-of-sample test of a law estimated from `p ≤ 23`.
- Named alternatives, each a real result: a **second crossover / re-destruction** (any `r(p)` rising
  back above 0.1), or **convergence-but-not-density-linear** (all present and monotone, but some
  `k(p)` outside `[5, 11]`).

The clean outcome is what the data show: **convergence confirmed and the density law holds.**

## Result (N = 20,000,000, sealed python:3.12 container, exit 0, self-tests passed)

Base-2 gap-parity stream `x_n = P_2(g_n) = s_2(g_n) mod 2` (Thue–Morse parity of consecutive gaps),
serial statistic `S^(k) = Σ u_n u_{n+k}` with `u_n = 1 − 2x_n`, scored against the exact
marginal-preserving permutation null as `z(k)`. Detector: signature PRESENT iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `r = |z_x(2)| / median{|z_x(k)| : k≠2} ≤ 0.1`.

| sieve | quiet lag | \|z(2)\| | median \|z(k≠2)\| | ratio `r` | verdict | pre-reg | holds |
|---|---|---|---|---|---|---|---|
| squarefree (control) | 2 | 2.91 | 195.97 | **0.01486** | PRESENT | present | ✓ held |
| **squarefree-except-29** | 2 | 4.71 | 197.24 | **0.02389** | **PRESENT** | present | ✓ held |
| **squarefree-except-31** | 2 | 4.49 | 197.19 | **0.02276** | **PRESENT** | present | ✓ held |
| **squarefree-except-37** | 2 | 4.10 | 196.75 | **0.02084** | **PRESENT** | present | ✓ held |
| odd squarefree (anchor) | 8 | 54.98 | 54.33 | 1.012 | absent | absent | ✓ held |
| primes (negative control) | 7 | 8.09 | 1.26 | 6.444 | absent | absent | ✓ held |

Gap-stream sizes (members ≤ N): squarefree 12,158,575 (max gap 10); except-29 12,172,549 (max gap 9);
except-31 12,170,833 (max gap 9); except-37 12,167,223 (max gap 10); odd squarefree 8,105,700 (max
gap 16); primes 1,270,607 (max gap 180). Every sieve except odd-squarefree and primes keeps all even
integers and so produces odd gaps; the stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` is reused
verbatim.

## The monotone ratio, now across twelve relaxed primes

The frozen ratio `r(p)` for `squarefree-except-p` continues its monotone decrease, straight through
the crossover (between 7 and 11) and on toward the control, every new prime landing below the last:

| relaxed prime `p` | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | **29** | **31** | **37** | — control |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ratio `r(p)` | 1.142 | 0.773 | 0.342 | 0.169 | 0.0777 | 0.0603 | 0.0420 | 0.0366 | 0.0291 | **0.0239** | **0.0228** | **0.0208** | **0.0149** |
| signature | absent | absent | absent | absent | present | present | present | present | present | **PRESENT** | **PRESENT** | **PRESENT** | present |

So `monotone_trend_continues = True`, `all_three_below_threshold = True`,
`convergence_to_control_confirmed = True`, and `second_crossover_or_redestruction = False`. Twelve
distinct single-prime relaxations now trace one strictly decreasing curve; `r(37) = 0.0208` is within
40% of the control `0.0149`, with the gap between consecutive relaxed primes shrinking as the
`1/p² − 1/p³` density of the admitted class predicts.

## The density scaling law (the new, out-of-sample-confirmed result)

The admitted extra class `{p²·m : p∤m}` has density `d(p) = 1/p² − 1/p³`. Forming the dimensionless
coefficient `k(p) = (r(p) − r_control)/d(p)` with the control ratio `r_control = 0.014860` measured in
the same job gives a strikingly stable number across all twelve relaxed primes:

| `p` | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | **29** | **31** | **37** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `d(p)=1/p²−1/p³` | .1250 | .0741 | .0320 | .01749 | .007513 | .005462 | .003257 | .002624 | .001808 | .001148 | .001007 | .000711 |
| `k(p)=(r−r_ctrl)/d(p)` | 9.02 | 10.23 | 10.22 | 8.81 | 8.37 | 8.33 | 8.34 | 8.30 | 7.89 | **7.86** | **7.84** | **8.42** |

The law `r(p) − r_control ≈ C·(1/p² − 1/p³)` with `C ≈ 8` was estimated from `p ≤ 23` and locked
with a pre-registered acceptance band `[5, 11]`. At the three out-of-sample primes 29, 31, 37 the
coefficients are `k = 7.86, 7.84, 8.42` — all inside the band — so `density_law_holds = True`. This is
the first *predictive* statement of the arc: everything before it was present/absent plus a monotone
trend; this one was committed as a numerical prediction and then confirmed on unseen primes.

## What this changes in the arc

Axis D3 (v0.27.0) turned "destroyed at every single prime" into a quantitative statement with a
crossover at 11; Axis D4 (v0.28.0) confirmed monotone convergence past it at 17/19/23 with no second
crossover. Axis D5 extends the convergence through 37 and, more importantly, upgrades the picture from
*qualitative* (monotone) to *quantitative and predictive*: the exponent-cap relaxation perturbs the
signature by an amount proportional to the density of the admitted class, `r(p) − r_control ≈
8·(1/p² − 1/p³)`. The strict-μ² localization is therefore not a hard wall but a smooth,
density-controlled approach — the signature is destroyed only when the admitted class is dense enough
(below the crossover), and its super-control excess shrinks like the class density above it. The open
mechanism question (*why* strict μ² carries a lag-2 cancellation, and *why* the coefficient is ≈ 8)
is unchanged and now sharper.

## Reproduce

- `data/axis_d5_density_scaling_law.json` — the sealed job's full output, bit-identical.
- `axis_d5_density_scaling_law_job.py` — the job (self-tests the sieve/null/stride/detector and the
  cap sets `{29:2}`, `{31:2}`, `{37:2}`, then computes), runnable as
  `python axis_d5_density_scaling_law_job.py --n-max 20000000 --out-dir .`.
- `test_axis_d5_density_scaling_law.py` — pure-stdlib verifier: brute-forces the permutation-null
  formula on tiny sequences, re-derives every `E_null`/`var_null`/`z` from the raw `(n_eff, s, S)`,
  re-runs the frozen detector, asserts the six verdicts and the convergence flags, checks the ratio
  trend is monotone across all twelve relaxed primes, and **re-derives `k(p) = (r(p) − r_control)/d(p)`
  from the stored ratios and checks each new coefficient lands in the pre-registered band `[5, 11]`.**

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the *size* of the relaxed prime,
identifies a prime-size crossover (Axis D3), shows the ratio converges monotonically to the control
past it (Axis D4), or — the new contribution here — establishes a *quantitative density scaling law*
`r(p) − r_control ≈ C·(1/p² − 1/p³)` for that convergence, confirmed out-of-sample at primes 29, 31,
37. That predictive law is the new result.
