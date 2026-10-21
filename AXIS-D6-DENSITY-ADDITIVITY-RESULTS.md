# Axis D6 — the density scaling law is additive: two-prime relaxations perturb the signature by the sum of the single-prime densities

**Headline (two pre-registered predictions, both held).** Relaxing the strict squarefree (μ²)
exponent cap at **two** primes at once — `squarefree-except-{p,q}` — keeps the base-2 gap-parity
lag-2 unique-quiet signature **present** at all three pairs {17,19}, {23,29}, {31,37} (no pair is
destroyed), and the two-prime super-control excess `r({p,q}) − r_control` equals the **sum** of the
single-prime excesses: the dimensionless coefficient `k({p,q}) = (r({p,q}) − r_control)/(d(p) + d(q))`
with `d(p) = 1/p² − 1/p³` lands inside the pre-registered band `[5, 11]` for all three pairs
(`k = 8.33, 7.91, 8.06`), clustered on the **same** `C ≈ 8` as the single-prime Axis-D5 law. Both the
**convergence** prediction and the new **additivity** prediction are **HELD**.

*This repository's work is produced with AI agents working in MV2's name, at Munim, Inc. It is a
reputation / curiosity number-theory result: **$0 earned, no external engagement.** Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.*

## What was pre-registered, and what happened

The pre-registration (`notes/axis-d6-density-additivity-strict-mu2-prereg.md`, committed before this
sealed job ran) locked six signature verdicts and a quantitative additivity test, with three named
outcomes.

- **Primary (signature verdicts).** `squarefree-except-{17,19}`, `-except-{23,29}`,
  `-except-{31,37}` all **PRESENT** (each pair's summed density keeps the ratio below the 0.1
  threshold); control **present**; odd-squarefree anchor and primes negative control **absent**.
- **Secondary (density additivity — the new quantitative test).** The two-prime coefficient
  `k({p,q}) = (r({p,q}) − r_control)/(d(p) + d(q))`, with `r_control` measured in the *same* job,
  lands inside the committed band `[5, 11]` — the band validated for the single-prime law in Axis D5
  — for each of the three pairs. A pass is positive evidence that two independent prime relaxations
  perturb the signature **additively in the admitted density**.
- Named alternatives, each a real result: **super-additive destruction** (any pair absent though each
  single relaxation is present), or **present-but-not-additive** (all present but some `k({p,q})`
  outside `[5, 11]`).

The clean outcome is what the data show: **convergence confirmed and the density law is additive.**

## Result (N = 20,000,000, sealed python:3.12 container, exit 0, self-tests passed)

Base-2 gap-parity stream `x_n = P_2(g_n) = s_2(g_n) mod 2` (Thue–Morse parity of consecutive gaps),
serial statistic `S^(k) = Σ u_n u_{n+k}` with `u_n = 1 − 2x_n`, scored against the exact
marginal-preserving permutation null as `z(k)`. Detector: signature PRESENT iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `r = |z_x(2)| / median{|z_x(k)| : k≠2} ≤ 0.1`.

| sieve | quiet lag | \|z(2)\| | median \|z(k≠2)\| | ratio `r` | verdict | pre-reg | holds |
|---|---|---|---|---|---|---|---|
| squarefree (control) | 2 | 2.91 | 195.97 | **0.01486** | PRESENT | present | ✓ held |
| **squarefree-except-{17,19}** | 2 | 12.93 | 202.43 | **0.06387** | **PRESENT** | present | ✓ held |
| **squarefree-except-{23,29}** | 2 | 7.62 | 199.24 | **0.03824** | **PRESENT** | present | ✓ held |
| **squarefree-except-{31,37}** | 2 | 5.68 | 197.95 | **0.02870** | **PRESENT** | present | ✓ held |
| odd squarefree (anchor) | 8 | 54.98 | 54.33 | 1.012 | absent | absent | ✓ held |
| primes (negative control) | 7 | 8.09 | 1.26 | 6.444 | absent | absent | ✓ held |

Gap-stream sizes (members ≤ N): squarefree 12,158,575 (max gap 10); except-{17,19} 12,230,413 (max
gap 10); except-{23,29} 12,194,594 (max gap 9); except-{31,37} 12,179,492 (max gap 9); odd squarefree
8,105,700 (max gap 16); primes 1,270,607 (max gap 180). Every pair sieve keeps all even integers and
so produces odd gaps; the stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` is reused verbatim.

## The additivity law (the new, out-of-sample-confirmed result)

Axis D5 established that a **single**-prime relaxation lifts the frozen ratio above the control by
`r(p) − r_control ≈ C·d(p)`, `d(p) = 1/p² − 1/p³`, `C ≈ 8`. If two independent relaxations compose,
relaxing at both `p` and `q` should lift the ratio by `C·(d(p) + d(q))` — the overlap class
`{p²q²·m}` has density `∼ 1/(p²q²) ∼ 10⁻⁵–10⁻⁶`, negligible against `d(p) ∼ 10⁻³`, so the union's
admitted density is `d(p) + d(q)` to first order. Forming the two-prime coefficient
`k({p,q}) = (r({p,q}) − r_control)/(d(p) + d(q))` with `r_control = 0.014860` measured in the same job:

| pair `{p,q}` | {17,19} | {23,29} | {31,37} |
|---|---|---|---|
| `d(p)+d(q)` | 0.0058810 | 0.0029562 | 0.0017177 |
| `r({p,q})` | 0.06387 | 0.03824 | 0.02870 |
| `r({p,q}) − r_control` | 0.04901 | 0.02338 | 0.01384 |
| `k({p,q})` | **8.33** | **7.91** | **8.06** |

All three coefficients lie inside the pre-registered band `[5, 11]`, so `additivity_law_holds = True`.
More tellingly, they sit right on the same `C ≈ 8` measured for the single primes in Axis D5
(`k(p) = 7.86–10.2` for `p = 2…37`, clustering near 8 for the larger primes). The two-prime excess is
the sum of the single-prime excesses: **the density scaling law composes linearly.** This was
committed as a numerical prediction and then confirmed on pairs never used to fit the law.

A back-of-envelope additivity check makes the point directly. Using the Axis-D5 single-prime law
`r(p) − r_control ≈ 8·d(p)`, the predicted pair excess is `8·(d(p)+d(q))`: for {17,19} that is
`8·0.0058810 = 0.04705` vs. the measured `0.04901`; for {23,29} `8·0.0029562 = 0.02365` vs. `0.02338`;
for {31,37} `8·0.0017177 = 0.01374` vs. `0.01384`. The measured two-prime excesses match the summed
single-prime prediction to within a few percent — no sign of destructive or constructive interaction
between the two relaxed primes.

## What this changes in the arc

Axis D3 (v0.27.0) turned "destroyed at every single prime" into a quantitative statement with a
crossover at 11; Axis D4 (v0.28.0) confirmed monotone convergence past it at 17/19/23; Axis D5
(v0.29.0) upgraded the convergence to a *predictive* single-prime density scaling law
`r(p) − r_control ≈ 8·(1/p² − 1/p³)`. Axis D6 shows that law is **additive**: two simultaneous
single-prime relaxations perturb the signature by the sum of their densities, with the same
coefficient. The strict-μ² localization is therefore not only a smooth, density-controlled approach
(D5) but a **linear** one in the admitted density — the perturbations from distinct relaxed primes
superpose. The open mechanism question (*why* strict μ² carries a lag-2 cancellation, and *why* the
coefficient is ≈ 8) is unchanged and now sharper: any mechanism must reproduce a per-prime excess that
is linear and additive in `1/p² − 1/p³`.

## Reproduce

- `data/axis_d6_density_additivity.json` — the sealed job's full output, bit-identical.
- `axis_d6_density_additivity_job.py` — the job (self-tests the sieve/null/stride/detector and the
  cap sets `{17:2,19:2}`, `{23:2,29:2}`, `{31:2,37:2}` and their unions, then computes), runnable as
  `python axis_d6_density_additivity_job.py --n-max 20000000 --out-dir .`.
- `test_axis_d6_density_additivity.py` — pure-stdlib verifier: brute-forces the permutation-null
  formula on tiny sequences, re-derives every `E_null`/`var_null`/`z` from the raw `(n_eff, s, S)`,
  re-runs the frozen detector, asserts the six verdicts and the convergence flags, and **re-derives
  `k({p,q}) = (r({p,q}) − r_control)/(d(p) + d(q))` from the stored ratios and checks each pair's
  coefficient lands in the pre-registered band `[5, 11]`.**

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the *size* of the relaxed prime
(Axis D3 crossover, Axis D4 convergence), establishes the single-prime density scaling law (Axis D5),
or — the new contribution here — shows that law is **additive** under a simultaneous two-prime
relaxation. That additivity result is new.
