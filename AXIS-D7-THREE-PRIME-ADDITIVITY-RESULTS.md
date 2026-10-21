# Axis D7 — the density scaling law stays additive for three primes: a three-prime relaxation perturbs the signature by the sum of the three single-prime densities

**Headline (two pre-registered predictions, both held).** Relaxing the strict squarefree (μ²)
exponent cap at **three** primes at once — `squarefree-except-{p,q,s}` — keeps the base-2 gap-parity
lag-2 unique-quiet signature **present** at all three triples {17,19,23}, {29,31,37}, {41,43,47} (no
triple is destroyed), and the three-prime super-control excess `r({p,q,s}) − r_control` equals the
**sum** of the three single-prime excesses: the dimensionless coefficient
`k({p,q,s}) = (r({p,q,s}) − r_control)/(d(p) + d(q) + d(s))` with `d(p) = 1/p² − 1/p³` lands inside
the pre-registered band `[5, 11]` for all three triples (`k = 8.29, 7.99, 8.31`), clustered on the
**same** `C ≈ 8` as the single-prime Axis-D5 law and the two-prime Axis-D6 law. The triple {41,43,47}
was never measured singly or in a pair, so it is a fully out-of-sample test of both convergence and
additivity. Both the **convergence** prediction and the **three-prime additivity** prediction are
**HELD**.

*This repository's work is produced with AI agents working in MV2's name, at Munim, Inc. It is a
reputation / curiosity number-theory result: **$0 earned, no external engagement.** Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.*

## What was pre-registered, and what happened

The pre-registration (`notes/axis-d7-three-prime-additivity-strict-mu2-prereg.md`, committed before
this sealed job ran) locked six signature verdicts and a quantitative three-prime additivity test,
with three named outcomes.

- **Primary (signature verdicts).** `squarefree-except-{17,19,23}`, `-except-{29,31,37}`,
  `-except-{41,43,47}` all **PRESENT** (each triple's summed density keeps the ratio below the 0.1
  threshold); control **present**; odd-squarefree anchor and primes negative control **absent**.
- **Secondary (three-prime additivity — the new quantitative test).** The three-prime coefficient
  `k({p,q,s}) = (r({p,q,s}) − r_control)/(d(p) + d(q) + d(s))`, with `r_control` measured in the *same*
  job, lands inside the committed band `[5, 11]` — the band validated for the single-prime law in
  Axis D5 and the two-prime law in Axis D6 — for each of the three triples. A pass is positive
  evidence that three independent prime relaxations perturb the signature **additively in the admitted
  density**, with no detectable third-order interaction term.
- Named alternatives, each a real result: **super-additive destruction** (any triple absent though
  each single and double relaxation at these primes is present), or **present-but-not-additive** (all
  present but some `k({p,q,s})` outside `[5, 11]`, i.e. a pairwise or triple cross term appears).

The clean outcome is what the data show: **convergence confirmed and the density law is additive to
third order.**

## Result (N = 20,000,000, sealed python:3.12 container, exit 0, self-tests passed)

Base-2 gap-parity stream `x_n = P_2(g_n) = s_2(g_n) mod 2` (Thue–Morse parity of consecutive gaps),
serial statistic `S^(k) = Σ u_n u_{n+k}` with `u_n = 1 − 2x_n`, scored against the exact
marginal-preserving permutation null as `z(k)`. Detector: signature PRESENT iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `r = |z_x(2)| / median{|z_x(k)| : k≠2} ≤ 0.1`.

| sieve | quiet lag | \|z(2)\| | median \|z(k≠2)\| | ratio `r` | verdict | pre-reg | holds |
|---|---|---|---|---|---|---|---|
| squarefree (control) | 2 | 2.91 | 195.97 | **0.01486** | PRESENT | present | ✓ held |
| **squarefree-except-{17,19,23}** | 2 | 16.07 | 204.40 | **0.07860** | **PRESENT** | present | ✓ held |
| **squarefree-except-{29,31,37}** | 2 | 7.52 | 199.23 | **0.03776** | **PRESENT** | present | ✓ held |
| **squarefree-except-{41,43,47}** | 2 | 5.48 | 197.55 | **0.02776** | **PRESENT** | present | ✓ held |
| odd squarefree (anchor) | 8 | 54.98 | 54.33 | 1.012 | absent | absent | ✓ held |
| primes (negative control) | 7 | 8.09 | 1.26 | 6.444 | absent | absent | ✓ held |

Gap-stream sizes (members ≤ N): squarefree 12,158,575 (max gap 10); except-{17,19,23} 12,252,558 (max
gap 10); except-{29,31,37} 12,193,494 (max gap 9); except-{41,43,47} 12,177,465 (max gap 10); odd
squarefree 8,105,700 (max gap 16); primes 1,270,607 (max gap 180). Every triple sieve keeps all even
integers and so produces odd gaps; the stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` is reused
verbatim.

## The three-prime additivity law (the new, out-of-sample-confirmed result)

Axis D5 established that a **single**-prime relaxation lifts the frozen ratio above the control by
`r(p) − r_control ≈ C·d(p)`, `d(p) = 1/p² − 1/p³`, `C ≈ 8`. Axis D6 showed two relaxations compose:
`r({p,q}) − r_control ≈ C·(d(p) + d(q))`. If a third independent relaxation remains independent,
relaxing at `p`, `q` **and** `s` should lift the ratio by `C·(d(p) + d(q) + d(s))` — the pairwise
overlap classes `{p²q²·m}` have density `∼ 1/(p²q²) ∼ 10⁻⁵–10⁻⁶` and the triple overlap `{p²q²s²·m}`
has density `∼ 10⁻⁹`, all negligible against `d(p) ∼ 10⁻³`, so the union's admitted density is
`d(p) + d(q) + d(s)` to first order. Forming the three-prime coefficient
`k({p,q,s}) = (r({p,q,s}) − r_control)/(d(p) + d(q) + d(s))` with `r_control = 0.014860` measured in
the same job:

| triple `{p,q,s}` | {17,19,23} | {29,31,37} | {41,43,47} |
|---|---|---|---|
| `d(p)+d(q)+d(s)` | 0.0076891 | 0.0028658 | 0.0015517 |
| `r({p,q,s})` | 0.07860 | 0.03776 | 0.02776 |
| `r({p,q,s}) − r_control` | 0.06374 | 0.02290 | 0.01290 |
| `k({p,q,s})` | **8.29** | **7.99** | **8.31** |

All three coefficients lie inside the pre-registered band `[5, 11]`, so `additivity_law_holds = True`.
More tellingly, they sit right on the same `C ≈ 8` measured for the single primes in Axis D5
(`k(p) = 7.86–10.2`) and the pairs in Axis D6 (`k({p,q}) = 8.33, 7.91, 8.06`). The three-prime excess
is the sum of the three single-prime excesses: **the density scaling law composes linearly to third
order.** This was committed as a numerical prediction and then confirmed on a triple — {41,43,47} —
whose primes were never measured singly or in a pair.

A back-of-envelope additivity check makes the point directly. Using the Axis-D5 single-prime law
`r(p) − r_control ≈ 8·d(p)`, the predicted triple excess is `8·(d(p)+d(q)+d(s))`: for {17,19,23} that
is `8·0.0076891 = 0.06151` vs. the measured `0.06374`; for {29,31,37} `8·0.0028658 = 0.02293` vs.
`0.02290`; for {41,43,47} `8·0.0015517 = 0.01241` vs. `0.01290`. The measured three-prime excesses
match the summed single-prime prediction to within a few percent — no sign of a pairwise or triple
interaction term at this order.

## What this changes in the arc

Axis D3 (v0.27.0) turned "destroyed at every single prime" into a quantitative statement with a
crossover at 11; Axis D4 (v0.28.0) confirmed monotone convergence past it at 17/19/23; Axis D5
(v0.29.0) upgraded the convergence to a *predictive* single-prime density scaling law
`r(p) − r_control ≈ 8·(1/p² − 1/p³)`; Axis D6 (v0.30.0) showed that law is additive across two
primes. Axis D7 shows the additivity persists to **three** primes: three simultaneous single-prime
relaxations perturb the signature by the sum of their densities, with the same coefficient and no
third-order cross term. The strict-μ² localization is therefore not only a smooth, density-controlled
approach (D5) but a **linear and additive** one in the admitted density, through at least third order
— the perturbations from distinct relaxed primes superpose. The open mechanism question (*why* strict
μ² carries a lag-2 cancellation, and *why* the coefficient is ≈ 8) is unchanged and now sharper
still: any mechanism must reproduce a per-prime excess that is linear and additive in `1/p² − 1/p³`
across one, two and three simultaneously relaxed primes.

## Reproduce

- `data/axis_d7_three_prime_additivity.json` — the sealed job's full output, bit-identical.
- `axis_d7_three_prime_additivity_job.py` — the job (self-tests the sieve/null/stride/detector and the
  cap sets `{17:2,19:2,23:2}`, `{29:2,31:2,37:2}`, `{41:2,43:2,47:2}` and their nine-prime union
  against full factorization over [1, 20000], then computes), runnable as
  `python axis_d7_three_prime_additivity_job.py --n-max 20000000 --out-dir .`.
- `test_axis_d7_three_prime_additivity.py` — pure-stdlib verifier: brute-forces the permutation-null
  formula on tiny sequences, re-derives every `E_null`/`var_null`/`z` from the raw `(n_eff, s, S)`,
  re-runs the frozen detector, asserts the six verdicts and the convergence flags, and **re-derives
  `k({p,q,s}) = (r({p,q,s}) − r_control)/(d(p) + d(q) + d(s))` from the stored ratios and checks each
  triple's coefficient lands in the pre-registered band `[5, 11]`.**

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-some-primes" indicator scales with the size of the relaxed primes
(Axis D3 crossover, Axis D4 convergence), establishes the single-prime density scaling law (Axis D5),
shows it is additive across two primes (Axis D6), or — the new contribution here — shows that
additivity persists under a simultaneous three-prime relaxation with no third-order interaction term.
That three-prime additivity result is new.
