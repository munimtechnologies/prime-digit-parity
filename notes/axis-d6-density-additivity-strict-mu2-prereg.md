# Pre-registration — Axis D6: is the density scaling law additive? A two-prime exponent-cap relaxation

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the experiment is
run. The work on this repository is produced with AI agents working in MV2's name. Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

Axis D5 (v0.29.0) established the first *quantitative, predictive* statement of the arc: relaxing the
strict squarefree (μ²) exponent cap at a **single** prime `p` — admitting `p²·m` with `p∤m` but
forbidding `p³|n` — perturbs the frozen ratio `r(p) = |z(2)| / median|z(k≠2)|` above the squarefree
control by an amount proportional to the admitted-class density `d(p) = 1/p² − 1/p³`:

    r(p) − r_control  ≈  C · d(p),   C ≈ 8,

with the dimensionless coefficient `k(p) = (r(p) − r_control)/d(p)` landing in `[5, 11]` for every
single prime measured so far (`p = 2 … 37`; `k` ranges ≈ 7.85–10.2, confirmed out-of-sample at 29,
31, 37).

That law describes one prime at a time. The obvious next question it raises is **structural: are two
independent single-prime relaxations independent?** If relaxing at `p` and relaxing at `q` perturb the
signature independently, then relaxing at **both** at once — `squarefree-except-{p,q}`, admitting the
`p²`-classes **and** the `q²`-classes — should perturb the ratio by the **sum** of the two
single-prime densities:

    r({p,q}) − r_control  ≈  C · (d(p) + d(q)).

The overlap class `{p²q²·m}` has density `≈ 1/(p²q²) ∼ 10⁻⁵–10⁻⁶`, negligible against
`d(p) ∼ 10⁻³`, so the admitted extra density of the union is `d(p)+d(q)` to first order. Axis D6 tests
this **additivity** prediction out-of-sample.

## Detector, null, observable (unchanged; graded base 2 only)

Identical to the rest of the arc and reused **verbatim** from the Axis D / D2 / D3 / D4 / D5 jobs:
consecutive gaps `g_n`, Thue–Morse parity `x_n = P_2(g_n) = s_2(g_n) mod 2`, character `u_n = 1 − 2x_n`,
serial statistic `S^(k) = Σ u_n u_{n+k}` scored against the exact marginal-preserving permutation null
as `z(k)`, base-4 stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`. Non-tautological only at base 2/4
(`notes/base2-even-admission-nontautology.md`), so graded in base 2 exclusively.

## Test sets (all to a common N = 20,000,000)

A multiplicative indicator fixed by the exponent cap `c_p` at each prime; `n` is a member iff
`v_p(n) ≤ c_p` for every prime `p`. The three pairs are disjoint and chosen so the summed density keeps
the ratio comfortably below the 0.1 detector threshold (so every pair should stay present):

- **S0 — full squarefree** (control). `c_p = 1` everywhere.
- **S-{17,19} — squarefree-except-{17,19}**. `c_17 = c_19 = 2`: admits `289|n` and `361|n` but not
  `4913|n` or `6859|n`. Summed density `d(17)+d(19) = 0.005881`.
- **S-{23,29} — squarefree-except-{23,29}**. `c_23 = c_29 = 2`: admits `529|n` and `841|n` but not
  `12167|n` or `24389|n`. Summed density `d(23)+d(29) = 0.002956`.
- **S-{31,37} — squarefree-except-{31,37}**. `c_31 = c_37 = 2`: admits `961|n` and `1369|n` but not
  `29791|n` or `50653|n`. Summed density `d(31)+d(37) = 0.001718`.
- **S-anchor — odd squarefree** (A056911). `c_2 = 0`.
- **S-neg — primes** (A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim from the graded Axis D … D5 jobs and re-self-tested (count / full gap stream / max gap against
full-factorization reference over [1, 20000] for every cap set used — including `{17:2,19:2}`,
`{23:2,29:2}`, `{31:2,37:2}` and their union — plus the mechanism contrast that each pair admits even
members and produces odd gaps, and the exact counts `Q(10⁶)=607926`, `π(10⁶)=78498`). Self-tests run
first; the job aborts non-zero on any failure, so only self-validated numbers are emitted.

## Pre-registered predictions

**Primary (signature verdicts), following the confirmed convergence / density-law reading:**

| label | set | pre-reg |
|---|---|---|
| **P-D6-sf** | full squarefree (control) | PRESENT |
| **P-D6-pair-17_19** | squarefree-except-{17,19} | **PRESENT** |
| **P-D6-pair-23_29** | squarefree-except-{23,29} | **PRESENT** |
| **P-D6-pair-31_37** | squarefree-except-{31,37} | **PRESENT** |
| **P-D6-odd** (anchor) | odd squarefree | ABSENT |
| **P-D6-primes** (neg control) | primes | ABSENT |

**Secondary (the new quantitative test — additivity of the density law).** Pre-registered: the
two-prime coefficient

    k({p,q}) = (r({p,q}) − r_control) / (d(p) + d(q)),   d(p) = 1/p² − 1/p³,

with `r_control` the control ratio measured in the *same* job, lands inside the committed band
**`[5, 11]`** — the same band validated for the single-prime law in Axis D5 — for each of the three
pairs {17,19}, {23,29}, {31,37}. The job records `density_additivity_law.additivity_law_holds`. A pass
is positive evidence that two independent prime relaxations perturb the signature **additively in the
admitted density**; it would be the first evidence that the density scaling law composes linearly.

**Decisive read / named alternatives (all are real results).**
- `convergence_holds` True **and** `additivity_law_holds` True → the clean outcome: two-prime
  relaxations stay present and their excess is the sum of the single-prime excesses (density-additive).
- **SUPER-ADDITIVE DESTRUCTION.** If any pair scores ABSENT (`r({p,q})` above 0.1) even though each
  single relaxation at the smaller primes is present, the two relaxations interact destructively;
  reported honestly and prominently.
- **PRESENT BUT NOT ADDITIVE.** If all three pairs stay present and below threshold but some `k({p,q})`
  falls outside `[5, 11]`, the convergence survives a double relaxation but the excess is not simply
  the sum of densities (sub- or super-additive); reported honestly.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the size of the relaxed prime
(Axis D3 crossover, Axis D4 convergence), establishes the single-prime density scaling law (Axis D5),
or — the new question here — tests whether that law is **additive** under a simultaneous relaxation at
two primes. That additivity test is the new contribution.
