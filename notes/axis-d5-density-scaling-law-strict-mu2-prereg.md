# Pre-registration — Axis D5: does the single-prime exponent-cap relaxation keep converging to the squarefree control at primes 29, 31, 37, and does the super-control excess obey a density scaling law?

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the experiment
is run. The work on this repository is produced with AI agents working in MV2's name. Code is
licensed Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

Axis D3 (v0.27.0) found a **crossover** at prime 11: relaxing the strict squarefree (μ²) exponent
cap at a single prime `p` — admitting `p²·m` with `p∤m` but forbidding `p³|n` — destroys the base-2
gap-parity **lag-2 unique-quiet** signature only *below* a prime-size crossover. Axis D4 (v0.28.0)
confirmed **convergence-to-control** at 17, 19, 23: all present, `r(p)` still falling, no second
crossover. The frozen ratio `r(p) = |z(2)| / median|z(k≠2)|` for `squarefree-except-p`, across the
nine relaxed primes measured so far:

| relaxed prime `p` | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | — control |
|---|---|---|---|---|---|---|---|---|---|---|
| ratio `r(p)` | 1.142 | 0.773 | 0.342 | 0.169 | 0.0777 | 0.0603 | 0.0420 | 0.0366 | 0.0291 | **0.0149** |

The admitted extra class `{p²·m : p∤m}` has density `d(p) = 1/p² − 1/p³`. Forming the dimensionless
coefficient `k(p) = (r(p) − r_control) / d(p)` over those nine primes gives (with
`r_control ≈ 0.0149`):

    k(2)≈9.0, k(3)≈10.2, k(5)≈10.2, k(7)≈8.8, k(11)≈8.36, k(13)≈8.31, k(17)≈8.32, k(19)≈8.27, k(23)≈7.85.

`k(p)` is strikingly stable (range ≈ 7.85–10.2, slowly declining), suggesting a **quantitative
density scaling law**: the super-control excess of the destruction ratio is roughly *proportional to
the density of the admitted class*,

    r(p) − r_control  ≈  C · (1/p² − 1/p³),   C ≈ 8.

This would be the first *predictive, quantitative* statement of the arc (everything before is
present/absent plus a monotone trend). Axis D5 tests it out-of-sample at the next three primes.

## Detector, null, observable (unchanged; graded base 2 only)

Identical to the rest of the arc and reused **verbatim** from the Axis D / D2 / D3 / D4 jobs:
consecutive gaps `g_n`, Thue–Morse parity `x_n = P_2(g_n) = s_2(g_n) mod 2`, character
`u_n = 1 − 2x_n`, serial statistic `S^(k) = Σ u_n u_{n+k}` scored against the exact
marginal-preserving permutation null as `z(k)`, base-4 stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`.
Non-tautological only at base 2/4 (`notes/base2-even-admission-nontautology.md`), so graded in base 2
exclusively.

## Test sets (all to a common N = 20,000,000)

A multiplicative indicator fixed by the exponent cap `c_p` at each prime; `n` is a member iff
`v_p(n) ≤ c_p` for every prime `p`.

- **S0 — full squarefree** (control). `c_p = 1` everywhere.
- **S-29 — squarefree-except-29**. `c_29 = 2`: admits `841|n` but not `24389|n`.
- **S-31 — squarefree-except-31**. `c_31 = 2`: admits `961|n` but not `29791|n`.
- **S-37 — squarefree-except-37**. `c_37 = 2`: admits `1369|n` but not `50653|n`.
- **S-anchor — odd squarefree** (A056911). `c_2 = 0`.
- **S-neg — primes** (A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim from the graded Axis D / D2 / D3 / D4 jobs and re-self-tested (count / full gap stream / max
gap against full-factorization reference over [1, 20000] for every cap set used — including `{29:2}`,
`{31:2}`, `{37:2}` and `{29:2,31:2,37:2}` — plus the mechanism contrast that each admits even members
and produces odd gaps, and the exact counts `Q(10⁶)=607926`, `π(10⁶)=78498`). Self-tests run first;
the job aborts non-zero on any failure, so only self-validated numbers are emitted.

## Pre-registered predictions

**Primary (signature verdicts), following the confirmed convergence-to-control reading:**

| label | set | pre-reg |
|---|---|---|
| **P-D5-sf** | full squarefree (control) | PRESENT |
| **P-D5-ex29** | squarefree-except-29 | **PRESENT** |
| **P-D5-ex31** | squarefree-except-31 | **PRESENT** |
| **P-D5-ex37** | squarefree-except-37 | **PRESENT** |
| **P-D5-odd** (anchor) | odd squarefree | ABSENT |
| **P-D5-primes** (neg control) | primes | ABSENT |

**Secondary (monotone trend).** The committed Axis-D4 ratio is `r(23)=0.0291`. Pre-registered: the
monotone decrease continues, `r(29) < 0.0291`, `r(31) < r(29)`, `r(37) < r(31)`, with all three
**below the 0.1 threshold** (all PRESENT) and all moving toward the control `0.0149`
(`monotone_trend_continues`, `all_three_below_threshold`, `convergence_to_control_confirmed`).

**Tertiary (the new quantitative test — density scaling law).** Pre-registered: the dimensionless
coefficient `k(p) = (r(p) − r_control) / d(p)` with `d(p) = 1/p² − 1/p³` and `r_control` the control
ratio measured in the *same* job lands, for each of `p = 29, 31, 37`, inside the committed band
**`[5, 11]`** (bracketing the historical range 7.85–10.2 with margin). The job records
`density_scaling_law.density_law_holds`. This is a genuine out-of-sample test of a law estimated from
`p ≤ 23`; a pass is positive evidence that the super-control excess is density-linear, a fail (with
the signatures still present) means the convergence is real but not simply proportional to `d(p)`.

**Decisive read / named alternatives (all are real results).**
- `convergence_to_control_confirmed` True **and** `density_law_holds` True → the clean outcome: the
  signature survives for every single-prime relaxation at least as large as 11, and its approach to
  the control is quantitatively density-linear.
- **SECOND CROSSOVER / RE-DESTRUCTION.** If any of `r(29), r(31), r(37)` rises back **above 0.1**
  (that sieve ABSENT), the monotone-convergence picture is wrong; reported honestly and prominently.
- **CONVERGES BUT NOT DENSITY-LINEAR.** If all three stay present and monotone but some `k(p)` falls
  outside `[5, 11]`, the signature converges but the simple density law is incomplete; reported
  honestly.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the *size* of the relaxed
prime. Axis D3 established a crossover at 11; Axis D4 confirmed monotone convergence at 17/19/23;
Axis D5 tests whether that convergence continues at 29/31/37 and, for the first time, whether the
super-control excess obeys a quantitative density scaling law `r(p) − r_control ≈ C·(1/p² − 1/p³)`.
