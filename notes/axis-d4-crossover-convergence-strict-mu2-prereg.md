# Pre-registration — Axis D4: past the crossover, does `r(p)` keep falling monotonically toward the squarefree control (convergence-to-control), or does the signature get destroyed again at some larger prime?

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the experiment
is run. The work on this repository is produced with AI agents working in MV2's name. Code is
licensed Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

Axis D3 (v0.27.0) found a **crossover**: relaxing the strict squarefree (μ²) exponent cap at a single
prime `p` — admitting `p²·m` with `p∤m` but forbidding `p³|n` — destroys the base-2 gap-parity
**lag-2 unique-quiet** signature only *below* a prime-size crossover. The frozen ratio
`r(p) = |z(2)| / median|z(k≠2)|` falls monotonically and crosses below the 0.1 detector threshold
between `p=7` and `p=11`:

| relaxed prime `p` | 2 | 3 | 5 | 7 | 11 | 13 | — (control) |
|---|---|---|---|---|---|---|---|
| ratio `r(p)` | 1.142 | 0.773 | 0.342 | 0.169 | **0.0777** | **0.0603** | **0.0149** |
| signature | absent | absent | absent | absent | **present** | **present** | present |

The detector scores the signature PRESENT iff `argmin_{k∈1..8}|z_x(k)| == 2` **and**
`r ≤ 0.1`. Axis D3 confirmed the *convergence-to-control* reading over the *every-prime-destroys*
reading at 11 and 13. The natural next question: the convergence reading predicts
`r(p) → 0.0149` (the control) as `p → ∞`, because the admitted class `{p²m : p∤m}` has density
`1/p² − 1/p³ → 0`. Does `r(p)` in fact keep falling monotonically, staying present, at the next
primes — or does something non-monotone happen (a second crossover, a plateau, re-destruction)?

## Detector, null, observable (unchanged; graded base 2 only)

Identical to the rest of the arc and reused **verbatim** from the Axis D / D2 / D3 jobs: consecutive
gaps `g_n`, Thue–Morse parity `x_n = P_2(g_n) = s_2(g_n) mod 2`, character `u_n = 1 − 2x_n`, serial
statistic `S^(k) = Σ u_n u_{n+k}` scored against the exact marginal-preserving permutation null as
`z(k)`, base-4 stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`. Non-tautological only at base 2/4
(`notes/base2-even-admission-nontautology.md`), so graded in base 2 exclusively.

## Test sets (all to a common N = 20,000,000)

A multiplicative indicator fixed by the exponent cap `c_p` at each prime; `n` is a member iff
`v_p(n) ≤ c_p` for every prime `p`.

- **S0 — full squarefree** (control). `c_p = 1` everywhere.
- **S-17 — squarefree-except-17**. `c_17 = 2`: admits `289|n` but not `4913|n`.
- **S-19 — squarefree-except-19**. `c_19 = 2`: admits `361|n` but not `6859|n`.
- **S-23 — squarefree-except-23**. `c_23 = 2`: admits `529|n` but not `12167|n`.
- **S-anchor — odd squarefree** (A056911). `c_2 = 0`.
- **S-neg — primes** (A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim from the graded Axis D / D2 / D3 jobs and re-self-tested (count / full gap stream / max gap
against full-factorization reference over [1, 20000] for every cap set used — including `{17:2}`,
`{19:2}`, `{23:2}` and `{17:2,19:2,23:2}` — plus the mechanism contrast that each admits even members
and produces odd gaps, and the exact counts `Q(10⁶)=607926`, `π(10⁶)=78498`). Self-tests run first;
the job aborts non-zero on any failure, so only self-validated numbers are emitted.

## Pre-registered predictions

**Primary (signature verdicts).** Following the Axis D3 convergence-to-control reading:

| label | set | pre-reg |
|---|---|---|
| **P-D4-sf** | full squarefree (control) | PRESENT |
| **P-D4-ex17** | squarefree-except-17 | **PRESENT** |
| **P-D4-ex19** | squarefree-except-19 | **PRESENT** |
| **P-D4-ex23** | squarefree-except-23 | **PRESENT** |
| **P-D4-odd** (anchor) | odd squarefree | ABSENT |
| **P-D4-primes** (neg control) | primes | ABSENT |

**Secondary (quantitative trend).** The committed Axis-D3 ratio is `r(13)=0.0603`. Pre-registered:
the monotone decrease continues, `r(17) < 0.0603`, `r(19) < r(17)`, `r(23) < r(19)`, with all three
**below the 0.1 threshold** (so all three PRESENT) and all moving toward the control `0.0149`. The
job records `monotone_trend_continues`, `all_three_below_threshold` and
`convergence_to_control_confirmed`.

**Decisive read / named alternatives (both are real results).**
- `convergence_to_control_confirmed` is True iff all three relaxed-prime sieves are PRESENT with the
  control PRESENT and anchor/neg-control ABSENT, and `r` continues to fall monotonically. This
  confirms the crossover at 11 is a true one-way boundary and the signature survives for every prime
  at least as large.
- **SECOND CROSSOVER / RE-DESTRUCTION (the surprise).** If any of `r(17), r(19), r(23)` rises back
  **above 0.1** (that sieve scores ABSENT), the simple monotone convergence picture is wrong: the
  signature can be destroyed again at a larger prime. Reported honestly and prominently.
- **NON-MONOTONE BUT STILL PRESENT.** If all three stay below 0.1 but the sequence is not strictly
  decreasing (e.g. `r(19) > r(17)`), the density-scaling picture is incomplete; reported honestly.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the *size* of the relaxed
prime. Axis D3 established a crossover at 11; Axis D4 tests whether the ratio then converges
monotonically to the squarefree control, or whether a second crossover / re-destruction occurs.
