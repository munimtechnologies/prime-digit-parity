# Pre-registration — Axis D3: is single-prime exponent-cap relaxation destructive at *every* prime, or does the lag-2 signature survive above a large-prime crossover?

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the experiment
is run. The work on this repository is produced with AI agents working in MV2's name. Code is
licensed Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc — and the tension that motivates it

Axis D (v0.24.0) and its generality replication Axis D2 (v0.26.0) showed that relaxing the strict
squarefree (μ²) exponent cap at a single prime `p` — admitting `p²·m` with `p∤m` but still
forbidding `p³|n` — destroys the base-2 gap-parity **lag-2 unique-quiet** signature, at every
prime tested so far: 2, 3, 5, 7. The scale-free detector scores the signature PRESENT iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`.

Axis D2 surfaced a **monotone trend** in the frozen, committed ratio `r(p) = |z(2)| / median|z(k≠2)|`
for `squarefree-except-p`:

| relaxed prime `p` | quiet lag | ratio `r(p)` | signature |
|---|---|---|---|
| 2 | 5 | 1.142 | absent |
| 3 | 8 | 0.773 | absent |
| 5 | 2 | 0.342 | absent |
| 7 | 2 | 0.169 | absent |
| — (strict squarefree control) | 2 | **0.015** | present |

`r(p)` falls monotonically as `p` grows. This creates a **sharp, testable tension**:

- **Strict-μ² "every prime" reading.** Relaxing the cap at *any* single prime destroys the
  signature, so `r(p) > 0.1` for all `p` — the ratio asymptotes *above* the 0.1 threshold.
- **Convergence-to-control reading.** As `p → ∞`, `squarefree-except-p → squarefree` (the admitted
  class `{p²m : p∤m}` has density `1/p² − 1/p³ → 0`), so `r(p) → 0.015`, the control value. Then
  `r(p)` must eventually cross *below* 0.1 and the detector would score the signature **PRESENT**
  for a large enough single prime.

These predict opposite things. Axis D3 tests which by pushing to the next two primes, **11 and
13**, where the trend (if it continues) would have `r` approaching the 0.1 boundary.

## Detector, null, observable (unchanged; graded base 2 only)

Identical to the rest of the arc and reused verbatim from the Axis D / D2 jobs: consecutive gaps
`g_n`, Thue–Morse parity `x_n = P_2(g_n) = s_2(g_n) mod 2`, character `u_n = 1 − 2x_n`, serial
statistic `S^(k) = Σ u_n u_{n+k}` scored against the exact marginal-preserving permutation null as
`z(k)`, base-4 stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`. Non-tautological only at base 2/4
(`notes/base2-even-admission-nontautology.md`), so graded in base 2 exclusively.

## Test sets (all to a common N = 20,000,000)

A multiplicative indicator fixed by the exponent cap `c_p` at each prime; `n` is a member iff
`v_p(n) ≤ c_p` for every prime `p`.

- **S0 — full squarefree** (control). `c_p = 1` everywhere.
- **S-11 — squarefree-except-11** (PRIMARY). `c_11 = 2`: admits `121|n` but not `1331|n`.
- **S-13 — squarefree-except-13** (replication). `c_13 = 2`: admits `169|n` but not `2197|n`.
- **S-anchor — odd squarefree** (A056911). `c_2 = 0`.
- **S-neg — primes** (A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim from the graded Axis D / D2 jobs and re-self-tested (count / full gap stream / max gap
against full-factorization reference over [1, 20000] for every cap set used — including `{11:2}`,
`{13:2}`, `{11:2,13:2}` — plus the mechanism contrast that S-11 contains 121 but not 1331 and S-13
contains 169 but not 2197, both admitting even members and producing odd gaps, and the exact
counts `Q(10⁶)=607926`, `π(10⁶)=78498`). Self-tests run first; the job aborts non-zero on any
failure, so only self-validated numbers are emitted.

## Pre-registered predictions

**Primary (signature verdicts).** Following the strict-μ² "every prime" reading:

| label | set | pre-reg |
|---|---|---|
| **P-D3-sf** | full squarefree (control) | PRESENT |
| **P-D3-ex11** (PRIMARY) | squarefree-except-11 | **ABSENT** |
| **P-D3-ex13** (replication) | squarefree-except-13 | **ABSENT** |
| **P-D3-odd** (anchor) | odd squarefree | ABSENT |
| **P-D3-primes** (neg control) | primes | ABSENT |

**Secondary (quantitative trend).** The committed Axis-D/D2 ratios are `r(5)=0.342`, `r(7)=0.169`.
Pre-registered: the monotone decrease continues, `r(11) < 0.169` and `r(13) < r(11)`, while **both
remain above 0.1** (so both are ABSENT). The job records `monotone_trend_continues` and
`both_remain_above_threshold`.

**Decisive read / named alternatives (both are real results).**
- `strict_mu2_destruction_persists_at_11_and_13` is True iff S-11 and S-13 are both ABSENT with
  the control PRESENT and anchor/neg-control ABSENT — the strict-μ² reading survives to 11/13.
- **CROSSOVER (the sharp surprise).** If `r(11) ≤ 0.1` or `r(13) ≤ 0.1`, that sieve scores the
  signature **PRESENT**. This would falsify the "relaxing the cap at every single prime destroys
  it" statement: a large-enough single-prime relaxation *preserves* the signature. It would force
  the arc's conclusion to be restated as "destroyed only below a prime-size crossover," and would
  be reported honestly and prominently as the headline.
- **TREND BREAK.** If `r(11) > 0.169` (non-monotone), the simple density-scaling picture is
  incomplete; reported honestly.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS
369, 2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2
gap-parity cancellation of a "squarefree-except-one-prime" indicator scales with the *size* of the
relaxed prime, nor tests whether a large-enough single-prime relaxation crosses back into carrying
the signature.
