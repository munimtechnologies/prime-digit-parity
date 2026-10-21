# Pre-registration — Axis D8: is the coefficient C prime-independent, or does it drift? A dense single-prime sweep across the convergent regime

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the experiment is
run. The work on this repository is produced with AI agents working in MV2's name. Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

Axis D5 (v0.29.0) established the single-prime density scaling law: relaxing the strict squarefree
(μ²) exponent cap at a **single** prime `p` — admitting `p²·m` with `p∤m` but forbidding `p³|n` —
perturbs the frozen ratio `r(p) = |z(2)| / median|z(k≠2)|` above the squarefree control by an amount
proportional to the admitted-class density `d(p) = 1/p² − 1/p³`:

    r(p) − r_control  ≈  C · d(p),   k(p) = (r(p) − r_control)/d(p) ∈ [5, 11]   (p = 2 … 37).

Axis D6 (v0.30.0) and Axis D7 (v0.31.0) showed that law is **additive** across two and three primes,
with the same coefficient `C ≈ 8` (`k({p,q}) = 8.33, 7.91, 8.06`; `k({p,q,s}) = 8.29, 7.99, 8.31`).

Across those runs the single-prime coefficients `k(p)` measured for `p = 2 … 37` scatter over roughly
`7.86 … 10.2`, visibly tighter and closer to 8 for the larger primes. That raises the natural
question this axis tests directly: **is `C` a single universal constant that the per-prime `k(p)`
converge to as `p` grows, or does `k(p)` drift systematically with prime size?** Axis D8 answers it
with a dense single-prime sweep across a block of primes **larger than any measured singly so far**
(41 … 97), all comfortably past the Axis-D3 crossover and inside the convergent regime.

## Detector, null, observable (unchanged; graded base 2 only)

Identical to the rest of the arc and reused **verbatim** from the Axis D / D2 / … / D7 jobs:
consecutive gaps `g_n`, Thue–Morse parity `x_n = P_2(g_n) = s_2(g_n) mod 2`, character
`u_n = 1 − 2x_n`, serial statistic `S^(k) = Σ u_n u_{n+k}` scored against the exact
marginal-preserving permutation null as `z(k)`, base-4 stride split
`P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`. Non-tautological only at base 2/4, so graded in base 2 exclusively.

## Test sets (all to a common N = 20,000,000)

A multiplicative indicator fixed by the exponent cap `c_p`; `n` is a member iff `v_p(n) ≤ c_p` for
every prime `p`. The sweep primes are the thirteen primes `41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
83, 89, 97` — **none measured singly before** (Axis D5 stopped at 37), so the whole sweep is
out-of-sample. For each sweep prime `p` the sieve `squarefree-except-{p}` sets `c_p = 2` (admits
`p²|n`, forbids `p³|n`), `c_q = 1` otherwise.

- **S0 — full squarefree** (control). `c_p = 1` everywhere.
- **S-{p}** for each `p ∈ {41,…,97}` — `c_p = 2`. Admitted extra density `d(p) = 1/p² − 1/p³`.
- **S-anchor — odd squarefree** (A056911). `c_2 = 0`.
- **S-neg — primes** (A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim from the graded Axis D … D7 jobs and re-self-tested (count / full gap stream / max gap
against full-factorization reference over [1, 20000] for the squarefree control, `{2:0}`, each sweep
prime cap `{p:2}` and the union of all thirteen sweep primes, plus the mechanism contrast that each
`{p:2}` sieve admits even members and produces odd gaps, and the exact counts `Q(10⁶)=607926`,
`π(10⁶)=78498`). Self-tests run first; the job aborts non-zero on any failure, so only self-validated
numbers are emitted.

## Pre-registered predictions

**Primary (signature verdicts), following the confirmed convergence / additive density-law reading:**

| label | set | pre-reg |
|---|---|---|
| **P-D8-sf** | full squarefree (control) | PRESENT |
| **P-D8-sweep-{p}** for every `p ∈ {41,…,97}` | squarefree-except-{p} | **PRESENT** |
| **P-D8-odd** (anchor) | odd squarefree | ABSENT |
| **P-D8-primes** (neg control) | primes | ABSENT |

Each sweep prime is past the crossover and deep in the convergent regime, so its ratio should stay
well below the 0.1 detector threshold and its coefficient `k(p) = (r(p) − r_control)/d(p)` should
land inside the committed band **`[5, 11]`** — the same band validated for the single-prime law
(Axis D5), the two-prime law (Axis D6) and the three-prime law (Axis D7), with `r_control` measured
in the *same* job.

**Secondary (the new quantitative test — universality vs drift of `C`).** The job fits an ordinary
least-squares line `k(p) = a + b·p` over the thirteen sweep primes and reports the slope `b`, its
standard error `SE(b)`, the sweep mean and median of `k(p)`, and the spread `max k − min k`. The
decisive read is pre-registered on the slope:

- **UNIVERSAL C (no drift).** If `|b| ≤ 2·SE(b)` (slope consistent with zero at ~2σ), the per-prime
  coefficient has no detectable dependence on prime size across the large-prime regime: `C` is a
  single constant `≈ 8` that `k(p)` has converged to. The job records `coefficient_is_universal = True`.
- **C DRIFTS WITH PRIME SIZE.** If `|b| > 2·SE(b)`, the coefficient drifts systematically with `p`;
  reported honestly with the sign and size of the slope, and with the extrapolated `k(∞) = a + b·∞`
  limit discussed (a negative `b` would mean `k(p)` is still decreasing toward a limit below the
  low-prime values, consistent with the Axis-D5 observation that larger primes cluster nearer 8).

Both are real results. The spread and mean/median are reported regardless, so even a null slope
quantifies how tightly the large-prime coefficients cluster on `C ≈ 8`.

**Named alternative (all are real results).**
- **ANY SWEEP PRIME ABSENT.** If any `squarefree-except-{p}` scores ABSENT (`r(p)` above 0.1) despite
  being well past the crossover, the monotone-convergence picture of the arc has a hole; reported
  prominently with the offending prime.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the size of the relaxed prime
(Axis D3 crossover, Axis D4 convergence), establishes the single-prime density scaling law (Axis D5),
shows it additive across two and three primes (Axis D6, D7), or — the new question here — pins down
whether the density-law coefficient `C` is a single prime-independent constant or drifts with prime
size via a dense out-of-sample single-prime sweep. That universality test is the new contribution.
