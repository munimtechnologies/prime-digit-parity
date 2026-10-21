# Pre-registration — Axis D9: does the density-law coefficient depend only on admitted density, or on the structure of the admitted class? A mixed exponent-cap (c=2 vs c=3) contrast

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the experiment is
run. The work on this repository is produced with AI agents working in MV2's name. Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

Axis D5 (v0.29.0) established the single-prime density scaling law: relaxing the strict squarefree
(μ²) exponent cap at a single prime `p` — admitting the extra class `{p²·m : p∤m}` but forbidding
`p³|n` — perturbs the frozen ratio `r(p) = |z(2)| / median|z(k≠2)|` above the squarefree control by
an amount proportional to the admitted-class **density** `d₂(p) = 1/p² − 1/p³`:

    r(p) − r_control  ≈  C · d₂(p),   k(p) = (r(p) − r_control)/d₂(p) ∈ [5, 11].

Axis D6 (v0.30.0) and Axis D7 (v0.31.0) showed that law is additive across two and three primes with
the same coefficient `C ≈ 8`. Axis D8 (queued) asks whether `C` is prime-independent across a large
single-prime sweep. **Every result so far pins the law to the *density* of the admitted class.** But
density is not the only thing that changes when you relax a cap: the *structure* of the admitted
integers changes too. Axis D9 isolates that question.

## The question

Relaxing prime `p` from cap 1 to cap **3** (admit `v_p(n) ≤ 3`: allow `p³|n`, forbid `p⁴|n`) admits
two disjoint extra classes relative to the squarefree control:

- the `v_p = 2` class `{p²m : p∤m}`, density `d₂(p) = 1/p² − 1/p³` — the *same* class the cap-2
  relaxation admits; and
- the `v_p = 3` class `{p³m : p∤m}`, density `e(p) = 1/p³ − 1/p⁴`.

So the total extra density for cap 3 is `d₃(p) = d₂(p) + e(p) = 1/p² − 1/p⁴`. If the coefficient is a
property of **density alone**, the same `C` multiplies both classes and the two super-control excesses
at the *same* prime stand in a fixed, parameter-free ratio:

    ρ(p) ≡ (r_cap3(p) − r_control) / (r_cap2(p) − r_control)  =  d₃(p)/d₂(p)  =  (p+1)/p.

Because this is a ratio of excesses at the *same* prime, the coefficient `C` **cancels exactly** — so
the prediction `ρ(p) = (p+1)/p` holds whether `C` is a universal constant or drifts with prime size
(Axis D8). That makes D9 well-posed independent of D8's outcome.

If instead the `v_p = 3` class contributes with its own effective coefficient `C₃ ≠ C₂` (coefficient
depends on class **structure**, e.g. on which power of `p` is admitted, not just density), then

    r_cap3 − r_control = C₂·d₂(p) + C₃·e(p),
    ρ(p) = 1 + (C₃/C₂)·(e(p)/d₂(p)) = 1 + (C₃/C₂)·(1/p),

since `e(p)/d₂(p) = (1/p³ − 1/p⁴)/(1/p² − 1/p³) = 1/p` exactly. This gives a clean, dimensionless
estimator of the structure ratio directly from the measured excess ratio:

    C₃/C₂  =  p · (ρ(p) − 1).

Density-only ⇔ `C₃/C₂ = 1` at every prime. The signal in `ρ` is of order `1/p`, so smaller present
primes carry the most resolution; the ensemble `C₃/C₂` across several primes is the robust test.

## Detector, null, observable (unchanged; graded base 2 only)

Identical to the rest of the arc and reused **verbatim** from the Axis D … D8 jobs: consecutive gaps
`g_n`, Thue–Morse parity `x_n = P_2(g_n) = s_2(g_n) mod 2`, character `u_n = 1 − 2x_n`, serial
statistic `S^(k) = Σ u_n u_{n+k}` scored against the exact marginal-preserving permutation null as
`z(k)`, base-4 stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`. Non-tautological only at base 2/4, so
graded in base 2 exclusively. Signature present ⇔ `argmin_k|z_x(k)| = 2` AND
`|z_x(2)| ≤ 0.1·median|z_x(k≠2)|`.

## Test sets (all to a common N = 40,000,000)

A multiplicative indicator fixed by the exponent cap `c_p`; `n` is a member iff `v_p(n) ≤ c_p` for
every prime `p`. The contrast primes are the seven primes `13, 17, 19, 23, 29, 31, 37`, all past the
Axis-D3 crossover and inside the convergent regime where `squarefree-except-{p}` (cap 2) is already
confirmed PRESENT. For each `p` we add the **cap-3** sieve not measured before.

- **S0 — full squarefree** (control). `c_p = 1` everywhere. Measured in the same job → `r_control`.
- **S2-{p}** for each `p` — `c_p = 2` (allow `p²|n`, forbid `p³|n`). Extra density `d₂(p)=1/p²−1/p³`.
- **S3-{p}** for each `p` — `c_p = 3` (allow `p³|n`, forbid `p⁴|n`). Extra density `d₃(p)=1/p²−1/p⁴`.
- **S-anchor — odd squarefree** (A056911). `c_2 = 0`.
- **S-neg — primes** (A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim and re-self-tested: count / full gap stream / max gap against a full-factorization reference
over [1, 20000] for the control, `{2:0}`, each `{p:2}` and `{p:3}` cap, and the union of all caps;
the mechanism contrast that each relaxed sieve admits even members and produces odd gaps; and the
exact counts `Q(10⁶)=607926` (A005117, Gegenbauer) and `π(10⁶)=78498`. The cap-3 membership checker
is self-tested to admit `p³` and reject `p⁴` for every contrast prime. Self-tests run first; the job
aborts non-zero on any failure, so only self-validated numbers are emitted.

## Pre-registered predictions

**Primary (signature verdicts), following the confirmed convergence / density-law reading:**

| label | set | pre-reg |
|---|---|---|
| **P-D9-sf** | full squarefree (control) | PRESENT |
| **P-D9-cap2-{p}** for every `p ∈ {13,…,37}` | squarefree-except-{p}, cap 2 | **PRESENT** |
| **P-D9-cap3-{p}** for every `p ∈ {13,…,37}` | squarefree-except-{p}, cap 3 | **PRESENT** |
| **P-D9-odd** (anchor) | odd squarefree | ABSENT |
| **P-D9-primes** (neg control) | primes | ABSENT |

Each cap-2 and cap-3 relaxation adds only a small density on top of the control, so both should stay
below the 0.1 detector threshold (PRESENT). As a secondary consistency check, each cap-2 coefficient
`k₂(p) = (r_cap2(p) − r_control)/d₂(p)` should again land in the committed band **`[5, 11]`**.

**Secondary (the new quantitative test — density vs structure).** For each contrast prime the job
reports the measured excess ratio `ρ(p)` and the implied structure ratio `C₃/C₂ = p·(ρ(p) − 1)`, plus
the ensemble mean, median and spread of `C₃/C₂`. The decisive read is pre-registered on the ensemble:

- **DENSITY-ONLY (structure-independent).** If the median `C₃/C₂ ∈ [0.5, 1.5]` and no individual
  prime has `|C₃/C₂ − 1| > 3`, the coefficient is consistent with depending on admitted density alone
  (both admitted powers of `p` contribute with the same `C`). The job records
  `coefficient_is_density_only = True`.
- **STRUCTURE-DEPENDENT.** If the per-prime `C₃/C₂` deviate from 1 **consistently in one direction**
  (all above 1.5 or all below 0.5) with a median outside `[0.5, 1.5]`, the `v_p = 3` class contributes
  with a different effective coefficient than the `v_p = 2` class; reported with the sign and size of
  `C₃/C₂ − 1` (positive ⇒ higher powers push harder per unit density, negative ⇒ weaker).

**Resolution caveat (pre-registered honestly).** The signal in `ρ` is `~1/p` (7.7 % at p=13 down to
2.7 % at p=37), and the only "noise" is finite-N systematic drift of the ratios (the observable is a
deterministic count, with no sampling randomness). If the per-prime `C₃/C₂` scatter widely without a
consistent sign — i.e. the drift is larger than the structure signal at N = 4·10⁷ — the result is
reported as **inconclusive at this cutoff**, with the spread quantified and a larger-N follow-up noted.
An inconclusive outcome is a real, honestly-reported result.

**Named alternative (all are real results).**
- **ANY CONTRAST SIEVE ABSENT.** If any cap-2 or cap-3 sieve scores ABSENT (`r` above 0.1) despite the
  prime being past the crossover, the monotone-convergence picture has a hole; reported prominently
  with the offending prime and cap.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work studies the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator at all, let alone its density scaling law
(Axis D5), its additivity (D6, D7), or — the new question here — whether that law's coefficient is a
function of admitted **density alone** or also of the **structure** (admitted exponent) of the extra
class, isolated by a parameter-free `C`-cancelling cap-2 vs cap-3 excess ratio. That density-vs-
structure contrast is the new contribution.
