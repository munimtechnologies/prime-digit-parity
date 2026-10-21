# Pre-registration — Axis D7: does the density law stay additive for three primes? A three-prime exponent-cap relaxation

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the experiment is
run. The work on this repository is produced with AI agents working in MV2's name. Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

Axis D5 (v0.29.0) established the single-prime density scaling law: relaxing the strict squarefree
(μ²) exponent cap at a **single** prime `p` — admitting `p²·m` with `p∤m` but forbidding `p³|n` —
perturbs the frozen ratio `r(p) = |z(2)| / median|z(k≠2)|` above the squarefree control by an amount
proportional to the admitted-class density `d(p) = 1/p² − 1/p³`:

    r(p) − r_control  ≈  C · d(p),   C ≈ 8,   k(p) = (r(p) − r_control)/d(p) ∈ [5, 11]  (p = 2 … 37).

Axis D6 (v0.30.0) then showed that law is **additive across two primes**: relaxing at both primes of
a pair `{p,q}` at once perturbs the ratio by the *sum* of the two single-prime densities,

    r({p,q}) − r_control  ≈  C · (d(p) + d(q)),

with the two-prime coefficient `k({p,q})` landing in `[5, 11]` for all three pairs {17,19}, {23,29},
{31,37} (`k = 8.33, 7.91, 8.06`). Two independent single-prime relaxations perturb the signature
additively in the admitted density.

The obvious next structural question: **does additivity extend to three primes, or does a
higher-order (pairwise or triple) interaction term appear at third order?** If three independent
single-prime relaxations remain independent, then relaxing at `p`, `q` **and** `s` at once —
`squarefree-except-{p,q,s}`, admitting the `p²`-, `q²`- and `s²`-classes — should perturb the ratio
by the **sum** of the three single-prime densities:

    r({p,q,s}) − r_control  ≈  C · (d(p) + d(q) + d(s)).

The pairwise-overlap classes `{p²q²·m}` have density `≈ 1/(p²q²) ∼ 10⁻⁵–10⁻⁶` and the triple overlap
`{p²q²s²·m}` has density `∼ 10⁻⁹`, all negligible against `d(p) ∼ 10⁻³`, so the admitted extra
density of the union is `d(p)+d(q)+d(s)` to first order. Axis D7 tests this three-prime additivity
prediction out-of-sample.

## Detector, null, observable (unchanged; graded base 2 only)

Identical to the rest of the arc and reused **verbatim** from the Axis D / D2 / D3 / D4 / D5 / D6
jobs: consecutive gaps `g_n`, Thue–Morse parity `x_n = P_2(g_n) = s_2(g_n) mod 2`, character
`u_n = 1 − 2x_n`, serial statistic `S^(k) = Σ u_n u_{n+k}` scored against the exact
marginal-preserving permutation null as `z(k)`, base-4 stride split
`P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`. Non-tautological only at base 2/4
(`notes/base2-even-admission-nontautology.md`), so graded in base 2 exclusively.

## Test sets (all to a common N = 20,000,000)

A multiplicative indicator fixed by the exponent cap `c_p` at each prime; `n` is a member iff
`v_p(n) ≤ c_p` for every prime `p`. The three triples are disjoint and chosen so the summed density
keeps the ratio comfortably below the 0.1 detector threshold (so every triple should stay present):

- **S0 — full squarefree** (control). `c_p = 1` everywhere.
- **S-{17,19,23}** — `c_17 = c_19 = c_23 = 2`: admits `289|n`, `361|n`, `529|n` but not the cubes.
  Summed density `d(17)+d(19)+d(23) = 0.007689`. Predicted `r ≈ 0.0149 + 8·0.007689 ≈ 0.076`.
- **S-{29,31,37}** — `c_29 = c_31 = c_37 = 2`: admits `841|n`, `961|n`, `1369|n` but not the cubes.
  Summed density `d(29)+d(31)+d(37) = 0.002866`. Predicted `r ≈ 0.038`.
- **S-{41,43,47}** — `c_41 = c_43 = c_47 = 2`: admits `1681|n`, `1849|n`, `2209|n` but not the cubes.
  Summed density `d(41)+d(43)+d(47) = 0.001551`. Predicted `r ≈ 0.027`. These three primes were never
  measured singly, so this triple is a fully out-of-sample test of both convergence and additivity.
- **S-anchor — odd squarefree** (A056911). `c_2 = 0`.
- **S-neg — primes** (A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim from the graded Axis D … D6 jobs and re-self-tested (count / full gap stream / max gap
against full-factorization reference over [1, 20000] for every cap set used — including each triple
`{17:2,19:2,23:2}`, `{29:2,31:2,37:2}`, `{41:2,43:2,47:2}` and the union of all nine relaxed primes —
plus the mechanism contrast that each triple admits even members and produces odd gaps, and the exact
counts `Q(10⁶)=607926`, `π(10⁶)=78498`). Self-tests run first; the job aborts non-zero on any
failure, so only self-validated numbers are emitted.

## Pre-registered predictions

**Primary (signature verdicts), following the confirmed convergence / additive density-law reading:**

| label | set | pre-reg |
|---|---|---|
| **P-D7-sf** | full squarefree (control) | PRESENT |
| **P-D7-triple-17_19_23** | squarefree-except-{17,19,23} | **PRESENT** |
| **P-D7-triple-29_31_37** | squarefree-except-{29,31,37} | **PRESENT** |
| **P-D7-triple-41_43_47** | squarefree-except-{41,43,47} | **PRESENT** |
| **P-D7-odd** (anchor) | odd squarefree | ABSENT |
| **P-D7-primes** (neg control) | primes | ABSENT |

**Secondary (the new quantitative test — three-prime additivity of the density law).**
Pre-registered: the three-prime coefficient

    k({p,q,s}) = (r({p,q,s}) − r_control) / (d(p) + d(q) + d(s)),   d(p) = 1/p² − 1/p³,

with `r_control` the control ratio measured in the *same* job, lands inside the committed band
**`[5, 11]`** — the same band validated for the single-prime law (Axis D5) and the two-prime law
(Axis D6) — for each of the three triples. The job records
`density_additivity_law.additivity_law_holds`. A pass is positive evidence that three independent
prime relaxations perturb the signature **additively in the admitted density** with no detectable
third-order interaction term.

**Decisive read / named alternatives (all are real results).**
- `convergence_holds` True **and** `additivity_law_holds` True → the clean outcome: three-prime
  relaxations stay present and their excess is the sum of the three single-prime excesses
  (density-additive to third order).
- **SUPER-ADDITIVE DESTRUCTION.** If any triple scores ABSENT (`r({p,q,s})` above 0.1) even though
  each single and double relaxation at these primes is present, the three relaxations interact
  destructively; reported honestly and prominently.
- **PRESENT BUT NOT ADDITIVE.** If all three triples stay present and below threshold but some
  `k({p,q,s})` falls outside `[5, 11]`, the convergence survives a triple relaxation but the excess
  is not simply the sum of the three densities (a pairwise or triple cross term appears); reported
  honestly, with the size of the departure from additivity quantified.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-some-primes" indicator scales with the size of the relaxed
primes (Axis D3 crossover, Axis D4 convergence), establishes the single-prime density scaling law
(Axis D5), shows it is additive across two primes (Axis D6), or — the new question here — tests
whether that additivity persists for a *simultaneous relaxation at three primes* or breaks down into
a higher-order interaction. That three-prime additivity test is the new contribution.
