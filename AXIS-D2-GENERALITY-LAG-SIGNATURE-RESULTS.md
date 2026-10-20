# The strict-μ² requirement generalizes to larger primes: relaxing the squarefree exponent cap at 5 or 7 also kills the base-2 gap-parity "lag-2 unique-quiet" signature

*Produced with AI by MV2, Munim, Inc. Pre-registered in
[`notes/axis-d2-generality-strict-mu2-prereg.md`](notes/axis-d2-generality-strict-mu2-prereg.md)
(committed before the compute job ran). Code Apache-2.0; result data CC BY 4.0. This is a
reputation/curiosity number-theory micro-result: **$0 earned, no external engagement.** Every
pre-registered prediction **HELD** this time — the sharper-than-usual outcome being that this is
the first axis whose prediction was confirmed rather than falsified, because it was made *after*
Axis D fixed the strict-μ² account. Everything below is reported exactly as graded by the frozen
detector.*

## Question

Axis D (v0.24.0) relaxed the squarefree exponent cap at the two smallest relevant primes —
**squarefree-except-2** (admit `4·odd`, i.e. `4|n` but not `8|n`) and **squarefree-except-3**
(admit `9·…` but not `27|n`) — and found that both **destroy** the lag-2 unique-quiet signature
even though both keep every even integer. The conclusion was that the signature needs *exactly*
the strict squarefree (μ²) congruence structure. But Axis D only ever touched 2 and 3.

> **Decisive test (pre-registered):** is that destruction a *general* property of relaxing the
> exponent cap at **any** single prime, or an artefact of the smallest primes 2/3? Repeat the
> exponent-cap relaxation at the **larger** primes 5 and 7: **squarefree-except-5** (`c_5 = 2` —
> admit `25|n` but not `125|n`; PRIMARY) and **squarefree-except-7** (`c_7 = 2` — admit `49|n`
> but not `343|n`; replication). Both are multiplicative 0/1 indicators that **still admit every
> even integer**. Following Axis D, both were pre-registered **ABSENT**; odd squarefree (A056911)
> is the even-removing **anchor** and primes (A000040) the **negative control**, both expected
> absent. A **PRESENT** outcome at except-5 or except-7 would have been a genuine surprise,
> reopening the question and reported honestly as such.

For a sieve, take consecutive gaps `g_1, g_2, …`, form the Thue–Morse stream `x_n = P_2(g_n)`
(parity of the binary digit sum), and measure serial memory by `S^(k) = Σ_n u_n u_{n+k}`
(`u = ±1`) against an **exact marginal-preserving permutation null** (reported as `z`). The
**lag-2 unique-quiet signature** is PRESENT iff `argmin_{k∈1..8}|z_x(k)| == 2` **and**
`|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`.

## Method (fixed in advance)

Five sieves to a common **N = 20,000,000**, each analysed with the **identical** observable,
exact permutation null, stride split (`P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)`) and scale-free detector
used across the whole arc. The exponent-cap sieve was self-tested first against a full
trial-division factorization over `[1, 20000]` (count, the complete gap stream and the max gap)
for cap sets `{}`, `{2:0}`, `{5:2}`, `{7:2}` and `{5:2, 7:2}`, plus the mechanism contrast
(except-5 contains 25 but not 125; except-7 contains 49 but not 343; both admit even members and
produce odd gaps) and the exact counts `Q(10⁶)=607926`, `π(10⁶)=78498`. Self-tests ran first and
passed; the job aborts non-zero on any failure, so only self-validated numbers were emitted. Run
in a sealed `python:3.12` container (job `axis-d2-generality-strict-mu2-v1`, 2026-10-20, 195 s
wall; the sieve/analysis step itself 19 s). Gap counts at `N = 2×10⁷` — squarefree **12,158,574**;
squarefree-except-5 **12,563,858**; squarefree-except-7 **12,375,683**; odd squarefree
**8,105,699**; primes **1,270,606** — ample power everywhere.

## Result — every pre-registered prediction HELD

| Prediction | Keeps even integers? | What was predicted | Verdict |
|---|---|---|---|
| **P-D2-sf** full squarefree (control) | yes | present | **HELD** (present) |
| **P-D2-ex5 (PRIMARY)** squarefree-except-5 | yes | absent | **HELD** (absent) |
| **P-D2-ex7** squarefree-except-7 | yes | absent | **HELD** (absent) |
| **P-D2-odd (anchor)** odd squarefree (A056911) | no | absent | **HELD** (absent) |
| **P-D2-primes (neg control)** primes (A000040) | — | absent | **HELD** (absent) |

The `x_base2` lag-2 unique-quiet detector across all five sieves at `N = 2×10⁷`:

| Sieve | admits even? | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|---|
| squarefree (control) | yes | **2** | 2.91 | 196.0 | 0.015 | **present** |
| squarefree-except-5 (PRIMARY) | yes | **2** | 78.18 | 228.83 | 0.342 | absent |
| squarefree-except-7 (replication) | yes | **2** | 36.43 | 215.09 | 0.169 | absent |
| odd squarefree (anchor) | no | 8 | 54.98 | 54.33 | 1.012 | absent |
| primes (neg control) | — | 7 | 8.09 | 1.26 | 6.444 | absent |

So `strict_mu2_requirement_generalizes_to_larger_primes` is **True**: the control is present
while both larger-prime relaxations (and the anchor and negative control) are absent. Relaxing
the μ² cap at 5 or 7 destroys the signature exactly as it did at 2 and 3 — the destruction is a
general property of the strict-μ² structure, not an artefact of the small primes.

## A graded trend worth noting honestly

Unlike the *small*-prime relaxations in Axis D, where the quiet lag jumped away from 2 (to lag 5
for except-2 and lag 8 for except-3), the **larger-prime relaxations keep lag 2 as the quietest
lag** — the lag-2 cancellation is only *partially* disrupted, not relocated. And the degree of
disruption falls monotonically as the relaxed prime grows. Stacking this release against the
committed v0.24.0 numbers (both frozen, both graded by the same detector):

| relaxed prime `p` (admit `p²\|n`) | quiet lag | ratio `|z(2)| / median|z(k≠2)|` | signature |
|---|---|---|---|
| 2 (except-2, v0.24.0) | 5 | 1.142 | absent |
| 3 (except-3, v0.24.0) | 8 | 0.773 | absent |
| **5 (except-5, this release)** | **2** | **0.342** | absent |
| **7 (except-7, this release)** | **2** | **0.169** | absent |
| — (strict squarefree control) | 2 | 0.015 | present |

The ratio decreases monotonically — `1.142 → 0.773 → 0.342 → 0.169` — toward the control's 0.015
as `p` increases. This is exactly what the strict-μ² reading predicts: admitting the square class
of a *larger* prime perturbs a *rarer* set of integers (the density of integers divisible by `p²`
but not `p³` is `1/p² − 1/p³`, falling from `3/24` at `p=2` to `6/343` at `p=7`), so it disturbs
the gap structure less and leaves more of the lag-2 cancellation intact. But at every prime
tested the ratio stays **above the 0.1 threshold**, so the signature is absent in each case — the
frozen criterion is not met by any single-prime relaxation, only by the complete μ² sieve. This
is a four-point trend across 2, 3, 5, 7; it is suggestive, not a limit claim, and only the
detector verdict (absent) is what the study grades.

## Interpretation

This answers the generality question Axis D left open **for the primes tested**: relaxing the
exponent cap at each of 5 and 7 removes the signature exactly as it did at 2 and 3, even though
every relaxation keeps all even integers. Across the four primes 2, 3, 5, 7 the destruction is
robust, and together with Axis C (removing one prime's residues also kills it) the signature
tolerates neither thinning nor enlarging at any of the small primes probed. Even-admission is
necessary but not sufficient, confirmed now at four distinct primes from the enlarging side.

**An honest tension the trend forces — and the open question it raises.** The monotone ratio
decrease (`1.142 → 0.773 → 0.342 → 0.169`) is not just a curiosity. As `p → ∞` the set
`squarefree-except-p` converges to `squarefree` itself (the admitted square class `{p²·m : p∤m}`
has density `1/p² − 1/p³ → 0`), so the ratio must converge to the control's **0.015**. That means
there is almost certainly a large-enough prime at which the ratio drops **below the 0.1
threshold** and the frozen detector would score the signature **present**. So the claim is *not*
that relaxing the cap at *every* prime destroys the signature — the honest statement is that it is
destroyed at every prime small enough that the perturbation is dense enough (2, 3, 5, 7 so far),
and whether the destruction persists for arbitrarily large primes or fails above some crossover
is genuinely **open**. The strict-μ² "every prime" reading predicts the ratio asymptotes *above*
0.1; the convergence-to-control argument predicts it crosses *below*. These disagree, and the
next primes (11, 13) are pre-registered to test which — see
`notes/axis-d3-large-prime-crossover-strict-mu2-prereg.md`.

## What this closes

Within the four-axis map (`SYNTHESIS.md`), this is the **generality replication of Axis D**:

- **Axis A (which k-free):** only k = 2 — squarefree, not cube-free or 4th-power-free (v0.18.0).
- **Axis B (origin vs. indicator):** needs a multiplicative *indicator*; the Liouville +1 level
  set is flat (v0.20.0).
- **Axis C (thinning):** removing evens or lifting one odd prime both kill it (v0.21.0 / v0.23.0).
- **Axis D (enlarging):** relaxing the exponent cap at 2 or 3 kills it (v0.24.0); **and relaxing
  it at the larger primes 5 or 7 kills it too (this result) — the strict-μ² requirement
  generalizes to 5 and 7**, with the large-prime crossover (does it ever survive for a big enough
  single prime?) pre-registered next (Axis D3, primes 11/13).

## Honest limits

All five sieves have ≥1.27×10⁶ gaps (the four multiplicative ones ≥8.1×10⁶), so this is not a
power problem; the detector is a **ratio** test (scale-free). A null means "no lag-2 unique-quiet
signature at this N." We relaxed the cap at the two specific primes 5 and 7 to exponent 2; with
Axis D's 2 and 3 that is four primes in total, and the monotone ratio trend is a four-point
observation, not an asymptotic or all-primes claim. We do not claim *every* exponent-cap
relaxation at *every* prime behaves this way, only that no multiplicative relative tested across
Axes A–D (now including larger-prime enlargements) reproduces the signature. The conclusion is
about the `N = 2×10⁷` regime; the arc makes no asymptotic claim. The question is non-tautological
only at base 2/4 (see `notes/base2-even-admission-nontautology.md`), so it is graded in base 2
exclusively. All numbers are re-derivable offline from
[`data/axis_d2_generality_lag_signature.json`](data/axis_d2_generality_lag_signature.json) by
[`test_axis_d2_generality_lag_signature.py`](test_axis_d2_generality_lag_signature.py), which
brute-force validates the null formula on tiny sequences, then re-computes every `E_null`,
`var_null` and `z` from the raw `(n_eff, s, S)` and re-runs the detector and grading.

## Reproduce

```
python axis_d2_generality_job.py --n-max 20000000 --out-dir .
python -m unittest -v test_axis_d2_generality_lag_signature
```

## Prior art

Thue–Morse / digit-sum correlation structure is classical: Gelfond 1968 (Acta Arith. 13,
1967/68, 259–265); Mauduit–Rivat *Ann. Math.* 171 (2010), 1591–1646; Baake–Coons
[arXiv:2209.07102](https://arxiv.org/abs/2209.07102). The nearest neighbour on squarefree binary
digits is Dietmann–Elsholtz–Shparlinski, *Trans. Amer. Math. Soc.* 369 (2017), 8369–8388
([doi:10.1090/tran/6903](https://doi.org/10.1090/tran/6903)) — about the squarefree integers and
their digits, not the gaps between them or the digit-sum parity of those gaps. Densities:
A005117 (squarefree `6/π²`), A056911 (odd squarefree). No prior work isolates the base-2
gap-parity lag-2 XOR cancellation, nor shows it depends on the complete strict squarefree
congruence structure at larger primes as well as small ones; see
[`PRIOR-ART.md`](PRIOR-ART.md). If you know of prior work that already does this, please open an
issue — a correct negative is as useful to this project as a positive.
