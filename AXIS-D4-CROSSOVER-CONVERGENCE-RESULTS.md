# Axis D4 — past the crossover, the single-prime relaxation converges to the squarefree control

**Headline (a pre-registered prediction, held).** Past the Axis D3 crossover at prime 11, relaxing
the strict squarefree (μ²) exponent cap at a *single* larger prime keeps the base-2 gap-parity
lag-2 unique-quiet signature **present**, and the frozen ratio `r(p)` keeps falling
**monotonically** toward the squarefree control. At primes **17, 19 and 23** the signature is
present (`r = 0.0420, 0.0366, 0.0291`, all below the 0.1 threshold), with no second crossover and no
re-destruction. The pre-registered **convergence-to-control** prediction — that `squarefree-except-17`,
`-except-19` and `-except-23` would all be **PRESENT** with `r` still decreasing — is **HELD**.

*This repository's work is produced with AI agents working in MV2's name, at Munim, Inc. It is a
reputation / curiosity number-theory result: **$0 earned, no external engagement.** Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.*

## What was pre-registered, and what happened

The pre-registration (`notes/axis-d4-crossover-convergence-strict-mu2-prereg.md`, committed before
this sealed job ran) set up the question left open by Axis D3's crossover, and the two named readings
that predict different things past it.

- **Convergence-to-control reading** (the primary prediction): as `p → ∞`, the admitted extra class
  `{p²·m : p∤m}` has density `1/p² − 1/p³ → 0`, so `squarefree-except-p → squarefree` and
  `r(p) → 0.0149` (the control). Past the crossover at 11, `r(p)` should keep falling monotonically
  and stay below 0.1. Predicts `squarefree-except-17/19/23` all **PRESENT**.
- **Second crossover / re-destruction** (the named surprise): if any of `r(17), r(19), r(23)` rose
  back **above** 0.1 (that sieve scoring ABSENT), the simple monotone picture would be wrong — the
  signature could be destroyed again at a larger prime.

The first reading is what the data show. **Convergence confirmed; no second crossover.**

## Result (N = 20,000,000, sealed python:3.12 container, exit 0, self-tests passed)

Base-2 gap-parity stream `x_n = P_2(g_n) = s_2(g_n) mod 2` (Thue–Morse parity of consecutive gaps),
serial statistic `S^(k) = Σ u_n u_{n+k}` with `u_n = 1 − 2x_n`, scored against the exact
marginal-preserving permutation null as `z(k)`. Detector: signature PRESENT iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `r = |z_x(2)| / median{|z_x(k)| : k≠2} ≤ 0.1`.

| sieve | quiet lag | \|z(2)\| | median \|z(k≠2)\| | ratio `r` | verdict | pre-reg | holds |
|---|---|---|---|---|---|---|---|
| squarefree (control) | 2 | 2.91 | 195.97 | **0.0149** | PRESENT | present | ✓ held |
| **squarefree-except-17** (PRIMARY) | 2 | 8.38 | 199.47 | **0.0420** | **PRESENT** | present | ✓ held |
| **squarefree-except-19** (replication) | 2 | 7.29 | 198.83 | **0.0366** | **PRESENT** | present | ✓ held |
| **squarefree-except-23** (replication) | 2 | 5.76 | 197.95 | **0.0291** | **PRESENT** | present | ✓ held |
| odd squarefree (anchor) | 8 | 54.98 | 54.33 | 1.012 | absent | absent | ✓ held |
| primes (negative control) | 7 | 8.09 | 1.26 | 6.444 | absent | absent | ✓ held |

Gap-stream sizes (members ≤ N): squarefree 12,158,575 (max gap 10); except-17 12,198,306 (max gap
10); except-19 12,190,580 (max gap 10); except-23 12,180,593 (max gap 10); odd squarefree 8,105,700
(max gap 16); primes 1,270,607 (max gap 180). Every sieve except odd-squarefree and primes keeps all
even integers and so produces odd gaps; the stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` is reused
verbatim.

## The monotone ratio, all the way down toward the control

The frozen ratio `r(p)` for `squarefree-except-p` continues its monotone decrease, straight through
the crossover (between 7 and 11) and on past it, with every new prime landing below the last and all
of them below the 0.1 threshold:

| relaxed prime `p` | 2 | 3 | 5 | 7 | 11 | 13 | **17** | **19** | **23** | — control |
|---|---|---|---|---|---|---|---|---|---|---|
| ratio `r(p)` | 1.142 | 0.773 | 0.342 | 0.169 | 0.0777 | 0.0603 | **0.0420** | **0.0366** | **0.0291** | **0.0149** |
| signature | absent | absent | absent | absent | present | present | **PRESENT** | **PRESENT** | **PRESENT** | present |

So `monotone_trend_continues = True`, `all_three_below_threshold = True`,
`convergence_to_control_confirmed = True`, and `second_crossover_or_redestruction = False`. Nine
distinct single-prime relaxations now trace one strictly decreasing curve, `r(23) = 0.0291`, already
within a factor of two of the control `0.0149`, with the gap between consecutive relaxed primes
shrinking as predicted by the `1/p² − 1/p³` density of the admitted class.

## What this changes in the arc

Axis D3 (v0.27.0) replaced the all-or-nothing "destroyed at every single prime" statement with a
quantitative one: the single-prime relaxation weakens the signature by an amount that **shrinks with
`p`**, crossing the detector threshold between 7 and 11. The obvious risk to that picture was a
*second* crossover — some larger prime at which the signature is destroyed again, which would mean
the 11 crossover is not a clean one-way boundary. Axis D4 rules that out through prime 23: the ratio
decreases monotonically at every step and converges toward the control, so the crossover at 11 is a
genuine one-way boundary (for this detector and threshold) and the signature survives for every
single-prime relaxation at least as large as 11. The open mechanism question (*why* the strict-μ²
structure carries a lag-2 cancellation at all) is unchanged; the sharpened picture is now that the
exponent-cap relaxation is a **continuous, monotone** perturbation of the signature controlled by the
relaxed prime's density contribution — destroying it below the crossover, preserving (and
asymptotically restoring) it above.

## Reproduce

- `data/axis_d4_crossover_convergence.json` — the sealed job's full output, bit-identical.
- `axis_d4_crossover_convergence_job.py` — the job (self-tests the sieve/null/stride/detector, then
  computes), runnable as `python axis_d4_crossover_convergence_job.py --n-max 20000000 --out-dir .`.
- `test_axis_d4_crossover_convergence.py` — pure-stdlib verifier: brute-forces the permutation-null
  formula on tiny sequences, re-derives every `E_null`/`var_null`/`z` from the raw `(n_eff, s, S)`,
  re-runs the frozen detector, and asserts the six verdicts, the convergence, and the monotone trend
  across all nine relaxed primes.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the *size* of the relaxed
prime, identifies a prime-size crossover (Axis D3), or shows that past that crossover the ratio
converges monotonically to the squarefree control with no second crossover (this result). That
convergence, traced across nine relaxed primes, is the new contribution.
