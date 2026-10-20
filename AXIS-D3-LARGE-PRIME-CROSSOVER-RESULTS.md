# Axis D3 — the strict-μ² exponent-cap requirement has a large-prime crossover

**Headline (a pre-registered prediction, falsified).** Relaxing the strict squarefree (μ²)
exponent cap at a *single* prime does **not** destroy the base-2 gap-parity lag-2 unique-quiet
signature at *every* prime. At the larger primes **11 and 13** the signature **returns**. The
pre-registered PRIMARY prediction — that `squarefree-except-11` and `squarefree-except-13` would
both be **ABSENT**, following the "destroyed at every single prime" reading — is **FALSIFIED**. The
destruction holds only **below a prime-size crossover** that sits between 7 and 11.

*This repository's work is produced with AI agents working in MV2's name, at Munim, Inc. It is a
reputation / curiosity number-theory result: **$0 earned, no external engagement.** Code is licensed
Apache-2.0; result data is licensed CC BY 4.0.*

## What was pre-registered, and what happened

The pre-registration (`notes/axis-d3-large-prime-crossover-strict-mu2-prereg.md`, committed before
this sealed job ran) set up two readings that predict opposite things, and ran the single decisive
experiment that separates them.

- **Strict-μ² "every prime" reading** (the primary prediction): relaxing the cap at *any* single
  prime `p` — admitting `p²·m` with `p∤m` but still forbidding `p³|n` — destroys the signature, so
  the frozen ratio `r(p) = |z(2)| / median|z(k≠2)|` stays **above** the 0.1 detector threshold for
  all `p`. Predicts `squarefree-except-11` and `squarefree-except-13` **ABSENT**.
- **Convergence-to-control reading**: as `p → ∞`, the admitted class `{p²m : p∤m}` has density
  `1/p² − 1/p³ → 0`, so `squarefree-except-p → squarefree` and `r(p) → 0.015` (the control). Then
  `r(p)` must eventually cross **below** 0.1 and the signature would score **PRESENT**.

The second reading is what the data show. **A crossover was detected.**

## Result (N = 20,000,000, sealed python:3.12 container, exit 0, self-tests passed)

Base-2 gap-parity stream `x_n = P_2(g_n) = s_2(g_n) mod 2` (Thue–Morse parity of consecutive gaps),
serial statistic `S^(k) = Σ u_n u_{n+k}` with `u_n = 1 − 2x_n`, scored against the exact
marginal-preserving permutation null as `z(k)`. Detector: signature PRESENT iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `r = |z_x(2)| / median{|z_x(k)| : k≠2} ≤ 0.1`.

| sieve | quiet lag | \|z(2)\| | median \|z(k≠2)\| | ratio `r` | verdict | pre-reg | holds |
|---|---|---|---|---|---|---|---|
| squarefree (control) | 2 | 2.91 | 195.97 | **0.0149** | PRESENT | present | ✓ held |
| **squarefree-except-11** (PRIMARY) | 2 | 15.87 | 204.23 | **0.0777** | **PRESENT** | absent | **✗ falsified** |
| **squarefree-except-13** (replication) | 2 | 12.20 | 202.27 | **0.0603** | **PRESENT** | absent | **✗ falsified** |
| odd squarefree (anchor) | 8 | 54.98 | 54.33 | 1.012 | absent | absent | ✓ held |
| primes (negative control) | 7 | 8.09 | 1.26 | 6.444 | absent | absent | ✓ held |

Gap-stream sizes (members ≤ N): squarefree 12,158,575 (max gap 10); except-11 12,250,682 (max gap
9); except-13 12,225,375 (max gap 10); odd squarefree 8,105,700 (max gap 16); primes 1,270,607 (max
gap 180). Every sieve except odd-squarefree and primes keeps all even integers and so produces odd
gaps; the stride split `P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` is reused verbatim.

## The monotone ratio, all the way down through the crossover

The frozen ratio `r(p)` for `squarefree-except-p` continues its monotone decrease and crosses below
the 0.1 threshold between 7 and 11:

| relaxed prime `p` | 2 | 3 | 5 | 7 | **11** | **13** |
|---|---|---|---|---|---|---|
| ratio `r(p)` | 1.142 | 0.773 | 0.342 | 0.169 | **0.0777** | **0.0603** |
| signature | absent | absent | absent | absent | **PRESENT** | **PRESENT** |

So `monotone_trend_continues = True`, `both_remain_above_threshold = False`,
`crossover_detected = True`, and `strict_mu2_destruction_persists_at_11_and_13 = False`. The crossing
of the 0.1 line happens strictly between `p=7` (r=0.169, absent) and `p=11` (r=0.0777, present).

## What this changes in the arc

Axes A–D2 built a clean map in which the lag-2 unique-quiet signature appeared to need *exactly* the
strict squarefree congruence structure — tolerating neither thinning a residue class (Axis C) nor
enlarging the exponent cap at any prime (Axis D for 2/3, Axis D2 for 5/7). Axis D3 shows that last
clause was an artefact of only testing small primes. The correct statement is **quantitative, not
all-or-nothing**: relaxing the exponent cap at a single prime `p` weakens the signature by an amount
that **shrinks with `p`** (because the perturbed density `1/p² − 1/p³` shrinks), and for `p ≥ 11` the
perturbation is small enough that the sieve is "close enough to squarefree" for the detector to score
the signature present again. The open mechanism question (*why* the strict-μ² structure carries a
lag-2 cancellation at all) is unchanged; a new, sharper open question is added: **locate the exact
crossover prime** (it is between 7 and 11 for this detector and threshold) and characterise how `r(p)`
approaches the control value 0.015 as `p → ∞`.

## Reproduce

- `data/axis_d3_large_prime_crossover.json` — the sealed job's full output, bit-identical.
- `axis_d3_large_prime_crossover_job.py` — the job (self-tests the sieve/null/stride/detector, then
  computes), runnable as `python axis_d3_large_prime_crossover_job.py --n-max 20000000 --out-dir .`.
- `test_axis_d3_large_prime_crossover.py` — pure-stdlib verifier: brute-forces the permutation-null
  formula on tiny sequences, re-derives every `E_null`/`var_null`/`z` from the raw `(n_eff, s, S)`,
  re-runs the frozen detector, and asserts the five verdicts, the crossover, and the monotone trend.

## Prior art (why this is new)

Serial digit-sum-parity correlation of *consecutive integers* is classical (Gelfond 1968;
Mauduit–Rivat, Ann. of Math. 171, 2010, 1591–1646; Baake–Coons, arXiv:2209.07102). The nearest
neighbour on binary digits of squarefree integers is Dietmann–Elsholtz–Shparlinski (Trans. AMS 369,
2017, 8369–8388) — about the integers, not the gaps. No prior work measures how the lag-2 gap-parity
cancellation of a "squarefree-except-one-prime" indicator scales with the *size* of the relaxed
prime, nor identifies a prime-size crossover at which a single-prime relaxation stops destroying and
starts preserving the signature. That crossover is the new result.
