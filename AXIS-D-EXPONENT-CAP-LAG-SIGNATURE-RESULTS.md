# Relaxing the squarefree exponent cap at one prime also kills the signature: the base-2 gap-parity "lag-2 unique-quiet" signature needs the *strict* squarefree (μ²) congruence structure

*Produced with AI by MV2, Munim, Inc. Pre-registered in
[`notes/axis-d-exponent-cap-lag-signature-prereg.md`](notes/axis-d-exponent-cap-lag-signature-prereg.md)
(committed before the compute job ran). Code Apache-2.0; result data CC BY 4.0. This is a
reputation/curiosity number-theory micro-result: **$0 earned, no external engagement.** The
PRIMARY pre-registered prediction was **FALSIFIED** — the sharper outcome the pre-registration
flagged as "a real result either way." Everything below is reported exactly as graded by the
frozen detector.*

## Question

For a sieve, take consecutive gaps `g_1, g_2, …`, form the Thue–Morse stream
`x_n = P_2(g_n)` (parity of the binary digit sum), and measure serial memory by the lag-`k`
statistic `S^(k) = Σ_n u_n u_{n+k}` (`u = ±1`) against an **exact marginal-preserving
permutation null** (reported as `z`). The exact identity

    P_2(g) = P_4(g) XOR P_4(floor(g/2))

splits `x` into base-4 stride streams `a_n = P_4(g_n)` and `b_n = P_4(⌊g_n/2⌋)`.

The **lag-2 unique-quiet signature** (present iff `argmin_{k∈1..8}|z_x(k)| == 2` **and**
`|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) was pinned to **squarefree (A005117)** at v0.18.0.
Axis C established that the signature needs the full squarefree congruence structure from the
*thinning* side: removing even integers (odd squarefree, v0.21.0) **and** keeping evens while
lifting one odd prime (squarefree coprime to 3/5/15, v0.23.0) both destroy it. So
even-admission is necessary but not sufficient.

Axis D asks the **complementary, enlarging** question Axis C cannot:

> **Decisive test (pre-registered):** instead of removing residue classes, *add some back*.
> Relax the μ² exponent cap at **one** small prime while keeping every other prime squarefree:
> **squarefree-except-2** (`A_2 = {0,1,2}` — admit 4·odd like 4, 12, 20, 36 … but not 8|n;
> PRIMARY) and **squarefree-except-3** (`A_3 = {0,1,2}` — admit 9, 18, 45 … but not 27|n;
> replication at an *odd* prime). Both are multiplicative 0/1 indicators that **still admit
> every even integer**, so odd gaps and the half-shifted base-4 stride phase survive. On the
> "even-admission is the driver" reading these should keep the signature. Odd squarefree
> (A056911) is the even-removing **anchor** (expected absent); primes (A000040) the **negative
> control** (expected absent).

## Method (fixed in advance)

Five sieves to a common **N = 20,000,000**, each analysed with the **identical** observable,
exact permutation null, stride split and scale-free detector used across the arc. Self-tests
ran first and passed before any number was emitted. Run in a sealed `python:3.12` container
(job `axis-d-exponent-cap-lag-signature-v1`, 2026-10-20, 66 s wall; the sieve/analysis step
itself 20 s). Gap counts at N = 2·10⁷ — squarefree **12,158,574**; squarefree-except-2
**14,184,990**; squarefree-except-3 **13,171,778**; odd squarefree **8,105,699**; primes
**1,270,606** — ample power everywhere. Natural densities (prereg): squarefree `6/π² ≈ 0.6079`,
squarefree-except-2 `7/π² ≈ 0.7094`, squarefree-except-3 `13/(2π²) ≈ 0.6586`.

## Result — the PRIMARY prediction was FALSIFIED

| Prediction | Keeps even integers? | What was predicted | Verdict |
|---|---|---|---|
| **P-D-sf** full squarefree (control) | yes | present | **HELD** (present) |
| **P-D-ex2 (PRIMARY)** squarefree-except-2 | yes | present | **FALSIFIED** (absent) |
| **P-D-ex3** squarefree-except-3 | yes | present | **FALSIFIED** (absent) |
| **P-D-odd (anchor)** odd squarefree (A056911) | no | absent | **HELD** (absent) |
| **P-D-primes (neg control)** primes (A000040) | — | absent | **HELD** (absent) |

The `x_base2` lag-2 unique-quiet detector across all five sieves at N = 2·10⁷:

| Sieve | admits even? | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|---|
| squarefree | yes | **2** | 2.91 | 196.0 | 0.015 | **present** |
| squarefree-except-2 | yes | 5 | 189.20 | 165.63 | 1.142 | absent |
| squarefree-except-3 | yes | 8 | 57.59 | 74.55 | 0.773 | absent |
| odd squarefree (A056911) | no | 8 | 54.98 | 54.33 | 1.012 | absent |
| primes (A000040) | — | 7 | 8.09 | 1.26 | 6.444 | absent |

## Interpretation

Relaxing the squarefree exponent cap at even a **single** prime — admitting `4·odd`
(squarefree-except-2) or `9·(coprime-3)` (squarefree-except-3) — **destroys** the lag-2
unique-quiet signature. The quiet lag jumps to 5 and 8, and the lag-2 `|z|` is of the same
order as the other lags (ratio 1.14 and 0.77, far from the 0.1 threshold). Both relaxations
**keep every even integer** and still produce odd gaps, so this fails the "even-admission is
sufficient" reading from the *enlarging* side — exactly as Axis C failed it from the *thinning*
side.

Put together with Axis C, this gives a clean, symmetric conclusion:

- **The signature tolerates neither thinning nor enlarging.** Remove one prime's residues
  (Axis C: coprime 3/5/15, odd squarefree) → gone. Add back one prime's square class (Axis D:
  squarefree-except-2/-3) → gone. Only the **exact strict squarefree (μ²) sieve** — every prime
  exponent ≤ 1, no prime's class removed — carries the lag-2 unique-quiet cancellation.
- **Even-admission is necessary but not sufficient, confirmed from both sides.** The anchor
  (odd squarefree, no evens) is absent, so admitting evens is necessary; but every
  even-admitting perturbation tested (Axis C coprime-k; Axis D except-2/-3) is also absent, so
  it is not sufficient. The decisive grading flag
  `signature_tolerates_exponent_cap_relaxation` is **False**.
- The *shape* of each absence is informative, as before. The exponent-cap relaxations are
  **loud and non-uniform** (large, lag-dependent `|z|`, no near-null island at lag 2): loosening
  the cap does not quiet the serial structure, it removes the specific destructive lag-2
  alignment. Only the complete μ² sieve produces that near-perfect lag-2 cancellation.

## What this closes

This closes **Axis D** and, with it, the four-axis map of what the signature depends on:

- **Axis A (which k-free):** only k = 2 — squarefree, not cube-free or 4th-power-free (v0.18.0).
- **Axis B (multiplicative origin vs. indicator):** needs a multiplicative *indicator*; the
  Liouville +1 level set (multiplicative origin, non-multiplicative indicator) is flat (v0.20.0).
- **Axis C (thinning):** removing evens or lifting one odd prime both kill it (v0.21.0 / v0.23.0).
- **Axis D (enlarging):** relaxing the exponent cap at one prime kills it too (this result).

The signature localizes to **exactly the strict squarefree congruence structure** — a
multiplicative indicator, even-admission, *and* strict μ² at every prime, none of the three
alone and no proper multiplicative relative tested reproducing it.

## Honest limits

All five sieves have ≥1.27·10⁶ gaps (the four multiplicative ones ≥8.1·10⁶), so this is not a
power problem; the detector is a **ratio** test (scale-free). A null means "no lag-2
unique-quiet signature at this N"; here the quiet lag is 5 / 8 / 8 / 7, clearly not lag 2, so
each absence is unambiguous. We relaxed the cap at two specific small primes (2 and 3) to
exponent 2; we do not claim *every* exponent-cap relaxation at *every* prime behaves this way,
only that no multiplicative relative tested across Axes A–D reproduces the signature. The
conclusion is about the N = 2·10⁷ regime; the arc makes no asymptotic claim. The question is
non-tautological only at base 2/4 (see `notes/base2-even-admission-nontautology.md`), so it is
graded in base 2 exclusively. All numbers are re-derivable offline from
[`data/axis_d_exponent_cap_lag_signature.json`](data/axis_d_exponent_cap_lag_signature.json) by
[`test_axis_d_exponent_cap_lag_signature.py`](test_axis_d_exponent_cap_lag_signature.py), which
brute-force validates the null formula on tiny sequences, then re-computes every `E_null`,
`var_null` and `z` from the raw `(n_eff, s, S)` and re-runs the detector and grading.

## Reproduce

```
python axis_d_exponent_cap_lag_signature_job.py --n-max 20000000 --out-dir .
python -m unittest -v test_axis_d_exponent_cap_lag_signature
```

## Prior art

Thue–Morse / digit-sum correlation structure is classical: Gelfond 1968; Mauduit–Rivat
*Ann. Math.* 171 (2010), 1591–1646; Baake–Coons [arXiv:2209.07102](https://arxiv.org/abs/2209.07102).
Squarefree gaps: A076259; squarefree density `6/π²` (A005117); odd squarefree `3/π²` (A056911).
No prior work isolates the base-2 gap-parity lag-2 XOR cancellation, nor shows it depends on the
complete, strict squarefree congruence structure rather than on even-admission or a multiplicative
indicator alone; see [`PRIOR-ART.md`](PRIOR-ART.md).
