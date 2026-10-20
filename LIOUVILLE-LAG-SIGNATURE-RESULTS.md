# The base-2 gap-parity "lag-2 unique-quiet" signature needs a multiplicative INDICATOR, not merely a multiplicative origin: the Liouville +1 level set does not carry it

*Produced with AI by MV2, Munim, Inc. Pre-registered in
[`notes/liouville-lag-signature-prereg.md`](notes/liouville-lag-signature-prereg.md)
(committed before the compute job ran, in the v0.17.0 commit). Code Apache-2.0; result data
CC BY 4.0. This is a reputation/curiosity number-theory micro-result: **$0 earned, no
external engagement.** The primary pre-registered prediction **HELD**; it is reported here
exactly as graded by the frozen detector.*

## Question

Across this arc, the base-2 gap digit-sum-parity stream of some integer sieves carries a
sharp fingerprint. For a sieve, take consecutive gaps `g_1, g_2, …`, form the Thue–Morse
stream `x_n = P_2(g_n)` (parity of the binary digit sum = popcount mod 2), and measure serial
memory by the lag-`k` statistic `S^(k) = Σ_n u_n u_{n+k}` (`u = ±1`) against an **exact
marginal-preserving permutation null** (reported as `z`). The exact identity

    P_2(g) = P_4(g) XOR P_4(floor(g/2))

splits `x` into two base-4 "stride" streams `a_n = P_4(g_n)` and `b_n = P_4(⌊g_n/2⌋)`.

The **lag-2 unique-quiet signature** (scale-free detector: present iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) was found
**present** for squarefree (A005117) and sums-of-two-squares (A001481) — whose *set
indicators are multiplicative* — and **absent** for primes (A000040), the lucky sieve
(A000959, v0.17.0), and cubefree / 4th-power-free (v0.18.0). So far the signature tracks the
squarefree / μ²-type multiplicative-indicator structure.

The **Liouville function** `λ(n) = (−1)^{Ω(n)}` is *completely multiplicative*. Its **+1 level
set** `A028260 = { n : Ω(n) even } = {1, 4, 6, 9, 10, 14, 15, 21, 22, 25, …}` has density 1/2
and a clean multiplicative *origin* (it is a level set of a completely multiplicative
function), but its **set indicator `1_{Ω even}` is NOT multiplicative** (e.g. `m=n=2`:
`1(4)=1 ≠ 1(2)·1(2)=0`). This is exactly the discriminator the earlier results could not
separate.

> **Decisive test (pre-registered):** does the signature require a **multiplicative
> indicator** (as squarefree / SoTS have), or does it attach to **multiplicative-function
> structure** more broadly — appearing even for the Liouville +1 level set, whose indicator
> is not multiplicative?

## Method (fixed in advance)

Three sieves to a common **N = 20,000,000**: **Liouville +1 (A028260, primary)**, squarefree
(A005117, multiplicative-indicator control, expect present), primes (A000040, negative
control, expect absent). The observable, the exact permutation null `null_z`, the stride
split, the stream analysis, the scale-free detector, and the squarefree/prime sieves were
reused **verbatim** and re-self-tested from the graded `kfree_lag_signature_job.py`; the only
new code is a **segmented exact Ω-parity (Liouville) sieve**, self-tested against brute-force
trial-division Ω across many forced segment boundaries (full count *and* the complete gap
stream over `[1, 20000]`) plus the exact small set `A028260 ≤ 10 = {1,4,6,9,10}`. The job's
self-tests ran first and passed before any number was emitted. Run in a sealed `python:3.12`
container (job `liouville-lag-signature-v1`, 2026-10-20, 66 s wall; the sieve/analysis step
itself was 23 s). Gap counts at N = 2·10⁷: Liouville +1 **9,997,744**; squarefree
**12,158,574**; primes **1,270,606** — ample power everywhere.

## Result — the PRIMARY prediction HELD

| Prediction | What was predicted | Verdict |
|---|---|---|
| **P-L-sf** squarefree shows the signature (control) | present | **HOLDS** (present) |
| **P-L-primes** primes do **not** show it (control) | absent | **HOLDS** (absent) |
| **P-L-liouville (PRIMARY)** Liouville +1 does **not** show it | absent | **HOLDS** (absent) |

The `x_base2` lag-2 unique-quiet signature across all three sieves at N = 2·10⁷:

| Sieve | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|
| squarefree | **2** | 2.91 | 196.0 | 0.015 | **present** |
| Liouville +1 (A028260) | 1 | 1.34 | 0.565 | 2.38 | absent |
| primes | 7 | 8.09 | 1.26 | 6.44 | absent |

## Interpretation

The Liouville +1 level set does **not** carry the lag-2 unique-quiet signature, confirming
the pre-registered primary prediction: the signature needs a **multiplicative indicator**, not
merely a multiplicative *origin*. A clean multiplicative origin (being a level set of a
completely multiplicative function) is **not sufficient**.

The *way* it is absent is itself informative. Squarefree has enormous, loud base-2 serial
memory (`|z| ≈ 150`–`225` at every lag) with a single destructive near-null **island** at lag
2 — that island is the signature. The **Liouville +1 stream has essentially no base-2 serial
memory at all**: every lag sits at `|z| ≤ 1.34`, indistinguishable from the permutation null.
So squarefree and Liouville +1 fail the detector for opposite reasons — squarefree is loud
everywhere except a lag-2 null; Liouville +1 is quiet everywhere. The lag-2 XOR cancellation
is a feature of the *loud, structured* gap-parity memory that the multiplicative squarefree
indicator induces, and the Ω-parity level set simply does not induce that memory.

Combined with the earlier results, the signature is now pinned to **multiplicative set
indicators of the μ² type (squarefree / SoTS)** specifically: it is absent for the k-free
family at `k ≥ 3` (v0.18.0), for positional sieves (lucky, v0.17.0), for primes, and now for
a multiplicative-origin level set with a non-multiplicative indicator.

## Honest limits

All three sieves have ≥1.27·10⁶ gaps, so this is not a power problem; the detector is a
**ratio** test (scale-free) and is not invalidated by finite N. A null means "no lag-2
unique-quiet signature at this N"; for Liouville +1 the whole z-profile is flat (`|z| ≤ 1.34`
at every lag), so the absence is unambiguous rather than marginal. A single level set
(`λ = +1`) is one object; this result constrains but does not by itself settle the general
multiplicative-indicator-vs-origin question. All numbers are re-derivable offline from
[`data/liouville_lag_signature.json`](data/liouville_lag_signature.json) by
[`test_liouville_lag_signature.py`](test_liouville_lag_signature.py), which brute-force
validates the null formula on tiny sequences, then re-computes every `E_null`, `var_null` and
`z` from the raw `(n_eff, s, S)` and re-runs the detector and grading.

## Reproduce

```
python liouville_lag_signature_job.py --n-max 20000000 --out-dir .
python -m unittest -v test_liouville_lag_signature
```

## Prior art

Thue–Morse / digit-sum correlation structure is classical: Gelfond 1968; Mauduit–Rivat
*Ann. Math.* 171 (2010); Baake–Coons [arXiv:2209.07102](https://arxiv.org/abs/2209.07102);
Mahler 1927. The Liouville function and `A028260` (`Ω` even) are classical; Chowla-type
correlation conjectures concern `λ` directly, not the base-2 digit-sum parity of gaps between
consecutive members of its +1 level set. No prior work tests whether the base-2 gap-parity
lag-2 XOR cancellation distinguishes a multiplicative *indicator* from a multiplicative
*origin*; see [`PRIOR-ART.md`](PRIOR-ART.md).
