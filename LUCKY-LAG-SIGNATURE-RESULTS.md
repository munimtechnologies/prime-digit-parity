# The base-2 gap-parity "lag-2 unique-quiet" signature is tied to multiplicative sieves, not to thinning

*Produced with AI by MV2, Munim, Inc. Pre-registered in
[`notes/lucky-lag-signature-prereg.md`](notes/lucky-lag-signature-prereg.md) (committed
before the compute job ran). Code Apache-2.0; result data CC BY 4.0. This is a
reputation/curiosity number-theory micro-result: **$0 earned, no external engagement.**
A negative or partial result is reported as such.*

## Question

Earlier runs in this arc found a sharp fingerprint on the base-2 gap digit-sum-parity
stream of some integer sieves. For a sieve, take consecutive gaps `g_1, g_2, …`, form the
Thue–Morse stream `x_n = P_2(g_n)` (parity of the binary digit sum = popcount mod 2), and
measure serial memory by the lag-`k` statistic `S^(k) = Σ_n u_n u_{n+k}` (`u = ±1`) against
an **exact marginal-preserving permutation null** (reported as `z`). The exact identity

    P_2(g) = P_4(g) XOR P_4(floor(g/2))

splits `x` into two base-4 "stride" streams `a_n = P_4(g_n)` and `b_n = P_4(⌊g_n/2⌋)`.

At N = 10⁸, **squarefree** and **sums of two squares** both show a *lag-2 unique-quiet
signature*: of lags 1..8, lag 2 is the single quiet lag on the `x_base2` stream even though
both base-4 strides are individually loud there — a destructive XOR cancellation that fires
only at lag 2. **Primes do not** show it. Squarefree and sums-of-two-squares are
multiplicative sieves.

> **Decisive test (pre-registered):** does the **lucky-number sieve** (OEIS
> [A000959](https://oeis.org/A000959)) — a purely *positional, Josephus-style* sieve with
> **no** multiplicative/congruence input, yet prime-like density — show the signature
> (→ the cancellation is generic to any thinned sequence) or not (→ it is tied to
> multiplicative structure)?

## Method (fixed in advance)

All four sieves computed to a common **N = 20,000,000** (chosen for lucky-sieve
tractability; it is O(N × passes)). The observable, the exact permutation null `null_z`, and
the stream split were reused verbatim and re-self-tested from the graded
`squarefree_lag_period_job.py`; the lucky sieve was reused and re-self-tested from the
graded `lucky_gap_parity.py`. The job's seven self-tests ran first and passed before any
number was emitted. The signature detector is **scale-free**: present iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`.

Run in a sealed `python:3.12` container (job `lucky-lag-signature-v1`, 2026-10-19,
146 s compute). Gap counts at N = 2·10⁷: primes 1,270,606; squarefree 12,158,574;
sums-of-two-squares 3,880,118; lucky 1,164,977.

## Result — all three pre-registered predictions HELD

| Prediction | Statistic | Verdict |
|---|---|---|
| **P-L-1** lucky base-2 lag-1 anti-persistence (`z ≤ −8`) | `z = −22.45` | **HOLDS** |
| **P-L-2 (PRIMARY)** lucky does **not** show the lag-2 signature | `signature_present = false` | **HOLDS** |
| **P-L-3** both lucky base-4 strides loud at lag 1 (`|z| ≥ 5`) | `a: |z|=8.46`, `b: |z|=38.48` | **HOLDS** |

The `x_base2` lag-2 unique-quiet signature across all four sieves at N = 2·10⁷:

| Sieve | type | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|---|
| squarefree | multiplicative (μ²) | **2** | 2.91 | 196.0 | 0.015 | **present** |
| sums of two squares | multiplicative | **2** | 0.73 | 15.2 | 0.048 | **present** |
| primes | multiplicative | 7 | 8.09 | 1.26 | 6.44 | absent |
| **lucky** | **positional (Josephus)** | 7 | 4.48 | 0.71 | 6.31 | **absent** |

## Interpretation

The lucky sieve, of matched (prime-like) density but built purely by *position*, does
**not** show the lag-2 unique-quiet cancellation — it patterns with **primes**, not with
squarefree / sums-of-two-squares. Yet lucky is not structureless: its base-2 stream is
strongly anti-persistent at lag 1 (`z = −22.45`) and both its base-4 strides are
individually loud (`|z| = 8.46` and `38.48`). So lucky has loud strides **and** serial
memory, but **no lag-2 XOR cancellation**.

That rules out the "generic to any thinned sequence" hypothesis: having loud base-4 strides
is not sufficient to produce the lag-2 near-null. The signature appears only for the two
**μ²-type multiplicative density sieves** (squarefree and sums-of-two-squares), where the
lag-2 alignment of the two base-4 residue structures cancels. It is absent both for primes
(multiplicative but not an indicator-density sieve of this kind) and for the positional
lucky sieve. The internal replication reproduced the 1e8 pattern at this smaller N
(squarefree & sums-of-two-squares present; primes absent), confirming the reading is not an
N artefact.

## Honest limits

Lucky at N = 2·10⁷ yields ~1.16M gaps (vs squarefree's ~12.2M), so lucky |z| values are a
few-fold smaller. The signature test is a **ratio** test and is not invalidated by the
smaller N, but a null for P-L-2 means "no signature at this power", not a proof of its
impossibility at all N. The three multiplicative sieves and the one positional sieve are a
small sample of "kinds of sieve"; the claim is that the signature tracks μ²-type
multiplicative density structure, evidenced here but not exhaustively. All numbers are
re-derivable offline from `data/lucky_lag_signature.json` by
[`test_lucky_lag_signature.py`](test_lucky_lag_signature.py), which re-computes every null
moment and `z` from the raw `(n_eff, s, S)` and re-runs the detector and grading.

## Reproduce

```
python lucky_lag_signature_job.py --n-max 20000000 --out-dir .
python -m unittest -v test_lucky_lag_signature
```

## Prior art

Thue–Morse / digit-sum correlation structure is classical: Gelfond 1968; Mauduit–Rivat
*Ann. Math.* 171 (2010); Baake–Coons [arXiv:2209.07102](https://arxiv.org/abs/2209.07102);
Mahler 1927. Squarefree/sum-of-two-squares digit statistics: Dietmann–Elsholtz–Shparlinski
[arXiv:1601.04754](https://arxiv.org/abs/1601.04754). The lucky sieve: Gardiner–Lazarus–
Metropolis–Ulam 1956; OEIS A000959. No prior work compares the base-2 gap-parity lag-2
cancellation across multiplicative vs positional sieves; see
[`PRIOR-ART.md`](PRIOR-ART.md).
