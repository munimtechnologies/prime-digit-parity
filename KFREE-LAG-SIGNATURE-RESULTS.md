# The base-2 gap-parity "lag-2 unique-quiet" signature is specific to squarefree (k=2), not a k-free family property

*Produced with AI by MV2, Munim, Inc. Pre-registered in
[`notes/kfree-lag-signature-prereg.md`](notes/kfree-lag-signature-prereg.md) (committed
before the compute job ran, in the v0.17.0 commit). Code Apache-2.0; result data CC BY 4.0.
This is a reputation/curiosity number-theory micro-result: **$0 earned, no external
engagement.** The primary pre-registered prediction was **falsified**; it is reported here
exactly as such.*

## Question

Across this arc, the base-2 gap digit-sum-parity stream of some integer sieves carries a
sharp fingerprint. For a sieve, take consecutive gaps `g_1, g_2, …`, form the Thue–Morse
stream `x_n = P_2(g_n)` (parity of the binary digit sum = popcount mod 2), and measure
serial memory by the lag-`k` statistic `S^(k) = Σ_n u_n u_{n+k}` (`u = ±1`) against an
**exact marginal-preserving permutation null** (reported as `z`). The exact identity

    P_2(g) = P_4(g) XOR P_4(floor(g/2))

splits `x` into two base-4 "stride" streams `a_n = P_4(g_n)` and `b_n = P_4(⌊g_n/2⌋)`.

The **lag-2 unique-quiet signature** (scale-free detector: present iff
`argmin_{k∈1..8}|z_x(k)| == 2` **and** `|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) was found
**present** for squarefree (A005117) and sums-of-two-squares (A001481), and **absent** for
primes (A000040) and the positional lucky sieve (A000959, graded v0.17.0). The reading so
far: the lag-2 XOR cancellation tracks μ²-type multiplicative density structure, not
thinning.

Squarefree is the **k = 2** member of the **k-free family**: `n` is *k-free* iff no prime
power `p^k` divides `n`, with multiplicative indicator `∏_p [v_p(n) < k]` for every `k ≥ 2`.

> **Decisive test (pre-registered):** is the signature a property of the whole **k-free
> family** (so **cubefree** k=3 and **4th-power-free** k=4 also show it → family property),
> or is it **specific to squarefree (k=2)**?

## Method (fixed in advance)

Four sieves to a common **N = 20,000,000**: squarefree (k=2, control), **cubefree (k=3,
primary)**, 4th-power-free (k=4), and primes (negative control). The k-free densities are
`1/ζ(k)` ≈ 0.608, 0.832, 0.924, so each gives >1.2·10⁷ gaps — ample power. The observable,
the exact permutation null `null_z`, the stride split, the stream analysis and the
scale-free signature detector were reused **verbatim** and re-self-tested from the graded
`lucky_lag_signature_job.py`; the only new code is a generalized segmented k-free sieve,
self-tested against exact small counts (k-free in `[1,100]` = 61 / 85 / 93; squarefree
count `Q(10⁶) = 607926`) and `π(10⁶) = 78498`. The job's seven self-tests ran first and
passed before any number was emitted. Run in a sealed `python:3.12` container (job
`kfree-lag-signature-v1`, 2026-10-19, 104 s wall; the sieve/analysis step itself was 17 s).
Gap counts at N = 2·10⁷: primes 1,270,606; squarefree 12,158,574; cubefree 16,638,148;
4th-power-free 18,478,767.

## Result — the PRIMARY prediction was FALSIFIED

| Prediction | What was predicted | Verdict |
|---|---|---|
| **P-K-2** squarefree shows the signature (control) | present | **HOLDS** (present) |
| **P-K-3 (PRIMARY)** cubefree shows the signature | present | **FALSIFIED** (absent) |
| **P-K-4** 4th-power-free shows the signature | present | **FALSIFIED** (absent) |
| **P-K-primes** primes do **not** show the signature (control) | absent | **HOLDS** (absent) |

The `x_base2` lag-2 unique-quiet signature across all four sieves at N = 2·10⁷:

| Sieve | k | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|---|
| squarefree | 2 | **2** | 2.91 | 196.0 | 0.015 | **present** |
| cubefree | 3 | 5 | 59.65 | 59.18 | 1.008 | absent |
| 4th-power-free | 4 | 3 | 8.556 | 8.557 | 1.000 | absent |
| primes | — | 7 | 8.09 | 1.26 | 6.44 | absent |

## Interpretation

The lag-2 unique-quiet cancellation does **not** survive raising `k`. Only squarefree
(k = 2) shows it; cubefree and 4th-power-free do not. The signature is therefore **specific
to squarefree**, not a property of the k-free family — exactly the sharper outcome the
pre-registration named as the falsification branch. Combined with the earlier controls
(sums-of-two-squares present; primes and lucky absent), the signature is now pinned to the
**k = 2 / μ² structure specifically**, not to "multiplicative indicator density" in general:
cubefree and 4th-power-free have multiplicative indicators of even higher density and still
lack it.

The two higher-`k` sieves are not quiet — they carry *different* structure:

- **Cubefree** has a strong, uniform base-2 memory (`|z| ≈ 59` at every lag, no lag singled
  out), and its **base-4 strides show a sharp loud peak at lag 7** (`a`: `z ≈ +1787`; `b`:
  `z ≈ +1872`), far louder than squarefree's strides. So cubefree has *more* serial
  structure than squarefree, yet the base-2 XOR produces **no near-null at lag 2** — the
  destructive lag-2 alignment is a squarefree-specific coincidence, not a generic k-free one.
- **4th-power-free** is uniformly weakly loud in base 2 (`|z| ≈ 8.56` at every lag) with no
  distinguished lag; its base-4 strides are uniformly loud (`|z| ≈ 247`–`259`).

So the picture is: the base-2 lag-2 unique-quiet XOR cancellation is a fingerprint of the
**squarefree (and sums-of-two-squares) density structure in particular**, and it degrades
into ordinary (non-cancelling) serial memory as soon as `k ≥ 3`.

## Honest limits

All four sieves here have ≥1.27·10⁶ gaps, so this is not a power problem; the signature test
is a **ratio** test (scale-free) and is not invalidated by finite N. A null for P-K-3/P-K-4
means "no lag-2 unique-quiet signature at this N", and at this N the quiet lag is clearly
*not* lag 2 (it is lag 5 and lag 3 respectively, each far from a null), so the falsification
is unambiguous rather than marginal. The claim is that the signature tracks the specific
k = 2 (squarefree / SoTS) structure, evidenced across k = 2, 3, 4 plus the prime and lucky
controls — not proven for all multiplicative sieves. All numbers are re-derivable offline
from `data/kfree_lag_signature.json` by
[`test_kfree_lag_signature.py`](test_kfree_lag_signature.py), which brute-force-validates
the null formula on tiny sequences, then re-computes every `E_null`, `var_null` and `z` from
the raw `(n_eff, s, S)` and re-runs the detector and grading.

## Reproduce

```
python kfree_lag_signature_job.py --n-max 20000000 --out-dir .
python -m unittest -v test_kfree_lag_signature
```

## Prior art

Thue–Morse / digit-sum correlation structure is classical: Gelfond 1968; Mauduit–Rivat
*Ann. Math.* 171 (2010); Baake–Coons [arXiv:2209.07102](https://arxiv.org/abs/2209.07102);
Mahler 1927. Squarefree digit statistics: Dietmann–Elsholtz–Shparlinski
[arXiv:1601.04754](https://arxiv.org/abs/1601.04754). k-free densities `1/ζ(k)` are
classical (Gegenbauer). No prior work isolates the base-2 gap-parity lag-2 XOR cancellation
as squarefree-specific against the k-free family; see [`PRIOR-ART.md`](PRIOR-ART.md).
