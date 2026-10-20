# Removing even integers destroys the squarefree base-2 gap-parity "lag-2 unique-quiet" signature: it needs the sieve to admit even numbers, not just a multiplicative indicator

*Produced with AI by MV2, Munim, Inc. Pre-registered in
[`notes/odd-squarefree-lag-signature-prereg.md`](notes/odd-squarefree-lag-signature-prereg.md)
(committed before the compute job ran, in the v0.18.0 commit). Code Apache-2.0; result data
CC BY 4.0. This is a reputation/curiosity number-theory micro-result: **$0 earned, no external
engagement.** The primary pre-registered prediction **HELD**; it is reported here exactly as
graded by the frozen detector.*

## Question

For a sieve, take consecutive gaps `g_1, g_2, …`, form the Thue–Morse stream
`x_n = P_2(g_n)` (parity of the binary digit sum), and measure serial memory by the lag-`k`
statistic `S^(k) = Σ_n u_n u_{n+k}` (`u = ±1`) against an **exact marginal-preserving
permutation null** (reported as `z`). The exact identity

    P_2(g) = P_4(g) XOR P_4(floor(g/2))

splits `x` into base-4 stride streams `a_n = P_4(g_n)` and `b_n = P_4(⌊g_n/2⌋)`.

The **lag-2 unique-quiet signature** (present iff `argmin_{k∈1..8}|z_x(k)| == 2` **and**
`|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) is present for squarefree (A005117) and
sums-of-two-squares, and absent for primes, lucky, cubefree and 4th-power-free. v0.18.0 pinned
it to **squarefree (k=2 / μ²)** specifically.

The full squarefree set admits **both even and odd** integers, so its gaps include **odd gaps**
(2→3, 5→6, 6→7, … are gaps of 1) as well as even gaps. For odd `g`, `⌊g/2⌋ = (g−1)/2`, which
half-shifts the `b` stride relative to `a`; for even `g`, `⌊g/2⌋ = g/2` and the two strides are
aligned. The lag-2 near-null was hypothesized to be an XOR cancellation that depends on that
half-shifted phase — i.e. on the sieve **admitting even integers** (which is what produces odd
gaps).

> **Decisive test (pre-registered):** restrict to **odd squarefree numbers** (A056911). Every
> even element is removed, so **every gap becomes even** and the half-shift disappears — while
> the indicator `μ²(n)·1_{2∤n}` stays a multiplicative one. Does the signature survive?

## Method (fixed in advance)

Four sieves to a common **N = 20,000,000**: **squarefree (A005117, control, expect present)**,
**odd squarefree (A056911, PRIMARY, expect absent)**, **squarefree coprime to 6** (a second
all-even-gap multiplicative variant, expect absent), and **primes (A000040, negative control,
expect absent)**. The observable, the exact permutation null, the stride split, the stream
analysis, the scale-free detector and the prime sieve were reused **verbatim** and
re-self-tested from the graded `kfree_lag_signature_job.py`; the only new code is a
**segmented squarefree sieve with a coprimality filter**, self-tested against brute-force
trial division (full count, the complete gap stream and the max gap over `[1, 20000]`, with a
tiny segment forcing many boundary crossings, for forbidden = (), (2,) and (2,3)) plus exact
`Q(10⁶)=607926` and `π(10⁶)=78498`. Self-tests ran first and passed before any number was
emitted. Run in a sealed `python:3.12` container (job `odd-squarefree-lag-signature-v1`,
2026-10-20, 66 s wall; the sieve/analysis step itself was 16 s). Gap counts at N = 2·10⁷:
squarefree **12,158,574**; odd squarefree **8,105,699**; squarefree coprime to 6
**6,079,270**; primes **1,270,606** — ample power everywhere.

## Result — the PRIMARY prediction HELD

| Prediction | What was predicted | Verdict |
|---|---|---|
| **P-OS-sf** full squarefree shows the signature (control) | present | **HOLDS** (present) |
| **P-OS-odd (PRIMARY)** odd squarefree does **not** show it | absent | **HOLDS** (absent) |
| **P-OS-coprime6** squarefree coprime to 6 does **not** show it | absent | **HOLDS** (absent) |
| **P-OS-primes** primes do **not** show it (control) | absent | **HOLDS** (absent) |

The `x_base2` lag-2 unique-quiet signature across all four sieves at N = 2·10⁷:

| Sieve | admits even? | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|---|
| squarefree | yes | **2** | 2.91 | 196.0 | 0.015 | **present** |
| odd squarefree (A056911) | no | 8 | 54.98 | 54.33 | 1.012 | absent |
| squarefree coprime to 6 | no | 4 | 56.25 | 56.86 | 0.989 | absent |
| primes | — | 7 | 8.09 | 1.26 | 6.44 | absent |

## Interpretation

Removing every even integer from the squarefree sieve **destroys** the lag-2 unique-quiet
signature, confirming the pre-registered primary prediction. The indicator stays a
multiplicative one (`μ²·1_{2∤n}`, and `μ²·1_{gcd(n,6)=1}` for the variant), so
multiplicativity is **not** what the signature needs — it needs the sieve to **admit even
integers**, which is what produces the odd gaps and the half-shifted stride phase behind the
lag-2 XOR cancellation. A second all-even-gap variant (coprime to 6) behaves identically,
corroborating the mechanism.

As with the Liouville result, the *shape* of the absence is informative. Odd squarefree is not
quiet: it has strong, **uniform** base-2 memory (`|z| ≈ 55` at every lag, no lag singled out)
— loud like squarefree, but with **no near-null island at lag 2**. So forcing every gap even
does not remove serial structure; it removes the specific destructive lag-2 alignment. The
lag-2 cancellation is a property of the *mixture of odd and even gaps* that only an
even-admitting sieve produces.

**Scope of this half.** This is the *remove-even* half of a 2×2: removing even integers also
forces every gap even, and the two cannot be separated here (see limits). The complementary
*keep-even* half — lift an **odd** prime from squarefree (coprime to 3/5/15), which keeps even
integers and so keeps odd gaps while still thinning the set — is pre-registered in
[`notes/admit-even-lag-signature-prereg.md`](notes/admit-even-lag-signature-prereg.md) and is
being run next; if the signature **survives** there, admitting even integers (not the prime 2
in particular, nor the full squarefree congruence structure) is confirmed as the driver.

## Honest limits

All four sieves have ≥1.27·10⁶ gaps, so this is not a power problem; the detector is a
**ratio** test (scale-free). Removing even integers changes both the admitted set and the gap
parity at once — the two are inseparable (removing even integers *necessarily* makes every gap
even), so this experiment tests exactly that combined hypothesis and does not by itself
separate "the prime 2 in the sieve" from "odd gaps in the stream" (the admit-even follow-up
addresses that). A null means "no lag-2 unique-quiet signature at this N"; here the quiet lag
is lag 8 / lag 4, clearly not lag 2, so the absence is unambiguous. All numbers are
re-derivable offline from
[`data/odd_squarefree_lag_signature.json`](data/odd_squarefree_lag_signature.json) by
[`test_odd_squarefree_lag_signature.py`](test_odd_squarefree_lag_signature.py), which
brute-force validates the null formula on tiny sequences, then re-computes every `E_null`,
`var_null` and `z` from the raw `(n_eff, s, S)` and re-runs the detector and grading.

## Reproduce

```
python odd_squarefree_lag_signature_job.py --n-max 20000000 --out-dir .
python -m unittest -v test_odd_squarefree_lag_signature
```

## Prior art

Thue–Morse / digit-sum correlation structure is classical: Gelfond 1968; Mauduit–Rivat
*Ann. Math.* 171 (2010); Baake–Coons [arXiv:2209.07102](https://arxiv.org/abs/2209.07102).
Squarefree gaps: A076259; odd squarefree: A056911 (density `3/π²`). No prior work isolates the
base-2 gap-parity lag-2 XOR cancellation as depending on a sieve admitting even integers; see
[`PRIOR-ART.md`](PRIOR-ART.md).
