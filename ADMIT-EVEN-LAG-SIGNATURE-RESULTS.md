# Admitting even integers is necessary but NOT sufficient: the squarefree base-2 gap-parity "lag-2 unique-quiet" signature needs the full squarefree congruence structure

*Produced with AI by MV2, Munim, Inc. Pre-registered in
[`notes/admit-even-lag-signature-prereg.md`](notes/admit-even-lag-signature-prereg.md)
(committed before the compute job ran). Code Apache-2.0; result data CC BY 4.0. This is a
reputation/curiosity number-theory micro-result: **$0 earned, no external engagement.** The
PRIMARY pre-registered prediction was **FALSIFIED** — the sharper outcome the pre-registration
flagged as "equally publishable." Everything below is reported exactly as graded by the frozen
detector.*

## Question

For a sieve, take consecutive gaps `g_1, g_2, …`, form the Thue–Morse stream
`x_n = P_2(g_n)` (parity of the binary digit sum), and measure serial memory by the lag-`k`
statistic `S^(k) = Σ_n u_n u_{n+k}` (`u = ±1`) against an **exact marginal-preserving
permutation null** (reported as `z`). The exact identity

    P_2(g) = P_4(g) XOR P_4(floor(g/2))

splits `x` into base-4 stride streams `a_n = P_4(g_n)` and `b_n = P_4(⌊g_n/2⌋)`.

The **lag-2 unique-quiet signature** (present iff `argmin_{k∈1..8}|z_x(k)| == 2` **and**
`|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) was pinned to **squarefree (A005117)** at v0.18.0.
The companion *remove-even* study (v0.21.0, odd squarefree A056911) found it **absent** once
every even integer is removed, and read that as: the signature needs the sieve to **admit even
integers** (which produces the odd gaps and the half-shifted base-4 stride phase behind the
lag-2 XOR cancellation), not just a multiplicative indicator.

That experiment could not separate two things: removing the even integers also forces **every
gap even**. This study fills the complementary cell.

> **Decisive test (pre-registered):** lift an **odd** prime from the squarefree sieve —
> squarefree **coprime to 3** (PRIMARY), and to **5** and to **15** as replications. These
> **keep every even integer**, so odd gaps and the half-shifted stride phase survive, yet the
> sieve is no longer the full squarefree set. Prediction: the signature **survives**. Odd
> squarefree (A056911) is included as the even-removing **anchor** (expected absent). If the
> odd-prime lifts keep the signature while the anchor loses it, "admitting even integers" is
> confirmed as the driver.

## Method (fixed in advance)

Five sieves to a common **N = 20,000,000**, each analysed with the **identical** observable,
exact permutation null, stride split, scale-free detector and segmented squarefree-coprime
sieve reused **verbatim** (and re-self-tested) from `odd_squarefree_lag_signature_job.py`. The
self-test additionally asserts the mechanism contrast: every odd-prime-lift sieve contains even
members **and at least one odd gap**, whereas the prime-2 lift (odd squarefree) has all even
gaps; plus exact `Q(10⁶)=607926` and `π(10⁶)=78498`. Self-tests ran first and passed before any
number was emitted. Run in a sealed `python:3.12` container (job
`admit-even-lag-signature-v1`, 2026-10-20, 90 s wall; the sieve/analysis step itself 27 s).
Gap counts at N = 2·10⁷: squarefree **12,158,574**; coprime-3 **9,118,929**; coprime-5
**10,132,122**; coprime-15 **7,599,091**; odd squarefree **8,105,699** — ample power everywhere.

## Result — the PRIMARY prediction was FALSIFIED

| Prediction | Keeps even integers? | What was predicted | Verdict |
|---|---|---|---|
| **P-AE-sf** full squarefree (control) | yes | present | **HELD** (present) |
| **P-AE-c3 (PRIMARY)** squarefree coprime to 3 | yes | present | **FALSIFIED** (absent) |
| **P-AE-c5** squarefree coprime to 5 | yes | present | **FALSIFIED** (absent) |
| **P-AE-c15** squarefree coprime to 15 | yes | present | **FALSIFIED** (absent) |
| **P-AE-odd (anchor)** odd squarefree (A056911) | no | absent | **HELD** (absent) |

The `x_base2` lag-2 unique-quiet detector across all five sieves at N = 2·10⁷:

| Sieve | admits even? | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|---|
| squarefree | yes | **2** | 2.91 | 196.0 | 0.015 | **present** |
| squarefree coprime 3 | yes | 8 | 352.04 | 427.08 | 0.824 | absent |
| squarefree coprime 5 | yes | 4 | 253.89 | 274.63 | 0.925 | absent |
| squarefree coprime 15 | yes | 3 | 308.64 | 445.32 | 0.693 | absent |
| odd squarefree (A056911) | no | 8 | 54.98 | 54.33 | 1.012 | absent |

## Interpretation

Keeping every even integer while removing a **single odd prime's** multiples (coprime to 3, or
5, or 15) **destroys** the lag-2 unique-quiet signature — the quiet lag jumps to 8 / 4 / 3 and
the lag-2 `|z|` is of the same order as the other lags (ratio 0.69–0.92, nowhere near the 0.1
threshold). So the prediction failed, and it failed the same way in all three independent lifts.

This is **sharper** than the hypothesis it refutes:

- **Admitting even integers is necessary but NOT sufficient.** The odd-prime lifts keep the
  even integers, the odd gaps and the half-shifted stride phase (the self-test proves each has
  even members and odd gaps), yet the signature is gone. So the v0.21.0 reading — "the driver
  is admitting even integers" — is **incomplete**: even-admission alone does not buy the
  signature.
- **The signature needs the full squarefree congruence structure.** It appears on A005117 and
  on *none* of the proper multiplicative sub-sieves tested — whether they admit even integers
  (coprime 3/5/15) or not (odd squarefree, and coprime-6 from v0.21.0). Removing any one prime's
  square-class residues is enough to break it. This **strengthens** the v0.18.0 result that the
  effect is specific to squarefree: it is specific to the *complete* squarefree sieve, not
  inherited by any thinned multiplicative relative we have tried.
- The decisive grading flag `signature_tracks_even_admission_not_full_mu2` is **False**: the
  odd-prime lifts did not keep the signature, so the data reject the even-admission account and
  point to the full-congruence account.

As in the earlier halves, the *shape* of each absence is informative. The odd-prime lifts are
**loud and non-uniform** (large, lag-dependent `|z|`, no near-null island at lag 2) — the
thinning does not quiet the serial structure, it removes the specific destructive lag-2
alignment. Only the full squarefree set produces that near-perfect lag-2 cancellation.

## What this closes, and what it opens

This closes Axis C of the arc as a clean 2×2: *remove-even* (v0.21.0, absent) and *keep-even,
lift-odd-prime* (this study, absent) both kill the signature, so neither even-admission nor a
multiplicative indicator is sufficient — the full μ² congruence structure is required. It
sharpens the open question for **Axis D** (allow a squared factor only at the small primes 2
and 3, pre-registered in
[`notes/axis-d-exponent-cap-lag-signature-prereg.md`](notes/axis-d-exponent-cap-lag-signature-prereg.md)):
given how fragile the signature is to *removing* one prime's residues, does *adding back* the
`4|n` / `9|n` integers destroy it too? That is the natural next perturbation.

## Honest limits

All five sieves have ≥7.6·10⁶ gaps, so this is not a power problem; the detector is a **ratio**
test (scale-free). A null means "no lag-2 unique-quiet signature at this N"; here the quiet lag
is 8 / 4 / 3, clearly not lag 2, so each absence is unambiguous. We tested three specific odd
lifts (3, 5, 15); we do not claim *every* proper sub-sieve behaves this way, only that no
multiplicative relative tested across Axes A–C reproduces the signature. The conclusion is
about the N = 2·10⁷ regime; the arc makes no asymptotic claim. All numbers are re-derivable
offline from [`data/admit_even_lag_signature.json`](data/admit_even_lag_signature.json) by
[`test_admit_even_lag_signature.py`](test_admit_even_lag_signature.py), which brute-force
validates the null formula on tiny sequences, then re-computes every `E_null`, `var_null` and
`z` from the raw `(n_eff, s, S)` and re-runs the detector and grading.

## Reproduce

```
python admit_even_lag_signature_job.py --n-max 20000000 --out-dir .
python -m unittest -v test_admit_even_lag_signature
```

## Prior art

Thue–Morse / digit-sum correlation structure is classical: Gelfond 1968; Mauduit–Rivat
*Ann. Math.* 171 (2010); Baake–Coons [arXiv:2209.07102](https://arxiv.org/abs/2209.07102).
Squarefree gaps: A076259; squarefree density `6/π²` (A005117); odd squarefree `3/π²` (A056911).
No prior work isolates the base-2 gap-parity lag-2 XOR cancellation, nor shows it depends on the
complete squarefree congruence structure rather than on even-admission alone; see
[`PRIOR-ART.md`](PRIOR-ART.md).
