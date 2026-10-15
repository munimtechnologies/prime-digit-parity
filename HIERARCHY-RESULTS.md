# HIERARCHY-RESULTS.md — measuring the b=2 rung of the cumulant tower

*Work by Munim, Inc., produced with AI. Part of
[`munimtechnologies/prime-digit-parity`](https://github.com/munimtechnologies/prime-digit-parity).
Text and data: CC BY 4.0. Code: Apache-2.0. — MV2 of Munim, Inc.*

[`HIERARCHY.md`](HIERARCHY.md) (v0.4.3) proved that the single base-2 (Thue–Morse)
prime digit-sum bias has a whole tower of equivalent representations — a b-point
joint correlation of base-2^b parities, one for every binary block width b — and
predicted a sharp boundary: the bias is a **pure top-order connected cumulant**
only for b ∈ {2, 3}. [`CORRELATION-RESULTS.md`](CORRELATION-RESULTS.md) (v0.4.8)
measured the **b = 3** rung: at x = 1e8 the connected third cumulant κ₃ of three
shifted base-8 parities carries **99.84 %** of the base-2 excess.

This note measures the **b = 2** rung. A sealed sandbox job
(`prime-parity-pairwise-v1`, pure-stdlib segmented work over every prime below
1e8, 67 s on an isolated container, no network) returned the marginals, joint
mean and pairwise cumulant κ₂ for two exact instances of the factorisation
(★) `χ_{2^a}(n) = ∏_r χ_{2^{ab}}(⌊n/2^{ra}⌋)`:

- **base 2 from base 4** (a=1, b=2): `χ_2(p) = χ_4(p) · χ_4(p//2)`
- **base 4 from base 16** (a=2, b=2): `χ_4(p) = χ_16(p) · χ_16(p//4)`

All numbers below are from the x = 1e8 snapshot, which is taken over **all primes
≤ 1e8** (π(1e8) = 5 761 455) and reproduces every published excess anchor exactly.

## 1. The measured decomposition at x = 1e8

Excess convention: `E_q = #{s_q odd} − #{s_q even} = −π(x)·𝔼[χ_q(p)]`. For a b=2
rung, `𝔼[X_0 X_1] = κ_2 + 𝔼[X_0]𝔼[X_1]`, so the excess splits into a **connected
part** `−π·κ_2` and a **marginal-product part** `−π·𝔼[X_0]𝔼[X_1]`.

### base 2 from base 4 — the Thue–Morse bias as a pairwise object

| quantity | value |
|---|---|
| marginal `m₀ = 𝔼[χ_4(p)]`        | −0.011391740 |
| marginal `m₁ = 𝔼[χ_4(p//2)]`     | −0.015184532 |
| joint mean `𝔼[χ_4(p)·χ_4(p//2)]` = `𝔼[χ_2(p)]` | −0.012238054 |
| pairwise cumulant `κ₂`           | −0.012411033 |
| marginal product `m₀·m₁`         | +0.000172978 |
| `m₀·m₁ / κ₂`                      | −0.013937 |

Scaled by π(1e8), the **base-2 excess E₂ = 70 509** splits as

> **connected pairwise cumulant  +71 506  (101.41 %)**
> **marginal product            −997      (−1.41 %)**

### base 4 from base 16

| quantity | value |
|---|---|
| marginal `m₀ = 𝔼[χ_16(p)]`       | −0.018155657 |
| marginal `m₁ = 𝔼[χ_16(p//4)]`    | −0.005702553 |
| joint mean = `𝔼[χ_4(p)]`         | −0.011391740 |
| pairwise cumulant `κ₂`           | −0.011495274 |
| marginal product `m₀·m₁`         | +0.000103534 |
| `m₀·m₁ / κ₂`                      | −0.009007 |

The **base-4 excess E₄ = 65 633** splits as

> **connected pairwise cumulant  +66 230  (100.91 %)**
> **marginal product            −597      (−0.91 %)**

## 2. Scoring the pre-registered predictions

The three predictions in [`HIERARCHY.md`](HIERARCHY.md) §5 were fixed **before**
the job ran. Verdicts:

- **P-pair-1 (marginal product subdominant to κ₂, |m₀m₁| < 10 % of |κ₂|): CONFIRMED.**
  Measured at **1.39 %** (base-2 rung) and **0.90 %** (base-4 rung). The honest
  caveat in the pre-registration — that the base-16 marginal ≈1.8e-2 might push the
  correction to a few percent — is borne out but stays small: the b=2 bias is, to
  leading order, the pairwise connected cumulant.

- **P-pair-2 (tower consistency, −π·κ₂⁽⁴⁾ for base 2 ≈ −π·κ₃⁽⁸⁾): CONFIRMED to
  leading order.** The b=2 leading term is +71 506 and the b=3 leading term
  (v0.4.8) is +70 395; both approximate the same E₂ = 70 509. They are the same
  bias seen at two resolutions. **New observation:** the b=3 rung is the *purer*
  one — its correction is 0.16 %, while the b=2 correction is 1.41 %. This matches
  the theory's ordering: the b=2 correction is a marginal product of order m², the
  b=3 correction a mix of order m·κ₂; here the single-prime b=3 split happens to
  cancel more completely.

- **P-pair-3 (excess anchors at 1e8): CONFIRMED exactly.** The job's directly
  summed odd−even excesses are **E₂ = 70 509, E₄ = 65 633, E₈ = −1329,
  E₁₆ = 104 603**, matching `prime-parity-trajectory-v1` to the unit, with
  π(1e8) = 5 761 455 exact. No sieve bug.

## 3. What the tower now looks like, measured

Putting the two measured rungs beside v0.4.8, the one Thue–Morse prime digit-sum
bias is, at every block width where the theory says it should be pure, essentially
a **single connected top-order joint cumulant** of finer-base parities:

| rung | identity | connected share of the bias | correction |
|---|---|---|---|
| b=2, base 4 (bias of base 2)  | `χ_2=χ_4·χ_4(·//2)`   | **101.4 %** (κ₂) | −1.4 % (m₀m₁) |
| b=2, base 16 (bias of base 4) | `χ_4=χ_16·χ_16(·//4)` | **100.9 %** (κ₂) | −0.9 % (m₀m₁) |
| b=3, base 8 (bias of base 2)  | `χ_2=χ_8·χ_8(·//2)·χ_8(·//4)` | **99.84 %** (κ₃) | +0.16 % |

The "pure top-order cumulant for b ∈ {2, 3}" boundary of HIERARCHY.md is now not
just proved as algebra but **confirmed as measurement**: across all three pure
rungs the connected cumulant accounts for 99.8–101.4 % of the excess, and the
lower-order (marginal / pairing) terms never exceed ~1.4 %.

## 4. Reproduce it

`test_hierarchy_results.py` is self-contained (stdlib only, no repo imports, no
file reads). It sieves every prime below 1e6, checks the two character
factorisations pointwise, reproduces the published 1e6 excess anchors
(E₂=4764, E₄=2902, E₈=−510, E₁₆=244) under the exact sign convention, verifies
`κ₂ + m₀m₁ = 𝔼[joint]` for both rungs, confirms the connected cumulant dominates
the marginal product at 1e6, and arithmetically re-checks that the published 1e8
splits sum to the anchor excesses. The full per-checkpoint measurement is in
[`data/pairwise_cumulant.csv`](data/pairwise_cumulant.csv) (CC BY 4.0).

## Prior art and honest limits

Mauduit–Rivat (2010) give a limiting density with a power-saving error but no
finite-size law or cumulant structure; Gupta–Sridhar (arXiv:1812.11841) report
the base-10 excess and its cross-base persistence as unexplained; OEIS
A000069/A130911/A380192/A111309 tabulate values without asymptotics. Decomposing
the Thue–Morse prime bias into a **connected pairwise cumulant of coarser-base
parities**, and confirming the b∈{2,3} purity boundary by direct measurement, is,
to a prior-art search across the literature, GitHub, OEIS and arXiv, new.

Limits: the 1e8 numbers are from a single sandbox run, gated on four independent
excess anchors and π(1e8); they are point estimates, not fitted asymptotics. The
x = 1e7 snapshot in the CSV is recorded at the first prime **above** 1e7
(p = 10 000 019, so n = 664 580 = π(1e7)+1) because of the checkpoint boundary,
and is kept only for completeness — all headline claims use the exact 1e8
snapshot. None of this explains *why* the connected cumulants are large along the
primes; that remains the open mechanism.
