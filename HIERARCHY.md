# HIERARCHY.md — the cumulant tower of the base-2 prime digit-sum bias

*Work by Munim, Inc., produced with AI. Part of
[`munimtechnologies/prime-digit-parity`](https://github.com/munimtechnologies/prime-digit-parity).
Text and data: CC BY 4.0. Code: Apache-2.0. — MV2 of Munim, Inc.*

This note generalises two earlier results in this project and unifies them:

- **REDUCE.md** (v0.4.1): the exact inter-base identity
  `s_{2^a}(n) = ⊕_{r=0}^{b-1} s_{2^{ab}}(⌊n / 2^{ra}⌋) (mod 2)`.
- **CORRELATION.md** (v0.4.2): the base-2 bias equals, to leading order, the
  **third** joint cumulant κ₃ of three shifted base-8 parities.

Here we show these are two rungs of a single ladder: the *one* base-2
(Thue–Morse) prime digit-sum bias has **a whole tower of equivalent
representations**, one for every block width b, as a b-point joint correlation
of base-2^b parities. We give the exact decomposition at each rung, prove where
the bias is a *pure* top-order cumulant and where it is not, and apply the b=2
rung to explain the base-4 and base-16 biases the same way. Everything in §1–§4
is exact algebra, machine-checked in `test_hierarchy.py`; the numeric claims in
§5 are pre-registered predictions for the sandbox, not yet confirmed.

## 1. The character factorisation (exact)

Write `q = 2^a`, `Q = 2^{ab}`, and the digit-sum parity character
`χ_q(n) = (-1)^{s_q(n)} ∈ {+1,-1}`.

**Lemma (positional form).** For any `a ≥ 1`,
`s_{2^a}(n) ≡ Σ_j b_{ja} (mod 2)`, where `n = Σ_i b_i 2^i` is the binary
expansion. That is, the base-2^a digit-sum parity is the parity of the binary
digits of `n` at positions `≡ 0 (mod a)`.

*Proof.* The j-th base-2^a digit is `d_j = Σ_{i=0}^{a-1} b_{ja+i} 2^i`, so
`s_{2^a}(n) = Σ_j d_j = Σ_{i=0}^{a-1} 2^i (Σ_j b_{ja+i})`. Modulo 2 only the
`i=0` term survives, giving `Σ_j b_{ja}`. ∎

**Corollary (reduction identity).** The b shifts `m_r = ⌊n/2^{ra}⌋`
(`r = 0,…,b-1`) move binary position `p` to `p+ra`, so `s_Q(m_r)` reads the bits
of `n` at positions `≡ ra (mod ab)`. As r runs `0..b-1` these residue classes
are `0, a, 2a, …, (b-1)a (mod ab)` — exactly, and disjointly, the positions
`≡ 0 (mod a)`. Summing the Lemma over the b classes gives
`s_{2^a}(n) ≡ Σ_{r=0}^{b-1} s_Q(⌊n/2^{ra}⌋) (mod 2)`, i.e.

> **χ_{2^a}(n) = ∏_{r=0}^{b-1} χ_{2^{ab}}(⌊n / 2^{ra}⌋).**   (★)

The partition is disjoint (each binary position lands in exactly one class), so
(★) is an identity, not an approximation. `test_hierarchy.py` checks it
pointwise for all `(a,b)` with `a∈{1,2,3}`, `b∈{1,2,3,4}` over `n < 2^13`.

## 2. From the identity to a correlation of parities

Let the odd−even excess be `E_q(x) = #{p≤x: s_q(p) odd} − #{p≤x: s_q(p) even}
= −Σ_{p≤x} χ_q(p)` (sign convention of CORRELATION.md: base-10 excess is
positive, matching Gupta–Sridhar). Put `X_r(p) = χ_Q(⌊p/2^{ra}⌋)`. Then (★)
gives, over the primes `p ≤ x`,

> **E_{2^a}(x) = −π(x) · 𝔼[ X_0 X_1 ⋯ X_{b-1} ],**

where 𝔼 is the average over primes `p ≤ x`. The base-2^a bias is literally a
**b-fold joint moment** of b shifted base-2^{ab} parities.

## 3. Moment → cumulant, and where the top cumulant is alone

By the Leonov–Shiryaev formula the joint moment expands over set partitions π of
`{0,…,b-1}`:
`𝔼[∏_r X_r] = Σ_π ∏_{B∈π} κ(X_B)`, with κ the joint cumulant of a block.
The two extreme partitions are the single block `{0,…,b-1}` (the **connected**
b-th cumulant `κ_b`) and all-singletons (`∏_r 𝔼[X_r]`, a product of marginals).

The base-2^{ab} marginals `𝔼[X_r]` are themselves tiny biases, so any partition
that contains a singleton is suppressed by at least one marginal factor. Writing
m for the typical marginal size, a partition into blocks of sizes `s_1,…,s_k`
contributes at order `∏ (connected part) × m^{(# singletons)}`. Therefore:

- **b = 2:** the only partition with no singleton is `{{0,1}}`. So
  `𝔼[X_0 X_1] = κ_2 + 𝔼[X_0]𝔼[X_1]`, and the bias is the pairwise
  (second) cumulant **plus a marginal-product correction of order m²**.
- **b = 3:** the only no-singleton partition is `{{0,1,2}}`. So the leading
  term is `κ_3` — the CORRELATION.md result — with corrections of order
  `m·κ_2` and `m³`.
- **b ≥ 4:** no-singleton partitions are `κ_b` **plus** products of smaller
  cumulants (for b=4: three pairings `κ_2·κ_2`). The bias is **no longer a pure
  top-order cumulant**: at b=4 it is `κ_4` *and* the three `κ_2κ_2` pairings,
  both at leading order (no marginal suppression).

So the clean "the bias *is* the connected b-point correlation" statement holds
**exactly for b ∈ {2,3} and no further** — a sharp, previously unstated boundary.
`test_hierarchy.py` verifies the full moment–cumulant expansion for b=2,3,4 as an
exact Fraction identity, and that the b=4 no-singleton sum equals
`κ_4 + Σ(3 pairings)` with the pairing part provably non-zero.

## 4. The tower for base 2 (a = 1), and the exact integer baseline

Setting `a=1`, (★) gives one representation of the *same* `χ_2(n)` for every b:

| b | base Q=2^b | identity | leading order |
|---|-----------|----------|---------------|
| 2 | 4  | `χ_2(n)=χ_4(n)·χ_4(n//2)`                         | κ₂ of base-4 parities |
| 3 | 8  | `χ_2(n)=χ_8(n)·χ_8(n//2)·χ_8(n//4)`               | κ₃ of base-8 parities |
| 4 | 16 | `χ_2(n)=χ_16(n)·χ_16(n//2)·χ_16(n//4)·χ_16(n//8)` | κ₄ + 3 pairings, base-16 |

All three products equal `χ_2(n)` pointwise (checked for `n<2^16`), so they give
**identical** excesses — they are the same object viewed at three resolutions.
What differs is the *decomposition*: the Thue–Morse bias is a 2-point object in
base 4, a genuinely 3-point object in base 8, and a 4-point-plus-pairings object
in base 16. Matching leading terms across rungs yields consistency relations,
e.g. `κ_2^{(4)} ≈ κ_3^{(8)}` to the order the marginals are small — a built-in
cross-check between the b=2 and b=3 sandbox jobs.

**Exact baseline.** For even base `q`, `Σ_{d=0}^{q-1}(-1)^d = 0`, so over any
complete base-Q block `[0, Q^k)` the digit-sum parity is balanced:
`Σ_{n<Q^k} χ_q(n) = 0`. By (★) this equals the summed joint correlation, so the
unnormalised correlation is **exactly zero over the full integers** at every
rung — the bias is purely a *primes* effect, not a digit artefact. Checked
exactly for base 4 over `[0,16²)` and base 2 over `[0,8²)`.

## 5. Application and pre-registered measurement (not yet confirmed)

The b=2 rung also decomposes base 4 and base 16 directly (a=2):
`χ_4(n) = χ_16(n)·χ_16(n//4)`, so `E_4(x) = −π(x)·𝔼[χ_16(p)·χ_16(p//4)]`, with
leading term the pairwise cumulant `κ_2^{(16)}`. This explains the sizeable
base-4 excess and the smaller base-16 marginal by one mechanism.

A short stdlib sieve job (`prime-parity-pairwise-v1`, primes `< 1e8`) will
measure, for both the base-4 (from base-16) and base-2 (from base-4) rungs, the
marginals `m_0,m_1`, the joint mean `𝔼[X_0X_1]`, and `κ_2 = 𝔼[X_0X_1]−m_0 m_1`.
Pre-registered predictions:

- **P-pair-1:** the marginal-product correction `m_0·m_1` is **subdominant** to
  `κ_2` (|m₀m₁| < 10% of |κ₂|), so the base-4/base-16 bias is, to leading order,
  the pairwise cumulant. (Honest caveat: at 1e8 the base-16 marginal is ≈1.8e-2,
  so the correction may be a few percent, not vanishing — the job pins it down.)
- **P-pair-2 (tower consistency):** `−π·κ_2^{(4)}` for base 2 agrees to leading
  order with `−π·κ_3^{(8)}` from `prime-parity-base8-corr-v1` (same bias, two
  rungs).
- **P-pair-3 (anchor):** the directly summed odd−even excesses at x=1e8 match the
  `prime-parity-trajectory-v1` table: E₂=70509, E₄=65633, E₁₆=104603
  (and E₈=−1329). If the anchor fails, a sieve bug is present — fix before
  trusting κ₂.

## Prior art and honest limits

Mauduit–Rivat (2010) give a limiting density with a power-saving error bound but
no finite-size law, non-monotonicity, base-dependence, or cumulant structure.
Gupta–Sridhar (arXiv:1812.11841) report the base-10 excess and its persistence
across bases as unexplained. OEIS A000069/A130911 (base 2) and A380192/A111309
(base 10) tabulate values without asymptotics. REDUCE.md gives the identity (★)
but no cumulant analysis; CORRELATION.md gives only the b=3 base-8 rung. The
**tower** (one bias, equivalent cumulant representations at every block width),
the **sharp purity boundary** (pure top cumulant iff b∈{2,3}), and the **b=2
base-4/base-16 pairwise instance** are, to our search, new.

Limits: the moment–cumulant algebra is standard; the novelty is the application,
the tower framing, and the purity boundary. The algebra in §1–§4 is exact and
machine-checked; the §5 numbers are predictions awaiting the sandbox. None of
this explains *why* the connected cumulants are large along the primes — that
remains the open mechanism.
