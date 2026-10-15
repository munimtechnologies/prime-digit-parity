# The general-b structure of the prime digit-sum parity correlation

*Work by Munim, Inc., produced with AI. MV2 of Munim, Inc. Code: Apache-2.0.
Text and data: CC BY 4.0.*

## Summary

`HIERARCHY.md` and `MOMENTS.md` wrote the base-2 (Thue-Morse) prime digit-sum
bias as a *b*-point correlation of shifted base-2^*b* digit-sum parities, and
showed by hand for *b* = 2, 3, 4 that this correlation is a single "top
cumulant" only for *b* ∈ {2, 3}, with the first extra terms appearing at *b* = 4
as three κ₂·κ₂ pairings. This note gives the **general-*b*** statement and its
**closed form**:

> **Lemma.** To leading order the *b*-point correlation
> E[A₀A₁···A_{b-1}] over the primes is the sum, over every set partition of
> {0,…,b−1} whose blocks all have size ≥ 2, of the product of the joint
> cumulants of the blocks. The number of such *singleton-free* partitions is
> **A000296(b)**, the number of partitions of a *b*-set into blocks of size ≥ 2,
> with the closed form
>
> &nbsp;&nbsp;&nbsp;&nbsp; a(b) = Σ_{k=0}^{b} (−1)^k C(b,k) B_{b−k}
>
> (B = Bell numbers). The correlation is a *single* top cumulant iff a(b) = 1,
> and a(b) = 1 **iff b ∈ {2, 3}**. So the purity boundary MOMENTS.md observed at
> b ≤ 4 is a corollary that holds for every b.

Everything below is proved by hand and machine-checked exactly (integers and
`Fraction`, no floating point) in `test_general.py`, including a full b = 5
moment→cumulant round trip run on the primes below 10⁶.

## 1. The character factorisation, for every b

For an even base the digit-sum parity is the mod-2 shadow of Legendre's
digit-sum formula (see `EVENBASE.md`). For base 2^*b* it specialises to

&nbsp;&nbsp;&nbsp;&nbsp; chi_2(n) = ∏_{r=0}^{b−1} chi_{2^b}(n ≫ r),

because chi_{2^b}(m) reads the bits of *m* at positions ≡ 0 (mod *b*), and the
shift by *r* reads positions ≡ *r* (mod *b*); the *b* shifts partition all bit
positions, so the product reads every bit — i.e. chi_2(n) = (−1)^{s₂(n)}. Write
A_r(n) = chi_{2^b}(n ≫ r). `test_general.py` checks this pointwise for
b = 2, 3, 4, 5 with the independent base-2^*b* digit-sum definition.

The base-2 excess is therefore the *b*-point correlation

&nbsp;&nbsp;&nbsp;&nbsp; −E_2(X) = Σ_{p ≤ X} chi_2(p)
= Σ_{p ≤ X} A₀(p)···A_{b-1}(p) = N · E[A₀···A_{b-1}],

with the sign convention of the repository (E_q = #odd − #even = −Σ chi_q).

## 2. Moments, cumulants, and the singleton-free sum

The joint moment of a product of random variables is the sum over set
partitions of products of joint cumulants:

&nbsp;&nbsp;&nbsp;&nbsp; E[A₀···A_{b-1}] = Σ_{π ⊢ {0,…,b−1}} ∏_{B ∈ π} κ(A_i : i ∈ B).

A block of size 1 contributes κ₁(A_i) = E[A_i], the **marginal mean** of a
single shifted base-2^*b* parity over the primes. Those marginals are small:
measured on the primes below 10⁶ each |E[A_i]| ≈ 10⁻³ (the single-character
excess at 10⁶ is in the hundreds, against N = 78 498). So any partition
carrying *s* singletons is suppressed by ≈ (10⁻³)^s relative to a singleton-free
partition of the same data. At leading order only the **singleton-free**
partitions survive — proving the Lemma's first half. This is an
order-of-magnitude statement, not a convergence theorem; the honest limits are
in §5.

## 3. Counting the survivors: OEIS A000296

The number of partitions of a *b*-set into blocks of size ≥ 2 is the integer
sequence **A000296**:

| b | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|----|----|
| a(b) | 1 | 1 | 4 | 11 | 41 | 162 | 715 |

`test_general.py` reproduces a(b) for b = 2…8 three independent ways: brute
enumeration of all set partitions; the Bell inclusion-exclusion closed form
a(b) = Σ_k (−1)^k C(b,k) B_{b−k}; and the explicit block-size tables of §4.

The closed form makes the purity boundary immediate. a(b) counts the full set
plus everything else; a(b) = 1 means *only* the full block survives, i.e. the
correlation is a single top cumulant κ_b. Since a(2) = a(3) = 1 and a(b) > 1 for
all b ≥ 4, the correlation is pure **iff b ∈ {2, 3}** — the boundary
`HIERARCHY.md`/`MOMENTS.md` found at b ≤ 4, now established for all b.

## 4. The structure just past the boundary (b = 5, b = 6)

Grouping the singleton-free partitions by the sizes of their blocks gives the
exact cumulant content of each correlation. The first two cases beyond the
boundary:

- **b = 5** (a = 11): κ₅ **+** 10·(κ₃·κ₂). The ten come from the C(5,3) = 10
  ways to pick the size-3 block.
- **b = 6** (a = 41): κ₆ **+** 15·(κ₄·κ₂) **+** 10·(κ₃·κ₃) **+** 15·(κ₂·κ₂·κ₂),
  with 15 = C(6,2), 10 = C(6,3)/2, 15 = 6!/(2!³·3!); total 1+15+10+15 = 41.

`test_general.py` verifies both tables exactly, and — on the primes below 10⁶ —
performs the full **b = 5 round trip**: it builds all 32 joint moments of
A₀…A₄, Möbius-inverts them to the joint cumulants, and reconstructs the base-2
excess from the Bell(5) = 52 partition sum, recovering the anchor exactly
(N·E[A₀…A₄] = −4764 = −E_2(10⁶)). Of those 52 partitions exactly 11 are
singleton-free (one κ₅, ten κ₃·κ₂); the other 41 each carry at least one small
marginal and are the leading-order-suppressed remainder.

## 5. What is verified, and the honest limits

Verified exactly (integers / `Fraction`, stdlib only, in `test_general.py`):
the character factorisation for b = 2,3,4,5; the anchor excesses
E_b(10⁶) = 4764, 2902, −510, 244 for the power-of-two bases 2, 4, 8, 16 under
the repo's sign convention; A000296(b) for b = 2…8 by three independent methods; the b = 5 and
b = 6 block-size tables; and the full b = 5 moment→cumulant round trip on the
primes below 10⁶ reproducing −4764.

Limits, stated plainly: the "only singleton-free partitions survive" step is a
finite-N order-of-magnitude argument from the measured ≈10⁻³ marginals, not a
proof that the suppressed terms vanish asymptotically; a rigorous bound would
need uniform control of the single-character prime sums (Mauduit-Rivat gives
such control for the asymptotic density, but not a term-by-term cumulant bound
here). The cumulant *magnitudes* at a given b still require the shifted
base-2^*b* counts from the sandbox (two such jobs are queued); this note settles
the *combinatorial structure* — which terms exist and how many — which is
compute-free.

## 6. Prior art

Gupta & Sridhar (arXiv:1812.11841) report one base-10 cutoff; Mauduit-Rivat
(Annals 2010) and Drmota give asymptotic balance with a power-saving error, an
upper bound only. Neither they, OEIS's digit-sum-parity entries (A000069,
A130911, A380192, A111309), nor this repository's own `HIERARCHY.md`/`MOMENTS.md`
state the general-*b* singleton-free decomposition, identify the survivor count
with **A000296**, give its Bell closed form, derive the purity boundary as a
corollary for all *b*, or give the worked b = 5 / b = 6 cumulant tables and the
b = 5 prime-data round trip. A prior-art search (web, arXiv, OEIS, GitHub) found
no prior statement of these; A000296 is a known sequence but has not, to our
search, been connected to the prime digit-sum parity correlation.

---

*Reproduce:* `python -m unittest -v test_general.py` (Python 3, standard library
only). Reference implementation in `general.py`.
