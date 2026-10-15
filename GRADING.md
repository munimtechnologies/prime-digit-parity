# The singleton grading of the prime digit-sum parity correlation

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. Code in this
repository is Apache-2.0; this text and the data are CC BY 4.0.*

This note grades the block-parity cumulant expansion of the prime digit-sum
parity bias by the **number of singleton blocks**. The grading is a small-
marginal expansion: its order-0 term is exactly the singleton-free structure
counted by OEIS A000296 (the leading-order law of [`GENERAL.md`](GENERAL.md)),
its full sum is the Bell-number moment expansion of [`MOMENTS.md`](MOMENTS.md),
and the two are joined by a single bivariate generating function. A structural
corollary falls out: the next-to-last order of every such correlation vanishes
identically. Everything here is exact (integers and `Fraction`) and is checked
both combinatorially and on the primes below 10^6 in
[`test_grading.py`](test_grading.py).

## Setup (recap)

For block width `b`, the base-2 (Thue-Morse) digit-sum character factorises
exactly as a product of `b` shifted base-2^b parities:

```
chi_2(n) = prod_{r=0}^{b-1} chi_{2^b}(n >> r) =: prod_{r=1}^{b} X_r(n),
```

(see [`EVENBASE.md`](EVENBASE.md); the parity of a base-2^b digit sum is the
parity of the bits at positions that are multiples of `b`). The base-2 excess
`E_2 = #odd - #even = -sum_p chi_2(p)` is therefore the `b`-point correlation of
`X_1, ..., X_b` over the primes. By the moment-cumulant (cluster) formula,

```
E[ X_1 ... X_b ] = sum over set partitions pi of {1..b}  prod_{B in pi} kappa(B),
```

where `kappa(B)` is the joint cumulant of the blocks' variables. A **singleton**
block `{r}` contributes the marginal mean `E[X_r] = delta_r`. Empirically each
`delta_r` is small -- at 10^6, `E[chi_8(p)] = -E_8/pi(10^6) = 510/78498 ~ 6.5e-3`
-- so a partition with `k` singletons carries `k` factors of a small marginal and
is of order `delta^k`.

## The grading law

Classify the set partitions of `{1..b}` by their number of singletons `k`. The
count is

```
P(b, k) = C(b, k) * A000296(b - k).
```

*Proof.* Choose the `k` singleton positions (`C(b,k)` ways); the remaining
`b - k` elements must form a partition with **no** singleton, and those number
`A000296(b-k)` by definition. ∎

Two boundary rows recover the results we already published:

* `k = 0` (no singletons): `P(b,0) = A000296(b)` -- the leading-order,
  all-connected term of [`GENERAL.md`](GENERAL.md).
* summed over `k`: `sum_k P(b,k) = Bell(b)` -- the full moment expansion of
  [`MOMENTS.md`](MOMENTS.md), via the classical `Bell(b) = sum_j C(b,j) A000296(j)`.

So the singleton grading is exactly the interpolation between GENERAL (the
leading piece) and MOMENTS (the whole thing), organised by powers of the small
marginal.

### Triangle P(b,k)

```
b \ k :    0    1    2    3    4    5    6  | Bell   A000296(b)
    0 :    1                                |    1         1
    1 :    0    1                           |    1         0
    2 :    1    0    1                       |    2         1
    3 :    1    3    0    1                  |    5         1
    4 :    4    4    6    0    1             |   15         4
    5 :   11   20   10   10    0    1        |   52        11
    6 :   41   66   60   20   15    0    1   |  203        41
```

(Each row b>=1 has a zero in column k=b-1 -- the vanishing next-to-last order;
see the corollary.)

## Generating function and recurrence

Marking singletons with a variable `u`, a set partition's exponential
generating function splits blocks of size 1 (contributing `u*x`) from blocks of
size >= 2 (contributing `e^x - 1 - x`):

```
sum_{b>=0}  ( sum_k P(b,k) u^k )  x^b / b!  =  exp( u*x + e^x - 1 - x ).
```

At `u = 0` this is `exp(e^x - 1 - x)`, the EGF of A000296 (leading order); at
`u = 1` it is `exp(e^x - 1)`, the Bell EGF (full expansion). Extracting the
`u^k` coefficient reproduces `P(b,k) = C(b,k) A000296(b-k)`.

Differentiating `A(x) = exp(e^x - 1 - x)` gives `A'(x) = (e^x - 1) A(x)`, i.e.
the recurrence

```
A000296(n+1) = sum_{m=0}^{n-1} C(n, m) A000296(m).
```

This is the full Bell recurrence `Bell(n+1) = sum_{m=0}^{n} C(n,m) Bell(m)` with
the top term dropped -- the combinatorial signature of forbidding singletons.
The test verifies it against brute-force enumeration through `A000296(9) = 3425`,
alongside the closed form `A000296(n) = sum_j (-1)^{n-j} C(n,j) Bell(j)` and the
Bernhart identity `Bell(n) = A000296(n) + A000296(n+1)`.

## Corollary: the next-to-last order vanishes

Because `A000296(1) = 0` (a one-element set has no singleton-free partition),

```
P(b, b-1) = C(b, b-1) * A000296(1) = 0   for every b >= 1.
```

No set partition has exactly `b-1` singletons -- if `b-1` points are singletons
the last one is forced to be a singleton too. Hence **the order-`delta^{b-1}`
term of every block-parity correlation is identically zero**: there is no
"one connected pair plus all-else-marginal" contribution at any width. For
`b = 2` this says the linear-in-marginal term is absent (the pair correlation is
a connected cumulant plus the pure product `delta_1 delta_2`, nothing between);
for `b = 3` the `delta^2` order is empty, as the row shows.

## Numerical confirmation on the primes below 10^6

`test_grading.py` sieves the 78,498 primes `p <= 10^6` and verifies, exactly in
`Fraction` arithmetic:

* the tower factorisation pointwise: `chi_4(p) chi_4(p>>1) = chi_2(p)` and
  `chi_8(p) chi_8(p>>1) chi_8(p>>2) = chi_2(p)` for every prime, and that both
  towers reproduce the anchor `E_2 = 4764` (with `E_4 = 2902`, `E_8 = -510`);
* **b = 2** (base-4): `E[AB] = kappa_2(A,B) + E[A]E[B]` with the middle order
  empty (`P(2,1)=0`), and the all-singleton term `E[A]E[B]` strictly below the
  correlation in magnitude;
* **b = 3** (base-8): the graded reconstruction
  `E[XYZ] = kappa_3 + (E[X]kappa_{YZ}+E[Y]kappa_{XZ}+E[Z]kappa_{XY}) + E[X]E[Y]E[Z]`
  holds exactly, the `delta^2` order is empty (`P(3,2)=0`), and the highest-order
  all-singleton term `E[X]E[Y]E[Z]` is orders of magnitude below both the signal
  and the leading connected cumulant `kappa_3` -- the grading in action.

This matches the earlier measurement that at 10^8 the connected triple cumulant
carries 99.84% of the base-8-to-base-2 excess
([`CORRELATION-RESULTS.md`](CORRELATION-RESULTS.md)): the `k=0` term dominates
and the singleton orders are the small corrections the grading predicts.

## Prior art and what is new

The combinatorial facts are classical: A000296 is the sequence of singleton-free
set partitions ("associated Bell numbers"), with EGF `exp(e^x-1-x)`, the
Bernhart relations `Bell(n) = A000296(n)+A000296(n+1)`, and the associated-Stirling
framing surveyed by Coquereaux (*J. Integer Sequences* 27, 2024; arXiv:2305.01100).
What is new, to the best of a prior-art search, is assembling the **singleton
grading into a small-marginal expansion of the prime digit-sum block-parity
correlation**: the explicit interpolation `P(b,k) = C(b,k)A000296(b-k)` between
the leading-order law (A000296, GENERAL) and the full moment expansion (Bell,
MOMENTS), its bivariate generating function, and the resulting
`order-delta^{b-1} = 0` corollary, all verified on the primes. No source ties
A000296 or this grading to the prime digit-sum parity bias.

### Limits

The dominance of the `k=0` term and the smallness of the singleton orders are
empirical (the marginals `delta_r` are small but not proven to vanish); the
grading law, generating function, recurrence and corollary are exact theorems.
The numerical check is at 10^6; the 10^8 dominance figure is quoted from the
earlier confirmed run, not recomputed here.

---
*Reproduce:* `python -m unittest -v test_grading.py` (pure standard library).
*Part of* `munimtechnologies/prime-digit-parity`. *— MV2, Munim, Inc.*
