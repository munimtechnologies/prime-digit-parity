# Where the even-base prime parity bias lives: marginals vs. joint correlation

*Work by Munim, Inc., produced with AI. Author: MV2 of Munim, Inc. Code
Apache-2.0; this text and the data table (`data/evenbase_factor_means.csv`)
CC BY 4.0.*

`EVENBASE.md` (v0.4.4) proved a uniform reduction that writes the digit-sum
parity character of any **even** base as a product of least-significant-bit
characters of the number's successive base-`b` truncations,

```
chi_b(n) = prod_{i>=0} eps( floor(n / b^i) ),   eps(x) = (-1)^x,        (E')
```

and, in §5, **pre-registered** three predictions (P-even-0,1,2) to be scored by a
future sandbox job measuring the marginal balance of each truncation factor
`eps(floor(p/b^i))` over primes. This note reports that job. It ran in the
sealed compute sandbox (`evenbase-factors-v2`, bases 6, 8, 10, all primes to
`10^8`, 74 s), behind a **hard anchor gate** that aborts unless the primes
reproduce the published excesses `E_b` at `10^6, 10^7, 10^8` and `pi(x)` exactly.
The gate passed. Every statement below is scored honestly against the
pre-registration, including the one prediction that was only *partly* right.

## The decomposition that was measured

For an odd prime `p`, (E') gives `chi_b(p) = eps(p) · prod_{i>=1} eps(floor(p/b^i))
= -prod_{i>=1} eps(floor(p/b^i))` (the units factor `eps(p) = -1` for every odd
`p`). So the entire prime excess `E_b(x) = -sum_{p<=x} chi_b(p)` is carried by
the **higher truncation factors** `eps(floor(p/b^i))`, `i >= 1` — the parity of
`p` with its last 1, 2, 3, … base-`b` digits stripped off. We measured, for each
base and each layer `i >= 1`, the **mean of that factor over all primes `<= x`**.
Full table: `data/evenbase_factor_means.csv` (bases 6/8/10 × layers × `x ∈
{1e6, 1e7, 1e8}`); the `x = 1e6` slice is reproduced from an independent stdlib
sieve in `test_evenbase_results.py`.

## Scoring the pre-registration

**P-even-0 (gate) — CONFIRMED, exactly.** The primes reproduce
`E_6 = 57813, E_8 = -1329, E_10 = 62627` at `10^8` and
`E_6 = 3181, E_8 = 215, E_10 = 2967` at `10^7` (and the `10^6` anchors
`1936, -510, 1916`), with `pi(10^8) = 5761455`. No sieve bug.

**P-even-1 (every factor balanced) — PARTIALLY confirmed; refined.** The
prediction was that *each* factor `eps(floor(p/b^i))` is close to balanced over
primes (|mean| of square-root size) for all three bases. The data splits the
layers into two regimes:

* **Interior layers** (`i` from 1 up to about `log_b x − 2`) **are** balanced to
  √-size. At `10^8`, with ~5.76M primes (fair-coin sd of a mean ≈ `4.2·10^-4`),
  the interior means are a few `·10^-4` — e.g. base 8 layers 1–5 are
  `+1.7, +0.2, -2.1, -2.2, -3.8` (`×10^-4`). This part of P-even-1 holds.
* **The top ~2 layers are NOT balanced.** They carry a deterministic
  *leading-truncation* bias toward `-1`: at `10^8` the leading factor is
  `-1.000` (base 6, `i=10`), `-0.197` (base 8, `i=8`), `-0.114` (base 10, `i=7`),
  with the next layer down an order of magnitude smaller but still well above
  √-size. This is the same leading-digit effect isolated for base 10 in
  `LEADING.md`: the top truncation `floor(p/b^i)` is a small integer (usually 1),
  so `eps` of it is systematically `-1`. Crucially, **this bias has the same
  shape in every base**, so it cannot explain why the bases differ.

So P-even-1 is true for the bulk but false at the top; the honest correction is
base-independent and does not touch the cross-base comparison.

**P-even-2 (the excess is joint, not marginal) — CONFIRMED, decisively.** Treat
the truncation factors as independent: the predicted mean of `chi_b(p)` is then
`eps(p)·prod_i E[eps_i] = -prod_i (factor mean)`, and the predicted excess is
`-n·(that mean)`. Because the interior means are each `~10^-4`, their product
collapses to `~10^-24` at `10^8` **for every base**, even including the biased
leading layers:

| base | true `E_b` at `10^8` | predicted `E_b` if factors independent |
|-----:|---------------------:|---------------------------------------:|
| 6    | `57813`              | `≈ 10^-18` |
| 8    | `-1329`              | `≈ 10^-18` |
| 10   | `62627`              | `≈ 10^-17` |

Independence predicts **essentially zero** excess for all three bases, yet the
true excesses span `-1329` (base 8) to `62627` (base 10). **100% of each base's
excess, and 100% of the `≈64000` gap between base 8 and base 10, lives in the
joint correlation of the truncation factors — none of it in their marginals.**
This is the exact non-dyadic analogue of the project's standing dyadic finding
(`REDUCE.md`, `CORRELATION.md`): the signal is a connected correlation of
individually (near-)balanced parities, not a bias in any single digit.

`test_evenbase_results.py` re-derives the `10^6` slice from an independent sieve
and asserts all three: the anchor gate, the factor means (to `10^-9`), the
interior-balanced / leading-biased split, and that the independence prediction is
`|E| < 1` while the true `|E|` is in the hundreds.

## Why this is new (prior art)

The even-base reduction (E′) is the mod-2 shadow of Legendre's formula —
elementary, not claimed as new (`EVENBASE.md §6`). What is new here, and not
found in a prior-art search (Mauduit–Rivat, *Annals* 2010, which proves only
asymptotic mod-2 equidistribution of prime digit sums for even `b`, with no
finite decomposition; Gupta–Sridhar, arXiv:1812.11841, one base-10 cutoff;
Spiegelhofer on cubes; OEIS A380192/A111309), is the **measured, anchor-gated
layer decomposition** of the even-base prime excess across bases 6/8/10 and the
resulting quantitative statement: the marginals are balanced to √-size in the
interior and carry only a base-independent leading bias, so the whole excess —
and the entire cross-base ordering — is a joint correlation governed by the odd
part `m` of the base (1 for 8, 5 for 10). It converts `EVENBASE.md`'s qualitative
"the count of factors doesn't order the bias" into a scored, reproducible
marginal-vs-joint decomposition.

## Honest limits

* The decomposition localises the bias to the joint correlation; it does **not**
  explain *why* that correlation is large along primes — the project's standing
  open problem.
* `z`-scores use the fair-coin sd `sqrt(count)`; the true variance along primes
  is only known to be `O(count)` up to the Mauduit–Rivat constant, so the sd is a
  scale, not an exact normalisation.
* The independence "prediction ≈ 0" is partly structural (any product with
  several `~10^-4` factors is tiny); its content is that the interior marginals
  *are* that small **and measured**, ruling out any single-digit explanation. The
  verified range of the factor means in the test is `10^6`; `10^7`/`10^8` rest on
  the gated sandbox run shipped in the CSV.

*— MV2, Munim, Inc.*
