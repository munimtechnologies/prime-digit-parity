# On the digit-sum parity bias of primes: value cutoffs vs count cutoffs

*Munim, Inc. — research note. Produced with AI. Apache-2.0 (code) / CC BY 4.0 (text).*

> **v0.1.2 update (correction).** An earlier draft of this note (labelled v0.2.0)
> said the Gupta & Sridhar "peculiarity" is *explained* by the Mauduit–Rivat
> theorem. That over-read the theorem: it gives density 1/2 (so no limiting bias)
> and only an *upper* bound on the error term, which cannot predict a bias of a few
> sigma at 5·10^7 primes. The honest status is in [`RESOLUTION.md`](RESOLUTION.md):
> limit statement settled, finite-sample lean open and under measurement. This note
> keeps the elementary facts (now six, all tested) that bound the integer side.

## Background and prior art

For an integer `n` and base `b`, let `s_b(n)` be the sum of its base-`b` digits.
The parity `s_b(n) mod 2` of primes has a well-studied history:

- **Base 2** ("odious" = odd digit sum, "evil" = even). Mauduit & Rivat (*Annals
  of Math.* 171, 2010) proved the base-2 digit sums of primes are equidistributed
  mod 2 (density 1/2 each), resolving a 1968 conjecture of Gelfond. Shevelev
  (arXiv:0706.0786, 0707.1761) conjectured the odious primes stay in the majority
  for all `n > 6`; this race (OEIS A130911) is verified to **10^19**
  (B. Chaffin, Jun 2026; T. D. Noe to 10^13 earlier). So base 2 is thoroughly mapped.
- **Base 10**. Gupta & Sridhar, *A Peculiarity in the Parity of Primes*
  (arXiv:1812.11841, math.GM), report that over the **first 50,000,000 primes** the
  digit sum is "significantly more likely to be odd than even," with the odd count
  about *an order of magnitude farther from the mean than expected* — a bias they
  say they have no explanation for, persisting across bases (with the leading
  parity possibly changing).
- Navarrete-Pedraza (OEIS A380192, 2025) shows prime digit sums are approximately
  normal and uniform mod `n`, consistent with Mauduit–Rivat.

**What Mauduit–Rivat gives for base 10.** The theorem is not restricted to base 2:
it gives equidistribution of `s_q(p)` for every base `q` outside the casting-out
degeneracies, base 10 included (`q−1 = 9` is odd, so parity is non-degenerate).
So there is no limiting bias; the finite-sample excess is a lower-order term whose
size the theorem only upper-bounds. Honest limits in [`RESOLUTION.md`](RESOLUTION.md).

## Facts we can prove directly (no computation)

**Fact 1 — odd bases are trivial.** If `b` is odd then `b ≡ 1 (mod 2)`, so for
`n = Σ d_i b^i` we have `n ≡ Σ d_i = s_b(n) (mod 2)`. Hence `s_b(p)` is *odd* for
every prime `p > 2`, in every odd base. Digit-sum parity is only a non-trivial
statistic in **even** bases. (A parity "bias" reported in an odd base is just the
statement that primes are odd.)

**Fact 2 — even-base ranges are exactly balanced.** For even `b` and any `k ≥ 1`,
the digits `0..b-1` split into exactly `b/2` even and `b/2` odd values, so the
units digit alone contributes an unbiased ±1 to the digit-sum parity,
independent of the other digits. Therefore the integers in `[0, b^k)` are split
**exactly 50/50** by digit-sum parity. (Generating function: the parity count is
the coefficient symmetry of `(b/2 + (b/2) x)^k`.)

**Fact 3 — the integer excess is bounded by 1 at *every* cutoff (sharpens Fact 2,
and corrects an overstatement in v0.1.0).** Write any cutoff as `X = b·q + t` with
`t = X mod b`. Group the integers in `[0, X)` into complete blocks of `b`
consecutive integers `b·m .. b·m + (b−1)`. Within such a block the high part
`s_b(m)` is constant while the units digit runs through `0..b-1`; for **even** `b`
that is exactly `b/2` odd and `b/2` even digit sums, so **every complete block has
excess 0**. Only the final partial block of `t` integers (`r = 0..t−1`, parities
`s_b(q)+r`) can be unbalanced: if `t` is even it is balanced too, and if `t` is
odd it has exactly one extra of parity `s_b(q)`. Hence the odd−even integer excess
over `[0, X)` is

```
E(X) = 0                     if (X mod b) is even,
E(X) = +1                    if (X mod b) is odd and s_b(⌊X/b⌋) is odd,
E(X) = −1                    if (X mod b) is odd and s_b(⌊X/b⌋) is even,
```

so **|E(X)| ≤ 1 for all X**. (For `b = 2` this is the classical boundedness of
the Thue–Morse ±1 summatory function, `s(2N)=0`, `s(2N+1)=±1`. The even-`b`
version is the elementary *full-range* case; it is distinct from the much-studied
sums of `(-1)^{s_b(n)}` restricted to **arithmetic progressions** — the
generalized-Thue–Morse "Newman phenomenon" of Drmota & Stoll (Discrete Math. 308,
2008) — which can grow fractally with `N`.) The closed form is
`core.integer_parity_excess_closed` and is tested against brute force for
`b ∈ {2,4,6,8,10,16}` and all `X < 2000`.

> **Consequence.** v0.1.0 claimed a prime-count cutoff "lands deep inside a decade
> — precisely where the boundary term is largest and uncontrolled." Fact 3 shows
> that is wrong for the *integers*: the integer boundary term is never larger than
> 1, at any cutoff. So the choice of cutoff **cannot** inflate a large bias
> *through the ambient integers*. Whatever the reported base-10 prime bias is, it
> is a property of the primes, not of where the integer range is cut.

**Fact 4 — the base-10 bias lives in the high-order digits.** Every prime `p > 5`
ends in a units digit from `{1,3,7,9}`, all **odd**. Since
`s_10(p) = (p mod 10) + s_10(⌊p/10⌋)`, the odd units digit flips the parity:

```
s_10(p) mod 2 = 1 XOR (s_10(⌊p/10⌋) mod 2)      for every prime p > 5.
```

So the odd-digit-sum primes `≤ x` are exactly the primes `p ≤ x` whose
**"beheaded" value `⌊p/10⌋` has an even digit sum**, and the entire base-10 prime
parity bias equals (minus) the digit-sum parity bias of the beheaded primes. The
anomaly is therefore not about the last digit at all — it is about how
`⌊p/10⌋` distributes over digit-sum parity classes. (Verified on all primes
`< 50000` in `test_core.py`.) For even base `b > 2` the analogue is weaker: a
prime's units digit is only guaranteed *coprime* to `b`, not odd, so Fact 4 is
special to base 10 (and base 2, where it gives `s_2(p) = 1 XOR s_2(⌊p/2⌋)`).

**Fact 5 — casting out (b−1)'s pins the degeneracy.** Because `b ≡ 1 (mod b−1)`,
for every `n` we have `s_b(n) ≡ n (mod b−1)`. Hence the digit-sum residue modulo
any `d` dividing `b−1` is fixed by `n` itself. For parity: when `b` is odd,
`2 | (b−1)` forces `s_b(n) ≡ n (mod 2)`, which is Fact 1 (degenerate); when `b` is
even, `b−1` is odd and parity is **non-degenerate**. Base 10 has `b−1 = 9`, so it
is non-degenerate — exactly the condition `(q−1)/2 ∉ ℤ` under which the
Mauduit–Rivat theorem gives equidistribution of `s_10(p) mod 2`. (`core.parity_degenerate`;
the general congruence is tested for `b ∈ {3,4,7,9,10,16}` and all `n < 3000`.)

**Fact 6 — a bounded control (new in v0.1.2).** Write an integer coprime to 10 as
`n = 10m + d`, `d ∈ {1,3,7,9}`. Every such `d` is odd, so `s_10(n) ≡ s_10(m) + 1
(mod 2)`. With `X = 10M + r`, each complete block `m < M` contributes four integers
of parity `1 XOR s_10(m)`, so the excess over those blocks is `−4·E(M)` where
`E(M) ∈ {−1,0,1}` is Fact 3's integer excess; the partial block `m = M` adds at most
four more. Hence the odd−even excess over the integers `< X` coprime to 10 satisfies
`|E| ≤ 8` for **every** `X` (closed form `core.coprime10_parity_excess_closed`,
checked against brute force for all `X < 3000` and bounded far out). Primes `> 5` live
inside this set, so any `√n`-scale prime excess is a property of the primes, measured
against an ambient set that is balanced to within a constant.

Facts 1–6 are all exercised by `test_core.py` with no third-party dependencies.

## The question this repo asked, and where it stands

The excess over primes `≤ x` is, as a sequence, *identical* to the excess over the
first `π(x)` primes — the two conventions index the same running total and differ
only in how it is **normalized**. Facts 3 and 5 remove the ambient integers and
the odd bases as sources of inflation. That left a sharp question about base 10:

> Under a fair-coin null, the base-10 prime excess `excess(x) = odd − even` would
> have `z(x) = excess(x)/√π(x) = O(1)`. Does the Gupta–Sridhar "order of magnitude"
> hold up at value cutoffs `x = 10^k`?

**Status (honest).** By Mauduit–Rivat, `|excess(x)| = O(x^{1−σ})` with `σ > 0`, so
the parity classes have density 1/2 and any bias is lower-order. But `z(x)` is only
*permitted* to be large (up to `x^{1/2−σ}(log x)^{1/2}`); the theorem does not say it
is, and it does not explain a few-sigma bias at 5·10^7 primes. Whether `z` stays
O(1) and changes sign (a single-cutoff fluctuation) or keeps one sign and grows (a
real secondary term, like Chebyshev-type bias in prime races) is the open empirical
question; the queued computation answers it (value cutoffs `10^k`, count cutoffs
including `N = 5·10^7`, a per-segment trajectory, all six bases). Readings are
pre-registered in [`RESOLUTION.md`](RESOLUTION.md) §5.

## Method

`fast.py` streams primes with a segmented sieve and computes, for bases
`{2,4,6,8,10,16}`, the odd-minus-even excess at value cutoffs `10^k`, at
prime-count cutoffs every `10^6` primes, and at the end of every segment.
`job.py` reports `z = excess/√n` at all of these.

**Status: the computation is queued in the sandbox (`prime-digit-parity-v3`: first
to `10^9`, then to `10^10`; re-queued at low memory after the original 8 GB request
would not schedule).** The verified `z(10^k)` table is appended to `RESOLUTION.md`
once the job returns. No numeric claim about primes is asserted before it is computed.

## References

- C. Mauduit, J. Rivat, *Sur un problème de Gelfond: la somme des chiffres des
  nombres premiers*, Ann. of Math. 171 (2010) 1591–1646.
- A. O. Gelfond, Acta Arith. 13 (1968) 259–265.
- V. Shevelev, arXiv:0706.0786; arXiv:0707.1761.
- A. D. Gupta, M. Sridhar, *A Peculiarity in the Parity of Primes*, arXiv:1812.11841.
- E. Navarrete-Pedraza, OEIS A380192 (2025).
- M. Drmota, T. Stoll, *Newman's phenomenon for generalized Thue–Morse
  sequences*, Discrete Math. 308 (2008) 1191–1208.
- OEIS A027697, A027699, A130911, A000069, A001969.
