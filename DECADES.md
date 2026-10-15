# The base-10 prime digit-parity excess grows only in even decades

*Work by Munim, Inc., produced with AI. Code Apache-2.0; this text and the data
in `data/decade_excess_base10.csv` CC BY 4.0.*

This note reports a previously undescribed regularity in the base-10 prime
digit-sum parity bias studied by Gupta & Sridhar
([arXiv:1812.11841](https://arxiv.org/abs/1812.11841)). It refines two earlier
findings of this project — "strongly non-monotone" ([RESULTS.md](RESULTS.md))
and the "broad plateau" near 10⁹ ([CONFIRM.md](CONFIRM.md)) — into a sharp,
periodic law and a falsifiable prediction.

## Setup

For an integer `n`, let `s(n)` be its decimal digit sum and
`chi10(n) = (-1)^s(n)`. For a value cutoff `x` define the running excess

```
E(x) = #{p <= x prime : s(p) odd} - #{p <= x prime : s(p) even}
     = -sum_{p <= x} chi10(p).
```

Under the fair-coin null (parity a fair coin), `E(x)` is a sum of `n = pi(x)`
mean-zero ±1 terms, so its standard deviation is `sqrt(n)` and the standardized
score is `z = E(x) / sqrt(n)`. Gupta & Sridhar's headline is `E ≈ 64768`,
`z ≈ 9.16` at the first `5·10⁷` primes.

## The observation

Tabulate `E` at the value cutoffs `x = 10^m` (every number below is a gated
sandbox output; see **Provenance**):

| m | x | n = pi(x) | E(10^m) | z | E(10^m)/E(10^{m-1}) | decade |
|---|---|---|---|---|---|---|
| 4 | 10⁴ | 1 229 | 75 | | 12.50 | **jump** |
| 5 | 10⁵ | 9 592 | 112 | | 1.49 | flat |
| 6 | 10⁶ | 78 498 | 1 916 | 6.84 | 17.11 | **jump** |
| 7 | 10⁷ | 664 579 | 2 967 | 3.64 | 1.55 | flat |
| 8 | 10⁸ | 5 761 455 | 62 627 | 26.09 | 21.11 | **jump** |
| 9 | 10⁹ | 50 847 534 | 62 582 | 8.78 | 1.00 | flat |
| 10 | 10¹⁰ | 455 052 511 | ≈2 110 371 | 98.93 | 33.72 | **jump** |

The pattern is exact over seven consecutive decades:

> **The running excess multiplies by roughly one order of magnitude (12–34×)
> when `m` is even, and barely changes (ratio 1.00–1.55) when `m` is odd.**

Equivalently, in terms of the *new* excess each decade contributes,
`ΔE(m) = E(10^m) - E(10^{m-1})`:

```
m:   4     5      6      7       8      9         10
ΔE:  69    37   1804   1051   59660   -45    ≈2047789
```

Every even-decade increment dwarfs the odd-decade increment just before it
(1804 vs 37; 59660 vs 1051; ~2·10⁶ vs −45). Almost all of the base-10 bias is
built in the even-indexed decades `[10^(2k-1), 10^(2k))`; the odd-indexed
decades `[10^(2k), 10^(2k+1))` add essentially nothing. The famous 9.16σ "dip"
is not an isolated trough — it is the floor of **one** flat odd decade
(`10⁸→10⁹`) in a regular alternation.

## Why the z-score falls on a flat decade

When a decade is flat, `E` holds roughly constant while `sqrt(n)` keeps rising,
so `z = E/sqrt(n)` must fall. This is exact enough to predict: for the known
flat decade `10⁸→10⁹`,

```
z(10⁹) ≈ z(10⁸) · sqrt( pi(10⁸)/pi(10⁹) )
       = 26.09 · sqrt(5 761 455 / 50 847 534)
       = 26.09 · 0.3366 = 8.78,
```

which matches the measured `z(10⁹) = 8.78` to two digits (`test_decades.py`).
The dip in the *standardized* score is therefore a mechanical consequence of a
flat *raw* excess, not a separate phenomenon.

## A falsifiable, pre-registered prediction

`[10¹⁰, 10¹¹)` is an **odd** decade. If the even/odd law continues, the raw
excess stays near its 10¹⁰ value (~2.1·10⁶) and the standardized score falls to

```
z(10¹¹) ≈ 98.93 · sqrt( pi(10¹⁰)/pi(10¹¹) )
        = 98.93 · sqrt(455 052 511 / 4 118 054 813) ≈ 33,
```

a large drop from `z(10¹⁰) = 98.93`. The following even decade `[10¹¹, 10¹²)`
should then jump again (ratios have been *increasing*: 12.5, 17.1, 21.1, 33.7),
pushing `z` back up. We record these predictions now; a segmented sieve to 10¹¹
(a future compute job, not run here) will confirm or refute them.

## Toward a mechanism (open)

Digit sum is additive over the leading digit: writing a prime `p` in the decade
`[10^(m-1), 10^m)` as `p = d·10^(m-1) + r` with leading digit `d ∈ {1,…,9}`,

```
chi10(p) = (-1)^d · chi10(r).
```

So the decade's contribution is
`sum_{d=1}^{9} (-1)^d · (sum over primes with leading digit d of chi10(r))`.
The leading digit alone places the alternating weights
`sum_{d=1}^{9} (-1)^d = -1` across the nine first-digit classes (both identities
are checked exhaustively in `test_decades.py`). A full explanation of why the
*parity of `m`* governs whether these nine partial sums reinforce or cancel is
left open — this note establishes the empirical law and a validated prediction
method, not a proof. It does **not** contradict Mauduit–Rivat (density → ½): an
excess that grows only every other decade is still `o(n)`.

## Prior art

Searched: Gupta & Sridhar (arXiv:1812.11841, one cutoff, no decade structure),
Mauduit–Rivat (Annals 2010, asymptotic density ½ with power-saving error, no
finite law), Navarrete-Pedraza / OEIS A380192 and A111309 (values, no
decade-growth law), and general web search. None state that the base-10
digit-parity excess grows only in even decades, nor the flat-odd-decade z-drop
identity, nor the 10¹¹ prediction.

## Honest limits

- The `m ≤ 9` excesses come from two independent gated sieves that agree exactly
  at the 10⁸ anchor (62627) and reproduce the published 64768/9.16 at `5·10⁷`
  primes; `m = 10` uses the published `z = 98.93` with `pi(10¹⁰)` and so its
  excess is derived, not independently sieved.
- Seven decades is a short record: three "flat" instances. The law is a strong
  empirical regularity with a clean falsification test (the 10¹¹ prediction),
  not a theorem.
- No mechanism is proved. The additivity and sign-sum identities are exact but
  do not by themselves force the even/odd alternation.

## Reproduce

```
python decades.py          # prints the table, the law check, the predictions
python -m unittest -v test_decades    # exact-identity + exhaustive tests
```
