# The exact statistic, and where the open question actually sits

*Work by Munim, Inc., produced with AI. Code Apache-2.0; this text CC BY 4.0.*

This note pins down the single number the whole project is about, states the
fair-coin baseline exactly, and records a prior-art check (done 2026-10-14) of
what the literature does and does not settle. No new numeric claim is made
about the asymptotics; that waits for the computation (see "Pending", below).

## The headline number, stated exactly

Gupta & Sridhar, *A Peculiarity in the Parity of Primes* (arXiv:1812.11841),
report that among the **first 50,000,000 primes**, exactly **25,032,384** have
an **odd** base-10 digit sum. Everything else follows by arithmetic:

| quantity | value |
|---|---|
| primes counted, `n` | 50,000,000 |
| odd digit sum | 25,032,384 |
| even digit sum | 24,967,616 |
| excess (odd − even) | **64,768** |
| fair-coin std, `sqrt(n)` | 7071.0678… |
| z-score, excess / sqrt(n) | **≈ 9.16 σ** |

So "about 10× farther from the mean than a fair coin" is, precisely, **9.16
standard deviations**. A fair coin essentially never does this (two-sided
p ≈ 10⁻²⁰). The bias at `n = 5×10⁷` is unmistakably real. `stats.py`
reproduces every cell above from the raw count; `test_stats.py` checks it.

## What is settled, and what is not (prior-art check, 2026-10-14)

- **Settled (limit):** Mauduit & Rivat, *Sur un problème de Gelfond: la somme
  des chiffres des nombres premiers* (Ann. of Math. 171 (2010)), prove the
  base-10 prime digit sum is equidistributed mod 2 — asymptotic density 1/2.
  So the *relative* bias → 0. This is an asymptotic density with an upper
  bound on the error; it does **not** predict a finite-sample z of ~9.
- **Not what we need (2026 work):** two recent papers on digit sums of primes
  do not touch this. arXiv:2606.04677 (*An Explicit Surjectivity Threshold for
  Digit Sums of Primes*) makes the Drmota–Mauduit–Rivat constants explicit to
  bound where every residue class coprime to 9 is hit — a coverage result, not
  a parity-bias result. arXiv:2607.10654 (*The Prime Digit Distribution
  Conjecture*) proves an **averaged individual-digit** frequency → 1/10 and
  explicitly flags the digit-**sum** parity question as logically distinct and
  out of its scope. Neither cites arXiv:1812.11841.
- **Open, quantitatively:** no work we found predicts or explains the *size* of
  the finite-sample excess, nor says whether `z(x)` for primes up to `x` stays
  near 9, grows (a genuine one-signed secondary term), or shrinks (the 2018
  figure being partly an artifact of stopping at a fixed prime count). That is
  the question this repository measures.

## An honest credit, and what our Fact 4 adds

Gupta & Sridhar already note the bias **survives removing the last digit** of
each prime. Our **Fact 4** (see `core.py`, `NOTES.md`) is a proof of *why*:
every prime `p > 5` ends in a digit from {1,3,7,9}, all odd, so the parity of
`p`'s digit sum is exactly `1 XOR` the parity of the digit sum of `p // 10`.
The entire base-10 bias is therefore carried by the "beheaded" integers
`floor(p/10)`. We claim the clean identity and its consequences, not the
empirical observation itself, which is theirs.

## Pending computation (no verdict yet)

Two jobs feed the next version:

1. `prime-parity-trajectory-v1` — a fast, pure-stdlib probe to `x = 10⁸`:
   odd-minus-even excess and `z = excess / sqrt(π(x))` at value cutoffs
   `x = 10^1 … 10^8` for even bases {2,4,6,8,10,12,16}; an empirical check of
   Fact 4 on every prime; and an empirical check of Fact 6 (integers coprime to
   10 have base-10 parity excess bounded by 8). Base 2 must show the odious
   (more-1s) lead as a sanity control (OEIS A130911).
2. `prime-digit-parity-v3` — the larger run to `x = 10⁹…10¹⁰`, covering the
   actual `5×10⁷`-prime regime and beyond.

**Pre-registered reading** (from `RESOLUTION.md` §5): if `z(10^k)` hovers O(1)
around the ~9 level with sign structure tied to the cutoff convention, the 2018
figure is partly a single-cutoff effect; if `z(10^k)` is one-signed and grows,
there is a real secondary term worth isolating (compared against a
coprime-to-10 random null in the same range). We will write `RESULTS.md` and
tag a numeric release only once the sandbox returns these tables.
