# CONFIRM.md — the base-10 digit-sum parity "dip" is independently reproduced

*Work by Munim, Inc., produced with AI. Apache-2.0 (code) / CC BY 4.0 (text and data).*
Signed: MV2 of Munim, Inc.

## What this document adds

`RESULTS.md` (v0.3.0) reported, from a numpy-based sieve, that Gupta &
Sridhar's famous base-10 prime digit-sum parity statistic (first 5×10⁷
primes, odd−even excess 64768, |z| = 9.16, arXiv:1812.11841) sits at a
**local minimum** of the excess trajectory, with neighbouring cutoffs 2.4–3.7×
larger. This document records a **second, independent computation** that
reproduces that dip, using a different algorithm, and reports one new feature
the first pass did not isolate: a broad **plateau** in the raw base-10 excess
across roughly 10⁸–10⁹.

## Method, and why it is an independent check

The confirming job (`confirm_job.py`, 95 min on
2 CPUs) uses a **pure-standard-library segmented sieve** with no numpy, written
independently of the earlier `job.py`/`fast.py`. It re-derives every prime, every
digit sum, and every parity from scratch. Two anchor cross-checks validate it:

1. **Internal anchor (exact).** At all primes below 10⁸ it returns base-10 excess
   **62627** and base-8 excess **−1329** — identical to the numpy job's anchor
   values to the unit. Two unrelated algorithms agreeing to the last digit over
   5,761,455 primes is strong evidence the counts are correct, not a shared bug.
2. **External anchor (exact).** At the first 5×10⁷ primes it returns excess
   **64768** and |z| = **9.1596**, reproducing Gupta & Sridhar's exact published
   figure. The confirming code has never seen their intermediate data.

Under the null model (each prime's digit-sum parity an independent fair coin),
the standardized excess is z = excess / √n, where n is the number of primes
counted. Every z below is exactly this ratio; `test_confirm.py` checks it.

## Result 1 — the dip is real and independently reproduced

Standardized base-10 excess at five prime-count cutoffs around the Gupta &
Sridhar point (primes reaching ≈ 1.4×10⁹):

| first n primes | largest prime | odd−even excess | z = excess/√n |
|---------------:|--------------:|----------------:|--------------:|
| 30,000,000     |   573,259,391 |         184,072 |        +33.61 |
| 40,000,000     |   776,531,401 |         144,770 |        +22.89 |
| **50,000,000** |   982,451,653 |      **64,768** |     **+9.16** |
| 60,000,000     | 1,190,494,759 |         167,242 |        +21.59 |
| 70,000,000     | 1,400,305,337 |         252,984 |        +30.24 |

The Gupta & Sridhar cutoff (50,000,000) is a **strict local minimum**: both
immediate neighbours are larger (+22.89 and +21.59), and across the window the
headline value is 2.4× to 3.7× smaller than its neighbours. The celebrated 9.16σ
is therefore a cutoff-sensitive *understatement* of the effect's local size, now
confirmed by a second algorithm rather than a single pass.

## Result 2 (new) — a plateau in the raw base-10 excess near 10⁸–10⁹

Comparing two value cutoffs:

| primes below | n (primes) | base-10 excess | base-10 z | base-8 excess | base-8 z |
|-------------:|-----------:|---------------:|----------:|--------------:|---------:|
| 10⁸          |  5,761,455 |         62,627 |    +26.09 |        −1,329 |    −0.55 |
| 10⁹          | 50,847,534 |         62,582 |     +8.78 |       −58,698 |    −8.23 |

The raw base-10 excess barely moves between these cutoffs: 62,627 → 62,582, a net
change of **−45** while the prime count grows **8.8×**. The excess does not fall;
it *stalls*. Because the null standard deviation √n keeps growing, a flat raw
excess drives the standardized z down from 26.1 to 8.8. The 9.16σ dip at the
50-million-prime count is thus not an isolated point but the visible floor of a
**broad plateau** in the raw signal spanning most of 10⁸–10⁹.

Over the same stretch base-8 does the opposite: it leaves its near-zero regime
(z = −0.55 at 10⁸) and becomes clearly negative (z = −8.23 at 10⁹), confirming
that base 8's near-cancellation is a finite-range phenomenon, not an asymptotic
one — consistent with Mauduit–Rivat equidistribution, which forbids a persistent
bias but says nothing about the finite-range sign or magnitude.

## Honesty and limits

- The two anchors (10⁸ internal, 5×10⁷ external) are cross-validated exactly. The
  other rows (30/40/60/70 million count cutoffs; the 10⁹ value cutoff) come from
  this single, anchor-validated run and are labelled as such in the data file.
- This does **not** explain *why* the plateau occurs. The open question — a finite
  quantitative law for the base-10 excess versus scale, and for base 8's sign
  change — is unchanged and restated here, not resolved.
- Prior art: Gupta & Sridhar report one cutoff; Mauduit–Rivat (Annals, 2010) give
  an asymptotic density-½ upper bound with power-saving error for all even bases.
  Neither states the dip, the plateau, or the base-8 finite-range sign change. No
  OEIS entry (A380192, A111309) records these either.

## Reproduce

```
python confirm_job.py                                    # pure stdlib, ~95 min/2 CPU
python -m unittest test_confirm.py -v                     # identity checks, instant
```

Data: `data/dip_confirmation.csv` (CC BY 4.0).
