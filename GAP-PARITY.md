# Digit-sum parity of the prime-GAP sequence

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.*

## The question

Every published result on digit-sum parity and the primes studies the parity of
the digit sum **of the primes themselves** — Gelfond's problem, proved for the
average by Mauduit–Rivat (*Annals* 2010): primes split evenly between even and
odd decimal digit sum. A separate thread studies the *last digits* of
consecutive primes (Lemke Oliver–Soundararajan, arXiv:1604.02443).

Here we ask a question that, as far as three fresh prior-art searches (2026-10-15,
2026-10-16 ×2) could find, **nobody has asked**: what is the digit-sum parity of
the prime *gaps* g_n = p_{n+1} − p_n? We compute, for bases b ∈ {2,4,6,8,10,16},

    E_b^gap(M) = #{n<M : χ_b(g_n) = −1} − #{n<M : χ_b(g_n) = +1} = − Σ_{n<M} χ_b(g_n),
    χ_b(m) = (−1)^{s_b(m)},   s_b = base-b digit sum,

over all gaps between consecutive primes ≤ 10^6, 10^7, 10^8, with z = E/√(n_gaps).
χ_2 is the Thue–Morse value of the gap.

## Why this is new

- Mauduit–Rivat (primes 2010, squares 2009) and Spiegelhofer (cubes,
  arXiv:2308.09498) evaluate the digit sum **along** a subsequence of the
  integers; the gap sequence g_n is a *derived* sequence, not a subsequence.
- arXiv:2102.01018 studies gaps *within the Thue–Morse word* — a different object
  (positions, not values).
- Rudin–Shapiro along primes (arXiv:2306.10629) and TM-along-primes
  (arXiv:2001.07017, 1803.01689) stay on the primes themselves.
- This repo's own `cross-sequence` result (v0.6.0) compares primes vs squarefree
  vs semiprime vs all-integers — again subsequences of ℕ, not gaps.
- **OEIS cross-check:** the closest catalogued objects are all distinct
  from χ_b(g_n): **A001223** (the gaps g_n themselves), **A028334** (half-gaps
  g_n/2 — relevant only to our base-2 identity below), **A255311** (the decimal
  *concatenation* of gap sizes, a constant, not a parity statistic), **A372429**
  (binary digit indices of prime(n), i.e. of the primes, not the gaps), and
  **A010060** (the Thue–Morse sequence itself). No OEIS entry gives the running
  parity excess of digit sums of prime gaps in any base.

No source computes χ_b(g_n). → novel.

## Anchors (gate — all must pass or the number is a bug, not a result)

Since there are no published numeric anchors for this new sequence, every run is
gated on *exact, hand-derived* identities that validate the method:

1. **π(x)** = 78498 / 664579 / 5761455 at 10^6/10^7/10^8 (sieve correct).
2. **Telescoping:** Σ g_n = p_M − 2 at each cutoff (gap extraction correct).
3. **Single odd gap:** exactly one odd gap exists (g_1 = p_2 − p_1 = 1); all
   later gaps are even. So Σ_n [g_n odd] = 1 in every prefix.
4. **Half-gap Thue–Morse identity:** for even g, s_2(g) = s_2(g/2), so
   χ_2(g_n) = χ_2(g_n/2). A second, independent E_2 tally over half-gaps
   (the A028334 sequence) must equal the direct one.
5. **Two-routine digit guard:** repeated-division parity == explicit-digit parity
   for the first 200 000 gaps, every base.

`test_gap_parity.py` re-derives anchors 1–5 on primes < 10^6 with its own sieve,
stdlib only. The measurement job reported, before writing any number:

> `gate PASSED (pi(x), telescoping, single odd gap, half-gap identity, two-routine guard).`

The half-gap identity (anchor 4) held exactly: the direct base-2 excess and the
independent half-gap tally agree to the digit at every cutoff
(−7 563 / −77 744 / −756 388).

## Result

Gate: **PASSED** (quoted above; from the sealed CPU-sandbox run log).

E_b^gap and z = E/√(n_gaps), by base and cutoff:

| base | E @10^6 | z @10^6 | E @10^7 | z @10^7 | E @10^8 | z @10^8 |
|-----:|--------:|--------:|--------:|--------:|--------:|--------:|
| 2  | −7 563 | −26.99 | −77 744 | −95.37 | −756 388 | −315.12 |
| 4  | 16 607 | +59.27 | 132 234 | +162.21 | 1 051 924 | +438.25 |
| 6  | 3 031 | +10.82 | 31 814 | +39.03 | 341 536 | +142.29 |
| 8  | −11 497 | −41.04 | −69 764 | −85.58 | −437 046 | −182.08 |
| 10 | −15 033 | −53.66 | −96 440 | −118.30 | −665 788 | −277.38 |
| 16 | −38 231 | −136.46 | −264 548 | −324.51 | −1 901 898 | −792.36 |

(n_gaps = π(cutoff) − 1 = 78 497 / 664 578 / 5 761 454.)

## Which hypothesis the data support (decided from the z-column; prereg run 25)

Three hypotheses were registered **before** the numbers were seen
(`notes/prime-gap-parity-prereg.md`):

- **H0 (null):** gap parities behave like a fair ±1 walk, |z_b| = O(1) for all b.
  Registered as the *most likely* outcome, on the reasoning that the half-gaps
  g_n/2 have no obvious parity bias.
- **H1 (prime-specific valley):** the base-8 near-cancellation that the primes
  show (E_8 along the primes at 10^8 = −1329, z ≈ −0.55, three to four orders of
  magnitude below every other base) is a property of the *primes*; on gaps it
  disappears or moves.
- **H2 (2^k-subsample mechanism):** the base-8 feature reappears on gaps ⇒ it is
  about power-of-two subsampling of a structured even sequence, not primeness.

**Verdict: H1 — and H0 is decisively rejected.** The prereg's expected null is
wrong: the gap sequence carries a *large and growing* digit-sum-parity bias in
every base (|z| climbing from tens at 10^6 to 142–792 at 10^8, all monotone in the
cutoff — the opposite of an O(1) fair walk). Within that, base 8 is **not**
special: its |z| at 10^8 is 182.08, mid-pack and in fact larger than base 6
(142.29), while the extreme base is 16 (|z| = 792.36). Nothing resembles the
primes' near-zero base-8 valley (|z_8| ≈ 0.55 along the primes). So the base-8
cancellation **does not transfer** to the gaps — it is a property of the primes
themselves (H1), not of 2^k subsampling of any structured even sequence (H2,
rejected: base 8 is not the anomalously small base; base 6 is the smallest and no
base cancels).

## Why the gaps are so biased (mechanism, honestly)

The large signal is not mysterious, and it is weaker evidence of deep structure
than the prime anomaly it was designed to probe. All gaps past g_1 are even, and
the gap-size distribution is sharply non-uniform — concentrated on small even
values and multiples of 6 (g ∈ {2, 4, 6} alone cover a large fraction, and 6 is
the single most common gap in this range). Those favored sizes have digit-sum
parities that do **not** cancel: e.g. in base 10, gaps 2/4/6/8 have *even* digit
sum while 10/12/14/16/18 have *odd*, and the former are far more frequent, which
is exactly the direction of the measured negative base-10 excess. So the base-b
excess of the gaps is, to first order, a weighted readout of the gap-size
histogram, not an independent arithmetic coincidence. The scientific content is
the **contrast**: the primes' hallmark base-8 cancellation is absent on this
derived sequence, confirming it is prime-specific rather than an artifact of the
base-8 digit machinery acting on any structured even sequence.

## Honest limits

- The magnitude of E_b^gap largely reflects the (well-known, non-uniform)
  gap-size distribution; we claim the *measurement* and the *H1/H0 contrast*, not
  that the gaps hide a new arithmetic law beyond their size histogram.
- H0 was pre-registered as most likely and is reported here as decisively
  **wrong** — a pre-registered prediction that failed, stated plainly.
- Measured to 10^8 (5.76 M gaps); the trend beyond is not claimed. z grows with
  the cutoff, consistent with a non-zero mean χ_b(g_n) rather than diffusive
  cancellation, but an asymptotic rate is not asserted.
- Single machine, single sieve; the two-routine guard + half-gap identity are the
  internal cross-checks that stand in for an external anchor on this new sequence.

## Reproduce

```
python -m unittest -v test_gap_parity   # anchors on primes < 1e6, stdlib only
```

Full measurement: `job.py --out-dir out --limit 100000000` (segmented stdlib
sieve; CPU-only; 67 s on 4 cores in the sealed sandbox). Data:
`data/gap_parity_excess.csv` (CC BY 4.0). Code Apache-2.0.

— MV2, Munim, Inc.
