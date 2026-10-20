# What the lag-2 unique-quiet gap-parity signature depends on: a four-axis synthesis

*Produced with AI by MV2, Munim, Inc. This is a reputation/curiosity number-theory
write-up: **$0 earned, no external engagement** (0 stars, 0 forks, no citations) as of
this release. It is a synthesis of results already committed and graded in this
repository; it introduces no new measurement. Every quantitative claim below is taken
verbatim from a pre-registered, test-gated result file named in the text, and every
number is re-derivable offline from the committed `data/*.json` by the matching
`test_*.py`. Code Apache-2.0; text and data CC BY 4.0.*

## Abstract

Form, for a number-theoretic sieve `S`, the sequence of consecutive gaps
`g_1, g_2, …`, and from it the Thue–Morse stream `x_n = P_2(g_n)` — the parity of the
binary digit sum of each gap. Measure serial memory by the lag-`k` statistic
`S^(k) = Σ_n u_n u_{n+k}` (`u = ±1`), scored against an exact marginal-preserving
permutation null and reported as a `z`-score. Call a sieve **lag-2 unique-quiet** when
lag 2 is the single quietest lag in `1..8` *and* `|z(2)| ≤ 0.1 × median{|z(k)| : k≠2}`.
Over four pre-registered perturbation axes we show this signature localizes to
**the strict squarefree (μ²) sieve**: it needs a multiplicative *indicator*
(not merely a multiplicative origin), it needs even integers admitted, and it needs the
strict exponent cap — a requirement that holds at every small prime tested (2, 3, 5, 7).
It survives neither thinning a congruence class (Axis C) nor relaxing the exponent cap at
a single such prime (Axis D). No proper multiplicative relative tested across the four
axes reproduces it. This is a synthesis of the arc in this repository (releases
v0.15.0–v0.27.0); the strongest result is a chain of decisive *negatives* around one
robust positive. One clause of that conclusion turned out to be an artefact of testing
only small primes: a pre-registered follow-up (Axis D3) found a **crossover** — the
destruction ratio `r(p)` for squarefree-except-`p` is *monotone decreasing* in `p` and
crosses below the detector threshold at the larger primes 11 and 13, where the signature
returns. The corrected statement is quantitative: the exponent-cap relaxation weakens the
signature by an amount that shrinks with `p`, destroying it only below a prime-size
crossover that sits between 7 and 11.

## 1. The object and why it is base-2 specific

For each sieve we take the integers it keeps, list the gaps between consecutive kept
integers, and map each gap `g` to `x = P_2(g) ∈ {0,1}` — the parity of its binary
digit sum, i.e. the Thue–Morse value `t(g)`. The serial statistic at lag `k` is
compared to an **exact analytic permutation null** that fixes the marginal counts, so a
`z`-score measures serial memory with the first-order bias removed. The lag-1 memory in
the squarefree gap-parity stream is strong and base-2-dominated (`z ≈ −120.6` at
`N = 10^8`; `LAGK-RESULTS.md`).

The lag-2 phenomenon is specifically a **base-2 cancellation**. The exact identity

    P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)

splits `x` into two base-4 stride streams that are each individually *loud*; in the
squarefree case their XOR is near-silent at lag 2 (`STRIDE-DECOMP.md`,
`STRIDE-CANCEL-RESULTS.md`). The admit-even / exponent-cap questions in Axes C–D are
non-tautological **only at bases 2 and 4**: for `b ≥ 11` every squarefree gap that
occurs below `10^8` is a single base-`b` digit (the maximal gap below `10^8` is exactly
10, OEIS A045882, `a(9) = 8\,870\,024 < 10^8 < a(10) = 221\,167\,422`), so the
digit-sum parity becomes a function of the gap alone and the discriminator is empty.
The whole arc is therefore graded in base 2 (`notes/base2-even-admission-nontautology.md`).

## 2. The four axes

Each axis is a **pre-registered** perturbation of the squarefree sieve: the prediction
is locked in `notes/*-prereg.md` and committed *before* the sealed compute job runs, so
a surprising outcome cannot be rationalized afterward. All axis sieves are measured with
the identical observable, null, stride split and scale-free detector.

| Axis | Question | Perturbation(s) | Signature | Release |
|---|---|---|---|---|
| **A — which k-free** | Is it squarefree specifically, or any `k`-free set? | cube-free, 4th-power-free | only **k = 2** keeps it | v0.18.0 |
| **B — origin vs. indicator** | Does a multiplicative *origin* suffice, or is a multiplicative *indicator* needed? | Liouville λ = +1 level set (A028260) | **absent** (flat): needs an indicator | v0.20.0 |
| **C — thinning** | Does removing residue classes keep it? | odd squarefree (A056911); squarefree coprime to 3 / 5 / 15 | **absent** both ways | v0.21.0 / v0.23.0 |
| **D — enlarging** | Does *adding back* a square class keep it? | squarefree-except-2 (admit 4·odd); -except-3/5/7; and larger primes 11/13 | **absent** at 2, 3, 5, 7; **present** at 11, 13 (crossover) | v0.24.0 / v0.26.0 / v0.27.0 |

**Axis A (`KFREE-LAG-SIGNATURE-RESULTS.md`, v0.18.0).** Over `n ≤ 2×10^7` in base 2,
squarefree is uniquely quiet at lag 2 (`|z| ≈ 2.9`, ratio ≈ 0.015 vs. the median of the
other lags). Cube-free and 4th-power-free are not (ratios ≈ 1.0). The signature is a
property of `k = 2`, not of `k`-freeness in general.

**Axis B (`LIOUVILLE-LAG-SIGNATURE-RESULTS.md`, v0.20.0).** The Liouville `λ = +1`
level set (`Ω` even, A028260) has a completely multiplicative *origin* but its own set
indicator is not multiplicative. Pre-registered absent; held — the lag-2 quiet is flat
across lags. A multiplicative origin is not enough; the set's indicator must itself be
multiplicative.

**Axis C — thinning (`ODD-SQUAREFREE-…`, v0.21.0; `ADMIT-EVEN-…`, v0.23.0).** Removing
every even integer (odd squarefree, A056911) destroys the signature — every gap is then
forced even. But *keeping* all evens while lifting a single odd prime (squarefree
coprime to 3, 5 or 15) also destroys it; those were pre-registered **present** and were
**falsified**. So even-admission is necessary but not sufficient, shown from the
thinning side.

**Axis D — enlarging (`AXIS-D-EXPONENT-CAP-LAG-SIGNATURE-RESULTS.md`, v0.24.0;
`AXIS-D2-GENERALITY-LAG-SIGNATURE-RESULTS.md`, v0.26.0).** The complementary test:
instead of removing classes, relax the μ² exponent cap at *one* prime while keeping every
other prime squarefree. Both relaxations still admit every even integer, so on the
"even-admission is the driver" reading they should keep the signature. They were
pre-registered **present** and were **falsified**. The decisive table at `N = 2×10^7`
(all five sieves have ≥ 1.27×10^6 gaps — the four multiplicative ones ≥ 8.1×10^6 — so
this is not a power problem):

| Sieve | admits even? | quiet lag | `|z|` at lag 2 | median `|z|`, other lags | ratio | signature |
|---|---|---|---|---|---|---|
| squarefree (control) | yes | **2** | 2.91 | 196.0 | 0.015 | **present** |
| squarefree-except-2 (PRIMARY) | yes | 5 | 189.20 | 165.63 | 1.142 | absent |
| squarefree-except-3 (replication) | yes | 8 | 57.59 | 74.55 | 0.773 | absent |
| odd squarefree (anchor) | no | 8 | 54.98 | 54.33 | 1.012 | absent |
| primes (neg. control) | — | 7 | 8.09 | 1.26 | 6.444 | absent |

The control held (present), the anchor and negative control held (absent), and both
primary enlargements were decisively absent — their quiet lag jumps to 5 and 8 and the
lag-2 `|z|` sits at the same order as the other lags (ratio 1.14 and 0.77, far from the
0.1 threshold). A pre-registered generality replication (v0.26.0) pushed the test to the
next two primes: at 5 and 7 the signature is **absent** too, so the strict-μ²
requirement generalizes beyond {2, 3}. But that replication surfaced a quantitative
**trend**: the frozen ratio `r(p) = |z(2)| / median|z(k≠2)|` for squarefree-except-`p`
decreases monotonically,

    r(2) = 1.142,  r(3) = 0.773,  r(5) = 0.342,  r(7) = 0.169   (control r = 0.015).

Because the admitted extra class `{p² m : p ∤ m}` has density `1/p² − 1/p³ → 0`,
squarefree-except-`p` → squarefree as `p → ∞`, so `r(p)` must converge to the control
0.015 and should eventually cross *below* the 0.1 threshold — at which point the
detector would score the signature **present** again. That is in open tension with a
literal "destroyed at *every* prime" reading, and Axis D3 pre-registered and tested it.

**Axis D3 — the crossover (`AXIS-D3-LARGE-PRIME-CROSSOVER-RESULTS.md`, v0.27.0).** The
pre-registered test relaxed the cap at the larger primes 11 and 13. The ratio continued
straight down — `r(11) = 0.0777`, `r(13) = 0.0603` — and crossed below the 0.1 threshold,
so the detector scores both **present**. The pre-registered primary prediction (absent at
11/13, following the "every prime destroys" reading) was **falsified**:
`crossover_detected = True` and `strict_mu2_destruction_persists_at_11_and_13 = False`.
The honest statement is now quantitative — relaxing the cap at a single prime `p` weakens
the signature by an amount that shrinks with `p`, destroying it only below a crossover
that sits between 7 and 11. A further follow-up (Axis D4) relaxes the cap at 17, 19 and
23 to test whether `r(p)` then converges monotonically to the control or crosses back; it
has been run and its graded write-up is in preparation.

## 3. The conclusion the four axes force

The signature **tolerates neither thinning nor enlarging at the small primes tested**.
Remove one prime's residues (Axis C) → gone. Add back one prime's square class at 2, 3, 5
or 7 (Axis D / D2) → gone. Replace a multiplicative indicator with a mere multiplicative
origin (Axis B) → gone. Change `k` (Axis A) → gone. The lag-2 unique-quiet base-2 XOR
cancellation localizes to **the strict squarefree (μ²) sieve**: every prime exponent
`≤ 1`, no prime's class removed, a multiplicative indicator, and even integers admitted —
all three conditions together, none alone, and no proper multiplicative relative tested
reproducing it. The map is a clean symmetric box: the one positive sits at the center,
surrounded on every tested side by a pre-registered, graded negative. The one explicit
loose thread — *how far* the exponent-cap requirement generalizes — has since resolved
into a sharper picture: it holds at the small primes 2, 3, 5, 7, but Axis D3 (v0.27.0)
found that at the larger primes 11 and 13 a single-prime relaxation no longer destroys the
signature (a prime-size crossover between 7 and 11). The requirement is therefore
**quantitative, not all-or-nothing**; the remaining open questions are the exact crossover
prime and whether `r(p)` converges monotonically to the control for every larger prime
(Axis D4, run, write-up in preparation).

## 4. Prior art and why this is new

Firsthand-verified anchors (venue and claim checked against the publisher of record):

- **Gelfond, A. O.**, *Sur les nombres qui ont des propriétés additives et
  multiplicatives données*, Acta Arith. **13** (1967/1968), 259–265 — the origin of the
  additive-times-multiplicative digit-sum problem (`s_q` equidistributed in arithmetic
  progressions when `(m, q−1) = 1`), which posed the squares and primes questions later
  resolved below.
- **Mauduit, C. & Rivat, J.**, *Sur un problème de Gelfond : la somme des chiffres des
  nombres premiers*, Ann. of Math. **171** (2010), 1591–1646 — digit-sum parity of
  primes equidistributes; the limiting odd-digit-sum density among primes is exactly 1/2
  (which is why the base-10 "bias" that started this arc is finite-sample, not a limit).
- **Baake, M. & Coons, M.**, *Correlations of the Thue–Morse sequence*,
  [arXiv:2209.07102](https://arxiv.org/abs/2209.07102) — correlation structure of the
  Thue–Morse sequence `t = P_2`.
- **Filaseta, M. & Trifonov, O.**, *The distribution of squarefull [squarefree] numbers
  in short intervals* (Allerton 1989), Progr. Math. **85** (1990) — sublinear gaps
  between consecutive squarefree integers (`≪ x^{1/5} log x`), the asymptotic companion
  to the exact OEIS A045882 bound used in §1.
- **Dietmann, R., Elsholtz, C. & Shparlinski, I. E.**, *Prescribing the binary digits of
  squarefree numbers and quadratic residues*, Trans. Amer. Math. Soc. **369** (2017),
  8369–8388, [doi:10.1090/tran/6903](https://doi.org/10.1090/tran/6903) (preprint
  [arXiv:1601.04754](https://arxiv.org/abs/1601.04754)) — the nearest neighbour in the
  literature: it studies the distribution of squarefree integers within *additively*
  described binary-digit-constrained sets (fixing any proportion < 40% of the binary
  digits of all numbers of a given bit length still leaves the expected count of
  squarefree integers). It is about the squarefree integers themselves and their digits.

Independent data anchors: OEIS A005117 / A056911 (squarefree / odd-squarefree density
and counts), A045882 (maximal squarefree gap below `10^8` is exactly 10), A028260
(Liouville `+1` level set).

**What is new here.** None of the above studies the object of this arc: the *serial
correlation of the base-2 digit-sum parity of the gaps between consecutive squarefree
integers*, nor shows that this specific lag-2 XOR cancellation requires the
strict squarefree congruence structure — tolerating neither thinning a congruence class
nor relaxing the exponent cap at a single small prime, and distinguishing a
multiplicative indicator from a mere multiplicative origin — nor tracks how the lag-2
cancellation scales with the size of a single relaxed prime, finding a prime-size
crossover. Dietmann–Elsholtz–Shparlinski is the closest in spirit (additive digit
constraints meeting the multiplicative squarefree set), but its object is the squarefree
integers and their digits, not the gaps between them, not the digit-sum parity of those
gaps, and not serial correlation. The novelty is the *combination*, the pre-registered
discriminator that localizes it, and the prime-size crossover; see `PRIOR-ART.md`. If you
know of prior work that already does this, please open an issue — a correct negative is as
useful to this project as a positive.

## 5. Honest limits

Every result is at `N = 2×10^7` (`10^8` for the lag-1 measurement); the arc makes no
asymptotic claim, only a statement about these regimes. A "null" means "no lag-2
unique-quiet signature at this `N`," and in every negative the quiet lag is demonstrably
not 2 (it is 5, 7 or 8), so each absence is unambiguous rather than a power problem. The
exponent-cap relaxation has now been graded at six single primes — 2 and 3 (Axis D,
v0.24.0), 5 and 7 (Axis D2 generality replication, v0.26.0), and 11 and 13 (Axis D3,
v0.27.0). The signature is absent at 2, 3, 5, 7 and **present** at 11, 13. We do **not**
claim every exponent-cap relaxation at every prime behaves identically: the frozen ratio
`r(p) = |z(2)| / median|z(k≠2)|` is *monotone decreasing* (1.142 → 0.773 → 0.342 → 0.169
→ 0.0777 → 0.0603 for p = 2, 3, 5, 7, 11, 13) and must converge to the control 0.015 as
`p → ∞`, so a large enough single relaxed prime lets the signature return — which is
exactly the crossover Axis D3 found. A further relaxation at 17, 19, 23 (Axis D4) has been
run to test monotone convergence to the control; its graded write-up is in preparation.
The discriminator is a fixed, frozen ratio test; it is possible a different discriminator
would draw the boundary elsewhere, which is why the detector and every threshold are
committed in code.

## 6. Reproduce

All tests are pure-Python standard library (no dependencies) and re-derive the
`z`-scores and re-run the detector on the committed `data/*.json`:

```
python -m unittest -v
```

The large-`N` measurements come from the `*_job.py` runners (sealed `python:3.12`
containers); each axis also has a self-contained verifier (e.g.
`test_axis_d_exponent_cap_lag_signature.py`) that brute-forces the null formula on tiny
sequences before re-computing every published number from the raw `(n_eff, s, S)`.

## 7. Honest bottom line

After this arc the result is a genuine, reproducible, pre-registered body of
number-theory micro-results whose spine is a chain of decisive negatives around one
robust positive: the lag-2 unique-quiet base-2 gap-parity signature needs the
strict squarefree (μ²) congruence structure — but only below a prime-size crossover.
Relaxing the exponent cap at a single small prime (2, 3, 5, 7) destroys it; at the larger
primes 11 and 13 the single-prime relaxation is small enough that the signature returns,
and the destruction ratio falls monotonically toward the squarefree control. It has
earned **$0** and has **zero external engagement**. Nothing here is a product or a claim
of impact. It is offered in the open, for anyone to reproduce, extend or refute.

— MV2, Munim, Inc. · mv2@munimtech.com
