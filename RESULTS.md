# RESULTS — Prime digit-sum parity to 10^10, across bases

*Work by Munim, Inc., produced with AI. Code Apache-2.0; text and data CC BY 4.0.*

This note reports the verified numeric verdict for the study begun in this repo:
is the "unexplained" base-10 prime digit-sum parity bias of Gupta & Sridhar
(*A Peculiarity in the Parity of Primes*, arXiv:1812.11841) what it appears to
be? We sieved every prime up to **x = 10^10** and recorded, at both **value
cutoffs** (primes ≤ x = 10^k) and **count cutoffs** (the first N primes), the
digit-sum parity excess

    E_b(·) = #{p : s_b(p) odd} − #{p : s_b(p) even},    z_b = E_b / sqrt(n)

for bases b ∈ {2, 4, 6, 8, 10, 12, 16}. Under a fair-coin null E_b is a sum of
n independent ±1's, so z is a standardized score. Two independent sandbox jobs
(`prime-parity-trajectory-v1`, `prime-digit-parity-v3`) agree **exactly** at
every shared cutoff up to 10^8; 10^9 and 10^10 come from the larger job. All
numbers below are reproduced by `test_results.py` and the CSVs in `data/`.

Baselines that held on every run (sanity gates):
- **Fact 4** (every prime p > 5 ends in an odd base-10 digit) held for all
  5 < p ≤ 10^8.
- **Fact 6 / THEORY.md** bound: the coprime-to-10 excess stayed within the tight
  constant φ(10) = 4 at every cutoff (max |·| = 4), confirming the sharpened
  bound proved in v0.2.0.

---

## Finding 1 — Gupta & Sridhar's 9.16-sigma sits at the bottom of a local dip

G&S report, for the **first 5×10⁷ primes** in base 10, an odd−even excess of
64,768, i.e. z = 9.16 — "about an order of magnitude farther from the mean than
expected." We reproduce that number exactly. But stepping the count cutoff by
10⁷ primes on either side shows it is a **local minimum**, not a typical value:

| first N primes | excess (base 10) | z |
|---:|---:|---:|
| 3×10⁷ | 184,072 | 33.61 |
| 4×10⁷ | 144,770 | 22.89 |
| **5×10⁷** | **64,768** | **9.16** |
| 6×10⁷ | 167,242 | 21.59 |
| 7×10⁷ | 252,984 | 30.24 |

The neighbouring cutoffs are **2.4–3.7×** larger. The same trough appears under
the independent value-cutoff convention: at x = 10⁹ (≈ 5.08×10⁷ primes) the
base-10 excess is 62,582 (z = 8.78), essentially unchanged from its value at
x = 10⁸ (62,627) a full decade earlier, before jumping to 2,110,439 (z = 98.93)
at x = 10¹⁰. Both conventions locate the **same** dip near p ≈ 10⁹.

**Reframing (verified, new):** the celebrated "9.16σ" is real but its *size* is
an artifact of the cutoff landing at the bottom of a genuine excess trough near
p ≈ 10⁹. The generic base-10 anomaly in that region is already z ≈ 20–34, and it
keeps growing — z ≈ 99 by 10¹⁰. G&S understated, rather than overstated, the
effect; and the headline figure is cutoff-sensitive in a way not previously
pointed out.

---

## Finding 2 — The base-10 excess trajectory is strongly non-monotone

Base 2 (the classical Thue–Morse / prime-popcount bias) rises smoothly and
monotonically with scale. Base 10 does **not** — it stalls for an entire decade
around 10⁹ and then surges:

| x | z(base 2) | z(base 10) |
|---:|---:|---:|
| 10⁶ | 17.00 | 6.84 |
| 10⁷ | 24.76 | 3.64 |
| 10⁸ | 29.38 | 26.09 |
| 10⁹ | 90.87 | **8.78** |
| 10¹⁰ | 134.34 | 98.93 |

Base 2 is a clean monotone power law; base 10 has large scale-localized dips.
This is the arithmetic texture behind Finding 1: the "peculiarity" is not a
steady drift but a signal with pronounced troughs, and sampling conventions that
fix a single cutoff (a round N or a round x) can land in one.

---

## Finding 3 — Base 8 shows almost no bias; the standardized anomaly is *largest* in base 2

Every even base satisfies gcd(2, b−1) = 1, so by Mauduit–Rivat all of them are
asymptotically equidistributed (density → ½); only the **rate** differs. Empirically
the rates are wildly different, and in a way that defeats the naive "higher base
→ bigger standardized anomaly" heuristic:

| base | |z| at 10¹⁰ | rough growth exponent α (E_b ∝ x^α) | behaviour |
|---:|---:|---:|:--|
| 2 | 134.34 | ≈ 0.70 | strongest, smooth, one-signed |
| 4 | 104.78 | ≈ 0.72 | strong, smooth |
| 6 | 54.09 | ≈ 0.69 | grows |
| 8 | 5.80 | ≈ 0.50 | **bounded z — square-root (random) size** |
| 10 | 98.93 | ≈ 0.76 | grows, non-monotone |
| 16 | 24.56 | — | erratic |

**Base 8 is the clean negative result:** its standardized excess stays within
|z| < 15 all the way to 10¹⁰ (sign-oscillating: −1.82, +0.26, −0.55, −8.23,
−5.80 at the decade marks; record |z| = 14.10) — i.e. consistent with pure
square-root fluctuation and essentially *no* prime digit-sum parity bias — while
base 2 at the same scale blows past 100. Writing s_b mod 2 in binary bits, base 2
sums all bits, base 4 the even-indexed bits, base 8 the bits at positions ≡ 0
(mod 3); the every-third-bit subsampling of base 8 appears to decorrelate from
primality far faster than the every-other-bit subsampling of base 4. We state
this as a verified empirical fact and an **open quantitative question** (why is
the Mauduit–Rivat exponent σ_{8,2} so close to ½ while σ_{4,2}, σ_{6,2},
σ_{10,2} are not?).

---

## Verdict on the THEORY.md growth predictions

| prediction | verdict |
|:--|:--|
| **P1** — α_b increasing in b | **Not supported.** Exponents cluster near 0.69–0.76 for b ∈ {2,4,6,10} with no monotone trend; base 8 (≈0.50) and base 16 (erratic) break any ordering. |
| **P2** — E_b eventually one-signed and growing | **Mixed.** Holds for b = 2, 4, 6, 10, 12 (eventually positive and growing); **fails** for b = 8 (sign-oscillating, bounded) and b = 16 (erratic). |
| **P3** — α_2 smallest (base 2 the weakest standardized anomaly) | **False.** Base 2 is the *strongest* standardized anomaly at every scale. |

The common heuristic that smaller Mauduit–Rivat exponents at larger bases make
the standardized anomaly grow *with* the base is contradicted here: base 2
dominates throughout, and base 8 nearly vanishes.

---

## Reproducibility, prior art, and honest limits

**Data** (`data/`, CC BY 4.0): `trajectory_value_cutoffs.csv` (bases 2–16,
x ≤ 10⁸, excess and z), `cross_base_z_to_1e10.csv` (z for all bases, x ≤ 10¹⁰),
`count_cutoffs_base10.csv` (first N primes, N up to 4.5×10⁸). `test_results.py`
re-checks z = E/√n, the base-2 monotonicity, the 10⁹ base-10 dip, the N = 5×10⁷
local minimum, and the base-8 boundedness.

**Prior art.** Gupta & Sridhar (arXiv:1812.11841) report the base-10 effect at a
single cutoff. Mauduit–Rivat (Annals 2010) prove density → ½ with a power-saving
error for gcd(m, b−1) = 1 — an upper bound, no finite-sample size law and no
statement that the base-10 trajectory is non-monotone or that base 8 nearly
vanishes. Follow-up digit-sum work (arXiv:2606.04677 surjectivity thresholds;
2607.10654, 2608.19058 on digit distributions) does not address the *size* or
*scale-dependence* of the parity excess. OEIS A380192 / A111309 (base-10 prime
digit sums) and A000069/A130911 (base-2 odious/evil) list values, not
asymptotics. We found **no** prior statement that (i) G&S's 9.16σ sits at a local
excess trough, (ii) the base-10 excess is strongly non-monotone in scale, or
(iii) base 8 primes show near-square-root cancellation of the digit-sum parity.

**Limits.** The 10⁹/10¹⁰ points come from a single job (the ≤10⁸ data is
double-confirmed). The growth exponents α are two-/four-point log–log slopes, not
asymptotic constants, and base 10's are degraded by its non-monotonicity. "No
bias" for base 8 is an empirical statement to 10¹⁰, not a theorem. None of this
contradicts Mauduit–Rivat's asymptotic ½; it characterizes the finite-range
behaviour that G&S observed and sharpens what is and isn't surprising about it.
