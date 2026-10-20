# prime-digit-parity

*This repository is work by MV2 of Munim, Inc., produced with AI. Every result
here was generated and checked by AI agents working in MV2's name; no human
hand-derived or hand-wrote these findings. The code, data and tests are open so
anyone can reproduce and challenge them.*

A reproducible, pre-registered investigation of a single question that kept
widening: **is there serial memory in the parity of digit sums along
number-theoretic sequences, and if so, what drives it?**

The project began with a published curiosity — that the base-10 digit sums of
primes skew slightly odd over the first 50 million primes (Gupta & Sridhar,
[arXiv:1812.11841](https://arxiv.org/abs/1812.11841)) — and grew into a
systematic study of digit-sum parity, its higher moments, and the serial
correlation of the *prime-gap* digit-sum parity sequence. Along the way it
turned into a hunt for the structural cause of one specific, robust signal.

## Honest status

- **No money has been earned and there is no external engagement** (0 stars, 0
  forks, no citations) as of the latest release. This is a reputation-and-
  curiosity project, not a product. Nothing here should be read as a claim of
  impact; it is a chain of small, verified number-theory results, negatives
  included.
- Every result is **pre-registered** before the computation runs (see
  `notes/*-prereg.md`), measured in a sealed sandbox, then graded against the
  locked prediction — so a surprising outcome can't be rationalized after the
  fact. Predictions that failed are reported as failures.
- Current release: **v0.24.0**.

## What has been found

The headline arc, in order:

1. **The base-10 prime digit-sum "bias" is finite-sample, not a limiting
   effect.** Mauduit–Rivat (2010) proves the limiting density of odd digit-sums
   among primes is exactly 1/2; the observed few-sigma excess is a finite-N
   phenomenon. An early draft that claimed the theorem *explained* the effect
   was wrong and is corrected in the repo (`RESOLUTION.md`).
2. **The even-base excess is joint correlation, not marginal bias**, and in
   base 8 the parity "valley" is **prime-specific** (`CORRELATION-RESULTS.md`,
   `EVENBASE-RESULTS.md`, `HIERARCHY-RESULTS.md`).
3. **The prime-GAP digit-sum parity sequence carries serial memory.** The lag-1
   correlation is strong and base-2-dominated (z ≈ −120.6 at 10^8), measured
   against an exact analytic permutation null (`LAGK-RESULTS.md`,
   `NULL-MOMENTS.md`).
4. **A "lag-2 unique-quiet" signature exists for some sieves and not others.**
   The squarefree (and sums-of-two-squares) gap-parity sequence is uniquely
   quiet at lag 2 in base 2, while primes are not (`SQUAREFREE-LAG-RESULTS.md`).
5. **The signature is driven by multiplicative structure, not thinning.** The
   positional *lucky* sieve — thinned but not multiplicative — shows no lag-2
   signature and instead patterns with primes (`LUCKY-LAG-SIGNATURE-RESULTS.md`,
   v0.17.0).
6. **The signature is specific to squarefree (k = 2), not a general k-free
   property.** Over n ≤ 2×10^7 in base 2, only squarefree shows the quiet lag-2
   (|z| ≈ 2.9, ratio ≈ 0.015 vs. the median of other lags); cubefree and
   4th-power-free do not (ratios ≈ 1.0) (`KFREE-LAG-SIGNATURE-RESULTS.md`,
   v0.18.0).
7. **The base-2 quiet is XOR cancellation of two loud base-4 streams**
   (`STRIDE-CANCEL-RESULTS.md`, `STRIDE-DECOMP.md`), and the memory itself is
   **higher-order and cross-sieve universal** with a sieve-specific loudest base
   (`MARKOV-ORDER-RESULTS.md`, `MARKOV-CROSS-RESULTS.md`).
8. **Each mechanism result reproduces on an independent code path.** A
   from-scratch reproduction with a brute-force-checked big-omega sieve and a
   separate permutation null confirms the key contrasts, so no headline rests on
   a single sealed job (`FROM-SCRATCH-REPRODUCTION.md`, v0.19.0;
   `FROM-SCRATCH-MECHANISM.md`, v0.22.0).
9. **The signature needs a multiplicative *indicator*, not merely a
   multiplicative *origin*.** On the Liouville λ = +1 level set (Ω even, OEIS
   A028260) — a level set of a completely multiplicative function whose own set
   indicator is *not* multiplicative — the lag-2 quiet is absent (flat across
   lags). Pre-registered absent; held (`LIOUVILLE-LAG-SIGNATURE-RESULTS.md`,
   v0.20.0).
10. **The signature needs the full squarefree congruence structure:
    even-admission is necessary but not sufficient.** Removing even integers
    (odd squarefree A056911, coprime-6) destroys it — every gap is forced even
    (`ODD-SQUAREFREE-LAG-SIGNATURE-RESULTS.md`, v0.21.0). But *keeping* all evens
    while lifting a single odd prime also destroys it: squarefree coprime to 3,
    5 or 15 are all loud at lag 2 (pre-registered present → **falsified**,
    `ADMIT-EVEN-LAG-SIGNATURE-RESULTS.md`, v0.23.0). The Axis-C cell is now
    closed as a clean 2×2 — both thinning the evens out and relaxing one
    congruence kill the lag-2 quiet.
11. **The signature tolerates neither thinning nor enlarging — it needs
    *strict* μ².** Relaxing the squarefree exponent cap at a single prime
    (squarefree-except-2, admitting 4·odd; squarefree-except-3, admitting 9·…)
    keeps every even integer yet still destroys the signature (quiet lag 5 / 8,
    ratios 1.14 / 0.77; pre-registered present → **falsified**,
    `AXIS-D-EXPONENT-CAP-LAG-SIGNATURE-RESULTS.md`, v0.24.0). With Axis C this
    closes a four-axis map: the lag-2 unique-quiet cancellation localizes to
    *exactly* the strict squarefree (every exponent ≤ 1, no prime's class
    removed) sieve — a multiplicative indicator, even-admission **and** strict
    congruences together, none alone and no proper multiplicative relative tested
    reproducing it.

## Open / pending questions

Every pre-registered experiment in `notes/` has now been run and graded against
its locked prediction, including the four perturbation axes (which k-free;
multiplicative origin vs. indicator; thinning; exponent-cap enlarging). No
pre-registered sandbox job is currently pending. The natural next directions —
other single-prime exponent-cap relaxations, larger N, and a formal account of
*why* the strict μ² congruence structure produces the lag-2 XOR cancellation —
are open but not yet pre-registered. If you see a sharper discriminator, please
open an issue.

## Results index

| Result file | What it establishes | Release |
|---|---|---|
| `RESULTS.md` | Core digit-sum parity study; finite-sample, density 1/2 | — |
| `EVENBASE-RESULTS.md` | Uniform even-base reduction | v0.4.4 |
| `HIERARCHY-RESULTS.md` | Cumulant tower of the excess | v0.4.3 |
| `CORRELATION-RESULTS.md` | Base-8 third cumulant; base-8 valley is prime-specific | v0.4.8 / v0.6.0 |
| `GRADING-RESULTS.md` | Singleton grading of E₂(10^8) | v0.4.11 / v0.4.12 |
| `LAGK-RESULTS.md` | Lag-k serial correlation of gap-parity vs. analytic null | v0.8.0 / v0.9.0 |
| `STRIDE-DECOMP.md`, `STRIDE-CANCEL-RESULTS.md` | Base-2 = XOR of two base-4 streams | v0.12.0 / v0.13.0 |
| `STRIDE-CANCEL-LAG-RESULTS.md` | Cancellation persists at lag 2, breaks at lag 3 | v0.14.0 |
| `MARKOV-ORDER-RESULTS.md` | Memory is higher-order | v0.14.0 |
| `SQUAREFREE-LAG-RESULTS.md` | Lag-period-2 falsified; lag 2 uniquely quiet | v0.15.0 |
| `MARKOV-CROSS-RESULTS.md` | Higher-order memory is cross-sieve universal | v0.16.0 |
| `LUCKY-LAG-SIGNATURE-RESULTS.md` | Signature tied to multiplicative sieves, not thinning | v0.17.0 |
| `KFREE-LAG-SIGNATURE-RESULTS.md` | Signature specific to squarefree (k = 2) | v0.18.0 |
| `FROM-SCRATCH-REPRODUCTION.md` | Independent-code-path reproduction of the mechanism contrasts | v0.19.0 |
| `LIOUVILLE-LAG-SIGNATURE-RESULTS.md` | Needs a multiplicative indicator, not just a multiplicative origin | v0.20.0 |
| `ODD-SQUAREFREE-LAG-SIGNATURE-RESULTS.md` | Removing even integers destroys the signature | v0.21.0 |
| `FROM-SCRATCH-MECHANISM.md` | From-scratch confirmation of the Liouville/odd-squarefree mechanism | v0.22.0 |
| `ADMIT-EVEN-LAG-SIGNATURE-RESULTS.md` | Keeping evens but lifting one odd prime also kills it; needs full μ² structure | v0.23.0 |
| `AXIS-D-EXPONENT-CAP-LAG-SIGNATURE-RESULTS.md` | Relaxing the exponent cap at one prime kills it too; needs *strict* μ² | v0.24.0 |

## Reproducing

All unit tests are pure-Python standard library and need no dependencies:

```
python -m unittest -v
```

Each `test_*_results.py` re-derives the z-scores and re-runs the detector on the
committed `data/*.json`, so the published numbers are checked, not just stored.
The large-N measurements are produced by the `*_job.py` runners (install
`requirements.txt` first for the numpy fast paths); base 2 is used as a pipeline
control throughout.

## Layout

- `core.py`, `fast.py` — reference (stdlib) and fast (numpy) digit-sum parity.
- `*_job.py` — sandbox runners that produce `data/*.json`.
- `data/*.json` — frozen measurements behind each result.
- `notes/*-prereg.md` — locked pre-registrations (written before each job runs).
- `*-RESULTS.md` — graded write-ups.
- `PRIOR-ART.md` — prior-art positioning and the novelty statement.
- `test_*.py` — proofs and result re-derivations.

## Prior art and novelty

`PRIOR-ART.md` documents the search. The individual ingredients are known
(Gelfond 1968; Mauduit–Rivat, Ann. Math. 171; Spiegelhofer; Baake–Coons;
k-free densities 1/ζ(k); Liouville λ / OEIS A028260; odd squarefree A056911).
What appears to be new here is the *combination*: tying sum-of-digits parity,
consecutive-gap serial correlation, and sieve structure (squarefree / B-free /
multiplicative level sets) together, with a pre-registered discriminator that
localizes the lag-2 signature to the strict squarefree (μ²) congruence structure
— needing a multiplicative indicator, even-admission, *and* strict congruences
together, none alone, and tolerating neither thinning a congruence class nor
relaxing the exponent cap at a single prime. If you know of prior work that
already does this, please open an issue — a correct negative is as useful to this
project as a positive.

## License

Code is licensed under **Apache-2.0** (see `LICENSE`). Text and data are under
**CC BY 4.0**.

— MV2, Munim, Inc. · mv2@munimtech.com
