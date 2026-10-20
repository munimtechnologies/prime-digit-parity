# Verifying the odd-squarefree density claim behind the Axis-C control

*This note was produced with AI by MV2 of Munim, Inc. It verifies a load-bearing
claim with independent sources and a self-contained derivation; it reports only what
the sources and the arithmetic support.*

## The claim being verified

The data-verification block of `odd_squarefree_lag_signature_job.py` (the job behind
`ODD-SQUAREFREE-LAG-SIGNATURE-RESULTS.md`, released v0.21.0) anchors its sieve with an
*exact* self-check and a density comment:

```
assert qsf == 607926, qsf        # squarefree count <= 10^6
assert pc  == 78498,  pc         # prime count     <= 10^6
# odd squarefree are 2/3 of squarefree in density;
# here just re-check the full squarefree count exactly.
```

Two facts carry weight here and were, until now, asserted without an external citation:

1. the exact anchor **607926** (squarefree integers ≤ 10⁶), and
2. the density statement **"odd squarefree are 2/3 of squarefree."**

Both are confirmed below, independently of this repository's own code.

## 1. The exact anchor: squarefree ≤ 10⁶ = 607,926

OEIS **A071172** ("Number of squarefree integers ≤ 10ⁿ") lists, for n = 0,1,2,…:

```
1, 7, 61, 608, 6083, 60794, 607926, 6079291, 60792694, 607927124, …
```

So A071172(6) = **607926**, matching the job's assertion exactly, and A071172(7) =
6,079,291. The prime-count anchor π(10⁶) = **78498** is the standard value (OEIS A006880).
Both self-checks are therefore correct against external references, not just internally.

As a consistency check on the sieve's calibration, 607926 / 10⁶ = 0.607926, against the
squarefree density 6/π² = 0.60792710…; the two agree to ~10⁻⁶, well inside the proven
error term O(√x)/x = O(x^(−1/2)) ≈ 10⁻³ at x = 10⁶.

## 2. The density statement: odd squarefree are 2/3 of squarefree

Let d_sf be the natural density of squarefree integers and d_odd that of **odd**
squarefree integers. Writing each as an Euler product over the local conditions:

- Squarefree: for every prime p, exclude p² | n →  d_sf = ∏_p (1 − 1/p²) = 1/ζ(2) = **6/π²**.
- Odd squarefree: at p = 2 require n odd (local factor 1 − 1/2 = 1/2, not 1 − 1/4);
  at every odd prime p require p² ∤ n (local factor 1 − 1/p²). Hence

  d_odd = (1/2) · ∏_{p odd} (1 − 1/p²).

The odd-prime product follows by removing the p = 2 factor from 6/π²:

  ∏_{p odd} (1 − 1/p²) = (6/π²) / (1 − 1/4) = (6/π²) / (3/4) = **8/π²**,

so

  d_odd = (1/2) · (8/π²) = **4/π²** ≈ 0.405285.

Therefore

  d_odd / d_sf = (4/π²) / (6/π²) = **2/3**, exactly.

This is independently corroborated by OEIS **A056911** ("Odd squarefree numbers"), whose
asymptotic a(n) ~ n·π²/4 is the inverse-density form of d_odd = 4/π², and by D. Forgues'
comment there that among squarefree numbers the fraction p/(p+1) = 2/3 (at p = 2) are odd.
So the job's "2/3 in density" comment is exact, not approximate.

## 3. Why this matters for the result: data volume is ample

The odd-squarefree arm of Axis C removes even integers, which could, in principle, thin
the data enough to weaken the lag-signature test. The density pins down exactly how much
data survives:

- below 2×10⁶:  4/π² · 2×10⁶ ≈ **810,569** odd squarefree integers;
- below 2×10⁷ (the job's default `--n-max`):  4/π² · 2×10⁷ ≈ **8,105,695**.

Each leaves millions of consecutive-gap observations — far above any floor needed for the
standardized lag statistics in `ODD-SQUAREFREE-LAG-SIGNATURE-RESULTS.md`. Removing the
even integers costs exactly one third of the squarefree sample (6/π² → 4/π²), a known and
modest reduction, not an order-of-magnitude loss. The reported "loud + uniform" (signature
absent) outcome for odd squarefree is therefore not an artifact of data starvation.

## Limits

This note verifies the two *stated* anchors and the density ratio; it does not re-derive
the lag-signature result itself, which rests on the job's own (queue-bound) compute run.
The densities are asymptotic; at finite x the counts carry the standard O(√x) error, which
is negligible at the scales used here but is not zero.

## Sources

- OEIS A071172 — Number of squarefree integers ≤ 10ⁿ (A071172(6) = 607926). https://oeis.org/A071172
- OEIS A056911 — Odd squarefree numbers; a(n) ~ n·π²/4, density 4/π². https://oeis.org/A056911
- OEIS A006880 — π(10ⁿ); π(10⁶) = 78498. https://oeis.org/A006880
- Standard: density of squarefree integers = 1/ζ(2) = 6/π² (e.g. Hardy & Wright, *An
  Introduction to the Theory of Numbers*, §18.6).

— MV2, Munim, Inc., mv2@munimtech.com
