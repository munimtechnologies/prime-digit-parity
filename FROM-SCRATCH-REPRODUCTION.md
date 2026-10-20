# From-scratch reproduction of the squarefree lag-2 gap-parity signature

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. Code is Apache-2.0;
this text is CC BY 4.0.*

## What this adds

Every other statistical test in this repository checks a result by re-deriving its
numbers from a committed JSON file that a sealed compute job produced at N = 2×10⁷.
That is good for auditing a specific run, but it trusts the committed data.

`test_from_scratch_signature.py` trusts **no committed data**. With the Python
standard library only, it:

1. sieves the squarefree integers in [1, N] and the primes in [2, N] from scratch
   (slice-assignment sieves),
2. forms the consecutive gaps of each sequence,
3. maps every gap `g` to its base-2 digit-sum parity `t(g)` (the Thue-Morse bit),
4. recomputes, for lags 1…12, the agreement z-score of the parity stream under an
   **exact marginal-preserving permutation null**, and
5. runs the published "lag-2 unique-quiet" detector.

A reader can run `python -m unittest test_from_scratch_signature` and watch the
signature appear in a few seconds, with no downloads and no reliance on any earlier
result in this repo.

## The permutation null is proven correct, not assumed

For a fixed bit multiset (n₁ ones among L bits), the lag-k agreement count T_k has
an exact mean and variance under a uniformly random permutation. The lag-k position
pairs split into k chains by residue mod k; consecutive edges of a chain share a
vertex, every other edge pair is vertex-disjoint, and the covariance depends only on
which case holds. This gives a closed form in falling factorials of n₁ and n₀.

`TestNullMomentFormula` checks that closed form against **exhaustive brute force**
over all permutations for six small cases, so the z-scores used in the reproduction
rest on a formula that is verified, not merely stated.

## Result

At N = 2,000,000 (about 1.2 million squarefree gaps and about 150 thousand prime
gaps):

- **Squarefree:** lag 2 is the uniquely quietest lag, and its |z| is at most 10% of
  the median |z| of the other lags. The detector fires. This is the same
  "lag-2 unique-quiet" signature reported at N = 2×10⁷ in `KFREE-LAG-SIGNATURE-RESULTS.md`.
- **Primes (negative control):** the detector does not fire.

The headline project result was established at N = 2×10⁷; this artifact confirms the
same qualitative signature, and the same squarefree-vs-primes contrast, holds at a
far smaller N that an ordinary test run can recompute end to end. It is a
reproducibility artifact, not a new claim: the novelty here is that the signature is
now independently reproducible from first principles in seconds, with the null
distribution itself checked against brute force.

## Honest limits

- The absolute z magnitudes scale with √N, so they are smaller than at N = 2×10⁷;
  the *detector verdict* (which is scale-invariant in form) is what is reproduced.
- This confirms, it does not extend, the published finding. The open mechanism
  questions (even-admission vs the squarefree exponent cap) remain the subject of the
  pre-registered compute jobs described in the `notes/` prereg files.
