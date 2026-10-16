# An analytic permutation null for the lag-1 gap-parity correlation

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.*

## What this adds

The forthcoming `GAP-CORR.md` result measures the lag-1 serial correlation of the
digit-sum parity of prime gaps,

    S_b = sum_{n=1}^{N-1} chi_b(g_n) * chi_b(g_{n+1}),   chi_b(m) = (-1)^{s_b(m)},

against a **marginal-preserving null**: the exact multiset of the N parity values
{chi_b(g_n)} is held fixed and reshuffled, so only the *order* is randomised. Until
now that null was estimated empirically, by averaging S_b over K random shuffles and
taking the sample mean and sample standard deviation to form a z-score.

This note replaces the empirical null with the **exact** first two moments of S_b
under the uniform permutation of a +/-1 multiset. The connected signal can then be
reported as an analytic z, `z_b = (S_b - E[S_b]) / sqrt(Var[S_b])`, with no Monte
Carlo error, and the empirical shuffle mean/std become an independent cross-check.

## The closed form

Let `(x_1, ..., x_N)` be a uniform random permutation of a fixed multiset of values
in `{+1, -1}`, and write `s = sum_i x_i`. Because every `x_i^2 = 1`, we have
`sum_i x_i^2 = N`, `sum_i x_i^3 = s`, and `sum_i x_i^4 = N`. For the linear
(path, non-wrapping) statistic `S = sum_{i=1}^{N-1} x_i x_{i+1}`:

> **E[S] = (s^2 - N) / N**
>
> **E[S^2] = (N - 1) + (s^4 - 4 N s^2 + 4 s^2 + N^2 - 2 N) / (N (N - 1))**
>
> **Var[S] = E[S^2] - E[S]^2**

With `s = -E_b^gap` (the gap-parity excess of `GAP-PARITY.md`), `E[S]` is exactly the
null mean `(E_b^2 - N)/N` already used in `GAP-CORR`; the new content is the variance.

### Derivation

For distinct ordered position pairs the symmetry of a uniform permutation gives
`rho_2 := E[x_i x_j] = (sum_{k != l} x_k x_l)/(N(N-1)) = (s^2 - N)/(N(N-1))`, so
`E[S] = (N-1) rho_2 = (s^2 - N)/N`.

Expanding `S^2 = sum_{e,f} (edge_e)(edge_f)` over the `N-1` path edges splits by how
many positions the two edges share:

- **same edge** (`N-1` terms): `E[x_i^2 x_{i+1}^2] = 1`;
- **edges sharing one vertex** (`2(N-2)` ordered terms): `E[x_i x_{i+1}^2 x_{i+2}]
  = E[x_i x_{i+2}] = rho_2`;
- **disjoint edges** (`(N-2)(N-3)` ordered terms): `E[x_i x_{i+1} x_j x_{j+1}]
  = rho_4`, the mean product over four distinct positions.

Newton's identities with `p_1=s, p_2=N, p_3=s, p_4=N` give the fourth elementary
symmetric sum, hence `M_4 = sum_{distinct 4-tuples} = s^4 - 6 N s^2 + 8 s^2 + 3 N^2
- 6 N` and `rho_4 = M_4 / (N(N-1)(N-2)(N-3))`. Collecting the three terms,

    E[S^2] = (N-1) + [2(N-2)(s^2 - N) + M_4] / (N(N-1))
           = (N-1) + (s^4 - 4 N s^2 + 4 s^2 + N^2 - 2 N) / (N(N-1)),

and `Var[S] = E[S^2] - ((s^2-N)/N)^2`.

### Worked checks (also in the test)

| multiset | N | s | E[S] | Var[S] |
|---|---:|---:|---:|---:|
| {+,+,−} | 3 | 1 | −2/3 | 8/9 |
| {+,−,−,−} | 4 | −2 | 0 | 1 |
| {+,+,+,+} | 4 | 4 | 3 | 0 |

## Why it is not a numerical artefact

`test_null_moments.py` brute-forces **every** permutation of eight small `+/-1`
multisets (including the one-edge path `N=2`, two one-sign multisets whose variance
must be exactly 0, and an `N=7` set over all 5040 permutations) and checks the
closed-form mean, `E[S^2]` and variance against the exact enumeration with
`fractions.Fraction` — no floating point, no sampling. The identities are exact.

## Prior art and honest attribution

The permutation (randomisation) moments of a serial-correlation statistic are
classical: **Wald & Wolfowitz, "An exact test for randomness in the non-parametric
case based on serial correlation", Ann. Math. Statist. 14 (1943) 378–388**, and the
general permutation-test framework of **E. J. G. Pitman (1937)**. Those treat the
**circular** statistic (`N` edges with wrap-around) and general real data. The
contribution here is only the explicit **linear/path, `+/-1`-multiset**
specialisation that yields a clean analytic null for this particular new observable
(the digit-sum parity of prime gaps), together with an exhaustive exact check. No
novelty is claimed for the permutation-moment method itself.

## How GAP-CORR uses it

When the sealed measurement of `S_b` returns, `GAP-CORR.md` will report, per base and
cutoff, both the empirical K-shuffle `z` and the analytic `z_b = (S_b - E[S_b]) /
sqrt(Var[S_b])` from `null_moments.perm_mean`/`perm_var`, and will note their
agreement as a cross-check. The same machinery extends to the lag-k statistic (a
union of paths), deferred to the lag-k write-up.

## Reproduce

```
python -m unittest -v test_null_moments   # exact, stdlib only, no network
```

Code Apache-2.0. — MV2, Munim, Inc.
