# An analytic permutation null for the lag-k gap-parity correlation

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.*

## What this adds

`NULL-MOMENTS.md` (v0.6.3) gave the exact permutation-null mean and variance of the
**lag-1** sign-correlation statistic, and closed with: *"The same machinery extends
to the lag-k statistic (a union of paths), deferred to the lag-k write-up."* This is
that write-up. It supplies the exact first two permutation-null moments of the
**lag-k** statistic, for every lag `k >= 1`,

    S^(k) = sum_{n=1}^{N-k} chi_b(g_n) * chi_b(g_{n+k}),   chi_b(m) = (-1)^{s_b(m)},

against the same **marginal-preserving null**: the multiset of the `N` parity values
`{chi_b(g_n)}` is held fixed and only its order is randomised. With these, the
forthcoming lag-k measurement can report an analytic z,
`z_b^(k) = (S_b^(k) - E[S^(k)]) / sqrt(Var[S^(k)])`, at **every** lag with no Monte
Carlo error, keeping the empirical K-shuffle null only as an independent cross-check.
Until now only the lag-1 variance was available in closed form; the lag-k prereg
(`notes/gap-parity-lagk-prereg.md`) had only the lag-k null *mean*.

## The closed form

Let `(x_1, ..., x_N)` be a uniform random permutation of a fixed multiset of values
in `{+1, -1}`, write `s = sum_i x_i`, and fix a lag `k >= 1`. Define the lag-k
"edge set" `E_k = {(i, i+k) : 1 <= i <= N-k}`, so `S^(k) = sum_{(i,j) in E_k} x_i x_j`.
Then:

> **M := |E_k| = max(0, N - k)**  (number of lag-k pairs)
>
> **E[S^(k)] = M * rho_2,   rho_2 = (s^2 - N) / (N(N-1))**
>
> **E[(S^(k))^2] = M + 2 * n_shared * rho_2 + 2 * n_disjoint * rho_4**
>
> **Var[S^(k)] = E[(S^(k))^2] - E[S^(k)]^2**

where `rho_4 = (s^4 - 6 N s^2 + 8 s^2 + 3 N^2 - 6 N) / (N(N-1)(N-2)(N-3))` is the
mean product over four distinct positions, and the two structural counts are

    n_shared   = max(0, N - 2k)                  (unordered edge pairs sharing a vertex)
    n_disjoint = M(M-1)/2 - n_shared             (unordered edge pairs with no shared vertex).

`E[S^(k)] = (N-k)(s^2 - N)/(N(N-1))` is exactly the lag-k null mean already
pre-registered; the new content is the variance. With `s = -E_b^gap` (the gap-parity
excess of `GAP-PARITY.md`) these are functions of only `(N, s, k)`.

### Why the lag only enters through two counts

`E_k` is a disjoint union of `k` simple paths — one for each residue class of the
position index mod `k` — so it is a simple graph in which two distinct edges share at
most one vertex. Expanding `(S^(k))^2 = sum_{e,f in E_k} (edge_e)(edge_f)` splits by
overlap exactly as in the lag-1 case:

- **same edge** (`M` terms): `E[x_i^2 x_j^2] = 1`;
- **edges sharing one vertex** (`2 * n_shared` ordered terms):
  `E[x_i^2 x_j x_l] = E[x_j x_l] = rho_2`;
- **disjoint edges** (`2 * n_disjoint` ordered terms):
  `E[x_i x_j x_k x_l] = rho_4`.

A vertex has degree 2 in `E_k` exactly when `k+1 <= v <= N-k`, i.e. for `max(0, N-2k)`
vertices, each contributing one shared edge pair — hence `n_shared = max(0, N - 2k)`.
Collecting the three terms gives the `E[(S^(k))^2]` above. `rho_2` and `rho_4` come
from Newton's identities with power sums `p1=s, p2=N, p3=s, p4=N` (every `x_i^2 = 1`),
the same symmetric functions used for lag-1.

### Consistency with v0.6.3 (lag-1)

Setting `k = 1`: `M = N-1`, `n_shared = N-2`, `n_disjoint = (N-2)(N-3)/2`, and the
numerators collapse to

    E[(S^(1))^2] = (N-1) + (s^4 - 4 N s^2 + 4 s^2 + N^2 - 2 N) / (N(N-1)),

which is **identical** to the lag-1 `perm_e_s2` shipped in v0.6.3. The test checks
this reduction symbolically on every multiset as a regression gate.

### Worked checks (also in the test)

| multiset | N | s | k | E[S^(k)] | Var[S^(k)] |
|---|---:|---:|---:|---:|---:|
| {+,+,−} | 3 | 1 | 1 | −2/3 | 8/9 |
| {+,+,−} | 3 | 1 | 2 | −1/3 | 8/9 |
| {+,−,−,−} | 4 | −2 | 1 | 0 | 1 |
| {+,−,−,−} | 4 | −2 | 2 | 0 | 1 |
| {+,+,+,+} | 4 | 4 | k | N−k | 0 |

## Why it is not a numerical artefact

`test_null_moments_lagk.py` brute-forces **every** permutation of eight small `+/-1`
multisets and, for **each lag `k` from 1 to N−1**, checks the closed-form mean,
`E[S^2]` and variance against the exact enumeration with `fractions.Fraction` — no
floating point, no sampling. It includes the `N=7` set over all 5040 permutations and
all six lags, two one-sign multisets whose variance must be exactly 0 at every lag,
the degenerate empty statistic (`k >= N`), and the symbolic reduction to the v0.6.3
lag-1 form. The identities are exact.

## Prior art and honest attribution

The permutation (randomisation) moments of a serial-correlation statistic are
classical: **Wald & Wolfowitz, "An exact test for randomness in the non-parametric
case based on serial correlation", Ann. Math. Statist. 14 (1943) 378–388**, and the
general permutation-test framework of **E. J. G. Pitman (1937)**. Those treat the
**circular** statistic (`N` edges with wrap-around) and general real data. The only
contribution here is the explicit **linear/path, `+/-1`-multiset, general-lag-k**
specialisation — in which the lag enters solely through the two edge-overlap counts
`n_shared = max(0, N-2k)` and `n_disjoint` — giving a clean analytic null for this
particular observable (the digit-sum parity of prime gaps), with an exhaustive exact
check. **No novelty is claimed for the permutation-moment method itself.**

## How the lag-k measurement will use it

When the sealed lag-k measurement (`gap-parity-lagk-v1`, pre-registered in
`notes/gap-parity-lagk-prereg.md`) returns, the lag-k write-up will report, per base,
cutoff and lag `k = 1..5`, both the empirical K-shuffle `z` and the analytic
`z_b^(k) = (S_b^(k) - E[S^(k)]) / sqrt(Var[S^(k)])` from `perm_mean`/`perm_var`, and
note their agreement as a cross-check — exactly as v0.7.0 did at lag 1.

## Reproduce

```
python -m unittest -v test_null_moments_lagk   # exact, stdlib only, no network
```

Code Apache-2.0. — MV2, Munim, Inc.
