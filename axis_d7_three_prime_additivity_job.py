#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Produced with AI by MV2 of Munim, Inc.
"""Axis D7: does the density scaling law stay ADDITIVE for THREE primes? Relax the strict squarefree
(mu^2) exponent cap at THREE primes at once and test whether the super-control excess of the ratio is
the SUM of the three single-prime excesses, or whether a higher-order (pairwise/triple) interaction
term appears.

Axis D5 (v0.29.0) established the single-prime density scaling law: squarefree-except-p perturbs the
frozen ratio r(p) = |z(2)| / median|z(k!=2)| above the squarefree control by an amount proportional
to the admitted-class density d(p) = 1/p^2 - 1/p^3:
    r(p) - r_control  ~  C * d(p),   C ~ 8   (k(p) = (r(p)-r_control)/d(p) in [5,11] for p = 2..37).
Axis D6 (v0.30.0) showed that law is ADDITIVE across TWO primes:
    r({p,q}) - r_control  ~  C * (d(p) + d(q)),   k({p,q}) in [5,11] for {17,19},{23,29},{31,37}.
If three independent single-prime relaxations remain independent, then relaxing at p, q AND s at once
(squarefree-except-{p,q,s}, admitting the p^2-, q^2- and s^2-classes) should perturb the ratio by the
SUM of the three single-prime densities:
    r({p,q,s}) - r_control  ~  C * (d(p) + d(q) + d(s)).
(The pairwise-overlap classes {p^2 q^2 m} have density ~1/(p^2 q^2) ~ 1e-5-1e-6 and the triple overlap
~1e-9, all negligible vs d(p) ~ 1e-3, so the admitted extra density of the union is d(p)+d(q)+d(s) to
first order.) This job tests that THREE-prime additivity prediction out-of-sample at three disjoint
prime triples chosen so the summed density keeps the ratio comfortably below the 0.1 detector
threshold (so all three should stay PRESENT):
    {17,19,23}: d=0.007689   {29,31,37}: d=0.002866   {41,43,47}: d=0.001551.
The last triple {41,43,47} was never measured singly, so it is a fully out-of-sample test of both
convergence and additivity.

Six sets to a common cutoff N:
  * squarefree             (control)      : c_p = 1 all p                      (A005117, 6/pi^2)
  * squarefree_except_17_19_23            : c_17=c_19=c_23=2, c_p=1 otherwise  (allow 289|n,361|n,529|n)
  * squarefree_except_29_31_37            : c_29=c_31=c_37=2, c_p=1 otherwise  (allow 841|n,961|n,1369|n)
  * squarefree_except_41_43_47            : c_41=c_43=c_47=2, c_p=1 otherwise  (allow 1681|n,1849|n,2209|n)
  * odd_squarefree         (anchor)       : c_2 = 0, c_p = 1 odd p             (A056911, all gaps even)
  * primes                 (neg control)  : A000040

Pre-registered in prime-digit-parity/notes/axis-d7-three-prime-additivity-strict-mu2-prereg.md
(committed before this job runs). Primary predictions: P-D7-sf PRESENT, P-D7-triple-{17_19_23,
29_31_37,41_43_47} PRESENT, P-D7-odd ABSENT, P-D7-primes ABSENT. Secondary (the new test -- THREE-prime
ADDITIVITY of the density law): additivity_law_holds, i.e. the three-prime coefficient
    k({p,q,s}) = (r({p,q,s}) - r_control) / (d(p) + d(q) + d(s))
lands in the committed band [5, 11] for each triple. Named alternatives (reported honestly): a triple
scoring ABSENT (r > 0.1) even though each single and double relaxation at these primes is present
would mean the three relaxations interact destructively (super-additive destruction); all present but
some k outside [5,11] would mean convergence survives but the excess is not simply the sum of the three
densities (a pairwise or triple cross term appears).

The observable, exact permutation null `null_z`, base-4 stride split, stream analysis, scale-free
detector, general exponent-cap sieve and prime-gap collector are reused VERBATIM (and re-self-tested)
from axis_d6_density_additivity_job.py. The exponent-cap sieve is self-tested against a full
factorization over [1, 20000] (count, the complete gap stream and the max gap, a tiny segment forcing
boundary crossings) for caps {}, {2:0}, {17:2,19:2,23:2}, {29:2,31:2,37:2}, {41:2,43:2,47:2} and the
union of all nine relaxed primes, plus the mechanism contrast (each triple admits even members and
produces odd gaps) and the exact counts Q(10^6)=607926, pi(10^6)=78498. Self-tests run first; the job
exits non-zero on any failure.
Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.
"""
import argparse
import json
import math
import os
import statistics
import sys
from array import array
from fractions import Fraction
from itertools import permutations

LAGS = [1, 2, 3, 4, 5, 6, 7, 8]

# The three disjoint prime triples relaxed together (each prime in a triple gets cap 2).
TRIPLES = {
    "17_19_23": (17, 19, 23),
    "29_31_37": (29, 31, 37),
    "41_43_47": (41, 43, 47),
}
THRESHOLD = 0.1
# Pre-registered band for the additivity coefficient k({p,q,s}) = (r - r_control)/(d(p)+d(q)+d(s)),
# the SAME band validated for the single-prime law (Axis D5) and the two-prime law (Axis D6).
ADDITIVITY_LAW_BAND = (5.0, 11.0)


# ------------------------------------------------- observable + null (verbatim from arc)
def digit_sum_parity(n, b):
    s = 0
    while n:
        s += n % b
        n //= b
    return s & 1


def lag_stats(pb):
    n_eff = len(pb)
    ones = pb.count(1)
    s = n_eff - 2 * ones
    out = {}
    for k in LAGS:
        if k >= n_eff:
            out[k] = None
            continue
        a = int.from_bytes(pb[: n_eff - k], "big")
        c = int.from_bytes(pb[k:], "big")
        d_k = (a ^ c).bit_count()
        m = n_eff - k
        out[k] = m - 2 * d_k
    return n_eff, s, out


def null_z(n_eff, s, S, k):
    N_ = n_eff
    M = N_ - k
    n_shared = max(0, N_ - 2 * k)
    n_disjoint = M * (M - 1) // 2 - n_shared
    rho2 = Fraction(s * s - N_, N_ * (N_ - 1))
    rho4 = Fraction(
        s ** 4 - 6 * N_ * s * s + 8 * s * s + 3 * N_ * N_ - 6 * N_,
        N_ * (N_ - 1) * (N_ - 2) * (N_ - 3),
    )
    eS = M * rho2
    eS2 = M + 2 * n_shared * rho2 + 2 * n_disjoint * rho4
    var = eS2 - eS * eS
    varf = float(var)
    z = (S - float(eS)) / math.sqrt(varf) if varf > 0 else 0.0
    return float(eS), varf, z


# ------------------------------------------------- sieves
def _small_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    return [i for i in range(2, limit + 1) if sieve[i]]


def collect_prime_gaps(N):  # verbatim
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    gaps = array("H")
    prev = None
    maxg = 0
    cnt = 0
    for n in range(2, N + 1):
        if sieve[n]:
            cnt += 1
            if prev is not None:
                g = n - prev
                if g > maxg:
                    maxg = g
                gaps.append(g)
            prev = n
    del sieve
    return gaps, cnt, maxg


def collect_exponent_cap_gaps(N, caps=None, seg_size=1 << 22):
    """Segmented sieve of {n <= N : v_p(n) <= caps.get(p, 1) for every prime p}, i.e. n not divisible
    by p^(cap_p + 1) for any p. caps={} -> squarefree (A005117); {17:2,19:2,23:2} ->
    squarefree-except-{17,19,23} (allow 289|n,361|n,529|n, forbid 4913|n,6859|n,12167|n); {2:0} -> odd
    squarefree (A056911). We mark multiples of each forbidden prime power q = p^(cap_p + 1) <= N; a
    surviving position is a member."""
    caps = dict(caps or {})
    root = int(round(N ** 0.5)) + 2
    primes = _small_primes(root)
    qs = []
    for p in primes:
        cap = caps.get(p, 1)
        q = p ** (cap + 1)
        if q <= N:
            qs.append(q)
    gaps = array("H")
    cnt = 0
    prev = None
    maxg = 0
    lo = 1
    while lo <= N:
        hi = min(lo + seg_size, N + 1)
        length = hi - lo
        seg = bytearray(b"\x01") * length
        for q in qs:
            if q > hi - 1:
                continue
            start = ((lo + q - 1) // q) * q
            for m in range(start, hi, q):
                seg[m - lo] = 0
        for off in range(length):
            if seg[off]:
                v = lo + off
                cnt += 1
                if prev is not None:
                    g = v - prev
                    if g > maxg:
                        maxg = g
                    gaps.append(g)
                prev = v
        lo = hi
    return gaps, cnt, maxg


# ------------------------------------------------- per-sieve stream analysis (verbatim)
def stream_bit(v, which):
    if which == "a_base4_even":
        return digit_sum_parity(v, 4)
    if which == "b_base4_odd":
        return digit_sum_parity(v >> 1, 4)
    if which == "x_base2":
        return digit_sum_parity(v, 2)
    raise ValueError(which)


STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]


def analyze_sieve(gaps, max_gap):
    out = {}
    for name in STREAMS:
        tbl = bytes(stream_bit(v, name) for v in range(max_gap + 1))
        pb = bytes(tbl[g] for g in gaps)
        n_eff, s, S = lag_stats(pb)
        frac_odd = (n_eff - s) / 2 / n_eff if n_eff else 0.0
        per_lag = {}
        for k in LAGS:
            if S[k] is None:
                continue
            eS, var, z = null_z(n_eff, s, S[k], k)
            per_lag[str(k)] = {"S": S[k], "E_null": eS, "var_null": var, "z": z}
        out[name] = {"n_eff": n_eff, "s": s, "marginal_frac_odd": frac_odd,
                     "lags": per_lag}
    return out


# ------------------------------------------------- scale-free signature detector (verbatim)
def lag2_signature(x_lags):
    az = {k: abs(x_lags[str(k)]["z"]) for k in LAGS}
    quiet_lag = min(az, key=az.get)
    others = [az[k] for k in LAGS if k != 2]
    med_others = statistics.median(others)
    ratio = (az[2] / med_others) if med_others > 0 else float("inf")
    present = bool(quiet_lag == 2 and az[2] <= THRESHOLD * med_others)
    return {
        "abs_z_by_lag": {str(k): az[k] for k in LAGS},
        "quiet_lag": quiet_lag,
        "abs_z_lag2": az[2],
        "median_abs_z_other_lags": med_others,
        "ratio_lag2_over_median_others": ratio,
        "signature_present": present,
    }


# ------------------------------------------------- density scaling law
def density(p):
    """Density of the single-prime admitted extra class {p^2 m : p does not divide m}: 1/p^2 - 1/p^3."""
    return 1.0 / (p * p) - 1.0 / (p * p * p)


# ------------------------------------------------- self tests
def digit_parity_tests():
    assert digit_sum_parity(19, 10) == 0
    assert digit_sum_parity(7, 2) == 1
    assert digit_sum_parity(8, 8) == 1
    assert digit_sum_parity(255, 16) == 0


def test_stride_identity():
    for g in range(0, 4096):
        a = digit_sum_parity(g, 4)
        b = digit_sum_parity(g >> 1, 4)
        assert (a ^ b) == digit_sum_parity(g, 2), g
        assert digit_sum_parity(g, 2) == bin(g).count("1") % 2, g


def test_lag_stats():
    pb = bytes([0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1])
    n_eff, s, S = lag_stats(pb)
    assert s == (pb.count(0) - pb.count(1))
    for k in LAGS:
        if k >= n_eff:
            continue
        disag = sum(1 for i in range(len(pb) - k) if pb[i] != pb[i + k])
        assert S[k] == (len(pb) - k) - 2 * disag, (k, S[k])


def test_null_formula():
    for seq in ([1, 1, 1, -1, -1, -1, 1, -1], [1, 1, 1, 1, -1, -1, 1]):
        N_ = len(seq)
        s = sum(seq)
        for k in (1, 2, 3, 4):
            tot = tot2 = cnt = 0
            for perm in permutations(seq):
                Sk = sum(perm[i] * perm[i + k] for i in range(N_ - k))
                tot += Sk
                tot2 += Sk * Sk
                cnt += 1
            emp_mean = Fraction(tot, cnt)
            emp_m2 = Fraction(tot2, cnt)
            M = N_ - k
            n_shared = max(0, N_ - 2 * k)
            n_disjoint = M * (M - 1) // 2 - n_shared
            rho2 = Fraction(s * s - N_, N_ * (N_ - 1))
            rho4 = Fraction(
                s ** 4 - 6 * N_ * s * s + 8 * s * s + 3 * N_ * N_ - 6 * N_,
                N_ * (N_ - 1) * (N_ - 2) * (N_ - 3),
            )
            assert emp_mean == M * rho2, (k, emp_mean, M * rho2)
            assert emp_m2 == M + 2 * n_shared * rho2 + 2 * n_disjoint * rho4, k


def _in_set(n, caps):
    """Obviously-correct membership by full trial-division factorization: n is a member iff
    v_p(n) <= caps.get(p, 1) for every prime p."""
    if n < 1:
        return False
    m = n
    d = 2
    while d * d <= m:
        if m % d == 0:
            e = 0
            while m % d == 0:
                m //= d
                e += 1
            if e > caps.get(d, 1):
                return False
        d += 1
    if m > 1 and 1 > caps.get(m, 1):
        return False
    return True


def _ref_members(limit, caps):
    return [v for v in range(1, limit + 1) if _in_set(v, caps)]


def test_exponent_cap_sieve():
    # Membership checker sanity against hand facts for the three-prime relaxations.
    assert _in_set(1, {}) and _in_set(2, {}) and _in_set(30, {})
    assert not _in_set(4, {}) and not _in_set(289, {})          # 2^2, 17^2
    # {17,19,23}: squares allowed, cubes forbidden.
    T1 = {17: 2, 19: 2, 23: 2}
    assert _in_set(289, T1) and _in_set(361, T1) and _in_set(529, T1)       # 17^2, 19^2, 23^2 allowed
    assert not _in_set(4913, T1) and not _in_set(6859, T1) and not _in_set(12167, T1)  # 17^3,19^3,23^3
    assert not _in_set(4, T1)                                   # unrelated square still forbidden
    # {29,31,37}
    T2 = {29: 2, 31: 2, 37: 2}
    assert _in_set(841, T2) and _in_set(961, T2) and _in_set(1369, T2)      # 29^2, 31^2, 37^2 allowed
    assert not _in_set(24389, T2) and not _in_set(29791, T2) and not _in_set(50653, T2)  # 29^3,31^3,37^3
    # {41,43,47}
    T3 = {41: 2, 43: 2, 47: 2}
    assert _in_set(1681, T3) and _in_set(1849, T3) and _in_set(2209, T3)    # 41^2, 43^2, 47^2 allowed
    assert not _in_set(68921, T3) and not _in_set(79507, T3) and not _in_set(103823, T3)  # 41^3,43^3,47^3
    assert not _in_set(2, {2: 0}) and _in_set(3, {2: 0})        # odd squarefree excludes evens
    # Exact reference match (count, FULL gap stream, max gap) over [1, 20000] for each cap set, with a
    # deliberately tiny segment so many segment boundaries are crossed.
    LIMIT = 20000
    cap_sets = ({}, {2: 0}, T1, T2, T3,
                {17: 2, 19: 2, 23: 2, 29: 2, 31: 2, 37: 2, 41: 2, 43: 2, 47: 2})
    for caps in cap_sets:
        members = _ref_members(LIMIT, caps)
        ref_gaps = [members[i + 1] - members[i] for i in range(len(members) - 1)]
        gaps, cnt, maxg = collect_exponent_cap_gaps(LIMIT, caps, seg_size=97)
        assert cnt == len(members), (caps, cnt, len(members))
        assert list(gaps) == ref_gaps, caps
        assert maxg == (max(ref_gaps) if ref_gaps else 0), (caps, maxg)
    # MECHANISM CONTRAST: each three-prime set admits even integers AND produces at least one odd gap,
    # unlike the even-removing anchor.
    for caps in (T1, T2, T3):
        members = _ref_members(20000, caps)
        gaps, _, _ = collect_exponent_cap_gaps(20000, caps, seg_size=97)
        assert any(v % 2 == 0 for v in members), ("no even member", caps)
        assert any(g % 2 == 1 for g in gaps), ("no odd gap", caps)
    g_odd, _, _ = collect_exponent_cap_gaps(20000, {2: 0}, seg_size=97)
    assert all(x % 2 == 0 for x in g_odd), "odd squarefree must have all even gaps"
    # A005117 count Q(10^6) = 607926 (Gegenbauer); re-check the squarefree control exactly.
    _, qsf, _ = collect_exponent_cap_gaps(10 ** 6, {}, seg_size=1 << 20)
    assert qsf == 607926, qsf


def test_primes():
    _, pc, _ = collect_prime_gaps(10 ** 6)
    assert pc == 78498, pc


def test_signature_detector():
    def mk(zs):
        return {str(k): {"z": zs[k - 1]} for k in range(1, 9)}
    sq_like = mk([-502, -7, -496, -426, -441, -439, -427, -334])
    assert lag2_signature(sq_like)["signature_present"] is True
    pr_like = mk([112.8, 11.4, 0.62, 0.95, 0.19, 2.6, 2.2, 0.64])
    assert lag2_signature(pr_like)["signature_present"] is False
    weak = mk([10, 3, 12, 11, 9, 13, 8, 10])
    assert lag2_signature(weak)["signature_present"] is False


def test_density_helper():
    # d(p) = 1/p^2 - 1/p^3, a positive decreasing function of p; triple density adds the three.
    assert abs(density(2) - (0.25 - 0.125)) < 1e-15
    assert density(3) < density(2) and density(47) < density(41)
    assert abs((density(17) + density(19) + density(23)) - 0.0076891) < 1e-5  # summed triple sanity
    assert abs((density(29) + density(31) + density(37)) - 0.0028658) < 1e-5
    assert abs((density(41) + density(43) + density(47)) - 0.0015517) < 1e-5


def run_selftests():
    digit_parity_tests()
    test_stride_identity()
    test_lag_stats()
    test_null_formula()
    test_exponent_cap_sieve()
    test_primes()
    test_signature_detector()
    test_density_helper()


# ------------------------------------------------- grading
def grade(res):
    sig = {s: lag2_signature(res[s]["x_base2"]["lags"]) for s in res}

    p_sf = {"squarefree_signature_present": sig["squarefree"]["signature_present"],
            "prediction": "present (control replication)",
            "holds": bool(sig["squarefree"]["signature_present"]),
            "detail": sig["squarefree"]}

    triple_pred = {}
    for tag, primes in TRIPLES.items():
        key = f"squarefree_except_{tag}"
        triple_pred[f"P_D7_triple_{tag}"] = {
            f"{key}_signature_present": sig[key]["signature_present"],
            "prediction": f"present (three-prime relaxation; r expected below 0.1 for triple {primes})",
            "holds": bool(sig[key]["signature_present"]),
            "detail": sig[key],
        }

    p_odd = {"odd_squarefree_signature_present": sig["odd_squarefree"]["signature_present"],
             "prediction": "absent (removes even integers -> all gaps even, anchor)",
             "holds": bool(not sig["odd_squarefree"]["signature_present"]),
             "detail": sig["odd_squarefree"]}

    p_pr = {"primes_signature_present": sig["primes"]["signature_present"],
            "prediction": "absent (negative control)",
            "holds": bool(not sig["primes"]["signature_present"]),
            "detail": sig["primes"]}

    r_control = sig["squarefree"]["ratio_lag2_over_median_others"]

    all_triples_present = bool(
        all(sig[f"squarefree_except_{tag}"]["signature_present"] for tag in TRIPLES))
    all_triples_below_threshold = bool(
        all(sig[f"squarefree_except_{tag}"]["ratio_lag2_over_median_others"] <= THRESHOLD
            for tag in TRIPLES))
    convergence_holds = bool(
        sig["squarefree"]["signature_present"]
        and all_triples_present
        and not sig["odd_squarefree"]["signature_present"]
        and not sig["primes"]["signature_present"]
    )
    some_triple_destroyed = bool(
        any(not sig[f"squarefree_except_{tag}"]["signature_present"] for tag in TRIPLES))

    # ---- THREE-prime ADDITIVITY: k({p,q,s}) = (r({p,q,s}) - r_control) / (d(p)+d(q)+d(s)).
    ratio_by_triple = {}
    triple_density = {}
    k_by_triple = {}
    for tag, primes in TRIPLES.items():
        r_triple = sig[f"squarefree_except_{tag}"]["ratio_lag2_over_median_others"]
        d_triple = sum(density(p) for p in primes)
        ratio_by_triple[tag] = r_triple
        triple_density[tag] = d_triple
        k_by_triple[tag] = (r_triple - r_control) / d_triple
    lo, hi = ADDITIVITY_LAW_BAND
    additivity_law_holds = bool(all(lo <= k_by_triple[tag] <= hi for tag in TRIPLES))

    out = {
        "P_D7_sf_control": p_sf,
    }
    out.update(triple_pred)
    out.update({
        "P_D7_odd_anchor": p_odd,
        "P_D7_primes_negative_control": p_pr,
        "convergence_holds": convergence_holds,
        "some_triple_destroyed": some_triple_destroyed,
        "all_triples_present": all_triples_present,
        "all_triples_below_threshold": all_triples_below_threshold,
        "threshold": THRESHOLD,
        "density_additivity_law": {
            "model": "r({p,q,s}) - r_control ~ C * (d(p)+d(q)+d(s)), d(p) = 1/p^2 - 1/p^3",
            "r_control_measured": r_control,
            "triples": {tag: list(TRIPLES[tag]) for tag in TRIPLES},
            "ratio_by_triple": ratio_by_triple,
            "summed_density_by_triple": triple_density,
            "coefficient_k_by_triple": k_by_triple,
            "preregistered_band": list(ADDITIVITY_LAW_BAND),
            "additivity_law_holds": additivity_law_holds,
        },
        "signature_all_sieves": sig,
    })
    return out


# ------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=20_000_000)
    ap.add_argument("--out-dir", default="/out")
    args = ap.parse_args()
    N = args.n_max
    os.makedirs(args.out_dir, exist_ok=True)

    print("[self-test] running ...", flush=True)
    try:
        run_selftests()
    except AssertionError as e:
        json.dump({"self_tests": {"all_passed": False, "error": str(e)}},
                  open(os.path.join(args.out_dir, "results.json"), "w"), indent=2)
        print("SELF-TEST FAILED:", e, file=sys.stderr)
        sys.exit(1)
    print("[self-test] passed", flush=True)

    res = {}
    counts = {}
    collectors = (
        ("squarefree", lambda n: collect_exponent_cap_gaps(n, {})),
        ("squarefree_except_17_19_23", lambda n: collect_exponent_cap_gaps(n, {17: 2, 19: 2, 23: 2})),
        ("squarefree_except_29_31_37", lambda n: collect_exponent_cap_gaps(n, {29: 2, 31: 2, 37: 2})),
        ("squarefree_except_41_43_47", lambda n: collect_exponent_cap_gaps(n, {41: 2, 43: 2, 47: 2})),
        ("odd_squarefree", lambda n: collect_exponent_cap_gaps(n, {2: 0})),
        ("primes", lambda n: collect_prime_gaps(n)),
    )
    for key, fn in collectors:
        print(f"[sieve] {key} up to {N} ...", flush=True)
        gaps, cnt, maxg = fn(N)
        print(f"[sieve] {key}: count={cnt} gaps={len(gaps)} max_gap={maxg}", flush=True)
        res[key] = analyze_sieve(gaps, maxg)
        counts[key] = {"count": cnt, "n_gaps": len(gaps), "max_gap": maxg}
        for name in STREAMS:
            zs = " ".join(
                f"z{k}={res[key][name]['lags'].get(str(k), {}).get('z', float('nan')):+9.2f}"
                for k in LAGS)
            print(f"    {key:28s} {name:14s} {zs}", flush=True)
        del gaps

    grading = grade(res)
    out = {
        "title": "Axis D7: does the Axis-D5/D6 density scaling law stay ADDITIVE for THREE primes? "
                 "Relaxing the strict squarefree exponent cap at three primes {p,q,s} at once should "
                 "perturb the ratio by the SUM of the single-prime densities, "
                 "r({p,q,s}) - r_control ~ C*(d(p)+d(q)+d(s)).",
        "n_max": N,
        "lags": LAGS,
        "identity": "P_2(g) = P_4(g) XOR P_4(floor(g/2))",
        "sieves": {
            "squarefree": "A005117 (control, c_p=1; expect present)",
            "squarefree_except_17_19_23": "c_17=c_19=c_23=2 (allow 289|n,361|n,529|n; expect present)",
            "squarefree_except_29_31_37": "c_29=c_31=c_37=2 (allow 841|n,961|n,1369|n; expect present)",
            "squarefree_except_41_43_47": "c_41=c_43=c_47=2 (allow 1681|n,1849|n,2209|n; expect present)",
            "odd_squarefree": "A056911, c_2=0 (anchor: removes even, expect absent)",
            "primes": "A000040 (negative control, expect absent)",
        },
        "null": "exact marginal-preserving permutation null of S^(k), reported as z",
        "signature": "lag-2 unique-quiet: argmin_k|z_x(k)| == 2 AND |z_x(2)| <= 0.1 x "
                     "median|z_x(k!=2)| (scale-free ratio test)",
        "counts": counts,
        "per_sieve": res,
        "grading": grading,
        "self_tests": {"all_passed": True},
        "produced_by": "AI agents working as MV2, Munim, Inc.",
    }
    with open(os.path.join(args.out_dir, "results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("[grading]", json.dumps(grading, indent=2), flush=True)
    print("[done] wrote results.json", flush=True)


if __name__ == "__main__":
    main()
