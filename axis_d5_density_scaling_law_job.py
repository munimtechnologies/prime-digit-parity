#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Produced with AI by MV2 of Munim, Inc.
"""Axis D5: does the single-prime exponent-cap relaxation keep converging to the squarefree control
at the next primes 29, 31, 37 -- and does the super-control excess obey a quantitative DENSITY
SCALING LAW?

Axis D3 (v0.27.0) found a crossover at prime 11; Axis D4 (v0.28.0) confirmed convergence-to-control
at 17, 19, 23 with no second crossover. Across the nine relaxed primes measured so far the frozen
ratio r(p) = |z(2)| / median|z(k!=2)| for squarefree-except-p is:
    r(2)=1.142, r(3)=0.773, r(5)=0.342, r(7)=0.169, r(11)=0.0777, r(13)=0.0603,
    r(17)=0.0420, r(19)=0.0366, r(23)=0.0291   (control r_control ~ 0.0149).
The admitted extra class {p^2 m : p does not divide m} has density d(p) = 1/p^2 - 1/p^3. Computing
the dimensionless coefficient k(p) = (r(p) - r_control) / d(p) over those nine primes gives a
strikingly stable number (~8-10, slowly declining): a candidate QUANTITATIVE LAW
    r(p) - r_control  ~  C * d(p),   C ~ 8.
This job pushes to the next three primes 29, 31, 37 to test TWO pre-registered things at once:
  (1) convergence continues -- all three PRESENT, r still falling monotonically, no second crossover;
  (2) the density scaling law holds out-of-sample -- k(29), k(31), k(37) each land in the committed
      band [5, 11] (consistent with the historical k range 7.85-10.2).

Six sets to a common cutoff N:
  * squarefree             (control)      : c_p = 1 all p                 (A005117, 6/pi^2)
  * squarefree_except_29                  : c_29 = 2, c_p = 1 otherwise   (allows 841|n, not 24389|n)
  * squarefree_except_31                  : c_31 = 2, c_p = 1 otherwise   (allows 961|n, not 29791|n)
  * squarefree_except_37                  : c_37 = 2, c_p = 1 otherwise   (allows 1369|n, not 50653|n)
  * odd_squarefree         (anchor)       : c_2 = 0, c_p = 1 odd p        (A056911, all gaps even)
  * primes                 (neg control)  : A000040

Pre-registered in prime-digit-parity/notes/axis-d5-density-scaling-law-strict-mu2-prereg.md
(committed before this job runs). Primary predictions: P-D5-sf PRESENT, P-D5-ex29/31/37 PRESENT,
P-D5-odd ABSENT, P-D5-primes ABSENT; monotone continuation r(29) < 0.0291, r(31) < r(29),
r(37) < r(31). Secondary (first quantitative/predictive test of the arc): density_law_holds, i.e.
k(29), k(31), k(37) all in [5, 11]. Named alternatives (reported honestly): a SECOND CROSSOVER /
RE-DESTRUCTION (some r rises above 0.1), or the density law FAILING (some k outside [5, 11]) while
the signatures still hold -- which would mean the convergence is real but not simply density-linear.

The observable, exact permutation null `null_z`, base-4 stride split, stream analysis, scale-free
detector, general exponent-cap sieve and prime-gap collector are reused VERBATIM (and
re-self-tested) from axis_d4_crossover_convergence_job.py. The exponent-cap sieve is self-tested
against a full factorization over [1, 20000] (count, the complete gap stream and the max gap, a tiny
segment forcing boundary crossings) for caps {}, {2:0}, {29:2}, {31:2}, {37:2} and
{29:2, 31:2, 37:2}, plus the mechanism contrast (except-p contains p^2 but not p^3; all admit even
members and produce odd gaps) and the exact counts Q(10^6)=607926, pi(10^6)=78498. Self-tests run
first; the job exits non-zero on any failure.
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

# Committed frozen ratios from Axis D (v0.24.0), D2 (v0.26.0), D3 (v0.27.0) and D4 (v0.28.0), for the
# pre-registered monotone-trend and density-law checks. These are published numbers, not recomputed.
PRIOR_RATIOS = {2: 1.142, 3: 0.773, 5: 0.342, 7: 0.169,
                11: 0.07770937284253937, 13: 0.06033990564490313,
                17: 0.042025554336370244, 19: 0.03664851374134167,
                23: 0.029122410406888552}
THRESHOLD = 0.1
# Pre-registered band for the density-law coefficient k(p) = (r(p) - r_control) / (1/p^2 - 1/p^3),
# chosen to bracket the historical range (7.85-10.2) with margin.
DENSITY_LAW_BAND = (5.0, 11.0)


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
    """Segmented sieve of {n <= N : v_p(n) <= caps.get(p, 1) for every prime p}, i.e. n not
    divisible by p^(cap_p + 1) for any p. caps={} -> squarefree (A005117); {29:2} ->
    squarefree-except-29 (allow 841|n, forbid 24389|n); {2:0} -> odd squarefree (A056911). We mark
    multiples of each forbidden prime power q = p^(cap_p + 1) <= N; a surviving position is a
    member."""
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
    # Membership checker sanity against hand facts.
    assert _in_set(1, {}) and _in_set(2, {}) and _in_set(30, {})
    assert not _in_set(4, {}) and not _in_set(841, {})          # 2^2, 29^2
    assert _in_set(841, {29: 2}) and _in_set(1682, {29: 2})     # 841|n now allowed
    assert not _in_set(24389, {29: 2})                          # but 24389 = 29^3 still forbidden
    assert not _in_set(4, {29: 2})                              # unrelated square still forbidden
    assert _in_set(961, {31: 2}) and _in_set(1922, {31: 2})     # 961|n allowed
    assert not _in_set(29791, {31: 2}) and not _in_set(4, {31: 2})  # 29791=31^3 and 4 still forbidden
    assert _in_set(1369, {37: 2}) and _in_set(2738, {37: 2})    # 1369|n allowed
    assert not _in_set(50653, {37: 2}) and not _in_set(4, {37: 2})  # 50653=37^3 and 4 still forbidden
    assert not _in_set(2, {2: 0}) and _in_set(3, {2: 0})        # odd squarefree excludes evens
    # Exact reference match (count, FULL gap stream, max gap) over [1, 20000] for each cap set,
    # with a deliberately tiny segment so many segment boundaries are crossed.
    LIMIT = 20000
    for caps in ({}, {2: 0}, {29: 2}, {31: 2}, {37: 2}, {29: 2, 31: 2, 37: 2}):
        members = _ref_members(LIMIT, caps)
        ref_gaps = [members[i + 1] - members[i] for i in range(len(members) - 1)]
        gaps, cnt, maxg = collect_exponent_cap_gaps(LIMIT, caps, seg_size=97)
        assert cnt == len(members), (caps, cnt, len(members))
        assert list(gaps) == ref_gaps, caps
        assert maxg == (max(ref_gaps) if ref_gaps else 0), (caps, maxg)
    # MECHANISM CONTRAST: the enlarging sets admit even integers AND produce at least one odd
    # gap (half-shifted stride phase survives), unlike the even-removing anchor.
    for caps in ({29: 2}, {31: 2}, {37: 2}):
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
    # d(p) = 1/p^2 - 1/p^3, a positive decreasing function of p.
    assert abs(density(2) - (0.25 - 0.125)) < 1e-15
    assert density(3) < density(2) and density(37) < density(29)


def run_selftests():
    digit_parity_tests()
    test_stride_identity()
    test_lag_stats()
    test_null_formula()
    test_exponent_cap_sieve()
    test_primes()
    test_signature_detector()
    test_density_helper()


# ------------------------------------------------- density scaling law
def density(p):
    """Density of the admitted extra class {p^2 m : p does not divide m}: 1/p^2 - 1/p^3."""
    return 1.0 / (p * p) - 1.0 / (p * p * p)


# ------------------------------------------------- grading
def grade(res):
    sieves = list(res.keys())
    sig = {s: lag2_signature(res[s]["x_base2"]["lags"]) for s in sieves}

    p_sf = {"squarefree_signature_present": sig["squarefree"]["signature_present"],
            "prediction": "present (control replication)",
            "holds": bool(sig["squarefree"]["signature_present"]),
            "detail": sig["squarefree"]}

    # Following Axis D3/D4 convergence-to-control, each larger-prime relaxation is predicted PRESENT.
    p_ex29 = {"squarefree_except_29_signature_present": sig["squarefree_except_29"]["signature_present"],
              "prediction": "present (convergence-to-control; r(29) expected below 0.1)",
              "holds": bool(sig["squarefree_except_29"]["signature_present"]),
              "detail": sig["squarefree_except_29"]}

    p_ex31 = {"squarefree_except_31_signature_present": sig["squarefree_except_31"]["signature_present"],
              "prediction": "present (convergence-to-control; r(31) expected below 0.1)",
              "holds": bool(sig["squarefree_except_31"]["signature_present"]),
              "detail": sig["squarefree_except_31"]}

    p_ex37 = {"squarefree_except_37_signature_present": sig["squarefree_except_37"]["signature_present"],
              "prediction": "present (convergence-to-control; r(37) expected below 0.1)",
              "holds": bool(sig["squarefree_except_37"]["signature_present"]),
              "detail": sig["squarefree_except_37"]}

    p_odd = {"odd_squarefree_signature_present": sig["odd_squarefree"]["signature_present"],
             "prediction": "absent (removes even integers -> all gaps even, anchor)",
             "holds": bool(not sig["odd_squarefree"]["signature_present"]),
             "detail": sig["odd_squarefree"]}

    p_pr = {"primes_signature_present": sig["primes"]["signature_present"],
            "prediction": "absent (negative control)",
            "holds": bool(not sig["primes"]["signature_present"]),
            "detail": sig["primes"]}

    r_control = sig["squarefree"]["ratio_lag2_over_median_others"]
    r29 = sig["squarefree_except_29"]["ratio_lag2_over_median_others"]
    r31 = sig["squarefree_except_31"]["ratio_lag2_over_median_others"]
    r37 = sig["squarefree_except_37"]["ratio_lag2_over_median_others"]

    all_three_below_threshold = bool(r29 <= THRESHOLD and r31 <= THRESHOLD and r37 <= THRESHOLD)
    # Monotone continuation vs the committed r(23)=0.0291.
    monotone_trend_continues = bool(r29 < PRIOR_RATIOS[23] and r31 < r29 and r37 < r31)
    convergence_to_control_confirmed = bool(
        sig["squarefree"]["signature_present"]
        and sig["squarefree_except_29"]["signature_present"]
        and sig["squarefree_except_31"]["signature_present"]
        and sig["squarefree_except_37"]["signature_present"]
        and not sig["odd_squarefree"]["signature_present"]
        and not sig["primes"]["signature_present"]
        and monotone_trend_continues
    )
    second_crossover_or_redestruction = bool(
        not sig["squarefree_except_29"]["signature_present"]
        or not sig["squarefree_except_31"]["signature_present"]
        or not sig["squarefree_except_37"]["signature_present"]
    )

    # ---- DENSITY SCALING LAW: k(p) = (r(p) - r_control) / d(p), d(p) = 1/p^2 - 1/p^3.
    ratio_by_prime = {**PRIOR_RATIOS, 29: r29, 31: r31, 37: r37}
    k_by_prime = {p: (ratio_by_prime[p] - r_control) / density(p) for p in ratio_by_prime}
    k29, k31, k37 = k_by_prime[29], k_by_prime[31], k_by_prime[37]
    lo, hi = DENSITY_LAW_BAND
    density_law_holds = bool(lo <= k29 <= hi and lo <= k31 <= hi and lo <= k37 <= hi)

    return {
        "P_D5_sf_control": p_sf,
        "P_D5_ex29": p_ex29,
        "P_D5_ex31": p_ex31,
        "P_D5_ex37": p_ex37,
        "P_D5_odd_anchor": p_odd,
        "P_D5_primes_negative_control": p_pr,
        "convergence_to_control_confirmed": convergence_to_control_confirmed,
        "second_crossover_or_redestruction": second_crossover_or_redestruction,
        "ratio_by_relaxed_prime": ratio_by_prime,
        "monotone_trend_continues": monotone_trend_continues,
        "all_three_below_threshold": all_three_below_threshold,
        "threshold": THRESHOLD,
        "density_scaling_law": {
            "model": "r(p) - r_control ~ C * (1/p^2 - 1/p^3)",
            "r_control_measured": r_control,
            "density_by_prime": {str(p): density(p) for p in ratio_by_prime},
            "coefficient_k_by_prime": {str(p): k_by_prime[p] for p in ratio_by_prime},
            "preregistered_band": list(DENSITY_LAW_BAND),
            "k_29_31_37": [k29, k31, k37],
            "density_law_holds": density_law_holds,
        },
        "signature_all_sieves": sig,
    }


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
        ("squarefree_except_29", lambda n: collect_exponent_cap_gaps(n, {29: 2})),
        ("squarefree_except_31", lambda n: collect_exponent_cap_gaps(n, {31: 2})),
        ("squarefree_except_37", lambda n: collect_exponent_cap_gaps(n, {37: 2})),
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
            print(f"    {key:22s} {name:14s} {zs}", flush=True)
        del gaps

    grading = grade(res)
    out = {
        "title": "Axis D5: does the single-prime exponent-cap relaxation keep converging to the "
                 "squarefree control at primes 29, 31, 37, and does the super-control excess "
                 "r(p) - r_control obey a density scaling law ~ C*(1/p^2 - 1/p^3)?",
        "n_max": N,
        "lags": LAGS,
        "identity": "P_2(g) = P_4(g) XOR P_4(floor(g/2))",
        "sieves": {
            "squarefree": "A005117 (control, c_p=1; expect present)",
            "squarefree_except_29": "c_29=2, c_p=1 else (allow 841|n not 24389|n; expect present)",
            "squarefree_except_31": "c_31=2, c_p=1 else (allow 961|n not 29791|n; expect present)",
            "squarefree_except_37": "c_37=2, c_p=1 else (allow 1369|n not 50653|n; expect present)",
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
