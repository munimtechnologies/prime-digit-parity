#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Produced with AI by MV2 of Munim, Inc.
"""Axis D3: is single-prime exponent-cap relaxation destructive at EVERY prime, or does the
base-2 gap-parity "lag-2 unique-quiet" signature survive above a large-prime crossover?

Axis D (v0.24.0) and its generality replication Axis D2 (v0.26.0) showed that relaxing the
squarefree exponent cap at a single prime p (admit p^2|n but not p^3|n) destroys the signature at
every prime tested: 2, 3, 5, 7. Axis D2 also surfaced a MONOTONE trend in the committed frozen
ratio r(p) = |z(2)| / median|z(k!=2)| for squarefree-except-p:
    r(2)=1.142, r(3)=0.773, r(5)=0.342, r(7)=0.169   (control squarefree r=0.015, PRESENT).
As p grows, squarefree-except-p -> squarefree (admitted class {p^2 m : p does not divide m} has
density 1/p^2 - 1/p^3 -> 0), so r(p) must converge to the control 0.015 and should eventually
cross BELOW the 0.1 threshold, at which point the detector would score the signature PRESENT. The
strict-mu^2 "every prime" reading instead predicts r(p) asymptotes ABOVE 0.1. These disagree.

This job pushes to the next two primes, 11 and 13, to test which. Five sets to a common cutoff N:
  * squarefree             (control)      : c_p = 1 all p                 (A005117, 6/pi^2)
  * squarefree_except_11   (PRIMARY)      : c_11 = 2, c_p = 1 otherwise   (allows 121|n, not 1331|n)
  * squarefree_except_13   (replication)  : c_13 = 2, c_p = 1 otherwise   (allows 169|n, not 2197|n)
  * odd_squarefree         (anchor)       : c_2 = 0, c_p = 1 odd p        (A056911, all gaps even)
  * primes                 (neg control)  : A000040

Pre-registered in prime-digit-parity/notes/axis-d3-large-prime-crossover-strict-mu2-prereg.md
(committed before this job runs). Primary predictions (strict-mu^2 "every prime" reading):
P-D3-sf PRESENT, P-D3-ex11 (PRIMARY) ABSENT, P-D3-ex13 ABSENT, P-D3-odd ABSENT, P-D3-primes ABSENT.
Secondary quantitative: the monotone decrease continues (r(11) < 0.169, r(13) < r(11)) while both
stay ABOVE 0.1. The sharp alternative CROSSOVER: if r(11) <= 0.1 or r(13) <= 0.1, that sieve
scores PRESENT, falsifying "every single prime destroys it" -- reported honestly as the headline.

The observable, exact permutation null `null_z`, base-4 stride split, stream analysis, scale-free
detector, general exponent-cap sieve and prime-gap collector are reused VERBATIM (and
re-self-tested) from axis_d2_generality_job.py / axis_d_exponent_cap_lag_signature_job.py. The
exponent-cap sieve is self-tested against a full factorization over [1, 20000] (count, the complete
gap stream and the max gap, tiny segment forcing boundary crossings) for caps {}, {2:0}, {11:2},
{13:2} and {11:2, 13:2}, plus the mechanism contrast (except-11 contains 121 but not 1331;
except-13 contains 169 but not 2197; both admit even members and produce odd gaps) and the exact
counts Q(10^6)=607926, pi(10^6)=78498. Self-tests run first; the job exits non-zero on any failure.
Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.

v2: identical to v1 except one cosmetic fix in grade() -- the v1 run completed the full
computation and all self-tests but crashed assembling the grading dict, because
`dict(PRIOR_RATIOS, **{11: r11, 13: r13})` requires string keyword keys; replaced with a dict
literal `{**PRIOR_RATIOS, 11: r11, 13: r13}`, which allows the integer keys. No measurement,
sieve, null or detector logic changed.
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

# Committed frozen ratios from Axis D (v0.24.0) and Axis D2 (v0.26.0), for the pre-registered
# monotone-trend check. These are published numbers, not recomputed here.
PRIOR_RATIOS = {2: 1.142, 3: 0.773, 5: 0.342, 7: 0.169}
THRESHOLD = 0.1


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
    divisible by p^(cap_p + 1) for any p. caps={} -> squarefree (A005117); {11:2} ->
    squarefree-except-11 (allow 121|n, forbid 1331|n); {13:2} -> squarefree-except-13 (allow
    169|n, forbid 2197|n); {2:0} -> odd squarefree (A056911). We mark multiples of each forbidden
    prime power q = p^(cap_p + 1) <= N; a surviving position is a member."""
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
    assert not _in_set(4, {}) and not _in_set(121, {})          # 2^2, 11^2
    assert _in_set(121, {11: 2}) and _in_set(242, {11: 2})      # 121|n now allowed
    assert not _in_set(1331, {11: 2})                           # but 1331 = 11^3 still forbidden
    assert not _in_set(4, {11: 2})                              # unrelated square still forbidden
    assert _in_set(169, {13: 2}) and _in_set(338, {13: 2})      # 169|n allowed
    assert not _in_set(2197, {13: 2}) and not _in_set(4, {13: 2})  # 2197 and 4 still forbidden
    assert not _in_set(2, {2: 0}) and _in_set(3, {2: 0})        # odd squarefree excludes evens
    # Exact reference match (count, FULL gap stream, max gap) over [1, 20000] for each cap set,
    # with a deliberately tiny segment so many segment boundaries are crossed.
    LIMIT = 20000
    for caps in ({}, {2: 0}, {11: 2}, {13: 2}, {11: 2, 13: 2}):
        members = _ref_members(LIMIT, caps)
        ref_gaps = [members[i + 1] - members[i] for i in range(len(members) - 1)]
        gaps, cnt, maxg = collect_exponent_cap_gaps(LIMIT, caps, seg_size=97)
        assert cnt == len(members), (caps, cnt, len(members))
        assert list(gaps) == ref_gaps, caps
        assert maxg == (max(ref_gaps) if ref_gaps else 0), (caps, maxg)
    # MECHANISM CONTRAST: the enlarging sets admit even integers AND produce at least one odd
    # gap (half-shifted stride phase survives), unlike the even-removing anchor.
    for caps in ({11: 2}, {13: 2}):
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


def run_selftests():
    digit_parity_tests()
    test_stride_identity()
    test_lag_stats()
    test_null_formula()
    test_exponent_cap_sieve()
    test_primes()
    test_signature_detector()


# ------------------------------------------------- grading
def grade(res):
    sieves = list(res.keys())
    sig = {s: lag2_signature(res[s]["x_base2"]["lags"]) for s in sieves}

    p_sf = {"squarefree_signature_present": sig["squarefree"]["signature_present"],
            "prediction": "present (control replication)",
            "holds": bool(sig["squarefree"]["signature_present"]),
            "detail": sig["squarefree"]}

    # P-D3-ex11 (PRIMARY): following the strict-mu^2 reading, relaxing the cap at 11 is predicted
    # to DESTROY it (ABSENT); the CROSSOVER alternative is PRESENT (r(11) <= 0.1).
    p_ex11 = {"squarefree_except_11_signature_present": sig["squarefree_except_11"]["signature_present"],
              "prediction": "absent (relax exponent cap at 11 to {0,1,2}; strict-mu^2 predicted)",
              "holds": bool(not sig["squarefree_except_11"]["signature_present"]),
              "detail": sig["squarefree_except_11"]}

    # P-D3-ex13 (replication): relaxing the cap at 13 is predicted to DESTROY it.
    p_ex13 = {"squarefree_except_13_signature_present": sig["squarefree_except_13"]["signature_present"],
              "prediction": "absent (relax exponent cap at 13 to {0,1,2}; strict-mu^2 predicted)",
              "holds": bool(not sig["squarefree_except_13"]["signature_present"]),
              "detail": sig["squarefree_except_13"]}

    p_odd = {"odd_squarefree_signature_present": sig["odd_squarefree"]["signature_present"],
             "prediction": "absent (removes even integers -> all gaps even, anchor)",
             "holds": bool(not sig["odd_squarefree"]["signature_present"]),
             "detail": sig["odd_squarefree"]}

    p_pr = {"primes_signature_present": sig["primes"]["signature_present"],
            "prediction": "absent (negative control)",
            "holds": bool(not sig["primes"]["signature_present"]),
            "detail": sig["primes"]}

    r11 = sig["squarefree_except_11"]["ratio_lag2_over_median_others"]
    r13 = sig["squarefree_except_13"]["ratio_lag2_over_median_others"]

    # Decisive: strict-mu^2 destruction PERSISTS at 11 and 13 iff control present while both
    # large-prime relaxations (and anchor and neg control) are absent.
    destruction_persists = bool(
        sig["squarefree"]["signature_present"]
        and not sig["squarefree_except_11"]["signature_present"]
        and not sig["squarefree_except_13"]["signature_present"]
        and not sig["odd_squarefree"]["signature_present"]
        and not sig["primes"]["signature_present"]
    )
    # CROSSOVER: the sharp surprise -- a large single-prime relaxation that KEEPS the signature.
    crossover_detected = bool(
        sig["squarefree_except_11"]["signature_present"]
        or sig["squarefree_except_13"]["signature_present"]
    )
    # Secondary quantitative trend vs the committed r(7)=0.169.
    monotone_trend_continues = bool(r11 < PRIOR_RATIOS[7] and r13 < r11)
    both_remain_above_threshold = bool(r11 > THRESHOLD and r13 > THRESHOLD)

    return {
        "P_D3_sf_control": p_sf,
        "P_D3_ex11_primary": p_ex11,
        "P_D3_ex13_replication": p_ex13,
        "P_D3_odd_anchor": p_odd,
        "P_D3_primes_negative_control": p_pr,
        "strict_mu2_destruction_persists_at_11_and_13": destruction_persists,
        "crossover_detected": crossover_detected,
        "ratio_by_relaxed_prime": {**PRIOR_RATIOS, 11: r11, 13: r13},
        "monotone_trend_continues": monotone_trend_continues,
        "both_remain_above_threshold": both_remain_above_threshold,
        "threshold": THRESHOLD,
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
        ("squarefree_except_11", lambda n: collect_exponent_cap_gaps(n, {11: 2})),
        ("squarefree_except_13", lambda n: collect_exponent_cap_gaps(n, {13: 2})),
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
        "title": "Axis D3: does single-prime exponent-cap relaxation destroy the base-2 gap-parity "
                 "lag-2 unique-quiet signature at the larger primes 11 and 13, or does the ratio "
                 "cross below the 0.1 threshold (signature survives)?",
        "n_max": N,
        "lags": LAGS,
        "identity": "P_2(g) = P_4(g) XOR P_4(floor(g/2))",
        "sieves": {
            "squarefree": "A005117 (control, c_p=1; expect present)",
            "squarefree_except_11": "c_11=2, c_p=1 else (PRIMARY: allow 121|n not 1331|n; expect absent)",
            "squarefree_except_13": "c_13=2, c_p=1 else (allow 169|n not 2197|n; expect absent)",
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
