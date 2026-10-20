#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Produced with AI by MV2 of Munim, Inc.
"""Axis D generality: does the STRICT-squarefree (mu^2) requirement for the base-2 gap-parity
"lag-2 unique-quiet" signature hold when the exponent cap is relaxed at a LARGER prime?

Axis D (v0.24.0) showed that relaxing the squarefree exponent cap at prime 2 (squarefree-except-2,
allow 4|n) or 3 (squarefree-except-3, allow 9|n) DESTROYS the signature even though both sets keep
every even integer. Conclusion: the signature needs the strict squarefree congruence structure --
it tolerates neither thinning a residue class (Axis C) nor enlarging the exponent cap (Axis D).

Axis D used only the two smallest odd-relevant primes (2, 3). This replication asks whether the
destruction is specific to those primes or is a general property of relaxing the cap at ANY single
prime, by repeating the exponent-cap relaxation at the larger primes 5 and 7. Five sets to a
common cutoff N:
  * squarefree            (control)      : c_p = 1 all p                 (A005117, 6/pi^2)
  * squarefree_except_5   (PRIMARY)      : c_5 = 2, c_p = 1 otherwise    (allows 25|n, not 125|n)
  * squarefree_except_7   (replication)  : c_7 = 2, c_p = 1 otherwise    (allows 49|n, not 343|n)
  * odd_squarefree        (anchor)       : c_2 = 0, c_p = 1 odd p        (A056911, all gaps even)
  * primes                (neg control)  : A000040

Pre-registered in prime-digit-parity/notes/axis-d2-generality-strict-mu2-prereg.md (committed
before this job runs). Pre-registered predictions, following the Axis D result: P-D2-sf PRESENT,
P-D2-ex5 (PRIMARY) ABSENT, P-D2-ex7 ABSENT, P-D2-odd ABSENT, P-D2-primes ABSENT -- i.e. the
strict-mu^2 requirement generalizes and relaxing the cap at 5 or 7 destroys the signature just as
it did at 2 and 3. A PRESENT outcome at except-5 and/or except-7 would be a genuine surprise,
reopening the question of whether the destruction is specific to the small primes 2/3, and would
be reported honestly as such.

The observable (base-2 gap digit-sum parity = Thue-Morse), the exact marginal-preserving
permutation null `null_z`, the base-4 stride split (P_2(g)=P_4(g) XOR P_4(g>>1)), the stream
analysis, the scale-free lag-2 signature detector, the general exponent-cap sieve and the
prime-gap collector are reused VERBATIM (and re-self-tested) from axis_d_exponent_cap_lag_signature_job.py.
The exponent-cap sieve is self-tested against an obviously-correct full factorization over
[1, 20000] (count, the complete gap stream and the max gap, with a tiny segment forcing many
boundary crossings) for caps {}, {2:0}, {5:2}, {7:2} and {5:2, 7:2}, plus the mechanism contrast
(except-5/except-7 admit even members and odd gaps; except-5 contains 25 but not 125; except-7
contains 49 but not 343) and the exact counts Q(10^6)=607926 and pi(10^6)=78498. Self-tests run
first; the job exits non-zero if any fails, so only self-validated numbers are emitted.
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
    divisible by p^(cap_p + 1) for any p. caps={} -> squarefree (A005117); {5:2} ->
    squarefree-except-5 (allow 25|n, forbid 125|n); {7:2} -> squarefree-except-7 (allow 49|n,
    forbid 343|n); {2:0} -> odd squarefree (A056911). We mark multiples of each forbidden prime
    power q = p^(cap_p + 1) <= N; a surviving position is a member."""
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
    present = bool(quiet_lag == 2 and az[2] <= 0.1 * med_others)
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
    assert not _in_set(4, {}) and not _in_set(25, {})          # 2^2, 5^2
    assert _in_set(25, {5: 2}) and _in_set(50, {5: 2})         # 25|n now allowed
    assert not _in_set(125, {5: 2})                            # but 125 = 5^3 still forbidden
    assert not _in_set(4, {5: 2})                              # unrelated square still forbidden
    assert _in_set(49, {7: 2}) and _in_set(98, {7: 2})         # 49|n allowed
    assert not _in_set(343, {7: 2}) and not _in_set(4, {7: 2})  # 343 and 4 still forbidden
    assert not _in_set(2, {2: 0}) and _in_set(3, {2: 0})       # odd squarefree excludes evens
    # Exact reference match (count, FULL gap stream, max gap) over [1, 20000] for each cap set,
    # with a deliberately tiny segment so many segment boundaries are crossed.
    LIMIT = 20000
    for caps in ({}, {2: 0}, {5: 2}, {7: 2}, {5: 2, 7: 2}):
        members = _ref_members(LIMIT, caps)
        ref_gaps = [members[i + 1] - members[i] for i in range(len(members) - 1)]
        gaps, cnt, maxg = collect_exponent_cap_gaps(LIMIT, caps, seg_size=97)
        assert cnt == len(members), (caps, cnt, len(members))
        assert list(gaps) == ref_gaps, caps
        assert maxg == (max(ref_gaps) if ref_gaps else 0), (caps, maxg)
    # MECHANISM CONTRAST: the enlarging sets admit even integers AND produce at least one odd
    # gap (half-shifted stride phase survives), unlike the even-removing anchor.
    for caps in ({5: 2}, {7: 2}):
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

    # P-D2-ex5 (PRIMARY): following Axis D, relaxing the cap at 5 is predicted to DESTROY it.
    p_ex5 = {"squarefree_except_5_signature_present": sig["squarefree_except_5"]["signature_present"],
             "prediction": "absent (relax exponent cap at 5 to {0,1,2}; strict-mu^2 predicted)",
             "holds": bool(not sig["squarefree_except_5"]["signature_present"]),
             "detail": sig["squarefree_except_5"]}

    # P-D2-ex7 (replication): relaxing the cap at 7 is predicted to DESTROY it.
    p_ex7 = {"squarefree_except_7_signature_present": sig["squarefree_except_7"]["signature_present"],
             "prediction": "absent (relax exponent cap at 7 to {0,1,2}; strict-mu^2 predicted)",
             "holds": bool(not sig["squarefree_except_7"]["signature_present"]),
             "detail": sig["squarefree_except_7"]}

    # P-D2-odd (anchor): odd squarefree does NOT show the signature.
    p_odd = {"odd_squarefree_signature_present": sig["odd_squarefree"]["signature_present"],
             "prediction": "absent (removes even integers -> all gaps even, anchor)",
             "holds": bool(not sig["odd_squarefree"]["signature_present"]),
             "detail": sig["odd_squarefree"]}

    # P-D2-primes (negative control): primes do NOT show it.
    p_pr = {"primes_signature_present": sig["primes"]["signature_present"],
            "prediction": "absent (negative control)",
            "holds": bool(not sig["primes"]["signature_present"]),
            "detail": sig["primes"]}

    # Decisive read: the strict-mu^2 requirement GENERALIZES to larger primes iff the control is
    # present while BOTH large-prime exponent relaxations (and the anchor and neg control) are
    # absent. If instead either except-5 or except-7 is PRESENT, the destruction was specific to
    # the small primes 2/3 and the strict-mu^2 account does not generalize (a surprise).
    strict_mu2_generalizes = bool(
        sig["squarefree"]["signature_present"]
        and not sig["squarefree_except_5"]["signature_present"]
        and not sig["squarefree_except_7"]["signature_present"]
        and not sig["odd_squarefree"]["signature_present"]
        and not sig["primes"]["signature_present"]
    )

    return {
        "P_D2_sf_control": p_sf,
        "P_D2_ex5_primary": p_ex5,
        "P_D2_ex7_replication": p_ex7,
        "P_D2_odd_anchor": p_odd,
        "P_D2_primes_negative_control": p_pr,
        "strict_mu2_requirement_generalizes_to_larger_primes": strict_mu2_generalizes,
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
        ("squarefree_except_5", lambda n: collect_exponent_cap_gaps(n, {5: 2})),
        ("squarefree_except_7", lambda n: collect_exponent_cap_gaps(n, {7: 2})),
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
        "title": "Axis D generality: does relaxing the squarefree exponent cap at a LARGER prime "
                 "(squarefree-except-5, squarefree-except-7) destroy the base-2 gap-parity "
                 "lag-2 unique-quiet signature, as it does at 2 and 3?",
        "n_max": N,
        "lags": LAGS,
        "identity": "P_2(g) = P_4(g) XOR P_4(floor(g/2))",
        "sieves": {
            "squarefree": "A005117 (control, c_p=1; expect present)",
            "squarefree_except_5": "c_5=2, c_p=1 else (PRIMARY: allow 25|n not 125|n; expect absent)",
            "squarefree_except_7": "c_7=2, c_p=1 else (allow 49|n not 343|n; expect absent)",
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
