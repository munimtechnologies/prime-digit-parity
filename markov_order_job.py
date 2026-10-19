"""Sandbox job: is the prime-gap digit-sum parity sequence first-order Markov?

For every prime gap up to N_MAX, and for each base in {2,4,8,10,12,16}, we form
the base-b digit-sum parity bit of the gap and ask whether the resulting binary
sequence is first-order Markov or carries genuine higher-order memory, using the
triple conditional-independence test H0: s_{n+1} _||_ s_{n-1} | s_n.

Two readings are produced per base:
  * global  : one 2x2x2 table over the whole sequence, df = 2, exact p = e^{-G2/2};
  * blocked : the sequence is cut into `blocks` contiguous pieces and the G2 are
              summed (df = 2*blocks), controlling for the slow drift of the gap
              distribution (gaps grow like log x), so a rejection cannot be a mere
              non-stationarity artifact.
We also report r1, r2 and the first-order prediction r1^2 as an intuitive
companion (a first-order Markov sequence satisfies r2 = r1^2).

Self-tests run first and the job exits non-zero if any fails, so only
self-validated numbers are ever emitted. Produced with AI by MV2, Munim, Inc.
Apache-2.0.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

import markov_order as mo


def sieve_np(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def np_triple_counts(bits):
    b = bits.astype(np.int64)
    code = 4 * b[:-2] + 2 * b[1:-1] + b[2:]
    return np.bincount(code, minlength=8)


def np_lag_k(bits, k):
    u = 1.0 - 2.0 * bits.astype(np.float64)
    u -= u.mean()
    denom = float(np.dot(u, u))
    if denom == 0.0:
        return 0.0
    return float(np.dot(u[:-k], u[k:]) / denom)


def blocked_g2(bits, blocks):
    """Sum of the per-block triple G2 statistics; df = 2 * (#usable blocks)."""
    n = bits.size
    edges = np.linspace(0, n, blocks + 1, dtype=np.int64)
    g2_sum = 0.0
    used = 0
    for i in range(blocks):
        seg = bits[edges[i]:edges[i + 1]]
        if seg.size < 3:
            continue
        res = mo.g_test_triples(np_triple_counts(seg).tolist())
        g2_sum += res["G2"]
        used += 1
    df = 2 * used
    p = mo.chi2_sf(g2_sum, df) if df > 0 else 1.0
    return {"G2": g2_sum, "df": df, "p": p, "blocks_used": used}


def gen_markov1(length, p1, rng):
    """First-order binary Markov chain. p1[b] = P(next=1 | current=b)."""
    x = np.empty(length, dtype=np.int8)
    x[0] = 0
    r = rng.random(length)
    p0, p1_ = p1[0], p1[1]
    for i in range(1, length):
        x[i] = 1 if r[i] < (p1_ if x[i - 1] else p0) else 0
    return x


def gen_markov2(length, q, rng):
    """Second-order binary chain. q[(a,b)] = P(next=1 | prev2=a, prev1=b)."""
    x = np.empty(length, dtype=np.int8)
    x[0], x[1] = 0, 1
    r = rng.random(length)
    for i in range(2, length):
        x[i] = 1 if r[i] < q[(int(x[i - 2]), int(x[i - 1]))] else 0
    return x


def run_selftests():
    checks = {}
    # 1. Sieve vs known prime counts.
    checks["pi_1e6"] = int(sieve_np(1_000_000).size) == 78498
    checks["pi_1e4"] = int(sieve_np(10_000).size) == 1229
    # 2. numpy triple counts match the pure-Python reference.
    rng = np.random.default_rng(12345)
    rb = rng.integers(0, 2, size=5000).astype(np.int8)
    checks["triple_crosscheck"] = (
        np_triple_counts(rb).tolist() == mo.triple_counts(rb.tolist()))
    # 3. numpy autocorr matches the pure-Python reference.
    checks["lag1_crosscheck"] = abs(
        np_lag_k(rb, 1) - mo.lag_k_autocorr(rb.tolist(), 1)) < 1e-9
    checks["lag2_crosscheck"] = abs(
        np_lag_k(rb, 2) - mo.lag_k_autocorr(rb.tolist(), 2)) < 1e-9
    # 4. chi2 survival closed forms.
    checks["chi2_df2"] = all(
        abs(mo.chi2_sf(x, 2) - np.exp(-x / 2.0)) < 1e-9
        for x in (1.0, 5.0, 30.0))
    checks["chi2_df200"] = 0.0 < mo.chi2_sf(220.0, 200) < 1.0
    # 5. FALSE-POSITIVE gate: a genuinely first-order chain must NOT be flagged.
    m1 = gen_markov1(500_000, {0: 0.3, 1: 0.7}, np.random.default_rng(7))
    g1 = mo.g_test_triples(np_triple_counts(m1).tolist())
    checks["markov1_not_rejected"] = g1["G2"] < 40.0
    # 6. POWER gate: a genuine second-order chain must be flagged.
    q = {(0, 0): 0.2, (0, 1): 0.5, (1, 0): 0.5, (1, 1): 0.8}
    m2 = gen_markov2(500_000, q, np.random.default_rng(9))
    g2 = mo.g_test_triples(np_triple_counts(m2).tolist())
    checks["markov2_rejected"] = g2["G2"] > 100.0
    checks["markov1_G2"] = g1["G2"]
    checks["markov2_G2"] = g2["G2"]
    passed = all(v is True for k, v in checks.items() if not k.endswith("_G2"))
    checks["all_passed"] = passed
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=100_000_000)
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--blocks", type=int, default=100)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, "results.json")

    t0 = time.time()
    selftests = run_selftests()
    if not selftests["all_passed"]:
        json.dump({"self_tests": selftests, "error": "self-tests failed"},
                  open(out_path, "w"), indent=2)
        print("SELF-TESTS FAILED:", selftests)
        sys.exit(1)

    primes = sieve_np(args.n_max)
    gaps = np.diff(primes).astype(np.int64)
    max_gap = int(gaps.max())
    bases = [2, 4, 8, 10, 12, 16]

    per_base = {}
    for base in bases:
        tbl = np.array(mo.parity_table(max_gap, base), dtype=np.int8)
        bits = tbl[gaps]
        counts = np_triple_counts(bits)
        glob = mo.g_test_triples(counts.tolist())
        blk = blocked_g2(bits, args.blocks)
        r1 = np_lag_k(bits, 1)
        r2 = np_lag_k(bits, 2)
        r1sq = r1 * r1
        per_base[str(base)] = {
            "triple_counts": counts.tolist(),
            "n_triples": int(bits.size - 2),
            "global": glob,
            "blocked": blk,
            "r1": r1,
            "r2": r2,
            "r1_squared": r1sq,
            "r2_minus_r1sq": r2 - r1sq,
            "frac_ones": float(bits.mean()),
        }

    g2_by_base = {b: per_base[b]["global"]["G2"] for b in per_base}
    argmax_base = max(g2_by_base, key=g2_by_base.get)

    results = {
        "title": "Markov-order control: is prime-gap digit-sum parity first-order Markov?",
        "n_max": args.n_max,
        "n_primes": int(primes.size),
        "n_gaps": int(gaps.size),
        "max_gap": max_gap,
        "bases": bases,
        "blocks": args.blocks,
        "per_base": per_base,
        "global_G2_by_base": g2_by_base,
        "base_with_largest_global_G2": argmax_base,
        "self_tests": selftests,
        "elapsed_sec": round(time.time() - t0, 2),
        "produced_by": "AI agents working as MV2, Munim, Inc.",
    }
    json.dump(results, open(out_path, "w"), indent=2)
    print(json.dumps({k: results[k] for k in
                      ("n_max", "n_primes", "n_gaps", "max_gap",
                       "global_G2_by_base", "base_with_largest_global_G2",
                       "elapsed_sec")}, indent=2))


if __name__ == "__main__":
    main()
