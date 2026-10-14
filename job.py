"""Run the prime digit-sum parity study and write results to --out-dir.

Usage:
    python job.py --out-dir OUT [--limit 1000000000]

Outputs in OUT:
    results.json   full checkpoint data (value cutoffs and count cutoffs)
    summary.json   derived z-scores and the value-vs-count comparison
    summary.txt    human-readable summary

The headline question (see NOTES.md / RESOLUTION.md): the finite-sample base-10
prime digit-sum parity bias reported by Gupta & Sridhar (arXiv:1812.11841) was
measured at a PRIME-COUNT cutoff (first 50,000,000 primes).  Integer digit-sum
parity is exactly balanced only at VALUE cutoffs x = b**k.  Does the prime bias
look the same at value cutoffs x = 10**k as at count cutoffs, and does the
fair-coin z grow (consistent with the Mauduit-Rivat power-saving error term)?
This job measures both so RESULTS can answer it with verified numbers.

Work by Munim, Inc.  Produced with AI (see README.md).  Apache-2.0.
"""

from __future__ import annotations

import argparse
import json
import math
import os

import core  # self-check reference (stdlib)
import fast


def self_check():
    """Cross-check the fast path against the stdlib reference on a small range."""
    import numpy as np

    pr = np.array(core.primes_upto(20000), dtype=np.int64)
    for b in (2, 4, 6, 8, 10, 16):
        fast_odd = int(fast.parity_base(pr, b).sum())
        ref_odd, _ = core.prime_parity_counts(20000, b)
        assert fast_odd == ref_odd, f"mismatch base {b}: {fast_odd} != {ref_odd}"
    # even-base exact balance on a complete range, via the reference
    for b in (2, 10, 16):
        odd, even = core.integer_parity_counts(b**4, b)
        assert odd == even, f"range not balanced base {b}"
    return True


def z_scores(rows, bases):
    """Add a 'z' field: excess / sqrt(n) (fair-coin standard deviation is
    sqrt(n), since excess = odd - even = 2*odd - n and a balanced coin gives
    Var(2*odd - n) = n)."""
    out = []
    for r in rows:
        n = r["n"]
        z = {}
        for b in bases:
            z[b] = (r["excess"][b] / math.sqrt(n)) if n > 0 else 0.0
        out.append({"x": r.get("x"), "n": n, "excess": r["excess"], "z": z})
    return out


def max_abs_z(rows, bases, min_n=10_000):
    """Per base, the checkpoint with the largest |z| among rows with n >= min_n."""
    out = {}
    for b in bases:
        best = None
        for r in rows:
            if r["n"] >= min_n:
                z = r["z"][b]
                if best is None or abs(z) > abs(best["z"]):
                    best = {"z": z, "n": r["n"], "x": r["x"]}
        out[b] = best
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--limit", type=int, default=1_000_000_000)
    ap.add_argument("--count-step", type=int, default=1_000_000)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    self_check()

    bases = [2, 4, 6, 8, 10, 16]
    # value checkpoints: powers of ten plus a few linear points
    vpts = []
    for k in range(1, 11):
        vpts.append(10**k)
    for m in (2, 3, 5, 7):
        for k in (7, 8, 9):
            vpts.append(m * 10**k)
    vpts = sorted(set(v for v in vpts if v <= args.limit))

    res = fast.run(args.limit, bases, vpts, args.count_step)

    with open(os.path.join(args.out_dir, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1)

    summary = {
        "limit": args.limit,
        "bases": bases,
        "value_cutoffs_z": z_scores(res["value_checkpoints"], bases),
        "count_cutoffs_z": z_scores(res["count_checkpoints"], bases),
        "trajectory_z": z_scores(res["segment_checkpoints"], bases),
        "max_abs_z": max_abs_z(
            z_scores(res["count_checkpoints"] + res["segment_checkpoints"], bases),
            bases,
        ),
        "final": res["final"],
        "note": (
            "excess = #(odd digit-sum primes) - #(even). z = excess/sqrt(n); "
            "|z| ~ O(1) would be fair-coin behaviour. Compare base-10 z at "
            "value cutoffs x=10^k vs at prime-count cutoffs."
        ),
    }
    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    lines = [f"prime digit-sum parity, primes up to {args.limit}", ""]
    lines.append("VALUE cutoffs x = 10^k  (excess=odd-even, z=excess/sqrt(n)):")
    for r in summary["value_cutoffs_z"]:
        if r["x"] in (10**k for k in range(1, 11)):
            b10 = r["excess"][10]
            z10 = r["z"][10]
            b2 = r["excess"][2]
            z2 = r["z"][2]
            lines.append(
                f"  x={r['x']:>13} n={r['n']:>10}  base10 excess={b10:>8} "
                f"z={z10:+.2f}  base2 excess={b2:>8} z={z2:+.2f}"
            )
    lines.append("")
    lines.append("VALUE cutoffs x = 10^k, z for every base (cross-base table):")
    lines.append("  " + "x".rjust(14) + "".join(f"  base{b:>3}" for b in bases))
    for r in summary["value_cutoffs_z"]:
        if r["x"] in (10**k for k in range(1, 11)):
            lines.append(
                "  " + str(r["x"]).rjust(14)
                + "".join(f"  {r['z'][b]:>+7.2f}" for b in bases)
            )
    lines.append("")
    lines.append("Largest |z| seen (count and segment checkpoints, n >= 10000):")
    for b in bases:
        m = summary["max_abs_z"][b]
        if m is not None:
            lines.append(f"  base {b:>2}: z={m['z']:+.2f} at n={m['n']} (x={m['x']})")
    lines.append("")
    lines.append("COUNT cutoffs (first N primes):")
    for r in summary["count_cutoffs_z"]:
        if r["n"] % 10_000_000 == 0:
            lines.append(
                f"  N={r['n']:>10} (p={r['x']:>13})  base10 excess="
                f"{r['excess'][10]:>8} z={r['z'][10]:+.2f}"
            )
    with open(os.path.join(args.out_dir, "summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
