#!/usr/bin/env python3
"""
compare_benchmarks.py - Compare Python vs R: precision and performance.

Usage:
    python3 benchmark_python.py   # step 1: generate data + run Python
    Rscript benchmark_R.R         # step 2: run R on same data
    python3 compare_benchmarks.py # step 3: compare
"""

import numpy as np
import json
import csv
import os
import sys


def load_pvals(path):
    """Load p-values from CSV (handles both plain and header formats)."""
    with open(path) as f:
        first = f.readline().strip()
    if first.startswith('"') or first.startswith("p") or first.startswith("f"):
        return np.loadtxt(path, delimiter=",", skiprows=1, usecols=-1)
    return np.loadtxt(path, delimiter=",")


def precision_report(py_vals, r_vals, label):
    """Print precision comparison for one set of values."""
    n = len(py_vals)
    abs_diff = np.abs(py_vals - r_vals)
    # Relative diff (avoid div by zero)
    denom = np.maximum(np.abs(r_vals), 1e-20)
    rel_diff = abs_diff / denom

    print(f"  {label} ({n} values)")
    print(f"    Max |diff|:       {abs_diff.max():.4e}")
    print(f"    Mean |diff|:      {abs_diff.mean():.4e}")
    print(f"    Median |diff|:    {np.median(abs_diff):.4e}")
    print(f"    Max rel diff:     {rel_diff.max() * 100:.6f}%")
    print(f"    Mean rel diff:    {rel_diff.mean() * 100:.6f}%")
    print(f"    Correlation:      {np.corrcoef(py_vals, r_vals)[0,1]:.15f}")
    print(f"    Exact match (12 dp): {np.sum(np.abs(py_vals - r_vals) < 1e-12)}/{n}")

    # Tail agreement
    for alpha in [0.05, 0.01, 0.001]:
        py_sig = py_vals < alpha
        r_sig = r_vals < alpha
        agree = np.sum(py_sig == r_sig)
        py_n = py_sig.sum()
        r_n = r_sig.sum()
        print(f"    Sig agreement (p<{alpha}): {agree}/{n}  "
              f"(Py={py_n}, R={r_n}, both={np.sum(py_sig & r_sig)})")
    print()


def main():
    # Check files exist
    required = {
        "bench_results_python.json": "Run: python3 benchmark_python.py",
        "bench_results_R.csv":       "Run: Rscript benchmark_R.R",
        "bench_bb_pvals_py.csv":     "Run: python3 benchmark_python.py",
        "bench_bb_pvals_R.csv":      "Run: Rscript benchmark_R.R",
    }
    for f, hint in required.items():
        if not os.path.exists(f):
            print(f"ERROR: {f} not found. {hint}")
            sys.exit(1)

    # ---- Load timing results ----
    with open("bench_results_python.json") as f:
        py_timing = json.load(f)

    r_timing = {}
    with open("bench_results_R.csv") as f:
        reader = csv.DictReader(f)
        row = next(reader)
        r_timing = {k: float(v) for k, v in row.items() if "sec" in k}
        r_timing.update({k: int(float(v)) for k, v in row.items() if "sig" in k})

    # ---- Load p-values ----
    bb_py  = load_pvals("bench_bb_pvals_py.csv")
    bb_r   = load_pvals("bench_bb_pvals_R.csv")
    ibb_py = load_pvals("bench_ibb_pvals_py.csv")
    ibb_r  = load_pvals("bench_ibb_pvals_R.csv")
    fc_py  = load_pvals("bench_ibb_fc_py.csv")
    fc_r   = load_pvals("bench_ibb_fc_R.csv")

    n_genes = len(bb_py)

    # ================================================================
    print("=" * 65)
    print(f"  BENCHMARK COMPARISON: Python vs R countdata")
    print(f"  {n_genes} genes x 200 samples, 48 threads")
    print("=" * 65)
    print()

    # ---- Performance ----
    print("-" * 65)
    print("PERFORMANCE")
    print("-" * 65)

    py_bb  = py_timing["bb_sec"]
    py_ibb = py_timing["ibb_sec"]
    r_bb   = r_timing["bb_sec"]
    r_ibb  = r_timing["ibb_sec"]

    print(f"  {'':15} {'Python':>10} {'R':>10} {'Ratio (R/Py)':>14}")
    print(f"  {'':15} {'------':>10} {'------':>10} {'------------':>14}")
    print(f"  {'bb.test':15} {py_bb:>9.3f}s {r_bb:>9.3f}s {r_bb/py_bb:>13.2f}x")
    print(f"  {'ibb.test':15} {py_ibb:>9.3f}s {r_ibb:>9.3f}s {r_ibb/py_ibb:>13.2f}x")
    print(f"  {'total':15} {py_bb+py_ibb:>9.3f}s {r_bb+r_ibb:>9.3f}s "
          f"{(r_bb+r_ibb)/(py_bb+py_ibb):>13.2f}x")
    print()

    ratio = (r_bb + r_ibb) / (py_bb + py_ibb)
    if ratio > 1.05:
        print(f"  -> Python is ~{ratio:.1f}x FASTER than R")
    elif ratio < 0.95:
        print(f"  -> R is ~{1/ratio:.1f}x FASTER than Python")
    else:
        print(f"  -> Similar performance (ratio={ratio:.2f}x)")
    print()

    # ---- Precision ----
    print("-" * 65)
    print("PRECISION")
    print("-" * 65)
    print()

    precision_report(bb_py, bb_r, "bb.test p-values")
    precision_report(ibb_py, ibb_r, "ibb.test p-values")
    precision_report(fc_py, fc_r, "ibb.test fold-change")

    # ---- Summary table ----
    print("-" * 65)
    print("SUMMARY")
    print("-" * 65)
    print(f"  bb.test:  Py {py_bb:.3f}s vs R {r_bb:.3f}s | "
          f"max p-val diff = {np.max(np.abs(bb_py - bb_r)):.2e} | "
          f"corr = {np.corrcoef(bb_py, bb_r)[0,1]:.12f}")
    print(f"  ibb.test: Py {py_ibb:.3f}s vs R {r_ibb:.3f}s | "
          f"max p-val diff = {np.max(np.abs(ibb_py - ibb_r)):.2e} | "
          f"corr = {np.corrcoef(ibb_py, ibb_r)[0,1]:.12f}")
    print()


if __name__ == "__main__":
    main()