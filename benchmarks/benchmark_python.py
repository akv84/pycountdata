#!/usr/bin/env python3
"""
benchmark_python.py - pycountdata benchmark on simulated RNA-seq data.
5000 genes x 200 samples (100/group), 48 threads.

Run this FIRST — it generates the shared input data.

Usage:
    conda activate countdata
    python3 benchmark_python.py
"""

import numpy as np
import time
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from countdata import bb_test, ibb_test

N_GENES = 2000
N_SAMPLES = 40
N_PER_GROUP = 20
N_THREADS = 48

np.random.seed(2024)

print(f"Simulating RNA-seq: {N_GENES} genes x {N_SAMPLES} samples ({N_PER_GROUP}/group)")
print(f"Threads: {N_THREADS}")
print()

# --- Simulate realistic RNA-seq counts ---
# Gene expression means: log-normal, typical RNA-seq range (1-2000 counts)
gene_means = np.exp(np.random.normal(2.5, 1.5, N_GENES)).clip(1, 2000)

# Library sizes: 10M-30M total reads per sample
library_sizes = np.random.uniform(1e7, 3e7, N_SAMPLES)

# Gene-level counts: NB with gene-specific mean and dispersion
# Counts represent reads mapping to each gene (typically 0 - few thousand)
x = np.zeros((N_GENES, N_SAMPLES), dtype=float)
for i in range(N_GENES):
    mu = gene_means[i] * (library_sizes / library_sizes.mean())
    # BCV^2 ~ 0.1 + 1/sqrt(mu) typical for RNA-seq
    disp = 0.1 + 10.0 / (1 + gene_means[i])
    r = 1.0 / disp
    p = np.clip(r / (r + mu), 1e-10, 1 - 1e-10)
    x[i] = np.random.negative_binomial(r, p)

# tx = library sizes (total reads per sample, same across genes)
# This ensures x[i,j] << tx[i,j] as in real RNA-seq
tx = np.tile(library_sizes, (N_GENES, 1))
group = np.array([1] * N_PER_GROUP + [2] * N_PER_GROUP)

# Save shared input data for R
np.savetxt("bench_x.csv", x, delimiter=",", fmt="%.0f")
np.savetxt("bench_tx.csv", tx, delimiter=",", fmt="%.2f")
np.savetxt("bench_group.csv", group.reshape(1, -1), delimiter=",", fmt="%d")

print(f"Data: bench_x.csv, bench_tx.csv, bench_group.csv")
print(f"Counts: min={x.min():.0f}, median={np.median(x):.0f}, max={x.max():.0f}")
print()

# --- bb.test ---
print("=" * 50)
print("bb.test (unpaired)")
print("=" * 50)
t0 = time.perf_counter()
bb_res = bb_test(x, tx, group, n_threads=N_THREADS, verbose=True)
bb_time = time.perf_counter() - t0

bb_p = bb_res["p.value"]
print(f"Time:       {bb_time:.3f} s")
print(f"Sig p<0.05: {np.sum(bb_p < 0.05)}/{N_GENES}")
print(f"Sig p<0.01: {np.sum(bb_p < 0.01)}/{N_GENES}")
print(f"First 5 p:  {bb_p[:5]}")
print()

# --- ibb.test ---
print("=" * 50)
print("ibb.test (paired)")
print("=" * 50)
t0 = time.perf_counter()
ibb_res = ibb_test(x, tx, group, n_threads=N_THREADS, verbose=True)
ibb_time = time.perf_counter() - t0

ibb_p = ibb_res["p.value"]
ibb_fc = ibb_res["fc"]
print(f"Time:       {ibb_time:.3f} s")
print(f"Sig p<0.05: {np.sum(ibb_p < 0.05)}/{N_GENES}")
print(f"Sig p<0.01: {np.sum(ibb_p < 0.01)}/{N_GENES}")
print(f"First 5 p:  {ibb_p[:5]}")
print(f"First 5 fc: {ibb_fc[:5]}")
print()

# --- Save p-values for comparison ---
np.savetxt("bench_bb_pvals_py.csv", bb_p, fmt="%.15e")
np.savetxt("bench_ibb_pvals_py.csv", ibb_p, fmt="%.15e")
np.savetxt("bench_ibb_fc_py.csv", ibb_fc, fmt="%.15e")

with open("bench_results_python.json", "w") as f:
    json.dump({"bb_sec": round(bb_time, 4), "ibb_sec": round(ibb_time, 4),
               "bb_sig005": int(np.sum(bb_p < 0.05)),
               "ibb_sig005": int(np.sum(ibb_p < 0.05))}, f, indent=2)

print(f"Summary: bb={bb_time:.3f}s  ibb={ibb_time:.3f}s")