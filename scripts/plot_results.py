#!/usr/bin/env python3
"""Create the performance graphs for the CT-488 CCA from results/summary.txt."""
import os
import sys
import matplotlib
matplotlib.use("Agg")  # save files directly; no window needed
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "..")
SUMMARY = os.path.join(ROOT, "results", "summary.txt")
FIG_DIR = os.path.join(ROOT, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

SIZES = [100, 500, 1000, 5000, 10000]
PROCS = [1, 2, 4, 8]

# Read the summary table: sequential medians and MPI medians per (N, P)
seq = {}
mpi = {}
with open(SUMMARY) as f:
    for line in f:
        parts = line.split()
        if len(parts) != 7 or parts[0] == "N":
            continue
        n = int(parts[0])
        if parts[1] == "seq":
            seq[n] = float(parts[2])
        else:
            p = int(parts[1])
            mpi[(n, p)] = {
                "total": float(parts[2]),
                "speedup": float(parts[3]),
                "eff": float(parts[4]),
                "comp": float(parts[5]),
                "comm": float(parts[6]),
            }

def save(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("Saved", path)

# Figure 1: speedup vs. processes
fig, ax = plt.subplots(figsize=(7, 5))
for n in SIZES:
    ps = [p for p in PROCS if (n, p) in mpi]
    ax.plot(ps, [mpi[(n, p)]["speedup"] for p in ps], marker="o", label=f"N={n}")
ax.plot(PROCS, PROCS, "k--", label="Ideal (speedup = P)")
ax.set_xlabel("Number of MPI processes (P)")
ax.set_ylabel("Speedup over sequential baseline")
ax.set_title("Speedup vs. number of processes")
ax.set_xticks(PROCS)
ax.grid(True, alpha=0.3)
ax.legend()
save(fig, "speedup_vs_processes.png")

# Figure 2: efficiency vs. processes
fig, ax = plt.subplots(figsize=(7, 5))
for n in SIZES:
    ps = [p for p in PROCS if (n, p) in mpi]
    ax.plot(ps, [mpi[(n, p)]["eff"] for p in ps], marker="o", label=f"N={n}")
ax.axhline(1.0, color="k", linestyle="--", label="Ideal (efficiency = 1)")
ax.set_xlabel("Number of MPI processes (P)")
ax.set_ylabel("Parallel efficiency")
ax.set_title("Parallel efficiency vs. number of processes")
ax.set_xticks(PROCS)
ax.grid(True, alpha=0.3)
ax.legend()
save(fig, "efficiency_vs_processes.png")

# Figure 3: computation vs. communication at N = 10000
n = 10000
labels = [f"P={p}" for p in PROCS]
comp = [mpi[(n, p)]["comp"] for p in PROCS]
comm = [mpi[(n, p)]["comm"] for p in PROCS]
fig, ax = plt.subplots(figsize=(7, 5))
ax.bar(labels, comp, label="Computation (slowest process)")
ax.bar(labels, comm, bottom=comp, label="Communication (slowest process)")
ax.set_ylabel("Time (s)")
ax.set_title(f"Computation vs. communication time, N={n}")
ax.grid(axis="y", alpha=0.3)
ax.legend()
save(fig, "computation_vs_communication_n10000.png")

# Figure 4: sequential runtime vs. N, with an O(N^2) reference line
ns = sorted(seq)
times = [seq[k] for k in ns]
ref = [seq[1000] * (k / 1000) ** 2 for k in ns]
fig, ax = plt.subplots(figsize=(7, 5))
ax.loglog(ns, times, "o-", label="Measured (sequential)")
ax.loglog(ns, ref, "k--", label="O(N^2) reference")
ax.set_xlabel("Number of bodies (N)")
ax.set_ylabel("Time (s)")
ax.set_title("Sequential runtime vs. N (log-log)")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
save(fig, "sequential_scaling_loglog.png")
