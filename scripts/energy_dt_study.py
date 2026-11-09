#!/usr/bin/env python3
"""Energy drift versus time step for N=100, all runs covering simulated time 1.0."""
import csv
import math
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "..")
EXE = os.path.join(ROOT, "nbody_seq")
G = 1.0
SOFTENING = 1e-3
N = 100
OUT_DIR = os.path.join(ROOT, "results")

def total_energy(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    pos = [(float(r["x"]), float(r["y"]), float(r["z"])) for r in rows]
    vel = [(float(r["vx"]), float(r["vy"]), float(r["vz"])) for r in rows]
    mass = [float(r["m"]) for r in rows]
    n = len(mass)
    kinetic = sum(0.5 * mass[i] * (vel[i][0]**2 + vel[i][1]**2 + vel[i][2]**2)
                  for i in range(n))
    potential = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dx = pos[j][0] - pos[i][0]
            dy = pos[j][1] - pos[i][1]
            dz = pos[j][2] - pos[i][2]
            r = math.sqrt(dx * dx + dy * dy + dz * dz + SOFTENING * SOFTENING)
            potential -= G * mass[i] * mass[j] / r
    return kinetic + potential

def run(steps, dt, name):
    out = os.path.join(OUT_DIR, name)
    subprocess.run([EXE, str(N), str(steps), str(dt), out],
                   check=True, capture_output=True)
    return out

e_start = total_energy(run(0, 0.01, "dt_study_init.csv"))
print(f"Initial total energy: {e_start:.6f}")
print(f"{'dt':>10} {'steps':>8} {'final energy':>14} {'relative drift':>16}")
for dt in [0.01, 0.001, 0.0001]:
    steps = round(1.0 / dt)
    e_end = total_energy(run(steps, dt, f"dt_study_{dt:g}.csv"))
    drift = (e_end - e_start) / abs(e_start)
    print(f"{dt:>10g} {steps:>8d} {e_end:>14.6f} {drift:>16.3e}")
