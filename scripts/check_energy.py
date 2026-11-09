#!/usr/bin/env python3
"""Check energy conservation of the sequential N-body simulation (N = 100).

Compares total energy (kinetic + potential) at the start with the energy
after 100 and after 1000 steps, using the state files written by nbody_seq.
"""
import csv
import math
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "..")
G = 1.0
SOFTENING = 1e-3

def load(path):
    with open(os.path.join(ROOT, path)) as f:
        return list(csv.DictReader(f))

def total_energy(rows):
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
    return kinetic, potential, kinetic + potential

e0 = total_energy(load("results/seq_init.csv"))
for label, path in [("after 100 steps", "results/seq_n100.csv"),
                    ("after 1000 steps", "results/traj_final_n100.csv")]:
    ek, ep, et = total_energy(load(path))
    drift = (et - e0[2]) / abs(e0[2])
    print(f"{label}: kinetic={ek:.6f} potential={ep:.6f} total={et:.6f} "
          f"relative drift={drift:.3e}")
print(f"initial: kinetic={e0[0]:.6f} potential={e0[1]:.6f} total={e0[2]:.6f}")
