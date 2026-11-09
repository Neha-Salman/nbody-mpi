#!/usr/bin/env python3
"""Plot body trajectories (x-y projection) from results/trajectory_n100.csv."""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")  # save the figure to a file; no window needed
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "..")
SOURCE = os.path.join(ROOT, "results", "trajectory_n100.csv")
OUTPUT = os.path.join(ROOT, "figures", "trajectories_n100.png")
SHOW = 10  # number of bodies to draw

# Collect (x, y) positions for the first SHOW bodies, in step order
tracks = {}
with open(SOURCE) as f:
    for row in csv.DictReader(f):
        i = int(row["id"])
        if i < SHOW:
            tracks.setdefault(i, []).append((float(row["x"]), float(row["y"])))

fig, ax = plt.subplots(figsize=(7, 7))
for i, points in sorted(tracks.items()):
    xs, ys = zip(*points)
    line, = ax.plot(xs, ys, linewidth=1, label=f"body {i}")
    colour = line.get_color()
    ax.plot(xs[0], ys[0], "o", color=colour)   # start position
    ax.plot(xs[-1], ys[-1], "s", color=colour)  # end position

ax.set_xlabel("x (simulation units)")
ax.set_ylabel("y (simulation units)")
ax.set_title(f"Trajectories of {SHOW} bodies, N=100, 1000 steps (circle = start, square = end)")
ax.set_aspect("equal")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=7, loc="upper right")
fig.tight_layout()
fig.savefig(OUTPUT, dpi=150)
plt.close(fig)
print("Saved", OUTPUT)
