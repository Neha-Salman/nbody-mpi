#!/bin/bash
# Reproducible sequential baseline experiments for the CT-488 CCA.
# Compiles the program, then runs it for each required body count.
set -e

# Move to the project root, whatever folder the script is run from
cd "$(dirname "$0")/.."

STEPS=100
DT=0.01
LOG=results/sequential_timings.log

gcc -O2 -Wall -Wextra -o nbody_seq src/nbody_seq.c -lm

: > "$LOG"   # clear the log before a fresh run
for N in 100 500 1000 5000 10000; do
    ./nbody_seq "$N" "$STEPS" "$DT" "results/seq_n${N}.csv" | tee -a "$LOG"
done
echo "Sequential experiments finished. Log saved to $LOG"
