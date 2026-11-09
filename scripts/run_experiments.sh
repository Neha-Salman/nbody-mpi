#!/bin/bash
# Reproducible experiments for the CT-488 CCA.
# Runs the sequential baseline and the MPI version for every body count
# and process count. Each configuration is repeated REPEATS times.
set -e

# Move to the project root, whatever folder the script is run from
cd "$(dirname "$0")/.."

STEPS=100
DT=0.0001
REPEATS=3
LOG=results/experiments.log

gcc -O2 -Wall -Wextra -o nbody_seq src/nbody_seq.c -lm
mpicc -O2 -Wall -Wextra -o nbody_mpi src/nbody_mpi.c -lm

: > "$LOG"   # clear the log before a fresh run

for N in 100 500 1000 5000 10000; do
    for rep in $(seq 1 "$REPEATS"); do
        echo "SEQ N=$N rep=$rep $(./nbody_seq "$N" "$STEPS" "$DT")" | tee -a "$LOG"
    done
    for P in 1 2 4 8; do
        for rep in $(seq 1 "$REPEATS"); do
            echo "MPI P=$P N=$N rep=$rep $(mpirun --oversubscribe -np "$P" ./nbody_mpi "$N" "$STEPS" "$DT")" | tee -a "$LOG"
        done
    done
done

echo "Experiments finished. Results saved to $LOG"
