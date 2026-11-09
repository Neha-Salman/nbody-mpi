# Distributed N-Body Simulation using MPI (CT-488 CCA)

This project implements a direct O(N²) gravitational N-body simulation in two versions:

- src/nbody_seq.c: sequential baseline, used for correctness and performance comparison
- src/nbody_mpi.c: distributed version using MPI, with block data decomposition

The project measures computation time, communication time, total time, speedup, parallel efficiency, and scalability for body counts from 100 to 10,000 and process counts from 1 to 8.

## Requirements

Tested on Ubuntu 26.04 LTS running in VMware Workstation, with 4 virtual CPUs.

~~~bash
sudo apt update
sudo apt install -y build-essential openmpi-bin libopenmpi-dev git python3-matplotlib
~

## Project structure

~
src/        C source code (sequential and MPI)
scripts/    Scripts to build, run experiments, analyse results, and plot
results/    Experiment log, summary table, and verification outputs
figures/    Performance graphs (PNG)
docs/       Report and presentation
~

## Build

~~~bash
gcc -O2 -Wall -Wextra -o nbody_seq src/nbody_seq.c -lm
mpicc -O2 -Wall -Wextra -o nbody_mpi src/nbody_mpi.c -lm
~

## Running

Sequential run:

~~~bash
./nbody_seq N STEPS DT [output.csv]
~

MPI run with P processes:

~~~bash
mpirun -np P ./nbody_mpi N STEPS DT [output.csv]
~

Example: 100 bodies, 100 steps, time step 0.01, on 4 processes:

~~~bash
mpirun -np 4 ./nbody_mpi 100 100 0.01 results/mpi_p4_n100.csv
~

## Parameters

| Parameter | Meaning |
|---|---|
| N | Number of bodies |
| STEPS | Number of time steps |
| DT | Size of each time step |
| P | Number of MPI processes (must be at most N) |
| output.csv | Optional file for the final state of each body |

Fixed constants in the code: gravitational constant G = 1.0, softening 1e-3, and random seed 42. The seed makes every run start from identical initial conditions.

## Reproducing the experiments

Run the full experiment (75 measured runs, about 15 minutes):

~~~bash
./scripts/run_experiments.sh
~

This writes every run to results/experiments.log.

Summarise the results (median of 3 repeats, speedup, efficiency):

~~~bash
./scripts/analyze.sh
~

Output goes to results/summary.txt.

Create the performance graphs:

~~~bash
python3 scripts/plot_results.py
~

Output goes to figures/.

Run the sequential baseline alone:

~~~bash
./scripts/run_sequential.sh
~

## Verification

The MPI results were compared with the sequential results for N = 100 after 100 steps. The maximum absolute difference was 0 for both 1 and 4 processes, within the stated tolerance of 1e-10.

~~~bash
paste -d, results/seq_n100.csv results/mpi_p4_n100.csv | awk -F, 'NR>1 {for (k=2; k<=7; k++) {d=$k-$(k+8); if (d<0) d=-d; if (d>m) m=d}} END {printf "max abs difference = %.3e\n", m}'
~

## Notes on the measurements

- The VM has 4 virtual CPUs. The host processor (Intel i7-7600U) has 2 physical cores with hyperthreading, so runs with more than 2 processes share physical cores and speedup levels off.
- Runs with 8 processes oversubscribe the 4 virtual CPUs and use --oversubscribe.
- Timings on a virtual machine are noisy. Each configuration is repeated 3 times and the median is reported.
- The sequential baseline stores bodies as structs, while the MPI version uses packed arrays. This makes the MPI version slightly faster on one process, so small-N speedups are not reliable.

## Licence

Academic project for CT-488 Distributed Computing, NED University of Engineering & Technology.
