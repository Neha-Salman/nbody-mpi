/*
 * nbody_mpi.c - Distributed direct O(N^2) N-body simulation using MPI.
 *
 * Usage: mpirun -np P ./nbody_mpi N STEPS DT [output.csv]
 *   P          number of MPI processes (must be <= N)
 *   N          number of bodies
 *   STEPS      number of time steps
 *   DT         size of each time step
 *   output.csv optional file that receives the final state of every body
 *
 * Decomposition: bodies are split into contiguous blocks, one per process.
 * The first (N mod P) processes receive one extra body, so block sizes
 * differ by at most one.
 *
 * Each time step:
 *   1. Each process computes forces on its own bodies, using all positions.
 *   2. Each process updates the velocities and positions of its own bodies.
 *   3. MPI_Allgatherv shares the updated positions with every process.
 * Masses never change, so every process generates them locally and they
 * are never communicated.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define G 1.0            /* gravitational constant (simulation units) */
#define SOFTENING 1e-3   /* avoids division by zero when two bodies get very close */
#define SEED 42ULL       /* same seed as nbody_seq.c, so initial conditions match */

/* Deterministic pseudo-random generator (LCG). Returns a value in [0, 1). */
static unsigned long long rng_state;

static double rng_uniform(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / 9007199254740992.0;
}

/*
 * Fill positions, velocities, and masses for all N bodies.
 * Uses the same generator and the same draw order as nbody_seq.c,
 * so every process starts from identical initial conditions.
 * Layout: pos and vel hold x,y,z for body i at indices 3*i, 3*i+1, 3*i+2.
 */
static void init_state(double *pos, double *vel, double *mass, int n) {
    rng_state = SEED;
    for (int i = 0; i < n; i++) {
        pos[3 * i + 0] = 2.0 * rng_uniform() - 1.0;
        pos[3 * i + 1] = 2.0 * rng_uniform() - 1.0;
        pos[3 * i + 2] = 2.0 * rng_uniform() - 1.0;
        vel[3 * i + 0] = 0.1 * (2.0 * rng_uniform() - 1.0);
        vel[3 * i + 1] = 0.1 * (2.0 * rng_uniform() - 1.0);
        vel[3 * i + 2] = 0.1 * (2.0 * rng_uniform() - 1.0);
        mass[i] = 1.0 / n;
    }
}

/*
 * Accelerations for the bodies this process owns (indices start .. start+count-1).
 * Every body interacts with every other body, using the positions of all N bodies.
 * The results go into acc, which holds x,y,z for each local body.
 */
static void compute_local_accelerations(const double *pos, const double *mass,
                                        int n, int start, int count, double *acc) {
    for (int li = 0; li < count; li++) {
        int i = start + li;
        double axi = 0.0, ayi = 0.0, azi = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double dx = pos[3 * j + 0] - pos[3 * i + 0];
            double dy = pos[3 * j + 1] - pos[3 * i + 1];
            double dz = pos[3 * j + 2] - pos[3 * i + 2];
            double r2 = dx * dx + dy * dy + dz * dz + SOFTENING * SOFTENING;
            double inv_r3 = 1.0 / (r2 * sqrt(r2));
            double f = G * mass[j] * inv_r3;
            axi += f * dx;
            ayi += f * dy;
            azi += f * dz;
        }
        acc[3 * li + 0] = axi;
        acc[3 * li + 1] = ayi;
        acc[3 * li + 2] = azi;
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Argument checking: every process runs this, so all of them agree */
    if (argc < 4) {
        if (rank == 0)
            fprintf(stderr, "Usage: mpirun -np P %s N STEPS DT [output.csv]\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int n = atoi(argv[1]);
    int steps = atoi(argv[2]);
    double dt = atof(argv[3]);

    if (n <= 0 || steps < 0 || dt <= 0.0 || nprocs > n) {
        if (rank == 0)
            fprintf(stderr, "Error: need N > 0, STEPS >= 0, DT > 0, and processes <= N\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Data decomposition: contiguous blocks, first (n mod nprocs) ranks get one extra */
    int *counts = malloc((size_t)nprocs * sizeof(int));
    int *displs = malloc((size_t)nprocs * sizeof(int));
    if (!counts || !displs) {
        fprintf(stderr, "Rank %d: memory allocation failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int base = n / nprocs;
    int rem = n % nprocs;
    int offset = 0;
    for (int r = 0; r < nprocs; r++) {
        counts[r] = base + (r < rem ? 1 : 0);
        displs[r] = offset;
        offset += counts[r];
    }
    int start = displs[rank];
    int local_n = counts[rank];

    /* Counts and displacements in doubles for the position exchange (3 values per body) */
    int *counts3 = malloc((size_t)nprocs * sizeof(int));
    int *displs3 = malloc((size_t)nprocs * sizeof(int));
    /* Counts and displacements for the final gather (6 values per body) */
    int *counts6 = malloc((size_t)nprocs * sizeof(int));
    int *displs6 = malloc((size_t)nprocs * sizeof(int));
    for (int r = 0; r < nprocs; r++) {
        counts3[r] = 3 * counts[r];
        displs3[r] = 3 * displs[r];
        counts6[r] = 6 * counts[r];
        displs6[r] = 6 * displs[r];
    }

    double *pos = malloc((size_t)3 * n * sizeof(double));
    double *vel = malloc((size_t)3 * n * sizeof(double));
    double *mass = malloc((size_t)n * sizeof(double));
    double *acc = malloc((size_t)3 * local_n * sizeof(double));
    if (!pos || !vel || !mass || !acc) {
        fprintf(stderr, "Rank %d: memory allocation failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    init_state(pos, vel, mass, n);

    /* Timing: total, computation, and communication are measured separately */
    double t_comp = 0.0, t_comm = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    for (int s = 0; s < steps; s++) {
        double t0 = MPI_Wtime();

        /* Forces on local bodies, then update local velocities and positions */
        compute_local_accelerations(pos, mass, n, start, local_n, acc);
        for (int li = 0; li < local_n; li++) {
            int i = start + li;
            vel[3 * i + 0] += acc[3 * li + 0] * dt;
            vel[3 * i + 1] += acc[3 * li + 1] * dt;
            vel[3 * i + 2] += acc[3 * li + 2] * dt;
            pos[3 * i + 0] += vel[3 * i + 0] * dt;
            pos[3 * i + 1] += vel[3 * i + 1] * dt;
            pos[3 * i + 2] += vel[3 * i + 2] * dt;
        }

        double t1 = MPI_Wtime();

        /* Share updated positions: each process contributes its own block */
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       pos, counts3, displs3, MPI_DOUBLE, MPI_COMM_WORLD);

        double t2 = MPI_Wtime();
        t_comp += t1 - t0;
        t_comm += t2 - t1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t_total = MPI_Wtime() - t_start;

    /* Reduce timings across processes: slowest total, slowest and fastest compute, slowest comm */
    double max_total, max_comp, min_comp, max_comm;
    MPI_Reduce(&t_total, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_comp, &max_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_comp, &min_comp, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_comm, &max_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("P=%d N=%d steps=%d dt=%g total=%.6f s comp_max=%.6f s comp_min=%.6f s comm_max=%.6f s\n",
               nprocs, n, steps, dt, max_total, max_comp, min_comp, max_comm);
    }

    /* Optional output: gather final positions and velocities on rank 0 */
    if (argc >= 5) {
        double *local_buf = malloc((size_t)6 * local_n * sizeof(double));
        double *all = NULL;
        if (rank == 0) all = malloc((size_t)6 * n * sizeof(double));
        if (!local_buf || (rank == 0 && !all)) {
            fprintf(stderr, "Rank %d: memory allocation failed\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int li = 0; li < local_n; li++) {
            int i = start + li;
            local_buf[6 * li + 0] = pos[3 * i + 0];
            local_buf[6 * li + 1] = pos[3 * i + 1];
            local_buf[6 * li + 2] = pos[3 * i + 2];
            local_buf[6 * li + 3] = vel[3 * i + 0];
            local_buf[6 * li + 4] = vel[3 * i + 1];
            local_buf[6 * li + 5] = vel[3 * i + 2];
        }

        MPI_Gatherv(local_buf, 6 * local_n, MPI_DOUBLE,
                    all, counts6, displs6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            FILE *fp = fopen(argv[4], "w");
            if (!fp) {
                perror(argv[4]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            fprintf(fp, "id,x,y,z,vx,vy,vz,m\n");
            for (int i = 0; i < n; i++) {
                const double *b = &all[6 * i];
                fprintf(fp, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                        i, b[0], b[1], b[2], b[3], b[4], b[5], mass[i]);
            }
            fclose(fp);
            printf("Final state written to %s\n", argv[4]);
        }

        free(local_buf);
        free(all);
    }

    free(counts); free(displs); free(counts3); free(displs3);
    free(counts6); free(displs6);
    free(pos); free(vel); free(mass); free(acc);

    MPI_Finalize();
    return 0;
}
