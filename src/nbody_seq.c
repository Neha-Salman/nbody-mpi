/*
 * nbody_seq.c - Sequential direct O(N^2) N-body simulation (baseline).
 *
 * Usage: ./nbody_seq N STEPS DT [output.csv]
 *   N          number of bodies
 *   STEPS      number of time steps
 *   DT         size of each time step
 *   output.csv optional file that receives the final state of every body
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define G 1.0            /* gravitational constant (simulation units) */
#define SOFTENING 1e-3   /* avoids division by zero when two bodies get very close */
#define SEED 42ULL       /* fixed seed so every run uses identical initial conditions */

typedef struct {
    double x, y, z;      /* position */
    double vx, vy, vz;   /* velocity */
    double m;            /* mass */
} Body;

/* Deterministic pseudo-random generator (LCG). Returns a value in [0, 1). */
static unsigned long long rng_state;

static double rng_uniform(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / 9007199254740992.0;
}

/* Fill the bodies with reproducible random positions, velocities, and masses. */
static void init_bodies(Body *b, int n) {
    rng_state = SEED;
    for (int i = 0; i < n; i++) {
        b[i].x  = 2.0 * rng_uniform() - 1.0;
        b[i].y  = 2.0 * rng_uniform() - 1.0;
        b[i].z  = 2.0 * rng_uniform() - 1.0;
        b[i].vx = 0.1 * (2.0 * rng_uniform() - 1.0);
        b[i].vy = 0.1 * (2.0 * rng_uniform() - 1.0);
        b[i].vz = 0.1 * (2.0 * rng_uniform() - 1.0);
        b[i].m  = 1.0 / n;
    }
}

/* Direct O(N^2) computation: acceleration of every body due to all others. */
static void compute_accelerations(const Body *b, int n,
                                  double *ax, double *ay, double *az) {
    for (int i = 0; i < n; i++) {
        double axi = 0.0, ayi = 0.0, azi = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double dx = b[j].x - b[i].x;
            double dy = b[j].y - b[i].y;
            double dz = b[j].z - b[i].z;
            double r2 = dx * dx + dy * dy + dz * dz + SOFTENING * SOFTENING;
            double inv_r3 = 1.0 / (r2 * sqrt(r2));
            double f = G * b[j].m * inv_r3;
            axi += f * dx;
            ayi += f * dy;
            azi += f * dz;
        }
        ax[i] = axi;
        ay[i] = ayi;
        az[i] = azi;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s N STEPS DT [output.csv]\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int steps = atoi(argv[2]);
    double dt = atof(argv[3]);

    if (n <= 0 || steps < 0 || dt <= 0.0) {
        fprintf(stderr, "Error: N must be > 0, STEPS >= 0, and DT > 0\n");
        return 1;
    }

    Body *bodies = malloc((size_t)n * sizeof(Body));
    double *ax = malloc((size_t)n * sizeof(double));
    double *ay = malloc((size_t)n * sizeof(double));
    double *az = malloc((size_t)n * sizeof(double));
    if (!bodies || !ax || !ay || !az) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(bodies); free(ax); free(ay); free(az);
        return 1;
    }

    init_bodies(bodies, n);

    /* Time only the simulation loop */
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int s = 0; s < steps; s++) {
        compute_accelerations(bodies, n, ax, ay, az);
        for (int i = 0; i < n; i++) {
            bodies[i].vx += ax[i] * dt;
            bodies[i].vy += ay[i] * dt;
            bodies[i].vz += az[i] * dt;
            bodies[i].x  += bodies[i].vx * dt;
            bodies[i].y  += bodies[i].vy * dt;
            bodies[i].z  += bodies[i].vz * dt;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    printf("N=%d steps=%d dt=%g time=%.6f s\n", n, steps, dt, elapsed);

    /* Optionally save the final state for later verification */
    if (argc >= 5) {
        FILE *fp = fopen(argv[4], "w");
        if (!fp) {
            perror(argv[4]);
            free(bodies); free(ax); free(ay); free(az);
            return 1;
        }
        fprintf(fp, "id,x,y,z,vx,vy,vz,m\n");
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    i, bodies[i].x, bodies[i].y, bodies[i].z,
                    bodies[i].vx, bodies[i].vy, bodies[i].vz, bodies[i].m);
        }
        fclose(fp);
        printf("Final state written to %s\n", argv[4]);
    }

    free(bodies); free(ax); free(ay); free(az);
    return 0;
}
