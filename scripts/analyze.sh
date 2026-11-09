#!/bin/bash
# Summarises results/experiments.log: median of the repeated runs,
# speedup relative to the sequential baseline, and parallel efficiency.
set -e

cd "$(dirname "$0")/.."

LOG=results/experiments.log
OUT=results/summary.txt

awk '
function getval(key, s) {
    if (match(s, key "=[0-9.eE+-]+"))
        return substr(s, RSTART + length(key) + 1, RLENGTH - length(key) - 1) + 0
    return -1
}
function med3(x, y, z,   mn, mx) {
    mn = x; mx = x
    if (y < mn) mn = y
    if (z < mn) mn = z
    if (y > mx) mx = y
    if (z > mx) mx = z
    return x + y + z - mn - mx
}
/^SEQ/ {
    n = getval("N", $0)
    k = "SEQ " n
    c[k]++
    t[k, c[k]] = getval("time", $0)
}
/^MPI/ {
    p = getval("P", $0)
    n = getval("N", $0)
    k = "MPI " p " " n
    c[k]++
    tot[k, c[k]] = getval("total", $0)
    cp[k, c[k]]  = getval("comp_max", $0)
    cm[k, c[k]]  = getval("comm_max", $0)
}
END {
    ns = split("100 500 1000 5000 10000", Ns, " ")
    ps = split("1 2 4 8", Ps, " ")
    printf "%-7s %-3s %-14s %-10s %-11s %-14s %-14s\n", "N", "P", "Median_total_s", "Speedup", "Efficiency", "Median_comp_s", "Median_comm_s"
    for (i = 1; i <= ns; i++) {
        n = Ns[i]
        sk = "SEQ " n
        seq = med3(t[sk, 1], t[sk, 2], t[sk, 3])
        printf "%-7s %-3s %-14.6f %-10s %-11s %-14s %-14s\n", n, "seq", seq, "1.000", "-", "-", "-"
        for (j = 1; j <= ps; j++) {
            p = Ps[j]
            k = "MPI " p " " n
            mt = med3(tot[k, 1], tot[k, 2], tot[k, 3])
            mc = med3(cp[k, 1], cp[k, 2], cp[k, 3])
            mm = med3(cm[k, 1], cm[k, 2], cm[k, 3])
            sp = seq / mt
            printf "%-7s %-3s %-14.6f %-10.3f %-11.3f %-14.6f %-14.6f\n", n, p, mt, sp, sp / p, mc, mm
        }
    }
}' "$LOG" | tee "$OUT"

echo "Summary saved to $OUT"
