"""
Topology dependence of the Liouvillian degeneracy structure
============================================================
Compares d_real(k) and d_total(k) across Chain, Star, Ring, Complete
topologies for N=3..6 to determine whether the inner degeneracy
structure (k >= 2) depends on the graph.

Output: simulations/results/topology_degeneracy_comparison.txt
"""

import numpy as np
from pathlib import Path
import sys
import os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
GAMMA_GRID = 0.1  # 2 * gamma, Lindblad convention
TOL = 1e-8
TOPOLOGIES = ["chain", "star", "ring", "complete"]

output_lines = []
def log(msg=""):
    print(msg)
    output_lines.append(msg)


def load_eigenvalues(topo, n):
    if topo == "chain":
        path = RESULTS_DIR / f"rmt_eigenvalues_N{n}.csv"
    else:
        path = RESULTS_DIR / f"rmt_eigenvalues_{topo}_N{n}.csv"
    if not path.exists():
        return None
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    return data[:, 0] + 1j * data[:, 1]


def degeneracy_sequences(eigs, n):
    """Compute d_real(k) and d_total(k) for k=0..N."""
    d_real = []
    d_total = []
    for k in range(n + 1):
        target = -k * GAMMA_GRID
        # All eigenvalues at this grid position
        on_grid = np.abs(eigs.real - target) < TOL
        d_total.append(int(np.sum(on_grid)))
        # Purely real eigenvalues
        real_on_grid = on_grid & (np.abs(eigs.imag) < TOL)
        d_real.append(int(np.sum(real_on_grid)))
    return d_real, d_total


log("=" * 75)
log("TOPOLOGY DEPENDENCE OF DEGENERACY STRUCTURE")
log("=" * 75)
log(f"Topologies: {', '.join(TOPOLOGIES)}")
log(f"Grid spacing: {GAMMA_GRID} (= 2*gamma, gamma=0.05)")
log()

all_data = {}  # (topo, N) -> (d_real, d_total)
verification_pass = True

for N in range(3, 7):
    log("=" * 75)
    log(f"N = {N}")
    log("=" * 75)

    # Load and compute for each topology
    for topo in TOPOLOGIES:
        eigs = load_eigenvalues(topo, N)
        if eigs is None:
            log(f"  {topo:10s}: NO DATA")
            continue

        d_real, d_total = degeneracy_sequences(eigs, N)
        all_data[(topo, N)] = (d_real, d_total)

        is_palindrome_real = d_real == d_real[::-1]
        is_palindrome_total = d_total == d_total[::-1]

        log(f"  {topo:10s}: d_real = {d_real}  palindrome={'✓' if is_palindrome_real else '✗'}")
        log(f"  {'':10s}  d_total= {d_total}  palindrome={'✓' if is_palindrome_total else '✗'}")

        # Verification: k=0 must be N+1, k=1 must be >= 2N, palindrome
        if d_real[0] != N + 1:
            log(f"  *** BUG: d_real(0) = {d_real[0]} != {N+1} ***")
            verification_pass = False
        if d_real[1] < 2 * N:
            log(f"  *** BUG: d_real(1) = {d_real[1]} < {2*N} (proof guarantees >= 2N) ***")
            verification_pass = False
        elif d_real[1] > 2 * N:
            log(f"  NOTE: d_real(1) = {d_real[1]} > {2*N} (extra modes from multi-weight mixing)")
        if not is_palindrome_real:
            log(f"  *** BUG: d_real not palindromic ***")
            verification_pass = False

    # Comparison table
    log(f"\n  d_real comparison (N={N}):")
    header = f"  {'Topology':10s}"
    for k in range(N + 1):
        header += f"  k={k:d}"
    log(header)
    log("  " + "-" * (12 + 6 * (N + 1)))

    chain_dreal = all_data.get(("chain", N), (None, None))[0]

    for topo in TOPOLOGIES:
        if (topo, N) not in all_data:
            continue
        d_real, _ = all_data[(topo, N)]
        row = f"  {topo:10s}"
        for k in range(N + 1):
            marker = ""
            if chain_dreal is not None and d_real[k] != chain_dreal[k]:
                marker = "*"
            row += f"  {d_real[k]:3d}{marker}"
        log(row)

    # Check: are all d_real sequences identical?
    d_reals = [all_data[(t, N)][0] for t in TOPOLOGIES if (t, N) in all_data]
    all_same_real = all(d == d_reals[0] for d in d_reals)
    log(f"\n  d_real identical across topologies: {'✓ YES' if all_same_real else '✗ NO — TOPOLOGY-DEPENDENT'}")

    # Same for d_total
    log(f"\n  d_total comparison (N={N}):")
    header = f"  {'Topology':10s}"
    for k in range(N + 1):
        header += f"  k={k:d}"
    log(header)
    log("  " + "-" * (12 + 8 * (N + 1)))

    chain_dtotal = all_data.get(("chain", N), (None, None))[1]

    for topo in TOPOLOGIES:
        if (topo, N) not in all_data:
            continue
        _, d_total = all_data[(topo, N)]
        row = f"  {topo:10s}"
        for k in range(N + 1):
            marker = ""
            if chain_dtotal is not None and d_total[k] != chain_dtotal[k]:
                marker = "*"
            row += f"  {d_total[k]:5d}{marker}"
        log(row)

    d_totals = [all_data[(t, N)][1] for t in TOPOLOGIES if (t, N) in all_data]
    all_same_total = all(d == d_totals[0] for d in d_totals)
    log(f"\n  d_total identical across topologies: {'✓ YES' if all_same_total else '✗ NO — TOPOLOGY-DEPENDENT'}")

    log()

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────

log("=" * 75)
log("SUMMARY")
log("=" * 75)
log()

# Collect results per N
log("d_real(k) comparison across all N:")
log()

any_difference_real = False
any_difference_total = False

for N in range(3, 7):
    d_reals = {t: all_data[(t, N)][0] for t in TOPOLOGIES if (t, N) in all_data}
    d_totals = {t: all_data[(t, N)][1] for t in TOPOLOGIES if (t, N) in all_data}

    same_real = len(set(tuple(v) for v in d_reals.values())) == 1
    same_total = len(set(tuple(v) for v in d_totals.values())) == 1

    if not same_real:
        any_difference_real = True
        log(f"  N={N}: d_real DIFFERS across topologies")
        for t, d in d_reals.items():
            log(f"    {t:10s}: {d}")
    else:
        log(f"  N={N}: d_real identical = {list(d_reals.values())[0]}")

    if not same_total:
        any_difference_total = True

log()

if any_difference_real:
    log("RESULT A: d_real(k) is TOPOLOGY-DEPENDENT for k >= 2.")
    log("The SWAP-invariance proof covers only k=0 and k=1.")
    log("Inner degeneracy structure depends on the bond graph.")
else:
    log("RESULT B: d_real(k) is TOPOLOGY-INDEPENDENT for all k tested.")
    log("Strong evidence for a universal formula d_real(k, N).")

log()
log(f"d_total topology-dependent: {'YES' if any_difference_total else 'NO'}")
log(f"Verification (k=0, k=1, palindrome): {'PASS' if verification_pass else 'FAIL'}")

# Save
out_path = RESULTS_DIR / "topology_degeneracy_comparison.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(output_lines))
print(f"\n>>> Results saved to: {out_path}")
