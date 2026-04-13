#!/usr/bin/env python3
"""
CROSS-TERM TOPOLOGY DEPENDENCE

Scientific question:
    The anti-commutator {L_H, L_Dc} (where L_Dc = L_D + Sg*I) vanishes
    exactly at N=2 (Pythagorean decomposition) but has relative magnitude
    1/sqrt(48) at N=3 for the Heisenberg chain. Is this constant
    chain-specific, or universal across topologies? Is there a closed-form
    relation between graph topology and the cross-term magnitude?

Methodology:
    For N=3 and N=4, build the Heisenberg Hamiltonian on each of
    {chain, ring, star, complete} with uniform Z-dephasing. Compute:
      1. L_H (Hamiltonian superoperator)
      2. L_Dc = L_D + Sg*I (centered dissipator), Sg = N*gamma
      3. {L_H, L_Dc} (anti-commutator)
      4. ||{L_H, L_Dc}|| / (||L_H|| * ||L_Dc||) (relative orthogonality)
    Verify gamma-independence by sweeping gamma in {0.001, 0.01, 0.05, 0.1, 0.5}.

Outputs:
    simulations/results/cross_term_topology/
        cross_term_topology.txt      (run log)
        cross_term_topology.json     (structured data)
        cross_term_topology.png      (bar chart)

Task file: ClaudeTasks/TASK_CROSS_TERM_TOPOLOGY.md
Parent: experiments/PRIMORDIAL_QUBIT_ALGEBRA.md ("What to compute next" item 1)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json
import os
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================
N_VALUES = [3, 4]
GAMMA = 0.05
J = 1.0
TOPOLOGIES = ["chain", "star", "ring", "complete"]
GAMMA_SWEEP = [0.001, 0.01, 0.05, 0.1, 0.5]

REPO = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
OUT_DIR = os.path.join(REPO, "simulations", "results", "cross_term_topology")
LOG_PATH = os.path.join(OUT_DIR, "cross_term_topology.txt")
JSON_PATH = os.path.join(OUT_DIR, "cross_term_topology.json")
PLOT_PATH = os.path.join(OUT_DIR, "cross_term_topology.png")

os.makedirs(OUT_DIR, exist_ok=True)

_log_file = open(LOG_PATH, "w", encoding="utf-8", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    _log_file.write(msg + "\n")
    _log_file.flush()


# ============================================================
# PAULI HELPERS (self-contained, no imports from other scripts)
# ============================================================
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, k, N):
    """Pauli op on site k, identity on the others."""
    ops = [I2] * N
    ops[k] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result


# ============================================================
# TOPOLOGY GENERATORS
# ============================================================
def chain_bonds(N):
    """Linear chain: 0-1-2-..-(N-1). Returns N-1 bonds."""
    return [(i, i + 1) for i in range(N - 1)]


def star_bonds(N):
    """Star with center at site 0. Returns N-1 bonds."""
    return [(0, i) for i in range(1, N)]


def ring_bonds(N):
    """Ring: 0-1-2-..-(N-1)-0. Returns N bonds."""
    return [(i, (i + 1) % N) for i in range(N)]


def complete_bonds(N):
    """Complete graph: all pairs. Returns N*(N-1)/2 bonds."""
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


BOND_GENERATORS = {
    "chain": chain_bonds,
    "star": star_bonds,
    "ring": ring_bonds,
    "complete": complete_bonds,
}


def bonds_str(bonds):
    """Human-readable bond list."""
    return ", ".join(f"({a},{b})" for a, b in bonds)


def degree_sequence(N, bonds):
    """Sorted degree sequence of the graph."""
    deg = [0] * N
    for a, b in bonds:
        deg[a] += 1
        deg[b] += 1
    return sorted(deg)


# ============================================================
# HAMILTONIAN AND LIOUVILLIAN BUILDERS
# ============================================================
def build_hamiltonian(N, bonds, J):
    """Heisenberg XXX Hamiltonian from a bond list."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P in (SX, SY, SZ):
            H += J * site_op(P, i, N) @ site_op(P, j, N)
    return H


def build_LH(H):
    """L_H = -1j * (H kron I - I kron H^T)."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def build_LD(N, gamma):
    """Dissipator for uniform Z-dephasing at rate gamma per site."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_D = np.zeros((d * d, d * d), dtype=complex)
    for k in range(N):
        Lk = site_op(SZ, k, N)
        LkdLk = Lk.conj().T @ Lk
        L_D += gamma * (np.kron(Lk, Lk.conj())
                        - 0.5 * (np.kron(LkdLk, Id)
                                 + np.kron(Id, LkdLk.T)))
    return L_D


def compute_cross_term(L_H, L_D, N, gamma):
    """
    Compute ||{L_H, L_Dc}|| / (||L_H|| * ||L_Dc||)
    where L_Dc = L_D + Sg*I, Sg = N*gamma.
    Returns (norm_cross, rel_ortho, norm_LH, norm_LDc).
    """
    d2 = L_H.shape[0]
    Sg = N * gamma
    L_Dc = L_D + Sg * np.eye(d2, dtype=complex)
    anticomm = L_H @ L_Dc + L_Dc @ L_H
    norm_cross = np.linalg.norm(anticomm)
    norm_LH = np.linalg.norm(L_H)
    norm_LDc = np.linalg.norm(L_Dc)
    denom = norm_LH * norm_LDc
    rel = norm_cross / denom if denom > 0 else 0.0
    return norm_cross, rel, norm_LH, norm_LDc


# ============================================================
# CLOSED-FORM SEARCH
# ============================================================
def try_identify(val, tol=1e-6):
    """Try to identify val as 1/sqrt(k) for small integer k, or 0."""
    if val < 1e-12:
        return "0"
    k_candidate = 1.0 / (val ** 2)
    k_int = round(k_candidate)
    if k_int > 0 and abs(k_candidate - k_int) < 0.5:
        check = 1.0 / np.sqrt(k_int)
        if abs(val - check) < tol:
            return f"1/sqrt({k_int})"
    # Try p/sqrt(k) for small p
    for p in range(2, 6):
        k2 = p ** 2 / (val ** 2)
        k2_int = round(k2)
        if k2_int > 0 and abs(k2 - k2_int) < 0.5:
            check = p / np.sqrt(k2_int)
            if abs(val - check) < tol:
                return f"{p}/sqrt({k2_int})"
    return None


# ============================================================
# MAIN
# ============================================================
log("=" * 90)
log("CROSS-TERM TOPOLOGY DEPENDENCE")
log(f"Date: {datetime.now()}")
log(f"N = {N_VALUES}, gamma = {GAMMA}, J = {J}")
log(f"Topologies: {TOPOLOGIES}")
log("=" * 90)

# ============================================================
# ANCHOR: verify N=3 chain = 1/sqrt(48)
# ============================================================
log()
log("=" * 90)
log("ANCHOR: N=3 chain Z-dephasing, expected 1/sqrt(48) = 0.144338...")
log("=" * 90)

H_anchor = build_hamiltonian(3, chain_bonds(3), J)
L_H_anchor = build_LH(H_anchor)
L_D_anchor = build_LD(3, GAMMA)
_, rel_anchor, _, _ = compute_cross_term(L_H_anchor, L_D_anchor, 3, GAMMA)
expected_chain3 = 1.0 / np.sqrt(48)

log(f"  Measured: {rel_anchor:.10f}")
log(f"  Expected: {expected_chain3:.10f}")
log(f"  Deviation: {abs(rel_anchor - expected_chain3):.2e}")

if abs(rel_anchor - expected_chain3) < 1e-8:
    log("  ANCHOR: PASS")
else:
    log("  ANCHOR: FAIL")
    _log_file.close()
    raise SystemExit(1)

# ============================================================
# TOPOLOGY EQUIVALENCE CHECK
# ============================================================
log()
log("=" * 90)
log("TOPOLOGY EQUIVALENCE CHECK")
log("=" * 90)

for N in N_VALUES:
    log(f"\n  N = {N}:")
    for topo in TOPOLOGIES:
        bonds = BOND_GENERATORS[topo](N)
        degs = degree_sequence(N, bonds)
        log(f"    {topo:>10}: {len(bonds)} bonds, edges = [{bonds_str(bonds)}], "
            f"degrees = {degs}")

    if N == 3:
        ring_set = set(tuple(sorted(b)) for b in ring_bonds(3))
        comp_set = set(tuple(sorted(b)) for b in complete_bonds(3))
        if ring_set == comp_set:
            log("    NOTE: ring == complete at N=3 (triangle = K_3)")

        chain_deg = degree_sequence(3, chain_bonds(3))
        star_deg = degree_sequence(3, star_bonds(3))
        if chain_deg == star_deg:
            log("    NOTE: chain and star are isomorphic at N=3 "
                "(both are the path P_3, degree sequence [1, 1, 2])")

# ============================================================
# MAIN COMPUTATION: cross term per (N, topology)
# ============================================================
log()
log("=" * 90)
log(f"CROSS-TERM RELATIVE ORTHOGONALITY (gamma = {GAMMA})")
log("=" * 90)

results = []

for N in N_VALUES:
    L_D = build_LD(N, GAMMA)
    log(f"\n  N = {N}, superoperator dim = {4**N}")
    log(f"  {'Topology':>10}  {'Edges':>5}  {'||{L_H,L_Dc}||':>16}  "
        f"{'rel_ortho':>12}  {'||L_H||':>12}  {'||L_Dc||':>12}")
    log(f"  {'-' * 80}")

    for topo in TOPOLOGIES:
        bonds = BOND_GENERATORS[topo](N)
        H = build_hamiltonian(N, bonds, J)
        L_H = build_LH(H)
        norm_cross, rel, norm_LH, norm_LDc = compute_cross_term(
            L_H, L_D, N, GAMMA)

        results.append({
            "N": N,
            "topology": topo,
            "n_edges": len(bonds),
            "bonds": [(a, b) for a, b in bonds],
            "norm_cross": float(norm_cross),
            "rel_ortho": float(rel),
            "norm_LH": float(norm_LH),
            "norm_LDc": float(norm_LDc),
            "gamma": GAMMA,
        })

        log(f"  {topo:>10}  {len(bonds):>5}  {norm_cross:>16.8f}  "
            f"{rel:>12.8f}  {norm_LH:>12.4f}  {norm_LDc:>12.4f}")

# ============================================================
# GAMMA SWEEP: verify gamma-independence
# ============================================================
log()
log("=" * 90)
log("GAMMA SWEEP (verify gamma-independence)")
log(f"Values: {GAMMA_SWEEP}")
log("=" * 90)

gamma_sweeps = []
gamma_independent = True

for topo in TOPOLOGIES:
    for N in N_VALUES:
        bonds = BOND_GENERATORS[topo](N)
        H = build_hamiltonian(N, bonds, J)
        L_H = build_LH(H)

        log(f"\n  {topo} N={N}:")
        log(f"    {'gamma':>10}  {'rel_ortho':>14}")
        log(f"    {'-' * 28}")

        sweep_vals = []
        for gamma in GAMMA_SWEEP:
            L_D = build_LD(N, gamma)
            _, rel, _, _ = compute_cross_term(L_H, L_D, N, gamma)
            sweep_vals.append(rel)
            gamma_sweeps.append({
                "topology": topo,
                "N": N,
                "gamma": gamma,
                "rel_ortho": float(rel),
            })
            log(f"    {gamma:>10.3f}  {rel:>14.10f}")

        spread = max(sweep_vals) - min(sweep_vals)
        status = "CONSTANT" if spread < 1e-8 else "VARIES"
        if spread >= 1e-8:
            gamma_independent = False
        log(f"    spread = {spread:.2e} ({status})")

log(f"\n  Overall gamma-independence: "
    f"{'CONFIRMED' if gamma_independent else 'FAILED'}")

# ============================================================
# ANALYSIS: closed-form identification
# ============================================================
log()
log("=" * 90)
log("CLOSED-FORM ANALYSIS")
log("=" * 90)

for N in N_VALUES:
    log(f"\n  N = {N}:")
    n_results = [r for r in results if r["N"] == N]
    for r in n_results:
        val = r["rel_ortho"]
        form = try_identify(val)
        r["closed_form"] = form
        if form:
            log(f"    {r['topology']:>10}: {val:.10f} = {form}")
        else:
            log(f"    {r['topology']:>10}: {val:.10f} (no simple closed form)")

# Check secondary hypothesis
log()
log("SECONDARY HYPOTHESIS: complete graph N=3 vanishing cross term")
comp_n3 = [r for r in results if r["topology"] == "complete" and r["N"] == 3]
if comp_n3:
    val = comp_n3[0]["rel_ortho"]
    if val < 1e-12:
        log(f"  CONFIRMED: complete graph N=3 has zero cross term ({val:.2e})")
    else:
        log(f"  REFUTED: complete graph N=3 has nonzero cross term ({val:.10f})")

# Edge-count vs cross term
log()
log("EDGE-COUNT CORRELATION:")
log(f"  {'N':>3}  {'Topology':>10}  {'Edges':>5}  {'rel_ortho':>14}  {'Closed form':>16}")
log(f"  {'-' * 60}")
for r in results:
    form = r.get("closed_form", "?")
    log(f"  {r['N']:>3}  {r['topology']:>10}  {r['n_edges']:>5}  "
        f"{r['rel_ortho']:>14.10f}  {form if form else '---':>16}")

# Check if same edges -> same cross term
log()
log("GRAPH ISOMORPHISM CHECK:")
for N in N_VALUES:
    n_res = [r for r in results if r["N"] == N]
    for i, r1 in enumerate(n_res):
        for r2 in n_res[i + 1:]:
            d1 = degree_sequence(N, r1["bonds"])
            d2 = degree_sequence(N, r2["bonds"])
            same_degs = d1 == d2
            same_val = abs(r1["rel_ortho"] - r2["rel_ortho"]) < 1e-8
            if same_degs or same_val:
                log(f"  N={N}: {r1['topology']} vs {r2['topology']}: "
                    f"same_degrees={same_degs}, same_cross_term={same_val}")

# ============================================================
# SAVE JSON
# ============================================================
output_data = {
    "metadata": {
        "date": str(datetime.now()),
        "J": J,
        "gamma_default": GAMMA,
        "gamma_sweep": GAMMA_SWEEP,
        "N_values": N_VALUES,
        "topologies": TOPOLOGIES,
    },
    "results": [{k: v for k, v in r.items()} for r in results],
    "gamma_sweeps": gamma_sweeps,
}

with open(JSON_PATH, "w", encoding="utf-8") as jf:
    json.dump(output_data, jf, indent=2)
log(f"\nSaved JSON: {JSON_PATH}")

# ============================================================
# PLOT: bar chart grouped by N, colored by topology
# ============================================================
log()
log("=" * 90)
log("GENERATING PLOT")
log("=" * 90)

fig, ax = plt.subplots(figsize=(10, 6))

topo_colors = {
    "chain": "#1f77b4",
    "star": "#ff7f0e",
    "ring": "#2ca02c",
    "complete": "#d62728",
}

n_topos = len(TOPOLOGIES)
bar_width = 0.18
x_base = np.arange(len(N_VALUES))

for idx, topo in enumerate(TOPOLOGIES):
    vals = []
    for N in N_VALUES:
        match = [r for r in results if r["N"] == N and r["topology"] == topo]
        vals.append(match[0]["rel_ortho"] if match else 0)
    x_pos = x_base + (idx - n_topos / 2 + 0.5) * bar_width
    bars = ax.bar(x_pos, vals, bar_width, label=topo,
                  color=topo_colors[topo], edgecolor="black", linewidth=0.5)
    for bar, val in zip(bars, vals):
        if val > 0.001:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.003,
                    f"{val:.4f}", ha="center", va="bottom", fontsize=8)

ax.axhline(expected_chain3, color="red", linestyle="--", alpha=0.4,
           label=f"1/\u221a48 \u2248 {expected_chain3:.4f} (N=3 chain)")

ax.set_xlabel("System size N")
ax.set_ylabel(r"Relative orthogonality $||\{L_H, L_{D,c}\}||"
              r" / (||L_H||\cdot||L_{D,c}||)$")
ax.set_title("Cross-term topology dependence, Heisenberg XXX, "
             "uniform Z-dephasing")
ax.set_xticks(x_base)
ax.set_xticklabels([f"N={N}" for N in N_VALUES])
ax.legend(loc="upper left", fontsize=9)
ax.grid(True, axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(PLOT_PATH, dpi=150, bbox_inches="tight")
log(f"Saved plot: {PLOT_PATH}")
plt.close()

# ============================================================
# DONE
# ============================================================
log()
log("=" * 90)
log("RUN COMPLETE")
log(f"Log:  {LOG_PATH}")
log(f"JSON: {JSON_PATH}")
log(f"Plot: {PLOT_PATH}")
log("=" * 90)
_log_file.close()
