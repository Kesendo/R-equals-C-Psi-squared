#!/usr/bin/env python3
"""
CROSS-TERM SHADOW-CROSSING FORMULA (EQ-012)

Conjecture: for shadow-crossing couplings (one Pauli in {X,Y}, one in
{I,Z}), the cross-term formula becomes:

    R^2(N, crossing) = (N-1) / (N * 4^(N-1))

replacing N-2 -> N-1 from the shadow-balanced case, because the bond
sites contribute a variance of 1 (instead of 0).

Tests: XZ, YZ, ZX, ZY couplings at N=3,4,5,6. Topology check at N=4.
Gamma independence check. Bond-site variance enumeration.

Outputs: simulations/results/cross_term_crossing/

Task: ClaudeTasks/TASK_CROSS_TERM_CROSSING.md
Closes: EQ-012
"""

import numpy as np
import json
import os
from datetime import datetime

REPO = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
OUT_DIR = os.path.join(REPO, "simulations", "results", "cross_term_crossing")
LOG_PATH = os.path.join(OUT_DIR, "cross_term_crossing.txt")
JSON_PATH = os.path.join(OUT_DIR, "cross_term_crossing.json")
os.makedirs(OUT_DIR, exist_ok=True)

_log = open(LOG_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _log.write(msg + "\n")

# ============================================================
# PAULI / LIOUVILLIAN HELPERS
# ============================================================
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

def site_op(op, k, N):
    ops = [I2] * N
    ops[k] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]

def complete_bonds(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]

def build_H(N, bonds, terms):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for coeff, Pi, Pj in terms:
            H += coeff * site_op(Pi, i, N) @ site_op(Pj, j, N)
    return H

def build_LH(H):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))

def build_LD(N, gamma):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_D = np.zeros((d * d, d * d), dtype=complex)
    for k in range(N):
        Lk = site_op(SZ, k, N)
        LkdLk = Lk.conj().T @ Lk
        L_D += gamma * (np.kron(Lk, Lk.conj())
                        - 0.5 * (np.kron(LkdLk, Id) + np.kron(Id, LkdLk.T)))
    return L_D

def compute_rel_ortho(L_H, L_D, N, gamma):
    d2 = L_H.shape[0]
    L_Dc = L_D + N * gamma * np.eye(d2, dtype=complex)
    anti = L_H @ L_Dc + L_Dc @ L_H
    nA = np.linalg.norm(anti)
    nH = np.linalg.norm(L_H)
    nD = np.linalg.norm(L_Dc)
    return nA / (nH * nD) if (nH * nD) > 0 else 0.0

def formula_balanced(N):
    return np.sqrt((N - 2) / (N * 4.0 ** (N - 1))) if N >= 2 else 0.0

def formula_crossing(N):
    return np.sqrt((N - 1) / (N * 4.0 ** (N - 1))) if N >= 2 else 0.0

# ============================================================
# COUPLING TYPES
# ============================================================
CROSSING_COUPLINGS = {
    "XZ": [(1, SX, SZ)],
    "YZ": [(1, SY, SZ)],
    "ZX": [(1, SZ, SX)],
    "ZY": [(1, SZ, SY)],
    "XZ+YZ": [(1, SX, SZ), (1, SY, SZ)],
}

GAMMA = 0.05

# ============================================================
# MAIN
# ============================================================
log("=" * 70)
log("CROSS-TERM SHADOW-CROSSING FORMULA (EQ-012)")
log(f"Date: {datetime.now()}")
log("=" * 70)

# ============================================================
# PART 1: Verify conjecture at N=3,4,5,6
# ============================================================
log()
log("=" * 70)
log("PART 1: CONJECTURE TEST")
log(f"  Conjecture: R^2 = (N-1) / (N * 4^(N-1))")
log("=" * 70)

results = []

for N in [3, 4, 5, 6]:
    bonds = chain_bonds(N)
    L_D = build_LD(N, GAMMA)
    pred_bal = formula_balanced(N)
    pred_cross = formula_crossing(N)

    log(f"\n  N = {N}, dim = {4**N}")
    log(f"  Predicted balanced:  {pred_bal:.10f}")
    log(f"  Predicted crossing:  {pred_cross:.10f}")
    log(f"  {'Coupling':>8}  {'Measured':>14}  {'Pred cross':>14}  "
        f"{'Dev':>10}  {'Match':>6}")
    log(f"  {'-' * 60}")

    for name, terms in CROSSING_COUPLINGS.items():
        H = build_H(N, bonds, terms)
        L_H = build_LH(H)
        if np.linalg.norm(L_H) < 1e-15:
            log(f"  {name:>8}  {'(L_H=0)':>14}")
            continue
        rel = compute_rel_ortho(L_H, L_D, N, GAMMA)
        dev = abs(rel - pred_cross)
        match = "YES" if dev < 1e-8 else "NO"
        log(f"  {name:>8}  {rel:>14.10f}  {pred_cross:>14.10f}  "
            f"{dev:>10.2e}  {match:>6}")
        results.append({"N": N, "coupling": name, "topology": "chain",
                         "rel_ortho": float(rel), "predicted": float(pred_cross),
                         "deviation": float(dev), "match": dev < 1e-8})

# ============================================================
# PART 2: Topology independence at N=4
# ============================================================
log()
log("=" * 70)
log("PART 2: TOPOLOGY INDEPENDENCE (N=4, chain vs complete)")
log("=" * 70)

N = 4
L_D = build_LD(N, GAMMA)
pred = formula_crossing(N)

for topo_name, bond_fn in [("chain", chain_bonds), ("complete", complete_bonds)]:
    bonds = bond_fn(N)
    terms = [(1, SX, SZ)]
    H = build_H(N, bonds, terms)
    L_H = build_LH(H)
    rel = compute_rel_ortho(L_H, L_D, N, GAMMA)
    dev = abs(rel - pred)
    log(f"  {topo_name:>10} ({len(bonds)} bonds): {rel:.10f}, "
        f"dev = {dev:.2e}")

# ============================================================
# PART 3: Gamma independence
# ============================================================
log()
log("=" * 70)
log("PART 3: GAMMA INDEPENDENCE (N=3, XZ coupling)")
log("=" * 70)

N = 3
bonds = chain_bonds(N)
H = build_H(N, bonds, [(1, SX, SZ)])
L_H = build_LH(H)

for gamma in [0.001, 0.01, 0.05, 0.1, 0.5]:
    L_D = build_LD(N, gamma)
    rel = compute_rel_ortho(L_H, L_D, N, gamma)
    log(f"  gamma = {gamma:.3f}: R = {rel:.10f}")

# ============================================================
# PART 4: Bond-site variance enumeration
# ============================================================
log()
log("=" * 70)
log("PART 4: BOND-SITE VARIANCE (XZ coupling, 2-site Pauli basis)")
log("=" * 70)

paulis = [I2, SX, SY, SZ]
pauli_names = ["I", "X", "Y", "Z"]
w_xy = [0, 1, 1, 0]

H2 = np.kron(SX, SZ)  # XZ coupling at 2 sites
d = 4

log(f"  {'Source':>6}  {'Target':>6}  {'w_src':>5}  {'w_tgt':>5}  "
    f"{'sum':>5}  {'s=sum-2':>7}  {'|M|^2':>8}")
log(f"  {'-' * 55}")

s_values = []
M_sq_values = []

for a_idx in range(16):
    a0, a1 = a_idx // 4, a_idx % 4
    Pa = np.kron(paulis[a0], paulis[a1])
    comm = H2 @ Pa - Pa @ H2
    if np.linalg.norm(comm) < 1e-12:
        continue
    for b_idx in range(16):
        b0, b1 = b_idx // 4, b_idx % 4
        Pb = np.kron(paulis[b0], paulis[b1])
        elem = -1j * np.trace(Pb.conj().T @ comm) / d
        if abs(elem) < 1e-12:
            continue
        w_src = w_xy[a0] + w_xy[a1]
        w_tgt = w_xy[b0] + w_xy[b1]
        s = w_src + w_tgt - 2
        M_sq = abs(elem) ** 2
        s_values.append(s)
        M_sq_values.append(M_sq)
        log(f"  {pauli_names[a0]+pauli_names[a1]:>6}  "
            f"{pauli_names[b0]+pauli_names[b1]:>6}  "
            f"{w_src:>5}  {w_tgt:>5}  {w_src+w_tgt:>5}  {s:>7}  "
            f"{M_sq:>8.1f}")

s_arr = np.array(s_values)
M_arr = np.array(M_sq_values)
mean_s = np.average(s_arr, weights=M_arr)
mean_s2 = np.average(s_arr ** 2, weights=M_arr)
log(f"\n  Weighted <s>  = {mean_s:.6f} (should be 0)")
log(f"  Weighted <s^2> = {mean_s2:.6f} (should be 1)")
log(f"  Bond-site variance = {mean_s2 - mean_s**2:.6f}")

# ============================================================
# VERDICT
# ============================================================
log()
log("=" * 70)
log("VERDICT")
log("=" * 70)

all_match = all(r["match"] for r in results)
log(f"\n  All measurements match conjecture: {'YES' if all_match else 'NO'}")
log(f"  Max deviation: {max(r['deviation'] for r in results):.2e}")

if all_match:
    log(f"\n  CONFIRMED: R^2(N, crossing) = (N-1) / (N * 4^(N-1))")
    log(f"  The bond-site variance is 1 (from <s^2> = 1, <s> = 0).")
    log(f"  Total variance = spectator (N-2) + bond (1) = N-1.")
    log(f"  The shadow-crossing formula is the balanced formula with N-2 -> N-1.")

# ============================================================
# SAVE
# ============================================================
output = {
    "metadata": {"date": str(datetime.now()), "gamma": GAMMA},
    "results": results,
    "bond_variance": {"mean_s": float(mean_s), "mean_s2": float(mean_s2)},
    "all_match": True if all_match else False,
}
with open(JSON_PATH, "w", encoding="utf-8") as jf:
    json.dump(output, jf, indent=2)
log(f"\nSaved: {JSON_PATH}")
log("=" * 70)
_log.close()
