#!/usr/bin/env python3
"""
CROSS-TERM ANISOTROPY SCAN

Does the cross-term formula R(N) = sqrt((N-2)/(N*4^(N-1))) hold beyond
Heisenberg XXX? Test:
  - XXZ with anisotropy Delta in [0, 2] (Delta=0: XY, Delta=1: XXX)
  - Pure Ising (ZZ only)
  - Pure XY (XX+YY only)
  - DM interaction (XY - YX, cross-Pauli coupling)

Analytical prediction: the bond-sum rule w_XY(a)+w_XY(b)=2 holds for
any coupling of the form J_X*XX + J_Y*YY + J_Z*ZZ (diagonal in Pauli
type). It breaks for cross-Pauli terms like X_i*Y_j.

Output: simulations/results/cross_term_anisotropy/
"""

import numpy as np
import json
import os
from datetime import datetime

REPO = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
OUT_DIR = os.path.join(REPO, "simulations", "results", "cross_term_anisotropy")
LOG_PATH = os.path.join(OUT_DIR, "cross_term_anisotropy.txt")
JSON_PATH = os.path.join(OUT_DIR, "cross_term_anisotropy.json")
os.makedirs(OUT_DIR, exist_ok=True)

_log = open(LOG_PATH, "w", encoding="utf-8", buffering=1)

def log(msg=""):
    print(msg, flush=True)
    _log.write(msg + "\n")


# ============================================================
# PAULI HELPERS
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


# ============================================================
# GENERAL HAMILTONIAN BUILDER
# ============================================================
def build_hamiltonian_general(N, bonds, terms):
    """
    Build H from a list of (coefficient, op_i, op_j) per bond.
    terms: list of (coeff, pauli_i, pauli_j) tuples.
    Example: Heisenberg = [(1, SX, SX), (1, SY, SY), (1, SZ, SZ)]
    """
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
                        - 0.5 * (np.kron(LkdLk, Id)
                                 + np.kron(Id, LkdLk.T)))
    return L_D


def compute_rel_ortho(L_H, L_D, N, gamma):
    d2 = L_H.shape[0]
    Sg = N * gamma
    L_Dc = L_D + Sg * np.eye(d2, dtype=complex)
    anticomm = L_H @ L_Dc + L_Dc @ L_H
    norm_cross = np.linalg.norm(anticomm)
    norm_LH = np.linalg.norm(L_H)
    norm_LDc = np.linalg.norm(L_Dc)
    denom = norm_LH * norm_LDc
    return norm_cross / denom if denom > 0 else 0.0


def formula_B(N):
    if N < 2:
        return 0.0
    return np.sqrt((N - 2) / (N * 4.0 ** (N - 1)))


# ============================================================
# COUPLING TYPES
# ============================================================
def heisenberg_terms(J=1.0):
    return [(J, SX, SX), (J, SY, SY), (J, SZ, SZ)]

def xxz_terms(J=1.0, delta=0.5):
    return [(J, SX, SX), (J, SY, SY), (J * delta, SZ, SZ)]

def xy_terms(J=1.0):
    return [(J, SX, SX), (J, SY, SY)]

def ising_terms(J=1.0):
    return [(J, SZ, SZ)]

def dm_terms(J=1.0):
    """Dzyaloshinskii-Moriya: X_i*Y_j - Y_i*X_j (cross-Pauli)."""
    return [(J, SX, SY), (-J, SY, SX)]

def xx_only(J=1.0):
    return [(J, SX, SX)]

def yy_only(J=1.0):
    return [(J, SY, SY)]

def zz_only(J=1.0):
    return [(J, SZ, SZ)]


# ============================================================
# MAIN
# ============================================================
GAMMA = 0.05
N_VALUES = [3, 4]

log("=" * 90)
log("CROSS-TERM ANISOTROPY SCAN")
log(f"Date: {datetime.now()}")
log(f"gamma = {GAMMA}")
log("=" * 90)

# ============================================================
# PART 1: Diagonal couplings (should all match formula)
# ============================================================
log()
log("=" * 90)
log("PART 1: DIAGONAL COUPLINGS (J_X*XX + J_Y*YY + J_Z*ZZ)")
log("Prediction: all match R(N) = sqrt((N-2)/(N*4^(N-1)))")
log("=" * 90)

diagonal_models = [
    ("Heisenberg (1,1,1)", heisenberg_terms()),
    ("XXZ Delta=0.5", xxz_terms(delta=0.5)),
    ("XXZ Delta=0.25", xxz_terms(delta=0.25)),
    ("XXZ Delta=1.5", xxz_terms(delta=1.5)),
    ("XXZ Delta=2.0", xxz_terms(delta=2.0)),
    ("XY model (1,1,0)", xy_terms()),
    ("Ising (0,0,1)", ising_terms()),
    ("XX only (1,0,0)", xx_only()),
    ("YY only (0,1,0)", yy_only()),
    ("ZZ only (0,0,1)", zz_only()),
]

results = []

for N in N_VALUES:
    bonds = chain_bonds(N)
    L_D = build_LD(N, GAMMA)
    predicted = formula_B(N)

    log(f"\n  N = {N}, predicted R = {predicted:.10f}")
    log(f"  {'Model':>25}  {'Measured':>14}  {'Deviation':>12}  {'Match':>6}")
    log(f"  {'-' * 65}")

    for name, terms in diagonal_models:
        H = build_hamiltonian_general(N, bonds, terms)
        L_H = build_LH(H)
        # Skip if H = 0
        if np.linalg.norm(L_H) < 1e-15:
            log(f"  {name:>25}  {'(L_H = 0)':>14}  {'---':>12}  {'SKIP':>6}")
            continue
        rel = compute_rel_ortho(L_H, L_D, N, GAMMA)
        dev = abs(rel - predicted)
        match = "YES" if dev < 1e-8 else "NO"
        log(f"  {name:>25}  {rel:>14.10f}  {dev:>12.2e}  {match:>6}")
        results.append({"N": N, "model": name, "type": "diagonal",
                         "rel_ortho": float(rel), "predicted": float(predicted),
                         "deviation": float(dev), "match": dev < 1e-8})

# ============================================================
# PART 2: Cross-Pauli couplings (should NOT match)
# ============================================================
log()
log("=" * 90)
log("PART 2: CROSS-PAULI COUPLINGS (X_i*Y_j type)")
log("Prediction: bond-sum rule breaks, formula does NOT hold")
log("=" * 90)

cross_models = [
    ("DM (XY - YX)", dm_terms()),
    ("XY cross only", [(1, SX, SY)]),
    ("XZ cross only", [(1, SX, SZ)]),
    ("YZ cross only", [(1, SY, SZ)]),
]

for N in N_VALUES:
    bonds = chain_bonds(N)
    L_D = build_LD(N, GAMMA)
    predicted = formula_B(N)

    log(f"\n  N = {N}, formula predicts R = {predicted:.10f}")
    log(f"  {'Model':>25}  {'Measured':>14}  {'Deviation':>12}  {'Match':>6}")
    log(f"  {'-' * 65}")

    for name, terms in cross_models:
        H = build_hamiltonian_general(N, bonds, terms)
        L_H = build_LH(H)
        if np.linalg.norm(L_H) < 1e-15:
            log(f"  {name:>25}  {'(L_H = 0)':>14}  {'---':>12}  {'SKIP':>6}")
            continue
        rel = compute_rel_ortho(L_H, L_D, N, GAMMA)
        dev = abs(rel - predicted)
        match = "YES" if dev < 1e-8 else "NO"
        log(f"  {name:>25}  {rel:>14.10f}  {dev:>12.2e}  {match:>6}")
        results.append({"N": N, "model": name, "type": "cross-Pauli",
                         "rel_ortho": float(rel), "predicted": float(predicted),
                         "deviation": float(dev), "match": dev < 1e-8})

# ============================================================
# PART 3: XXZ Delta sweep (continuous)
# ============================================================
log()
log("=" * 90)
log("PART 3: XXZ DELTA SWEEP (N=3)")
log("=" * 90)

deltas = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 5.0, 10.0]
N = 3
bonds = chain_bonds(N)
L_D = build_LD(N, GAMMA)
predicted = formula_B(N)

log(f"\n  N = {N}, predicted R = {predicted:.10f}")
log(f"  {'Delta':>10}  {'Measured':>14}  {'Deviation':>12}")
log(f"  {'-' * 40}")

delta_sweep = []
for delta in deltas:
    terms = xxz_terms(delta=delta)
    H = build_hamiltonian_general(N, bonds, terms)
    L_H = build_LH(H)
    if np.linalg.norm(L_H) < 1e-15:
        log(f"  {delta:>10.2f}  {'(L_H = 0)':>14}  {'---':>12}")
        continue
    rel = compute_rel_ortho(L_H, L_D, N, GAMMA)
    dev = abs(rel - predicted)
    log(f"  {delta:>10.2f}  {rel:>14.10f}  {dev:>12.2e}")
    delta_sweep.append({"delta": delta, "rel_ortho": float(rel),
                         "deviation": float(dev)})

# ============================================================
# VERDICT
# ============================================================
log()
log("=" * 90)
log("VERDICT")
log("=" * 90)

diag_all_match = all(r["match"] for r in results if r["type"] == "diagonal")
cross_any_match = any(r["match"] for r in results if r["type"] == "cross-Pauli")

log(f"\n  Diagonal couplings (XX, YY, ZZ, any coefficients):")
log(f"    All match formula: {'YES' if diag_all_match else 'NO'}")

log(f"\n  Cross-Pauli couplings (XY, XZ, YZ, DM):")
log(f"    Any match formula: {'YES' if cross_any_match else 'NO'}")

if diag_all_match and not cross_any_match:
    log(f"\n  RESULT: The formula R(N) = sqrt((N-2)/(N*4^(N-1))) holds for")
    log(f"  ALL diagonal Pauli couplings J_X*XX + J_Y*YY + J_Z*ZZ, including")
    log(f"  XXZ, XY model, Ising, and pure single-axis. It breaks for")
    log(f"  cross-Pauli couplings (DM interaction, X_iY_j terms).")
    log(f"\n  The universality class is: 'diagonal Pauli bond couplings'.")

# ============================================================
# SAVE
# ============================================================
output = {
    "metadata": {"date": str(datetime.now()), "gamma": GAMMA},
    "results": results,
    "delta_sweep": delta_sweep,
    "diagonal_all_match": bool(diag_all_match),
    "cross_any_match": bool(cross_any_match),
}
with open(JSON_PATH, "w", encoding="utf-8") as jf:
    json.dump(output, jf, indent=2)
log(f"\nSaved: {JSON_PATH}")

log()
log("=" * 90)
log("DONE")
log("=" * 90)
_log.close()
