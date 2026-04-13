#!/usr/bin/env python3
"""
CROSS-TERM FORMULA CONFIRMATION AND PROOF (EQ-011)

Part 1: Numerical confirmation at N=5 (chain + complete) and N=6 (chain).
Part 2: Analytical proof via Pauli-basis decomposition.

Two candidate formulas:
  (A) Conjecture from CROSS_TERM_TOPOLOGY:  R = 1/sqrt(N * 2^(N+1))
  (B) Derived from single-bond proof:       R = sqrt((N-2) / (N * 4^(N-1)))

Both agree at N=3,4 but diverge at N>=5. This script decides which is correct.

Outputs:
    simulations/results/cross_term_formula/
        cross_term_formula.txt      (run log)
        cross_term_formula.json     (structured data)

Task file: ClaudeTasks/TASK_CROSS_TERM_FORMULA.md
Parent: experiments/CROSS_TERM_TOPOLOGY.md
"""

import numpy as np
import json
import os
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================
GAMMA = 0.05
J = 1.0

REPO = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
OUT_DIR = os.path.join(REPO, "simulations", "results", "cross_term_formula")
LOG_PATH = os.path.join(OUT_DIR, "cross_term_formula.txt")
JSON_PATH = os.path.join(OUT_DIR, "cross_term_formula.json")

os.makedirs(OUT_DIR, exist_ok=True)

_log_file = open(LOG_PATH, "w", encoding="utf-8", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    _log_file.write(msg + "\n")
    _log_file.flush()


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


def complete_bonds(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


def build_hamiltonian(N, bonds, J):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P in (SX, SY, SZ):
            H += J * site_op(P, i, N) @ site_op(P, j, N)
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


def compute_cross_term(L_H, L_D, N, gamma):
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
# FORMULA PREDICTIONS
# ============================================================
def formula_A(N):
    """Old conjecture: 1/sqrt(N * 2^(N+1))."""
    return 1.0 / np.sqrt(N * 2 ** (N + 1))


def formula_B(N):
    """Derived formula: sqrt((N-2) / (N * 4^(N-1)))."""
    if N < 2:
        return 0.0
    return np.sqrt((N - 2) / (N * 4.0 ** (N - 1)))


# ============================================================
# MAIN
# ============================================================
log("=" * 90)
log("CROSS-TERM FORMULA CONFIRMATION AND PROOF (EQ-011)")
log(f"Date: {datetime.now()}")
log(f"gamma = {GAMMA}, J = {J}")
log("=" * 90)

# ============================================================
# FORMULA COMPARISON TABLE
# ============================================================
log()
log("=" * 90)
log("FORMULA PREDICTIONS")
log("  (A) Old conjecture: R = 1/sqrt(N * 2^(N+1))")
log("  (B) Derived formula: R = sqrt((N-2) / (N * 4^(N-1)))")
log("=" * 90)
log()
log(f"  {'N':>3}  {'Formula A':>14}  {'Formula B':>14}  {'Match':>8}")
log(f"  {'-' * 45}")
for N in range(2, 8):
    a = formula_A(N)
    b = formula_B(N)
    match = "YES" if abs(a - b) < 1e-12 else "NO"
    log(f"  {N:>3}  {a:>14.10f}  {b:>14.10f}  {match:>8}")

log()
log("  N=2: Formula A gives 1/4, Formula B gives 0. Correct answer is 0.")
log("  N=3,4: Both formulas agree.")
log("  N>=5: Formulas diverge. This computation decides which is correct.")

# ============================================================
# ANCHOR: N=3 and N=4 (must match known values)
# ============================================================
log()
log("=" * 90)
log("ANCHOR: N=3 and N=4 (known values from CROSS_TERM_TOPOLOGY)")
log("=" * 90)

for N in [3, 4]:
    bonds = chain_bonds(N)
    H = build_hamiltonian(N, bonds, J)
    L_H = build_LH(H)
    L_D = build_LD(N, GAMMA)
    _, rel, _, _ = compute_cross_term(L_H, L_D, N, GAMMA)
    expected = formula_B(N)
    dev = abs(rel - expected)
    log(f"  N={N}: measured={rel:.10f}, predicted={expected:.10f}, "
        f"deviation={dev:.2e} ({'PASS' if dev < 1e-8 else 'FAIL'})")

# ============================================================
# PART 1: N=5 (chain + complete)
# ============================================================
log()
log("=" * 90)
log("PART 1a: N=5 (chain and complete)")
log(f"  Superoperator dim = {4**5} = 1024")
log(f"  Formula A predicts: {formula_A(5):.10f}")
log(f"  Formula B predicts: {formula_B(5):.10f}")
log("=" * 90)

results = []

for topo_name, bond_fn in [("chain", chain_bonds), ("complete", complete_bonds)]:
    bonds = bond_fn(5)
    H = build_hamiltonian(5, bonds, J)
    L_H = build_LH(H)
    L_D = build_LD(5, GAMMA)
    norm_cross, rel, norm_LH, norm_LDc = compute_cross_term(
        L_H, L_D, 5, GAMMA)

    dev_A = abs(rel - formula_A(5))
    dev_B = abs(rel - formula_B(5))

    log(f"\n  N=5 {topo_name} ({len(bonds)} bonds):")
    log(f"    measured rel_ortho   = {rel:.12f}")
    log(f"    Formula A prediction = {formula_A(5):.12f}  "
        f"deviation = {dev_A:.4e}")
    log(f"    Formula B prediction = {formula_B(5):.12f}  "
        f"deviation = {dev_B:.4e}")
    log(f"    Winner: {'A' if dev_A < dev_B else 'B'}")

    results.append({
        "N": 5, "topology": topo_name, "n_edges": len(bonds),
        "rel_ortho": float(rel),
        "norm_cross": float(norm_cross),
        "norm_LH": float(norm_LH),
        "norm_LDc": float(norm_LDc),
        "dev_formula_A": float(dev_A),
        "dev_formula_B": float(dev_B),
    })

# Topology independence at N=5
n5_vals = [r["rel_ortho"] for r in results if r["N"] == 5]
spread_5 = max(n5_vals) - min(n5_vals) if len(n5_vals) > 1 else 0
log(f"\n  Topology independence at N=5: spread = {spread_5:.2e} "
    f"({'CONFIRMED' if spread_5 < 1e-8 else 'BROKEN'})")

# ============================================================
# PART 1b: N=6 (chain only)
# ============================================================
log()
log("=" * 90)
log("PART 1b: N=6 (chain only)")
log(f"  Superoperator dim = {4**6} = 4096")
log(f"  Matrix memory: ~256 MB (dense complex128)")
log(f"  Formula A predicts: {formula_A(6):.10f}")
log(f"  Formula B predicts: {formula_B(6):.10f}")
log("=" * 90)

bonds_6 = chain_bonds(6)
H_6 = build_hamiltonian(6, bonds_6, J)
L_H_6 = build_LH(H_6)
L_D_6 = build_LD(6, GAMMA)
norm_cross_6, rel_6, norm_LH_6, norm_LDc_6 = compute_cross_term(
    L_H_6, L_D_6, 6, GAMMA)

dev_A_6 = abs(rel_6 - formula_A(6))
dev_B_6 = abs(rel_6 - formula_B(6))

log(f"\n  N=6 chain ({len(bonds_6)} bonds):")
log(f"    measured rel_ortho   = {rel_6:.12f}")
log(f"    Formula A prediction = {formula_A(6):.12f}  "
    f"deviation = {dev_A_6:.4e}")
log(f"    Formula B prediction = {formula_B(6):.12f}  "
    f"deviation = {dev_B_6:.4e}")
log(f"    Winner: {'A' if dev_A_6 < dev_B_6 else 'B'}")

results.append({
    "N": 6, "topology": "chain", "n_edges": len(bonds_6),
    "rel_ortho": float(rel_6),
    "norm_cross": float(norm_cross_6),
    "norm_LH": float(norm_LH_6),
    "norm_LDc": float(norm_LDc_6),
    "dev_formula_A": float(dev_A_6),
    "dev_formula_B": float(dev_B_6),
})

# Free N=6 memory
del H_6, L_H_6, L_D_6
log("  (N=6 matrices freed)")

# ============================================================
# FORMULA VERDICT
# ============================================================
log()
log("=" * 90)
log("FORMULA VERDICT")
log("=" * 90)

all_dev_B = [r["dev_formula_B"] for r in results]
all_dev_A = [r["dev_formula_A"] for r in results]

max_dev_B = max(all_dev_B)
max_dev_A = max(all_dev_A)

log(f"  Max deviation from Formula A: {max_dev_A:.4e}")
log(f"  Max deviation from Formula B: {max_dev_B:.4e}")

if max_dev_B < 1e-8 and max_dev_A > 1e-4:
    verdict = "FORMULA B CONFIRMED, FORMULA A REFUTED"
    log(f"\n  VERDICT: {verdict}")
    log("  The correct formula is:")
    log("    R(N) = sqrt((N-2) / (N * 4^(N-1)))  for all N >= 2")
    log("  The old conjecture 1/sqrt(N * 2^(N+1)) is refuted at N=5.")
elif max_dev_A < 1e-8 and max_dev_B > 1e-4:
    verdict = "FORMULA A CONFIRMED, FORMULA B REFUTED"
    log(f"\n  VERDICT: {verdict}")
elif max_dev_A < 1e-8 and max_dev_B < 1e-8:
    verdict = "BOTH FORMULAS CONFIRMED (unexpected)"
    log(f"\n  VERDICT: {verdict}")
else:
    verdict = "NEITHER FORMULA MATCHES"
    log(f"\n  VERDICT: {verdict}")

# ============================================================
# PART 2: ANALYTICAL PROOF VERIFICATION
# ============================================================
log()
log("=" * 90)
log("PART 2: ANALYTICAL PROOF VERIFICATION")
log("=" * 90)

log("""
  The proof proceeds in four steps:

  Step 1: ||L_Dc||^2 = gamma^2 * 4^N * N
    L_Dc is diagonal in the Pauli basis with eigenvalue gamma*(N - 2*w_XY).
    Sum of squares: Sum_a (N - 2*w_a)^2 = 4^N * N.
    (Each site contributes independently; cross-site terms vanish because
    Sum_{4 Paulis} (1 - 2*delta_XY) = 0.)

  Step 2: Every Heisenberg bond transition satisfies w_XY(a) + w_XY(b) = 2
    at the bond sites (the "bond-sum rule").
    This is the same property that makes the Pythagorean decomposition
    exact at N=2.

  Step 3: The spectator variance is N-2
    For any bond (i,j), the non-bond sites are spectators: unchanged by
    the transition, contributing (N-2-2*w_rest)^2 to the anti-commutator.
    The average over all 4^(N-2) spectator configs is:
      <(N-2-2*w_rest)^2> = N-2.

  Step 4: Assembly
    ||{L_H, L_Dc}||^2 = 4*gamma^2 * (N-2) * ||L_H||^2
    R^2 = 4*gamma^2*(N-2)*||L_H||^2 / (||L_H||^2 * gamma^2 * 4^N * N)
        = 4*(N-2) / (N * 4^N)
        = (N-2) / (N * 4^(N-1))
""")

# Verify Step 1: ||L_Dc||^2 = gamma^2 * 4^N * N
log("  Step 1 verification: ||L_Dc||^2 = gamma^2 * 4^N * N")
log(f"  {'N':>3}  {'||L_Dc||^2 (measured)':>22}  {'gamma^2 * 4^N * N':>22}  {'Match':>8}")
log(f"  {'-' * 60}")

for N in [3, 4, 5]:
    L_D = build_LD(N, GAMMA)
    Sg = N * GAMMA
    L_Dc = L_D + Sg * np.eye(4 ** N, dtype=complex)
    norm_sq = np.linalg.norm(L_Dc) ** 2
    predicted = GAMMA ** 2 * 4 ** N * N
    dev = abs(norm_sq - predicted)
    log(f"  {N:>3}  {norm_sq:>22.6f}  {predicted:>22.6f}  "
        f"{'PASS' if dev < 1e-8 else 'FAIL':>8}")

# Verify Step 2: bond-sum rule at each N
log()
log("  Step 2 verification: bond-sum rule (w_XY(a) + w_XY(b) = 2 at bond sites)")

paulis_1site = [I2, SX, SY, SZ]
pauli_names = ["I", "X", "Y", "Z"]
w_xy = [0, 1, 1, 0]  # XY weight for I, X, Y, Z

for N in [3, 4, 5]:
    # Build single-bond L_H in Pauli basis for bond (0,1)
    bonds = [(0, 1)]
    H = build_hamiltonian(N, bonds, J)
    L_H_single = build_LH(H)

    # Build Pauli basis for N sites
    pauli_strs = []
    for idx in range(4 ** N):
        digits = []
        tmp = idx
        for _ in range(N):
            digits.append(tmp % 4)
            tmp //= 4
        pauli_strs.append(digits)  # digits[k] is Pauli index at site k

    # Compute Pauli-basis L_H matrix elements
    d = 2 ** N
    pauli_mats = []
    for ps in pauli_strs:
        mat = paulis_1site[ps[0]]
        for k in range(1, N):
            mat = np.kron(mat, paulis_1site[ps[k]])
        pauli_mats.append(mat)

    violations = 0
    total_nonzero = 0
    for a_idx, (pa, psa) in enumerate(zip(pauli_mats, pauli_strs)):
        for b_idx, (pb, psb) in enumerate(zip(pauli_mats, pauli_strs)):
            # L_H in Pauli basis: L_ab = Tr(pa^dag [H, pb]) * (-1j) / d
            # Actually: Tr(pa^dag L_H(pb)) / d where L_H(pb) = -i[H, pb]
            comm = H @ pb - pb @ H
            elem = -1j * np.trace(pa.conj().T @ comm) / d
            if abs(elem) > 1e-12:
                total_nonzero += 1
                # Check bond-sum rule at sites 0,1
                w01_a = w_xy[psa[0]] + w_xy[psa[1]]
                w01_b = w_xy[psb[0]] + w_xy[psb[1]]
                if abs(w01_a + w01_b - 2) > 0.01:
                    violations += 1

    log(f"    N={N}: {total_nonzero} nonzero entries, "
        f"{violations} violations of bond-sum rule "
        f"({'PASS' if violations == 0 else 'FAIL'})")

# Verify Step 3+4: ||{L_H, L_Dc}||^2 = 4*gamma^2*(N-2)*||L_H||^2
log()
log("  Step 3+4 verification: ||{L_H, L_Dc}||^2 = 4*gamma^2*(N-2)*||L_H||^2")
log(f"  {'N':>3}  {'Topology':>10}  {'||anti||^2 (meas)':>18}  "
    f"{'4*g^2*(N-2)*||L_H||^2':>22}  {'ratio':>10}")
log(f"  {'-' * 70}")

for N in [2, 3, 4, 5]:
    for topo_name, bond_fn in [("chain", chain_bonds), ("complete", complete_bonds)]:
        bonds = bond_fn(N)
        H = build_hamiltonian(N, bonds, J)
        L_H = build_LH(H)
        L_D = build_LD(N, GAMMA)
        Sg = N * GAMMA
        L_Dc = L_D + Sg * np.eye(4 ** N, dtype=complex)
        anticomm = L_H @ L_Dc + L_Dc @ L_H
        anti_sq = np.linalg.norm(anticomm) ** 2
        predicted = 4 * GAMMA ** 2 * (N - 2) * np.linalg.norm(L_H) ** 2
        ratio = anti_sq / predicted if predicted > 0 else float("inf")
        log(f"  {N:>3}  {topo_name:>10}  {anti_sq:>18.6f}  "
            f"{predicted:>22.6f}  {ratio:>10.6f}")

# ============================================================
# SUMMARY TABLE
# ============================================================
log()
log("=" * 90)
log("COMPLETE DATA TABLE")
log("=" * 90)
log()
log(f"  {'N':>3}  {'Topology':>10}  {'Measured':>14}  {'Formula B':>14}  "
    f"{'Dev B':>10}  {'Formula A':>14}  {'Dev A':>10}")
log(f"  {'-' * 85}")

# Include known values from CROSS_TERM_TOPOLOGY (N=3,4)
all_data = []
for N in [2, 3, 4]:
    for topo_name, bond_fn in [("chain", chain_bonds)]:
        bonds = bond_fn(N)
        H = build_hamiltonian(N, bonds, J)
        L_H = build_LH(H)
        L_D = build_LD(N, GAMMA)
        _, rel, _, _ = compute_cross_term(L_H, L_D, N, GAMMA)
        fa = formula_A(N)
        fb = formula_B(N)
        da = abs(rel - fa)
        db = abs(rel - fb)
        log(f"  {N:>3}  {topo_name:>10}  {rel:>14.10f}  {fb:>14.10f}  "
            f"{db:>10.2e}  {fa:>14.10f}  {da:>10.2e}")
        all_data.append({"N": N, "topology": topo_name, "rel_ortho": float(rel),
                         "formula_A": float(fa), "formula_B": float(fb)})

for r in results:
    fa = formula_A(r["N"])
    fb = formula_B(r["N"])
    log(f"  {r['N']:>3}  {r['topology']:>10}  {r['rel_ortho']:>14.10f}  "
        f"{fb:>14.10f}  {r['dev_formula_B']:>10.2e}  "
        f"{fa:>14.10f}  {r['dev_formula_A']:>10.2e}")

# ============================================================
# SAVE JSON
# ============================================================
output_data = {
    "metadata": {
        "date": str(datetime.now()),
        "gamma": GAMMA, "J": J,
        "formula_A": "1/sqrt(N * 2^(N+1))",
        "formula_B": "sqrt((N-2) / (N * 4^(N-1)))",
    },
    "verdict": verdict,
    "results": results,
    "proof_identity": "||{L_H, L_Dc}||^2 = 4*gamma^2*(N-2)*||L_H||^2",
}

with open(JSON_PATH, "w", encoding="utf-8") as jf:
    json.dump(output_data, jf, indent=2)
log(f"\nSaved JSON: {JSON_PATH}")

# ============================================================
# DONE
# ============================================================
log()
log("=" * 90)
log(f"RUN COMPLETE. Verdict: {verdict}")
log(f"Log:  {LOG_PATH}")
log(f"JSON: {JSON_PATH}")
log("=" * 90)
_log_file.close()
