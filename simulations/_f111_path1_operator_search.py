"""
F111 Path 1: Search for a per-site Pi_y_par_D-flip operator that:
  1. Satisfies dissipator conjugation: M D[D] M^-1 = -D[D] - 2 gamma I.
  2. Achieves operator-level palindrome residual = 0 on every off-y_par(D)
     k=4 single-term Hamiltonian in the diagonal Klein cell.
  3. (Hopefully) does NOT achieve palindrome on on-y_par(D) k=4 templates,
     confirming the structural separation.

Per-site M structure: (permutation of {I, X, Y, Z}) x (phase per label).
We consider per-site permutations that satisfy the dissipator constraint
(swap immune-set with damped-set as two 2-cycles).

Per D:
  D = Z: permutations (I<->X, Y<->Z) and (I<->Y, X<->Z)
  D = X: permutations (I<->Z, X<->Y) and (I<->Y, X<->Z)
  D = Y: permutations (I<->X, Y<->Z) and (I<->Z, X<->Y)

Phase variants: 4 phase choices per label = 256 phase combos per permutation,
but reduced to {+1, -1, +i, -i} per label. We additionally enforce M^2 = +/- I_4
diagonal (so the 2-cycle phase products are +/-1).

The diagonal cell of D under D-dephasing has off-y_par(D) sector with:
  D=Z: #X odd, #Y odd, #Z even   (sigma = -2 in p^2 etc.)
  D=X: #Y odd, #Z odd, #X even
  D=Y: #X odd, #Z odd, #Y even

In each, ONE letter (the dephase letter D itself) is the "even-count axis".

Strategy:
  - Enumerate per-site M satisfying dissipator constraint.
  - For each, compute Q = M^otimes N and the operator-residual on each
    off-y_par k=4 H. Record which M achieves residual = 0 on ALL off-y_par.
  - If a winning M exists, examine its structure for closed-form interpretation.
"""

import os
from datetime import datetime
from itertools import product as iprod
from collections import Counter

import numpy as np

RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
os.makedirs(RESULTS_DIR, exist_ok=True)
OUT_LOG = os.path.join(RESULTS_DIR, "f111_path1_operator_search.txt")
f_log = open(OUT_LOG, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f_log.write(msg + "\n")
    f_log.flush()


# ---- Pauli matrices ----
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": sx, "Y": sy, "Z": sz}
LABELS = ["I", "X", "Y", "Z"]


def bit_a(letter): return 1 if letter in "XY" else 0
def bit_b(letter): return 1 if letter in "YZ" else 0


def term_klein(letters):
    return (sum(bit_a(c) for c in letters) % 2,
            sum(bit_b(c) for c in letters) % 2)


def term_y_par(letters):
    return sum(1 for c in letters if c == "Y") % 2


def diagonal_klein(dephase):
    return (bit_a(dephase), bit_b(dephase))


def y_par_dephase(dephase):
    return bit_a(dephase) & bit_b(dephase)


def build_pauli_op_full(letters):
    mat = PAULI[letters[0]]
    for k in range(1, len(letters)):
        mat = np.kron(mat, PAULI[letters[k]])
    return mat


def build_chain_k_body(N, template):
    k = len(template)
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    if k > N:
        return H
    for shift in range(N - k + 1):
        full = ["I"] * N
        for i in range(k):
            full[shift + i] = template[i]
        H = H + build_pauli_op_full(full)
    return H


def build_L(H, gamma, N, dephase):
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    D_op = PAULI[dephase]
    for k in range(N):
        ops = [I2] * N
        ops[k] = D_op
        Dk = ops[0]
        for o in ops[1:]:
            Dk = np.kron(Dk, o)
        L += gamma * (np.kron(Dk, Dk.conj()) - np.eye(d2))
    return L


def build_Q_Nsite(per_site, N):
    d = 2 ** N
    d2 = d * d
    basis_N = list(iprod(LABELS, repeat=N))
    label_to_idx = {bl: i for i, bl in enumerate(basis_N)}
    mats_N = []
    for bl in basis_N:
        mat = PAULI[bl[0]]
        for k in range(1, N):
            mat = np.kron(mat, PAULI[bl[k]])
        mats_N.append(mat)
    vecs_N = [m.flatten() for m in mats_N]

    Q = np.zeros((d2, d2), dtype=complex)
    for idx, bl in enumerate(basis_N):
        phase = 1.0
        tgt_labels = []
        for sl in bl:
            ph, t = per_site[sl]
            phase *= ph
            tgt_labels.append(t)
        tgt_idx = label_to_idx[tuple(tgt_labels)]
        Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
    return Q


def operator_palindrome_residual(Q, L, sigma):
    Q_inv = np.linalg.inv(Q)
    d2 = L.shape[0]
    M = Q @ L @ Q_inv + L + 2.0 * sigma * np.eye(d2)
    return np.linalg.norm(M)


def enumerate_k_terms_in_cell(k, klein, y_par):
    terms = []
    for s in iprod(LABELS, repeat=k):
        if all(c == "I" for c in s):
            continue
        seq = "".join(s)
        if term_klein(seq) == klein and term_y_par(seq) == y_par:
            terms.append(seq)
    return terms


# ---- Per-site permutations satisfying dissipator constraint per D ----
def dissipator_permutations(dephase):
    """For each dephase letter, return permutations of {I, X, Y, Z} that swap
    immune-set with damped-set as two 2-cycles."""
    if dephase == "Z":
        # immune={I,Z}, damped={X,Y}.
        # Pair (I,Z) with (X,Y): two options for pairing.
        # Option A: I<->X, Z<->Y  (i.e., (I,X) 2-cycle, (Y,Z) 2-cycle).
        # Option B: I<->Y, Z<->X  (i.e., (I,Y) 2-cycle, (X,Z) 2-cycle).
        return [
            {"I": "X", "X": "I", "Y": "Z", "Z": "Y"},
            {"I": "Y", "Y": "I", "X": "Z", "Z": "X"},
        ]
    elif dephase == "X":
        # immune={I,X}, damped={Y,Z}.
        return [
            {"I": "Y", "Y": "I", "X": "Z", "Z": "X"},
            {"I": "Z", "Z": "I", "X": "Y", "Y": "X"},
        ]
    elif dephase == "Y":
        # immune={I,Y}, damped={X,Z}.
        return [
            {"I": "X", "X": "I", "Y": "Z", "Z": "Y"},
            {"I": "Z", "Z": "I", "Y": "X", "X": "Y"},
        ]
    raise ValueError(dephase)


# ---- Enumerate phase variants for a given permutation ----
PHASES = [1.0, -1.0, 1j, -1j]


def enumerate_phase_variants(perm):
    """For each label, choose a phase from {+1, -1, +i, -i}. Return all M
    that have unitary 4x4 structure (each column has one non-zero unit entry)."""
    labels = ["I", "X", "Y", "Z"]
    for phase_tuple in iprod(PHASES, repeat=4):
        M = {label: (phase_tuple[idx], perm[label]) for idx, label in enumerate(labels)}
        yield M


def dissipator_check(M, dephase, gamma=0.05):
    """Verify M D[D] M^-1 = -D[D] - 2 gamma I as a 4x4 sanity check on 1 qubit."""
    Q1 = build_Q_Nsite(M, 1)
    Q1_inv = np.linalg.inv(Q1)
    d_op = PAULI[dephase]
    d2 = 4
    D_op = gamma * (np.kron(d_op, d_op.conj()) - np.eye(d2))
    lhs = Q1 @ D_op @ Q1_inv
    rhs = -D_op - 2 * gamma * np.eye(d2)
    return np.linalg.norm(lhs - rhs) < 1e-10


def hamiltonian_check_single(M, N, template, dephase, gamma=0.05):
    """Compute residual ||Q L_H Q^-1 - (-L_H)||_F for H = chain k-body of template."""
    sigma = N * gamma
    H = build_chain_k_body(N, template)
    L = build_L(H, gamma, N, dephase)
    Q = build_Q_Nsite(M, N)
    Q_inv = np.linalg.inv(Q)
    d2 = L.shape[0]
    diff = Q @ L @ Q_inv + L + 2 * sigma * np.eye(d2)
    return np.linalg.norm(diff)


def main():
    log("=" * 88)
    log("F111 Path 1: Search for per-site Pi_y_par_D-flip operator")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    N = 4
    k = 4
    gamma = 0.05

    log(f"\nSetup: N={N}, k={k}, gamma={gamma}")
    log("Goal: find per-site M such that Q = M^otimes 4 achieves palindrome")
    log("residual = 0 on every off-y_par(D) k=4 single-term Hamiltonian.")

    for dephase in ["Z", "X", "Y"]:
        log("\n" + "=" * 88)
        log(f"DEPHASE {dephase} (diagonal Klein = {diagonal_klein(dephase)}, "
            f"y_par(D) = {y_par_dephase(dephase)})")
        log("=" * 88)

        diag = diagonal_klein(dephase)
        off_y_par = 1 - y_par_dephase(dephase)
        on_y_par = y_par_dephase(dephase)

        off_templates = enumerate_k_terms_in_cell(k, diag, off_y_par)
        on_templates = enumerate_k_terms_in_cell(k, diag, on_y_par)
        log(f"\nOff-y_par templates: {len(off_templates)}")
        log(f"On-y_par templates: {len(on_templates)}")

        perms = dissipator_permutations(dephase)
        log(f"\nDissipator-valid permutations: {len(perms)}")

        # Enumerate phase variants for each permutation
        winners_off_only = []  # M achieving residual=0 on ALL off-y_par AND NON-zero on ALL on-y_par
        winners_off_all = []  # M achieving residual=0 on ALL off-y_par (regardless of on-y_par)
        all_off_zero_off_on_count = []  # (M, on_zero_count) for M with all off=0
        scan_count = 0
        for perm_idx, perm in enumerate(perms):
            log(f"\n-- Permutation {perm_idx + 1}: {perm} --")
            for M in enumerate_phase_variants(perm):
                scan_count += 1
                # Dissipator check (cheap, 1-qubit)
                if not dissipator_check(M, dephase, gamma):
                    continue
                # Hamiltonian check on all off-y_par
                off_residuals = []
                early_break = False
                for template in off_templates:
                    r = hamiltonian_check_single(M, N, template, dephase, gamma)
                    off_residuals.append(r)
                    if r > 1e-8:
                        early_break = True
                        break
                if early_break:
                    continue
                # All off_y_par = 0 — check on_y_par
                on_residuals = [hamiltonian_check_single(M, N, t, dephase, gamma)
                                for t in on_templates]
                on_zero_count = sum(1 for r in on_residuals if r < 1e-8)
                all_off_zero_off_on_count.append((M, on_zero_count, max(on_residuals)))
                # Strict winner: separates off from on
                if on_zero_count == 0:
                    winners_off_only.append(M)
                winners_off_all.append(M)
        log(f"\n  Scanned {scan_count} M variants for dephase {dephase}")
        log(f"  Winners (off all 0, on all > 0): {len(winners_off_only)}")
        log(f"  Winners (off all 0, regardless of on): {len(winners_off_all)}")

        if winners_off_only:
            log(f"\n  STRICT WINNERS (achieve palindrome on off-y_par ONLY):")
            for w_idx, w in enumerate(winners_off_only[:10]):
                log(f"    {w_idx + 1}. {format_per_site(w)}")
        elif winners_off_all:
            log(f"\n  WIDE WINNERS (zero on off-y_par; also some on-y_par):")
            for w_idx, (w, on_n, on_max) in enumerate(all_off_zero_off_on_count[:10]):
                log(f"    {w_idx + 1}. {format_per_site(w)} -- on_zero_count={on_n}/{len(on_templates)}, on_max_residual={on_max:.3e}")
        else:
            log(f"\n  NO M ACHIEVED ALL-OFF-ZERO; closest M analysis:")
            # Find M with minimum max-off residual
            min_max_off = np.inf
            best_M = None
            for perm in perms:
                for M in enumerate_phase_variants(perm):
                    if not dissipator_check(M, dephase, gamma):
                        continue
                    rs = [hamiltonian_check_single(M, N, t, dephase, gamma)
                          for t in off_templates]
                    if max(rs) < min_max_off:
                        min_max_off = max(rs)
                        best_M = M
            log(f"    Best M (min of max-off-residual): {format_per_site(best_M)}")
            log(f"    Max off-residual under best M: {min_max_off:.3e}")

    log("\n" + "=" * 88)
    log(f"Done: {datetime.now()}")
    log(f"Log: {OUT_LOG}")
    log("=" * 88)


def format_per_site(M):
    parts = []
    for L in ["I", "X", "Y", "Z"]:
        ph, tgt = M[L]
        ph_str = format_phase(ph)
        parts.append(f"{L}->{ph_str}{tgt}")
    return ", ".join(parts)


def format_phase(ph):
    if abs(ph - 1) < 1e-10: return "+"
    if abs(ph + 1) < 1e-10: return "-"
    if abs(ph - 1j) < 1e-10: return "+i"
    if abs(ph + 1j) < 1e-10: return "-i"
    return f"{ph:.3f}*"


if __name__ == "__main__":
    main()
