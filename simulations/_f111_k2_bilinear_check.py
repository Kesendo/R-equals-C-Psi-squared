"""
F111: classify each Π²-D-ODD bilinear in the diagonal cell of D, single-term, at k=2 N varied.
Specifically:
  Z-deph: diagonal Klein (0,1) Π²-Z-odd bilinears = {XY, YX}. y_par(XY)=y_par(YX)=1 (off-y_par(Z)).
  X-deph: diagonal Klein (1,0) Π²-X-odd bilinears = {YZ, ZY}. y_par=1 (off-y_par(X)).
  Y-deph: diagonal Klein (1,1) Π²-Y-odd bilinears = {XZ, ZX}. y_par=0 (off-y_par(Y)).

For each, classify at N=2, 3, 4, 5 under chain dephasing.

If these bilinears are ALWAYS SOFT (never hard), then the off-y_par single-term
k=4 N=4 case follows by structural reduction.

For PAIR Hamiltonians: check whether the pair sum is also always soft.
"""

import os
from datetime import datetime
from itertools import product as iprod
from collections import Counter, defaultdict

import numpy as np

RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
os.makedirs(RESULTS_DIR, exist_ok=True)
OUT_LOG = os.path.join(RESULTS_DIR, "f111_k2_bilinear_check.txt")
f_log = open(OUT_LOG, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f_log.write(msg + "\n")
    f_log.flush()


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": sx, "Y": sy, "Z": sz}
LABELS = ["I", "X", "Y", "Z"]


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


def canonical_pi_per_site(dephase):
    if dephase == "Z":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}
    elif dephase == "X":
        return {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}
    elif dephase == "Y":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}
    raise ValueError(dephase)


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


def spec_pairs(evals, sigma, tol):
    used = [False] * len(evals)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma
        best_j, best_d = -1, np.inf
        for j in range(len(evals)):
            if used[j]:
                continue
            d = abs(evals[j] - target)
            if d < best_d:
                best_d = d
                best_j = j
        if best_j < 0:
            return False, np.inf
        used[i] = True
        if best_j != i:
            used[best_j] = True
        if best_d > max_err:
            max_err = best_d
    return max_err < tol, max_err


def classify(N, template, dephase, gamma=0.05,
             op_tol=1e-10, spec_tol=1e-6):
    H = build_chain_k_body(N, template)
    L = build_L(H, gamma, N, dephase)
    sigma = N * gamma
    Pi = build_Q_Nsite(canonical_pi_per_site(dephase), N)
    Pi_inv = np.linalg.inv(Pi)
    M = Pi @ L @ Pi_inv + L + 2 * sigma * np.eye(L.shape[0])
    op_norm = np.linalg.norm(M)
    if op_norm < op_tol:
        return "truly", op_norm, 0.0
    evals = np.linalg.eigvals(L)
    paired, max_err = spec_pairs(evals, sigma, spec_tol)
    cls = "soft" if paired else "hard"
    return cls, op_norm, max_err


def classify_pair(N, templates, dephase, gamma=0.05,
                  op_tol=1e-10, spec_tol=1e-6):
    H = sum(build_chain_k_body(N, t) for t in templates)
    L = build_L(H, gamma, N, dephase)
    sigma = N * gamma
    Pi = build_Q_Nsite(canonical_pi_per_site(dephase), N)
    Pi_inv = np.linalg.inv(Pi)
    M = Pi @ L @ Pi_inv + L + 2 * sigma * np.eye(L.shape[0])
    op_norm = np.linalg.norm(M)
    if op_norm < op_tol:
        return "truly", op_norm, 0.0
    evals = np.linalg.eigvals(L)
    paired, max_err = spec_pairs(evals, sigma, spec_tol)
    cls = "soft" if paired else "hard"
    return cls, op_norm, max_err


def main():
    log("=" * 88)
    log("F111: Pi2-D-odd bilinear classification (off-y_par(D) at k=2)")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    gamma = 0.05
    N_range = [2, 3, 4, 5]

    off_bilinears = {
        "Z": ["XY", "YX"],
        "X": ["YZ", "ZY"],
        "Y": ["XZ", "ZX"],
    }

    log("\nSingle-bilinear classification across N:")
    for dephase, bilinears in off_bilinears.items():
        log(f"\nDephase {dephase}, off-y_par bilinears: {bilinears}")
        for bil in bilinears:
            for N in N_range:
                cls, op_norm, spec_err = classify(N, bil, dephase, gamma)
                log(f"  {bil} at N={N}: cls={cls:>5s}, "
                    f"||M||={op_norm:.3e}, spec_err={spec_err:.3e}")

    log("\nPair classification (pairs of off-y_par bilinears):")
    for dephase, bilinears in off_bilinears.items():
        log(f"\nDephase {dephase}, off-y_par bilinear pairs:")
        for b1 in bilinears:
            for b2 in bilinears:
                if b1 > b2:
                    continue
                for N in N_range:
                    cls, op_norm, spec_err = classify_pair(N, [b1, b2], dephase, gamma)
                    log(f"  ({b1}, {b2}) at N={N}: cls={cls:>5s}, "
                        f"||M||={op_norm:.3e}, spec_err={spec_err:.3e}")

    log("\n" + "=" * 88)
    log(f"Done: {datetime.now()}")
    log(f"Log: {OUT_LOG}")
    log("=" * 88)


if __name__ == "__main__":
    main()
