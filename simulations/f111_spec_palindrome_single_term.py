"""
F111 verification: classify each k=4 single-term Hamiltonian in the diagonal
Klein cell (split by y_par sector) as truly/soft/hard. Confirms the structural
claim at the single-term level: off-y_par(D) is ALL SOFT (never hard) at k=4
N=4.

This is faster than the 2-term pair enumeration (only 64 templates per dephase
vs 1500+ pairs).
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
OUT_LOG = os.path.join(RESULTS_DIR, "f111_spec_palindrome_single_term.txt")
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


def build_pi_canonical(per_site, N):
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


def canonical_pi_per_site(dephase):
    if dephase == "Z":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}
    elif dephase == "X":
        return {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}
    elif dephase == "Y":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}
    raise ValueError(dephase)


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
    Pi = build_pi_canonical(canonical_pi_per_site(dephase), N)
    Pi_inv = np.linalg.inv(Pi)
    M = Pi @ L @ Pi_inv + L + 2 * sigma * np.eye(L.shape[0])
    op_norm = np.linalg.norm(M)
    if op_norm < op_tol:
        return "truly", op_norm, 0.0
    evals = np.linalg.eigvals(L)
    paired, max_err = spec_pairs(evals, sigma, spec_tol)
    cls = "soft" if paired else "hard"
    return cls, op_norm, max_err


def enumerate_k_terms_in_cell(k, klein, y_par):
    terms = []
    for s in iprod(LABELS, repeat=k):
        if all(c == "I" for c in s):
            continue
        seq = "".join(s)
        if term_klein(seq) == klein and term_y_par(seq) == y_par:
            terms.append(seq)
    return terms


def main():
    log("=" * 88)
    log("F111 single-term verification: classify each k=4 H in diagonal cell")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    N = 4
    k = 4
    gamma = 0.05

    log(f"\nSetup: N={N}, k={k}, gamma={gamma}")
    log("Classify each k=4 single-term Hamiltonian in the diagonal Klein cell,")
    log("split by y_par(template). The F111 claim implies: off-y_par(D) sector")
    log("should be all SOFT (no hard).")

    for dephase in ["Z", "X", "Y"]:
        log("\n" + "=" * 88)
        log(f"DEPHASE {dephase} (diagonal Klein = {diagonal_klein(dephase)}, "
            f"y_par(D) = {y_par_dephase(dephase)})")
        log("=" * 88)

        diag = diagonal_klein(dephase)
        on_y_par = y_par_dephase(dephase)
        off_y_par = 1 - on_y_par

        for y_par_label, y_par_val in [("ON-y_par(D)", on_y_par),
                                        ("OFF-y_par(D)", off_y_par)]:
            templates = enumerate_k_terms_in_cell(k, diag, y_par_val)
            log(f"\n-- {y_par_label} (y_par={y_par_val}): {len(templates)} templates --")
            classes = Counter()
            sample_records = defaultdict(list)
            for template in templates:
                cls, op_norm, spec_err = classify(N, template, dephase, gamma)
                classes[cls] += 1
                sample_records[cls].append((template, op_norm, spec_err))
            for cls in ["truly", "soft", "hard"]:
                cnt = classes[cls]
                log(f"    {cls:>5s}: {cnt:>4d}")
            for cls in ["hard", "soft"]:
                if sample_records[cls]:
                    log(f"    Sample {cls} templates:")
                    for t, on, se in sample_records[cls][:5]:
                        log(f"      {t}: ||M||={on:.3e}, max_spec_err={se:.3e}")

    log("\n" + "=" * 88)
    log(f"Done: {datetime.now()}")
    log(f"Log: {OUT_LOG}")
    log("=" * 88)


if __name__ == "__main__":
    main()
