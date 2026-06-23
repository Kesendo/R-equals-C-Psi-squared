"""
F111 Path 2 test: Does the existing F108 Pi_5bilinear (per dephase) achieve
operator-level palindrome on k=4 off-y_par(D) Hamiltonians in the diagonal
Klein cell?

If YES on a representative sample (say all 2-term pairs in the off-y_par sector
of one diagonal cell at N=4), then the closed-form would follow as a corollary
of F108: the same Pi_5bilinear that handles Pi^2-D-even bilinears also handles
the off-y_par(D) Pi^2-D-odd Hamiltonians.

We will test:
  - Single-term Hamiltonians: each off-y_par(D) k=4 template, placed as chain
    k-body H, with the matching D-dephasing.
  - Two-term pair Hamiltonians: representative pairs from the off-y_par sector.

If even single-term H yields residual = 0 for off-y_par(D) templates,
that immediately closes F111: the existing F108 Pi_5bilinear ALSO handles
the off-y_par sector of the diagonal Klein cell.
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
OUT_LOG = os.path.join(RESULTS_DIR, "f111_path2_pi5bilinear_test.txt")
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


def bit_a(letter):
    return 1 if letter in "XY" else 0


def bit_b(letter):
    return 1 if letter in "YZ" else 0


def term_klein(letters):
    return (sum(bit_a(c) for c in letters) % 2,
            sum(bit_b(c) for c in letters) % 2)


def term_y_par(letters):
    return sum(1 for c in letters if c == "Y") % 2


def diagonal_klein(dephase):
    return (bit_a(dephase), bit_b(dephase))


def y_par_dephase(dephase):
    return bit_a(dephase) & bit_b(dephase)


# ---- Pi_5bilinear per-site action per dephase (F108 Part 1+2+3) ----
def pi5bilinear_per_site(dephase):
    if dephase == "Z":
        # Z-deph: I -> +X, X -> -I, Y -> +iZ, Z -> -iY
        return {"I": (1, "X"), "X": (-1, "I"), "Y": (1j, "Z"), "Z": (-1j, "Y")}
    elif dephase == "X":
        # X-deph: I -> +Z, Z -> -I, X -> -iY, Y -> +iX
        return {"I": (1, "Z"), "Z": (-1, "I"), "X": (-1j, "Y"), "Y": (1j, "X")}
    elif dephase == "Y":
        # Y-deph: I -> +X, X -> -I, Y -> -iZ, Z -> +iY
        return {"I": (1, "X"), "X": (-1, "I"), "Y": (-1j, "Z"), "Z": (1j, "Y")}
    raise ValueError(dephase)


def canonical_pi_per_site(dephase):
    """Canonical PiOperator (NOT 5bilinear) for sanity comparisons."""
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


def operator_palindrome_residual(Q, L, sigma):
    Q_inv = np.linalg.inv(Q)
    d2 = L.shape[0]
    M = Q @ L @ Q_inv + L + 2.0 * sigma * np.eye(d2)
    return np.linalg.norm(M)


def enumerate_k_terms_in_cell(k, klein, y_par):
    """All k-letter strings (non-trivial) with given Klein index and y_par."""
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
    log("F111 Path 2: Pi_5bilinear residual on k=4 off-y_par(D) Hamiltonians")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    N = 4
    k = 4
    gamma = 0.05
    sigma = N * gamma

    log(f"\nSetup: N={N}, k={k}, gamma={gamma}, sigma={sigma}")
    log("\nF111 claim: at k=4 N=4, every F87-hard pair in diagonal Klein cell")
    log("has y_par = y_par(D). Equivalently: off-y_par(D) sector contains NO hard pairs.")
    log("\nIf Pi_5bilinear achieves residual = 0 on off-y_par(D) k=4 Hamiltonians,")
    log("then by Pi-operator-level palindrome -> spec(L) palindromic -> NOT F87-hard.")
    log("\nPath 2 = closed-form via existing F108 operator without new construction.")

    for dephase in ["Z", "X", "Y"]:
        log("\n" + "=" * 88)
        log(f"DEPHASE {dephase} (diagonal Klein = {diagonal_klein(dephase)}, "
            f"y_par(D) = {y_par_dephase(dephase)})")
        log("=" * 88)

        diag = diagonal_klein(dephase)
        off_y_par = 1 - y_par_dephase(dephase)
        on_y_par = y_par_dephase(dephase)

        # Build operators once
        per_site_5bi = pi5bilinear_per_site(dephase)
        per_site_can = canonical_pi_per_site(dephase)
        Q_5bi = build_Q_Nsite(per_site_5bi, N)
        Q_can = build_Q_Nsite(per_site_can, N)

        # ----- Test 1: single-term k=4 Hamiltonian, OFF-y_par sector -----
        log(f"\n-- Test 1: single-term k=4 H, OFF-y_par sector "
            f"(Klein={diag}, y_par={off_y_par}) --")
        off_terms = enumerate_k_terms_in_cell(k, diag, off_y_par)
        log(f"  Number of off-y_par templates: {len(off_terms)}")
        residuals_off_5bi = []
        residuals_off_can = []
        for term in off_terms:
            H = build_chain_k_body(N, term)
            L = build_L(H, gamma, N, dephase)
            r_5bi = operator_palindrome_residual(Q_5bi, L, sigma)
            r_can = operator_palindrome_residual(Q_can, L, sigma)
            residuals_off_5bi.append(r_5bi)
            residuals_off_can.append(r_can)
        log(f"  Pi_5bilinear residual: min={min(residuals_off_5bi):.3e}, "
            f"max={max(residuals_off_5bi):.3e}, median={np.median(residuals_off_5bi):.3e}")
        log(f"  Pi_canonical residual: min={min(residuals_off_can):.3e}, "
            f"max={max(residuals_off_can):.3e}, median={np.median(residuals_off_can):.3e}")
        zero_count_5bi = sum(1 for r in residuals_off_5bi if r < 1e-8)
        log(f"  Pi_5bilinear EXACT (residual<1e-8) count: {zero_count_5bi} / "
            f"{len(residuals_off_5bi)}")

        # Show first 5 off-y_par templates with their residuals
        log(f"  First 5 off-y_par templates:")
        for term, r5, rc in list(zip(off_terms, residuals_off_5bi, residuals_off_can))[:5]:
            log(f"    {term}: 5bi={r5:.3e}, can={rc:.3e}")

        # ----- Test 2: single-term k=4 H, ON-y_par sector for comparison -----
        log(f"\n-- Test 2: single-term k=4 H, ON-y_par sector "
            f"(Klein={diag}, y_par={on_y_par}) for comparison --")
        on_terms = enumerate_k_terms_in_cell(k, diag, on_y_par)
        log(f"  Number of on-y_par templates: {len(on_terms)}")
        residuals_on_5bi = []
        residuals_on_can = []
        for term in on_terms:
            H = build_chain_k_body(N, term)
            L = build_L(H, gamma, N, dephase)
            r_5bi = operator_palindrome_residual(Q_5bi, L, sigma)
            r_can = operator_palindrome_residual(Q_can, L, sigma)
            residuals_on_5bi.append(r_5bi)
            residuals_on_can.append(r_can)
        log(f"  Pi_5bilinear residual: min={min(residuals_on_5bi):.3e}, "
            f"max={max(residuals_on_5bi):.3e}, median={np.median(residuals_on_5bi):.3e}")
        log(f"  Pi_canonical residual: min={min(residuals_on_can):.3e}, "
            f"max={max(residuals_on_can):.3e}, median={np.median(residuals_on_can):.3e}")
        zero_count_5bi_on = sum(1 for r in residuals_on_5bi if r < 1e-8)
        log(f"  Pi_5bilinear EXACT (residual<1e-8) count: {zero_count_5bi_on} / "
            f"{len(residuals_on_5bi)}")
        log(f"  First 5 on-y_par templates:")
        for term, r5, rc in list(zip(on_terms, residuals_on_5bi, residuals_on_can))[:5]:
            log(f"    {term}: 5bi={r5:.3e}, can={rc:.3e}")

        # ----- Summary for this dephase -----
        log(f"\n-- Summary for {dephase}-dephase --")
        log(f"  Pi_5bilinear on off-y_par single-term: "
            f"{zero_count_5bi}/{len(residuals_off_5bi)} EXACT")
        log(f"  Pi_5bilinear on on-y_par  single-term: "
            f"{zero_count_5bi_on}/{len(residuals_on_5bi)} EXACT")
        if (zero_count_5bi == len(residuals_off_5bi) and
                zero_count_5bi_on < len(residuals_on_5bi)):
            log(f"  => Path 2 SUCCESS for {dephase}: Pi_5bilinear separates "
                f"off-y_par (EXACT) from on-y_par (NON-zero)")
        elif zero_count_5bi == len(residuals_off_5bi):
            log(f"  => Pi_5bilinear EXACT on off-y_par; need check on-y_par as well")
        else:
            log(f"  => Path 2 FAILS for {dephase}: Pi_5bilinear NOT EXACT on all off-y_par")

    log("\n" + "=" * 88)
    log(f"Done: {datetime.now()}")
    log(f"Log: {OUT_LOG}")
    log("=" * 88)


if __name__ == "__main__":
    main()
