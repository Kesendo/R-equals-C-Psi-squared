"""
F111: combine Y^N (Hilbert-space conjugation) with various Pi operators
to construct a palindrome operator for off-y_par(D) sector.

Strategy:
- For off-y_par(D), find a Hilbert-space unitary V such that V H V^-1 = -H.
- V is one of {X^N, Y^N, Z^N} or product.
- Combine with Pi (canonical, 5bilinear, etc.) to get a Liouville-space
  operator R such that R L R^-1 = -L - 2 sigma I.

For Z-dephasing, off-y_par(Z) sector at k=4 N=4: V = Y^N gives V H V^-1 = -H.
Goal: find Pi such that Pi acts as "dissipator inverter only" so that
(V-superop) * Pi = palindrome.

Q_V L Q_V^-1 = -L_H + L_D (Hamiltonian flipped, dissipator unchanged).
We need Q_W L Q_W^-1 = -L - 2 sigma I.

So we want Q_W = Q_V * Pi where Pi achieves:
  Pi L_H Pi^-1 = L_H  (commute, since the -L_H comes from Q_V)
  Pi L_D Pi^-1 = -L_D - 2 sigma I  (the F1 dissipator side)

For Pi^2-D-even H, Pi_5bi anti-commutes ([H,.]). For Pi^2-D-odd H, it
neither commutes nor anti-commutes (residual 8 for both checks).

So a per-site Pi_5bi or canonical Pi won't commute with [H,.] for off-y_par.
We need a DIFFERENT operator P such that P [H,.] P^-1 = +[H,.] for off-y_par H
AND P L_D P^-1 = -L_D - 2 sigma I.

Test: composition (Q_V) * (P where P is one of: canonical Pi, Pi_5bi, identity).

For each composition, compute residual ||R L R^-1 + L + 2 sigma I|| on:
  (a) each off-y_par(D) single-term H at k=4 N=4
  (b) one representative off-y_par pair

If a composition achieves residual = 0 universally on off-y_par, that
gives the closed-form.
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
OUT_LOG = os.path.join(RESULTS_DIR, "f111_combined_operator_search.txt")
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


def pi5bilinear_per_site(dephase):
    if dephase == "Z":
        return {"I": (1, "X"), "X": (-1, "I"), "Y": (1j, "Z"), "Z": (-1j, "Y")}
    elif dephase == "X":
        return {"I": (1, "Z"), "Z": (-1, "I"), "X": (-1j, "Y"), "Y": (1j, "X")}
    elif dephase == "Y":
        return {"I": (1, "X"), "X": (-1, "I"), "Y": (-1j, "Z"), "Z": (1j, "Y")}
    raise ValueError(dephase)


def canonical_pi_per_site(dephase):
    if dephase == "Z":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}
    elif dephase == "X":
        return {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}
    elif dephase == "Y":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}
    raise ValueError(dephase)


def build_Q_V_superop(V, N):
    """Q_V = V ⊗ V* in vec basis (conjugation by Hilbert-space unitary V).
    V is a d x d unitary on the full N-qubit Hilbert space."""
    return np.kron(V, V.conj())


def enumerate_k_terms_in_cell(k, klein, y_par):
    terms = []
    for s in iprod(LABELS, repeat=k):
        if all(c == "I" for c in s):
            continue
        seq = "".join(s)
        if term_klein(seq) == klein and term_y_par(seq) == y_par:
            terms.append(seq)
    return terms


def residual_check(R, L, sigma):
    R_inv = np.linalg.inv(R)
    d2 = L.shape[0]
    M = R @ L @ R_inv + L + 2 * sigma * np.eye(d2)
    return np.linalg.norm(M)


def main():
    log("=" * 88)
    log("F111 combined operator search: Q_V * Pi compositions")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    N = 4
    k = 4
    gamma = 0.05
    sigma = N * gamma

    log(f"\nSetup: N={N}, k={k}, gamma={gamma}, sigma={sigma}")

    # Build flip operators V = X^N, Y^N, Z^N
    V_X = build_pauli_op_full("X" * N)
    V_Y = build_pauli_op_full("Y" * N)
    V_Z = build_pauli_op_full("Z" * N)
    Q_VX = build_Q_V_superop(V_X, N)
    Q_VY = build_Q_V_superop(V_Y, N)
    Q_VZ = build_Q_V_superop(V_Z, N)

    flip_letters = {
        "X": (V_X, Q_VX),
        "Y": (V_Y, Q_VY),
        "Z": (V_Z, Q_VZ),
    }

    for dephase in ["Z", "X", "Y"]:
        log("\n" + "=" * 88)
        log(f"DEPHASE {dephase}")
        log("=" * 88)

        diag = diagonal_klein(dephase)
        off_y_par = 1 - y_par_dephase(dephase)
        off_templates = enumerate_k_terms_in_cell(k, diag, off_y_par)
        log(f"\nOff-y_par(D)={off_y_par} templates: {len(off_templates)}")

        # Verify flip: for each flip letter V, check Q_V H Q_V^-1 = -H (or +H) on off-y_par
        log(f"\n-- Step 1: verify Q_V acts as H -> -H on off-y_par(D) Hamiltonians --")
        for V_letter, (V, Q_V) in flip_letters.items():
            n_flip = 0
            n_preserve = 0
            for t in off_templates:
                H = build_chain_k_body(N, t)
                H_conj = V @ H @ np.linalg.inv(V)
                if np.allclose(H_conj, -H):
                    n_flip += 1
                elif np.allclose(H_conj, H):
                    n_preserve += 1
            log(f"  V={V_letter}^{N}: flip={n_flip}/{len(off_templates)}, "
                f"preserve={n_preserve}/{len(off_templates)}")

        # Now: composition R = Q_V * Pi_X. Try Pi_canonical (D-deph) and Pi_5bi (D-deph).
        Q_can = build_Q_Nsite(canonical_pi_per_site(dephase), N)
        Q_5bi = build_Q_Nsite(pi5bilinear_per_site(dephase), N)

        # Find Q_V that flips ALL off-y_par templates
        log(f"\n-- Step 2: identify universal flip letter for off-y_par(D)={off_y_par} --")
        universal_flips = []
        for V_letter, (V, Q_V) in flip_letters.items():
            ok = True
            for t in off_templates:
                H = build_chain_k_body(N, t)
                H_conj = V @ H @ np.linalg.inv(V)
                if not np.allclose(H_conj, -H):
                    ok = False
                    break
            if ok:
                universal_flips.append(V_letter)
        log(f"  Universal flip letters: {universal_flips}")
        if not universal_flips:
            log(f"  -> no universal flip letter; F111 closed form would require ON-y_par Pi")
            continue

        # Try R = Q_V * Pi for V in universal_flips, Pi in {canonical, 5bi}
        log(f"\n-- Step 3: test R = Q_V * Pi on off-y_par templates --")
        compositions = []
        for V_letter in universal_flips:
            V, Q_V = flip_letters[V_letter]
            for pi_name, Pi in [("Pi_can", Q_can), ("Pi_5bi", Q_5bi),
                                 ("Pi_can_inv", np.linalg.inv(Q_can)),
                                 ("Pi_5bi_inv", np.linalg.inv(Q_5bi))]:
                R = Q_V @ Pi
                compositions.append((f"Q_V({V_letter}^{N}) * {pi_name}", R))
                R2 = Pi @ Q_V
                compositions.append((f"{pi_name} * Q_V({V_letter}^{N})", R2))

        for comp_name, R in compositions:
            try:
                residuals = []
                for t in off_templates:
                    H = build_chain_k_body(N, t)
                    L = build_L(H, gamma, N, dephase)
                    res = residual_check(R, L, sigma)
                    residuals.append(res)
                n_zero = sum(1 for r in residuals if r < 1e-8)
                log(f"  {comp_name}: n_zero={n_zero}/{len(off_templates)}, "
                    f"max={max(residuals):.3e}, median={np.median(residuals):.3e}")
            except Exception as e:
                log(f"  {comp_name}: ERROR {e}")

    log("\n" + "=" * 88)
    log(f"Done: {datetime.now()}")
    log(f"Log: {OUT_LOG}")
    log("=" * 88)


if __name__ == "__main__":
    main()
