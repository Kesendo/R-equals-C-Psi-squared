"""F108 Part 2: Pi-family scan for pure-Pi^2_X-even non-truly Hamiltonians + X-dephasing.

Mirrors f108_part1_pi_family_scan.py for the BitA-axis case. The conjecture:
the Z-deph Pi_5bilinear's mechanism lifts to X-deph by applying the X<->Z label
swap, giving Pi_5bilinear_X with per-site action

  I -> +Z,  Z -> -I,  X -> -iY,  Y -> +iX

(same I<->Z, X<->Y permutation as canonical X-deph Pi, with phase flips on Z->I
and Y->X back-arrows). Pi^2_X-even bilinears = bit_a-even pairs (excluding
I-containing trivial cases) = {ZZ, YY, YX, XY, XX}. Pure-Pi^2_X-even NON-truly
(F85 X-deph truly = #X even AND #Y even) bilinears = {YX, XY} (single-site
non-truly).

Scan checks if Pi_5bilinear_X gives EXACT operator-level palindrome
Pi * L * Pi^-1 + L + 2*sigma*I = 0 for every pure-Pi^2_X-even non-truly pair at
N=3,4,5 under X-dephasing.

Output: console + simulations/results/f108_part2_x_dephasing_scan.txt
"""

import os
from datetime import datetime
from itertools import combinations_with_replacement, product as iprod

import numpy as np

# ============================================================
# Output setup
# ============================================================
RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
OUT = os.path.join(RESULTS_DIR, "f108_part2_x_dephasing_scan.txt")
os.makedirs(RESULTS_DIR, exist_ok=True)
f_log = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f_log.write(msg + "\n")
    f_log.flush()


# ============================================================
# Pauli matrices
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": sx, "Y": sy, "Z": sz}
LABELS = ["I", "X", "Y", "Z"]


# ============================================================
# bit_a parity + F85 truly classification under X-dephasing
# ============================================================
def bit_a(bilinear):
    """Pi^2_X parity bit: #X + #Y mod 2."""
    return sum(1 for c in bilinear if c in "XY") % 2


def is_truly_xdeph(bilinear):
    """F85 X-deph truly: #X even AND #Y even."""
    nx = sum(1 for c in bilinear if c == "X")
    ny = sum(1 for c in bilinear if c == "Y")
    return (nx % 2 == 0) and (ny % 2 == 0)


# ============================================================
# Site operator builder
# ============================================================
def site_op(op, s, N):
    ops = [I2] * N
    ops[s] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H_chain(N, bonds, comps):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for c, J in comps.items():
            if J == 0:
                continue
            H += J * site_op(PAULI[c[0]], i, N) @ site_op(PAULI[c[1]], j, N)
    return H


def build_L_xdeph(H, gamma, N):
    """L = -i [H, .] + sum_k gamma * D[X_k] in vec (row-major) basis."""
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Xk = site_op(sx, k, N)
        L += gamma * (np.kron(Xk, Xk.conj()) - np.eye(d2))
    return L


# ============================================================
# Per-site Pi map as dict: label -> (phase, target_label)
# ============================================================
def build_Q_Nsite(M, N):
    """Q = M^{otimes N} as d^2 x d^2 superoperator in standard vec basis."""
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
            ph, t = M[sl]
            phase *= ph
            tgt_labels.append(t)
        tgt_idx = label_to_idx[tuple(tgt_labels)]
        Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
    return Q


# ============================================================
# Pi family library
# ============================================================
def Pi_X_5bilinear_predicted():
    """Predicted X-deph 5-bilinear variant by X<->Z swap from Z-deph variant.
    I -> +Z, Z -> -I, X -> -iY, Y -> +iX
    Same I<->Z, X<->Y permutation as canonical X-deph Pi; phase flips on Z->I
    and Y->X back arrows (mirror of Z-deph variant's X->-I and Z->-iY flips)."""
    return {"I": (1, "Z"), "Z": (-1, "I"), "X": (-1j, "Y"), "Y": (1j, "X")}


def Pi_X_canonical():
    """Canonical X-deph P1 per PiOperator.cs:
    I -> +Z, Z -> +I, X -> -iY, Y -> -iX"""
    return {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}


def Pi_X_variant_v2():
    """Alternative variant: Z<->I flipped, X<->Y normal."""
    return {"I": (1, "Z"), "Z": (-1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}


def Pi_X_variant_v3():
    """Alternative: Z<->I normal, X<->Y flipped."""
    return {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (1j, "X")}


def Pi_X_variant_v4():
    """Alternative: both flipped opposite way."""
    return {"I": (1, "Z"), "Z": (-1, "I"), "X": (1j, "Y"), "Y": (-1j, "X")}


# ============================================================
# Palindrome checks
# ============================================================
def spec_palindrome_score(L, sigma, tol=1e-5):
    """Fraction of eigenvalues with palindromic partner."""
    evals = np.linalg.eigvals(L)
    target_sum = -2.0 * sigma
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    paired = 0
    for i in range(n):
        if used[i]:
            continue
        partner = target_sum - evals[i]
        if abs(evals[i].real - target_sum / 2) < tol and abs(evals[i].imag) < tol:
            used[i] = True
            paired += 1
            continue
        best_j, best_d = -1, np.inf
        for j in range(n):
            if used[j] or j == i:
                continue
            dd = abs(evals[j] - partner)
            if dd < best_d:
                best_d = dd
                best_j = j
        if best_j >= 0 and best_d < tol:
            used[i] = True
            used[best_j] = True
            paired += 2
    return paired / n


def operator_palindrome_residual(Q, L, sigma):
    """||Q L Q^{-1} + L + 2 sigma I||_F."""
    try:
        Q_inv = np.linalg.inv(Q)
    except np.linalg.LinAlgError:
        return np.inf
    d2 = L.shape[0]
    M = Q @ L @ Q_inv + L + 2.0 * sigma * np.eye(d2)
    return np.linalg.norm(M)


# ============================================================
# Bilinear classification under X-dephasing
# ============================================================
# Pi^2_X-even = bit_a=0 = #X + #Y in bilinear is even
# Allowed pairs of letters: both in {I, Z} (bit_a 0) or both in {X, Y} (bit_a 1)
# Excluding I-containing trivial: {ZZ, XX, XY, YX, YY}
PI2X_EVEN_BILINEARS = ["ZZ", "XX", "XY", "YX", "YY"]
# Truly under X-deph: #X even AND #Y even -> {ZZ, XX, YY} (XX: #X=2 even, #Y=0; YY: #Y=2 even, #X=0)
TRULY_BILINEARS_XDEPH = ["ZZ", "XX", "YY"]
# Non-truly Pi^2_X-even: {XY, YX} (both have #X=1, #Y=1, both odd)
NONTRULY_PI2XEVEN_BILINEARS = ["XY", "YX"]


def pair_is_pure_pi2x_even_non_truly(t1, t2):
    """Both Pi^2_X-even, at least one NON-truly under X-deph."""
    if bit_a(t1) != 0 or bit_a(t2) != 0:
        return False
    if is_truly_xdeph(t1) and is_truly_xdeph(t2):
        return False
    return True


# ============================================================
# MAIN
# ============================================================
def main():
    log("=" * 88)
    log("F108 Part 2: Pi-family scan for pure-Pi^2_X-even non-truly H + X-dephasing")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    log("\nPi^2_X-even (bit_a=0) 2-site bilinears: {ZZ, XX, XY, YX, YY}")
    log("Truly under X-deph (F85): #X even AND #Y even -> {ZZ, XX, YY}")
    log("Pure-Pi^2_X-even NON-truly bilinears: {XY, YX}")

    # ============================================================
    # Phase 1: enumerate the pairs
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 1: Enumerate pure-Pi^2_X-even non-truly pairs")
    log("=" * 88)

    pairs = []

    # Single bilinears: only non-truly ones (XY, YX)
    for b in NONTRULY_PI2XEVEN_BILINEARS:
        pairs.append((b, {b: 1.0}))

    # Two-term distinct combos: both Pi^2_X-even, at least one non-truly
    for t1, t2 in combinations_with_replacement(PI2X_EVEN_BILINEARS, 2):
        if t1 == t2:
            continue
        if pair_is_pure_pi2x_even_non_truly(t1, t2):
            label = f"{t1}+{t2}"
            pairs.append((label, {t1: 1.0, t2: 1.0}))

    log(f"\n  Total pairs: {len(pairs)}")
    for label, _ in pairs:
        log(f"    {label}")

    # ============================================================
    # Phase 2: For each pair, scan Pi families
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 2: Pi-family operator-level palindrome scan (X-dephasing)")
    log("=" * 88)

    gamma = 0.05
    N_list = [3, 4, 5]

    families = {
        "X_canonical": Pi_X_canonical(),
        "X_5bilinear_predicted": Pi_X_5bilinear_predicted(),
        "X_variant_v2": Pi_X_variant_v2(),
        "X_variant_v3": Pi_X_variant_v3(),
        "X_variant_v4": Pi_X_variant_v4(),
    }

    results = {}

    for N in N_list:
        log(f"\n--- N = {N} (chain, open BC, gamma = {gamma}) ---")
        bonds = [(i, i + 1) for i in range(N - 1)]
        sigma = N * gamma

        Q_per_family = {}
        for fname, M in families.items():
            Q_per_family[fname] = build_Q_Nsite(M, N)

        header = f"  {'Pair':<10} {'SpecPal':>8}  " + \
                 "  ".join(f"{fname[:21]:>21}" for fname in families.keys())
        log(header)
        log(f"  {'-' * 10} {'-' * 8}  " + "  ".join("-" * 21 for _ in families))

        for label, comps in pairs:
            H = build_H_chain(N, bonds, comps)
            L = build_L_xdeph(H, gamma, N)
            spec_score = spec_palindrome_score(L, sigma, tol=1e-5)
            row = {}
            for fname in families:
                resid = operator_palindrome_residual(Q_per_family[fname], L, sigma)
                row[fname] = resid
            results[(label, N)] = (spec_score, row)

            spec_str = f"{spec_score*100:>6.1f}%"
            resid_strs = []
            for fname in families:
                r = row[fname]
                if r < 1e-8:
                    resid_strs.append(f"{'EXACT':>21}")
                elif r > 1e10:
                    resid_strs.append(f"{'SING':>21}")
                else:
                    resid_strs.append(f"{r:>21.3e}")
            log(f"  {label:<10} {spec_str:>8}  " + "  ".join(resid_strs))

    # ============================================================
    # Phase 3: Pattern extraction
    # ============================================================
    log("\n" + "=" * 88)
    log(f"Phase 3: Single-family coverage test (X-dephasing)")
    log("=" * 88)

    for fname in families:
        not_covered = []
        for label, _ in pairs:
            for N in N_list:
                resid = results[(label, N)][1][fname]
                if resid >= 1e-8:
                    not_covered.append((label, N, resid))
                    break
        if not not_covered:
            log(f"  {fname:<24} : COVERS ALL {len(pairs)} pairs at N={N_list}")
        else:
            log(f"  {fname:<24} : misses {len(not_covered)} cases:")
            for lab, N, r in not_covered[:5]:
                log(f"    {lab} at N={N}: residual {r:.3e}")

    # ============================================================
    # Phase 4: 2-qubit anti-commutation algebra for winner
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 4: 2-qubit anti-commutation check for predicted X-deph variant")
    log("=" * 88)

    M_predicted = Pi_X_5bilinear_predicted()
    Q2 = build_Q_Nsite(M_predicted, 2)
    BILINEARS = ["XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"]
    Id = np.eye(4)
    log()
    log(f"  {'Bilinear':<10} {'bit_a':<8} {'#X,#Y even (truly_X)':<22} {'{Q, [B,.]} norm':<20}")
    log("  " + "-" * 64)
    for term in BILINEARS:
        H2 = np.kron(PAULI[term[0]], PAULI[term[1]])
        C = -1j * (np.kron(H2, Id) - np.kron(Id, H2.T))
        err = np.linalg.norm(Q2 @ C + C @ Q2)
        ok = "ANTI-COMM" if err < 1e-8 else f"err={err:.2e}"
        b_tag = f"Pi2_X={'+' if bit_a(term) == 0 else '-'}"
        truly_tag = "truly_X" if is_truly_xdeph(term) else "non-truly_X"
        log(f"  {term:<10} {b_tag:<8} {truly_tag:<22} {ok}")

    # Per-site dissipator: M D[X] M^-1 = -D[X] - 2 gamma I ?
    log()
    log("  Per-site dissipator check (1-qubit):")
    sigma1 = 0.05
    D_X = sigma1 * (np.kron(sx, sx.conj()) - np.eye(4))
    Q1 = build_Q_Nsite(M_predicted, 1)
    Q1_inv = np.linalg.inv(Q1)
    conj = Q1 @ D_X @ Q1_inv
    target = -D_X - 2 * sigma1 * np.eye(4)
    err = np.linalg.norm(conj - target)
    log(f"    ||Q D[X] Q^-1 - (-D[X] - 2 gamma I)|| = {err:.4e}")
    log(f"    (Should be ~0 if X-deph variant works at dissipator level)")

    # ============================================================
    # Phase 5: Robustness
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 5: Robustness over random Pi^2_X-even non-truly H + X-deph")
    log("=" * 88)

    rng = np.random.default_rng(2026)
    M_winner = Pi_X_5bilinear_predicted()
    log()
    log(f"  {'Trial':>5} {'N':>3}  {'spec_pal':>9}  {'op_resid':>12}")
    log(f"  {'-' * 5} {'-' * 3}  {'-' * 9}  {'-' * 12}")
    for trial in range(5):
        for N in [3, 4, 5]:
            bonds = [(i, i + 1) for i in range(N - 1)]
            sigma = N * gamma
            d = 2 ** N
            H = np.zeros((d, d), dtype=complex)
            for i, j in bonds:
                for term in PI2X_EVEN_BILINEARS:
                    J = rng.normal()
                    if abs(J) < 1e-12:
                        continue
                    H += J * site_op(PAULI[term[0]], i, N) @ site_op(PAULI[term[1]], j, N)
            H = (H + H.conj().T) / 2.0
            L = build_L_xdeph(H, gamma, N)
            spec_score = spec_palindrome_score(L, sigma, tol=1e-5)
            Q = build_Q_Nsite(M_winner, N)
            resid = operator_palindrome_residual(Q, L, sigma)
            log(f"  {trial+1:>5} {N:>3}  {spec_score*100:>7.1f}%  {resid:>12.3e}")

    log()
    log(f"  Completed: {datetime.now()}")
    log("=" * 88)


if __name__ == "__main__":
    try:
        main()
    finally:
        f_log.close()
        print(f"\n>>> Results saved to: {OUT}")
