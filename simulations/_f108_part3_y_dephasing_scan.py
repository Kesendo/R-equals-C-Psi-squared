"""F108 Part 3: Pi-family scan for pure-Pi^2_Y-even non-truly Hamiltonians + Y-dephasing.

Mirrors _f108_part1_pi_family_scan.py and _f108_part2_x_dephasing_scan.py. The
conjecture: the Z-deph Pi_5bilinear's mechanism lifts to Y-deph by using the
Y-dephase-appropriate phase factors, giving Pi_5bilinear_Y with per-site action

  I -> +X,  X -> -I,  Y -> -iZ,  Z -> +iY

(same I<->X, Y<->Z permutation as canonical Y-deph Pi; sign flips on X->I and
Z->Y back-arrows mirror the Z-deph variant's flips but with Y-deph's -i phase
on the Y/Z 2-cycle instead of Z-deph's +i). Pi^2_Y-even bilinears = bit_b-even
pairs = same as Pi^2_Z-even: {XX, YY, YZ, ZY, ZZ}. F85 Y-deph truly = #Y even
AND #Z even (same as Z-deph).

Scan checks if Pi_5bilinear_Y gives EXACT operator-level palindrome
Pi * L * Pi^-1 + L + 2*sigma*I = 0 for every pure-Pi^2_Y-even non-truly pair at
N=3,4,5 under Y-dephasing.

Output: console + simulations/results/f108_part3_y_dephasing_scan.txt
"""

import os
from datetime import datetime
from itertools import combinations_with_replacement, product as iprod

import numpy as np

RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
OUT = os.path.join(RESULTS_DIR, "f108_part3_y_dephasing_scan.txt")
os.makedirs(RESULTS_DIR, exist_ok=True)
f_log = open(OUT, "w", buffering=1)


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


def bit_b(bilinear):
    """Pi^2_Y parity bit (same as Pi^2_Z): #Y + #Z mod 2."""
    return sum(1 for c in bilinear if c in "YZ") % 2


def is_truly_ydeph(bilinear):
    """F85 Y-deph truly: #Y even AND #Z even (same form as Z-deph per F107)."""
    ny = sum(1 for c in bilinear if c == "Y")
    nz = sum(1 for c in bilinear if c == "Z")
    return (ny % 2 == 0) and (nz % 2 == 0)


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


def build_L_ydeph(H, gamma, N):
    """L = -i [H, .] + sum_k gamma * D[Y_k] in vec (row-major) basis."""
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Yk = site_op(sy, k, N)
        L += gamma * (np.kron(Yk, Yk.conj()) - np.eye(d2))
    return L


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


def Pi_Y_5bilinear_predicted():
    """Predicted Y-deph 5-bilinear variant: same I<->X, Y<->Z permutation as
    canonical Y-deph Pi; sign flips on X->I and Z->Y back-arrows. Y-deph phase
    on Y/Z 2-cycle is -i (vs Z-deph's +i).
    I -> +X, X -> -I, Y -> -iZ, Z -> +iY"""
    return {"I": (1, "X"), "X": (-1, "I"), "Y": (-1j, "Z"), "Z": (1j, "Y")}


def Pi_Y_canonical():
    """Canonical Y-deph P1 per PiOperator.cs:
    I -> +X, X -> +I, Y -> -iZ, Z -> -iY"""
    return {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}


def Pi_Y_variant_v2():
    """Alternative: only X<->I sign-flipped, Y<->Z normal."""
    return {"I": (1, "X"), "X": (-1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}


def Pi_Y_variant_v3():
    """Alternative: only Y<->Z sign-flipped, X<->I normal."""
    return {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (1j, "Y")}


def Pi_Y_variant_v4():
    """Alternative: opposite Y/Z sign flip."""
    return {"I": (1, "X"), "X": (-1, "I"), "Y": (1j, "Z"), "Z": (-1j, "Y")}


def spec_palindrome_score(L, sigma, tol=1e-5):
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
    try:
        Q_inv = np.linalg.inv(Q)
    except np.linalg.LinAlgError:
        return np.inf
    d2 = L.shape[0]
    M = Q @ L @ Q_inv + L + 2.0 * sigma * np.eye(d2)
    return np.linalg.norm(M)


# Pi^2_Y-even = bit_b=0 = same set as Pi^2_Z-even (excluding I-containing)
PI2Y_EVEN_BILINEARS = ["XX", "YY", "YZ", "ZY", "ZZ"]
TRULY_BILINEARS_YDEPH = ["XX", "YY", "ZZ"]
NONTRULY_PI2YEVEN_BILINEARS = ["YZ", "ZY"]


def pair_is_pure_pi2y_even_non_truly(t1, t2):
    if bit_b(t1) != 0 or bit_b(t2) != 0:
        return False
    if is_truly_ydeph(t1) and is_truly_ydeph(t2):
        return False
    return True


def main():
    log("=" * 88)
    log("F108 Part 3: Pi-family scan for pure-Pi^2_Y-even non-truly H + Y-dephasing")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    log("\nPi^2_Y-even (bit_b=0) 2-site bilinears: {XX, YY, YZ, ZY, ZZ} (same set as Pi^2_Z-even)")
    log("Truly under Y-deph (F85): #Y even AND #Z even -> {XX, YY, ZZ}")
    log("Pure-Pi^2_Y-even NON-truly bilinears: {YZ, ZY}")

    log("\n" + "=" * 88)
    log("Phase 1: Enumerate pure-Pi^2_Y-even non-truly pairs")
    log("=" * 88)

    pairs = []
    for b in NONTRULY_PI2YEVEN_BILINEARS:
        pairs.append((b, {b: 1.0}))
    for t1, t2 in combinations_with_replacement(PI2Y_EVEN_BILINEARS, 2):
        if t1 == t2:
            continue
        if pair_is_pure_pi2y_even_non_truly(t1, t2):
            label = f"{t1}+{t2}"
            pairs.append((label, {t1: 1.0, t2: 1.0}))

    log(f"\n  Total pairs: {len(pairs)}")
    for label, _ in pairs:
        log(f"    {label}")

    log("\n" + "=" * 88)
    log("Phase 2: Pi-family operator-level palindrome scan (Y-dephasing)")
    log("=" * 88)

    gamma = 0.05
    N_list = [3, 4, 5]

    families = {
        "Y_canonical": Pi_Y_canonical(),
        "Y_5bilinear_predicted": Pi_Y_5bilinear_predicted(),
        "Y_variant_v2": Pi_Y_variant_v2(),
        "Y_variant_v3": Pi_Y_variant_v3(),
        "Y_variant_v4": Pi_Y_variant_v4(),
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
            L = build_L_ydeph(H, gamma, N)
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

    log("\n" + "=" * 88)
    log(f"Phase 3: Single-family coverage test (Y-dephasing)")
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

    log("\n" + "=" * 88)
    log("Phase 4: 2-qubit anti-commutation check for predicted Y-deph variant")
    log("=" * 88)

    M_predicted = Pi_Y_5bilinear_predicted()
    Q2 = build_Q_Nsite(M_predicted, 2)
    BILINEARS = ["XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"]
    Id = np.eye(4)
    log()
    log(f"  {'Bilinear':<10} {'bit_b':<8} {'#Y,#Z even (truly_Y)':<22} {'{Q, [B,.]} norm':<20}")
    log("  " + "-" * 64)
    for term in BILINEARS:
        H2 = np.kron(PAULI[term[0]], PAULI[term[1]])
        C = -1j * (np.kron(H2, Id) - np.kron(Id, H2.T))
        err = np.linalg.norm(Q2 @ C + C @ Q2)
        ok = "ANTI-COMM" if err < 1e-8 else f"err={err:.2e}"
        b_tag = f"Pi2_Y={'+' if bit_b(term) == 0 else '-'}"
        truly_tag = "truly_Y" if is_truly_ydeph(term) else "non-truly_Y"
        log(f"  {term:<10} {b_tag:<8} {truly_tag:<22} {ok}")

    log()
    log("  Per-site dissipator check (1-qubit):")
    sigma1 = 0.05
    D_Y = sigma1 * (np.kron(sy, sy.conj()) - np.eye(4))
    Q1 = build_Q_Nsite(M_predicted, 1)
    Q1_inv = np.linalg.inv(Q1)
    conj = Q1 @ D_Y @ Q1_inv
    target = -D_Y - 2 * sigma1 * np.eye(4)
    err = np.linalg.norm(conj - target)
    log(f"    ||Q D[Y] Q^-1 - (-D[Y] - 2 gamma I)|| = {err:.4e}")
    log(f"    (Should be ~0 if Y-deph variant works at dissipator level)")

    log("\n" + "=" * 88)
    log("Phase 5: Robustness over random Pi^2_Y-even non-truly H + Y-deph")
    log("=" * 88)

    rng = np.random.default_rng(2026)
    M_winner = Pi_Y_5bilinear_predicted()
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
                for term in PI2Y_EVEN_BILINEARS:
                    J = rng.normal()
                    if abs(J) < 1e-12:
                        continue
                    H += J * site_op(PAULI[term[0]], i, N) @ site_op(PAULI[term[1]], j, N)
            H = (H + H.conj().T) / 2.0
            L = build_L_ydeph(H, gamma, N)
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
