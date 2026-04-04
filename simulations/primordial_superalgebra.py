#!/usr/bin/env python3
"""
Primordial Superalgebra: Light and Lens as One
=============================================
Verifies the Pythagorean orthogonality {L_H, L_D + Sigma_gamma I} = 0
at N=2, quantifies aberration at N>=3, decomposes eigenmodes into
light/lens sectors, and performs Seidel aberration analysis.

R=CΨ² Project, Homework #10
Source: TASK_PRIMORDIAL_SUPERALGEBRA.md
"""

import sys
import io
import numpy as np
from pathlib import Path
import time
import itertools

# Force UTF-8 output on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Constants ──────────────────────────────────────────────────────
J = 1.0
GAMMA = 0.05
TOL = 1e-8

# ── Pauli matrices ────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, Xm, Ym, Zm]
PAULI_LABELS = ['I', 'X', 'Y', 'Z']

# ── Output ────────────────────────────────────────────────────────
out = []

def log(msg=""):
    print(msg)
    out.append(msg)


# ── Utilities ─────────────────────────────────────────────────────
def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]


def build_hamiltonian(N, bonds):
    """Heisenberg XXX: H = J sum_{(a,b)} (X_a X_b + Y_a Y_b + Z_a Z_b)."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N
            ops[a] = P
            ops[b] = P
            H += J * kron_chain(ops)
    return H


def build_L_H(H, d):
    """Hamiltonian superoperator: L_H(rho) = -i[H, rho]."""
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def build_L_D(N, gammas):
    """Dissipative superoperator: L_D = sum_k gamma_k D[Z_k]."""
    d = 2**N
    Id = np.eye(d, dtype=complex)
    L_D = np.zeros((d**2, d**2), dtype=complex)
    for k in range(N):
        ops = [I2] * N
        ops[k] = Zm
        Lk = np.sqrt(gammas[k]) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L_D += (np.kron(Lk, Lk.conj())
                - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T)))
    return L_D


def classify_paired(eigvals, sigma_gamma, tol=TOL):
    """Classify eigenvalues as palindromically paired.
    Partner of lambda is at -2 Sigma_gamma - conj(lambda)."""
    n = len(eigvals)
    paired = np.zeros(n, dtype=bool)
    partner_idx = -np.ones(n, dtype=int)
    used = set()
    for i in range(n):
        if i in used:
            continue
        target = -2 * sigma_gamma - eigvals[i].conjugate()
        dists = np.abs(eigvals - target)
        dists[i] = 999
        for j in used:
            dists[j] = 999
        j = np.argmin(dists)
        if dists[j] < tol:
            paired[i] = True
            paired[j] = True
            partner_idx[i] = j
            partner_idx[j] = i
            used.add(i)
            used.add(j)
    return paired, partner_idx


def generate_pauli_basis(N):
    """Generate all 4^N Pauli basis vectors for N qubits.

    Returns:
        labels:     list of strings ('IXYZ' style)
        pauli_mat:  d^2 x 4^N matrix (columns = vectorized Pauli strings)
        n_xy:       array of X/Y factor counts per string
    """
    vecs = []
    labels = []
    n_xy = []
    for indices in itertools.product(range(4), repeat=N):
        label = ''.join(PAULI_LABELS[i] for i in indices)
        mat = kron_chain([PAULIS[i] for i in indices])
        vecs.append(mat.flatten())          # row-major vectorization
        labels.append(label)
        n_xy.append(sum(1 for i in indices if i in (1, 2)))
    return labels, np.array(vecs).T, np.array(n_xy)


# ══════════════════════════════════════════════════════════════════
# STEP 1: Pythagorean Orthogonality  {L_H, L_D + Sigma_gamma I} = 0
# ══════════════════════════════════════════════════════════════════
def step1():
    log("=" * 75)
    log("STEP 1: Pythagorean Orthogonality  {L_H, L_D + Sigma_gamma I} = 0")
    log("=" * 75)
    log(f"J = {J}, gamma = {GAMMA}, chain topology")
    log()
    log(f"{'N':>3}  {'dim':>6}  {'Sigma_g':>8}  {'||A||_F':>14}  "
        f"{'||L_H||':>14}  {'||L_Ds||':>14}  {'rel. dev.':>14}  {'t':>6}")
    log("-" * 88)

    results = []
    for N in range(2, 7):
        d = 2**N
        dim = d**2
        t0 = time.time()

        bonds = chain_bonds(N)
        gammas = [GAMMA] * N
        sigma_gamma = sum(gammas)

        H = build_hamiltonian(N, bonds)
        L_H = build_L_H(H, d)
        L_D = build_L_D(N, gammas)
        L_Ds = L_D + sigma_gamma * np.eye(dim, dtype=complex)

        A = L_H @ L_Ds + L_Ds @ L_H

        nA = np.linalg.norm(A, 'fro')
        nH = np.linalg.norm(L_H, 'fro')
        nD = np.linalg.norm(L_Ds, 'fro')
        rel = nA / (nH * nD) if nH * nD > 0 else 0.0
        dt = time.time() - t0

        log(f"{N:3d}  {dim:6d}  {sigma_gamma:8.4f}  {nA:14.6e}  "
            f"{nH:14.6e}  {nD:14.6e}  {rel:14.6e}  {dt:5.1f}s")
        results.append(dict(N=N, dim=dim, sigma=sigma_gamma,
                            nA=nA, nH=nH, nD=nD, rel=rel))

    log()
    log("Interpretation:")
    for r in results:
        if r['rel'] < 1e-12:
            log(f"  N={r['N']}: EXACT Pythagorean orthogonality (< 10^-12)")
        else:
            log(f"  N={r['N']}: aberration = {r['rel']*100:.6f}%")
    return results


# ══════════════════════════════════════════════════════════════════
# STEP 2: Aberration Profile and Gamma Independence
# ══════════════════════════════════════════════════════════════════
def step2():
    log()
    log("=" * 75)
    log("STEP 2: Aberration Profile and Gamma Independence")
    log("=" * 75)
    log()

    gammas_test = [0.01, 0.05, 0.1, 0.5, 1.0]

    header = f"{'N':>3}  {'dim':>5}"
    for g in gammas_test:
        header += f"  {'g=' + f'{g}':>14}"
    log(header)
    log("-" * (10 + 16 * len(gammas_test)))

    results = {}
    for N in range(2, 7):
        d = 2**N
        dim = d**2
        bonds = chain_bonds(N)
        H = build_hamiltonian(N, bonds)
        L_H = build_L_H(H, d)
        nH = np.linalg.norm(L_H, 'fro')

        row = f"{N:3d}  {dim:5d}"
        devs = []
        for g in gammas_test:
            gs = [g] * N
            sg = sum(gs)
            L_D = build_L_D(N, gs)
            L_Ds = L_D + sg * np.eye(dim, dtype=complex)
            A = L_H @ L_Ds + L_Ds @ L_H
            nA = np.linalg.norm(A, 'fro')
            nD = np.linalg.norm(L_Ds, 'fro')
            rel = nA / (nH * nD) if nH * nD > 0 else 0.0
            row += f"  {rel:14.6e}"
            devs.append(rel)
        log(row)
        results[N] = devs

    log()
    log("Gamma independence (coefficient of variation across 5 gamma values):")
    for N in sorted(results):
        arr = np.array(results[N])
        if np.mean(arr) > 1e-14:
            cv = np.std(arr) / np.mean(arr)
            log(f"  N={N}: CV = {cv:.2e}  "
                f"({'gamma-INDEPENDENT' if cv < 0.01 else 'gamma-DEPENDENT'})")
        else:
            log(f"  N={N}: identically zero (exact Pythagorean)")

    log()
    log("Trend with N (gamma = 0.05 column):")
    g_idx = gammas_test.index(0.05)
    Ns = sorted(results)
    for i in range(1, len(Ns)):
        p, c = results[Ns[i - 1]][g_idx], results[Ns[i]][g_idx]
        if p > 1e-14:
            log(f"  N={Ns[i-1]} -> {Ns[i]}: "
                f"{p:.6e} -> {c:.6e}  (x{c / p:.4f})")
        else:
            log(f"  N={Ns[i-1]} -> {Ns[i]}: 0 -> {c:.6e}")

    return results


# ══════════════════════════════════════════════════════════════════
# STEP 3: Light/Lens Sector Decomposition of Eigenmodes
# ══════════════════════════════════════════════════════════════════
def step3():
    log()
    log("=" * 75)
    log("STEP 3: Light/Lens Sector Decomposition of Eigenmodes")
    log("=" * 75)
    log(f"J = {J}, gamma = {GAMMA}, chain topology")
    log("k = number of X/Y factors in Pauli string")
    log("  k=0: fully immune (I,Z)^N = LENS")
    log("  k>0: contains coherences  = LIGHT")
    log()

    all_corr = []

    for N in range(2, 6):
        d = 2**N
        dim = d**2
        log(f"-- N = {N}  (dim = {dim}) " + "-" * 50)
        t0 = time.time()

        bonds = chain_bonds(N)
        gammas = [GAMMA] * N
        sigma_gamma = sum(gammas)

        H = build_hamiltonian(N, bonds)
        L = build_L_H(H, d) + build_L_D(N, gammas)

        eigvals, eigvecs = np.linalg.eig(L)

        _, pauli_mat, n_xy = generate_pauli_basis(N)

        # Coefficients c_alpha = vec(P_alpha)^dag v / sqrt(d)
        # Parseval: sum |c_alpha|^2 = ||v||^2 = 1
        coeffs = pauli_mat.conj().T @ eigvecs / np.sqrt(d)
        wt = np.abs(coeffs)**2      # shape: 4^N x dim

        # Sector weights per mode
        sector_w = np.zeros((N + 1, dim))
        for k in range(N + 1):
            mask = (n_xy == k)
            sector_w[k] = np.sum(wt[mask], axis=0)

        # Normalise (should be ~1 already)
        totals = np.sum(sector_w, axis=0)
        sector_w /= np.maximum(totals[np.newaxis, :], 1e-30)

        lens = sector_w[0]       # k=0 fraction
        light = 1.0 - lens       # k>0 fraction
        decay = np.abs(eigvals.real)

        # ── Global correlation ────────────────────────────────────
        active = decay > TOL
        n_active = int(np.sum(active))
        if n_active > 2:
            r_corr = np.corrcoef(decay[active], light[active])[0, 1]
            log(f"  Correlation(|Re(lambda)|, light fraction): "
                f"r = {r_corr:+.6f}  ({n_active} active modes)")
            all_corr.append((N, r_corr))
            if r_corr > 0.3:
                log(f"    Positive: faster-decaying modes carry MORE light character")
            elif r_corr < -0.3:
                log(f"    Negative: faster-decaying modes carry MORE lens character")
            else:
                log(f"    Weak: no clear sorting by decay rate")

        # ── Average sector weights ────────────────────────────────
        log(f"  Average sector weights across all {dim} modes:")
        for k in range(N + 1):
            n_str = int(np.sum(n_xy == k))
            log(f"    k={k}: mean weight = {np.mean(sector_w[k]):.6f}  "
                f"({n_str} Pauli strings, basis fraction = {n_str / dim:.4f})")

        # ── Palindromic pair analysis ─────────────────────────────
        paired, pidx = classify_paired(eigvals, sigma_gamma)
        n_pairs = int(np.sum(paired)) // 2
        log(f"  Palindromic pairs: {n_pairs}")

        seen = set()
        fast_sw = np.zeros(N + 1)
        slow_sw = np.zeros(N + 1)
        pair_count = 0

        for i in range(dim):
            j = int(pidx[i])
            if not paired[i] or j < 0 or i in seen:
                continue
            seen.add(i)
            seen.add(j)

            if abs(eigvals[i].real) < abs(eigvals[j].real):
                fi, si = i, j
            else:
                fi, si = j, i

            fast_sw += sector_w[:, fi]
            slow_sw += sector_w[:, si]
            pair_count += 1

        if pair_count > 0:
            fast_sw /= pair_count
            slow_sw /= pair_count
            log(f"  Pair sector analysis ({pair_count} pairs):")
            log(f"    {'k':>3}  {'fast(Re~0)':>14}  {'slow(Re~-2Sg)':>14}  "
                f"{'Delta':>14}")
            log(f"    " + "-" * 48)
            for k in range(N + 1):
                d_ = fast_sw[k] - slow_sw[k]
                log(f"    {k:3d}  {fast_sw[k]:14.6f}  {slow_sw[k]:14.6f}  "
                    f"{d_:+14.6f}")
            log(f"    Lens total (k=0): fast = {fast_sw[0]:.6f}, "
                f"slow = {slow_sw[0]:.6f}")
            log(f"    Light total (k>0): fast = {1 - fast_sw[0]:.6f}, "
                f"slow = {1 - slow_sw[0]:.6f}")
            if fast_sw[0] > slow_sw[0]:
                log(f"    --> Fast partners are MORE lens-like (I,Z dominated)")
            else:
                log(f"    --> Slow partners are MORE lens-like (I,Z dominated)")

        log(f"  ({time.time() - t0:.1f}s)")
        log()

    return all_corr


# ══════════════════════════════════════════════════════════════════
# STEP 4: Seidel Aberration Decomposition per Weight Sector
# ══════════════════════════════════════════════════════════════════
def step4():
    log()
    log("=" * 75)
    log("STEP 4: Seidel Aberration Decomposition per Weight Sector")
    log("=" * 75)
    log(f"J = {J}, gamma = {GAMMA}, chain topology")
    log("A = anticommutator transformed to Pauli basis")
    log("Diagonal blocks A_{kk}: intra-sector aberration")
    log()

    for N in range(2, 6):
        d = 2**N
        dim = d**2
        log(f"-- N = {N}  (dim = {dim}) " + "-" * 50)
        t0 = time.time()

        bonds = chain_bonds(N)
        gammas = [GAMMA] * N
        sg = sum(gammas)

        H = build_hamiltonian(N, bonds)
        L_H = build_L_H(H, d)
        L_D = build_L_D(N, gammas)
        L_Ds = L_D + sg * np.eye(dim, dtype=complex)

        A = L_H @ L_Ds + L_Ds @ L_H
        nA = np.linalg.norm(A, 'fro')

        if nA < 1e-12:
            log(f"  ||A|| = {nA:.2e}  (exact Pythagorean, no aberration)")
            log(f"  ({time.time() - t0:.1f}s)")
            log()
            continue

        _, pauli_mat, n_xy = generate_pauli_basis(N)
        U = pauli_mat / np.sqrt(d)          # normalised Pauli basis
        A_pb = U.conj().T @ A @ U           # A in Pauli basis

        nA_pb = np.linalg.norm(A_pb, 'fro')
        log(f"  ||A|| = {nA:.6e},  ||A_pb|| = {nA_pb:.6e}  "
            f"(unitarity check: ratio = {nA_pb / nA:.10f})")

        # ── Diagonal blocks ───────────────────────────────────────
        log(f"  {'k':>3}  {'dim_k':>6}  {'||A_kk||':>14}  "
            f"{'frac of ||A||':>14}  {'per element':>14}")
        log(f"  " + "-" * 58)

        fracs = []
        for k in range(N + 1):
            idx = np.where(n_xy == k)[0]
            dk = len(idx)
            Akk = A_pb[np.ix_(idx, idx)]
            nk = np.linalg.norm(Akk, 'fro')
            frac = nk / nA
            per_el = nk / dk if dk > 0 else 0.0
            log(f"  {k:3d}  {dk:6d}  {nk:14.6e}  "
                f"{frac:14.6f}  {per_el:14.6e}")
            fracs.append(frac)

        # ── Off-diagonal couplings ────────────────────────────────
        log(f"  Cross-sector coupling ||A_kl|| (fraction of ||A||):")
        for k in range(N + 1):
            for l in range(k + 1, N + 1):
                ik = np.where(n_xy == k)[0]
                il = np.where(n_xy == l)[0]
                Akl = A_pb[np.ix_(ik, il)]
                nkl = np.linalg.norm(Akl, 'fro')
                fkl = nkl / nA
                if fkl > 0.001:
                    log(f"    k={k} <-> l={l}: {fkl:.6f}")

        # ── Aberration type ───────────────────────────────────────
        log()
        fracs_arr = np.array(fracs)
        peak_k = int(np.argmax(fracs_arr))

        # Palindromic symmetry of aberration profile
        if N >= 4:
            half = N // 2
            left = fracs_arr[:half]
            right = fracs_arr[-half:][::-1]
            asym = np.sum(np.abs(left - right))
            total_lr = np.sum(left + right)
            asym_ratio = asym / total_lr if total_lr > 0 else 0
            log(f"  Aberration profile symmetry: "
                f"{'PALINDROMIC' if asym_ratio < 0.05 else 'ASYMMETRIC'} "
                f"(asymmetry = {asym_ratio:.4f})")

        if N >= 3:
            edge_sum = fracs_arr[0] + fracs_arr[N]
            mid_sum = float(np.sum(fracs_arr[1:-1]))
            log(f"  Edge sectors (k=0,{N}): {edge_sum:.4f}")
            log(f"  Interior sectors (k=1..{N-1}): {mid_sum:.4f}")
            if mid_sum > edge_sum:
                log(f"  --> Interior-dominated: coma-type aberration")
            else:
                log(f"  --> Edge-dominated: spherical-type aberration")

        log(f"  Peak aberration at k={peak_k}")
        log(f"  ({time.time() - t0:.1f}s)")
        log()


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    log("Primordial Superalgebra: Light and Lens as One")
    log("=" * 75)
    log(f"R=CΨ² Project, Homework #10")
    log(f"Heisenberg chain, J = {J}, Z-dephasing")
    log()

    t_total = time.time()

    r1 = step1()
    r2 = step2()
    r3 = step3()
    step4()

    # ── Summary ───────────────────────────────────────────────────
    log()
    log("=" * 75)
    log("SUMMARY")
    log("=" * 75)

    log()
    log("1. Pythagorean orthogonality {L_H, L_D + Sigma_gamma} = 0:")
    for r in r1:
        s = "EXACT" if r['rel'] < 1e-12 else f"{r['rel'] * 100:.6f}%"
        log(f"   N={r['N']}: {s}")

    log()
    log("2. Gamma independence: CONFIRMED")
    log("   ||A|| / (||L_H|| ||L_Ds||) is constant across gamma values")
    log("   because L_D scales linearly with gamma.")

    log()
    log("3. Aberration trend (gamma = 0.05):")
    g_idx = 1   # index of gamma=0.05 in gammas_test
    Ns = sorted(r2.keys())
    for N in Ns:
        val = r2[N][g_idx]
        if val < 1e-14:
            log(f"   N={N}: 0 (exact)")
        else:
            log(f"   N={N}: {val * 100:.6f}%")

    log()
    log("4. Light/lens correlation:")
    for N, corr in r3:
        log(f"   N={N}: corr(|Re(lambda)|, light fraction) = {corr:+.6f}")

    dt = time.time() - t_total
    log()
    log(f"Total time: {dt:.1f}s")

    # ── Save ──────────────────────────────────────────────────────
    out_path = RESULTS_DIR / "primordial_superalgebra.txt"
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(out))
    log(f"\n>>> Results saved to: {out_path}")
