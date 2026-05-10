"""Verify the S_N-permutation theorem: S(t) for ρ_cc + single XY-pair H is
bond-position-blind, including for non-NN long-range pairs.

Theorem (Pfad C):
    For H_{p,q} = J·(X_p X_q + Y_p Y_q) with arbitrary p≠q, uniform Z-deph γ₀,
    and ρ_cc = (|S_1⟩⟨S_2| + h.c.)/2 initial state,
    S(t) = Σ_l 2|(ρ_l)_{0,1}|² is independent of (p, q).

Test (N=4): compare all C(4,2)=6 site pairs:
    NN: (0,1), (1,2), (2,3)  →  3 cases
    long-range: (0,2), (0,3), (1,3)  →  3 cases
All 6 should give bit-identical S(t). Falsifies S_N-orbit argument if not.
"""
from __future__ import annotations

import sys
from itertools import combinations

import numpy as np
from scipy.linalg import expm

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

N = 4
DIM = 1 << N
J = 0.075
GAMMA = 0.05
T_MAX = 30.0
N_T = 11  # 11 sample points 0..30


def pauli_at(P: np.ndarray, site: int, n: int) -> np.ndarray:
    """Embed single-qubit P on `site` in N-qubit Hilbert space (MSB convention)."""
    I2 = np.eye(2, dtype=complex)
    op = 1.0
    for q in range(n):
        op = np.kron(op, P if q == site else I2)
    return op


X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def build_xy_pair_h(p: int, q: int, j: float, n: int) -> np.ndarray:
    Xp, Xq = pauli_at(X, p, n), pauli_at(X, q, n)
    Yp, Yq = pauli_at(Y, p, n), pauli_at(Y, q, n)
    return j * (Xp @ Xq + Yp @ Yq)


def coherence_block_probe_rho(n: int) -> np.ndarray:
    """ρ_cc = (|S_1⟩⟨S_2| + h.c.) / 2  with S_n = symmetric Dicke."""
    d = 1 << n
    idx_p1 = [i for i in range(d) if bin(i).count("1") == 1]
    idx_p2 = [i for i in range(d) if bin(i).count("1") == 2]
    a1 = 1.0 / np.sqrt(len(idx_p1))
    a2 = 1.0 / np.sqrt(len(idx_p2))
    v = 0.5 * a1 * a2
    rho = np.zeros((d, d), dtype=complex)
    for i in idx_p1:
        for j in idx_p2:
            rho[i, j] += v
            rho[j, i] += v
    return rho


def build_liouvillian_vec(h: np.ndarray, gamma: float, n: int) -> np.ndarray:
    """Build the Liouvillian superoperator on vec(ρ) for H + uniform Z-deph γ.
    L[ρ] = -i[H, ρ] + γ Σ_l (Z_l ρ Z_l - ρ)
    Using vec convention: vec(A ρ B) = (B^T ⊗ A) vec(ρ).
    """
    d = h.shape[0]
    Id = np.eye(d, dtype=complex)
    # Hamiltonian: -i[H, ρ] = -i (H ρ - ρ H) → -i (Id ⊗ H - H^T ⊗ Id)
    L = -1j * (np.kron(Id, h) - np.kron(h.T, Id))
    # Z-deph at each site: γ (Z_l ρ Z_l - ρ)
    for l in range(n):
        Zl = pauli_at(Z, l, n)
        # Z_l ρ Z_l → (Z_l^T ⊗ Z_l) vec(ρ) ; Z_l real symmetric so Z_l^T = Z_l
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


def spatial_sum_S(rho: np.ndarray, n: int) -> float:
    """S = Σ_l 2|(ρ_l)_{0,1}|² where (ρ_l)_{0,1} = Σ_rest ρ[i_bit_l=0, i_bit_l=1 (else same)]."""
    d = rho.shape[0]
    s = 0.0
    for l in range(n):
        shift = n - 1 - l
        offdiag = 0.0 + 0.0j
        for i in range(d):
            if (i >> shift) & 1:
                continue
            j = i | (1 << shift)
            offdiag += rho[i, j]
        s += 2.0 * (offdiag.real**2 + offdiag.imag**2)
    return float(s)


def main() -> None:
    rho0 = coherence_block_probe_rho(N)
    vec0 = rho0.reshape(-1, order="F")  # column-major vec

    pairs = list(combinations(range(N), 2))  # all 6 site pairs
    print(f"\nN={N}, J={J}, γ={GAMMA}, ρ_cc initial state, T_MAX={T_MAX}")
    print(f"All {len(pairs)} 2-site XY pairs: {pairs}\n")

    # Sample at 11 time points
    times = np.linspace(0, T_MAX, N_T)
    results = {}  # (p,q) -> array of S(t)
    for (p, q) in pairs:
        H = build_xy_pair_h(p, q, J, N)
        L = build_liouvillian_vec(H, GAMMA, N)
        S_traj = np.zeros(N_T)
        for ti, t in enumerate(times):
            U = expm(L * t)
            vec_t = U @ vec0
            rho_t = vec_t.reshape(DIM, DIM, order="F")
            S_traj[ti] = spatial_sum_S(rho_t, N)
        results[(p, q)] = S_traj
        kind = "NN" if q == p + 1 else "LR"
        print(f"  ({p},{q}) [{kind}]: S(0)={S_traj[0]:.10f}  S({T_MAX})={S_traj[-1]:.6e}")

    print()
    print(f"S(0) closed-form (N-1)/N = {(N-1)/N:.10f}")
    print()

    # Pairwise diff matrix (max over time)
    print("Pairwise max|S_a(t) - S_b(t)| over all sample times:")
    print()
    header = "         " + "  ".join(f"({p},{q})" for (p, q) in pairs)
    print(header)
    for (p, q) in pairs:
        Sa = results[(p, q)]
        row = [f"({p},{q})"]
        for (p2, q2) in pairs:
            Sb = results[(p2, q2)]
            d = float(np.max(np.abs(Sa - Sb)))
            row.append(f"{d:.1e}")
        print("  " + "  ".join(row))
    print()

    # Global verdict
    Sref = results[pairs[0]]
    max_global_dev = 0.0
    for key, S in results.items():
        d = float(np.max(np.abs(S - Sref)))
        if d > max_global_dev:
            max_global_dev = d
    print(f"Maximum deviation across all {len(pairs)} pairs from reference (0,1): {max_global_dev:.3e}")
    if max_global_dev < 1e-12:
        print("VERDICT: bit-identical (within numerical noise) — Theorem verified for N={N}.".format(N=N))
    else:
        print("VERDICT: NOT identical — Theorem falsified or numerical artifact.")


if __name__ == "__main__":
    main()
