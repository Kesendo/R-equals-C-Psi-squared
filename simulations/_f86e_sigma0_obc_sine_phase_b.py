"""F86e σ_0 closed-form Phase B: V_inter in OBC sine basis, bond-pair structure.

Phase A established: σ_0 reproduced for N=5..15, Aitken² limit 2.86287 ± 2e-5,
all previously-tested closed forms (2√2, √(41/5), √(8+π/16), √(2π+2)) ruled out.

Phase B goal: expose the structural composition of the top singular vector
U_top of V_inter in two complementary bases:

  Basis 1 (bond-pair, direction (b) from C2InterChannelAnalytical obstruction
  note): label HD=1 entries by (a, b) where a = the single site of the
  popcount-1 state, b = the extra site in the popcount-2 state. U_top
  becomes an (N × (N-1)) matrix; if the structure decomposes as a sum of
  rank-1 sine-mode products, that's the closed form.

  Basis 2 (sine-mode, full Slater basis): transform V_inter to (k; k_1, k_2)
  index where k = single-particle sine mode (1..N), (k_1, k_2) = ordered
  2-particle Slater modes. In this basis M_H_total is diagonal with
  Δ = E_k − E_{k_1} − E_{k_2}; the projector structure becomes the carrier
  of σ_0.

Diagnostic outputs at N=5, 6, 7, 8, 10:
  - U_top as (a, b) bond-pair matrix
  - U_top in sine basis: where on the k-axis is its support?
  - SVD rank of V_inter (single vs double-degeneracy from chain mirror R)
  - σ_0² as candidate closed-form sum-over-modes

Output: simulations/results/f86e_sigma0_obc_sine_phase_b/
"""
from __future__ import annotations

import math
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")
import framework as fw  # noqa: E402


def obc_sine_modes(N: int) -> np.ndarray:
    """OBC sine basis: ψ[k-1, i] = √(2/(N+1)) · sin(πk(i+1)/(N+1)) for k=1..N, i=0..N-1."""
    psi = np.zeros((N, N))
    for k in range(1, N + 1):
        for i in range(N):
            psi[k - 1, i] = math.sqrt(2.0 / (N + 1)) * math.sin(
                math.pi * k * (i + 1) / (N + 1))
    return psi


def obc_dispersion(N: int, J: float = 1.0) -> np.ndarray:
    """E_k = 2J cos(πk/(N+1)) for k=1..N."""
    return np.array([2 * J * math.cos(math.pi * k / (N + 1)) for k in range(1, N + 1)])


def sites_of_integer(integer: int, N: int) -> list[int]:
    """Bit positions set in integer, sorted ascending. Convention: bit b = site b."""
    return [b for b in range(N) if (integer >> b) & 1]


def hd_subspace_projector(N: int, n: int, hd_value: int) -> np.ndarray:
    """Orthonormal projector P (M × n_HD) onto HD-subspace of (n, n+1) block."""
    P_n = fw.popcount_states(N, n)
    P_np1 = fw.popcount_states(N, n + 1)
    Mnp1 = len(P_np1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}

    cols = []
    labels = []
    for p in P_n:
        for q in P_np1:
            if bin(p ^ q).count("1") == hd_value:
                idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
                v = np.zeros(len(P_n) * Mnp1, dtype=complex)
                v[idx] = 1.0
                cols.append(v)
                labels.append((p, q))
    return (np.column_stack(cols) if cols else np.zeros(
        (len(P_n) * Mnp1, 0), dtype=complex)), labels


def sigma_0_with_basis(N: int, n: int = 1, gamma_0: float = 0.05):
    """Compute σ_0 with full basis labels for U_top, V_top."""
    _, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H_total = sum(M_H_per_bond)
    P_HD1, hd1_labels = hd_subspace_projector(N, n, 1)
    P_HD3, hd3_labels = hd_subspace_projector(N, n, 3)
    V_inter = P_HD1.conj().T @ M_H_total @ P_HD3
    U, S, Vh = np.linalg.svd(V_inter, full_matrices=False)
    return {
        "N": N,
        "sigma_0": float(S[0]),
        "sigma_all": S,
        "U_top": U[:, 0],
        "V_top": Vh[0, :].conj(),
        "hd1_labels": hd1_labels,
        "hd3_labels": hd3_labels,
        "V_inter": V_inter,
        "degeneracy_top": int(np.sum(np.abs(S - S[0]) < 1e-10)),
    }


def u_top_as_bond_pair_matrix(result: dict, N: int) -> np.ndarray:
    """Reshape U_top (HD=1 subspace, c=2 n=1) into (a, b)-indexed N × (N-1) matrix.

    HD=1 entry (p, q) has p with one bit at site a, q = p ∪ {site b ≠ a}.
    Matrix M[a, idx_of_b] = U_top entry for (p={a}, q={a, b}).
    For each a, the (N−1) remaining sites b are enumerated in ascending order.
    """
    M = np.zeros((N, N - 1), dtype=complex)
    for (p, q), amp in zip(result["hd1_labels"], result["U_top"]):
        p_sites = sites_of_integer(p, N)
        q_sites = sites_of_integer(q, N)
        assert len(p_sites) == 1 and len(q_sites) == 2
        a = p_sites[0]
        b_candidates = [s for s in q_sites if s != a]
        assert len(b_candidates) == 1
        b = b_candidates[0]
        # idx_of_b among the (N-1) sites ≠ a, in ascending order
        b_list = sorted([s for s in range(N) if s != a])
        b_idx = b_list.index(b)
        M[a, b_idx] = amp
    return M


def slater_unitary_p2(N: int, psi: np.ndarray):
    """Build U_p2: (M_p2 × M_p2) Slater-determinant transformation.

    Columns labeled by mode pairs (k_1, k_2), k_1 < k_2, k=1..N.
    Rows labeled by popcount-2 computational states (in fw.popcount_states ordering),
    each interpreted as a fermionic 2-site occupation.
    """
    p2_states = fw.popcount_states(N, 2)
    M_p2 = len(p2_states)
    p2_sites = [sites_of_integer(i, N) for i in p2_states]
    mode_pairs = [(k1, k2) for k1 in range(1, N + 1) for k2 in range(k1 + 1, N + 1)]
    assert len(mode_pairs) == M_p2, f"Mode pair count mismatch at N={N}: {len(mode_pairs)} vs {M_p2}"

    U = np.zeros((M_p2, M_p2), dtype=complex)
    for r, sites in enumerate(p2_sites):
        i1, i2 = sites
        for c, (k1, k2) in enumerate(mode_pairs):
            slater = (psi[k1 - 1, i1] * psi[k2 - 1, i2]
                      - psi[k1 - 1, i2] * psi[k2 - 1, i1]) / math.sqrt(2.0)
            U[r, c] = slater
    return U, mode_pairs


def transform_v_inter_to_sine_basis(N: int, result: dict, gamma_0: float = 0.05):
    """Lift V_inter back to M × M and transform to sine basis.

    Returns V_sine, U_top_sine, V_top_sine where the _sine vectors are
    expressed in the sine operator basis |k⟩⟨k_1 k_2|, indexed by (k; k_1, k_2).
    """
    n = 1
    P_n = fw.popcount_states(N, n)
    P_np1 = fw.popcount_states(N, n + 1)
    M_p1 = len(P_n)
    M_p2 = len(P_np1)
    M = M_p1 * M_p2

    # Lift V_inter to full M × M operator basis
    P_HD1, _ = hd_subspace_projector(N, n, 1)
    P_HD3, _ = hd_subspace_projector(N, n, 3)
    V_full = P_HD1 @ result["V_inter"] @ P_HD3.conj().T

    # Build sine transformations
    psi = obc_sine_modes(N)
    U_p1 = psi.T  # (N, N): rows = computational sites (i), cols = sine modes (k)
    U_p2, mode_pairs = slater_unitary_p2(N, psi)  # (M_p2, M_p2): rows = comp, cols = sine pairs

    # Operator basis transformation: vec(|k⟩⟨k1k2|) in comp basis
    #   = sum_{i,j} ψ_k(i) · U_p2[j, (k1,k2)]* · vec(|i⟩⟨j|)
    # So T_op[(i,j), (k,k1k2)] = ψ_k(i) · U_p2[j, (k1,k2)]*
    # vec(V_sine) = T_op^† @ vec(V_full @ T_op) ... mess.
    # Simpler: V_sine[(k, k1k2), (k', k1'k2')] = <k|<k1k2| V_full[as 4-tensor] |k'> |k1'k2'>
    # Convert V_full from (M, M) → (M_p1, M_p2, M_p1, M_p2) → contract with U_p1, U_p2

    V_full_4 = V_full.reshape(M_p1, M_p2, M_p1, M_p2)
    # Transform first index: V[k, j_o, i, j] = sum_i U_p1[i, k]^* · V[i, j_o, i, j]
    V_k = np.einsum('ia,iJjK->aJjK', U_p1.conj(), V_full_4)
    V_kK = np.einsum('aJjK,JM->aMjK', V_k, U_p2.conj())
    V_kKi = np.einsum('aMjK,jb->aMbK', V_kK, U_p1)
    V_sine_4 = np.einsum('aMbK,KN->aMbN', V_kKi, U_p2)
    # V_sine_4[a, M, b, N] = <a|<M| V_sine |b>|N> where a, b are single-particle modes
    # and M, N are 2-particle Slater mode-pair indices.
    V_sine = V_sine_4.reshape(M_p1 * M_p2, M_p1 * M_p2)

    # SVD of V_sine should give same singular values as V_full (basis change is unitary)
    U_s, S_s, Vh_s = np.linalg.svd(V_sine, full_matrices=False)

    return V_sine_4, V_sine, U_s[:, 0], Vh_s[0, :].conj(), S_s, mode_pairs


def main() -> None:
    gamma_0 = 0.05
    n = 1

    print("=" * 100)
    print("F86e Phase B: V_inter in OBC sine basis, bond-pair structure, sweet-spot N=7 anchor")
    print("=" * 100)

    out_dir = Path("simulations/results/f86e_sigma0_obc_sine_phase_b")
    out_dir.mkdir(parents=True, exist_ok=True)

    for N in [5, 6, 7, 8, 10]:
        print()
        print(f"--- N = {N} ---")
        t0 = time.time()
        result = sigma_0_with_basis(N, n=n, gamma_0=gamma_0)
        elapsed = time.time() - t0
        print(f"  σ_0 = {result['sigma_0']:.10f}  "
              f"(degen={result['degeneracy_top']}, "
              f"V_inter shape ({len(result['hd1_labels'])}, {len(result['hd3_labels'])}), "
              f"{elapsed:.2f}s)")

        if N == 7:
            print(f"  σ_0 − 2√2 = {result['sigma_0'] - 2 * math.sqrt(2):+.2e}")

        M_bond = u_top_as_bond_pair_matrix(result, N)
        print(f"  U_top in (a, b) bond-pair basis (rows = site a, cols = b-index 0..N−2):")
        for a in range(N):
            row = [f"{M_bond[a, j].real:+.4f}" for j in range(N - 1)]
            print(f"    a={a}: [{' '.join(row)}]")
        print(f"  |U_top|² = {np.sum(np.abs(M_bond) ** 2):.6f} (should be 1.0)")

        rank_bond = np.linalg.matrix_rank(M_bond, tol=1e-10)
        u_bond, s_bond, vh_bond = np.linalg.svd(M_bond, full_matrices=False)
        print(f"  Bond-pair matrix SVD: rank={rank_bond}, singular values = "
              f"{[f'{s:.4f}' for s in s_bond[:min(5, len(s_bond))]]}")

        # Phase B core: sine-basis transform
        V_sine_4, V_sine, U_top_sine, V_top_sine, S_sine, mode_pairs = \
            transform_v_inter_to_sine_basis(N, result, gamma_0)
        print(f"  V_sine top-3 singular values: {[f'{s:.6f}' for s in S_sine[:3]]}")
        print(f"  (should match: σ_0 = {result['sigma_0']:.6f})")

        # Show |U_top_sine|² distribution over (k, k_1, k_2)
        M_p2 = len(mode_pairs)
        U_top_sine_4 = U_top_sine.reshape(N, M_p2)  # rows: single-particle k (1..N), cols: pair (k_1, k_2)
        prob_k = np.sum(np.abs(U_top_sine_4) ** 2, axis=1)  # marginal over single-particle k
        prob_pair = np.sum(np.abs(U_top_sine_4) ** 2, axis=0)  # marginal over Slater pair

        print(f"  U_top single-particle marginal |c_k|² (k=1..N):")
        print(f"    {' '.join([f'k={k}:{prob_k[k-1]:.4f}' for k in range(1, N+1)])}")

        # Top contributing pairs
        top_pairs_idx = np.argsort(prob_pair)[::-1][:5]
        print(f"  U_top top-5 Slater pair contributions:")
        for idx in top_pairs_idx:
            k1, k2 = mode_pairs[idx]
            print(f"    (k_1={k1}, k_2={k2}): |c|² = {prob_pair[idx]:.4f}")

        E_k = obc_dispersion(N, J=1.0)
        # Compute Δ_{k; k_1, k_2} = E_k - E_{k_1} - E_{k_2} for the top pair
        if len(top_pairs_idx) > 0:
            kt1, kt2 = mode_pairs[top_pairs_idx[0]]
            print(f"  E_k = {[f'{E_k[k-1]:+.4f}' for k in range(1, N+1)]}")
            print(f"  Top pair Δ = E_k − E_{{k1}} − E_{{k2}} for k=1..N (top Slater {kt1},{kt2}):")
            deltas = [E_k[k-1] - E_k[kt1-1] - E_k[kt2-1] for k in range(1, N+1)]
            print(f"    {[f'k={k}:Δ={d:+.4f}' for k, d in enumerate(deltas, 1)]}")

        np.savez(out_dir / f"N{N}_sine_basis.npz",
                 sigma_0=result["sigma_0"],
                 sigma_all=result["sigma_all"],
                 U_top_comp=result["U_top"],
                 M_bond_pair=M_bond,
                 U_top_sine=U_top_sine,
                 V_top_sine=V_top_sine,
                 S_sine=S_sine,
                 mode_pairs=np.array(mode_pairs),
                 E_k=E_k)

    print()
    print(f"Phase B data saved: {out_dir}")
    print()
    print("Next: identify the (k; k_1, k_2) → σ_0 integral kernel at large N.")


if __name__ == "__main__":
    main()
