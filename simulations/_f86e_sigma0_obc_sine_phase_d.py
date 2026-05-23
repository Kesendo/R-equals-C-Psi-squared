"""F86e σ_0 closed-form Phase D: spectral decomposition of σ_0 over M_H eigenmodes.

Phase C established: U_top's bond-pair matrix decomposes with comp0 u(a) = ψ_N
(band-bottom mode), comp1 u(a) = ψ_{N−1}. The EP-partner is a band-edge object.

Phase D exact identity: σ_0 = U_top† · M_H_total · V_top.

M_H_total on the c=2 block is real-symmetric (kron(H_p1, I) − kron(I, H_p2)
with H_p1, H_p2 symmetric XY block Hamiltonians), hence orthogonally
diagonalisable: M_H |ω⟩ = Δ_ω |ω⟩, eigenvalues Δ_ω = E_k − E_{k_1} − E_{k_2},
eigenvectors |ω⟩ = |ψ_k⟩⟨ψ_{k_1} ψ_{k_2}| (the sine-operator basis).

So  σ_0 = Σ_ω Δ_ω · ⟨U_top|ω⟩* · ⟨ω|V_top⟩  =  Σ_ω Δ_ω · a_ω* · b_ω.

This is the exact spectral sum. As N → ∞ the Δ_ω become dense and the sum
turns into an integral. Phase D:
  1. verify the spectral sum reproduces σ_0 bit-exactly
  2. map the contribution density c_ω = Δ_ω a_ω* b_ω over the Δ-axis
  3. find where σ_0 accumulates (cumulative sum)
  4. extract the N-scaling of the density → the N→∞ integral

Output: simulations/results/f86e_sigma0_obc_sine_phase_d/
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")
import framework as fw  # noqa: E402


def sites_of_integer(integer: int, N: int) -> list[int]:
    return [b for b in range(N) if (integer >> b) & 1]


def hd_subspace_projector(N: int, n: int, hd_value: int):
    P_n = fw.popcount_states(N, n)
    P_np1 = fw.popcount_states(N, n + 1)
    Mnp1 = len(P_np1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}
    cols = []
    for p in P_n:
        for q in P_np1:
            if bin(p ^ q).count("1") == hd_value:
                idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
                v = np.zeros(len(P_n) * Mnp1, dtype=complex)
                v[idx] = 1.0
                cols.append(v)
    return (np.column_stack(cols) if cols else
            np.zeros((len(P_n) * Mnp1, 0), dtype=complex))


def phase_d_at(N: int, n: int = 1, gamma_0: float = 0.05):
    """Spectral decomposition of σ_0 over M_H eigenmodes.

    fw's M_H_total = −i·C, anti-Hermitian (it carries the −i of the
    Liouvillian unitary part). C = [H, ·] is real-symmetric because the
    XY chain Hamiltonian has real matrix elements. So C = i·M_H_total is
    real-symmetric and C |ω⟩ = Δ_ω |ω⟩ with Δ_ω = E_k − E_{k_1} − E_{k_2}.

      σ_0 = U_top† M_H_total V_top = −i · U_top† C V_top
          = −i · Σ_ω Δ_ω · a_ω* · b_ω
    """
    _, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H_total = sum(M_H_per_bond)

    # C = i·M_H_total should be real-symmetric.
    C = 1j * M_H_total
    c_imag_residual = float(np.linalg.norm(C.imag))
    asym = float(np.linalg.norm(C.real - C.real.T))
    C_real = C.real

    # Eigendecomposition (real symmetric → eigh)
    Delta, omega = np.linalg.eigh(C_real)  # Delta ascending, omega orthonormal

    # V_inter SVD → U_top, V_top in full block coordinates
    P_HD1 = hd_subspace_projector(N, n, 1)
    P_HD3 = hd_subspace_projector(N, n, 3)
    V_inter = P_HD1.conj().T @ M_H_total @ P_HD3
    U, S, Vh = np.linalg.svd(V_inter, full_matrices=False)
    sigma_0 = float(S[0])
    U_top = P_HD1 @ U[:, 0]    # lift to full block
    V_top = P_HD3 @ Vh[0, :].conj()

    # Spectral projections: σ_0 = −i Σ_ω Δ_ω a_ω* b_ω
    a = omega.T @ U_top         # a_ω = ⟨ω|U_top⟩  (complex)
    b = omega.T @ V_top         # b_ω = ⟨ω|V_top⟩  (complex)
    c = -1j * Delta * np.conj(a) * b   # contribution per eigenmode
    sigma_0_check = float(np.sum(c).real)
    sigma_0_imag = float(np.sum(c).imag)

    return {
        "N": N, "sigma_0": sigma_0, "sigma_0_check": sigma_0_check,
        "sigma_0_imag": sigma_0_imag,
        "Delta": Delta, "a": a, "b": b, "c": c,
        "asym": asym, "imag": c_imag_residual,
        "block_dim": C_real.shape[0],
    }


def main() -> None:
    gamma_0 = 0.05
    n = 1
    out_dir = Path("simulations/results/f86e_sigma0_obc_sine_phase_d")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("F86e Phase D: spectral decomposition σ_0 = Σ_ω Δ_ω · a_ω* · b_ω")
    print("=" * 100)

    records = []
    for N in range(5, 16):
        t0 = time.time()
        r = phase_d_at(N, n, gamma_0)
        elapsed = time.time() - t0
        err = abs(r["sigma_0"] - r["sigma_0_check"])
        print()
        print(f"--- N={N}  σ_0={r['sigma_0']:.8f}  "
              f"spectral-sum={r['sigma_0_check']:.8f}  "
              f"err={err:.2e}  ({elapsed:.2f}s) ---")
        print(f"  C=i·M_H asym={r['asym']:.1e}, imag-residual={r['imag']:.1e}, "
              f"spectral-sum imag={r['sigma_0_imag']:.1e}, block_dim={r['block_dim']}")

        Delta = r["Delta"]
        c = r["c"].real

        # Where does σ_0 accumulate? Sort eigenmodes by |contribution| descending.
        order = np.argsort(np.abs(c))[::-1]
        cum = np.cumsum(c[order])
        # How many modes for 50%, 90%, 99% of σ_0?
        targets = [0.5, 0.9, 0.99]
        n_for = {}
        for t in targets:
            idx = np.searchsorted(np.abs(cum) / abs(r["sigma_0_check"]), t)
            n_for[t] = idx + 1
        print(f"  Modes carrying 50% / 90% / 99% of σ_0: "
              f"{n_for[0.5]} / {n_for[0.9]} / {n_for[0.99]} "
              f"(of {len(c)} total)")

        # Top-8 contributing eigenmodes: Δ value + contribution
        print(f"  Top-8 contributing eigenmodes (Δ, contribution c_ω):")
        for idx in order[:8]:
            print(f"    Δ={Delta[idx]:+.5f}  c_ω={c[idx]:+.6f}  "
                  f"|a|={abs(r['a'][idx]):.4f} |b|={abs(r['b'][idx]):.4f}")

        # Contribution binned by Δ
        # Positive vs negative Δ contributions
        pos_c = float(np.sum(c[Delta > 1e-9]))
        neg_c = float(np.sum(c[Delta < -1e-9]))
        zero_c = float(np.sum(c[np.abs(Delta) <= 1e-9]))
        print(f"  Σc over Δ>0: {pos_c:+.6f}   Δ<0: {neg_c:+.6f}   Δ≈0: {zero_c:+.6f}")

        # Δ range and where the dominant contributions sit
        dom_Delta = Delta[order[:n_for[0.9]]]
        print(f"  Δ-range of the 90%-carrying modes: "
              f"[{dom_Delta.min():+.4f}, {dom_Delta.max():+.4f}]  "
              f"(full Δ-range [{Delta.min():+.4f}, {Delta.max():+.4f}])")

        records.append(r)
        np.savez(out_dir / f"N{N}_phase_d.npz",
                 Delta=Delta, a=r["a"], b=r["b"], c=c,
                 sigma_0=r["sigma_0"])

    # Cross-N: how does the contribution structure scale?
    print()
    print("=" * 100)
    print("Cross-N contribution structure")
    print("=" * 100)
    print(f"{'N':>3}  {'σ_0':>10}  {'#modes/90%':>11}  {'frac modes':>11}  "
          f"{'Δ>0 sum':>10}  {'Δ<0 sum':>10}")
    for r in records:
        Delta = r["Delta"]
        c = r["c"].real
        order = np.argsort(np.abs(c))[::-1]
        cum = np.cumsum(c[order])
        idx90 = np.searchsorted(np.abs(cum) / abs(r["sigma_0_check"]), 0.9) + 1
        pos_c = float(np.sum(c[Delta > 1e-9]))
        neg_c = float(np.sum(c[Delta < -1e-9]))
        print(f"{r['N']:>3}  {r['sigma_0']:>10.6f}  {idx90:>11}  "
              f"{idx90 / len(c):>11.4f}  {pos_c:>+10.5f}  {neg_c:>+10.5f}")

    print()
    print("Reading guide:")
    print("  If '#modes/90%' grows ∝ N, the σ_0 sum is a genuine Riemann sum →")
    print("    integral. The Δ-range of the 90%-modes is the integration domain.")
    print("  If it stays O(1), σ_0 is carried by a fixed finite set of modes →")
    print("    σ_0(∞) is a finite closed-form expression in E_k limits.")
    print(f"\nData saved: {out_dir}")


if __name__ == "__main__":
    main()
