"""F86e σ_0 closed-form Phase C: rank-1 sine-mode decomposition of U_top.

Phase B established: U_top (the EP-partner mode of V_inter in the HD=1
subspace) reshapes into a bond-pair matrix M_full[a, b] (N×N, zero diagonal)
of rank 2 at odd N, rank 3 at even N. The rank-2 pieces at odd N are
chain-mirror (R) partners.

Phase C goal: SVD M_full per N, expand each dominant rank-1 component's
left vector u(a) and right vector v(b) in the OBC sine basis. If u, v are
single sine modes ψ_k, σ_0(N) becomes a finite sine sum and σ_0(N→∞) the
corresponding band-edge integral.

For each N=5..15:
  - M_full rank, chain-mirror R-parity, imag residual
  - top-2 rank-1 components: σ_i, and u(a)/v(b) sine-mode content
  - whether u, v are single modes (clean) or combinations
  - the implied σ_0(N) reconstruction from the sine content

Output: simulations/results/f86e_sigma0_obc_sine_phase_c/
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


def obc_sine_modes(N: int) -> np.ndarray:
    """psi[k-1, i] = √(2/(N+1))·sin(πk(i+1)/(N+1)), k=1..N, i=0..N-1."""
    psi = np.zeros((N, N))
    for k in range(1, N + 1):
        for i in range(N):
            psi[k - 1, i] = math.sqrt(2.0 / (N + 1)) * math.sin(
                math.pi * k * (i + 1) / (N + 1))
    return psi


def obc_dispersion(N: int, J: float = 1.0) -> np.ndarray:
    return np.array([2 * J * math.cos(math.pi * k / (N + 1))
                     for k in range(1, N + 1)])


def sites_of_integer(integer: int, N: int) -> list[int]:
    return [b for b in range(N) if (integer >> b) & 1]


def hd_subspace_projector(N: int, n: int, hd_value: int):
    P_n = fw.popcount_states(N, n)
    P_np1 = fw.popcount_states(N, n + 1)
    Mnp1 = len(P_np1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}
    cols, labels = [], []
    for p in P_n:
        for q in P_np1:
            if bin(p ^ q).count("1") == hd_value:
                idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
                v = np.zeros(len(P_n) * Mnp1, dtype=complex)
                v[idx] = 1.0
                cols.append(v)
                labels.append((p, q))
    return (np.column_stack(cols) if cols else
            np.zeros((len(P_n) * Mnp1, 0), dtype=complex)), labels


def u_top_full_matrix(N: int, n: int = 1, gamma_0: float = 0.05):
    """M_full[a, b] (N×N, zero diagonal) = U_top entry for HD-1 state ({a},{a,b})."""
    _, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H_total = sum(M_H_per_bond)
    P_HD1, hd1_labels = hd_subspace_projector(N, n, 1)
    P_HD3, _ = hd_subspace_projector(N, n, 3)
    V_inter = P_HD1.conj().T @ M_H_total @ P_HD3
    U, S, Vh = np.linalg.svd(V_inter, full_matrices=False)
    U_top = U[:, 0]

    M_full = np.zeros((N, N), dtype=complex)
    for (p, q), amp in zip(hd1_labels, U_top):
        a = sites_of_integer(p, N)[0]
        b_sites = [s for s in sites_of_integer(q, N) if s != a]
        M_full[a, b_sites[0]] = amp
    return M_full, float(S[0]), S


def identify_sine_content(vec: np.ndarray, psi: np.ndarray):
    """Expand vec (length N) in sine basis. Return (overlaps, dominant_k, purity)."""
    overlaps = psi @ vec  # overlaps[k-1] = Σ_a ψ_k(a)·vec(a)
    dominant = int(np.argmax(np.abs(overlaps)))
    purity = abs(overlaps[dominant]) ** 2 / np.sum(np.abs(overlaps) ** 2)
    return overlaps, dominant + 1, purity


def main() -> None:
    gamma_0 = 0.05
    n = 1
    out_dir = Path("simulations/results/f86e_sigma0_obc_sine_phase_c")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("F86e Phase C: rank-1 sine-mode decomposition of U_top bond-pair matrix")
    print("=" * 100)

    records = []
    for N in range(5, 16):
        t0 = time.time()
        psi = obc_sine_modes(N)
        M_full, sigma_0, S_all = u_top_full_matrix(N, n, gamma_0)

        # Fix global phase: largest |entry| → real positive
        idx_max = np.unravel_index(int(np.argmax(np.abs(M_full))), M_full.shape)
        phase = M_full[idx_max] / abs(M_full[idx_max])
        M_full = M_full / phase
        imag_norm = float(np.linalg.norm(M_full.imag))
        M_real = M_full.real

        # SVD
        u_svd, s_svd, vh_svd = np.linalg.svd(M_real)
        rank = int(np.sum(s_svd > 1e-9))

        # Chain-mirror R-parity
        R = np.fliplr(np.eye(N))
        M_reflected = R @ M_real @ R
        r_even = float(np.linalg.norm(M_reflected - M_real))
        r_odd = float(np.linalg.norm(M_reflected + M_real))
        parity = "even" if r_even < r_odd else "odd"

        elapsed = time.time() - t0
        print()
        print(f"--- N={N}  σ_0={sigma_0:.6f}  rank(M_full)={rank}  "
              f"R-{parity}  imag={imag_norm:.1e}  ({elapsed:.2f}s) ---")
        print(f"  M_full singular values: {[f'{s:.5f}' for s in s_svd[:5]]}")
        sv_sq_sum = float(np.sum(s_svd ** 2))
        print(f"  Σσ_i² = {sv_sq_sum:.6f} (should be 1, U_top normalised)")

        comp_records = []
        for comp in range(min(3, rank)):
            u_comp = u_svd[:, comp]
            v_comp = vh_svd[comp, :]
            sv = float(s_svd[comp])

            u_ov, u_k, u_pur = identify_sine_content(u_comp, psi)
            v_ov, v_k, v_pur = identify_sine_content(v_comp, psi)

            u_clean = "single" if u_pur > 0.98 else (
                "near" if u_pur > 0.85 else "mixed")
            v_clean = "single" if v_pur > 0.98 else (
                "near" if v_pur > 0.85 else "mixed")

            print(f"  comp {comp} (σ={sv:.5f}):")
            print(f"    u(a): dominant sine k={u_k} (purity {u_pur:.4f}, {u_clean})  "
                  + "  ".join([f"k{kk+1}:{u_ov[kk]:+.3f}"
                               for kk in np.argsort(np.abs(u_ov))[::-1][:4]]))
            print(f"    v(b): dominant sine k={v_k} (purity {v_pur:.4f}, {v_clean})  "
                  + "  ".join([f"k{kk+1}:{v_ov[kk]:+.3f}"
                               for kk in np.argsort(np.abs(v_ov))[::-1][:4]]))
            comp_records.append({
                "comp": comp, "sigma": sv,
                "u_k": u_k, "u_purity": u_pur,
                "v_k": v_k, "v_purity": v_pur,
            })

        records.append({
            "N": N, "sigma_0": sigma_0, "rank": rank, "parity": parity,
            "M_singular": s_svd, "components": comp_records,
        })
        np.savez(out_dir / f"N{N}_phase_c.npz",
                 M_full=M_full, M_singular=s_svd,
                 sigma_0=sigma_0, E_k=obc_dispersion(N))

    # Cross-N pattern in dominant sine modes
    print()
    print("=" * 100)
    print("Cross-N pattern: dominant sine modes (k) of the rank-1 components")
    print("=" * 100)
    print(f"{'N':>3}  {'σ_0':>9}  {'rank':>4}  {'parity':>6}  "
          f"{'comp0 (σ,u_k,v_k,pur)':>34}  {'comp1 (σ,u_k,v_k,pur)':>34}")
    for r in records:
        c = r["components"]
        c0 = c[0] if len(c) > 0 else None
        c1 = c[1] if len(c) > 1 else None
        c0s = (f"σ={c0['sigma']:.3f} u_k={c0['u_k']} v_k={c0['v_k']} "
               f"p={min(c0['u_purity'], c0['v_purity']):.2f}") if c0 else "-"
        c1s = (f"σ={c1['sigma']:.3f} u_k={c1['u_k']} v_k={c1['v_k']} "
               f"p={min(c1['u_purity'], c1['v_purity']):.2f}") if c1 else "-"
        print(f"{r['N']:>3}  {r['sigma_0']:>9.5f}  {r['rank']:>4}  "
              f"{r['parity']:>6}  {c0s:>34}  {c1s:>34}")

    print()
    print("Reading guide:")
    print("  If u_k, v_k are stable single modes across N → σ_0 has a finite sine sum.")
    print("  If purity is high (>0.98) the components ARE single sine modes.")
    print("  Phase D would then write σ_0(N) as the sum and take the N→∞ integral.")
    print(f"\nData saved: {out_dir}")


if __name__ == "__main__":
    main()
