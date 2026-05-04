"""F86 EP-resonance Q-scan on a 5-proton water chain (Grotthuss / Heisenberg +
Z-dephasing) at the popcount-(2, 3) coherence block (chromaticity c = 3).

Reports HWHM-/Q* per bond class against F86's Tier-1-candidate prediction
(Interior ≈ 0.756, Endpoint ≈ 0.770) and a per-proton Bloch + Π²-odd-fraction
state-level reading at the interior Q_peak. Block-restricted L via the framework
primitives in `simulations/framework/coherence_block.py`; per-bond K_CC_pr
observable

    K_b(Q, t) = 2 · Re ⟨ρ(t) | S_kernel | ∂ρ/∂J_b⟩

with ρ₀ the Dicke probe and ∂ρ/∂J_b from the Duhamel formula (the simple
coherence-block-norm observable is monotonic in γ; only the J-derivative
picks up the EP).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import framework as fw  # noqa: E402


# Pauli helpers for the state-level diagnostic. Mirror the same in
# `proton_chain_memory_reading.py`; lift to a shared module if a third caller appears.
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, SX, SY, SZ]
BIT_B = [0, 0, 1, 1]


def kron_n(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def site_op(P, k, N):
    return kron_n([P if i == k else I2 for i in range(N)])


def heisenberg_chain(N, J=1.0):
    H = np.zeros((2**N, 2**N), dtype=complex)
    for b in range(N - 1):
        for P in [SX, SY, SZ]:
            H += (J / 4.0) * site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def z_dephasing_lindblad(H, gamma, N):
    d = 2**N
    eye = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, eye) - np.kron(eye, H.T))
    for k in range(N):
        Zk = site_op(SZ, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.kron(eye, eye))
    return L


def diagonalise(L):
    lambdas, R = np.linalg.eig(L)
    return lambdas, R, np.linalg.inv(R)


def evolve(rho0, lambdas, R, R_inv, t):
    d = rho0.shape[0]
    vec0 = rho0.flatten(order='F')
    c = R_inv @ vec0
    vec_t = R @ (np.exp(lambdas * t) * c)
    rho_t = vec_t.reshape((d, d), order='F')
    return (rho_t + rho_t.conj().T) / 2.0


def per_proton_bloch_mag(rho, N):
    return [
        np.sqrt(
            np.trace(site_op(SX, k, N) @ rho).real ** 2
            + np.trace(site_op(SY, k, N) @ rho).real ** 2
            + np.trace(site_op(SZ, k, N) @ rho).real ** 2
        )
        for k in range(N)
    ]


def pi2_odd_fraction_in_memory(rho, N):
    """Frobenius-fraction of Π²-odd Pauli content in (rho - kernel projection)."""
    d = 2**N
    rho_d0 = np.zeros_like(rho)
    for n in range(N + 1):
        P_n = np.zeros((d, d), dtype=complex)
        for b in range(d):
            if bin(b).count('1') == n:
                P_n[b, b] = 1.0
        rank_n = int(np.real(np.trace(P_n @ P_n)))
        if rank_n == 0:
            continue
        coeff = np.real(np.trace(P_n @ rho)) / rank_n
        rho_d0 = rho_d0 + coeff * P_n
    rho_d2 = rho - rho_d0
    rho_odd = np.zeros_like(rho)
    inv = 1.0 / d
    for k in range(4**N):
        kk = k
        idxs = []
        for _ in range(N):
            idxs.append(kk & 3)
            kk >>= 2
        if (sum(BIT_B[i] for i in idxs) & 1) == 0:
            continue
        sigma = kron_n([PAULIS[i] for i in idxs])
        coeff = np.trace(sigma @ rho_d2) * inv
        if abs(coeff) > 1e-14:
            rho_odd = rho_odd + coeff * sigma
    norm_total = np.linalg.norm(rho, 'fro') ** 2
    norm_static = np.linalg.norm(rho_d0, 'fro') ** 2
    norm_mem = np.linalg.norm(rho_d2, 'fro') ** 2
    norm_odd = np.linalg.norm(rho_odd, 'fro') ** 2
    return (
        norm_static / norm_total if norm_total > 0 else 0.0,
        norm_mem / norm_total if norm_total > 0 else 0.0,
        norm_odd / norm_mem if norm_mem > 1e-12 else 0.0,
    )


# ---------- F86 K_CC_pr Q-scan (block-restricted, framework primitives) ----------

def per_bond_K_curve(N, n, gamma_0, J_grid, t_grid):
    """Per-bond |∂S/∂J_b|_max(Q) via the Duhamel formula on the (n, n+1) block.

    Mirrors `simulations/_eq022_b1_step_e_resonance_shape.py`; lift to a
    workflow if a third caller appears. Returns (K_max_over_t, t_at_max),
    each shape (n_bonds, len(J_grid)).
    """
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    rho0 = fw.dicke_block_probe(N, n)
    S_kernel = fw.spatial_sum_coherence_kernel(N, n)
    n_bonds = N - 1
    M_H_total = sum(M_H_per_bond)

    K_curves = np.zeros((n_bonds, len(J_grid)))
    t_curves = np.zeros((n_bonds, len(J_grid)))

    for i_J, J in enumerate(J_grid):
        L = D + J * M_H_total
        evals, R = np.linalg.eig(L)
        R_inv = np.linalg.inv(R)
        c0 = R_inv @ rho0
        X_b_list = [R_inv @ Mb @ R for Mb in M_H_per_bond]

        bond_K_max = np.zeros(n_bonds)
        bond_t_at_max = np.zeros(n_bonds)
        for t in t_grid:
            e = np.exp(evals * t)
            lam_j = evals[:, None]
            lam_k = evals[None, :]
            with np.errstate(divide='ignore', invalid='ignore'):
                I_mat = np.where(np.abs(lam_k - lam_j) > 1e-10,
                                 (e[None, :] - e[:, None]) / (lam_k - lam_j),
                                 t * e[:, None])
            rho_t = R @ (e * c0)
            for b in range(n_bonds):
                F_b = X_b_list[b] * I_mat
                drho_b = R @ (F_b @ c0)
                K = 2.0 * float(np.real(np.vdot(rho_t, S_kernel @ drho_b)))
                K_abs = abs(K)
                if K_abs > bond_K_max[b]:
                    bond_K_max[b] = K_abs
                    bond_t_at_max[b] = t

        K_curves[:, i_J] = bond_K_max
        t_curves[:, i_J] = bond_t_at_max

    return K_curves, t_curves


def find_peak_with_interp(Q_grid, K_curve):
    """Parabolic-interpolated peak Q* and |K|max, plus linear-interp HWHM."""
    i_max = int(np.argmax(K_curve))
    K_max = K_curve[i_max]

    if 0 < i_max < len(Q_grid) - 1:
        x = Q_grid[i_max - 1: i_max + 2]
        y = K_curve[i_max - 1: i_max + 2]
        coefs = np.polyfit(x, y, 2)
        if coefs[0] < 0:
            Q_star = -coefs[1] / (2 * coefs[0])
            K_max_interp = coefs[2] - coefs[1] ** 2 / (4 * coefs[0])
            if abs(Q_star - Q_grid[i_max]) <= (Q_grid[i_max + 1] - Q_grid[i_max - 1]):
                K_max = K_max_interp
            else:
                Q_star = Q_grid[i_max]
        else:
            Q_star = Q_grid[i_max]
    else:
        Q_star = Q_grid[i_max]

    half = K_max / 2.0
    hwhm_left = None
    hwhm_right = None
    for i in range(i_max, -1, -1):
        if K_curve[i] < half:
            x0, x1 = Q_grid[i], Q_grid[i + 1]
            y0, y1 = K_curve[i], K_curve[i + 1]
            x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_left = Q_star - x_half
            break
    for i in range(i_max, len(Q_grid)):
        if K_curve[i] < half:
            x0, x1 = Q_grid[i - 1], Q_grid[i]
            y0, y1 = K_curve[i - 1], K_curve[i]
            x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_right = x_half - Q_star
            break

    return Q_star, K_max, hwhm_left, hwhm_right


# ---------- main ----------

def coherence_state(N, p_bits, q_bits):
    psi = np.zeros(2**N, dtype=complex)
    psi[p_bits] = 1.0 / np.sqrt(2.0)
    psi[q_bits] = 1.0 / np.sqrt(2.0)
    return np.outer(psi, psi.conj())


def main():
    N = 5
    n = 2
    gamma_0 = 0.05
    Q_grid = np.arange(0.20, 6.01, 0.025)
    J_grid = Q_grid * gamma_0
    t_peak = 1.0 / (4.0 * gamma_0)
    t_grid = np.linspace(0.6 * t_peak, 1.6 * t_peak, 21)

    print(f"F86 EP-resonance on water chain N={N}, popcount-({n}, {n + 1}), c={fw.chromaticity(N, n)}")
    print("=" * 72)
    print(f"γ₀ = {gamma_0}, Q ∈ [{Q_grid[0]}, {Q_grid[-1]}] (dQ = {Q_grid[1] - Q_grid[0]:.3f})")
    print(f"J = Q·γ₀ ∈ [{J_grid[0]:.4f}, {J_grid[-1]:.4f}]")
    print(f"t_peak = 1/(4γ₀) = {t_peak:.2f}, t-grid 21 points across [0.6, 1.6]·t_peak")
    print()

    K_curves, t_curves = per_bond_K_curve(N, n, gamma_0, J_grid, t_grid)

    n_bonds = N - 1
    interior_b = 2

    print(f"{'bond':>5} {'class':>9} {'Q_peak':>9} {'|K|_max':>10} {'HWHM-/Q*':>10} {'HWHM+/Q*':>10}")
    print("-" * 60)
    summary = []
    for b in range(n_bonds):
        cls = "Interior" if 0 < b < n_bonds - 1 else "Endpoint"
        Q_star, K_max, hl, hr = find_peak_with_interp(Q_grid, K_curves[b])
        rl = hl / Q_star if hl is not None and Q_star > 0 else float('nan')
        rr = hr / Q_star if hr is not None and Q_star > 0 else float('nan')
        print(f"{b:>5} {cls:>9} {Q_star:>9.4f} {K_max:>10.4e} {rl:>10.4f} {rr:>10.4f}")
        summary.append((b, cls, Q_star, K_max, hl, hr, rl, rr))

    print()
    print("F86 Tier-1-candidate prediction (`project_q_peak_ep_structure`):")
    print("  Interior HWHM-/Q* ≈ 0.756 ± 0.005")
    print("  Endpoint HWHM-/Q* ≈ 0.770")
    print()

    interior = [s for s in summary if s[1] == "Interior"]
    endpoint = [s for s in summary if s[1] == "Endpoint"]
    if interior:
        avg_int = np.mean([s[6] for s in interior if s[6] == s[6]])
        print(f"Interior (this run, {len(interior)} bonds) HWHM-/Q* mean = {avg_int:.4f}  (Δ vs 0.756: {avg_int - 0.756:+.4f})")
    if endpoint:
        avg_end = np.mean([s[6] for s in endpoint if s[6] == s[6]])
        print(f"Endpoint (this run, {len(endpoint)} bonds) HWHM-/Q* mean = {avg_end:.4f}  (Δ vs 0.770: {avg_end - 0.770:+.4f})")
    print()

    print("=== Universal shape, Interior bond b={} ===".format(interior_b))
    print("(Q − Q_peak)/Q_peak       K(Q)/|K|_max")
    print("-" * 40)
    Q_star_int = summary[interior_b][2]
    K_max_int = summary[interior_b][3]
    show_idx = list(range(0, len(Q_grid), max(1, len(Q_grid) // 20)))
    for i in show_idx:
        x = (Q_grid[i] - Q_star_int) / Q_star_int
        y = K_curves[interior_b, i] / K_max_int
        bar = "#" * int(60 * y)
        print(f"  {x:+8.3f}             {y:.3f}  {bar}")
    print()

    print(f"=== State-level reading at γ_at_peak = J/Q* = {1.0 / Q_star_int:.4f} ===")
    print("Probe: |ψ⟩ = (|00011⟩ + |00111⟩)/√2, popcount-(2, 3) HD=1 coherence state.")
    print("F86 prediction: Π²-odd / memory = 10/19 = 0.5263 exactly at t = 0")
    print("(popcount-mirror at n_p+n_q=N; see PROTON_WATER_CHAIN.md §EP-Resonance Inheritance).")
    print()
    print(f"{'t':>6} {'static':>9} {'memory':>9} {'Π²-odd/mem':>12}  per-proton |r|")
    print("-" * 80)

    p_bits = 0b00011
    q_bits = 0b00111
    H = heisenberg_chain(N, J=1.0)
    gamma_at_peak = 1.0 / Q_star_int
    L_full = z_dephasing_lindblad(H, gamma_at_peak, N)
    lambdas, R, R_inv = diagonalise(L_full)
    rho0 = coherence_state(N, p_bits, q_bits)

    for t in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
        rho_t = evolve(rho0, lambdas, R, R_inv, t)
        s, m, o = pi2_odd_fraction_in_memory(rho_t, N)
        rs = per_proton_bloch_mag(rho_t, N)
        bs = "  ".join(f"q{k}={r:.3f}" for k, r in enumerate(rs))
        print(f"{t:>6.1f} {s:>9.4f} {m:>9.4f} {o:>12.4f}  {bs}")

    print()
    print("=== Reading ===")
    print()
    print(f"Per-bond Q_peak values are chain-specific (no closed form survives the")
    print(f"two retractions in `docs/proofs/PROOF_F86_QPEAK.md`). Universal shape is")
    print(f"the structural law that survived: K(Q)/|K|_max = f(Q/Q_peak) within bond")
    print(f"class. Water inherits the algebra cleanly: same Heisenberg + Z-dephasing.")
    print()
    print(f"State-level diagnostic at the popcount-(2,3) HD=1 coherence state shows")
    print(f"Π²-odd/memory = 10/19 at t=0 (popcount-mirror at n_p+n_q=N) and tracks how")
    print(f"Heisenberg + Z-dephasing migrate the proton chain's per-proton Bloch magnitude.")


if __name__ == "__main__":
    main()
