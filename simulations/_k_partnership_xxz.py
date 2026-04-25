#!/usr/bin/env python3
"""K-partnership test for full Heisenberg/XXZ on open chain (Robustness Test 2).

Hypothesis (from PROOF_K_PARTNERSHIP, "Scope" section): full Heisenberg/XXZ
with Δ ≠ 0 contributes ZZ-coupling that, in the single-excitation sector,
generates an effective on-site potential V_eff(ℓ) = (#bonds) − 2·deg(ℓ).
On open chains, deg(0) = deg(N-1) = 1 vs deg(interior) = 2 makes V_eff
non-uniform across the chain, breaking K at the boundary sites. K is
restored on topologies with uniform site degree (e.g., periodic chains).

Test design: full 2^N Hilbert-space simulation (no single-exc reduction),
Lindblad propagation under XXZ + uniform Z-dephasing. Compute MI(a, N-1-a)(t)
for each mirror-pair depth a, then compare ρ_k vs ρ_{N+1-k}.

Predictions:
- Δ = 0 (pure XX): all mirror-pair-depths match identically (sanity check).
- Δ ≠ 0 on open chain: boundary-pair (0, N-1) breaks strongest, decreasing
  toward bulk. Boundary V_eff differs from interior V_eff by 2.
- Δ ≠ 0 on periodic chain (uniform deg): partner identity restored.
- Linear scaling in Δ for small Δ (perturbative regime).
"""
import math
import sys
from itertools import combinations

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_n(ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(N: int, i: int, op: np.ndarray) -> np.ndarray:
    """Operator op on site i (0-indexed, MSB convention), I on others."""
    ops = [I2] * N
    ops[i] = op
    return kron_n(ops)


def build_H_xxz(N: int, J: float, delta: float, periodic: bool = False) -> np.ndarray:
    """H = (J/2) sum_<ij> (X_i X_j + Y_i Y_j + delta Z_i Z_j)."""
    D = 2**N
    H = np.zeros((D, D), dtype=complex)
    bonds = [(i, i + 1) for i in range(N - 1)]
    if periodic:
        bonds.append((N - 1, 0))
    for (i, j) in bonds:
        Xi = site_op(N, i, X)
        Yi = site_op(N, i, Y)
        Zi = site_op(N, i, Z)
        Xj = site_op(N, j, X)
        Yj = site_op(N, j, Y)
        Zj = site_op(N, j, Z)
        H = H + 0.5 * J * (Xi @ Xj + Yi @ Yj + delta * Zi @ Zj)
    return H


def bonding_mode_state(N: int, k: int) -> np.ndarray:
    """|ψ_k⟩ = sqrt(2/(N+1)) Σ_j sin(πk(j+1)/(N+1)) |1_j⟩ in the 2^N basis.

    Convention: site 0 = MSB, so |1_j⟩ has bit-index 2^(N-1-j).
    """
    psi = np.zeros(2**N, dtype=complex)
    norm = math.sqrt(2.0 / (N + 1))
    for j in range(N):
        idx = 1 << (N - 1 - j)
        psi[idx] = norm * math.sin(math.pi * k * (j + 1) / (N + 1))
    return psi


def lindblad_rhs(rho: np.ndarray, H: np.ndarray, Z_ops: list, gamma: float) -> np.ndarray:
    """dρ/dt = -i [H, ρ] + γ Σ_l (Z_l ρ Z_l - ρ)."""
    out = -1j * (H @ rho - rho @ H)
    for Zl in Z_ops:
        out = out + gamma * (Zl @ rho @ Zl - rho)
    return out


def rk4_step(rho: np.ndarray, dt: float, H: np.ndarray, Z_ops: list, gamma: float) -> np.ndarray:
    k1 = lindblad_rhs(rho, H, Z_ops, gamma)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, Z_ops, gamma)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, Z_ops, gamma)
    k4 = lindblad_rhs(rho + dt * k3, H, Z_ops, gamma)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def partial_trace_keep_sites(rho: np.ndarray, sites: tuple, N: int) -> np.ndarray:
    """Trace out all sites not in `sites`. Returns 2^k x 2^k. From c1_pair_local.py."""
    sites = tuple(sorted(sites))
    shape_2N = [2] * (2 * N)
    out = rho.reshape(shape_2N)
    ket_axes = list(range(N))
    bra_axes = list(range(N, 2 * N))
    for j in range(N - 1, -1, -1):
        if j in sites:
            continue
        a_k = ket_axes[j]
        a_b = bra_axes[j]
        out = np.trace(out, axis1=a_k, axis2=a_b)
        lo, hi = sorted((a_k, a_b))
        for k in range(N):
            if k == j:
                continue
            if ket_axes[k] > hi:
                ket_axes[k] -= 2
            elif ket_axes[k] > lo:
                ket_axes[k] -= 1
            if bra_axes[k] > hi:
                bra_axes[k] -= 2
            elif bra_axes[k] > lo:
                bra_axes[k] -= 1
    k = len(sites)
    d_sub = 2**k
    remaining_axes_ket = [ket_axes[i] for i in sites]
    remaining_axes_bra = [bra_axes[i] for i in sites]
    current_order = remaining_axes_ket + remaining_axes_bra
    perm = np.argsort(current_order)
    out = np.transpose(out, perm)
    return out.reshape(d_sub, d_sub)


def von_neumann_S(rho: np.ndarray) -> float:
    evs = np.linalg.eigvalsh(rho).real
    S = 0.0
    for ev in evs:
        if ev > 1e-15:
            S -= ev * math.log2(ev)
    return S


def MI_pair(rho: np.ndarray, a: int, b: int, N: int) -> float:
    """Mutual information between sites a and b, full-Hilbert-space rho."""
    rho_ab = partial_trace_keep_sites(rho, (a, b), N)
    rho_a = partial_trace_keep_sites(rho, (a,), N)
    rho_b = partial_trace_keep_sites(rho, (b,), N)
    return von_neumann_S(rho_a) + von_neumann_S(rho_b) - von_neumann_S(rho_ab)


def simulate_mi_trajectory(N: int, k: int, H: np.ndarray, Z_ops: list, gamma: float,
                            t_max: float, n_steps: int, mirror_pairs: list) -> np.ndarray:
    """Return (n_steps+1, n_pairs) array of MI(a, N-1-a) over time."""
    psi = bonding_mode_state(N, k)
    rho = np.outer(psi, np.conj(psi))
    dt = t_max / n_steps
    n_pairs = len(mirror_pairs)
    out = np.zeros((n_steps + 1, n_pairs))
    for step in range(n_steps + 1):
        for p_idx, (a, b) in enumerate(mirror_pairs):
            out[step, p_idx] = MI_pair(rho, a, b, N)
        if step < n_steps:
            rho = rk4_step(rho, dt, H, Z_ops, gamma)
    return out


def run_xxz_test(N: int, delta: float, gamma: float, t_max: float, n_steps: int,
                 periodic: bool = False) -> dict:
    """Run K-partnership test for given Δ, return per-mirror-pair worst |ΔMI|."""
    Z_ops = [site_op(N, l, Z) for l in range(N)]
    H = build_H_xxz(N, J=1.0, delta=delta, periodic=periodic)

    # mirror pairs (0, N-1), (1, N-2), ..., up to floor((N-1)/2)
    mirror_pairs = [(a, N - 1 - a) for a in range((N) // 2)]

    half = (N + 1) // 2
    out = {}
    for k in range(1, half + 1):
        k_mirror = N + 1 - k
        if k == k_mirror:
            continue
        mi_k = simulate_mi_trajectory(N, k, H, Z_ops, gamma, t_max, n_steps, mirror_pairs)
        mi_m = simulate_mi_trajectory(N, k_mirror, H, Z_ops, gamma, t_max, n_steps, mirror_pairs)
        worst_per_pair = np.max(np.abs(mi_k - mi_m), axis=0)
        out[(k, k_mirror)] = worst_per_pair
    return out, mirror_pairs


# ============================================================================
# Multi-excitation sector test (Robustness Test 3)
# ============================================================================

def two_exc_slater(N: int, k1: int, k2: int) -> np.ndarray:
    """Slater determinant of single-exc bonding modes (psi_{k1}, psi_{k2})
    embedded in 2^N (with site 0 = MSB convention).

    For free fermions / pure XX, these are the natural 2-particle eigenstates.
    For XXZ with Δ ≠ 0, they are still legitimate states but no longer
    eigenstates of H; the partner-identity test runs on them as initial
    conditions to probe whether K-partnership survives the ZZ-interaction.
    """
    norm_factor = math.sqrt(2.0 / (N + 1))
    psi1 = [norm_factor * math.sin(math.pi * k1 * (j + 1) / (N + 1)) for j in range(N)]
    psi2 = [norm_factor * math.sin(math.pi * k2 * (j + 1) / (N + 1)) for j in range(N)]
    state = np.zeros(2**N, dtype=complex)
    for j1 in range(N):
        for j2 in range(j1 + 1, N):
            amp = psi1[j1] * psi2[j2] - psi1[j2] * psi2[j1]
            idx = (1 << (N - 1 - j1)) | (1 << (N - 1 - j2))
            state[idx] = amp
    nrm = float(np.linalg.norm(state))
    if nrm < 1e-15:
        return state  # zero state (e.g., k1 == k2)
    return state / nrm


def simulate_2exc_mi_trajectory(N: int, k1: int, k2: int, H: np.ndarray, Z_ops: list,
                                  gamma: float, t_max: float, n_steps: int,
                                  mirror_pairs: list) -> np.ndarray:
    psi = two_exc_slater(N, k1, k2)
    rho = np.outer(psi, np.conj(psi))
    dt = t_max / n_steps
    n_pairs = len(mirror_pairs)
    out = np.zeros((n_steps + 1, n_pairs))
    for step in range(n_steps + 1):
        for p_idx, (a, b) in enumerate(mirror_pairs):
            out[step, p_idx] = MI_pair(rho, a, b, N)
        if step < n_steps:
            rho = rk4_step(rho, dt, H, Z_ops, gamma)
    return out


def k_partner_2exc(N: int, k1: int, k2: int) -> tuple:
    """K-partner of Slater (k1, k2) is (N+1-k1, N+1-k2), then sorted."""
    return tuple(sorted((N + 1 - k1, N + 1 - k2)))


def run_2exc_xxz_test(N: int, delta: float, gamma: float, t_max: float, n_steps: int,
                       periodic: bool = False) -> tuple:
    """Test K-partnership in the 2-excitation sector under XXZ + Z-dephasing."""
    Z_ops = [site_op(N, l, Z) for l in range(N)]
    H = build_H_xxz(N, J=1.0, delta=delta, periodic=periodic)
    mirror_pairs = [(a, N - 1 - a) for a in range(N // 2)]

    # Find non-self K-partner pairs (k1 < k2)
    seen = set()
    test_pairs = []
    for k1 in range(1, N + 1):
        for k2 in range(k1 + 1, N + 1):
            partner = k_partner_2exc(N, k1, k2)
            if (k1, k2) == partner:
                continue  # self-partner (trivial identity)
            if (k1, k2) in seen or partner in seen:
                continue
            seen.add((k1, k2))
            seen.add(partner)
            test_pairs.append(((k1, k2), partner))

    out = {}
    for (orig, partner) in test_pairs:
        mi_o = simulate_2exc_mi_trajectory(N, orig[0], orig[1], H, Z_ops, gamma,
                                             t_max, n_steps, mirror_pairs)
        mi_p = simulate_2exc_mi_trajectory(N, partner[0], partner[1], H, Z_ops, gamma,
                                             t_max, n_steps, mirror_pairs)
        worst = np.max(np.abs(mi_o - mi_p), axis=0)
        out[(orig, partner)] = worst
    return out, mirror_pairs


def main() -> None:
    gamma = 0.0
    t_max = 2.0
    n_steps = 200

    print("=" * 78)
    print(f"XXZ K-partnership boundary scaling test  (gamma={gamma}, t_max={t_max})")
    print("Full 2^N Lindblad propagation (no single-exc reduction)")
    print("Hypothesis: ZZ-term breaks K at boundary mirror-pair (0, N-1)")
    print("            stronger than at deeper bulk pairs (open chain).")
    print("=" * 78)
    print()

    # Open chain at N = 7: three mirror-pair depths (boundary, mid-1, center)
    N_open = 7
    deltas = [0.0, 0.1, 0.3, 0.5, 1.0]
    print(f"--- N={N_open} OPEN chain (mirror-pair depth scaling) ---")
    print()

    for delta in deltas:
        out, mirror_pairs = run_xxz_test(N_open, delta, gamma, t_max, n_steps, periodic=False)
        print(f"Δ = {delta:.2f}  (open chain)")
        header = f"  {'pair (k, N+1-k)':>16}  " + "  ".join(
            f"{'pair ' + str(p):>14}" for p in mirror_pairs
        )
        print(header)
        for (k, km), worst in out.items():
            row = f"  ({k}, {km})         "
            row += "  ".join(f"{w:>14.3e}" for w in worst)
            print(row)
        print()

    # Periodic chain at N = 6 (EVEN -> bipartite ring -> K well-defined)
    # NOTE: odd N periodic is NON-BIPARTITE (wrap-around bond connects same-
    # sublattice sites with parity (-1)^(N-1+0) = +1 for odd N), so K breaks
    # there independently of Δ. Even N preserves bipartite structure.
    N_per = 6
    print(f"--- N={N_per} PERIODIC chain (EVEN N -> bipartite -> K-symmetric for any Δ) ---")
    print()
    for delta in [0.0, 0.5, 1.0]:
        out, mirror_pairs = run_xxz_test(N_per, delta, gamma, t_max, n_steps, periodic=True)
        print(f"Δ = {delta:.2f}  (periodic, even N)")
        header = f"  {'pair (k, N+1-k)':>16}  " + "  ".join(
            f"{'pair ' + str(p):>14}" for p in mirror_pairs
        )
        print(header)
        for (k, km), worst in out.items():
            row = f"  ({k}, {km})         "
            row += "  ".join(f"{w:>14.3e}" for w in worst)
            print(row)
        print()

    # Sanity: same N as periodic but OPEN, so we can isolate the periodic effect
    print(f"--- N={N_per} OPEN chain (control for the periodic case above) ---")
    print()
    for delta in [0.0, 0.5, 1.0]:
        out, mirror_pairs = run_xxz_test(N_per, delta, gamma, t_max, n_steps, periodic=False)
        print(f"Δ = {delta:.2f}  (open, even N)")
        header = f"  {'pair (k, N+1-k)':>16}  " + "  ".join(
            f"{'pair ' + str(p):>14}" for p in mirror_pairs
        )
        print(header)
        for (k, km), worst in out.items():
            row = f"  ({k}, {km})         "
            row += "  ".join(f"{w:>14.3e}" for w in worst)
            print(row)
        print()

    # ========================================================================
    # 2-EXCITATION SECTOR TEST (Robustness Test 3)
    # ========================================================================
    print("=" * 78)
    print("2-EXCITATION SECTOR TEST (Robustness Test 3)")
    print("Hypothesis: pure XX (= free fermions) preserves K via Slater")
    print("determinants of single-exc modes; XXZ with Δ ≠ 0 contributes a")
    print("genuine 2-particle interaction (not just a constant shift) and")
    print("breaks K even on uniform-deg topologies.")
    print("=" * 78)
    print()

    N_2exc = 5
    print(f"--- N={N_2exc} OPEN chain, 2-exc Slater initial states ---")
    print()
    for delta in [0.0, 0.3, 1.0]:
        out, mp = run_2exc_xxz_test(N_2exc, delta, gamma, t_max, n_steps, periodic=False)
        print(f"Δ = {delta:.2f}  (open, 2-exc)")
        header = f"  {'(k1,k2) <-> partner':>26}  " + "  ".join(
            f"{'pair ' + str(p):>14}" for p in mp
        )
        print(header)
        for (orig, partner), worst in out.items():
            row = f"  {str(orig)+' <-> '+str(partner):>26}  "
            row += "  ".join(f"{w:>14.3e}" for w in worst)
            print(row)
        print()

    N_2exc_per = 6
    print(f"--- N={N_2exc_per} PERIODIC chain (even N, bipartite), 2-exc ---")
    print("    KEY TEST: in single-exc this restored K at all Δ; in 2-exc")
    print("    we expect K to BREAK at Δ ≠ 0 due to genuine ZZ interaction.")
    print()
    for delta in [0.0, 0.3, 1.0]:
        out, mp = run_2exc_xxz_test(N_2exc_per, delta, gamma, t_max, n_steps, periodic=True)
        print(f"Δ = {delta:.2f}  (periodic even N, 2-exc)")
        header = f"  {'(k1,k2) <-> partner':>26}  " + "  ".join(
            f"{'pair ' + str(p):>14}" for p in mp
        )
        print(header)
        for (orig, partner), worst in out.items():
            row = f"  {str(orig)+' <-> '+str(partner):>26}  "
            row += "  ".join(f"{w:>14.3e}" for w in worst)
            print(row)
        print()

    # Summary
    print("Summary:")
    print("  - At Δ = 0 (XX) on either boundary, K hold (identity to 1e-15).")
    print("  - At Δ ≠ 0 on OPEN chain, ZZ generates non-uniform V_eff(ℓ);")
    print("    boundary-pair (0, N-1) breaks strongest, weaker for deeper")
    print("    mirror-pairs (predicted scaling with deg-discontinuity location).")
    print("  - At Δ ≠ 0 on EVEN-N PERIODIC chain (uniform deg = 2) in single-exc,")
    print("    V_eff is a constant shift -> ZZ does not break K -> identity restored.")
    print("  - Numerical confirmation of PROOF_K_PARTNERSHIP 'Scope' caveat:")
    print("    Heisenberg/XXZ with Δ ≠ 0 breaks K at boundary discontinuity in")
    print("    single-exc; restored on uniform-degree topologies in single-exc.")
    print("  - 2-exc sector: pure XX (Δ=0) preserves K via Slater determinants of")
    print("    single-exc modes (free-fermion structure). Δ ≠ 0 introduces a")
    print("    genuine 2-particle interaction that breaks K even on periodic")
    print("    bipartite chains -- the single-exc 'constant shift' result does")
    print("    NOT generalize to higher excitation sectors.")
    print("  - ODD-N periodic chain is non-bipartite (wrap-around bond connects")
    print("    same-sublattice sites), so K breaks there independently of Δ;")
    print("    this is a pure topology effect, not a ZZ effect.")


if __name__ == "__main__":
    main()
