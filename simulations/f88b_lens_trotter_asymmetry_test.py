"""Trotter-asymmetry test for F88b-Lens hard category.

OPEN_THREAD_GAMMA0_INFORMATION update 2026-05-18 hypothesis: hard H = XX+XY has
non-commuting bond bilinears (`[X⊗X, X⊗Y] = 2i·X⊗Z`) while soft H = XY+YX has
commuting bond bilinears (`[X⊗Y, Y⊗X] = 0`). Hardware (Marrakesh 2026-04-26)
uses 3 Trotter steps for t=0.8. The unexplained −0.106 hard-Δ in the F88b-Lens
reading (ideal Π²-odd/mem 0.381 vs hardware 0.276) might be Trotter
discretization error acting asymmetrically across the trichotomy.

Test: simulate the same |+−+⟩ evolution under each of truly / soft / hard with
different Trotter step counts (3, 5, 10, 30, 100, ∞), with the same hardware-
calibrated baseline (γ_Z=0.1, h_y=0.05, γ_T1=0.046). Trotter is first-order
split-operator: per step Δt = t/n_trotter, apply each H-term unitary in sequence
then apply the continuous dissipator over Δt.

Expected pattern if hypothesis is right:
- truly: stable across n_trotter (XX,YY commute on the same bond)
- soft: stable across n_trotter (XY,YX commute on the same bond)
- hard: Π²-odd/mem decreases as n_trotter decreases, trending from continuous
  ideal (0.381) toward the hardware value (0.276)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))
from f88b_lens_ibm_framework_snapshots import (  # noqa: E402
    I2, PAULIS, SY, SZ,
    f88b_lens_2qubit,
    lindbladian_plus_minus_plus_rho_q0q2,
)

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def trotter_plus_minus_plus_rho_q0q2(
    t: float,
    H_terms: list[tuple[str, str, float]],
    gamma: float = 0.0,
    h_y: float = 0.0,
    gamma_t1: float = 0.0,
    n_trotter: int = 3,
    N: int = 3,
) -> np.ndarray:
    """First-order split-operator Trotter simulation of |+−+⟩ at N=3.

    Per Trotter step Δt = t/n_trotter:
    - Apply each H-component unitary exp(-i·H_i·Δt) in sequence (bond bilinears
      first by H_terms order, then h_y single-site terms).
    - Apply the continuous dissipator superoperator exp(L_dissipator·Δt) where
      L_dissipator covers Z-dephasing (γ) and T1 amplitude damping (γ_t1).

    Partial trace over q1 → 2-qubit reduced ρ on (q0, q2).
    """
    if N != 3:
        raise NotImplementedError(f"only N=3 supported, got N={N}")

    psi_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi_minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(psi_plus, psi_minus), psi_plus)
    rho_0 = np.outer(psi, psi.conj())

    def site_op(P, k):
        ops = [I2] * N
        ops[k] = P
        return np.kron(np.kron(ops[0], ops[1]), ops[2])

    d = 2 ** N
    I_d = np.eye(d, dtype=complex)

    # Build per-component H operators (bilinear bond terms + h_y single-site)
    H_components = []
    for (a_letter, b_letter, coeff) in H_terms:
        Pa = PAULIS[a_letter]
        Pb = PAULIS[b_letter]
        for b in range(N - 1):
            H_components.append(coeff * site_op(Pa, b) @ site_op(Pb, b + 1))
    if h_y != 0.0:
        for k in range(N):
            H_components.append(h_y * site_op(SY, k))

    # Build continuous dissipator superoperator
    L_diss_super = np.zeros((d * d, d * d), dtype=complex)
    if gamma > 0.0:
        for l in range(N):
            Z_l = site_op(SZ, l)
            L_diss_super += gamma * (np.kron(Z_l.T, Z_l) - np.kron(I_d, I_d))
    if gamma_t1 > 0.0:
        sigma_minus = np.array([[0, 1], [0, 0]], dtype=complex)
        P1_single = np.array([[0, 0], [0, 1]], dtype=complex)
        for l in range(N):
            sm_l = site_op(sigma_minus, l)
            P1_l = site_op(P1_single, l)
            L_diss_super += gamma_t1 * (
                np.kron(sm_l.conj(), sm_l)
                - 0.5 * np.kron(I_d, P1_l)
                - 0.5 * np.kron(P1_l.T, I_d)
            )

    # Single Trotter step: unitary product, then continuous dissipator
    dt = t / n_trotter
    U_step = np.eye(d, dtype=complex)
    for H_op in H_components:
        U_step = expm(-1j * H_op * dt) @ U_step
    # Super-op: vec(U·ρ·U†) = (U* ⊗ U) vec(ρ)
    U_super = np.kron(U_step.conj(), U_step)
    if gamma > 0.0 or gamma_t1 > 0.0:
        E_super = expm(L_diss_super * dt)
        step_super = E_super @ U_super
    else:
        step_super = U_super

    vec_rho = rho_0.flatten(order='F')
    full_super = np.linalg.matrix_power(step_super, n_trotter)
    vec_rho_t = full_super @ vec_rho
    rho_full = vec_rho_t.reshape((d, d), order='F')

    rho_q0q2 = np.zeros((4, 4), dtype=complex)
    for i0 in range(2):
        for i2 in range(2):
            for j0 in range(2):
                for j2 in range(2):
                    s = 0.0 + 0j
                    for q1 in range(2):
                        idx_i = (i0 << 2) | (q1 << 1) | i2
                        idx_j = (j0 << 2) | (q1 << 1) | j2
                        s += rho_full[idx_i, idx_j]
                    rho_q0q2[(i0 << 1) | i2, (j0 << 1) | j2] = s
    return rho_q0q2


HARD_H_TERMS = [("X", "X", 1.0), ("X", "Y", 1.0)]
SOFT_H_TERMS = [("X", "Y", 1.0), ("Y", "X", 1.0)]
TRULY_H_TERMS = [("X", "X", 1.0), ("Y", "Y", 1.0)]

GAMMA = 0.1
H_Y = 0.05
GAMMA_T1 = 0.046
T = 0.8

# Hardware values from Marrakesh 2026-04-26 framework_snapshots
HW_VALUES = {"truly": 0.0297, "soft": 0.7444, "hard": 0.2763}


def main() -> None:
    print("Trotter-asymmetry test for F88b-Lens hard category")
    print("=" * 78)
    print(f"Setup: N=3, |+−+⟩, t={T}, γ_Z={GAMMA}, h_y={H_Y}, γ_T1={GAMMA_T1}")
    print(f"Hardware (Marrakesh 2026-04-26, 3 Trotter steps from JSON):")
    print(f"  truly={HW_VALUES['truly']:.4f}  soft={HW_VALUES['soft']:.4f}  hard={HW_VALUES['hard']:.4f}")
    print()

    trotter_counts = [3, 5, 10, 30, 100]

    print(f"  {'category':>8s}  {'cont.':>7s}   " +
          "  ".join(f"n={n:>3d}" for n in trotter_counts) +
          f"   {'hw':>7s}   {'Δ(n=3 − cont.)':>14s}")
    print("  " + "-" * 84)

    for label, H_terms in [("truly", TRULY_H_TERMS), ("soft", SOFT_H_TERMS), ("hard", HARD_H_TERMS)]:
        rho_continuous = lindbladian_plus_minus_plus_rho_q0q2(
            t=T, H_terms=H_terms, gamma=GAMMA, h_y=H_Y, gamma_t1=GAMMA_T1,
        )
        val_continuous = f88b_lens_2qubit(rho_continuous)["pi2_odd_in_memory"]

        trotter_vals = []
        for n in trotter_counts:
            rho_trotter = trotter_plus_minus_plus_rho_q0q2(
                t=T, H_terms=H_terms, gamma=GAMMA, h_y=H_Y, gamma_t1=GAMMA_T1, n_trotter=n,
            )
            v = f88b_lens_2qubit(rho_trotter)["pi2_odd_in_memory"]
            trotter_vals.append(v)

        delta_n3_vs_cont = trotter_vals[0] - val_continuous
        print(f"  {label:>8s}  {val_continuous:>7.4f}   " +
              "  ".join(f"{v:>5.4f}" for v in trotter_vals) +
              f"   {HW_VALUES[label]:>7.4f}   {delta_n3_vs_cont:>+14.4f}")

    print()
    print("Reading rule:")
    print("  - If hard's n=3 Trotter value sits BELOW continuous AND trends toward hardware")
    print("    as n decreases, the −0.106 hardware gap is (at least partly) Trotter-asymmetry.")
    print("  - Truly and soft should stay close to continuous across all n (commuting bilinears).")
    print("  - n=100 should match continuous to within Trotter error O(1/n²) ~ 0.0001.")


if __name__ == "__main__":
    main()
