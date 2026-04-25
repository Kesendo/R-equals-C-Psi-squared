#!/usr/bin/env python3
"""K_dwell scaling vs T1 amplitude damping strength.

Hypothesis (from K-symmetry vs symmetry-breaking discussion 2026-04-25):
- Pure Z-dephasing on Bell+ gives K_dwell/δ = 1.0801 (theoretical, K-symmetric).
- T1 amplitude damping breaks the U(1) excitation-number symmetry. It is a
  "second decay channel" / "second force" beyond pure dephasing.
- Hypothesis: the additional symmetry-breaking accelerates the cusp crossing,
  reducing K_dwell/δ below 1.0801.
- Kingston Hardware (April 2026, data/ibm_cusp_slowing_april2026): K_dwell/δ ≈ 0.67.
  Two pairs at different T1/T2: 0.649 (qubits 124-125) and 0.694 (qubits 14-15).
- If the hypothesis holds, the K_dwell/δ vs (γ_T1 / γ_z) curve should pass
  through 0.67 at the effective Kingston ratio.

Test: scan γ_T1 / γ_z and compute K_dwell/δ for Bell+ via direct Lindblad
propagation. Compare to pure-Z theory (1.0801) and Kingston measurement (0.67).
"""
import math
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# Pauli matrices and ladder operators
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # |0⟩⟨1|
SIGMA_PLUS = SIGMA_MINUS.conj().T                          # |1⟩⟨0|


def kron_n(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(N: int, i: int, op: np.ndarray) -> np.ndarray:
    ops = [I2] * N
    ops[i] = op
    return kron_n(*ops)


def bellplus_2q() -> np.ndarray:
    """|Ψ+⟩ = (|01⟩ + |10⟩)/√2 as a 4×4 density matrix."""
    psi = np.zeros(4, dtype=complex)
    psi[1] = 1.0 / math.sqrt(2)  # |01⟩ at index 1
    psi[2] = 1.0 / math.sqrt(2)  # |10⟩ at index 2
    return np.outer(psi, psi.conj())


def cpsi(rho: np.ndarray) -> float:
    """CΨ = C · Ψ where C = Tr(ρ²), Ψ = l1(ρ)/(d−1).

    l1 is the Baumgratz l1-coherence: sum over OFF-DIAGONAL |ρ_ij| only.
    For Bell+: l1 = |ρ_12| + |ρ_21| = 1, Ψ = 1/3, CΨ(0) = 1·1/3 = 1/3
    (matches F25).
    """
    d = rho.shape[0]
    C = float(np.trace(rho @ rho).real)
    abs_rho = np.abs(rho)
    diag_sum = float(np.sum(np.diag(abs_rho)))
    l1 = float(np.sum(abs_rho)) - diag_sum
    Psi = l1 / (d - 1)
    return C * Psi


def lindblad_rhs(rho, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops):
    """L_D[ρ] = γ_z Σ (Z_l ρ Z_l − ρ)
             + γ_T1 Σ (σ⁻_l ρ σ⁺_l − ½ {σ⁺_l σ⁻_l, ρ}).
    Repo convention for γ_z (see PROOF_C1_MIRROR_SYMMETRY setup):
        L_D = γ · Σ_i (Z_i ρ Z_i − ρ)
    so 2-qubit cross-coherence ρ_{01,10} decays with rate 4γ_z under pure Z,
    matching F25 f(t) = e^{−4γt}.
    """
    out = np.zeros_like(rho)
    for Zl in Z_ops:
        out = out + gamma_z * (Zl @ rho @ Zl - rho)
    if gamma_t1 > 0:
        for sm, sp in zip(sm_ops, sp_ops):
            spsm = sp @ sm
            out = out + gamma_t1 * (sm @ rho @ sp - 0.5 * (spsm @ rho + rho @ spsm))
    return out


def rk4_step(rho, dt, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops):
    k1 = lindblad_rhs(rho, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    k4 = lindblad_rhs(rho + dt * k3, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def K_dwell_over_delta(gamma_z: float, gamma_t1: float,
                       delta: float = 1e-3,
                       t_max: float = 5.0,
                       n_steps: int = 50000) -> float:
    """Compute K_dwell/δ for Bell+ under given Z and T1 dephasing rates.

    K_dwell = γ_z · t_dwell where t_dwell is the time interval during which
    |CΨ(t) − 1/4| < δ. Returns K_dwell/δ. Uses γ_z as the time anchor (the
    "known" rate from Kingston's T2 calibration).
    """
    N = 2
    Z_ops = [site_op(N, i, Z) for i in range(N)]
    sm_ops = [site_op(N, i, SIGMA_MINUS) for i in range(N)]
    sp_ops = [site_op(N, i, SIGMA_PLUS) for i in range(N)]

    rho = bellplus_2q()
    dt = t_max / n_steps
    times = np.empty(n_steps + 1)
    cps = np.empty(n_steps + 1)
    for step in range(n_steps + 1):
        cps[step] = cpsi(rho)
        times[step] = step * dt
        if step < n_steps:
            rho = rk4_step(rho, dt, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)

    in_dwell = np.abs(cps - 0.25) < delta
    if not np.any(in_dwell):
        return float("nan")
    indices = np.where(in_dwell)[0]
    t_dwell = times[indices[-1]] - times[indices[0]]
    return (gamma_z * t_dwell) / delta


def main() -> None:
    print("=" * 76)
    print("K_dwell/δ vs T1 amplitude damping strength")
    print("=" * 76)
    print()
    print("Bell+ on 2 qubits, Z-dephasing γ_z + T1 amplitude damping γ_T1.")
    print("Theory (pure Z, F57): K_dwell/δ = 1.0801, γ-invariant.")
    print("Kingston (April 2026): K_dwell/δ ≈ 0.67 (mean of two pairs).")
    print()

    # γ_z=1 sets the time scale; t_cross for Bell+ under pure Z is at
    # K_cross = γ·t_cross = 0.03735, so t_cross ≈ 0.037. δ = 10⁻² gives a
    # dwell window t_dwell ≈ 1.08·δ/γ = 0.0108. Use t_max ≈ 0.5 to safely
    # bracket the crossing for any T1 contribution; n_steps to resolve dwell.
    gamma_z = 1.0
    delta = 1e-2
    t_max = 0.5
    n_steps = 50000  # dt = 10⁻⁵, ≈ 1000 points across dwell window

    ratios = [0.0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.16, 1.5, 2.0, 2.18, 3.0, 5.0]

    print(f"  {'γ_T1/γ_z':>10}  {'K_dwell/δ':>12}  {'note':>30}")
    print(f"  {'-' * 10}  {'-' * 12}  {'-' * 30}")
    for ratio in ratios:
        gamma_t1 = ratio * gamma_z
        val = K_dwell_over_delta(gamma_z, gamma_t1, delta, t_max, n_steps)
        note = ""
        if abs(ratio - 0.0) < 1e-6:
            note = "pure Z (theory: 1.0801)"
        elif abs(ratio - 1.16) < 1e-3:
            note = "Kingston Pair A (T1=259, T2=150 μs)"
        elif abs(ratio - 2.18) < 1e-3:
            note = "Kingston Pair B (T1=350, T2=381 μs)"
        print(f"  {ratio:>10.2f}  {val:>12.4f}  {note:>30}")
    print()

    print("Kingston measurements:")
    print("  Pair A: K_dwell/δ = 0.649")
    print("  Pair B: K_dwell/δ = 0.694")
    print()
    print("Compare scan above with Kingston values at the matching ratios.")


if __name__ == "__main__":
    main()
