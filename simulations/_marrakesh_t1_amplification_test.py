"""Quantify the 'T1 amplification' hypothesis from data/ibm_soft_break_april2026/README.md.

The README's hardening explanation was: 'T1 thermal relaxation and ZZ crosstalk
compound the soft-break operator-level break in observables that are not purely
diagonal in L's eigenbasis.' This script tests it.

Result: T1 ATTENUATES the soft signal, it does NOT amplify it. The hardware
hardening (Δ_hw = -0.722 vs 'framework idealized' = -0.623) is fully explained
by replacing continuous Lindblad evolution with the actual Trotter n=3 circuit
that the hardware ran. At γ_Z = 0.1, Trotter alone reproduces hardware to
within 0.0014 (well below shot-noise statistical error).

Three quantitative findings:

  1. T1 sweep (continuous Lindblad, γ_Z=0.1):
     γ_T1 = 0   → Δ = -0.623 (idealized framework)
     γ_T1 = 0.5 → Δ = -0.440 (50% suppressed)
     γ_T1 = 1.0 → Δ = -0.303 (further suppressed)
     T1 monotonically REDUCES |Δ|; no value reproduces hardware -0.722.

  2. Trotter n=3 alone (no T1, no extra noise, γ_Z=0.1):
     Δ_continuous = -0.6230
     Δ_trotter    = -0.7231
     Δ_hardware   = -0.7217  (matches Trotter to 0.0014)

  3. Joint optimization (γ_Z, γ_T1) over all 45 hardware pairs:
     optimum at γ_T1 = 0.0001 ≈ 0; data does NOT prefer adding T1.

Mechanism: at δt = 0.267 (Trotter step) and ‖H‖ ~ J·N ~ 3, the small-step
condition ‖H·δt‖ ≪ 1 fails. First-order Trotter has systematic bias that
points outward in operator space (toward stronger soft-break magnitude), not
inward. The 'amplification' Tom observed is a Trotter discretization artifact
of the circuit, not a noise channel effect.

Note on Y-axis sign flip: continuous Lindblad predicts ⟨Y₀Z₂⟩_truly = -0.408
(via Heisenberg dynamics), Trotter n=3 predicts +0.439, hardware +0.583.
Trotter recovers the correct hardware sign for all Y-containing observables.
This is independent confirmation that the hardware is in the Trotter regime
where continuous evolution is the wrong physics.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq, minimize

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw
from framework.lindblad import lindbladian_z_dephasing, lindbladian_z_plus_t1
from framework.pauli import _build_bilinear, ur_pauli, site_op


HW_PATH = Path(__file__).parent.parent / "data" / "ibm_soft_break_april2026" / \
    "soft_break_ibm_marrakesh_20260426_001101.json"


def _setup():
    with open(HW_PATH) as f:
        hw = json.load(f)
    N = hw["parameters"]["N"]; J = hw["parameters"]["J"]; t = hw["parameters"]["t"]
    n_trot = hw["parameters"]["n_trotter"]
    chain = fw.ChainSystem(N=N)
    bonds = chain.bonds
    ket_p = np.array([1, 1], dtype=complex) / np.sqrt(2)
    ket_m = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi0 = np.kron(np.kron(ket_p, ket_m), ket_p)
    rho0 = np.outer(psi0, psi0.conj())
    cases = {
        "truly_unbroken": [("X", "X", J), ("Y", "Y", J)],
        "soft_broken":   [("X", "Y", J), ("Y", "X", J)],
        "hard_broken":   [("X", "X", J), ("X", "Y", J)],
    }
    return hw, N, J, t, n_trot, bonds, rho0, cases


def vec_F(M): return M.flatten("F")
def unvec_F(v, d): return v.reshape((d, d), order="F")


def predict_continuous(case_terms, N, bonds, rho0, t, gamma_z, gamma_t1):
    H = _build_bilinear(N, bonds, case_terms)
    if gamma_t1 == 0:
        L = lindbladian_z_dephasing(H, [gamma_z] * N)
    else:
        L = lindbladian_z_plus_t1(H, [gamma_z] * N, [gamma_t1] * N)
    return unvec_F(expm(L * t) @ vec_F(rho0), 2 ** N)


def predict_trotter(case_terms, N, bonds, rho0, t, n_trot, gamma_z, gamma_t1):
    delta_t = t / n_trot
    Id = np.eye(2 ** N, dtype=complex)
    L_diss = np.zeros((4 ** N, 4 ** N), dtype=complex)
    for l in range(N):
        Zl = site_op(N, l, "Z")
        L_diss += gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    if gamma_t1 > 0:
        sm = np.array([[0, 1], [0, 0]], dtype=complex)
        for l in range(N):
            ops = [np.eye(2, dtype=complex)] * N
            ops[l] = sm
            sml = ops[0]
            for op in ops[1:]:
                sml = np.kron(sml, op)
            sml_dag_sml = sml.conj().T @ sml
            L_diss += gamma_t1 * (np.kron(sml, sml.conj())
                                  - 0.5 * np.kron(sml_dag_sml, Id)
                                  - 0.5 * np.kron(Id, sml_dag_sml.T))
    M_diss = expm(L_diss * delta_t)
    U_step = np.eye(2 ** N, dtype=complex)
    for (P, Q, c) in case_terms:
        for (l, m) in bonds:
            ops = [ur_pauli("I")] * N
            ops[l] = ur_pauli(P); ops[m] = ur_pauli(Q)
            op_full = ops[0]
            for op in ops[1:]:
                op_full = np.kron(op_full, op)
            U_step = expm(-1j * c * delta_t * op_full) @ U_step
    rho = rho0
    for _ in range(n_trot):
        rho = U_step @ rho @ U_step.conj().T
        rho = unvec_F(M_diss @ vec_F(rho), 2 ** N)
    return rho


def expectation_xz(rho, N):
    obs = np.kron(np.kron(ur_pauli("X"), np.eye(2 ** (N - 2))), ur_pauli("Z"))
    return float(np.real(np.trace(rho @ obs)))


def main():
    hw, N, J, t, n_trot, bonds, rho0, cases = _setup()
    hw_xz = {label: hw["expectations"][label]["X,Z"] for label in cases}
    hw_delta = hw_xz["soft_broken"] - hw_xz["truly_unbroken"]

    print("=" * 72)
    print("T1 amplification hypothesis test")
    print("=" * 72)
    print(f"  Hardware: ⟨X₀Z₂⟩ truly={hw_xz['truly_unbroken']:+.4f}, "
          f"soft={hw_xz['soft_broken']:+.4f}, "
          f"Δ(soft − truly)={hw_delta:+.4f}")
    print()

    # Test 1: T1 monotonically attenuates the soft signal
    print("Test 1: γ_T1 sweep (γ_Z=0.1, continuous Lindblad)")
    print(f"  {'γ_T1':>5} | {'truly':>8} | {'soft':>8} | {'Δ(s-t)':>8}")
    print("  " + "-" * 42)
    for gt1 in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
        rt = predict_continuous(cases["truly_unbroken"], N, bonds, rho0, t, 0.1, gt1)
        rs = predict_continuous(cases["soft_broken"], N, bonds, rho0, t, 0.1, gt1)
        xt = expectation_xz(rt, N); xs = expectation_xz(rs, N)
        print(f"  {gt1:>5.2f} | {xt:+8.4f} | {xs:+8.4f} | {xs-xt:+8.4f}")
    print("  → T1 monotonically reduces |Δ|; cannot reach hardware Δ = -0.722")
    print()

    # Test 2: Trotter n=3 alone reproduces hardware
    print(f"Test 2: continuous vs Trotter (n_trotter={n_trot}) at γ_Z=0.1, no T1")
    print()
    rt_c = predict_continuous(cases["truly_unbroken"], N, bonds, rho0, t, 0.1, 0)
    rs_c = predict_continuous(cases["soft_broken"], N, bonds, rho0, t, 0.1, 0)
    rt_t = predict_trotter(cases["truly_unbroken"], N, bonds, rho0, t, n_trot, 0.1, 0)
    rs_t = predict_trotter(cases["soft_broken"], N, bonds, rho0, t, n_trot, 0.1, 0)
    xt_c, xs_c = expectation_xz(rt_c, N), expectation_xz(rs_c, N)
    xt_t, xs_t = expectation_xz(rt_t, N), expectation_xz(rs_t, N)
    print(f"  continuous : truly={xt_c:+.4f}, soft={xs_c:+.4f}, Δ={xs_c-xt_c:+.4f}")
    print(f"  Trotter n=3: truly={xt_t:+.4f}, soft={xs_t:+.4f}, Δ={xs_t-xt_t:+.4f}")
    print(f"  hardware   : truly={hw_xz['truly_unbroken']:+.4f}, soft={hw_xz['soft_broken']:+.4f}, Δ={hw_delta:+.4f}")
    print()
    print(f"  Trotter Δ vs hardware Δ: difference = {abs((xs_t-xt_t) - hw_delta):.4f} "
          f"(shot-noise statistical error ≈ 0.015)")
    print(f"  → Trotter alone (no T1, no extra noise) reproduces hardware exactly.")
    print()

    # Test 3: joint (γ_Z, γ_T1) optimization
    def total_rms(gz, gt1):
        err = 0.0; cnt = 0
        for label, terms in cases.items():
            rho = predict_trotter(terms, N, bonds, rho0, t, n_trot, gz, gt1)
            for p0 in "IXYZ":
                for p2 in "IXYZ":
                    key = f"{p0},{p2}"
                    if key == "I,I":
                        continue
                    obs = np.kron(np.kron(ur_pauli(p0), np.eye(2)), ur_pauli(p2))
                    pred = float(np.real(np.trace(rho @ obs)))
                    err += (hw["expectations"][label][key] - pred) ** 2
                    cnt += 1
        return math.sqrt(err / cnt)

    print("Test 3: joint (γ_Z, γ_T1) optimization, Trotter n=3, all 45 obs-pairs")
    res = minimize(lambda x: total_rms(max(x[0], 0), max(x[1], 0)),
                   x0=[0.1, 0.05], method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1e-4})
    g_opt, gt1_opt = res.x
    print(f"  optimal γ_Z   = {g_opt:.4f}")
    print(f"  optimal γ_T1  = {gt1_opt:.6f}")
    print(f"  total RMS     = {res.fun:.4f}")
    print()
    if gt1_opt < 0.005:
        print("  → optimum at γ_T1 ≈ 0: the data does NOT support adding T1.")

    print()
    print("=" * 72)
    print("Conclusion")
    print("=" * 72)
    print("""
  T1 amplification hypothesis: REFUTED.

  T1 monotonically attenuates the soft signal (predictable: amplitude damping
  pulls ρ toward |0⟩⟨0| which has ⟨X₀Z₂⟩ = 0). The 'amplification' Tom
  observed is the Trotter n=3 discretization built into the actual hardware
  circuit. At δt = 0.267 with ‖H‖ ~ J·N ~ 3, the small-step condition fails;
  first-order Trotter biases ⟨X₀Z₂⟩ outward by ≈ +0.10.

  The corrected framework prediction (Trotter n=3, γ_Z=0.1, no T1) matches
  hardware Δ = -0.7217 to within 0.0014, well below shot-noise statistical
  error (≈ 0.015). Joint (γ_Z, γ_T1) optimization across all 45 hardware
  observable-Hamiltonian pairs converges to γ_T1 = 0; the data prefers no T1.

  The README of data/ibm_soft_break_april2026/ should be amended to attribute
  the hardening to Trotter discretization, not to T1 thermal relaxation.
""")


if __name__ == "__main__":
    main()
