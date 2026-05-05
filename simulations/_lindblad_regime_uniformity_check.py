"""Lindblad-only regime-uniformity check (no Trotter, no T1, no gates, no QPU).

Tom's question after the 2026-05-05 Kingston hardware confirmation: does
the regime-uniformity effect (mixed chains read 13.5× higher truly-baseline
than uniform chains) appear in pure noiseless Lindblad simulation, before
adding any hardware-specific noise mechanisms?

If yes: the effect is structural to per-site γ_l mixing alone. Aer simulation
is not needed; the analytical follow-up is "find a closed-form truly-baseline
as a function of (γ_0, γ_1, γ_2)".

If no: the effect is hardware-specific (T1 amplitude damping, gate errors,
crosstalk, Trotter discretization, or interaction between these). Aer with
controlled per-qubit T1/T2 profiles is the next step.

Method
------
At N = 3, build the F87 trichotomy Hamiltonians (truly XX+YY, soft XY+YX,
mixed XY+YZ, pi2_even_nontruly YZ+ZY). For each Hamiltonian, evolve the
initial state |+−+⟩ continuously under
    L = -i [H, ·] + Σ_l γ_l (Z_l ρ Z_l - ρ)
for t = 0.8 in J = 1 units, with three γ_l profiles representing the
three hardware regimes:
    uniform-classical:  γ_l = (γ_low, γ_low, γ_low)
    uniform-quantum:    γ_l = (γ_high, γ_high, γ_high)
    regime-mixed:       γ_l = (γ_high, γ_low, γ_low)

Then partial-trace over q1 to get reduced ρ on (q0, q2), apply F88-Lens,
read truly-baseline (Π²-odd memory fraction). Compare across regimes.

Hardware anchors (Marrakesh + Kingston, 2026-04-26 / 05-05):
    uniform-classical [48, 49, 50]   truly = 0.0013
    uniform-quantum   [43, 56, 63]   truly = 0.0022
    regime-mixed      [0, 1, 2]      truly = 0.0297
    ratio mixed / uniform-classical = 22.8×
    ratio mixed / uniform-quantum   = 13.5×
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _f88_lens_ibm_framework_snapshots import f88_lens_2qubit, PAULIS

N = 3
T_EVAL = 0.8
J = 1.0


def site_op(letter: str, site: int, n: int) -> np.ndarray:
    """Single-qubit Pauli on `site`, identity elsewhere. Sites 0..n-1 left to right."""
    op = np.array([[1.0]], dtype=complex)
    for k in range(n):
        op = np.kron(op, PAULIS[letter] if k == site else PAULIS["I"])
    return op


def bilinear_hamiltonian(terms: list[tuple[str, str]], j: float = J) -> np.ndarray:
    """H = (J/2) Σ_bond Σ_(α,β) σ^α_l σ^β_{l+1} for given letter-pair list.
    For N=3, two bonds: (0,1) and (1,2)."""
    h = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for bond in [(0, 1), (1, 2)]:
        l, m = bond
        for a, b in terms:
            h += (j / 2.0) * site_op(a, l, N) @ site_op(b, m, N)
    return h


def lindbladian_z(h: np.ndarray, gamma_l: list[float]) -> np.ndarray:
    """L = -i [H, ·] + Σ_l γ_l (Z_l ρ Z_l - ρ) as 4^N × 4^N superoperator on vec(ρ)."""
    d = h.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(h, Id) - np.kron(Id, h.T))
    for l, gamma in enumerate(gamma_l):
        Z_l = site_op("Z", l, N)
        L += gamma * (np.kron(Z_l, Z_l.conj()) - np.kron(Id, Id))
    return L


def evolve(rho0: np.ndarray, h: np.ndarray, gamma_l: list[float], t: float) -> np.ndarray:
    L = lindbladian_z(h, gamma_l)
    vec0 = rho0.flatten("F")
    vec_t = expm(L * t) @ vec0
    return vec_t.reshape(rho0.shape, order="F")


def evolve_trotter(rho0: np.ndarray, h: np.ndarray, gamma_l: list[float],
                   t: float, n_trotter: int) -> np.ndarray:
    """First-order Trotter: alternating unitary H-step and Lindblad dephasing-step.
    Each step has δt = t / n_trotter. n_trotter = 3 matches the run_soft_break.py
    hardware-runs Trotter slicing."""
    dt = t / n_trotter
    d = h.shape[0]
    Id = np.eye(d, dtype=complex)
    U = expm(-1j * h * dt)
    L_dephase = np.zeros((d * d, d * d), dtype=complex)
    for l, gamma in enumerate(gamma_l):
        Z_l = site_op("Z", l, N)
        L_dephase += gamma * (np.kron(Z_l, Z_l.conj()) - np.kron(Id, Id))
    E_dephase = expm(L_dephase * dt)
    rho = rho0.copy()
    for _ in range(n_trotter):
        rho = U @ rho @ U.conj().T
        vec = rho.flatten("F")
        vec = E_dephase @ vec
        rho = vec.reshape(rho.shape, order="F")
    return rho


def xneel_state(n: int) -> np.ndarray:
    """|+−+⟩ for N=3 (alternating |+⟩ and |−⟩ along the X axis)."""
    plus = np.array([1.0, 1.0]) / np.sqrt(2.0)
    minus = np.array([1.0, -1.0]) / np.sqrt(2.0)
    bits = [plus, minus, plus]
    psi = bits[0]
    for k in range(1, n):
        psi = np.kron(psi, bits[k])
    return np.outer(psi, psi.conj())


def partial_trace_q1(rho: np.ndarray) -> np.ndarray:
    """Trace out the middle qubit q1, keep (q0, q2) in standard tensor order."""
    rho_4 = rho.reshape(2, 2, 2, 2, 2, 2)
    return np.einsum("ijkilm->jklm", rho_4).reshape(4, 4)


CATEGORIES = {
    "truly_unbroken":      [("X", "X"), ("Y", "Y")],
    "pi2_odd_pure":        [("X", "Y"), ("Y", "X")],
    "pi2_even_nontruly":   [("Y", "Z"), ("Z", "Y")],
    "mixed_anti_one_sixth":[("X", "Y"), ("Y", "Z")],
}

GAMMA_LOW = 0.05
GAMMA_HIGH = 0.20

PROFILES = {
    "uniform-classical": [GAMMA_LOW] * 3,
    "uniform-quantum":   [GAMMA_HIGH] * 3,
    "regime-mixed":      [GAMMA_HIGH, GAMMA_LOW, GAMMA_LOW],
    "regime-mixed-q1":   [GAMMA_LOW, GAMMA_HIGH, GAMMA_LOW],
    "regime-mixed-q2":   [GAMMA_LOW, GAMMA_LOW, GAMMA_HIGH],
}


def main():
    print("Regime-uniformity check at two simulation levels (no Aer, no QPU)")
    print(f"N = {N}, J = {J}, t = {T_EVAL}, γ_low = {GAMMA_LOW}, γ_high = {GAMMA_HIGH}")
    print("=" * 78)

    rho0 = xneel_state(N)

    print()
    print("=" * 78)
    print("LEVEL 1: continuous Lindblad (exp(L·t), no Trotter, no T1, no gates)")
    print("=" * 78)
    print()
    print(f"  {'category':<22} {'profile':<22} {'Π²-odd / memory':>16}")
    print("  " + "-" * 70)

    results = {}
    for cat, terms in CATEGORIES.items():
        h = bilinear_hamiltonian(terms)
        for profile_name, gamma_l in PROFILES.items():
            rho_t = evolve(rho0, h, gamma_l, T_EVAL)
            rho_q02 = partial_trace_q1(rho_t)
            lens = f88_lens_2qubit(rho_q02)
            results[(cat, profile_name)] = lens["pi2_odd_in_memory"]
            print(f"  {cat:<22} {profile_name:<22} {lens['pi2_odd_in_memory']:>16.6f}")

    print()
    print("=" * 78)
    print("LEVEL 2: first-order Trotter, n_trotter = 3 (matches hardware run, still no T1/gates)")
    print("=" * 78)
    print()
    print(f"  {'category':<22} {'profile':<22} {'Π²-odd / memory':>16}")
    print("  " + "-" * 70)

    trotter_results = {}
    for cat, terms in CATEGORIES.items():
        h = bilinear_hamiltonian(terms)
        for profile_name, gamma_l in PROFILES.items():
            rho_t = evolve_trotter(rho0, h, gamma_l, T_EVAL, n_trotter=3)
            rho_q02 = partial_trace_q1(rho_t)
            lens = f88_lens_2qubit(rho_q02)
            trotter_results[(cat, profile_name)] = lens["pi2_odd_in_memory"]
            print(f"  {cat:<22} {profile_name:<22} {lens['pi2_odd_in_memory']:>16.6f}")

    print()
    print("=" * 78)
    print("READING")
    print()

    print(f"truly-baseline at Level 1 (continuous Lindblad):")
    for prof in ["uniform-classical", "uniform-quantum", "regime-mixed", "regime-mixed-q1", "regime-mixed-q2"]:
        v = results[("truly_unbroken", prof)]
        print(f"  {prof:<25} = {v:.6f}")
    print()
    print(f"truly-baseline at Level 2 (Trotter n=3):")
    for prof in ["uniform-classical", "uniform-quantum", "regime-mixed", "regime-mixed-q1", "regime-mixed-q2"]:
        v = trotter_results[("truly_unbroken", prof)]
        print(f"  {prof:<25} = {v:.6f}")
    print()
    print(f"  hardware anchors (Marrakesh + Kingston, 2026-04-26 / 05-05):")
    print(f"    uniform-classical [48, 49, 50] = 0.0013")
    print(f"    uniform-quantum   [43, 56, 63] = 0.0022")
    print(f"    regime-mixed      [0, 1, 2]    = 0.0297  (~22.8× uniform-classical)")
    print()

    t_uc = trotter_results[("truly_unbroken", "uniform-classical")]
    t_uq = trotter_results[("truly_unbroken", "uniform-quantum")]
    t_m0 = trotter_results[("truly_unbroken", "regime-mixed")]
    print(f"Level 2 mixing penalty (mixed-q0 / uniform-classical):")
    if t_uc > 1e-9:
        print(f"  ratio = {t_m0 / t_uc:.2f}× (hardware: 22.8×)")
    print()

    print("Reading:")
    print("  - Level 1 (continuous Lindblad) gives truly = 0 exactly for all γ profiles,")
    print("    matching the framework's idealised prediction (M = 0 for truly H, regardless")
    print("    of γ structure). Pure γ-mixing alone cannot produce the hardware effect.")
    print("  - Level 2 (Trotter n=3) introduces non-zero truly via the discretisation error.")
    print("    If the Trotter level reproduces the mixing penalty (mixed > uniform by factor")
    print("    of ~10× or more), the effect is Trotter+γ structural and Aer is not needed.")
    print("    If Trotter alone gives a uniform truly across γ profiles, the mixing penalty")
    print("    is hardware-specific (T1, gates, crosstalk, readout) and Aer is the next step.")


if __name__ == "__main__":
    main()
