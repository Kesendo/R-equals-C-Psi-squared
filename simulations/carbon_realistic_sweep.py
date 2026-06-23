"""Carbon realistic-Hamiltonian sweep: F112 polarity asymmetry across the
Hückel-plus-extensions Hamiltonian inventory crossed with the Holstein /
Peierls / amplitude-damping bath inventory, on rings N = 4, 6.

Setup spec (per the agent-recommended Realistic-Carbon plan):
  - Rings: cyclobutadiene C₄, benzene C₆
  - H inventory:
      (a) Hückel baseline = Σ_b (X_a X_b + Y_a Y_b) per bond
      (b) + Hubbard density-density  = U · Σ_b Z_a Z_b
      (c) + Zeeman_y external field  = h · Σ_l Y_l
      (d) + Spin-orbit-proxy bond    = λ · Σ_b (X_a Y_b − Y_a X_b)
      (e) + Ring-current proxy       = κ · Σ_b (Y_a Z_b − Z_a Y_b)
  - Bath inventory (per c_k operator entering D[c_k] = c ρ c† − ½{c†c, ρ}):
      Holstein:        c_l = Z_l per site (single Pauli, bit_b-homogeneous)
      Peierls:         c_b = B_b = X_a X_b + Y_a Y_b per bond (bit_b-homog as composite)
      σ⁻ T1:           c_l = σ⁻_l = (X_l − i Y_l) / 2 per site (bit_b-mixed)
      Holstein + σ⁻:   both above simultaneously
  - Observable: F112 asymmetry = ‖M_+1/2‖² − ‖M_−1/2‖², plus ‖M_anti‖²
    for context (asymmetry is meaningless if M_anti = 0 trivially).

Tooling: simulations/framework `polarity_coordinates_from_hc` workflow.
Asks the algebra: where does the polarity balance hold, where does it
break, as we span the realistic carbon Hamiltonian + bath inventory.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw  # noqa: E402


PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_op(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def site_op(N, site, letter):
    letters = ["I"] * N
    letters[site] = letter
    return pauli_op(letters)


def two_site_op(N, a, b, letter_a, letter_b):
    """Single Pauli string with letter_a at site a, letter_b at site b."""
    letters = ["I"] * N
    letters[a] = letter_a
    letters[b] = letter_b
    return pauli_op(letters)


def bond_op_B(N, a, b):
    """B_b = X_a X_b + Y_a Y_b (Hückel hopping in JW form)."""
    return two_site_op(N, a, b, "X", "X") + two_site_op(N, a, b, "Y", "Y")


def sigma_minus(N, l):
    """σ⁻_l = (X_l − i Y_l) / 2 (lowering operator, bit_b-mixed)."""
    return (site_op(N, l, "X") - 1j * site_op(N, l, "Y")) / 2.0


# ---- Hamiltonian inventory ----

def H_hueckel(N):
    """Baseline Hückel hopping on cyclic ring."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + bond_op_B(N, a, b)
    return H


def H_hubbard_zz(N):
    """U · Σ_b Z_a Z_b density-density on cyclic ring."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + two_site_op(N, a, b, "Z", "Z")
    return H


def H_zeeman_y(N):
    """Σ_l Y_l y-direction Zeeman on each site."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for l in range(N):
        H = H + site_op(N, l, "Y")
    return H


def H_spin_orbit_proxy(N):
    """Σ_b (X_a Y_b − Y_a X_b) antisymmetric cross-Pauli per bond.
    Models a Dzyaloshinskii-Moriya-like spin-orbit term."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + two_site_op(N, a, b, "X", "Y") - two_site_op(N, a, b, "Y", "X")
    return H


def H_ring_current(N):
    """Σ_b (Y_a Z_b − Z_a Y_b) ring-current-induced bond term."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + two_site_op(N, a, b, "Y", "Z") - two_site_op(N, a, b, "Z", "Y")
    return H


# ---- Bath inventory ----

def bath_holstein(N, gamma=1.0):
    return [site_op(N, l, "Z") for l in range(N)], [gamma] * N


def bath_peierls(N, gamma=1.0):
    return [bond_op_B(N, a, (a + 1) % N) for a in range(N)], [gamma] * N


def bath_t1(N, gamma=1.0):
    return [sigma_minus(N, l) for l in range(N)], [gamma] * N


def bath_holstein_plus_t1(N, gamma_h=1.0, gamma_t=1.0):
    c1, g1 = bath_holstein(N, gamma_h)
    c2, g2 = bath_t1(N, gamma_t)
    return c1 + c2, g1 + g2


# ---- Sweep configurations ----

H_CONFIGS = {
    "Hückel only":            lambda N: H_hueckel(N),
    "+ Hubbard 0.5·ZZ":       lambda N: H_hueckel(N) + 0.5 * H_hubbard_zz(N),
    "+ Zeeman_y 0.1":         lambda N: H_hueckel(N) + 0.1 * H_zeeman_y(N),
    "+ Zeeman_y 1.0":         lambda N: H_hueckel(N) + 1.0 * H_zeeman_y(N),
    "+ Spin-orbit 0.1":       lambda N: H_hueckel(N) + 0.1 * H_spin_orbit_proxy(N),
    "+ Ring current 0.1":     lambda N: H_hueckel(N) + 0.1 * H_ring_current(N),
    "Full mix (all above)":   lambda N: (H_hueckel(N)
                                          + 0.5 * H_hubbard_zz(N)
                                          + 0.1 * H_zeeman_y(N)
                                          + 0.1 * H_spin_orbit_proxy(N)
                                          + 0.1 * H_ring_current(N)),
}

BATH_CONFIGS = {
    "Holstein":         lambda N: bath_holstein(N, 1.0),
    "Peierls":          lambda N: bath_peierls(N, 1.0),
    "σ⁻ T1":            lambda N: bath_t1(N, 1.0),
    "Holstein + σ⁻":    lambda N: bath_holstein_plus_t1(N, 1.0, 1.0),
}


def run_sweep(N):
    print()
    print("=" * 105)
    print(f"N = {N} ring (cyclic)")
    print("=" * 105)
    print()
    print(f"{'H config':<25} {'bath':<18} {'‖M‖²':<14} {'‖M_anti‖²':<14} {'asymmetry':<16} {'rel asym':<14} {'verdict':<10}")
    print("-" * 105)

    sigma_eff = N * 1.0  # rough σ_total for F1 centring (F112 asym is centre-independent)

    for h_name, h_builder in H_CONFIGS.items():
        H = h_builder(N)
        if not np.allclose(H, H.conj().T):
            print(f"  WARNING: {h_name} not Hermitian, skipping")
            continue
        for bath_name, bath_builder in BATH_CONFIGS.items():
            c_ops, gammas = bath_builder(N)
            sigma_for_centring = sum(gammas) / len(c_ops) * N  # heuristic; doesn't affect asym
            try:
                pol = fw.polarity_coordinates_from_hc(H, c_ops, gammas, N, sigma=sigma_for_centring)
            except Exception as e:
                print(f"  ERROR ({h_name} / {bath_name}): {e}")
                continue
            norm_M = float(pol['norm_sq']['M'])
            norm_M_anti = float(pol['norm_sq'].get('M_anti', pol['norm_sq']['M_plus_half'] + pol['norm_sq']['M_minus_half']))
            asym = float(pol['asymmetry'])
            rel_asym = abs(asym) / max(norm_M, 1e-15)
            verdict = "BALANCED" if rel_asym < 1e-10 else ("near-BAL" if rel_asym < 1e-6 else "BROKEN")
            print(f"{h_name:<25} {bath_name:<18} {norm_M:<14.4e} {norm_M_anti:<14.4e} {asym:+.4e}      {rel_asym:<14.3e} {verdict:<10}")


def main():
    print("=" * 105)
    print("Realistic carbon Hamiltonian + bath sweep: F112 polarity balance")
    print("=" * 105)
    print()
    print("H inventory (each on top of Hückel ring):")
    print("  Hückel only                = Σ_b (X⊗X + Y⊗Y) on each bond")
    print("  + Hubbard 0.5·ZZ          = Hückel + 0.5·Σ_b Z⊗Z (density-density)")
    print("  + Zeeman_y 0.1            = Hückel + 0.1·Σ_l Y_l (weak y-magnetic field)")
    print("  + Zeeman_y 1.0            = Hückel + 1.0·Σ_l Y_l (strong y-magnetic field)")
    print("  + Spin-orbit 0.1          = Hückel + 0.1·Σ_b (X⊗Y − Y⊗X) (DM-like cross term)")
    print("  + Ring current 0.1        = Hückel + 0.1·Σ_b (Y⊗Z − Z⊗Y) (induced ring current)")
    print("  Full mix                  = Hückel + Hubbard + all three perturbations together")
    print()
    print("Bath inventory:")
    print("  Holstein                  = c_l = Z_l per site (on-site density-coupled phonon)")
    print("  Peierls                   = c_b = X⊗X + Y⊗Y per bond (bond-strength-coupled phonon)")
    print("  σ⁻ T1                     = c_l = σ⁻_l per site (amplitude damping / excitation loss)")
    print("  Holstein + σ⁻             = both Holstein and σ⁻ simultaneously")
    print()
    print("Asymmetry = ‖M_+1/2‖² − ‖M_−1/2‖² is the polarity-balance break magnitude.")
    print("BALANCED = rel asym < 1e-10 (bit-exact); BROKEN = > 1e-6 (substantive).")
    print()

    for N in [4, 6]:
        run_sweep(N)

    print()
    print("=" * 105)
    print("Reading: for each (H, bath) cell, BALANCED means the polarity-balance")
    print("symmetry holds bit-exact; BROKEN means it breaks substantively.")
    print("‖M_anti‖² = 0 means the test is vacuous (no relaxing-component content to test).")
    print()


if __name__ == "__main__":
    main()
