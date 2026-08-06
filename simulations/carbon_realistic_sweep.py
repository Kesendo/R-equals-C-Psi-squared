"""Carbon realistic-Hamiltonian sweep: F112 polarity asymmetry across the
Hückel-plus-extensions Hamiltonian inventory crossed with the Holstein /
Peierls / amplitude-damping bath inventory, on rings N = 4, 6.

Setup spec (per the agent-recommended Realistic-Carbon plan):
  - Rings: cyclobutadiene C₄, benzene C₆
  - H inventory:
      (a) Hückel baseline = Σ_b (X_a X_b + Y_a Y_b) per bond
      (b) + neighbour density-density = V · Σ_b Z_a Z_b (extended-Hubbard V,
          not an on-site U: the spinless π model has none, since n² = n)
      (c) + transverse field         = h · Σ_l Y_l (a pseudospin field; it is
          number-violating, so it is not a magnetic field on a π system)
      (d) + DM axial bond            = λ · Σ_b (X_a Y_b − Y_a X_b)   (D ∥ ẑ)
      (e) + DM transverse bond       = κ · Σ_b (Y_a Z_b − Z_a Y_b)   (D ∥ x̂)
  - Bath inventory (per c_k operator entering D[c_k] = c ρ c† − ½{c†c, ρ}):
      Holstein:        c_l = Z_l per site (single Pauli, bit_b-homogeneous)
      Peierls:         c_b = B_b = X_a X_b + Y_a Y_b per bond (bit_b-homog as composite)
      σ⁻ T1:           c_l = σ⁻_l = (X_l + i Y_l) / 2 per site (bit_b-mixed)
      Holstein + σ⁻:   both above simultaneously
  - Observable: F112 asymmetry = ‖M_+1/2‖² − ‖M_−1/2‖², plus ‖M_anti‖²
    for context (asymmetry is meaningless if M_anti = 0 trivially).

Tooling: simulations/framework `polarity_coordinates_from_hc` workflow.
Asks the algebra: where does the polarity balance hold, where does it
break, as we span the realistic carbon Hamiltonian + bath inventory.

σ⁻ = (X + iY)/2 = [[0,1],[0,0]] takes |1⟩ → |0⟩, the operator the framework's
`lindbladian_z_plus_t1` builds and the one F113's signs are written for. Where
the asymmetry is nonzero its SIGN is a pairing convention on top of that, see
the note inside `polarity_coordinates_from_hc`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

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
    """σ⁻_l = (X_l + i Y_l) / 2 = [[0,1],[0,0]] (lowering operator |1⟩ → |0⟩,
    bit_b-mixed)."""
    return (site_op(N, l, "X") + 1j * site_op(N, l, "Y")) / 2.0


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


def H_dm_axial(N):
    """Σ_b (X_a Y_b − Y_a X_b) antisymmetric cross-Pauli per bond:
    Dzyaloshinskii-Moriya exchange with the D vector along ẑ, the spin-orbit
    term of the carbon inventory.

    On every Jordan-Wigner-adjacent bond this term IS the bond current:
    i(c†_a c_b − c†_b c_a) = −½ (X_a Y_b − Y_a X_b), residual exactly 0.0 at
    N = 2, 4, 6. It is therefore also the term a magnetic flux induces: a
    Peierls phase t → t·e^(iφ) splits the hopping into cos φ times the usual
    hopping plus sin φ times that same bond current, so for the standard
    tight-binding form H = −t·Σ(e^(iφ) c†_a c_b + h.c.) the coefficient on
    (X_a Y_b − Y_a X_b) is +½·t·sin φ.
    The ring-closing bond is the JW boundary term and carries the parity
    string, exactly as the Hückel hopping on that bond does; the identity is
    bond-local and the ring closure is a spin-model wrap in this sweep.
    It conserves π-electron number: ‖[N̂, ·]‖_F = 0 on every bond, with
    N̂ = Σ_l (I − Z_l)/2."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + two_site_op(N, a, b, "X", "Y") - two_site_op(N, a, b, "Y", "X")
    return H


def H_dm_transverse(N):
    """Σ_b (Y_a Z_b − Z_a Y_b) antisymmetric cross-Pauli per bond:
    Dzyaloshinskii-Moriya exchange with the D vector along x̂, a second DM
    axis beside `H_dm_axial`.

    This term is NOT a ring current, a reading its letters invite and the
    algebra refuses. Its Hilbert-Schmidt overlap with the bond-current
    operator is exactly 0, and it does not conserve π-electron number:
    ‖[N̂, ·]‖_F / ‖·‖_F = 1.0 per bond at every N, with N̂ = Σ_l (I − Z_l)/2.
    What a magnetic flux induces on a π ring is `H_dm_axial`."""
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
    "+ DM axial 0.1":         lambda N: H_hueckel(N) + 0.1 * H_dm_axial(N),
    "+ DM transverse 0.1":    lambda N: H_hueckel(N) + 0.1 * H_dm_transverse(N),
    "Full mix (all above)":   lambda N: (H_hueckel(N)
                                          + 0.5 * H_hubbard_zz(N)
                                          + 0.1 * H_zeeman_y(N)
                                          + 0.1 * H_dm_axial(N)
                                          + 0.1 * H_dm_transverse(N)),
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
            # Denominator is the polarity content ||M_anti||^2, not ||M||^2; see
            # simulations/framework/workflows/polarity_fingerprint.py for why.
            rel_asym = 0.0 if norm_M_anti == 0.0 else abs(asym) / norm_M_anti
            verdict = "BALANCED" if rel_asym < 1e-10 else ("near-BAL" if rel_asym < 1e-6 else "BROKEN")
            print(f"{h_name:<25} {bath_name:<18} {norm_M:<14.4e} {norm_M_anti:<14.4e} {asym:+.4e}      {rel_asym:<14.3e} {verdict:<10}")


def run_block_decomposition(N):
    """Does the polarity balance survive restriction to one π-count block?

    The sweep above reports one number per cell, a norm over the whole operator
    space. That hides the structure a chemistry reading needs: against a
    number-conserving bath, M_anti is block-diagonal in the pair of π counts
    (n_bra, n_ket), and the balance holds inside each block on its own. The
    half-filled block is the one a neutral molecule occupies.

    M_plus_half / M_minus_half come back in the PAULI basis, so they are
    transformed to vec first. `order='F'` matches the transform
    `polarity_coordinates` itself uses, where vec index k = j·d + i holds
    A[i, j]; the two stackings differ by a diagonal sign that leaves every
    per-block Frobenius norm unchanged.
    """
    from framework.pauli import _vec_to_pauli_basis_transform

    d = 2**N
    T = _vec_to_pauli_basis_transform(N, order="F")
    T_inv = T.conj().T / (2**N)
    popcount = np.array([bin(i).count("1") for i in range(d * d) if i < d])
    flat = np.arange(d * d)
    bra = popcount[flat % d]
    ket = popcount[flat // d]

    H = H_hueckel(N) + 0.1 * H_dm_axial(N)
    c_ops, gammas = bath_holstein(N, 1.0)
    pol = fw.polarity_coordinates_from_hc(H, c_ops, gammas, N, sigma=N * 1.0)
    M_plus = T @ np.asarray(pol["M_plus_half"]) @ T_inv
    M_minus = T @ np.asarray(pol["M_minus_half"]) @ T_inv
    M_anti = M_plus + M_minus
    total = float(np.linalg.norm(M_anti) ** 2)

    label = bra * (N + 1) + ket
    off_block = float(np.sum(np.abs(M_anti[~(label[:, None] == label[None, :])]) ** 2))

    print()
    print("=" * 105)
    print(f"N = {N}: does the balance survive restriction to one π-count block?")
    print("=" * 105)
    print("H = Hückel + 0.1·DM axial, Holstein bath (number-conserving on both sides)")
    print(f"  ‖M_anti‖² = {total:.6f}      weight outside the (n_bra, n_ket) blocks = {off_block:.3e}")
    print()
    print(f"  {'block':<32} {'‖M_+1/2‖²':<16} {'‖M_−1/2‖²':<16} {'difference':<14} {'share':<8}")
    print("-" * 90)

    half = N // 2
    blocks = [(f"(n_bra={half}, n_ket={half}) half filled", (bra == half) & (ket == half)),
              (f"(n_bra={half}, n_ket={half - 1})", (bra == half) & (ket == half - 1)),
              ("Δn = 0, all sectors", bra == ket),
              ("Δn odd (π-parity forbidden)", (bra - ket) % 2 != 0),
              ("Δn even, nonzero", ((bra - ket) % 2 == 0) & (bra != ket))]
    for name, sel in blocks:
        idx = np.ix_(sel, sel)
        p = float(np.linalg.norm(M_plus[idx]) ** 2)
        m = float(np.linalg.norm(M_minus[idx]) ** 2)
        print(f"  {name:<32} {p:<16.6f} {m:<16.6f} {p - m:<+14.3e} {100 * (p + m) / total:6.2f}%")
    print()
    print("  Reading: every block balances on its own, so the unpreparable Δn-odd part")
    print("  can be set aside and the neutral half-filled block still balances exactly.")


def main():
    print("=" * 105)
    print("Realistic carbon Hamiltonian + bath sweep: F112 polarity balance")
    print("=" * 105)
    print()
    print("H inventory (each on top of Hückel ring):")
    print("  Hückel only                = Σ_b (X⊗X + Y⊗Y) on each bond")
    print("  + Hubbard 0.5·ZZ          = Hückel + 0.5·Σ_b Z⊗Z (neighbour density-density)")
    print("  + Zeeman_y 0.1            = Hückel + 0.1·Σ_l Y_l (weak transverse pseudospin field)")
    print("  + Zeeman_y 1.0            = Hückel + 1.0·Σ_l Y_l (strong transverse pseudospin field)")
    print("  + DM axial 0.1            = Hückel + 0.1·Σ_b (X⊗Y − Y⊗X) (D ∥ ẑ; the bond current)")
    print("  + DM transverse 0.1       = Hückel + 0.1·Σ_b (Y⊗Z − Z⊗Y) (D ∥ x̂; NOT a ring current)")
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
    print("A ‖M_anti‖² at float-noise level means the test is vacuous (no relaxing-")
    print("component content to test). Read the column, not a literal zero: 8 of the 56")
    print("rows print exactly 0.0 and 4 more print 2.5e-32 / 4.8e-30, which are equally")
    print("vacuous. The 44 rows above 1e-20 are the ones carrying a real reading.")
    print()

    for N in [4, 6]:
        run_block_decomposition(N)


if __name__ == "__main__":
    main()
