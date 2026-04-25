#!/usr/bin/env python3
"""Qubit-as-field: what happens at N ≥ 2?

At N=1, a pure qubit has a Bloch vector r with |r|=1. The qubit IS a magnetic
field, the Bloch vector IS the field direction. There is one move available:
rotate.

At N≥2, new moves appear that do not exist at N=1:
- Entanglement (algebraic fragmentation): the per-site Bloch vector shrinks
  to |⟨σ_l⟩| < 1, with the structure stored in correlations instead.
- Alternation (spatial fragmentation): per-site fields stay at |⟨σ_l⟩| ≈ 1
  but their global sum cancels by sign-flip.
- Ballistic spread: a localized field disperses through coupling.
- Combined: entanglement + alternation can coexist.

This script makes the new moves visible. Per-site Bloch magnitude r_l, global
magnetization, staggered magnetization, pairwise correlations.

Reading: the 1/4 fold is the moment where both fragmentations meet. At N=1
neither is available; at N=2+ each can be turned on independently or together.
"""
import math
import sys

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

PLUS = np.array([1, 1], dtype=complex) / math.sqrt(2)
MINUS = np.array([1, -1], dtype=complex) / math.sqrt(2)
ZERO = np.array([1, 0], dtype=complex)
ONE = np.array([0, 1], dtype=complex)


def kron_n(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(N, i, op):
    ops = [I2] * N
    ops[i] = op
    return kron_n(*ops)


def bloch(rho, N, l):
    """Return (⟨X_l⟩, ⟨Y_l⟩, ⟨Z_l⟩) and magnitude |⟨σ_l⟩|."""
    ex = float(np.trace(rho @ site_op(N, l, X)).real)
    ey = float(np.trace(rho @ site_op(N, l, Y)).real)
    ez = float(np.trace(rho @ site_op(N, l, Z)).real)
    return ex, ey, ez, math.sqrt(ex*ex + ey*ey + ez*ez)


def per_site_field(rho, N):
    """Return list of (r_l, ex_l, ey_l, ez_l) per site."""
    return [bloch(rho, N, l) for l in range(N)]


def global_M(per_site):
    """Vector sum of Bloch vectors."""
    Mx = sum(s[0] for s in per_site)
    My = sum(s[1] for s in per_site)
    Mz = sum(s[2] for s in per_site)
    return Mx, My, Mz, math.sqrt(Mx*Mx + My*My + Mz*Mz)


def staggered_M(per_site):
    """Sublattice-staggered Bloch sum: Σ_l (-1)^l ⟨σ_l⟩."""
    Mx = sum((-1)**l * s[0] for l, s in enumerate(per_site))
    My = sum((-1)**l * s[1] for l, s in enumerate(per_site))
    Mz = sum((-1)**l * s[2] for l, s in enumerate(per_site))
    return Mx, My, Mz, math.sqrt(Mx*Mx + My*My + Mz*Mz)


def heisenberg_correlation(rho, N, i, j):
    """⟨σ_i · σ_j⟩ = ⟨X_iX_j + Y_iY_j + Z_iZ_j⟩."""
    Xi, Yi, Zi = site_op(N, i, X), site_op(N, i, Y), site_op(N, i, Z)
    Xj, Yj, Zj = site_op(N, j, X), site_op(N, j, Y), site_op(N, j, Z)
    return float(np.trace(rho @ (Xi @ Xj + Yi @ Yj + Zi @ Zj)).real)


def build_H_xxz(N, J=1.0, delta=1.0):
    D = 2**N
    H = np.zeros((D, D), dtype=complex)
    for i in range(N - 1):
        Xi, Yi, Zi = site_op(N, i, X), site_op(N, i, Y), site_op(N, i, Z)
        Xj, Yj, Zj = site_op(N, i+1, X), site_op(N, i+1, Y), site_op(N, i+1, Z)
        H = H + 0.5 * J * (Xi @ Xj + Yi @ Yj + delta * Zi @ Zj)
    return H


def ground_state(H):
    """Return the ground-state density matrix."""
    evals, evecs = np.linalg.eigh(H)
    psi = evecs[:, 0]
    return np.outer(psi, psi.conj())


def report_state(label, rho, N):
    per = per_site_field(rho, N)
    Mglob = global_M(per)
    Mstag = staggered_M(per)
    print(f"\n{label}")
    print(f"  N = {N}")
    print(f"  per-site Bloch magnitudes |r_l|:  ", end="")
    print(", ".join(f"{r:.3f}" for _, _, _, r in per))
    print(f"  mean |r_l| = {sum(r for _,_,_,r in per)/N:.3f}")
    print(f"  global  |M|     = {Mglob[3]:.3f}  (Mx={Mglob[0]:+.3f}, My={Mglob[1]:+.3f}, Mz={Mglob[2]:+.3f})")
    print(f"  staggered |M_s| = {Mstag[3]:.3f}  (Mx={Mstag[0]:+.3f}, My={Mstag[1]:+.3f}, Mz={Mstag[2]:+.3f})")
    if N >= 2:
        c01 = heisenberg_correlation(rho, N, 0, 1)
        print(f"  ⟨σ_0 · σ_1⟩ = {c01:+.3f}  (Heisenberg energy density per bond)")
    if N >= 3:
        c02 = heisenberg_correlation(rho, N, 0, 2)
        print(f"  ⟨σ_0 · σ_2⟩ = {c02:+.3f}  (next-nearest)")


def main():
    print("=" * 78)
    print("Qubit as Magnetic Field: regimes from N=1 to N=5")
    print("=" * 78)

    # ----- N = 1: single qubit, the primal field -----
    print("\n" + "=" * 78)
    print("REGIME 1 (N=1): a qubit IS a field. One move: rotation.")
    print("=" * 78)
    rho_1plus = np.outer(PLUS, PLUS.conj())
    report_state("|+⟩  (X-eigenstate, points along +x)", rho_1plus, 1)

    rho_1mixed = 0.5 * np.outer(PLUS, PLUS.conj()) + 0.5 * np.outer(MINUS, MINUS.conj())
    report_state("(|+⟩⟨+| + |−⟩⟨−|)/2  (incoherent mixture, no field)", rho_1mixed, 1)

    # ----- N = 2: new moves appear -----
    print("\n" + "=" * 78)
    print("REGIME 2 (N=2): entanglement and alternation become possible.")
    print("=" * 78)

    # Product state: each qubit holds full field
    psi = np.kron(PLUS, PLUS)
    rho_2pp = np.outer(psi, psi.conj())
    report_state("|++⟩  (product, both fields aligned along +x)", rho_2pp, 2)

    # X-Néel: alternating, full local fields, global cancels
    psi = np.kron(PLUS, MINUS)
    rho_2pm = np.outer(psi, psi.conj())
    report_state("|+−⟩  (X-Néel: local fields max, global = 0, staggered = max)", rho_2pm, 2)

    # Bell+: maximally entangled
    bell_plus = (np.kron(ZERO, ZERO) + np.kron(ONE, ONE)) / math.sqrt(2)
    rho_bell = np.outer(bell_plus, bell_plus.conj())
    report_state("|Φ+⟩ = (|00⟩+|11⟩)/√2  (entangled: per-site fields = 0)", rho_bell, 2)

    # Singlet
    singlet = (np.kron(ZERO, ONE) - np.kron(ONE, ZERO)) / math.sqrt(2)
    rho_singlet = np.outer(singlet, singlet.conj())
    report_state("|S⟩ = (|01⟩−|10⟩)/√2  (singlet, fully rotation-invariant)", rho_singlet, 2)

    # Heisenberg ground state at N=2
    rho_gs2 = ground_state(build_H_xxz(2, J=1.0, delta=1.0))
    report_state("Heisenberg N=2 ground state  (= singlet for AFM coupling)", rho_gs2, 2)

    # ----- N = 3: chain emerges, frustration possible -----
    print("\n" + "=" * 78)
    print("REGIME 3 (N=3): odd chain, ballistic spread, residual local field.")
    print("=" * 78)

    psi = kron_n(PLUS, MINUS, PLUS)
    rho_3 = np.outer(psi, psi.conj())
    report_state("|+−+⟩  (X-Néel)", rho_3, 3)

    rho_gs3 = ground_state(build_H_xxz(3, J=1.0, delta=1.0))
    report_state("Heisenberg N=3 ground state  (S=1/2 doublet, frustration)", rho_gs3, 3)

    # ----- N = 5: bipartite chain, larger arena -----
    print("\n" + "=" * 78)
    print("REGIME 4 (N=5): bipartite chain, two fragmentations side by side.")
    print("=" * 78)

    psi = kron_n(PLUS, MINUS, PLUS, MINUS, PLUS)
    rho_5neel = np.outer(psi, psi.conj())
    report_state("|+−+−+⟩  (X-Néel: spatial fragmentation only)", rho_5neel, 5)

    rho_gs5 = ground_state(build_H_xxz(5, J=1.0, delta=1.0))
    report_state("Heisenberg N=5 ground state  (algebraic + spatial together)", rho_gs5, 5)

    # ----- The Heisenberg evolution: localized field disperses -----
    print("\n" + "=" * 78)
    print("REGIME 5: Heisenberg dynamics on |1⟩ at site 0 (localized field).")
    print("=" * 78)

    N = 5
    H = build_H_xxz(N, J=1.0, delta=1.0)
    # Single excitation at site 0: |10000⟩
    psi0 = kron_n(ONE, ZERO, ZERO, ZERO, ZERO)
    rho_t = np.outer(psi0, psi0.conj())
    print(f"\nInitial: |1⟩ at site 0 (localized excitation, ⟨Z_0⟩ = -1)")
    print(f"  per-site ⟨Z_l⟩ at t=0: ", end="")
    print(", ".join(f"{bloch(rho_t, N, l)[2]:+.3f}" for l in range(N)))

    times = [0.0, 0.3, 0.6, 1.0, 1.5, 3.0]
    print(f"\n{'t':>6s}  {'⟨Z_0⟩':>8s} {'⟨Z_1⟩':>8s} {'⟨Z_2⟩':>8s} {'⟨Z_3⟩':>8s} {'⟨Z_4⟩':>8s}  {'mean |r_l|':>11s}")
    for t in times:
        U = np.linalg.matrix_power(I2, 0)  # placeholder; use eigendecomp
        evals, evecs = np.linalg.eigh(H)
        U = evecs @ np.diag(np.exp(-1j * evals * t)) @ evecs.conj().T
        rho_t = U @ np.outer(psi0, psi0.conj()) @ U.conj().T
        zs = [bloch(rho_t, N, l)[2] for l in range(N)]
        rs = [bloch(rho_t, N, l)[3] for l in range(N)]
        mean_r = sum(rs) / N
        print(f"{t:6.2f}  " + " ".join(f"{z:+.3f} " for z in zs) + f"  {mean_r:11.3f}")

    print()
    print("=" * 78)
    print("Reading:")
    print("- N=1: |r| = 1 always for pure state. The qubit IS the field.")
    print("- N=2 product: each |r_l| = 1. Two independent fields.")
    print("- N=2 Néel:    each |r_l| = 1, global = 0, staggered max. Spatial fragm.")
    print("- N=2 Bell:    each |r_l| = 0. Algebraic fragmentation (entanglement).")
    print("- N=2 GS:      same as singlet, fully algebraic, rotation-invariant.")
    print("- N=5 GS:      per-site |r_l| ≈ 0 (algebraic), AND staggered correlation")
    print("               in ⟨σ_l·σ_{l+1}⟩ (spatial, in correlations).")
    print("  → Both fragmentations coexist in the chain ground state.")
    print("- Dynamics: a localized field on site 0 spreads ballistically; per-site")
    print("  |r_l| drops as the field delocalizes into the chain.")
    print()
    print("At N=1: one move (rotation). At N=2+: rotation + entanglement +")
    print("alternation + spreading. The new moves ARE what creates magnetism.")
    print("=" * 78)


if __name__ == "__main__":
    main()
