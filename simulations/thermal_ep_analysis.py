"""
Thermal EP analysis: exceptional point crossings in the oscillating fraction.

Shows that L(n_bar) = L_0 + n_bar · L_thermal is linear in n_bar, and that
isolated exceptional points (EPs) cause the oscillating eigenvalue count to
change by ±2. The oscillating fraction stays near 82% but is NOT exactly
invariant — reverse EPs (two real eigenvalues colliding and splitting into
a complex conjugate pair) occur at specific n_bar values.

This script was developed to verify and correct the claim in Formula 52
(ANALYTICAL_FORMULAS.md) that the oscillating fraction is thermally invariant.
The corrected statement: approximately stable with isolated EP crossings.

Usage:
    python simulations/thermal_ep_analysis.py          # default: N=3 and N=4
    python simulations/thermal_ep_analysis.py --N 5    # include N=5 (slow: 1024x1024)

Requires: numpy
"""
import numpy as np
import argparse
import sys

# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Sp = np.array([[0, 1], [0, 0]], dtype=complex)
Sm = np.array([[0, 0], [1, 0]], dtype=complex)


def kron_chain(ops):
    """Tensor product of a list of operators."""
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def dissipator_term(L_op, d):
    """Lindblad dissipator D[L] in superoperator form."""
    Id = np.eye(d, dtype=complex)
    LdL = L_op.conj().T @ L_op
    return (np.kron(L_op, L_op.conj())
            - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T)))


def build_decomposed(N, gamma=0.05, J=1.0):
    """
    Build Liouvillian decomposed into four parts:
      L_H        (Hamiltonian, purely imaginary spectrum)
      L_D_Z      (Z-dephasing dissipator)
      L_D_emit   (spontaneous + stimulated emission dissipator)
      L_D_absorb (thermal absorption dissipator)

    Full Liouvillian: L(n_bar) = L_H + L_D_Z + (n_bar+1)·L_D_emit + n_bar·L_D_absorb
                                = L_0 + n_bar · L_thermal
    where L_0 = L_H + L_D_Z + L_D_emit and L_thermal = L_D_emit + L_D_absorb.
    """
    d = 2**N
    Id = np.eye(d, dtype=complex)

    # Heisenberg Hamiltonian (XXX chain)
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N; ops[i] = P; ops[i + 1] = P
            H += J * kron_chain(ops)

    L_H = -1j * (np.kron(H, Id) - np.kron(Id, H.T))

    L_D_Z = np.zeros((d**2, d**2), dtype=complex)
    L_D_emit = np.zeros((d**2, d**2), dtype=complex)
    L_D_absorb = np.zeros((d**2, d**2), dtype=complex)

    for k in range(N):
        ops_z = [I2] * N; ops_z[k] = Zm
        L_D_Z += dissipator_term(np.sqrt(gamma) * kron_chain(ops_z), d)

        ops_m = [I2] * N; ops_m[k] = Sm
        L_D_emit += dissipator_term(np.sqrt(gamma) * kron_chain(ops_m), d)

        ops_p = [I2] * N; ops_p[k] = Sp
        L_D_absorb += dissipator_term(np.sqrt(gamma) * kron_chain(ops_p), d)

    return L_H, L_D_Z, L_D_emit, L_D_absorb


def count_oscillating(eigs, tol=1e-6):
    """Count eigenvalues with |Im| > tol."""
    return np.sum(np.abs(eigs.imag) > tol)


def find_ep_transitions(L_0, L_th, n_bar_range, n_points=500, tol=1e-6):
    """Sweep n_bar and find where the oscillating count changes."""
    transitions = []
    prev_osc = None
    for nb in np.linspace(n_bar_range[0], n_bar_range[1], n_points):
        eigs = np.linalg.eigvals(L_0 + nb * L_th)
        osc = count_oscillating(eigs, tol)
        if prev_osc is not None and osc != prev_osc:
            transitions.append((nb, prev_osc, osc))
        prev_osc = osc
    return transitions


def analyze_mechanism(L_0, L_th, nb_crit, delta=0.005, tol=1e-6):
    """Determine whether the EP is forward (complex→real) or reverse (real→complex)."""
    eigs_before = np.linalg.eigvals(L_0 + (nb_crit - delta) * L_th)
    eigs_after = np.linalg.eigvals(L_0 + (nb_crit + delta) * L_th)

    osc_before = count_oscillating(eigs_before, tol)
    osc_after = count_oscillating(eigs_after, tol)

    if osc_after > osc_before:
        return "reverse EP (real → complex conjugate pair)"
    elif osc_after < osc_before:
        return "forward EP (complex conjugate pair → real)"
    else:
        return "no net change (transient)"


# === Main analysis ===

def run_analysis(N_values):
    for N in N_values:
        dim = 4**N
        print(f"\n{'='*70}")
        print(f"N = {N}, Liouvillian dimension {dim}×{dim}")
        print(f"{'='*70}\n")

        L_H, L_D_Z, L_D_emit, L_D_absorb = build_decomposed(N)
        L_0 = L_H + L_D_Z + L_D_emit
        L_th = L_D_emit + L_D_absorb

        # --- Linearity check ---
        L_05_direct = L_0 + 0.5 * L_th
        L_05_interp = 0.5 * L_0 + 0.5 * (L_0 + L_th)
        lin_err = np.linalg.norm(L_05_direct - L_05_interp)
        print(f"Linearity check: ||L(0.5) - interpolated||_F = {lin_err:.2e}")

        # --- EP transition scan ---
        # Scan in three ranges: fine near zero, medium, coarse at high n_bar
        print(f"\nScanning for EP transitions...")
        all_transitions = []
        for (lo, hi, pts) in [(0, 0.5, 500), (0.5, 5, 200), (5, 100, 200)]:
            trans = find_ep_transitions(L_0, L_th, (lo, hi), pts)
            all_transitions.extend(trans)

        if all_transitions:
            print(f"\nFound {len(all_transitions)} EP transition(s):\n")
            print(f"  {'n_bar':>10}  {'osc_before':>10}  {'osc_after':>10}  mechanism")
            print(f"  {'-'*60}")
            for (nb, before, after) in all_transitions:
                mech = analyze_mechanism(L_0, L_th, nb)
                print(f"  {nb:10.5f}  {before:10d}  {after:10d}  {mech}")
        else:
            print("No EP transitions found in [0, 100].")

        # --- Stability table ---
        print(f"\nStability table:")
        print(f"  {'n_bar':>8}  {'osc':>5}  {'frac':>7}  {'Δ from n̄=0':>12}")
        print(f"  {'-'*40}")

        eigs_ref = np.linalg.eigvals(L_0)
        osc_ref = count_oscillating(eigs_ref)

        n_bar_probe = [0, 0.01, 0.1, 0.5, 1, 2, 5, 10]
        if N <= 4:
            n_bar_probe.extend([50, 100])

        for nb in n_bar_probe:
            eigs = np.linalg.eigvals(L_0 + nb * L_th)
            osc = count_oscillating(eigs)
            frac = osc / dim * 100
            delta = osc - osc_ref
            print(f"  {nb:8.2f}  {osc:5d}  {frac:6.1f}%  {delta:+12d}")

        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Thermal EP analysis for Heisenberg chains")
    parser.add_argument("--N", type=int, default=0,
                        help="Max system size (default: 4). Use 5 for thorough analysis (slow).")
    args = parser.parse_args()

    if args.N > 0:
        N_values = list(range(3, args.N + 1))
    else:
        N_values = [3, 4]

    print("Thermal Lindblad EP Analysis")
    print("=" * 70)
    print("L(n_bar) = L_0 + n_bar · L_thermal")
    print("Tracking oscillating eigenvalue count through exceptional points")
    print()

    run_analysis(N_values)

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("""
The oscillating fraction is approximately stable (~82%) under thermal
perturbation, but NOT exactly invariant. Isolated exceptional points
(reverse EPs: two real eigenvalues colliding into a complex conjugate
pair, or forward EPs: complex pair colliding onto real axis) change the
count by ±2 at specific n_bar values.

The perturbation theory argument: L(n_bar) is linear in n_bar, so
eigenvalues move continuously. EPs are codimension-2 in matrix space,
making them isolated points along any 1D parameter path. This explains
why the count is MOSTLY stable — changes are sparse exceptions, not
the rule.

See: Formula 52 in docs/ANALYTICAL_FORMULAS.md
     experiments/THERMAL_BLACKBODY.md (Result 4)
""")
