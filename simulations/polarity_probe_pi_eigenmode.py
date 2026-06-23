"""Polarity probe #6: L hand-engineered as pure Π-conjugation +i eigenmode.

This is the candidate-breaker the POLARITY_COORDINATES reflection asked for:
build L by a path that does not factor through any (kron(c, c.conj()) +- something)
form, with deliberate +i Pi-eigenvalue imbalance hand-engineered into M_anti.

Approach:
1. Diagonalize Pi (order-4 unitary, eigenvalues {+1, -1, +i, -i}).
2. In Pi's eigenbasis, the matrix element L[i,j] has Pi-conjugation eigenvalue
   eigvals[i] / eigvals[j].
3. Build L with support ONLY on (i, j) pairs whose ratio = +i. This makes L
   a pure Pi-conjugation +i eigenmode: Pi @ L @ Pi^-1 = i * L.
4. Set sigma = 0 (skip the dephasing shift; it would only add to M_zero).
5. Feed L to polarity_coordinates_from_L. If the formulas are honest,
   M = Pi L Pi^-1 + L = (1 + i) L is also Pi-conj +i eigenmode, so:
      M_anti = (M - i*M) / 2 = (1-i)/2 * M
      M_plus_half = (M_anti - i*Pi*M_anti*Pi^-1) / 2 = (1-i)/2 * M  (kept)
      M_minus_half = (M_anti + i*Pi*M_anti*Pi^-1) / 2 = 0          (killed)
   So predicted asymmetry = ||M||^2 / 2, NON-ZERO.

If observed asymmetry is non-zero: balance broken. The structural identity is
"the standard L construction channels (kron(c, c.conj())) all give M with
balanced +/-i Pi content; outside those channels the balance can break".

If observed asymmetry is zero anyway: the balance is deeper than even this
hand-engineered breaker can violate. That would suggest an algebraic identity
INSIDE polarity_coordinates_from_L itself, independent of L.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np

from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L


def main():
    N = 3
    d2 = 4 ** N

    Pi = build_pi_full(N)

    # Pi is unitary order-4. Diagonalize.
    eigvals_raw, U = np.linalg.eig(Pi)

    # Pi is unitary so eigvals are on the unit circle. Snap to {+1, -1, +i, -i}.
    eigvals = np.zeros_like(eigvals_raw)
    for k, lam in enumerate(eigvals_raw):
        candidates = np.array([1.0, -1.0, 1j, -1j])
        eigvals[k] = candidates[np.argmin(np.abs(candidates - lam))]

    # Pi eigenspace dimensions
    dims = {lam: int(np.sum(eigvals == lam)) for lam in [1.0, -1.0, 1j, -1j]}
    print(f"Pi eigenvalue multiplicities at N={N}, d^2={d2}:")
    print(f"  +1: {dims[1.0]:>4d}")
    print(f"  -1: {dims[-1.0]:>4d}")
    print(f"  +i: {dims[1j]:>4d}")
    print(f"  -i: {dims[-1j]:>4d}")
    print(f"  total: {sum(dims.values())} (expect {d2})")
    print()

    # Verify U diagonalizes Pi
    U_inv = np.linalg.inv(U)
    Pi_reconstructed = U @ np.diag(eigvals) @ U_inv
    diag_residual = np.max(np.abs(Pi - Pi_reconstructed))
    print(f"Pi reconstruction residual (should be ~ 1e-12): {diag_residual:.3e}")
    print()

    # Pi-conjugation eigenvalue matrix: ratios[i, j] = eigvals[i] / eigvals[j]
    ratios = eigvals[:, None] / eigvals[None, :]

    # Mask for pure Pi-conj +i eigenmode entries
    plus_i_mask = np.abs(ratios - 1j) < 1e-6
    minus_i_mask = np.abs(ratios + 1j) < 1e-6

    n_plus_i = int(plus_i_mask.sum())
    n_minus_i = int(minus_i_mask.sum())
    print(f"Pi-conj +i eigenspace size: {n_plus_i} (i.e., {n_plus_i} (i,j) pairs)")
    print(f"Pi-conj -i eigenspace size: {n_minus_i}")
    print(f"  (should be equal by Pi's order-4 symmetry, see Pi pair-count formula)")
    print()

    # Build L in U-basis with random complex values ONLY on +i entries
    rng = np.random.default_rng(seed=42)
    L_ubasis = np.zeros((d2, d2), dtype=complex)
    n_fill = n_plus_i
    L_ubasis[plus_i_mask] = (
        rng.normal(size=n_fill) + 1j * rng.normal(size=n_fill)
    )

    # Transform L back to Pauli basis (this is the L_pauli that polarity_coordinates_from_L expects)
    L_pauli = U @ L_ubasis @ U_inv

    # Sanity check: Pi @ L_pauli @ Pi^-1 should equal i * L_pauli (Pi-conj +i eigenmode)
    Pi_inv = Pi.conj().T  # unitary
    Pi_L_Pi_inv = Pi @ L_pauli @ Pi_inv
    eigenmode_residual = np.max(np.abs(Pi_L_Pi_inv - 1j * L_pauli))
    print(f"Pi-conj +i eigenmode residual ||Pi L Pi^-1 - i L||_max: {eigenmode_residual:.3e}")
    print(f"  (should be ~ 1e-12 if L is correctly built as pure +i eigenmode)")
    print()

    # Run polarity_coordinates_from_L with sigma = 0
    result = polarity_coordinates_from_L(L_pauli, N=N, sigma=0.0)

    print("=" * 70)
    print("Polarity decomposition of hand-engineered Pi-conj +i eigenmode L")
    print("=" * 70)
    print(f"  ||M||^2          = {result['norm_sq']['M']:.6f}")
    print(f"  ||M_zero||^2     = {result['norm_sq']['M_zero']:.6f}")
    print(f"  ||M_plus_half||^2  = {result['norm_sq']['M_plus_half']:.6f}")
    print(f"  ||M_minus_half||^2 = {result['norm_sq']['M_minus_half']:.6f}")
    print(f"  asymmetry        = {result['asymmetry']:+.6e}")
    print(f"  ortho residual   = {result['orthogonality_residual']:.3e}")
    print()

    predicted = result['norm_sq']['M'] / 2
    print(f"PREDICTION (if pure +i eigenmode): asymmetry = ||M||^2 / 2 = {predicted:.6f}")
    print(f"OBSERVED:                          asymmetry = {result['asymmetry']:+.6f}")
    print()

    if abs(result['asymmetry']) > 1e-6:
        print("*** BALANCE BROKEN ***")
        print("Hand-engineered Pi-conj +i eigenmode L produces asymmetric +/-1/2.")
        print("The standard L construction pipeline (kron(c, c.conj())) was the source")
        print("of the +/-1/2 balance; the balance is NOT an algebraic identity of")
        print("polarity_coordinates_from_L itself.")
    else:
        print("*** BALANCE PRESERVED EVEN HERE ***")
        print("Even hand-engineered Pi-conj +i eigenmode L gives asymmetry = 0.")
        print("The balance IS an algebraic identity of polarity_coordinates_from_L")
        print("or of the way Pi acts on Liouville space, independent of L.")
        print("Deeper structural investigation needed.")


if __name__ == '__main__':
    main()
