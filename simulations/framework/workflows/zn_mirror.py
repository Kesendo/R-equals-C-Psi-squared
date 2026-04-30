"""Z⊗N-Mirror Symmetrie-Test zwischen zwei Dichtematrizen (transverse-field detection)."""
from __future__ import annotations

import numpy as np


def zn_mirror_diagnostic(chain, rho_a, rho_b, tol=1e-6):
    """Z⊗N-Mirror Symmetrie-Test zwischen zwei Dichtematrizen.

    Wenn ρ_a und ρ_b Z⊗N-Partner sind (ρ_b = (Z⊗N)·ρ_a·(Z⊗N)), gilt für
    jede Pauli-String-Erwartung:

        ⟨P⟩_b = (−1)^n_XY(P) · ⟨P⟩_a

    wobei n_XY(P) = Anzahl X- oder Y-Buchstaben in P. Diese Identität
    scheitert in Anwesenheit von transverse fields h_l X_l oder h_l Y_l;
    sie gilt für XXZ, Z-Dephasing, T1 (σ⁻σ⁺ pairs), und non-uniform
    Z-Detuning (Mini-Magnetfeld δ_l Z_l).

    Liefert max-Abweichung über alle 4^N Pauli-Strings → 'preserved' wenn
    unter tol, sonst 'broken'.

    Args:
        rho_a: 2^N × 2^N Dichtematrix
        rho_b: 2^N × 2^N Dichtematrix (vermuteter Z⊗N-Partner von ρ_a)
        tol: Toleranz für 'preserved' Verdict.

    Returns:
        dict mit 'max_violation' (max Pauli-String-Abweichung),
        'verdict' ('preserved' | 'broken'), 'worst_string' (Pauli-Label
        wo die Abweichung am größten), 'worst_a', 'worst_b' (die zwei
        Erwartungen die nicht stimmen).
    """
    from ..pauli import _PAULI_MATRICES, _k_to_indices, _pauli_label, bit_a
    N = chain.N
    d = 2 ** N
    rho_a = np.asarray(rho_a, dtype=complex)
    rho_b = np.asarray(rho_b, dtype=complex)
    if rho_a.shape != (d, d) or rho_b.shape != (d, d):
        raise ValueError(
            f"rho_a, rho_b must be {d}×{d} for N={N}; got "
            f"{rho_a.shape}, {rho_b.shape}"
        )
    max_violation = 0.0
    worst_k = 0
    worst_a = 0.0
    worst_b = 0.0
    for k in range(4 ** N):
        indices = _k_to_indices(k, N)
        # Build Pauli string σ_indices
        P = _PAULI_MATRICES[indices[0]]
        for idx in indices[1:]:
            P = np.kron(P, _PAULI_MATRICES[idx])
        exp_a = float(np.real(np.trace(rho_a @ P)))
        exp_b = float(np.real(np.trace(rho_b @ P)))
        n_xy = sum(bit_a(idx) for idx in indices)  # X, Y have bit_a=1
        sign = (-1) ** n_xy
        violation = abs(exp_b - sign * exp_a)
        if violation > max_violation:
            max_violation = violation
            worst_k = k
            worst_a = exp_a
            worst_b = exp_b
    worst_label = _pauli_label(worst_k, N)
    return {
        'max_violation': max_violation,
        'verdict': 'preserved' if max_violation < tol else 'broken',
        'worst_string': worst_label,
        'worst_a': worst_a,
        'worst_b': worst_b,
        'tol': tol,
    }
