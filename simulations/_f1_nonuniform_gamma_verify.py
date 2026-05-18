"""Verify the F1 H-block residual is γ-independent under site-dependent Z-dephasing.

Closes the F1 OpenQuestion "non-uniform γ_i: site-dependent dephasing" by a NEGATIVE
result: the conjectured Σγ_l² replacement of (Σγ)² does not occur, because the
dissipator-block residual M_D collapses to zero PER PAULI STRING under the F1 σ-shift
`2σ·I = 2·Σ_l γ_l · I`, regardless of how γ is distributed across sites.

This script verifies the closure derived in
`docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md`:

(1) H-block γ-independence at N=3 with a soft (Π²-odd) Hamiltonian XY+YX, across
    three γ patterns (uniform, scaled, random): ‖M‖² is bit-exact identical.
(2) Dissipator-block per-Pauli-string vanishing: with H=0, ‖M‖² ≈ 0 at machine
    precision for non-uniform γ at N=3, 4, 5.
(3) Scaling against c_H · F(N, G): at N=3, 4, 5 with a fixed H = XX+YZ chain
    (main Hamiltonian class, c_H anchored at N=2), the residual matches
    c_H · (N−1) · 4^(N−2) for uniform and non-uniform γ alike.
(4) Cross-Hamiltonian invariance: at N=3, 4, 5 across multiple H classes
    (soft XY+YX, main mixed XX+YZ) and γ patterns, every γ pattern at fixed
    (N, H) gives the same ‖M‖².

All assertions are bit-exact (absolute deviation = 0.0 between γ patterns,
operator-level identity from M_D = 0 per Pauli string).
"""
from __future__ import annotations

import sys

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from framework.lindblad import (  # noqa: E402
    lindbladian_z_dephasing,
    palindrome_residual,
)
from framework.pauli import _build_bilinear  # noqa: E402


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #


def m_norm_squared(N: int, H: np.ndarray, gamma_l: list[float]) -> float:
    """‖M‖²_F for the F1 residual M = Π·L·Π⁻¹ + L + 2σ·I, with σ = Σ_l γ_l."""
    L = lindbladian_z_dephasing(H, gamma_l)
    M = palindrome_residual(L, sum(gamma_l), N)
    return float(np.real(np.trace(M.conj().T @ M)))


def soft_xy_yx_chain(N: int) -> np.ndarray:
    """Soft (Π²-odd) Hamiltonian XY+YX on a chain (non-zero ‖M_H‖²; tests γ-independence)."""
    bonds = [(b, b + 1) for b in range(N - 1)]
    return _build_bilinear(N, bonds, [("X", "Y", 1.0), ("Y", "X", 1.0)])


def main_xx_yz_chain(N: int) -> np.ndarray:
    """Main (Π²-mixed) Hamiltonian XX+YZ on a chain; matches the anchor of
    PalindromeResidualScalingClaim.Verify (XX+YZ at N=2 gives c_H)."""
    bonds = [(b, b + 1) for b in range(N - 1)]
    return _build_bilinear(N, bonds, [("X", "X", 1.0), ("Y", "Z", 1.0)])


def gamma_patterns(N: int, seed: int = 42) -> list[tuple[str, list[float]]]:
    """Three γ patterns at fixed N: uniform, scaled, random."""
    rng = np.random.RandomState(seed + N)
    return [
        ("uniform", [0.1] * N),
        ("scaled",  [0.05 * (l + 1) for l in range(N)]),
        ("random",  [float(g) for g in rng.uniform(0.01, 0.20, N)]),
    ]


# --------------------------------------------------------------------------- #
# Sections                                                                    #
# --------------------------------------------------------------------------- #


def section_1_h_block_gamma_independence() -> None:
    """At N=3 with soft H = XY+YX, ‖M‖² is bit-exact identical across γ patterns."""
    print("(1) H-block γ-independence: at N=3 with soft H = XY+YX (Π²-odd),")
    print("    ‖M‖² is bit-exact identical for uniform, scaled, and random γ.")
    N = 3
    H = soft_xy_yx_chain(N)
    reference: float | None = None
    for label, gamma in gamma_patterns(N):
        norm_sq = m_norm_squared(N, H, gamma)
        if reference is None:
            reference = norm_sq
        diff = abs(norm_sq - reference)
        print(f"  {label:8s} γ={['%.3f' % g for g in gamma]}: "
              f"‖M‖² = {norm_sq:.6f}  |Δ from uniform| = {diff:.3e}")
        # Bit-exact equality (operator-level identity, not floating-point coincidence).
        assert diff == 0.0, f"γ-pattern '{label}' deviated by {diff:.3e}; expected exact 0"


def section_2_dissipator_per_string_vanishing() -> None:
    """M_D = 0 per Pauli string: at H=0 the F1 residual is machine-zero for any γ."""
    print("\n(2) Dissipator-block per-Pauli-string vanishing: with H=0,")
    print("    ‖M‖² ≈ 0 to machine precision for non-uniform γ at N=3, 4, 5.")
    print("    (Asserts the operator-level identity M_D = 0, not just averaged-to-zero.)")
    for N in (3, 4, 5):
        H_zero = np.zeros((2 ** N, 2 ** N), dtype=complex)
        for label, gamma in gamma_patterns(N):
            norm_sq = m_norm_squared(N, H_zero, gamma)
            print(f"  N={N} {label:8s} γ=[{gamma[0]:.3f}, ..., {gamma[-1]:.3f}]: ‖M‖² = {norm_sq:.3e}")
            # Tolerance allows for 4^N · ε² floating-point accumulation; at N=5 (1024-dim) the
            # observed values are ~10⁻³¹ (squared machine epsilon scale), comfortably below 1e-20.
            assert norm_sq < 1e-20, f"N={N} {label}: ‖M_D‖² = {norm_sq:.3e} exceeds machine-zero"


def section_3_c_h_factor_scaling() -> None:
    """At N=3, 4, 5 with H = XX+YZ chain, ‖M‖² = c_H · F(N, chain) for any γ pattern."""
    print("\n(3) Scaling against c_H · F(N, G): with H = XX+YZ chain (main class),")
    print("    anchor c_H from N=2 then assert c_H · (N−1) · 4^(N−2) at N=3, 4, 5")
    print("    holds bit-exact for uniform and non-uniform γ.")
    # Anchor c_H at N=2 (one bond). F(2, chain) = (N-1)·4^(N-2) = 1, so c_H = ‖M(N=2)‖².
    N0 = 2
    H0 = main_xx_yz_chain(N0)
    c_H = m_norm_squared(N0, H0, [0.05, 0.05])  # uniform γ at anchor is fine
    print(f"  Anchor c_H from N=2 (XX+YZ, uniform γ=0.05) = {c_H}")
    for N in (3, 4, 5):
        H = main_xx_yz_chain(N)
        F = (N - 1) * 4 ** (N - 2)
        predicted = c_H * F
        for label, gamma in gamma_patterns(N):
            observed = m_norm_squared(N, H, gamma)
            diff = abs(observed - predicted)
            print(f"  N={N} {label:8s}: predicted c_H · F = {predicted}, observed = {observed}, "
                  f"|Δ| = {diff:.3e}")
            # Bit-exact: M_D = 0 means the entire F1 residual norm is c_H · F (the H-block
            # closed form), with no γ-dependent correction in any direction.
            assert diff == 0.0, f"N={N} {label}: ‖M‖² deviated from c_H · F by {diff:.3e}"


def section_4_cross_hamiltonian_invariance() -> None:
    """γ-pattern invariance holds across distinct Π²-classes of H."""
    print("\n(4) Cross-Hamiltonian invariance: γ-pattern independence holds across")
    print("    multiple H classes (soft XY+YX, main XX+YZ) at N=3, 4, 5.")
    for N in (3, 4, 5):
        for h_label, H_builder in [
            ("soft  XY+YX",  soft_xy_yx_chain),
            ("main  XX+YZ",  main_xx_yz_chain),
        ]:
            H = H_builder(N)
            reference: float | None = None
            for label, gamma in gamma_patterns(N):
                norm_sq = m_norm_squared(N, H, gamma)
                if reference is None:
                    reference = norm_sq
                diff = abs(norm_sq - reference)
                print(f"  N={N} H={h_label}  {label:8s}: ‖M‖² = {norm_sq:.6f}  "
                      f"|Δ from uniform| = {diff:.3e}")
                assert diff == 0.0, (
                    f"N={N} H={h_label} {label}: γ-pattern leaked into ‖M‖² by {diff:.3e}; "
                    f"expected exact 0 by the M_D = 0 per-Pauli-string identity")


def main() -> None:
    print("F1 non-uniform γ closure verification")
    print("=" * 78)
    section_1_h_block_gamma_independence()
    section_2_dissipator_per_string_vanishing()
    section_3_c_h_factor_scaling()
    section_4_cross_hamiltonian_invariance()
    print("\nAll sections complete. F1 H-block residual ‖M_H‖² = c_H · F(N, G) is")
    print("γ-independent for uniform AND non-uniform γ; the OpenQuestion is closed by a")
    print("negative result (the conjectured Σγ_l² scaling does not occur because M_D = 0).")


if __name__ == "__main__":
    main()
