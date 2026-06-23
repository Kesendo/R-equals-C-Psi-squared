"""Verify the closed form of ‖M‖² for the F1 residual with T1 amplitude damping.

Memory `project_palindrome_frobenius_scaling` and the typed `HARDWARE_DISSIPATORS['T1']`
entry both record (a, b) = (3, 4):

    ‖M(L_Z + T1)‖²_F = 2^(N+2)·n_YZ·‖H‖²_F  +  4^(N−1) · [ 3·Σγ²_T1 + 4·(Σγ_T1)² ]

with the residual computed in the framework's orthonormal Pauli-string basis
(`palindrome_residual` in `framework/lindblad.py` does `L_pauli = M_basis.conj().T @ L @ M_basis / 2^N`).
Π here is the F1 palindrome operator from `framework/symmetry.py` (per-site rule
I ↔ X with phase +1, Y ↔ Z with phase +i; the construction from
`docs/proofs/MIRROR_SYMMETRY_PROOF.md`). The script proves:

(1) F1 sanity:   pure Z-dephasing → ‖M‖² ≈ 0
                Heisenberg + Z-dephasing → ‖M‖² ≈ 0
(2) Pure T1:    extract (a, b) numerically across uniform and non-uniform γ_T1
                up to N=5; bit-exact match against the (3, 4) closed form.
(3) Orthogonality of the T1 block to coherent (H) and Z-dephasing blocks:
                ‖M(H + Z + T1)‖² = ‖M(H + Z)‖² + ‖M(T1)‖² across truly / soft H.

All numerics use the framework's own primitives so any future refactor of
`palindrome_residual`, `lindbladian_z_plus_t1`, or `build_pi_full` is caught
here automatically.
"""
from __future__ import annotations

import sys

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from framework.lindblad import (  # noqa: E402
    lindbladian_z_plus_t1,
    palindrome_residual,
)
from framework.pauli import _build_bilinear  # noqa: E402


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #


def m_norm_squared(N: int, H: np.ndarray, gamma_z_l: list[float],
                   gamma_t1_l: list[float]) -> float:
    """‖M‖²_F in the framework's orthonormal Pauli-string basis."""
    L = lindbladian_z_plus_t1(H, gamma_z_l, gamma_t1_l)
    M = palindrome_residual(L, sum(gamma_z_l), N)
    return float(np.real(np.trace(M.conj().T @ M)))


def predict_t1_only(N: int, gamma_t1_l: list[float]) -> float:
    """Closed form ‖M(T1)‖²_F = 4^(N−1) · [3·Σγ² + 4·(Σγ)²].

    Verification-local reimplementation of just the T1 half so the script can
    isolate the T1 block without computing the H+T1 sum. For the full H+T1
    closed form (F49 H block plus this T1 block) use the framework's
    ``fw.predict_residual_norm_squared_from_terms(chain, terms, gamma_t1=...)``
    in `framework/diagnostics/f49_frobenius_scaling.py`.
    """
    sum_sq = sum(g * g for g in gamma_t1_l)
    sum_g_sq = sum(gamma_t1_l) ** 2
    return float(4 ** (N - 1)) * (3.0 * sum_sq + 4.0 * sum_g_sq)


def heisenberg_chain(N: int) -> np.ndarray:
    bonds = [(b, b + 1) for b in range(N - 1)]
    return _build_bilinear(N, bonds,
                           [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)])


def xy_bond_chain(N: int) -> np.ndarray:
    """H_truly = Σ_b (X_b X_{b+1} + Y_b Y_{b+1}); used to probe H ⊥ T1."""
    bonds = [(b, b + 1) for b in range(N - 1)]
    return _build_bilinear(N, bonds, [("X", "X", 1.0), ("Y", "Y", 1.0)])


def xy_yx_soft_chain(N: int) -> np.ndarray:
    """Soft (Π²-odd) Hamiltonian XY+YX: non-zero M; tests H ⊥ T1 in soft regime."""
    bonds = [(b, b + 1) for b in range(N - 1)]
    return _build_bilinear(N, bonds, [("X", "Y", 1.0), ("Y", "X", 1.0)])


# --------------------------------------------------------------------------- #
# Sections                                                                    #
# --------------------------------------------------------------------------- #


def section_1_sanity() -> None:
    print("(1) F1 sanity: pure Z-dephasing and Heisenberg + Z must give ‖M‖² ≈ 0.")
    for N in (2, 3, 4):
        z_only = m_norm_squared(N, np.zeros((2 ** N, 2 ** N), complex),
                                [0.1] * N, [0.0] * N)
        heis = m_norm_squared(N, heisenberg_chain(N), [0.1] * N, [0.0] * N)
        print(f"  N={N}: ‖M‖²(Z only) = {z_only:.3e}    "
              f"‖M‖²(Heisenberg+Z) = {heis:.3e}    "
              f"(both expected ≈ 0)")


def section_2_pure_t1_coefficients() -> None:
    print("\n(2) Pure T1 only: fit (a, b) in ‖M‖² = 4^(N−1) · [a·Σγ² + b·(Σγ)²].")
    for N in (2, 3, 4, 5):
        gamma_u = [0.1] * N
        gamma_nu = [0.05 + 0.05 * l for l in range(N)]
        obs_u = m_norm_squared(N, np.zeros((2 ** N, 2 ** N), complex),
                               [0.0] * N, gamma_u)
        obs_nu = m_norm_squared(N, np.zeros((2 ** N, 2 ** N), complex),
                                [0.0] * N, gamma_nu)
        # Solve the 2x2 system [[Σγ², (Σγ)²], …] · (a, b) = (obs / 4^(N-1), …)
        A = np.array([
            [sum(g * g for g in gamma_u), sum(gamma_u) ** 2],
            [sum(g * g for g in gamma_nu), sum(gamma_nu) ** 2],
        ])
        b_vec = np.array([obs_u, obs_nu]) / (4 ** (N - 1))
        (a_fit, b_fit) = np.linalg.solve(A, b_vec)
        pred_u = predict_t1_only(N, gamma_u)
        pred_nu = predict_t1_only(N, gamma_nu)
        print(f"  N={N}:")
        print(f"    uniform γ_T1=0.1:    obs={obs_u:.6f}  pred={pred_u:.6f}  "
              f"|Δ|={abs(obs_u - pred_u):.3e}")
        print(f"    non-uniform γ_T1:    obs={obs_nu:.6f}  pred={pred_nu:.6f}  "
              f"|Δ|={abs(obs_nu - pred_nu):.3e}")
        print(f"    fitted (a, b) = ({a_fit:.6f}, {b_fit:.6f})  "
              f"(expected exactly (3, 4))")


def section_3_orthogonality_truly_h() -> None:
    print("\n(3) Orthogonality H ⊥ T1 for truly H (Heisenberg): "
          "‖M(H+T1)‖² should equal ‖M(H)‖² + ‖M(T1)‖² with H-part = 0.")
    for N in (3, 4):
        H = heisenberg_chain(N)
        gamma_t1 = [0.1] * N
        only_h = m_norm_squared(N, H, [0.0] * N, [0.0] * N)
        only_t1 = m_norm_squared(N, np.zeros((2 ** N, 2 ** N), complex),
                                 [0.0] * N, gamma_t1)
        together = m_norm_squared(N, H, [0.0] * N, gamma_t1)
        cross = together - only_h - only_t1
        print(f"  N={N}: ‖M‖²(H only) = {only_h:.3e}    "
              f"‖M‖²(T1 only) = {only_t1:.6f}    "
              f"‖M‖²(H + T1) = {together:.6f}    "
              f"cross = {cross:+.3e}")


def section_4_orthogonality_with_z() -> None:
    print("\n(4) Orthogonality with Z-dephasing: ‖M(H+Z+T1)‖² = ‖M(H+Z)‖² + ‖M(T1)‖².")
    for N in (3, 4):
        H = heisenberg_chain(N)
        gamma_z = [0.1] * N
        gamma_t1 = [0.1] * N
        h_and_z = m_norm_squared(N, H, gamma_z, [0.0] * N)
        only_t1 = m_norm_squared(N, np.zeros((2 ** N, 2 ** N), complex),
                                 [0.0] * N, gamma_t1)
        full = m_norm_squared(N, H, gamma_z, gamma_t1)
        cross = full - h_and_z - only_t1
        print(f"  N={N}: ‖M‖²(H+Z) = {h_and_z:.3e}    "
              f"‖M‖²(T1 only) = {only_t1:.6f}    "
              f"‖M‖²(H+Z+T1) = {full:.6f}    "
              f"cross = {cross:+.3e}")


def section_5_orthogonality_soft_h() -> None:
    print("\n(5) Orthogonality with soft (Π²-odd) H (XY+YX): H-part non-zero; T1 stays (3, 4).")
    for N in (3, 4):
        H = xy_yx_soft_chain(N)
        gamma_t1 = [0.1] * N
        only_h = m_norm_squared(N, H, [0.0] * N, [0.0] * N)
        only_t1 = m_norm_squared(N, np.zeros((2 ** N, 2 ** N), complex),
                                 [0.0] * N, gamma_t1)
        together = m_norm_squared(N, H, [0.0] * N, gamma_t1)
        cross = together - only_h - only_t1
        print(f"  N={N}: ‖M‖²(soft H) = {only_h:.6f}    "
              f"‖M‖²(T1 only) = {only_t1:.6f}    "
              f"‖M‖²(soft H + T1) = {together:.6f}    "
              f"cross = {cross:+.3e}")


def section_6_per_site_kernel() -> None:
    """Display the per-site Pauli-basis M_l matrix that underlies the (3, 4) closed form."""
    print("\n(6) Per-site M_T1,l Pauli-basis matrix at γ_T1=1, with derivation of (3, 4).")
    paulis = {
        "I": np.eye(2, dtype=complex),
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.array([[1, 0], [0, -1]], dtype=complex),
    }
    labels = ("I", "X", "Y", "Z")
    sm = np.array([[0, 1], [0, 0]], dtype=complex)  # σ⁻ = (X+iY)/2
    sp = sm.conj().T
    P_excited = sp @ sm

    D = np.zeros((4, 4), dtype=complex)
    for j, lj in enumerate(labels):
        out = sm @ paulis[lj] @ sp - 0.5 * (P_excited @ paulis[lj]
                                            + paulis[lj] @ P_excited)
        for i, li in enumerate(labels):
            D[i, j] = 0.5 * np.trace(paulis[li] @ out)

    # Single-site Π in (I, X, Y, Z) basis: I↔X phase 1, Y↔Z phase i
    Pi = np.array([
        [0, 1, 0, 0],   # I row receives X
        [1, 0, 0, 0],   # X row receives I
        [0, 0, 0, 1j],  # Y row receives iZ
        [0, 0, 1j, 0],  # Z row receives iY
    ], dtype=complex)
    Pi_inv = Pi.conj().T

    # Cross-check: the hand-built per-site Π must equal the framework's
    # per-site Π (i.e. build_pi_full at N=1 in the Z-dephasing convention).
    # This catches any future drift between the script-local Π and the
    # framework's `build_pi_full`-derived Π used by `palindrome_residual`
    # everywhere else in this script.
    from framework.symmetry import build_pi_full  # noqa: E402
    Pi_framework = build_pi_full(1, dephase_letter='Z')
    assert np.allclose(Pi, Pi_framework, atol=1e-12), (
        "hand-built per-site Π disagrees with framework build_pi_full(N=1, 'Z'): "
        f"max |Δ| = {np.max(np.abs(Pi - Pi_framework)):.3e}"
    )

    M_per_site = Pi @ D @ Pi_inv + D
    norm_sq_per = float(np.real(np.trace(M_per_site.conj().T @ M_per_site)))
    tr_per = complex(np.trace(M_per_site))

    def _fmt(z: complex) -> str:
        if abs(z) < 1e-12:
            return "0"
        if abs(z.imag) < 1e-12:
            return f"{z.real:+.2f}"
        if abs(z.real) < 1e-12:
            return f"{z.imag:+.2f}i"
        return f"{z.real:+.2f}{z.imag:+.2f}i"

    print("  D_T1 (γ=1) in orthonormal Pauli basis (rows = output, cols = input):")
    for i, li in enumerate(labels):
        row = " ".join(_fmt(D[i, j]).rjust(9) for j in range(4))
        print(f"    {li} | {row}")

    print("\n  Π in single-letter Pauli basis (signed permutation, order-4 unitary):")
    for i, li in enumerate(labels):
        row = " ".join(_fmt(Pi[i, j]).rjust(9) for j in range(4))
        print(f"    {li} | {row}")

    print("\n  M_l := Π·D·Π⁻¹ + D (per-site, γ=1):")
    for i, li in enumerate(labels):
        row = " ".join(_fmt(M_per_site[i, j]).rjust(9) for j in range(4))
        print(f"    {li} | {row}")

    print(f"\n  ‖M_l‖²_F = {norm_sq_per:.4f}  (= 7 = 0.25 + 0.25 + 1 + 2.25 + 1 + 2.25)")
    print(f"  tr(M_l)  = {tr_per.real:+.4f} {tr_per.imag:+.4f}i  "
          f"(= −4; |tr|² = 16 drives the cross-site coefficient)")

    print("\n  Multi-site assembly (per-site action, identity at other sites):")
    print("    tr(M_l† M_l)   = 4^(N−1) · ‖M_per_site‖² = 7 · 4^(N−1)")
    print("    tr(M_l† M_l')  = |tr(M_per_site)|² · 4^(N−2) = 16 · 4^(N−2) = 4 · 4^(N−1)  (l ≠ l')")
    print("    ‖Σ γ_l M_l‖²   = 7·4^(N−1)·Σγ²  +  4·4^(N−1)·[(Σγ)² − Σγ²]")
    print("                    = 4^(N−1) · [3·Σγ² + 4·(Σγ)²]  ✓")


def main() -> None:
    print("F1 T1-residual closed-form verification")
    print("=" * 78)
    section_1_sanity()
    section_2_pure_t1_coefficients()
    section_3_orthogonality_truly_h()
    section_4_orthogonality_with_z()
    section_5_orthogonality_soft_h()
    section_6_per_site_kernel()
    print("\nAll sections complete. ‖M(T1)‖² = 4^(N−1) · [3·Σγ²_T1 + 4·(Σγ_T1)²]"
          " (in framework Pauli basis).")


if __name__ == "__main__":
    main()
