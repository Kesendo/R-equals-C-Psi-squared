"""Verify the closed form of ‖M‖² for the F1 residual under depolarizing noise.

The closed form derived in `docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md`:

    ‖M(depol)‖²_F = 4^(N−1) · [ (16/9) · Σγ²  +  16 · (Σγ)² ]

with the residual computed in the framework's orthonormal Pauli-string basis
(`palindrome_residual` in `framework/lindblad.py`). Π is the F1 palindrome
operator from `framework/symmetry.py` (per-site rule I ↔ X with phase +1,
Y ↔ Z with phase +i).

This script verifies:

(1) F1 sanity:    pure Z-dephasing → ‖M‖² ≈ 0.
(2) Pure depol:   fit (a, b) numerically across uniform and non-uniform γ
                  up to N=5; bit-exact match against the (16/9, 16) closed form.
(3) Orthogonality of the depol block to coherent (truly H) blocks.
(4) Orthogonality of the depol block to combined H+Z blocks.
(5) Orthogonality of the depol block to soft (Π²-odd) H (XY+YX).
(6) Per-site M_l kernel: D_depol, Π, M_l matrix display + derivation of (16/9, 16).
(7) Π²-trivial split assertions: ‖M − Π·M·Π⁻¹‖_F < 1e-13 and
                                 ‖M − Π²·M·Π²⁻¹‖_F < 1e-13 (M_anti = 0 exactly).

Note on σ-shift convention: depol's F1 residual is the BARE residual
Π·L·Π⁻¹ + L (σ = 0), NOT the Z-dephasing-style Π·L·Π⁻¹ + L + 2Σγ·I.
M_l = diag(−4/3, −4/3, −8/3, −8/3) has two distinct diagonal values; no
constant 2σ·I can equalise them. F5's scalar (2/3)Σγ is a complementary
trace-projection diagnostic of the same broken palindrome.

Framework primitive note: `simulations/framework/lindblad.py` has no
`lindbladian_depolarizing` primitive. We use `lindbladian_general(H, c_ops)`
with 3 jump operators per site (c_{l, P} = √(γ_l/3)·P_l for P ∈ {X, Y, Z}).
"""
from __future__ import annotations

import sys

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from framework.lindblad import (  # noqa: E402
    lindbladian_general,
    lindbladian_z_dephasing,
    palindrome_residual,
)
from framework.pauli import _build_bilinear, site_op  # noqa: E402
from framework.symmetry import build_pi_full  # noqa: E402


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #


def depol_c_ops(N: int, gamma_l: list[float]) -> list[np.ndarray]:
    """Three Lindblad operators per site: c_{l, P} = √(γ_l / 3) · P_l for P ∈ {X, Y, Z}."""
    c_ops: list[np.ndarray] = []
    for l, g in enumerate(gamma_l):
        if g == 0:
            continue
        s = np.sqrt(g / 3.0)
        for letter in ("X", "Y", "Z"):
            c_ops.append(s * site_op(N, l, letter))
    return c_ops


def m_norm_squared_depol_only(N: int, gamma_l: list[float]) -> float:
    """‖M(depol)‖²_F: pure depolarizing channel, no H, σ-shift = 0."""
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    L = lindbladian_general(H, depol_c_ops(N, gamma_l))
    M = palindrome_residual(L, 0.0, N)  # σ = 0 for depol (NOT Σγ)
    return float(np.real(np.trace(M.conj().T @ M)))


def m_norm_squared_combined(N: int, H: np.ndarray,
                            gamma_z_l: list[float],
                            gamma_depol_l: list[float],
                            sigma: float) -> float:
    """‖M‖²_F for combined H + Z-dephasing + depolarizing at the given σ-shift.

    palindrome_residual is linear in L: M(L_z + L_depol; σ) = M(L_z; σ) + M(L_depol; 0).
    With σ = Σγ_Z the Z-block residual vanishes and the depol-block carries the entire
    Frobenius norm equal to the σ=0 depol prediction.
    """
    c_ops: list[np.ndarray] = []
    for l, g in enumerate(gamma_z_l):
        if g == 0:
            continue
        c_ops.append(np.sqrt(g) * site_op(N, l, "Z"))
    c_ops.extend(depol_c_ops(N, gamma_depol_l))
    L = lindbladian_general(H, c_ops)
    M = palindrome_residual(L, sigma, N)
    return float(np.real(np.trace(M.conj().T @ M)))


def predict_depol(N: int, gamma_l: list[float]) -> float:
    """Closed form ‖M(depol)‖²_F = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]."""
    sum_sq = sum(g * g for g in gamma_l)
    sum_g_sq = sum(gamma_l) ** 2
    return float(4 ** (N - 1)) * (16.0 / 9.0 * sum_sq + 16.0 * sum_g_sq)


def heisenberg_chain(N: int) -> np.ndarray:
    bonds = [(b, b + 1) for b in range(N - 1)]
    return _build_bilinear(N, bonds,
                           [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)])


def xy_yx_soft_chain(N: int) -> np.ndarray:
    """Soft (Π²-odd) Hamiltonian XY+YX: non-zero M; tests H ⊥ depol in soft regime."""
    bonds = [(b, b + 1) for b in range(N - 1)]
    return _build_bilinear(N, bonds, [("X", "Y", 1.0), ("Y", "X", 1.0)])


def m_matrix_depol_only(N: int, gamma_l: list[float]) -> np.ndarray:
    """The full residual M (not its squared norm); used by section 7 for Π-equivariance."""
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    L = lindbladian_general(H, depol_c_ops(N, gamma_l))
    return palindrome_residual(L, 0.0, N)


# --------------------------------------------------------------------------- #
# Sections                                                                    #
# --------------------------------------------------------------------------- #


def section_1_sanity() -> None:
    print("(1) F1 sanity: pure Z-dephasing must give ‖M‖² ≈ 0.")
    for N in (2, 3, 4):
        H = np.zeros((2 ** N, 2 ** N), dtype=complex)
        L = lindbladian_z_dephasing(H, [0.1] * N)
        M = palindrome_residual(L, sum([0.1] * N), N)
        norm_sq = float(np.real(np.trace(M.conj().T @ M)))
        heis_L = lindbladian_z_dephasing(heisenberg_chain(N), [0.1] * N)
        heis_M = palindrome_residual(heis_L, sum([0.1] * N), N)
        heis_norm_sq = float(np.real(np.trace(heis_M.conj().T @ heis_M)))
        print(f"  N={N}: ‖M‖²(Z only) = {norm_sq:.3e}    "
              f"‖M‖²(Heisenberg+Z) = {heis_norm_sq:.3e}    (both expected ≈ 0)")


def section_2_pure_depol_coefficients() -> None:
    print("\n(2) Pure depol only: fit (a, b) in ‖M‖² = 4^(N−1) · [a·Σγ² + b·(Σγ)²]; expect (16/9, 16).")
    for N in (2, 3, 4, 5):
        gamma_u = [0.1] * N
        gamma_nu = [0.05 * (l + 1) for l in range(N)]
        obs_u = m_norm_squared_depol_only(N, gamma_u)
        obs_nu = m_norm_squared_depol_only(N, gamma_nu)
        # Solve 2x2 system: [[Σγ², (Σγ)²], …] · (a, b) = (obs / 4^(N-1), …)
        A = np.array([
            [sum(g * g for g in gamma_u), sum(gamma_u) ** 2],
            [sum(g * g for g in gamma_nu), sum(gamma_nu) ** 2],
        ])
        b_vec = np.array([obs_u, obs_nu]) / (4 ** (N - 1))
        (a_fit, b_fit) = np.linalg.solve(A, b_vec)
        pred_u = predict_depol(N, gamma_u)
        pred_nu = predict_depol(N, gamma_nu)
        print(f"  N={N}:")
        print(f"    uniform γ=0.1:       obs={obs_u:.6f}  pred={pred_u:.6f}  "
              f"|Δ|={abs(obs_u - pred_u):.3e}")
        print(f"    non-uniform γ:       obs={obs_nu:.6f}  pred={pred_nu:.6f}  "
              f"|Δ|={abs(obs_nu - pred_nu):.3e}")
        print(f"    fitted (a, b) = ({a_fit:.6f}, {b_fit:.6f})  "
              f"(expected exactly (16/9 ≈ 1.777778, 16))")
        assert abs(obs_u - pred_u) < 1e-10, f"N={N} uniform mismatch: |Δ| = {abs(obs_u - pred_u):.3e}"
        assert abs(obs_nu - pred_nu) < 1e-10, f"N={N} non-uniform mismatch: |Δ| = {abs(obs_nu - pred_nu):.3e}"


def section_3_orthogonality_truly_h() -> None:
    print("\n(3) Orthogonality H ⊥ depol for truly H (Heisenberg): "
          "‖M(H+depol)‖² should equal ‖M(H)‖² + ‖M(depol)‖² with H-part = 0 (both at σ = 0).")
    for N in (3, 4):
        H = heisenberg_chain(N)
        gamma_depol = [0.1] * N
        # σ = 0 for both: probes the cross-term cleanly.
        H_only_L = lindbladian_general(H, [])
        H_only_M = palindrome_residual(H_only_L, 0.0, N)
        only_h = float(np.real(np.trace(H_only_M.conj().T @ H_only_M)))
        only_depol = m_norm_squared_depol_only(N, gamma_depol)
        together_L = lindbladian_general(H, depol_c_ops(N, gamma_depol))
        together_M = palindrome_residual(together_L, 0.0, N)
        together = float(np.real(np.trace(together_M.conj().T @ together_M)))
        cross = together - only_h - only_depol
        print(f"  N={N}: ‖M‖²(H only, σ=0) = {only_h:.3e}    "
              f"‖M‖²(depol only) = {only_depol:.6f}    "
              f"‖M‖²(H + depol) = {together:.6f}    "
              f"cross = {cross:+.3e}")
        # Heisenberg with σ = 0 has residual = 2·Σγ_Z·I = 0 since γ_Z = 0 here, so M(H only, σ=0) = 0.
        assert abs(cross) < 1e-10, f"N={N} cross-term too large: {cross:.3e}"


def section_4_orthogonality_with_z() -> None:
    print("\n(4) Orthogonality with Z-dephasing: ‖M(H+Z+depol; σ=Σγ_Z)‖² ≡ ‖M(depol)‖²_F.")
    print("    palindrome_residual is linear in L: at σ=Σγ_Z the H+Z block residual = 0 and")
    print("    the depol block residual stays at its bare σ=0 value. Cross-term vanishes by")
    print("    Frobenius-orthogonality of the H+Z and depol blocks (same Step 6 mechanism as F1T1).")
    for N in (3, 4):
        H = heisenberg_chain(N)
        gamma_z = [0.1] * N
        gamma_depol = [0.1] * N
        # H + Z only at σ = Σγ_Z: residual ≈ 0 (truly H, Z-dephasing satisfies F1 bit-exactly).
        hz_L = lindbladian_z_dephasing(H, gamma_z)
        hz_M = palindrome_residual(hz_L, sum(gamma_z), N)
        h_and_z = float(np.real(np.trace(hz_M.conj().T @ hz_M)))
        # depol only at σ = 0: the closed-form residual.
        only_depol = m_norm_squared_depol_only(N, gamma_depol)
        # full H+Z+depol at σ=Σγ_Z: by linearity equals M(H+Z; Σγ_Z) + M(depol; 0).
        full = m_norm_squared_combined(N, H, gamma_z, gamma_depol, sigma=sum(gamma_z))
        cross = full - h_and_z - only_depol
        print(f"  N={N}: ‖M‖²(H+Z; σ=Σγ_Z) = {h_and_z:.3e}    "
              f"‖M‖²(depol; σ=0) = {only_depol:.6f}    "
              f"‖M‖²(H+Z+depol; σ=Σγ_Z) = {full:.6f}    "
              f"cross = {cross:+.3e}")
        assert abs(cross) < 1e-10, f"N={N} cross-term too large: {cross:.3e}"


def section_5_orthogonality_soft_h() -> None:
    print("\n(5) Orthogonality with soft (Π²-odd) H (XY+YX): H-part non-zero; depol stays (16/9, 16).")
    for N in (3, 4):
        H = xy_yx_soft_chain(N)
        gamma_depol = [0.1] * N
        # σ = 0 for both, to probe cross-term cleanly.
        H_only_L = lindbladian_general(H, [])
        H_only_M = palindrome_residual(H_only_L, 0.0, N)
        only_h = float(np.real(np.trace(H_only_M.conj().T @ H_only_M)))
        only_depol = m_norm_squared_depol_only(N, gamma_depol)
        together_L = lindbladian_general(H, depol_c_ops(N, gamma_depol))
        together_M = palindrome_residual(together_L, 0.0, N)
        together = float(np.real(np.trace(together_M.conj().T @ together_M)))
        cross = together - only_h - only_depol
        print(f"  N={N}: ‖M‖²(soft H) = {only_h:.6f}    "
              f"‖M‖²(depol only) = {only_depol:.6f}    "
              f"‖M‖²(soft H + depol) = {together:.6f}    "
              f"cross = {cross:+.3e}")
        assert abs(cross) < 1e-9, f"N={N} cross-term too large: {cross:.3e}"


def section_6_per_site_kernel() -> None:
    """Display the per-site Pauli-basis M_l matrix that underlies the (16/9, 16) closed form."""
    print("\n(6) Per-site M_depol,l Pauli-basis matrix at γ=1, with derivation of (16/9, 16).")
    paulis = {
        "I": np.eye(2, dtype=complex),
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.array([[1, 0], [0, -1]], dtype=complex),
    }
    labels = ("I", "X", "Y", "Z")

    # D_depol(ρ) = (1/3) · Σ_{P ∈ {X,Y,Z}} (P ρ P − ρ) at γ = 1.
    D = np.zeros((4, 4), dtype=complex)
    for j, lj in enumerate(labels):
        rho = paulis[lj]
        out = np.zeros((2, 2), dtype=complex)
        for P_letter in ("X", "Y", "Z"):
            P = paulis[P_letter]
            out += (1.0 / 3.0) * (P @ rho @ P - rho)
        for i, li in enumerate(labels):
            D[i, j] = 0.5 * np.trace(paulis[li] @ out)

    # Per-site Π in (I, X, Y, Z) basis: I↔X phase 1, Y↔Z phase i.
    Pi = np.array([
        [0, 1, 0, 0],   # I row receives X
        [1, 0, 0, 0],   # X row receives I
        [0, 0, 0, 1j],  # Y row receives iZ
        [0, 0, 1j, 0],  # Z row receives iY
    ], dtype=complex)
    Pi_inv = Pi.conj().T

    # Cross-check: the hand-built per-site Π must equal the framework's
    # per-site Π (build_pi_full at N=1 in the Z-dephasing convention).
    Pi_framework = build_pi_full(1, dephase_letter="Z")
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
            return f"{z.real:+.3f}"
        if abs(z.real) < 1e-12:
            return f"{z.imag:+.3f}i"
        return f"{z.real:+.3f}{z.imag:+.3f}i"

    print("  D_depol (γ=1) in orthonormal Pauli basis (rows = output, cols = input):")
    for i, li in enumerate(labels):
        row = " ".join(_fmt(D[i, j]).rjust(10) for j in range(4))
        print(f"    {li} | {row}")

    print("\n  Π in single-letter Pauli basis (signed permutation, order-4 unitary):")
    for i, li in enumerate(labels):
        row = " ".join(_fmt(Pi[i, j]).rjust(10) for j in range(4))
        print(f"    {li} | {row}")

    print("\n  M_l := Π·D·Π⁻¹ + D (per-site, γ=1):")
    for i, li in enumerate(labels):
        row = " ".join(_fmt(M_per_site[i, j]).rjust(10) for j in range(4))
        print(f"    {li} | {row}")

    print(f"\n  ‖M_l‖²_F = {norm_sq_per:.6f}  (= 160/9 ≈ {160/9:.6f}; 2·(4/3)² + 2·(8/3)² = 2·(16/9 + 64/9))")
    print(f"  tr(M_l)  = {tr_per.real:+.6f} {tr_per.imag:+.6f}i  "
          f"(= −8; −4/3 − 4/3 − 8/3 − 8/3; |tr|² = 64 drives the cross-site coefficient)")
    # Bit-exact assertions on the per-site kernel.
    assert abs(norm_sq_per - 160 / 9) < 1e-12, f"per-site ‖M_l‖² mismatch: {norm_sq_per} vs {160/9}"
    assert abs(tr_per.real + 8) < 1e-12 and abs(tr_per.imag) < 1e-12, f"per-site tr mismatch: {tr_per}"

    print("\n  Multi-site assembly (per-site action, identity at other sites):")
    print("    tr(M_l† M_l)   = 4^(N−1) · ‖M_per_site‖² = (160/9) · 4^(N−1)")
    print("    tr(M_l† M_l')  = |tr(M_per_site)|² · 4^(N−2) = 64 · 4^(N−2) = 16 · 4^(N−1)  (l ≠ l')")
    print("    ‖Σ γ_l M_l‖²   = (160/9)·4^(N−1)·Σγ²  +  16·4^(N−1)·[(Σγ)² − Σγ²]")
    print("                    = 4^(N−1) · [(160/9 − 16)·Σγ² + 16·(Σγ)²]")
    print("                    = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]  ✓")


def section_7_pi2_trivial_split() -> None:
    """Assert Π and Π² are exact symmetries of M(depol): M_anti = 0 identically."""
    print("\n(7) Π²-trivial split: M(depol) is Pauli-basis-diagonal ⟹ Π·M·Π⁻¹ = M exactly.")
    print("    Stronger than Π²·M·Π²⁻¹ = M: full Π conjugation is identity. Hence M_anti = 0,")
    print("    distinguishing depol from T1 (where M_anti = D_{T1, odd} carries F82/F84 content).")
    for N in (2, 3, 4):
        gamma = [0.1] * N
        M = m_matrix_depol_only(N, gamma)
        Pi = build_pi_full(N, dephase_letter="Z")
        Pi_inv = Pi.conj().T
        Pi2 = Pi @ Pi
        Pi2_inv = Pi_inv @ Pi_inv
        diff_pi = float(np.linalg.norm(M - Pi @ M @ Pi_inv))
        diff_pi2 = float(np.linalg.norm(M - Pi2 @ M @ Pi2_inv))
        # M_anti = (M − Π·M·Π⁻¹) / 2, expected zero.
        M_anti = (M - Pi @ M @ Pi_inv) / 2.0
        anti_norm_sq = float(np.real(np.trace(M_anti.conj().T @ M_anti)))
        print(f"  N={N}: ‖M − Π·M·Π⁻¹‖_F  = {diff_pi:.3e}    "
              f"‖M − Π²·M·Π²⁻¹‖_F = {diff_pi2:.3e}    "
              f"‖M_anti‖²_F = {anti_norm_sq:.3e}")
        assert diff_pi < 1e-13, f"N={N} Π-equivariance violated: ‖M − Π·M·Π⁻¹‖_F = {diff_pi:.3e}"
        assert diff_pi2 < 1e-13, f"N={N} Π²-equivariance violated: ‖M − Π²·M·Π²⁻¹‖_F = {diff_pi2:.3e}"
        assert anti_norm_sq < 1e-26, f"N={N} M_anti should be zero: ‖M_anti‖² = {anti_norm_sq:.3e}"


def main() -> None:
    print("F1 depol-residual closed-form verification")
    print("=" * 78)
    section_1_sanity()
    section_2_pure_depol_coefficients()
    section_3_orthogonality_truly_h()
    section_4_orthogonality_with_z()
    section_5_orthogonality_soft_h()
    section_6_per_site_kernel()
    section_7_pi2_trivial_split()
    print("\nAll sections complete. "
          "‖M(depol)‖² = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]  (in framework Pauli basis).")
    print("Structural surprises confirmed: M_anti(depol) = 0 (Π²-trivial); σ-shift = 0 (not Σγ).")


if __name__ == "__main__":
    main()
