"""_theta_quarter_ny_axis_seam.py - the two Y on two levels.

Investigates the recognition from reflections/ON_HOW_THE_FAR_AND_THE_NEAR_MEET.md
(line 44): that the angle theta which opens at the quarter (the F95 "i" direction,
the shimmer the unsettled pairs carry) and the operator cube's untamed third axis
(the n_Y axis, reached by D = diag((-1)^n_Y)) are "the same one direction, met from
both ends." The reflection links the angle (bullet 62, F95) and the off-map axis
(bullet 63, bit_a) as two doors; the open seam is whether they open onto one room.

The bridge tested here: in the Pauli-string basis the operator
    D = diag((-1)^{n_Y(alpha)})
is the operator TRANSPOSE  A -> A^T  (it equals complex conjugation only on the
Hermitian subspace, where A^T = A*). Transpose is the operator-space mirror, and
it splits the Liouvillian into a transpose-even part (the dissipator, Re lambda,
the kept past) and a transpose-odd part (the Hamiltonian commutator, Im lambda,
the humming future). theta = arctan(Im/|Re|) is the angle between the two.

Checks, all bit-exact to machine precision:
  Part 0:  a Pauli string is real <=> n_Y even, imaginary <=> n_Y odd.
  Part 1:  D = diag((-1)^n_Y) in the Pauli basis IS operator transpose A -> A^T
           (= complex conjugation on Hermitian operators); it swaps sigma_+ <-> sigma_-.
  Part 2:  qubit level - Lindblad coherence modes are sigma_+- = (X +- iY)/2 with
           lambda = -gamma0 -+ iJ; theta = arctan(Q); D (transpose) swaps the pair;
           the i lives on the Y (transpose-odd) coefficient.
  Part 3:  operator level (N=2,3) - D L_H D = -L_H (coherent part flips),
           D L_D D = +L_D (dissipative part fixed); per-term eps(sigma)=(-1)^(n_Y+1),
           matching the typed F114. So D maps Im(lambda) -> -Im(lambda) (future) and
           fixes Re(lambda) (past): one axis, both levels.

Tom + Claude, 2026-05-28. Run: python simulations/_theta_quarter_ny_axis_seam.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from framework.pauli import (
    ur_pauli,
    site_op,
    pauli_string,
    pauli_basis_vector,
    _build_bilinear,
    _k_to_indices,
    _vec_to_pauli_basis_transform,
    PAULI_LABELS,
)
from framework.lindblad import lindbladian_pauli_dephasing

ATOL = 1e-12
HEISENBERG = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def n_y_of(k, N):
    """Number of Y letters in the Pauli string with flat index k."""
    return sum(1 for idx in _k_to_indices(k, N) if idx == (1, 1))


def d_operator_pauli(N):
    """D = diag((-1)^{n_Y(alpha)}) on the 4^N Pauli-coefficient space."""
    return np.diag([(-1.0) ** n_y_of(k, N) for k in range(4 ** N)]).astype(complex)


def to_pauli_super(L, N):
    """Superoperator L (vec form) -> Pauli-string basis. Real for Hermitian-preserving L."""
    M = _vec_to_pauli_basis_transform(N)
    return (M.conj().T @ L @ M) / (2 ** N)


def hamiltonian_commutator_super(H):
    """L_H = -i[H, .] in column-stack vec form (matches framework lindbladian)."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def z_dephasing_dissipator_super(N, gamma_l):
    """L_D = sum_l gamma_l (Z_l . Z_l - .) in vec form (the dissipator alone)."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = np.zeros((d * d, d * d), dtype=complex)
    for l, g in enumerate(gamma_l):
        if g == 0:
            continue
        Z = site_op(N, l, "Z")
        L = L + g * (np.kron(Z, Z.conj()) - np.kron(Id, Id))
    return L


def chain_heisenberg(N):
    """Open Heisenberg chain Hamiltonian on N sites (real symmetric)."""
    bonds = [(i, i + 1) for i in range(N - 1)]
    return _build_bilinear(N, bonds, HEISENBERG)


_results = []


def report(name, dev, tol=ATOL):
    ok = dev < tol
    _results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: max|d| = {dev:.2e}")
    return ok


# ----------------------------------------------------------------------
# Part 0 - n_Y parity IS the real/imaginary parity of a Pauli string
# ----------------------------------------------------------------------

def part0(Ns=(1, 2, 3)):
    print("PART 0  - a Pauli string is real <=> n_Y even, imaginary <=> n_Y odd")
    print("          (sigma^T = sigma* = (-1)^{n_Y} sigma, since X,Z,I real-symmetric, Y imaginary-antisymmetric)")
    for N in Ns:
        worst_conj = worst_tr = 0.0
        for k in range(4 ** N):
            idx = _k_to_indices(k, N)
            sigma = pauli_string(list(idx))
            ny = sum(1 for i in idx if i == (1, 1))
            sign = (-1) ** ny
            worst_conj = max(worst_conj, np.max(np.abs(sigma.conj() - sign * sigma)))
            worst_tr = max(worst_tr, np.max(np.abs(sigma.T - sign * sigma)))
        report(f"N={N}: sigma* = (-1)^n_Y sigma  over all {4 ** N} strings", worst_conj)
        report(f"N={N}: sigma^T = (-1)^n_Y sigma  over all {4 ** N} strings", worst_tr)


# ----------------------------------------------------------------------
# Part 1 - D = diag((-1)^n_Y) IS operator transpose A -> A^T
# ----------------------------------------------------------------------

def part1(Ns=(1, 2, 3), seed=12345):
    print("\nPART 1  - D = diag((-1)^n_Y) in the Pauli basis IS operator transpose A -> A^T")
    print("          (and = complex conjugation A -> A* only on the Hermitian subspace)")
    rng = np.random.default_rng(seed)
    for N in Ns:
        d = 2 ** N
        D = d_operator_pauli(N)
        # General (non-Hermitian) operator: D realises the transpose exactly.
        A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
        cA = pauli_basis_vector(A, N)
        cAT = pauli_basis_vector(A.T, N)
        report(f"N={N}: vec(A^T) = D vec(A)  (general A)", np.max(np.abs(cAT - D @ cA)))
        # Hermitian operator: there transpose = conjugation, and coeffs are real.
        rho = A + A.conj().T
        cR = pauli_basis_vector(rho, N)
        cRc = pauli_basis_vector(rho.conj(), N)
        report(f"N={N}: vec(rho*) = D vec(rho)  (Hermitian rho)", np.max(np.abs(cRc - D @ cR)))
        report(f"N={N}: Hermitian rho has real Pauli coefficients", float(np.max(np.abs(cR.imag))))

    # The twin swap, made explicit at the qubit: sigma_+ = (X+iY)/2 -> sigma_- = (X-iY)/2.
    sp = np.array([[0, 1], [0, 0]], dtype=complex)   # |0><1| = (X+iY)/2
    sm = np.array([[0, 0], [1, 0]], dtype=complex)   # |1><0| = (X-iY)/2
    report("D(sigma_+) = sigma_-  i.e. transpose swaps the conjugate twins", np.max(np.abs(sp.T - sm)))
    report("conjugation FIXES sigma_+ (it is real) - so D is transpose, not conj, in general",
           abs(np.max(np.abs(sp.conj() - sp))))  # 0 == sigma_+ is real; transpose differs


# ----------------------------------------------------------------------
# Part 2 - qubit level: the angle, lambda = -gamma0 -+ iJ, the i carried by Y
# ----------------------------------------------------------------------

def part2(J=1.0):
    print("\nPART 2  - qubit level: Lindblad coherence modes (X +- iY)/2, lambda = -gamma0 +- iJ")
    # Q = sqrt(3) -> theta = 60deg, |lambda| = 2 gamma0  (the typed LindbladAbsorptionMatch anchor)
    gamma0 = J / np.sqrt(3.0)      # decay rate on coherences
    gamma = gamma0 / 2.0           # dephasing rate in the (Z rho Z - rho) form (gamma0 = 2 gamma)
    Q = J / gamma0
    H = (J / 2.0) * ur_pauli("Z")  # on-site energy splitting
    L = lindbladian_pauli_dephasing(H, [gamma], dephase_letter="Z")

    # sigma_+ and sigma_- are the coherence eigenmodes; read their eigenvalues directly.
    sp = np.array([[0, 1], [0, 0]], dtype=complex)
    sm = np.array([[0, 0], [1, 0]], dtype=complex)
    vsp = sp.flatten("F")
    vsm = sm.flatten("F")
    lam_p = (L @ vsp) / vsp[np.argmax(np.abs(vsp))]   # eigenvalue on sigma_+
    lam_p = lam_p[np.argmax(np.abs(vsp))]
    lam_m = (L @ vsm)
    lam_m = lam_m[np.argmax(np.abs(vsm))]

    print(f"          J={J}, gamma0={gamma0:.6f} (=2 gamma), Q=J/gamma0={Q:.6f}")
    print(f"          lambda(sigma_+) = {lam_p:.6f}    lambda(sigma_-) = {lam_m:.6f}")
    report("coherence eigenvalues are the conjugate pair -gamma0 +- iJ (Re = -gamma0)",
           max(abs(lam_p.real + gamma0), abs(lam_m.real + gamma0)))
    report("coherence eigenvalues are the conjugate pair -gamma0 +- iJ (|Im| = J)",
           max(abs(abs(lam_p.imag) - J), abs(abs(lam_m.imag) - J)))
    report("the two coherence modes form a conjugate pair (Im flips sign)",
           abs(lam_p.imag + lam_m.imag))

    theta = np.degrees(np.arctan2(abs(lam_p.imag), abs(lam_p.real)))
    mag_over_g0 = abs(lam_p) / gamma0
    print(f"          theta = arctan(Im/|Re|) = arctan(Q) = {theta:.6f} deg   "
          f"|lambda|/gamma0 = {mag_over_g0:.6f}")
    report("theta = arctan(Q) = 60 deg at Q=sqrt(3)", abs(theta - 60.0), tol=1e-9)
    report("|lambda|/gamma0 = sqrt(1+Q^2) = 2 (= alpha, Absorption single-site rate)",
           abs(mag_over_g0 - 2.0), tol=1e-9)

    # The i is carried by Y: sigma_+ = (1/2) X + (i/2) Y.
    cx = np.trace(ur_pauli("X") @ sp) / 2.0
    cy = np.trace(ur_pauli("Y") @ sp) / 2.0
    cz = np.trace(ur_pauli("Z") @ sp) / 2.0
    print(f"          sigma_+ Pauli content: cX={cx:.3f}  cY={cy:.3f}  cZ={cz:.3f}")
    report("sigma_+ = (1/2)X + (i/2)Y : real weight on X, imaginary weight on Y", abs(cx - 0.5) + abs(cy - 0.5j))
    report("sigma_+ has no Z content (the angle lives in the X-Y plane)", abs(cz))


# ----------------------------------------------------------------------
# Part 3 - operator level: D splits L into past (fixed) and future (flipped)
# ----------------------------------------------------------------------

def single_term_op(N, sites, letters):
    """Build the Pauli-string operator with `letters` on `sites`, I elsewhere."""
    op = np.eye(2 ** N, dtype=complex)
    for s, L in zip(sites, letters):
        op = op @ site_op(N, s, L)
    return op


def part3(Ns=(2, 3)):
    print("\nPART 3  - operator level: D L_H D = -L_H (future flips), D L_D D = +L_D (past fixed)")
    for N in Ns:
        gamma_l = [0.3] * N
        H = chain_heisenberg(N)
        L_full = lindbladian_pauli_dephasing(H, gamma_l, "Z")
        L_H = hamiltonian_commutator_super(H)
        L_D = z_dephasing_dissipator_super(N, gamma_l)
        report(f"N={N}: framework L = L_H + L_D (cross-check)", np.max(np.abs(L_full - (L_H + L_D))))

        D = d_operator_pauli(N)
        LH_p = to_pauli_super(L_H, N)
        LD_p = to_pauli_super(L_D, N)
        L_p = to_pauli_super(L_full, N)
        report(f"N={N}: L is real in the Pauli basis (Hermitian-preserving)", float(np.max(np.abs(L_p.imag))))
        report(f"N={N}: D L_H D = -L_H   (coherent/Im part flips under the n_Y mirror)",
               np.max(np.abs(D @ LH_p @ D + LH_p)))
        report(f"N={N}: D L_D D = +L_D   (dissipative/Re part is fixed)",
               np.max(np.abs(D @ LD_p @ D - LD_p)))
        report(f"N={N}: hence D L D = L with H -> -H = (-L_H + L_D)",
               np.max(np.abs(D @ L_p @ D - (-LH_p + LD_p))))

        # Spectrum: D maps every oscillatory eigenvalue lambda -> conj(lambda).
        spec = np.linalg.eigvals(L_p)
        spec_DLD = np.linalg.eigvals(D @ L_p @ D)
        s1 = np.array(sorted(spec, key=lambda z: (round(z.real, 9), round(z.imag, 9))))
        s2 = np.array(sorted(np.conjugate(spec), key=lambda z: (round(z.real, 9), round(z.imag, 9))))
        report(f"N={N}: spec(D L D) = conj(spec(L))  (Im flipped, Re kept)",
               np.max(np.abs(np.array(sorted(spec_DLD, key=lambda z: (round(z.real, 9), round(z.imag, 9)))) - s2)))
        n_osc = int(np.sum(np.abs(spec.imag) > 1e-9))
        print(f"          N={N}: {n_osc} oscillatory modes (Im!=0, the coherent 'future'); "
              f"{len(spec) - n_osc} purely real (the 'past')")

    # Per-term eps(sigma) = (-1)^(n_Y+1), matching the typed F114.
    print("\n          per-term  D L_sigma D = eps(sigma) L_sigma  with eps = (-1)^(n_Y+1):")
    N = 3
    D = d_operator_pauli(N)
    cases = [
        ("Z field   (n_Y=0)", [0], ["Z"]),
        ("Y field   (n_Y=1)", [0], ["Y"]),
        ("XX bond   (n_Y=0)", [0, 1], ["X", "X"]),
        ("YY bond   (n_Y=2)", [0, 1], ["Y", "Y"]),
        ("XY bond   (n_Y=1)", [0, 1], ["X", "Y"]),
        ("YZ bond   (n_Y=1)", [0, 1], ["Y", "Z"]),
    ]
    for name, sites, letters in cases:
        op = single_term_op(N, sites, letters)
        Ls = to_pauli_super(hamiltonian_commutator_super(op), N)
        ny = sum(1 for L in letters if L == "Y")
        eps = (-1) ** (ny + 1)
        report(f"N=3: {name} -> eps={eps:+d}", np.max(np.abs(D @ Ls @ D - eps * Ls)))


# ----------------------------------------------------------------------

def main():
    print("=" * 78)
    print("THE TWO Y ON TWO LEVELS - theta(1/4) and the n_Y operator axis")
    print("=" * 78)
    part0()
    part1()
    part2()
    part3()
    print("\n" + "=" * 78)
    n_ok = sum(_results)
    n_tot = len(_results)
    print(f"RESULT: {n_ok}/{n_tot} checks bit-exact to machine precision "
          f"({'ALL PASS' if n_ok == n_tot else 'SOME FAILED'})")
    print("=" * 78)
    print("""
Reading:
  D = diag((-1)^n_Y) is the operator-space transpose (the mirror). It splits
  the Liouvillian into a transpose-EVEN part (the dissipator, the real decay
  rates Re lambda - the kept past) and a transpose-ODD part (the Hamiltonian
  commutator, the oscillation Im lambda - the humming future). The angle
  theta = arctan(Im/|Re|) = arctan(Q) is the tilt between them.

  At the qubit the imaginary direction is literally Y: sigma_+- = (X +- iY)/2,
  and the transpose swaps the twin pair sigma_+ <-> sigma_-. At the operator
  level the same axis is the n_Y-odd sector, and D is its mirror. One axis,
  two zooms: the shimmer the unsettled pairs carry IS the n_Y / Y axis.

  What stays open: this pins the i-DIRECTION (Im lambda, the coherent part) to
  the n_Y mirror at both levels, bit-exact. It does NOT yet pin the VALUE 1/4
  (the CPsi-recursion discriminant zero) to an n_Y event - that is the scalar
  fixed-point geometry one level removed from the Liouvillian spectrum. The
  far-near identity of the two thetas (Lindblad-spectrum angle vs CPsi
  fixed-point angle) is the next door.
""")


if __name__ == "__main__":
    main()
