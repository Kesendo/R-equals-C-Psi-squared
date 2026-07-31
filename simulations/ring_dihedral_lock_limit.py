#!/usr/bin/env python3
"""The ring dihedral-lock constant c_N: does N->inf give 1/sqrt(2)? No, it gives ln 2.

PROOF_RING_N4_DIHEDRAL_LOCK Section 5 shows the BOUND Im_max(L) <= DeltaE_max(H) at every N (the
dephasing dissipator is self-adjoint in the operator inner product, so it contributes nothing to
Im(lambda) = Im<v, -i[H,.]v>, which the Hamiltonian spread bounds; individual frequencies DO move
under the watching, they just cannot leave that interval). SATURATION -- that the bound is reached, so the
inequality is an equality -- is a separate step: proved at N=4 by the realising mode of Section 4
(K_{2,2}-specific), and MEASURED at N=6 (full Liouvillian, STAGE 1 below) and at N=8 (the Q=2 anchor
Im/sigma = 1.4128 = 0.7064*Q of F1_DISSIPATION_GAP_PATTERN). Where saturation holds,

    c_N = Im_max / (J*N) = (E_max - E_min) / (J*N) = 1/4 - E0/(J*N)

is a pure HAMILTONIAN question (2^N), not the full Liouvillian (4^N). E_max = J*N/4 is the
ferromagnet; E0 is the antiferromagnetic Heisenberg-ring ground state. For N >= 10 the table below is
therefore the Hamiltonian spread; reading it as Im_max carries the saturation assumption with it.
In the limit the per-bond ground energy is the Hulthen/Bethe value E0/(J*N) -> 1/4 - ln2, so

    c_inf = 1/4 - (1/4 - ln2) = ln 2 = 0.693147...   (NOT 1/sqrt(2) = 0.707107)

The EVEN sequence c_4=3/4, c_6=0.7171, c_8=0.7064 DECREASES through 1/sqrt(2) toward ln2.
1/sqrt(2) is only a value it passes through (the same red-herring lesson as s*=0.709). The ODD
rings are frustrated, so their AFM ground state sits higher and their spread is smaller: c_5=0.6236,
c_7=0.6579, c_9=0.6719 lie BELOW ln2 and RISE toward it. Two monotone branches closing on the same
limit from opposite sides, which is why reading the even rows alone as one descent is wrong. We
confirm by computing the AFM Heisenberg ring ground state E0(N) (sparse, minimal-S_z ground state)
and forming c_N.

  STAGE 0  the c_N sequence and its ln2 limit (the Hamiltonian side): the even branch N = 4..16
           and the odd branch N = 5..13, with the two-branch structure asserted.
  STAGE 1  saturation, i.e. the premise: max|Im| of the FULL Liouvillian vs DeltaE_max(H). N=4 by
           default (256x256), N=6 with --slow (4096x4096, a minute or two), each at two gammas, since
           saturation is what licenses reading the c_N table as Im_max at all.
  STAGE 2  the finite-N constants against their integer minimal polynomials in c (N = 6, 8, 10).

Run: python simulations/ring_dihedral_lock_limit.py [--slow]
"""
from __future__ import annotations

import sys
import math

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def heisenberg_ring_ground_energy(N: int) -> float:
    """E0 of H = (1/4) sum_bonds (X_iX_j + Y_iY_j + Z_iZ_j) on the N-site ring (J=1).
    Built sparsely by bit manipulation over the full 2^N space; ZZ diagonal, XY hops a
    differing pair with matrix element 1/2. eigsh finds the smallest (AFM ground) eigenvalue."""
    dim = 1 << N
    bonds = [(i, (i + 1) % N) for i in range(N)]
    rows, cols, vals = [], [], []
    diag = np.zeros(dim)
    for s in range(dim):
        d = 0.0
        for (i, j) in bonds:
            bi = (s >> i) & 1
            bj = (s >> j) & 1
            d += 0.25 if bi == bj else -0.25            # (1/4) Z_i Z_j
            if bi != bj:                                # (1/4)(X_iX_j+Y_iY_j) hops the differing pair
                t = s ^ (1 << i) ^ (1 << j)
                rows.append(t); cols.append(s); vals.append(0.5)
        diag[s] = d
    H = sp.coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()
    H = H + sp.diags(diag)
    e0 = eigsh(H, k=1, which="SA", return_eigenvectors=False, maxiter=5000)
    return float(e0[0])


def heisenberg_ring_dense(N: int) -> np.ndarray:
    """The same H = (1/4) sum_bonds (X_iX_j + Y_iY_j + Z_iZ_j), dense, for the Liouvillian check."""
    dim = 1 << N
    H = np.zeros((dim, dim))
    for s in range(dim):
        for (i, j) in [(k, (k + 1) % N) for k in range(N)]:
            bi, bj = (s >> i) & 1, (s >> j) & 1
            if bi == bj:
                H[s, s] += 0.25
            else:
                H[s, s] -= 0.25
                H[s ^ (1 << i) ^ (1 << j), s] += 0.5
    return H


def liouvillian_max_im(H: np.ndarray, N: int, g: float) -> float:
    """max|Im| over the whole Liouvillian spectrum: L = -i[H, .] + g*sum_l (Z_l . Z_l - .)."""
    dim = 1 << N
    Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    z = np.array([[1.0 if not ((s >> l) & 1) else -1.0 for s in range(dim)] for l in range(N)])
    for l in range(N):
        Zl = np.diag(z[l])
        L += g * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return float(np.abs(np.linalg.eigvals(L).imag).max())


def main() -> None:
    LN2 = math.log(2.0)
    INV_SQRT2 = 1.0 / math.sqrt(2.0)
    print("=" * 78)
    print("  STAGE 0 -- ring dihedral-lock constant c_N = 1/4 - E0/(J*N),  J = 1")
    print(f"  candidates: ln2 = {LN2:.6f}   1/sqrt(2) = {INV_SQRT2:.6f}")
    print("=" * 78)
    print(f"  {'N':>3} {'E0':>11} {'E0/N':>10} {'c_N':>10} {'c_N - ln2':>11} {'c_N - 1/√2':>11}")
    prev = None
    for N in range(4, 17, 2):
        e0 = heisenberg_ring_ground_energy(N)
        cN = 0.25 - e0 / N
        cross = ""
        if prev is not None and (prev - INV_SQRT2) * (cN - INV_SQRT2) < 0:
            cross = "  <- crosses 1/√2 here"
        print(f"  {N:>3} {e0:>11.5f} {e0/N:>10.5f} {cN:>10.5f} {cN-LN2:>+11.5f} {cN-INV_SQRT2:>+11.5f}{cross}")
        prev = cN
    print()
    print("  The ODD branch, which the even table above does not see (frustrated rings:")
    print("  a higher AFM ground state, so a SMALLER spread, and c_N BELOW ln2 rising):")
    prev_odd = None
    for N in range(5, 14, 2):
        e0 = heisenberg_ring_ground_energy(N)
        cN = 0.25 - e0 / N
        assert cN < LN2, f"STAGE 0 GATE FIRED at N={N}: odd-ring c_N {cN} is not below ln2"
        if prev_odd is not None:
            assert cN > prev_odd, f"STAGE 0 GATE FIRED at N={N}: odd branch is not rising"
        print(f"  {N:>3} {e0:>11.5f} {e0/N:>10.5f} {cN:>10.5f} {cN-LN2:>+11.5f} {cN-INV_SQRT2:>+11.5f}")
        prev_odd = cN
    print()
    print("  reading: TWO monotone branches, not one. The even c_N fall to ln2 = 0.6931 (the")
    print("  Hulthen per-bond AFM ground energy 1/4 - ln2) from above, passing THROUGH")
    print("  1/sqrt(2)=0.7071 on the way; the odd ones rise to it from below. The N->inf limit")
    print("  is ln2, not 1/sqrt(2); 1/sqrt(2) is a value the even branch crosses, like s*.")

    # ---------------------------------------------------------------- STAGE 1: the premise itself
    print()
    print("=" * 78)
    print("  STAGE 1 -- saturation: is max|Im|(L) really EQUAL to DeltaE_max(H), not just bounded?")
    print("  Section 5 of the proof gives the bound at every N; equality is proved at N=4 only.")
    print("=" * 78)
    sizes = (4, 6) if "--slow" in sys.argv else (4,)
    for N in sizes:
        H = heisenberg_ring_dense(N)
        e = np.linalg.eigvalsh(H)
        dE = float(e.max() - e.min())
        for g in (0.05, 0.5):
            mx = liouvillian_max_im(H, N, g)
            ok = abs(mx - dE) < 1e-6
            print(f"  N={N:>2}, gamma={g:>4}: max|Im|(L) = {mx:.6f}   DeltaE_max(H) = {dE:.6f}   "
                  f"c_N = {dE / N:.6f}   {'SATURATED' if ok else 'NOT saturated'}")
            assert ok, f"STAGE 1 GATE FIRED at N={N}, g={g}: max|Im| {mx} != DeltaE_max {dE}"
    if len(sizes) == 1:
        print("  (N=6 needs --slow: a 4096x4096 eigendecomposition, a minute or two.)")
    print("  STAGE 1 PASS: the bound is reached, so c_N reads as Im_max at the sizes checked here")
    print("  (N=4 proved + measured, N=6 measured with --slow, N=8 from the F1_DISSIPATION Q=2 anchor).")

    # ------------------------------------------------- STAGE 2: the finite-N constants are ALGEBRAIC
    # E0(N) is an eigenvalue of an integer matrix (4H has integer entries), so every finite-N c_N is an
    # algebraic number, not a transcendental one -- and at the sizes below it is a LOW-degree one with a
    # small integer minimal polynomial, obtained by factoring the exact characteristic polynomial of the
    # S_z=0 sector over the rationals. Only the LIMIT ln2 is transcendental.
    print()
    print("=" * 78)
    print("  STAGE 2 -- the finite-N constants in closed form (minimal polynomials over Z)")
    print("=" * 78)
    MINPOLY = {                                   # coefficients, highest power first
        6: ([12, -10, 1], "c_6 = (5 + sqrt(13))/12,  E0 = -(2 + sqrt(13))/2"),
        8: ([512, -640, 232, -25], "c_8 = the LARGEST root of 512c^3 - 640c^2 + 232c - 25"),
        10: ([200000, -420000, 342000, -136600, 27500, -2494, 67],
             "c_10 = the largest root of 200000c^6 - 420000c^5 + 342000c^4 - 136600c^3 + 27500c^2 - 2494c + 67"),
    }
    for N, (coeffs, label) in MINPOLY.items():
        cN = 0.25 - heisenberg_ring_ground_energy(N) / N
        roots = [r.real for r in np.roots(coeffs) if abs(r.imag) < 1e-9]
        near = min(roots, key=lambda r: abs(r - cN))
        print(f"  N={N:>2}: c_N = {cN:.15f}   nearest root = {near:.15f}   |diff| = {abs(near - cN):.2e}")
        print(f"        {label}")
        assert abs(near - cN) < 1e-9, f"STAGE 2 GATE FIRED at N={N}: c_N is not a root of its minimal polynomial"
        assert abs(max(roots) - cN) < 1e-9, (
            f"STAGE 2 GATE FIRED at N={N}: c_N is a root but not the LARGEST one, which is what "
            f"the documents state (roots {sorted(roots)})")
    c6 = 0.25 - heisenberg_ring_ground_energy(6) / 6
    assert abs(c6 - (5 + math.sqrt(13)) / 12) < 1e-12, "STAGE 2 GATE FIRED: c_6 != (5+sqrt13)/12"
    print("  STAGE 2 PASS: c_6 is a quadratic surd, c_8 a cubic and c_10 a sextic algebraic number.")
    print("  The per-N constants have closed forms; what needs Bethe/Hulthen is only the N->inf LIMIT.")


if __name__ == "__main__":
    main()
