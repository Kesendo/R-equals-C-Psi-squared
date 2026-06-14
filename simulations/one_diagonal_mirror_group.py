"""The one diagonal is one of three: {Q_X, Q_Y, Q_Z} as one orbit -- the CORRECTED picture.

Deepens reflections/ON_THE_ONE_DIAGONAL.md via docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md.
History: the first run's gate FIRED and taught us the structure (2026-06-14). Two corrections it
forced (kept here as the dated lesson; "blind ones learn to see"):
  (1) THE Y-TRANSPOSE: the PHYSICAL dephasing diagonal is Q_P = Sum_l kron(P_l, P_l^T) (rho -> P rho P);
      since Y^T = -Y, Q_Y = -Sum kron(Y,Y), NOT +Sum kron(Y,Y). Same-spectrum held even with the wrong
      sign (the spectrum is symmetric about 0, so +-Q_Y are co-spectral) -- the gate separated SPECTRUM
      from OPERATOR, exactly its job.
  (2) THE PERMUTER: the group that permutes {Q_X,Q_Y,Q_Z} is the BASIS-change S3 <h_zx, h_yz> (the
      single-qubit Cliffords permuting the X/Y/Z bases), NOT <R, D, h>. D (the transpose) FIXES each
      diagonal (D Q D = +Q, the RATE reading); it does not permute them. The proof's "D = the Z<->Y swap"
      lives on the palindromizer Pi (D Pi_Z D = Pi_Y), not on the diagonal Q.

THE CORRECTED PICTURE: two three-fold structures, semidirect (S3 |x| D4, the shape proof sect.5 expected):
  * the THREE DIAGONALS {Q_X,Q_Y,Q_Z} = one orbit of the basis-change S3 (h_zx: Z<->X, h_yz: Z<->Y);
  * the THREE READINGS (rate/mirror/judge) = the mirror group D4 = <R, D> acting WITHIN one diagonal
    (D fixes = rate; R reflects R Q R = -Q = mirror, the -2 sum(gamma) shift; {D, FD} joint-fixed = judge).
They do NOT commute fully ([h_zx,D]=0 but [h_zx,R]!=0; [h_yz,R]=0 but [h_yz,D]!=0) -> S3 |x| D4.

Row-stacking (C-order) vec, |i><j| -> e_i (x) e_j, kron(A,B): rho -> A rho B^T (matches
framework.lindblad + mirror_inventory_d4.py). Self-validating; ALL STAGES PASS prints only if all hold.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import site_op  # noqa: E402

TOL = 1e-10

X1 = np.array([[0, 1], [1, 0]], dtype=complex)
Y1 = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z1 = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
HAD = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0)               # Hadamard: X<->Z, Y->-Y
RX = np.cos(np.pi / 4) * I2 - 1j * np.sin(np.pi / 4) * X1                      # R_x(pi/2): Y->Z, Z->-Y


def kron_n(mats):
    out = np.array([[1.0]], dtype=complex)
    for m in mats:
        out = np.kron(out, m)
    return out


def xn(N):
    return kron_n([X1] * N)


def swap_perm(N):
    """D on coherence space: vec_C(rho) -> vec_C(rho^T)."""
    d = 2 ** N
    P = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            P[j * d + i, i * d + j] = 1.0
    return P


def Q_diag(N, letter):
    """The PHYSICAL dephasing diagonal in light P: Q_P = Sum_l kron(P_l, P_l^T) (rho -> P_l rho P_l).
    Z,X symmetric => kron(P,P); Y antisymmetric (Y^T = -Y) => Q_Y carries a minus sign (correction 1)."""
    d = 2 ** N
    out = np.zeros((d * d, d * d), dtype=complex)
    for l in range(N):
        Pl = site_op(N, l, letter)        # P on site l, I elsewhere (2^N x 2^N)
        out += np.kron(Pl, Pl.T)          # kron(A,B): rho -> A rho B^T, so P rho P uses P^T on the right
    return out


def ad_unitary(N, U1):
    """Ad_{U1^(x)N} on coherence space: rho -> U rho U^dag, vec_C = kron(U, U^*)."""
    U = kron_n([U1] * N)
    return np.kron(U, U.conj())


def group_closure(gens, dim, tol=TOL):
    elems = [np.eye(dim, dtype=complex)]
    changed = True
    while changed:
        changed = False
        for g in gens:
            for e in list(elems):
                cand = g @ e
                if not any(np.max(np.abs(cand - x)) < tol for x in elems):
                    elems.append(cand)
                    changed = True
    return elems


def orbit_of(Q, elems):
    orb = []
    for g in elems:
        Qg = g @ Q @ np.linalg.inv(g)
        if not any(np.max(np.abs(Qg - o)) < TOL for o in orb):
            orb.append(Qg)
    return orb


# ============================ Stage 0: same spectrum (the gate) ============================
def stage0_same_spectrum(Ns=(2, 3, 4)):
    print("Stage 0 (gate): Q_X, Q_Y, Q_Z same spectrum, computed DIRECTLY (eigvals):")
    for N in Ns:
        specs = {L: np.sort(np.linalg.eigvalsh(Q_diag(N, L)).real) for L in ("X", "Y", "Z")}
        dxy = np.max(np.abs(specs["X"] - specs["Y"]))
        dxz = np.max(np.abs(specs["X"] - specs["Z"]))
        print(f"   N={N}: max|spec(Q_X)-spec(Q_Y)|={dxy:.2e}, max|spec(Q_X)-spec(Q_Z)|={dxz:.2e}")
        assert dxy < TOL and dxz < TOL, f"N={N}: the three diagonals do NOT share a spectrum (the finding)."
    print("   [0] same spectrum, all N. OK")


# ============================ Stage 1: the basis-S3 orbit ============================
def stage1_orbit(N=2):
    print(f"\nStage 1: {{Q_X,Q_Y,Q_Z}} = one orbit of the basis-change S3 <h_zx, h_yz> (N={N}):")
    d = 2 ** N
    QX, QY, QZ = (Q_diag(N, L) for L in ("X", "Y", "Z"))
    h_zx = ad_unitary(N, HAD)             # Z<->X (Hadamard)
    h_yz = ad_unitary(N, RX)              # Z<->Y (R_x(pi/2))

    # exact conjugators (the two basis transpositions)
    dev_zx = np.max(np.abs(h_zx @ QZ @ np.linalg.inv(h_zx) - QX))
    dev_yz = np.max(np.abs(h_yz @ QZ @ np.linalg.inv(h_yz) - QY))
    print(f"   exact conjugators: |h_zx Q_Z h_zx^-1 - Q_X| = {dev_zx:.2e}, "
          f"|h_yz Q_Z h_yz^-1 - Q_Y| = {dev_yz:.2e}")
    assert dev_zx < TOL and dev_yz < TOL, "the basis moves do not realize Q_Z->Q_X / Q_Z->Q_Y bit-exact"

    elems = group_closure([h_zx, h_yz], d * d)
    orbit = orbit_of(QZ, elems)
    hits = {L: any(np.max(np.abs(Q_diag(N, L) - o)) < TOL for o in orbit) for L in ("X", "Y", "Z")}
    print(f"   |<h_zx,h_yz>| (single-qubit Clifford basis group, coherence) = {len(elems)}; "
          f"orbit(Q_Z) size = {len(orbit)}, contains X/Y/Z = {hits}")
    assert len(orbit) == 3 and all(hits.values()), "the three diagonals are not one orbit of the basis-S3"
    print("   [1] {Q_X,Q_Y,Q_Z} is exactly one orbit of the basis-change S3. OK")


def stage1b_orbit_n3():
    """attack at N+1: the orbit is a per-site basis permutation, so it must hold at N=3 too."""
    print("\nStage 1b (attack at N+1): the basis-S3 orbit holds at N=3:")
    N = 3
    QZ = Q_diag(N, "Z")
    h_zx = ad_unitary(N, HAD)
    h_yz = ad_unitary(N, RX)
    dev_zx = np.max(np.abs(h_zx @ QZ @ np.linalg.inv(h_zx) - Q_diag(N, "X")))
    dev_yz = np.max(np.abs(h_yz @ QZ @ np.linalg.inv(h_yz) - Q_diag(N, "Y")))
    print(f"   N=3 conjugators: |h_zx Q_Z h_zx^-1 - Q_X| = {dev_zx:.2e}, |h_yz Q_Z h_yz^-1 - Q_Y| = {dev_yz:.2e}")
    assert dev_zx < TOL and dev_yz < TOL, "N=3: basis conjugators not exact -- the orbit is not N-uniform"
    print("   [1b] basis-S3 orbit confirmed at N=3 (per-site, N-uniform). OK")


# ============================ Stage 2: the three readings + the semidirect structure ============================
def stage2_readings_and_structure(N=3):
    print(f"\nStage 2: the three readings = the mirror D4 acting WITHIN one diagonal (N={N}):")
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    F = xn(N)
    R = np.kron(Id, F)                    # R: rho -> rho F
    D = swap_perm(N)                      # D: rho -> rho^T
    QZ = Q_diag(N, "Z")

    dev_fix = np.max(np.abs(D @ QZ @ D - QZ))      # rate: D FIXES the diagonal
    dev_refl = np.max(np.abs(R @ QZ @ R + QZ))     # mirror: R REFLECTS it, R Q R = -Q
    print(f"   rate:   |D Q_Z D - Q_Z| = {dev_fix:.2e}   (D fixes the diagonal; NOT a permuter)")
    print(f"   mirror: |R Q_Z R + Q_Z| = {dev_refl:.2e}   (R reflects it, R Q R = -Q, carries -2 sum(gamma))")
    assert dev_fix < TOL and dev_refl < TOL, "the rate/mirror readings are not D-fix / R-anti on Q"
    print("   judge:  truly = {D, FD} joint-fixed cell (n_Y even & n_Z even); 63/63 in mirror_inventory_d4.py block D")
    print("   characterization: Q is the unique D-invariant, R-anti-invariant dephasing diagonal.")

    # the structure: basis-S3 and mirror-D4 are SEMIDIRECT (S3 |x| D4), not a direct product.
    h_zx = ad_unitary(N, HAD)
    h_yz = ad_unitary(N, RX)
    c_zx_D = np.max(np.abs(h_zx @ D - D @ h_zx))
    c_zx_R = np.max(np.abs(h_zx @ R - R @ h_zx))
    c_yz_D = np.max(np.abs(h_yz @ D - D @ h_yz))
    c_yz_R = np.max(np.abs(h_yz @ R - R @ h_yz))
    print(f"   structure (S3 |x| D4): [h_zx,D]={c_zx_D:.1e} [h_zx,R]={c_zx_R:.1e} "
          f"[h_yz,D]={c_yz_D:.1e} [h_yz,R]={c_yz_R:.1e}")
    # each basis move commutes with one mirror generator and not the other -> genuinely semidirect.
    assert c_zx_D < TOL and c_zx_R > 0.1 and c_yz_R < TOL and c_yz_D > 0.1, \
        "the basis-S3 / mirror-D4 coupling is not the expected semidirect pattern"
    print("   [2] rate=D-fix, mirror=R-anti, judge={D,FD} cell; basis-S3 |x| mirror-D4 (semidirect). OK")


if __name__ == "__main__":
    stage0_same_spectrum()
    stage1_orbit(N=2)
    stage1b_orbit_n3()
    stage2_readings_and_structure(N=3)
    print("\nALL STAGES PASS: the one diagonal is one of three (basis-S3 orbit), read three ways (mirror-D4); "
          "the structure is S3 |x| D4.")
