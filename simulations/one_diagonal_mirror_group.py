"""The one diagonal is one of three: {Q_X, Q_Y, Q_Z} as one orbit -- the CORRECTED picture.

Deepens reflections/ON_THE_ONE_DIAGONAL.md via docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md.
History: the first run's gate FIRED and taught us the structure (2026-06-14). Two corrections it
forced (kept here as the dated lesson; "blind ones learn to see"):
  (1) THE Y-TRANSPOSE: the PHYSICAL dephasing diagonal is Q_P = Sum_l kron(P_l, P_l^T) (rho -> P rho P);
      since Y^T = -Y, Q_Y = -Sum kron(Y,Y), NOT +Sum kron(Y,Y). Same-spectrum held even with the wrong
      sign (the spectrum is symmetric about 0, so +-Q_Y are co-spectral) -- the gate separated SPECTRUM
      from OPERATOR, exactly its job.
  (2) THE PERMUTER: the group that permutes {Q_X,Q_Y,Q_Z} is the LETTER group of single-qubit basis
      moves, NOT <R, D, h>. D (the transpose) FIXES each diagonal (D Q D = +Q, the RATE reading); it
      does not permute them. The proof's "D = the Z<->Y swap" lives on the palindromizer Pi
      (D Pi_Z D = Pi_Y), not on the diagonal Q.
  (3) THE GROUP NAME and THE STRUCTURE (added 2026-08-01, found by review, not by this gate):
      <h_zx, h_yz> is NOT S3. R_x(pi/2) is a quarter-turn on the letters, so h_yz has order 4 and the
      closure has order 24 (the single-qubit Clifford group mod phase). Stage 1 printed that 24 all
      along under an "S3" label. The genuine letter-S3 has order 6 and needs the INVOLUTIVE
      transposition t_yz = Ad of (Y+Z)/sqrt(2); both groups induce the same orbit of 3, the order-6
      one faithfully, the order-24 one through a kernel of order 4 (Q_P is quadratic in P, hence
      sign-blind). And there is NO semidirect product S3 |x| D4: the letter moves do not NORMALIZE
      D4 (stage 3 below). The old stage 2 read the commutator pattern as evidence for a semidirect
      product; that inference is invalid, and relations of that shape hold in any group containing
      both factors.

THE CORRECTED PICTURE: two three-fold structures that do NOT lock into one group:
  * the THREE DIAGONALS {Q_X,Q_Y,Q_Z} = one orbit of the letter moves (h_zx: Z<->X, h_yz/t_yz: Z<->Y);
  * the THREE READINGS (rate/mirror/judge) = the mirror group D4 = <R, D> acting WITHIN one diagonal
    (D fixes = rate; R reflects R Q R = -Q = mirror, the -2 sum(gamma) shift; {D, FD} joint-fixed = judge).
The coherence-space closure <R, D, h_zx, t_yz> has order 96*2^N, not 48; see
docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md sect.5 "Resolution of the S3 side" and its gate
simulations/linear_s3_mirror_closure.py.

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
T_YZ = (Y1 + Z1) / np.sqrt(2.0)                          # the INVOLUTIVE Y<->Z transposition Clifford


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
    print(f"\nStage 1: {{Q_X,Q_Y,Q_Z}} = one orbit of the letter moves (N={N}):")
    d = 2 ** N
    QX, QY, QZ = (Q_diag(N, L) for L in ("X", "Y", "Z"))
    h_zx = ad_unitary(N, HAD)             # Z<->X (Hadamard)
    h_yz = ad_unitary(N, RX)              # Z<->Y (R_x(pi/2))

    # exact conjugators (the two basis transpositions)
    dev_zx = np.max(np.abs(h_zx @ QZ @ np.linalg.inv(h_zx) - QX))
    dev_yz = np.max(np.abs(h_yz @ QZ @ np.linalg.inv(h_yz) - QY))
    print(f"   exact conjugators: |h_zx Q_Z h_zx^-1 - Q_X| = {dev_zx:.2e}, "
          f"|h_yz Q_Z h_yz^-1 - Q_Y| = {dev_yz:.2e}")
    assert dev_zx < TOL and dev_yz < TOL, \
        "the basis moves do not realize Q_Z->Q_X / Q_Z->Q_Y to tolerance"

    # The two candidate letter groups. <h_zx,h_yz> is the order-24 Clifford group (R_x(pi/2) is a
    # quarter-turn, order 4); the genuine order-6 letter-S3 uses the involutive t_yz instead.
    t_yz = ad_unitary(N, T_YZ)
    dev_t = np.max(np.abs(t_yz @ QZ @ np.linalg.inv(t_yz) - QY))
    assert dev_t < TOL, "the involutive Y<->Z transposition does not carry Q_Z to Q_Y"
    for label, gens, want in (("<h_zx,h_yz>  (Clifford, R_x(pi/2) is order 4)", [h_zx, h_yz], 24),
                              ("<h_zx,t_yz>  (the genuine letter-S3)        ", [h_zx, t_yz], 6)):
        elems = group_closure(gens, d * d)
        orbit = orbit_of(QZ, elems)
        hits = {L: any(np.max(np.abs(Q_diag(N, L) - o)) < TOL for o in orbit) for L in ("X", "Y", "Z")}
        print(f"   |{label}| = {len(elems)} (expect {want}); "
              f"orbit(Q_Z) size = {len(orbit)}, contains X/Y/Z = {hits}")
        assert len(elems) == want, f"{label}: group order is not {want} -- the group is misnamed"
        assert len(orbit) == 3 and all(hits.values()), \
            f"{label}: the three diagonals are not one orbit"
    print("   [1] {Q_X,Q_Y,Q_Z} is one orbit of size 3 under BOTH letter groups; "
          "only the order-6 one is S3. OK")


def stage1b_orbit_n3():
    """attack at N+1: the orbit is a per-site basis permutation, so it must hold at N=3 too."""
    print("\nStage 1b (attack at N+1): the letter orbit holds at N=3:")
    N = 3
    QZ = Q_diag(N, "Z")
    h_zx = ad_unitary(N, HAD)
    h_yz = ad_unitary(N, RX)
    dev_zx = np.max(np.abs(h_zx @ QZ @ np.linalg.inv(h_zx) - Q_diag(N, "X")))
    dev_yz = np.max(np.abs(h_yz @ QZ @ np.linalg.inv(h_yz) - Q_diag(N, "Y")))
    print(f"   N=3 conjugators: |h_zx Q_Z h_zx^-1 - Q_X| = {dev_zx:.2e}, |h_yz Q_Z h_yz^-1 - Q_Y| = {dev_yz:.2e}")
    assert dev_zx < TOL and dev_yz < TOL, "N=3: basis conjugators off -- the orbit is not N-uniform"
    # the ORBIT CLOSURE at N=3, not just the two conjugators: the doc claims orbit size 3 at N=2,3.
    orbit = orbit_of(QZ, group_closure([h_zx, ad_unitary(N, T_YZ)], (2 ** N) ** 2))
    hits = {L: any(np.max(np.abs(Q_diag(N, L) - o)) < TOL for o in orbit) for L in ("X", "Y", "Z")}
    print(f"   N=3 orbit(Q_Z) under the letter-S3: size = {len(orbit)}, contains X/Y/Z = {hits}")
    assert len(orbit) == 3 and all(hits.values()), "N=3: the orbit is not exactly the three diagonals"
    print("   [1b] letter orbit confirmed at N=3 (per-site, N-uniform). OK")


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

    # The commutator pattern: each letter move commutes with one mirror generator and not the other.
    # This is a FACT about the two factors and NOT evidence of a semidirect product (stage 3 tests that).
    h_zx = ad_unitary(N, HAD)
    h_yz = ad_unitary(N, RX)
    c_zx_D = np.max(np.abs(h_zx @ D - D @ h_zx))
    c_zx_R = np.max(np.abs(h_zx @ R - R @ h_zx))
    c_yz_D = np.max(np.abs(h_yz @ D - D @ h_yz))
    c_yz_R = np.max(np.abs(h_yz @ R - R @ h_yz))
    print(f"   commutators: [h_zx,D]={c_zx_D:.1e} [h_zx,R]={c_zx_R:.1e} "
          f"[h_yz,D]={c_yz_D:.1e} [h_yz,R]={c_yz_R:.1e}")
    assert c_zx_D < TOL and c_zx_R > 0.1 and c_yz_R < TOL and c_yz_D > 0.1, \
        "the letter / mirror commutator pattern is not the expected one"
    print("   [2] rate=D-fix, mirror=R-anti, judge={D,FD} cell; commutator pattern as expected. OK")


# ==================== Stage 3: the letter moves do NOT normalize D4 ====================
def stage3_no_semidirect_product(N=3):
    """The gate the old stage 2 was missing. Until 2026-08-01 this script read the commutator
    pattern of stage 2 as 'so the structure is S3 |x| D4'. A semidirect product needs the letter
    moves to NORMALIZE D4; non-commutation is not that, and here the normalizer condition FAILS."""
    print(f"\nStage 3: the letter moves do NOT normalize D4, so there is no S3 |x| D4 (N={N}):")
    d = 2 ** N
    R = np.kron(np.eye(d, dtype=complex), xn(N))
    D = swap_perm(N)
    d4 = group_closure([R, D], d * d)
    assert len(d4) == 8, f"|<R,D>| = {len(d4)}, expected the 8-element D4"

    h_zx = ad_unitary(N, HAD)
    conj_R = h_zx @ R @ np.linalg.inv(h_zx)
    inside = any(np.max(np.abs(conj_R - g)) < TOL for g in d4)
    z_side = np.kron(np.eye(d, dtype=complex), kron_n([Z1] * N))
    dev_z = np.max(np.abs(conj_R - z_side))
    print(f"   |<R,D>| = {len(d4)}; h_zx R h_zx^-1 in <R,D>? {inside}")
    print(f"   h_zx R h_zx^-1 = one-sided multiplication by Z^(x)N: dev = {dev_z:.2e}")
    assert not inside, "h_zx R h_zx^-1 landed inside <R,D> -- the normalizer claim would need revisiting"
    assert dev_z < TOL, "h_zx R h_zx^-1 is not the one-sided Z^(x)N multiplication"
    print("   [3] the normalizer condition FAILS -> no semidirect product S3 |x| D4. "
          "The full closure (order 96*2^N, not 48) is gated in linear_s3_mirror_closure.py. OK")


if __name__ == "__main__":
    stage0_same_spectrum()
    stage1_orbit(N=2)
    stage1b_orbit_n3()
    stage2_readings_and_structure(N=3)
    stage3_no_semidirect_product(N=3)
    print("\nALL STAGES PASS: the one diagonal is one of three (a letter orbit of size 3), read three "
          "ways (mirror-D4); the two three-folds do NOT form a semidirect product S3 |x| D4.")
