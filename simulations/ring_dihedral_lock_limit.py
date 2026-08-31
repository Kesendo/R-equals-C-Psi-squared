#!/usr/bin/env python3
"""The ring dihedral-lock constant c_N: does N->inf give 1/sqrt(2)? No, it gives ln 2.

THE LOCK IS THE HAMILTONIAN SPREAD, AT EVERY N. F148 ("the imaginary reach is the
Hamiltonian spread, on every graph", docs/ANALYTICAL_FORMULAS.md) gives

    max|Im lambda_L|  =  DeltaE_max(H)      exactly, at every N and every gamma,

for isotropic Heisenberg under uniform Z-dephasing on the full single-site set {Z_l}.
Both halves are Tier 1, and neither is N=4's or the ring's:

  BOUND       the Z_l are Hermitian, so D is self-adjoint in the operator inner product
              and contributes nothing to Im(lambda) = Im<v, -i[H,.]v>, which the
              Hamiltonian spread bounds. Individual frequencies DO move under the
              watching, they just cannot leave that interval. Survives non-uniform gamma.
  ATTAINMENT  |ferro> = |0...0> is a Z-product state AND the maximum of H, and
              [H, sum_l Z_l] = 0 lets an E_min eigenvector be chosen inside a single
              Hamming rung k*. Then every basis pair of |ferro><ground| carries the same
              popcount(i^j) = k*, uniform dephasing acts on it as the scalar -2*gamma*k*,
              and lambda = -2*gamma*k* - i*DeltaE_max is an eigenvalue EXACTLY. Needs
              uniform gamma over the full {Z_l}, and the polarised extreme.

So no saturation assumption is carried at any N, and the lock reduces to

    c_N  =  Im_max / (J*N)  =  (E_max - E_min) / (J*N)  =  1/4 - E0/(J*N),

a pure HAMILTONIAN question (2^N), not the full Liouvillian (4^N). E_max = J*N/4 is the
ferromagnet, exactly: every bond term S_i.S_j has largest eigenvalue 1/4, so E <= J*N/4 on
the N-bond ring, and |ferro> attains all N of them at once. E0 is the antiferromagnetic
Heisenberg-ring ground state. In the limit the per-bond ground energy is the Hulthen/Bethe
value E0/(J*N) -> 1/4 - ln2, so

    c_inf  =  1/4 - (1/4 - ln2)  =  ln 2  =  0.693147...   (NOT 1/sqrt(2) = 0.707107)

The EVEN sequence c_4=3/4, c_6=0.7171, c_8=0.7064 DECREASES through 1/sqrt(2) toward ln2;
1/sqrt(2) is only a value it passes through, never the limit (the same lesson as the
birth-canal s* = 0.709: assume the pretty constant at your peril, compute it). The ODD
rings are frustrated, so their AFM ground state sits higher PER SITE and their per-site
spread is smaller: c_5=0.6236, c_7=0.6579, c_9=0.6719 lie BELOW ln2 and RISE toward it.
Two monotone branches closing on the same limit from opposite sides, which is why reading
the even rows alone as one descent is wrong.

  STAGE 0  the c_N ladder, both branches gated, and the headline gated by SHAPE rather
           than by a threshold: on the even branch c_N - ln2 ~ (pi^2/12)/N^2, the periodic
           c=1, v=pi/2 CFT finite-size form. It pins pi^2/12 only to about +-1.2% (the run
           prints the accepted window), but the refuted 1/sqrt(2) reading fails it outright,
           its column changing sign at N=8 and diverging. Printed alongside.
  STAGE 1  F148 from below. (a) ATTAINMENT at EVERY N of the ladder, via the
           |ferro><ground| certificate, which needs no 4^N object at all; (b) the BOUND,
           the full Liouvillian spectrum at N=4, and N=6 under --slow, each at two gammas,
           confirming that nothing else reaches higher; (c) the archived N=8 full-spectrum
           run, the one place the ring's whole 4^8 spectrum has ever been formed.
  STAGE 2  the finite-N constants against their integer minimal polynomials (N = 5, 6, 8,
           10), an odd rung among them: the odd branch is algebraic too.

Sibling gate: simulations/ring_n4_lock_gate.py (the N=4 closed form, F148's four scope
fences, and the same certificate at N=6, 8, 10). Proof: docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md.

Run: python simulations/ring_dihedral_lock_limit.py [--slow]
"""
from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EPS = np.finfo(float).eps
RING_N8_METRICS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "results", "f1_n8_n9_metrics", "ring_N8.json")


def _rung_states(N: int, p: int) -> np.ndarray:
    """The computational basis states of Hamming weight p."""
    return np.array([s for s in range(1 << N) if bin(s).count("1") == p], dtype=np.int64)


def heisenberg_ring_rung(N: int, p: int):
    """H = (1/4) sum_bonds (X_iX_j + Y_iY_j + Z_iZ_j) on the N-site ring (J=1), restricted
    to the Hamming rung p. H conserves popcount, so each rung is a block on its own: this
    is cheaper than the full 2^N space (at N=16, C(16,8)=12870 against 65536) and it is
    the reason a ground vector taken from one rung has a DEFINITE Hamming weight, which is
    exactly F148's attainment hypothesis. ZZ is diagonal; XY hops a differing pair with
    matrix element 1/2 and cannot leave the rung."""
    states = _rung_states(N, p)
    where = {int(s): i for i, s in enumerate(states)}
    dim = len(states)
    rows, cols, vals = [], [], []
    diag = np.zeros(dim)
    for c, s in enumerate(states):
        s = int(s)
        d = 0.0
        for i in range(N):
            j = (i + 1) % N
            bi, bj = (s >> i) & 1, (s >> j) & 1
            d += 0.25 if bi == bj else -0.25            # (1/4) Z_i Z_j
            if bi != bj:                                # (1/4)(X_iX_j+Y_iY_j) hops the pair
                rows.append(where[s ^ (1 << i) ^ (1 << j)]); cols.append(c); vals.append(0.5)
        diag[c] = d
    return sp.coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr() + sp.diags(diag)


def ring_ground(N: int, seed: int = 20260604) -> tuple[float, int, np.ndarray]:
    """(E0, k*, ground vector inside rung k*), by scanning the rungs. The global spin flip
    X^(x)N commutes with the isotropic H and maps rung p to rung N-p, so E0(p) = E0(N-p)
    and only p <= N/2 has to be solved; at odd N that mirror is why BOTH rungs (N-1)/2 and
    (N+1)/2 host the ground doublet, giving two realisers (below). eigsh gets a fixed
    start vector so the run is reproducible."""
    best_e, best_p, best_v = np.inf, -1, None
    rng = np.random.default_rng(seed)
    for p in range(N // 2 + 1):
        H = heisenberg_ring_rung(N, p)
        n = H.shape[0]
        if n <= 64:
            w, V = np.linalg.eigh(H.toarray())
            e, v = float(w[0]), V[:, 0]
        else:
            w, V = eigsh(H, k=1, which="SA", v0=rng.standard_normal(n), maxiter=50000)
            e, v = float(w[0]), V[:, 0]
        if e < best_e:
            best_e, best_p, best_v = e, p, v
    return best_e, best_p, best_v


def heisenberg_ring_ground_energy(N: int) -> float:
    """E0 of the AFM Heisenberg ring at J=1, the rung-resolved minimum."""
    return ring_ground(N)[0]


def ferro_energy_exact(N: int) -> float:
    """E_max, read off the p=0 rung. That rung is 1x1: |ferro> = |0...0> agrees on all N
    bonds and no XY term can hop out of it, so the entry is N*(1/4) with no eigensolver and
    no rounding (0.25 and its partial sums are exact in binary). That N/4 is the MAXIMUM is
    one line and needs no computation: each S_i.S_j has largest eigenvalue 1/4."""
    return float(heisenberg_ring_rung(N, 0).toarray()[0, 0])


def ulp_steps_to_exact_dE8(x: float) -> int:
    """How many doubles separate x from the EXACT DeltaE_max(8), with no eigensolver anywhere.

    c_8 is the largest root of 512c^3 - 640c^2 + 232c - 25 (STAGE 2). DeltaE_max(8) = 8*c_8,
    and substituting c = x/8 clears that cubic to the MONIC integer cubic

        Q(x) = x^3 - 10x^2 + 29x - 25,

    whose largest root DeltaE_max(8) is. Q is increasing there, so sign(Q(d)) places any
    double d against the true algebraic number, and over Fraction that sign is exact: no
    rounding on either side of the comparison. Walks from x until the sign flips."""
    def Q(d: float) -> Fraction:
        f = Fraction(d)
        return f ** 3 - 10 * f ** 2 + 29 * f - 25
    start = Q(x) < 0
    steps, d = 0, x
    while (Q(d) < 0) == start:
        d = np.nextafter(d, np.inf if start else -np.inf)
        steps += 1
    return steps


def certificate_collapsed(N: int, p: int, e0: float, v: np.ndarray) -> tuple[float, float]:
    """F148's realiser, checked without building any 4^N object.

    For M = |ferro><g| every term of L[M] stays rank-1 and collapses analytically:
      H@M   = (H|ferro>)<g|      = (N/4)*M          |ferro> is an exact H eigenvector,
      M@H   = |ferro>(H|g>)^dag  = E0*M             up to the ground eigenvector residual,
      Z_l M Z_l = |ferro>(Z_l|g>)^dag               since Z_l|ferro> = +|ferro> for all l,
      sum_l Z_l M Z_l = (N - 2k*)*M                 since |g> sits in rung k*,
    so D[M] = -2*gamma*k* * M and L[M] = (-2*gamma*k* - i*(N/4 - E0)) * M EXACTLY. The one
    numerical ingredient left is ||H|g> - E0|g>||, so the certificate's residual IS the
    ground eigenvector residual, and its error model is the eigensolver's own, eps*||H||.
    Returned as (residual, residual / (eps*||H||_1)) so the RATIO can be read across N
    rather than a bare threshold being asserted."""
    H = heisenberg_ring_rung(N, p)
    residual = float(np.linalg.norm(H @ v - e0 * v))
    scale = float(abs(H).sum(axis=0).max()) * EPS
    return residual, residual / scale


def heisenberg_ring_dense(N: int) -> np.ndarray:
    """The same H, dense on the full 2^N space, for the Liouvillian bound check."""
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


def certificate_direct(N: int, k_star: int, e0: float, g_full: np.ndarray, g: float) -> float:
    """The same certificate done literally, to check the collapse above is not wishful:
    build M = |ferro><g| in the 2^N x 2^N space and apply L term by term. Frobenius, the
    same norm certificate_collapsed uses, so the two are directly comparable and must agree.
    Affordable only while 2^N x 2^N is; the collapsed form is what runs at every N."""
    dim = 1 << N
    H = heisenberg_ring_dense(N)
    M = np.zeros((dim, dim))
    M[0, :] = g_full
    lam = -2.0 * g * k_star - 1j * (N / 4 - e0)
    LM = -1j * (H @ M - M @ H)
    for l in range(N):
        z = np.array([1.0 if not ((s >> l) & 1) else -1.0 for s in range(dim)])
        LM = LM + g * (z[:, None] * M * z[None, :] - M)
    return float(np.linalg.norm(LM - lam * M))


def liouvillian(H: np.ndarray, N: int, g: float) -> np.ndarray:
    """L = -i[H, .] + g*sum_l (Z_l . Z_l - .), dense on the 4^N operator space."""
    dim = 1 << N
    Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        zl = np.diag([1.0 if not ((s >> l) & 1) else -1.0 for s in range(dim)])
        L += g * np.kron(zl, zl)
    L[np.diag_indices(dim * dim)] -= g * N          # the -N*Id half of sum_l (Z_l . Z_l - .)
    return L


def main() -> None:
    LN2 = math.log(2.0)
    INV_SQRT2 = 1.0 / math.sqrt(2.0)
    CFT = math.pi ** 2 / 12.0          # c = 1, v = pi/2 on the periodic chain
    EVEN = list(range(4, 17, 2))
    ODD = list(range(5, 14, 2))

    print("=" * 78)
    print("  STAGE 0 -- ring dihedral-lock constant c_N = 1/4 - E0/(J*N),  J = 1")
    print(f"  candidates: ln2 = {LN2:.6f}   1/sqrt(2) = {INV_SQRT2:.6f}")
    print("=" * 78)
    print(f"  {'N':>3} {'E0':>11} {'E0/N':>10} {'c_N':>10} {'c_N - ln2':>11} {'c_N - 1/√2':>11}")
    c_even: dict[int, float] = {}
    prev = None
    for N in EVEN:
        e0 = heisenberg_ring_ground_energy(N)
        cN = 0.25 - e0 / N
        c_even[N] = cN
        cross = ""
        if prev is not None and (prev - INV_SQRT2) * (cN - INV_SQRT2) < 0:
            cross = "  <- crosses 1/√2 here"
        assert cN > LN2, f"STAGE 0 GATE FIRED at N={N}: even-ring c_N {cN} is not above ln2"
        if prev is not None:
            assert cN < prev, f"STAGE 0 GATE FIRED at N={N}: even branch is not falling"
        print(f"  {N:>3} {e0:>11.5f} {e0/N:>10.5f} {cN:>10.5f} {cN-LN2:>+11.5f} {cN-INV_SQRT2:>+11.5f}{cross}")
        prev = cN
    # "no N sits on 1/sqrt(2)" needs no tolerance: it is the strict chain c_6 > 1/sqrt2 > c_8
    # together with the monotone descent already asserted above.
    assert c_even[6] > INV_SQRT2 > c_even[8], (
        f"STAGE 0 GATE FIRED: 1/√2 is not crossed between N=6 and N=8 "
        f"({c_even[6]} , {INV_SQRT2} , {c_even[8]})")
    # c_4 = 3/4 is EXACT as a theorem (F149, the K_{2,2} Casimir gap 3J), but this route
    # reaches it through an eigensolver and so arrives with the eigensolver's error, not
    # exactly: measured 1 ulp low. Gate the deviation against that model, not equality.
    assert abs(c_even[4] - 0.75) <= 8 * np.spacing(0.75), (
        f"STAGE 0 GATE FIRED: c_4 is {c_even[4]}, off 3/4 by "
        f"{abs(c_even[4] - 0.75) / np.spacing(0.75):.1f} ulp")

    print()
    print("  The ODD branch, which the even table above does not see (frustrated rings: a")
    print("  higher AFM ground state PER SITE, so a smaller per-site spread, c_N BELOW ln2):")
    prev_odd = None
    for N in ODD:
        e0 = heisenberg_ring_ground_energy(N)
        cN = 0.25 - e0 / N
        assert cN < LN2, f"STAGE 0 GATE FIRED at N={N}: odd-ring c_N {cN} is not below ln2"
        if prev_odd is not None:
            assert cN > prev_odd, f"STAGE 0 GATE FIRED at N={N}: odd branch is not rising"
        print(f"  {N:>3} {e0:>11.5f} {e0/N:>10.5f} {cN:>10.5f} {cN-LN2:>+11.5f} {cN-INV_SQRT2:>+11.5f}")
        prev_odd = cN

    # --------------------------------------------- the headline, gated as a law not a number
    print()
    print("  The headline is not 'it looks like it tends to ln2'. On the periodic chain the")
    print("  c=1, v=pi/2 CFT finite-size form gives E0(N)/N = e_inf - pi^2/(12 N^2), hence")
    print(f"  c_N - ln2 ~ (pi^2/12)/N^2 with pi^2/12 = {CFT:.9f}. The refuted 1/sqrt(2)")
    print("  reading cannot satisfy any such law: its column changes sign and diverges.")
    print(f"  {'N':>3} {'N²(c_N - ln2)':>15} {'falsifier: N²(c_N - 1/√2)':>27}")
    prev_law = None
    for N in EVEN:
        law = N * N * (c_even[N] - LN2)
        assert law > CFT, f"STAGE 0 GATE FIRED at N={N}: N²(c_N - ln2) = {law} is below pi²/12"
        if prev_law is not None:
            assert law < prev_law, f"STAGE 0 GATE FIRED at N={N}: N²(c_N - ln2) is not falling to pi²/12"
        print(f"  {N:>3} {law:>15.6f} {N*N*(c_even[N]-INV_SQRT2):>27.6f}")
        prev_law = law
    # How fast the gap to pi²/12 closes is a model too, not a wish. The Heisenberg chain's
    # marginal operator leaves a logarithmic correction: gap ~ 1/ln³N predicts a factor
    # (ln16/ln4)³ = 8.0 across this table, against a measured 10.14, while a pure 1/N² tail
    # would give 16 and a constant offset 1. Bracket two-sided, so a stalling sequence AND a
    # suspiciously fast one both fire.
    gap_first, gap_last = EVEN[0] ** 2 * (c_even[EVEN[0]] - LN2) - CFT, prev_law - CFT
    shrink = gap_first / gap_last
    assert 5.0 < shrink < 20.0, (
        f"STAGE 0 GATE FIRED: the gap to pi²/12 shrank by {shrink:.2f} across "
        f"N = {EVEN[0]}..{EVEN[-1]}; the 1/ln³N model predicts 8.0, a 1/N² tail 16")
    print(f"  the gap to pi²/12 shrinks {shrink:.2f}x across the table (1/ln³N predicts 8.0).")
    print("  Honest about what this resolves: the three gates above accept any constant in")
    print("  (0.8114, 0.8311), so they pin pi²/12 to about ±1.2%. It is a SHAPE test that the")
    print("  1/sqrt(2) reading fails outright, not a precision measurement of the constant.")
    print(f"  falling toward {CFT:.6f} from above, never below it. The falsifier column")
    print("  crosses zero at N=8 and runs away: 1/sqrt(2) is passed through, not approached.")
    print()
    print("  reading: TWO monotone branches, not one. The even c_N fall to ln2 = 0.6931 (the")
    print("  Hulthen per-bond AFM ground energy 1/4 - ln2) from above, passing THROUGH")
    print("  1/sqrt(2)=0.7071 on the way; the odd ones rise to it from below.")

    # ------------------------------------------------------- STAGE 1: F148 from below
    print()
    print("=" * 78)
    print("  STAGE 1 -- F148 from below: max|Im|(L) = DeltaE_max(H), attainment and bound")
    print("=" * 78)
    print("  (a) ATTAINMENT. |ferro><ground| is an L eigenoperator at lambda = -2*gamma*k*")
    print("      - i*DeltaE_max, exactly, at every N and every gamma. Both inputs are exact:")
    print("      E_max = N/4 off the 1x1 p=0 rung, and |ground> sits in ONE rung by")
    print("      [H, sum_l Z_l] = 0. The residual below therefore IS the ground eigenvector")
    print("      residual (see certificate_collapsed), so it is READ against the")
    print("      eigensolver's own error model eps*||H||, not gated at a bare number.")
    print("      That residual certifies EIGEN-ness; that E0 is the MINIMUM is a separate")
    print("      claim, carried by the independent-seed rung rescan beside it.")
    print(f"  {'N':>3} {'k*':>3} {'DeltaE_max':>12} {'lambda(gamma=0.5)':>26} {'resid':>10} {'/eps||H||':>10}")
    for N in EVEN + ODD:
        e0, k_star, v = ring_ground(N)
        # |ferro> is the 1x1 p=0 rung, so N*(1/4) is exact in binary: a construction
        # guard on the rung builder, not a rounding test. It cannot fail arithmetically.
        e_max = ferro_energy_exact(N)
        assert e_max == N / 4, f"STAGE 1 GATE FIRED at N={N}: E_max is {e_max}, not exactly N/4"
        # MINIMALITY IS A SEPARATE CLAIM. The residual below certifies that (e0, v) is
        # an eigenpair, which ANY eigenpair satisfies; F148 needs the MINIMUM. ARPACK
        # missing the ground state is start-vector dependent, so re-run the whole rung
        # scan from an independent seed. The ground-to-first-excited gap on these rungs
        # is ~0.5 while seed-to-seed noise is ~2e-14, so any threshold between them
        # separates the two; 1e-9 sits nine orders below the gap and five above the noise.
        e0b, k_starb, _ = ring_ground(N, seed=91)
        assert abs(e0 - e0b) < 1e-9 and k_star == k_starb, (
            f"STAGE 1 GATE FIRED at N={N}: the rung scan is seed-dependent, E0 {e0} (rung "
            f"{k_star}) against {e0b} (rung {k_starb}); one of them is not the ground state")
        dE = e_max - e0
        resid, ratio = certificate_collapsed(N, k_star, e0, v)
        # measured 0.9 to 4.7 over the whole ladder at this seed, and 2.0 to 6.6 at N=16
        # across twelve ARPACK start vectors, so 100 leaves ~15x on the seed-varying worst
        # case and still fires long before anything structural could hide in it.
        assert ratio < 100, (f"STAGE 1 GATE FIRED at N={N}: certificate residual {resid} is "
                             f"{ratio:.1f} x eps*||H||, outside the eigensolver's error model")
        lam = -2 * 0.5 * k_star - 1j * dE
        print(f"  {N:>3} {k_star:>3} {dE:>12.8f} {str(np.round(lam, 8)):>26} {resid:>10.2e} {ratio:>10.1f}")
    print("      Counted by DISTINCT Re lambda, not by realiser: |1...1><g| carries the same")
    print("      lambda as |0...0><g| throughout. At ODD N the ground doublet lives in BOTH")
    print("      rungs (N∓1)/2 (spin-flip images), giving TWO values, -sigma ± gamma with")
    print("      sigma = N*gamma, exchanged by the palindrome lambda -> -conj(lambda) - 2*sigma.")
    print("      At EVEN N there is one value, Re lambda = -sigma: it sits ON the palindrome centre.")

    print()
    print("      cross-check that the rank-1 collapse is real: build M in the full 2^N x 2^N")
    print("      space, apply L term by term, and compare with the collapsed value. If the")
    print("      collapse holds these are the SAME number, so the check is their RATIO, not a")
    print("      threshold. Both are ~1e-15 reached by different arithmetic, so they agree to")
    print("      a factor rather than to digits; a broken collapse would be orders out.")
    for N in (4, 5, 6, 8):
        e0, k_star, v = ring_ground(N)
        full = np.zeros(1 << N)
        full[_rung_states(N, k_star)] = v
        collapsed, _ = certificate_collapsed(N, k_star, e0, v)
        for g in (0.05, 0.5, 2.0):
            direct = certificate_direct(N, k_star, e0, full, g)
            ratio = direct / collapsed
            assert 0.5 < ratio < 2.0, (
                f"STAGE 1 GATE FIRED at N={N}, g={g}: the direct residual {direct} and the "
                f"collapsed {collapsed} differ by {ratio:.3f}x; the rank-1 collapse is wrong")
        print(f"      N={N:>2}: direct {direct:.2e} vs collapsed {collapsed:.2e}, ratio {ratio:.4f}")

    print()
    print("  (b) BOUND. Nothing else in the 4^N spectrum reaches higher. The full Liouvillian,")
    print("      at N=4 (256x256) and, under --slow, N=6 (4096x4096), each at two gammas.")
    sizes = (4, 6) if "--slow" in sys.argv else (4,)
    for N in sizes:
        H = heisenberg_ring_dense(N)
        e = np.linalg.eigvalsh(H)
        dE = float(e.max() - e.min())
        for g in (0.05, 0.5):
            L = liouvillian(H, N, g)
            mx = float(np.abs(np.linalg.eigvals(L).imag).max())
            model = EPS * float(np.linalg.norm(L, 2))      # the eigensolver's own scale
            excess = mx - dE
            # measured 4.3 to 19.9 x eps*||L|| at N=4 and N=6, both gammas; 200 leaves ~10x.
            assert abs(excess) < 200 * model, (
                f"STAGE 1 GATE FIRED at N={N}, g={g}: max|Im| {mx} departs from DeltaE_max "
                f"{dE} by {excess}, which is {abs(excess)/model:.1f} x eps*||L||")
            print(f"      N={N:>2}, gamma={g:>4}: max|Im|(L) = {mx:.12f}   DeltaE_max = {dE:.12f}   "
                  f"excess = {excess:+.1e} = {abs(excess)/model:>5.1f} x eps*||L||")
    if len(sizes) == 1:
        print("      (N=6 needs --slow: a 4096x4096 eigendecomposition, a minute or two. It is a")
        print("       cross-check on (a), not a premise: (a) already covers N=4..16 above.)")

    print()
    print("  (c) The archived N=8 run, the one place the ring's whole 4^8 spectrum has been")
    print("      formed (F1 SLOW_N8 sweep, 2026-05-18, 65536 eigenvalues). Here an EXACT route")
    print("      exists, so no eigensolver appears on either side: DeltaE_max(8) = 8*c_8, and")
    print("      c = x/8 clears STAGE 2's cubic to the MONIC integer Q(x) = x^3 - 10x^2 + 29x")
    print("      - 25, whose largest root it is. Q evaluated over Fraction places the archived")
    print("      double against the true algebraic number with no rounding at all.")
    ran_archive = False
    if os.path.exists(RING_N8_METRICS):
        with open(RING_N8_METRICS, encoding="utf-8-sig") as fh:
            m = json.load(fh)
        assert m["N"] == 8 and m["JValue"] == 1 and m["GammaValue"] == 0.5, (
            f"STAGE 1 GATE FIRED: the archive is not the N=8, J=1, gamma=0.5 ring run "
            f"(N={m['N']}, J={m['JValue']}, gamma={m['GammaValue']})")
        archive = float(m["MaxImag"])
        steps = ulp_steps_to_exact_dE8(archive)
        sigma = m["N"] * m["GammaValue"]
        print(f"      MaxImag(4^8 Liouvillian) = {archive!r}")
        print(f"      exact DeltaE_max(8)        lies {steps} doubles above it, and nowhere else:")
        print(f"      Q is negative at that double and positive {steps} ulp up.")
        print(f"      Im/sigma = {archive/sigma:.10f} = c_8*Q at Q = J/gamma = "
              f"{m['JValue']/m['GammaValue']:.0f}, which the documents")
        print("      have been publishing to four digits as 1.4128 = 0.7064*Q.")
        assert steps <= 4, (
            f"STAGE 1 GATE FIRED: the archived MaxImag sits {steps} doubles from the exact "
            f"DeltaE_max(8); the 4^8 run and the closed form have parted company")
        ran_archive = True
    else:
        print(f"      (absent, so (c) did NOT run: {RING_N8_METRICS})")
    print("  STAGE 1 PASS: attainment certified at N = 4..16 even and 5..13 odd; the bound")
    print(f"  cross-checked on the full spectrum at N = {', '.join(str(n) for n in sizes)}"
          f"{'; the N=8 archive compared exactly' if ran_archive else '; (c) SKIPPED, archive absent'}.")

    # ------------------------------------------------- STAGE 2: the finite-N constants are ALGEBRAIC
    # E0(N) is an eigenvalue of an integer matrix (4H has integer entries), so every finite-N c_N is an
    # algebraic number, not a transcendental one -- and at the sizes below it is a LOW-degree one with a
    # small integer minimal polynomial, obtained by factoring the exact characteristic polynomial of the
    # S_z=0 sector over the rationals. Only the LIMIT ln2 is transcendental.
    print()
    print("=" * 78)
    print("  STAGE 2 -- finite-N constants in closed form (minimal polynomials over Z), one")
    print("  odd rung among them: the odd branch is algebraic too, it is not the even branch")
    print("  that is special. Only the N -> inf LIMIT needs Bethe/Hulthen.")
    print("=" * 78)
    MINPOLY = {                                   # coefficients, highest power first
        5: ([100, -80, 11], "c_5 = (4 + sqrt(5))/10,  the odd branch is algebraic too"),
        6: ([12, -10, 1], "c_6 = (5 + sqrt(13))/12,  E0 = -(2 + sqrt(13))/2"),
        8: ([512, -640, 232, -25], "c_8 = the LARGEST root of 512c^3 - 640c^2 + 232c - 25"),
        10: ([200000, -420000, 342000, -136600, 27500, -2494, 67],
             "c_10 = the largest root of 200000c^6 - 420000c^5 + 342000c^4 - 136600c^3 + 27500c^2 - 2494c + 67"),
    }
    # Budget for the comparison below (no exact route: an eigensolver on one side, a
    # polynomial root-finder on the other). eigsh contributes ~eps*||H||; np.roots on the
    # N=10 sextic contributes 62 ulp = 6.9e-15, measured by bracketing the exact root with
    # Fraction arithmetic on the integer polynomial. 1e-11 leaves ~1400x on the larger of the two and would still catch a wrong root, which is what this
    # is for: the nearest OTHER real root of the sextic sits 0.2419233 from c_10, and the
    # tightest pair anywhere in its root set is 0.0078194 apart.
    ROOT_BUDGET = 1e-11
    for N, (coeffs, label) in MINPOLY.items():
        cN = 0.25 - heisenberg_ring_ground_energy(N) / N
        roots = [r.real for r in np.roots(coeffs) if r.imag == 0.0]
        assert roots, f"STAGE 2 GATE FIRED at N={N}: np.roots returned no real root"
        largest = max(roots)
        others = sorted(abs(r - cN) for r in roots if abs(r - cN) > ROOT_BUDGET)
        print(f"  N={N:>2}: c_N = {cN:.15f}   largest root = {largest:.15f}   |diff| = {abs(largest - cN):.2e}")
        print(f"        {label}")
        assert abs(largest - cN) < ROOT_BUDGET, (
            f"STAGE 2 GATE FIRED at N={N}: c_N is not the LARGEST root of its minimal "
            f"polynomial, which is what the documents state (roots {sorted(roots)})")
        assert not others or others[0] > 1e-3, (
            f"STAGE 2 GATE FIRED at N={N}: another root sits {others[0]} away, inside the "
            f"margin this comparison relies on")
    # One side of this one is exact ((5+sqrt13)/12 in doubles), so the whole budget is the
    # eigensolver's: measured 2 ulp, gated at 64.
    c6 = 0.25 - heisenberg_ring_ground_energy(6) / 6
    c6_exact = (5 + math.sqrt(13)) / 12
    assert abs(c6 - c6_exact) <= 64 * np.spacing(c6_exact), (
        f"STAGE 2 GATE FIRED: c_6 is {abs(c6 - c6_exact) / np.spacing(c6_exact):.0f} ulp from "
        f"(5+sqrt13)/12")
    print("  STAGE 2 PASS: c_5 and c_6 are quadratic surds, c_8 a cubic, c_10 a sextic.")
    print("  The per-N constants have closed forms; what needs Bethe/Hulthen is only the N->inf LIMIT.")


if __name__ == "__main__":
    main()
