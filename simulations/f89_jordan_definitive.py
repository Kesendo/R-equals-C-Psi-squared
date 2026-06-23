"""ADVERSARIAL CONFIRMATION: definitive Jordan-structure test for the
F89 path-3 octic EP (q_EP = sqrt((-1+sqrt13)/6)), independent of EpCharacter.

Route: generalized-eigenvector counting at the exact double-root lambda_EP.
  g1 = dim null(L - lam*I)          (geometric mult, via SVD, ABSOLUTE tol)
  g2 = dim null((L - lam*I)^2)
  DEFECTIVE (Jordan)  <=> g1 = 1 AND g2 = 2
  DIABOLIC            <=> g1 = 2 (two ordinary eigenvectors)

We build BOTH:
  (A) full un-symmetrized 24x24 (SE,DE) block  -- no S2 compression can hide a pair
  (B) 12x12 S2-symmetric sector               -- cross-check

Convention: gamma = 1, J = q (Q = q).  ASCII prints only.
"""

from __future__ import annotations

import sys
import numpy as np

# ASCII-only console
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

np.set_printoptions(precision=6, suppress=True, linewidth=160)


# --------------------------------------------------------------------------
# Independent construction of the (SE, DE) block.  Built from the physics
# (SE single-excitation hop matrix + DE two-excitation hop matrix on a 4-site
# chain), NOT copied verbatim, so a construction bug is not inherited.
# --------------------------------------------------------------------------
def build_full_block(J: float, gamma: float):
    """Return (L 24x24 complex, basis list of (i, (j,k))).

    Basis index = (SE site i in 0..3, DE pair jk among the 6 unordered pairs).
    SE block:  M_SE[a,b] = 2J if |a-b|==1   (nearest-neighbour hop on 4-chain)
    DE block:  hop one of the two excitations to an adjacent empty site.
    Liouvillian super-action: SE hop carries -i*M_SE, DE hop carries +i*M_DE
      (opposite signs = ket vs bra side of the commutator), plus a real
      diagonal -2gamma (SE site overlaps the DE pair) or -6gamma (no overlap).
    """
    de_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    basis = [(i, jk) for i in range(4) for jk in de_pairs]
    index = {b: n for n, b in enumerate(basis)}

    # SE nearest-neighbour hop on the 4-site chain
    M_SE = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J

    # DE two-excitation hop: move j or k to an adjacent empty site
    M_DE = np.zeros((6, 6))
    for col, (j, k) in enumerate(de_pairs):
        for nj in (j - 1, j + 1):
            if 0 <= nj <= 3 and nj != k:
                np_pair = tuple(sorted((nj, k)))
                if np_pair in de_pairs:
                    M_DE[de_pairs.index(np_pair), col] += 2 * J
        for nk in (k - 1, k + 1):
            if 0 <= nk <= 3 and nk != j:
                np_pair = tuple(sorted((j, nk)))
                if np_pair in de_pairs:
                    M_DE[de_pairs.index(np_pair), col] += 2 * J

    L = np.zeros((24, 24), dtype=complex)
    for col, (i, jk) in enumerate(basis):
        # SE side (ket): -i * M_SE acting on the SE index
        for i2 in range(4):
            if M_SE[i2, i] != 0:
                L[index[(i2, jk)], col] += -1j * M_SE[i2, i]
        # DE side (bra): +i * M_DE acting on the DE index
        jk_idx = de_pairs.index(jk)
        for jk2 in range(6):
            if M_DE[jk_idx, jk2] != 0:
                L[index[(i, de_pairs[jk2])], col] += 1j * M_DE[jk_idx, jk2]
        # dissipative diagonal
        L[col, col] += -2 * gamma if i in jk else -6 * gamma
    return L, basis


def reflection_projector(basis):
    """Return (P_plus 24xK+, P_minus 24xK-) onto the +-1 sectors of the
    chain reflection R: site s -> 3-s.  P-columns are orthonormal."""
    perm = {0: 3, 1: 2, 2: 1, 3: 0}
    refl_pair = lambda jk: tuple(sorted((perm[jk[0]], perm[jk[1]])))
    index = {b: n for n, b in enumerate(basis)}

    plus, minus = [], []
    handled = set()
    for n, (i, jk) in enumerate(basis):
        if n in handled:
            continue
        m = index[(perm[i], refl_pair(jk))]
        if n == m:
            v = np.zeros(24, dtype=complex)
            v[n] = 1.0
            plus.append(v)              # fixed point -> symmetric only
            handled.add(n)
        else:
            vp = np.zeros(24, dtype=complex)
            vp[n] = vp[m] = 1.0 / np.sqrt(2)
            vm = np.zeros(24, dtype=complex)
            vm[n] = 1.0 / np.sqrt(2)
            vm[m] = -1.0 / np.sqrt(2)
            plus.append(vp)
            minus.append(vm)
            handled.update((n, m))
    return np.column_stack(plus), np.column_stack(minus)


def sym_block(L, basis):
    """12x12 S2-symmetric (R = +1) sector, matching the locator's L_sym."""
    Pp, _ = reflection_projector(basis)
    return Pp.conj().T @ L @ Pp, Pp


# --------------------------------------------------------------------------
# Jordan test primitives
# --------------------------------------------------------------------------
def nullity(M, tol):
    """dim of numerical null space via SVD, ABSOLUTE singular-value tol."""
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s < tol)), s


def jordan_counts(L, lam, tol):
    n = L.shape[0]
    A = L - lam * np.eye(n)
    g1, s1 = nullity(A, tol)
    g2, s2 = nullity(A @ A, tol)
    return g1, g2, s1, s2


def null_basis(M, k, tol):
    """Return an orthonormal basis (columns) for the k smallest right
    singular vectors of M (the approximate null space)."""
    _, _, Vh = np.linalg.svd(M)
    return Vh.conj().T[:, -k:] if k > 0 else np.zeros((M.shape[1], 0), dtype=complex)


# --------------------------------------------------------------------------
# GATE: validate the g1/g2 test on KNOWN answers
# --------------------------------------------------------------------------
def gate():
    print("=" * 70)
    print("GATE  (validate g1/g2 nullity test on known 2x2 matrices)")
    print("=" * 70)
    lam = -4.0 + 1.3j
    tol = 1e-8
    Jb = np.array([[lam, 1.0], [0.0, lam]], dtype=complex)     # defective
    Db = np.array([[lam, 0.0], [0.0, lam]], dtype=complex)     # diabolic
    for name, M, exp in (("defective [[l,1],[0,l]]", Jb, (1, 2)),
                         ("diabolic  diag(l,l)", Db, (2, 2))):
        g1, g2, s1, s2 = jordan_counts(M, lam, tol)
        ok = (g1, g2) == exp
        print(f"  {name:24s}: g1={g1} g2={g2}  (expect g1={exp[0]} g2={exp[1]})  "
              f"{'PASS' if ok else 'FAIL'}")
        print(f"      sv(A)={np.round(s1,12)}   sv(A^2)={np.round(s2,12)}")
    print()


# --------------------------------------------------------------------------
def analyse(L, basis, label, q_ep, gamma=1.0):
    print("=" * 70)
    print(f"{label}   (dim = {L.shape[0]})")
    print("=" * 70)

    # reproduce EP: coalescing octic pair -> lam_EP
    eigvals = np.linalg.eigvals(L)
    rates = -eigvals.real / gamma
    is_at = (np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6)
    octic = list(eigvals[~is_at])
    pairs = sorted((abs(octic[i] - octic[j]), i, j)
                   for i in range(len(octic)) for j in range(i + 1, len(octic)))
    d_min, i0, j0 = pairs[0]
    lam_ep = 0.5 * (octic[i0] + octic[j0])
    print(f"  octic modes found        : {len(octic)}")
    print(f"  min coalescing pair dist : {d_min:.3e}")
    print(f"  lam_EP (numeric)         : {lam_ep:.10f}")
    print(f"  analytic -4g + 2iJ       : {complex(-4*gamma, 2*q_ep*gamma):.10f}")
    print(f"  |Re+4g|={abs(lam_ep.real+4*gamma):.2e}  |Im-2J|={abs(lam_ep.imag-2*q_ep*gamma):.2e}")

    # Use the ANALYTIC lambda for the Jordan test (exact double root), and
    # also the numeric midpoint -- report both so verdict is not lam-choice
    # dependent.
    lam_exact = complex(-4 * gamma, 2 * q_ep * gamma)
    for tag, lam in (("numeric-midpoint", lam_ep), ("analytic -4g+2iJ", lam_exact)):
        print(f"\n  -- Jordan test at lam = {lam:.8f}  [{tag}] --")
        # show several abs tolerances: nullity must be stable across a gap
        for tol in (1e-6, 1e-5, 1e-4, 1e-3):
            g1, g2, s1, s2 = jordan_counts(L, lam, tol)
            print(f"     tol={tol:.0e}:  g1={g1}  g2={g2}")
        # detailed singular values at a representative tol
        g1, g2, s1, s2 = jordan_counts(L, lam, 1e-5)
        print(f"     smallest sv of (L-lamI)   : {np.sort(s1)[:4]}")
        print(f"     smallest sv of (L-lamI)^2 : {np.sort(s2)[:4]}")
        verdict = ("DEFECTIVE (Jordan block)" if (g1 == 1 and g2 == 2)
                   else "DIABOLIC (g1=2, two eigenvectors)" if g1 >= 2
                   else f"AMBIGUOUS (g1={g1}, g2={g2})")
        print(f"     => {verdict}")

    # ---- eigenvector independence (Gram / rank / angle) at lam_EP ----
    A = L - lam_ep * np.eye(L.shape[0])
    V = null_basis(A, 2, 1e-5)           # 2 smallest right singular vectors
    # normalise
    V = V / np.linalg.norm(V, axis=0, keepdims=True)
    G = V.conj().T @ V
    sG = np.linalg.svd(G, compute_uv=False)
    rankG = int(np.sum(sG > 1e-6))
    # mutual angle of the two raw null vectors (before orthogonalising)
    v0, v1 = V[:, 0], V[:, 1]
    cosang = abs(np.vdot(v0, v1)) / (np.linalg.norm(v0) * np.linalg.norm(v1))
    print(f"\n  -- eigenvector independence at lam_EP --")
    print(f"     Gram singular values     : {np.round(sG, 6)}")
    print(f"     Gram rank (tol 1e-6)     : {rankG}  (2 => independent, 1 => coalesced)")
    print(f"     |<v0|v1>| of SVD null vec: {cosang:.4f}")

    # generalized-eigenvector existence: solve (L-lamI) x = v1 in least squares;
    # residual ~0 AND x not in null space  => generalized vector exists (defective)
    # If g1 already = 2, no generalized vector is needed; report the residual anyway.
    try:
        x, res, rk, sv = np.linalg.lstsq(A, v0, rcond=None)
        resid = np.linalg.norm(A @ x - v0)
        print(f"     lstsq resid (L-lamI)x=v0 : {resid:.3e}  (small => v0 in range(A) => chain exists)")
    except Exception as e:
        print(f"     lstsq failed: {e}")

    return lam_ep, octic, is_at, eigvals


def parity_sectors(L_full, basis, lam_ep, q_ep, gamma=1.0):
    """Confirm both coalescing eigenvectors sit in the SAME reflection sector."""
    print("=" * 70)
    print("PARITY: do the two coalescing eigenvectors share one R-sector?")
    print("=" * 70)
    Pp, Pm = reflection_projector(basis)
    # reflection operator as a matrix
    perm = {0: 3, 1: 2, 2: 1, 3: 0}
    refl_pair = lambda jk: tuple(sorted((perm[jk[0]], perm[jk[1]])))
    index = {b: n for n, b in enumerate(basis)}
    R = np.zeros((24, 24))
    for n, (i, jk) in enumerate(basis):
        R[index[(perm[i], refl_pair(jk))], n] = 1.0

    A = L_full - lam_ep * np.eye(24)
    V = null_basis(A, 2, 1e-5)
    for c in range(V.shape[1]):
        v = V[:, c] / np.linalg.norm(V[:, c])
        rv = R @ v
        # reflection eigenvalue ~ <v|R|v>
        par = np.vdot(v, rv).real
        proj_plus = np.linalg.norm(Pp @ (Pp.conj().T @ v))
        proj_minus = np.linalg.norm(Pm @ (Pm.conj().T @ v))
        print(f"  null vec {c}: <v|R|v> = {par:+.4f}   "
              f"||P+ v||={proj_plus:.4f}  ||P- v||={proj_minus:.4f}")
    print("  (both ~ +1 and concentrated in P+  => same symmetric sector)\n")


def bracket_scan(builder, q_ep, gamma=1.0):
    print("=" * 70)
    print("BRACKET SCAN: does g1 ever drop to 1 (defective) near q_EP?")
    print("=" * 70)
    print("   q          min-pair-dist   g1(@midpoint)  g2   verdict")
    for dq in (-0.01, -0.003, -0.001, -3e-4, 0.0, 3e-4, 0.001, 0.003, 0.01):
        q = q_ep + dq
        L, basis = builder(q * gamma, gamma)
        ev = np.linalg.eigvals(L)
        rates = -ev.real / gamma
        is_at = (np.abs(rates - 2.0) < 1e-6) | (np.abs(rates - 6.0) < 1e-6)
        oc = list(ev[~is_at])
        pr = sorted((abs(oc[i] - oc[j]), i, j) for i in range(len(oc)) for j in range(i+1, len(oc)))
        d, i0, j0 = pr[0]
        lam = 0.5 * (oc[i0] + oc[j0])
        g1, g2, _, _ = jordan_counts(L, lam, 1e-5)
        v = "DEFECTIVE" if (g1 == 1 and g2 == 2) else ("diabolic" if g1 >= 2 else "?")
        print(f"  {q:.6f}   {d:.3e}     g1={g1}        g2={g2}   {v}")
    print()


def main():
    q_ep = np.sqrt((-1 + np.sqrt(13)) / 6)
    gamma = 1.0
    print(f"q_EP = sqrt((-1+sqrt13)/6) = {q_ep:.10f}")
    print(f"check 3q^4+q^2-1 = {3*q_ep**4 + q_ep**2 - 1:.3e}  (should be 0)\n")

    gate()

    # (A) full un-symmetrized block
    L_full, basis = build_full_block(q_ep * gamma, gamma)
    lam_full, _, _, _ = analyse(L_full, basis, "(A) FULL un-symmetrized 24x24 block", q_ep, gamma)
    print()
    parity_sectors(L_full, basis, lam_full, q_ep, gamma)

    # (B) 12x12 symmetric sector
    L_sym, _ = sym_block(L_full, basis)
    analyse(L_sym, None, "(B) 12x12 S2-symmetric sector", q_ep, gamma)
    print()

    # bracket scans
    print("\n### FULL block bracket ###")
    bracket_scan(build_full_block, q_ep, gamma)

    print("### SYM block bracket ###")
    def sym_builder(J, g):
        Lf, bs = build_full_block(J, g)
        Ls, _ = sym_block(Lf, bs)
        return Ls, None
    bracket_scan(sym_builder, q_ep, gamma)


if __name__ == "__main__":
    main()
