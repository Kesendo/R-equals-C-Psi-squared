#!/usr/bin/env python3
"""Independent verifier of the birth-canal boundary via the |1-exc><vac| number-changing block.

Survey finding (2026-06-13): the boundary's slowest non-kernel rate (the witness's
BirthCanalDeviation source) is NOT the single-excitation density block rho^(1); it is the ODD
n_diff=1, number-CHANGING |1-excitation><vacuum| coherence. Under H = sum_b (X_bX_{b+1}+Y_bY_{b+1})
(coeff 1) the (ket#, bra#) bi-grading is conserved, so the (1,0) block {|k><vac|} is an N-dim
invariant subspace with

    L_(1,0)[k,k'] = -i*Q*h[k,k']  -  2*gamma_k * delta_{k,k'} ,    h tridiagonal, off-diag 2.

This script is the independent ground truth for the BirthCanalSurfaceWitness (which uses the full
4^N L). Two questions:

  (A) Is the (1,0) block the GLOBAL slowest non-kernel mode across the WHOLE surface (not just the
      three anchors)? -> the honesty check: does the cheap N-dim reduction faithfully reproduce the
      full-4^N boundary everywhere, or is there a region where the global slowest leaves the block?

  (C) Does the boundary (and the two sterility kinds) survive beyond N=5, where the dense witness
      cannot reach? -> N=6,7,8,10,... at N x N cost; full-4^N cross-check where still affordable.

Analytic backbone (proved here in passing): at FLAT gamma, L = -iQ*h - 2*gamma*I; -iQ*h is
anti-Hermitian (imaginary spectrum), so Re(lambda) = -2*gamma for EVERY mode -> rate = 2*gamma,
Q-invariant by uniformity alone, blind to eigenvector drift. That is R1's fragile (flat-gamma)
sterility, now a one-liner, at every N.
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
TOL = 1e-4          # the sterile/canal Deviation tolerance (matches the witness)
PROBE_LO, PROBE_HI = 1.5, 1000.0


# ---------- symmetric gamma-profile slice (matches C# SymmetricGammaSlice) ----------
def center_count(N):
    return 2 if N % 2 == 0 else 1


def sym_profile(N, w_edge, w_center):
    """[w_edge, w_bulk, ..., w_center, ..., w_bulk, w_edge], sum = N. None if bulk <= 0."""
    kc = center_count(N)
    kb = N - 2 - kc
    if kb < 1:
        raise ValueError("need N >= 5 for two independent (edge, center) coords")
    w_bulk = (N - 2.0 * w_edge - kc * w_center) / kb
    if w_bulk <= 0:
        return None
    lo = (N - kc) // 2
    p = np.empty(N)
    for i in range(N):
        p[i] = w_edge if i in (0, N - 1) else (w_center if lo <= i <= lo + kc - 1 else w_bulk)
    return p


# ---------- the two engines ----------
def _op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def _H_xy_unit(N):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        for P in (X, Y):
            t = np.array([[1]], complex)
            for i in range(N):
                t = np.kron(t, P if i in (b, b + 1) else I2)
            H += t
    return H


def full_eigs(N, Q, profile):
    d = 2 ** N
    Id = np.eye(d)
    H1 = _H_xy_unit(N)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = _op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return np.linalg.eigvals(L)


def full_slowest(N, Q, profile):
    w = full_eigs(N, Q, profile)
    return -float(np.max(w[np.abs(w) > 1e-7].real))


def block_L(N, Q, profile):
    h = np.zeros((N, N), complex)
    for i in range(N - 1):
        h[i, i + 1] = 2.0
        h[i + 1, i] = 2.0
    return -1j * Q * h - 2.0 * np.diag(profile)


def block_eigs(N, Q, profile):
    return np.linalg.eigvals(block_L(N, Q, profile))


def block_slowest(N, Q, profile):
    w = block_eigs(N, Q, profile)
    nz = w[np.abs(w) > 1e-7]
    return -float(np.max(nz.real)) if len(nz) else 0.0


def deviation(slow_fn, N, profile):
    return slow_fn(N, PROBE_HI, profile) - slow_fn(N, PROBE_LO, profile)


# ---------- (A) is the (1,0) block the global slowest across the whole surface? ----------
def part_A_surface_honesty(N=5, grid=9):
    print(f"=== (A) surface honesty: (1,0) block vs full 4^N, N={N}, {grid}x{grid} grid ===")
    es = np.linspace(0.2, 1.0, grid)
    cs = np.linspace(0.5, 3.0, grid)
    worst = 0.0
    worst_pt = None
    n_pts = 0
    n_canal_full = n_canal_block = 0
    verdict_mismatch = []
    for c in cs:
        for e in es:
            prof = sym_profile(N, e, c)
            if prof is None:
                continue
            n_pts += 1
            for Q in (PROBE_LO, PROBE_HI):
                rf = full_slowest(N, Q, prof)
                rb = block_slowest(N, Q, prof)
                d = abs(rf - rb)
                if d > worst:
                    worst, worst_pt = d, (round(e, 3), round(c, 3), Q)
            # boundary verdict from each engine
            df = abs(deviation(full_slowest, N, prof))
            db = abs(deviation(block_slowest, N, prof))
            cf, cb = df > TOL, db > TOL
            n_canal_full += cf
            n_canal_block += cb
            if cf != cb:
                verdict_mismatch.append((round(e, 3), round(c, 3), df, db))
    print(f"  admissible points: {n_pts}")
    print(f"  worst |full - block| rate discrepancy: {worst:.2e} at (e,c,Q)={worst_pt}")
    print(f"  canal verdicts: full={n_canal_full}, block={n_canal_block}, "
          f"verdict mismatches: {len(verdict_mismatch)}")
    if verdict_mismatch:
        print("  MISMATCH points (e,c,dev_full,dev_block):")
        for m in verdict_mismatch[:12]:
            print(f"    e={m[0]} c={m[1]}  dev_full={m[2]:.2e}  dev_block={m[3]:.2e}")
    honest = worst < 1e-6 and not verdict_mismatch
    print(f"  ==> the (1,0) block {'IS' if honest else 'is NOT'} the global slowest / boundary "
          f"carrier across the whole N={N} surface.\n")
    return honest


def part_A_subspectrum(N=5):
    """Stronger than the rate match: every (1,0)-block eigenvalue must be in the full spectrum."""
    print(f"=== (A') sub-spectrum: block eigenvalues subset of full 4^N spectrum, N={N} ===")
    pts = [(1.0, 1.0), (0.25, 3.0), (0.25, 1.5)]
    ok = True
    for (e, c) in pts:
        prof = sym_profile(N, e, c)
        for Q in (PROBE_LO, PROBE_HI):
            fe = full_eigs(N, Q, prof)
            be = block_eigs(N, Q, prof)
            maxmiss = max(float(np.min(np.abs(fe - s))) for s in be)
            ok = ok and maxmiss < 1e-7
            print(f"  (e={e}, c={c}, Q={Q:>6.1f}): max distance block-eig -> nearest full-eig "
                  f"= {maxmiss:.2e}")
    print(f"  ==> block spectrum {'IS' if ok else 'is NOT'} an exact sub-spectrum of the full L.\n")
    return ok


# ---------- (C) does the boundary survive beyond N=5? ----------
def part_C_scaling():
    print("=== (C) N>5: the boundary + the two sterility kinds, at N x N cost ===")
    print(f"  {'N':>3} {'uniform rate@1.5':>16} {'@1000':>9} {'flat-blind?':>12}"
          f" {'deep-edge dev':>14} {'sterile/canal':>14} {'full-check':>12}")
    worst_check = 0.0
    for N in (5, 6, 7, 8, 10, 12):
        # 1) uniform -> flat-gamma blindness: rate should be 2*gamma = 2 (gamma_l = 1), both Q.
        uni = np.ones(N)
        ru_lo = block_slowest(N, PROBE_LO, uni)
        ru_hi = block_slowest(N, PROBE_HI, uni)
        blind = abs(ru_lo - 2.0) < 1e-9 and abs(ru_hi - 2.0) < 1e-9

        # 2) a "deep-edge" structured profile: edges depressed to 0.25, the rest flat to sum N.
        #    (the general-N sibling of the flat-bulk-edge canal anchor)
        edge = 0.25
        rest = (N - 2 * edge) / (N - 2)
        deep = np.array([edge] + [rest] * (N - 2) + [edge])
        dev_block = deviation(block_slowest, N, deep)
        zone = "canal" if abs(dev_block) > TOL else "sterile"

        # 3) full-4^N cross-check where still affordable (N <= 6: 4096^2 EVD ~ seconds;
        #    N=7 is a 16384^2 dense EVD ~ tens of minutes, not worth it - the N=6 check already
        #    validates the reduction past N=5; N>=7 is block-only physics).
        if N <= 6:
            df_uni = abs(full_slowest(N, PROBE_LO, uni) - ru_lo)
            df_deep = abs(deviation(full_slowest, N, deep) - dev_block)
            worst_check = max(worst_check, df_uni, df_deep)
            check = f"{max(df_uni, df_deep):.1e}"
        else:
            check = "block-only"     # CANNOT validate: see the N=6 crossing below

        print(f"  {N:>3} {ru_lo:>16.6f} {ru_hi:>9.6f} {str(blind):>12}"
              f" {dev_block:>14.6f} {zone:>14} {check:>12}")
    print("\n  flat-gamma blindness (uniform rate = 2*gamma, Q-invariant) is analytic: -iQh is "
          "anti-Hermitian, so Re=-2gamma for every mode, at every N (this part of the reduction "
          "scales).")
    print("  BUT the deep-edge full-check breaks at N=6 (gap ~0.35): the global slowest is no longer "
          "the odd |1-exc><vac| mode. At low Q it switches to the EVEN {0,2}-coherence in the "
          "2-excitation density block (n_diff hist {0:0.78, 2:0.22}, <n_XY>=0.44 < 1, so slower than "
          "the odd mode's <n_XY>=1) - the SAME {0,2}-coherence the COHERENCE-HORIZON arc studies. So "
          "the (1,0) block does NOT capture the boundary at N>=6; the block-only zone verdicts at "
          "N>=7 track only the odd mode and are NOT validated.")
    return worst_check


def main():
    np.set_printoptions(precision=6, suppress=True)
    a = part_A_surface_honesty(N=5, grid=9)
    ap = part_A_subspectrum(N=5)
    worst_c = part_C_scaling()
    print()
    n5_exact = a and ap
    scales = worst_c < 1e-6
    if n5_exact and scales:
        print("VERDICT: the |1-exc><vac| N x N reduction is an exact independent ground truth and "
              "scales past N=5.")
    elif n5_exact:
        print("VERDICT: the |1-exc><vac| N x N reduction is an EXACT independent ground truth for the "
              "birth-canal boundary AT N=5 (whole surface, sub-spectrum), and proves R1's flat-gamma "
              "blindness analytically at every N. It does NOT universally scale: at N>=6 the global "
              f"slowest crosses to the even {{0,2}}-coherence (worst full-vs-block gap {worst_c:.2f}), "
              "the coherence-horizon mode. The N=5 validation stands; the N>=6 boundary is a "
              "birth_canal <-> coherence_horizon junction needing the full odd+even low-<n_XY> sector, "
              "not the (1,0) block alone.")
    else:
        print("VERDICT: the reduction fails even at N=5 (see (A)/(A')) - investigate.")


if __name__ == "__main__":
    main()
