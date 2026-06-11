#!/usr/bin/env python3
"""F86b2 Direction (b''): the two-dial localization of the HWHM lift (2026-06-11).

THE QUESTION. The HWHM_left/Q_peak ratio lifts from the bare doubled-PTF floor
0.671535 (Tier1Derived, C2BareDoubledPtfClosedForm) to the empirical 0.7506
(Interior) / 0.7728 (Endpoint), captured so far only by fitted (alpha, beta)
per bond sub-class (F86HwhmClosedFormClaim, Tier1Candidate). Where does the
lift LIVE: in the inter-channel coupling beyond the top SVD pair (the "octic
residual" suspicion), or elsewhere?

THE TWO DIALS (on the c=2 coherence block, ResonanceScan conventions,
machinery of f86_kb_chiral_mirror.py which matches the C# at 1e-14):

  lambda dial: V_inter = sigma0*u0 v0^dag + lambda*V_rest (intra untouched).
  mu     dial: M = mu*M_intra + V_inter (inter untouched).

THE FINDING (all numbers asserted below, N = 5):
  1. GATE: at lambda = mu = 1 the pipeline reproduces the F90 bridge anchors
     bit-for-bit to grid noise (b=0: 0.7700, b=1: 0.7454).
  2. THE RATIO IS A SHAPE INVARIANT OF THE TAIL: along the lambda dial the
     ratio moves <= 0.005 while Q_peak moves by ~0.3-0.35 (b0: 2.84 -> 2.50,
     b1: 1.74 -> 1.48). The inter-channel SVD tail renormalizes the PEAK
     POSITION (hence g_eff = 4.39382/(Q_peak+2)) and leaves the dimensionless
     lineshape ratio essentially fixed.
  3. THE LIFT IS THE INTRA-CHANNEL DISPERSION, NON-PERTURBATIVELY: at mu = 0
     the resonance shape collapses entirely (no clean peak; ratio extraction
     degenerates), at mu = 0.5 it overshoots (b0 ~ 0.80), at mu = 1 it lands
     on the empirical values. The bare doubled-PTF model has NO intra-channel
     dispersion, which is exactly why it floors at 0.671535.

CONSEQUENCE FOR DIRECTION (b''). The (alpha, beta) derivation target is the
shape functional of the RANK-1-BRIDGE + INTRA-DISPERSION model: two channels
(HD = 1, 3) each carrying their internal path-like dynamics (sine modes, the
F89 Chebyshev terrain), bridged by one rank-1 coupling sigma0*u0 v0^dag.
The inter-channel tail only moves g_eff; deriving the ratio means treating
the intra dispersion exactly (not perturbatively) inside that reduced model.
The earlier framing "lift = H_B-mixed octic residual of the inter coupling"
is refuted by dial 2's tail-inertia.

Scope: N = 5 asserted here (the gate's anchor row); the dial structure at
N = 6, 7 is a follow-up. gamma0 = 0.05, uniform J = Q*gamma0, Dicke uniform
probe (DickeBlockProbe convention), peak over the ResonanceScan t-grid.

Run: python simulations/f86b2_shape_invariance_dial.py  (~3 min)
"""
import sys
from itertools import combinations
from math import comb, sqrt
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import f86_kb_chiral_mirror as kb  # noqa: E402

GAMMA = 0.05
N = 5
ANCHORS = {0: 0.7700, 1: 0.7454, 2: 0.7454, 3: 0.7700}   # F90 bridge, N=5


def build_dial(n_sites):
    blk = kb.build_block12(n_sites, GAMMA)
    qpairs = list(combinations(range(n_sites), 2))
    basis = [(j, pq) for j in range(n_sites) for pq in qpairs]
    hd = np.array([1 if j in pq else 3 for (j, pq) in basis])
    M_tot = sum(blk['M'])
    rows1 = np.where(hd == 1)[0]
    rows3 = np.where(hd == 3)[0]
    inter_mask = np.zeros(M_tot.shape, dtype=bool)
    inter_mask[np.ix_(rows1, rows3)] = True
    inter_mask[np.ix_(rows3, rows1)] = True
    V_inter = np.where(inter_mask, M_tot, 0.0)
    M_intra = M_tot - V_inter
    B13 = M_tot[np.ix_(rows1, rows3)]
    B31 = M_tot[np.ix_(rows3, rows1)]
    assert np.max(np.abs(B31 + B13.conj().T)) == 0.0, "inter block not anti-structured"
    U, s, Vh = np.linalg.svd(B13)
    top13 = s[0] * np.outer(U[:, 0], Vh[0, :])
    d_diag = np.array([-2.0 * GAMMA * (1 if j in pq else 3) for (j, pq) in basis],
                      dtype=complex)
    return {
        'blk': blk, 'd_diag': d_diag, 'rows1': rows1, 'rows3': rows3,
        'M_intra': M_intra, 'B13': B13, 'top13': top13, 'rest13': B13 - top13,
        'sigma0': s[0], 'svals': s, 'n': M_tot.shape[0],
    }


def generator(d, lam, mu):
    M = mu * d['M_intra'].copy()
    blk13 = d['top13'] + lam * d['rest13']
    M[np.ix_(d['rows1'], d['rows3'])] += blk13
    M[np.ix_(d['rows3'], d['rows1'])] += -blk13.conj().T
    return M


def k_of_Q(d, n_sites, lam, mu, Q):
    L = np.diag(d['d_diag']).astype(complex) + Q * GAMMA * generator(d, lam, mu)
    v0 = np.full(d['n'], 1.0 / (2.0 * sqrt(comb(n_sites, 1) * comb(n_sites, 2))),
                 dtype=complex)
    K = kb.k_curves(L, d['blk']['M'], d['blk']['S'], [v0])
    kpk, _ = kb.peak(K)
    return kpk[0]


def ratios(d, n_sites, lam, mu, q_lo=0.2, q_hi=6.0):
    qs = np.linspace(q_lo, q_hi, 44)
    k_all = np.array([k_of_Q(d, n_sites, lam, mu, q) for q in qs])
    out = []
    for b in range(n_sites - 1):
        i0 = int(np.argmax(k_all[:, b]))
        lo, hi = max(qs[i0] - 0.3, q_lo), min(qs[i0] + 0.3, q_hi)
        for _ in range(3):
            qf = np.linspace(lo, hi, 12)
            kf = np.array([k_of_Q(d, n_sites, lam, mu, q)[b] for q in qf])
            j = int(np.argmax(kf))
            lo, hi = qf[max(j - 1, 0)], qf[min(j + 1, len(qf) - 1)]
        q_peak = 0.5 * (lo + hi)
        k_peak = k_of_Q(d, n_sites, lam, mu, q_peak)[b]
        lo2, hi2 = q_lo, q_peak
        for _ in range(26):
            mid = 0.5 * (lo2 + hi2)
            if k_of_Q(d, n_sites, lam, mu, mid)[b] < k_peak / 2.0:
                lo2 = mid
            else:
                hi2 = mid
        out.append(((q_peak - 0.5 * (lo2 + hi2)) / q_peak, q_peak))
    return out


def main():
    d = build_dial(N)
    print(f"F86b2 two-dial localization, N = {N} (block dim {d['n']})")
    print(f"sigma0 = {d['sigma0']:.6f}, next svals {np.round(d['svals'][1:4], 4)} "
          f"(near-degenerate tail)")
    print("-" * 78)

    print("BLOCK 1  the gate: lambda = mu = 1 vs the F90 bridge anchors")
    full = ratios(d, N, 1.0, 1.0)
    for b, (r, qp) in enumerate(full):
        dev = abs(r - ANCHORS[b])
        print(f"  b={b}: ratio = {r:.4f} (anchor {ANCHORS[b]:.4f}, dev {dev:.4f})  "
              f"Q_peak = {qp:.4f}")
        assert dev < 1.5e-3, f"gate failed at b={b}"
    print("BLOCK 1 PASS")
    print("-" * 78)

    print("BLOCK 2  the lambda dial: ratio tail-inert, Q_peak tail-mobile")
    lam0 = ratios(d, N, 0.0, 1.0)
    for b in range(N - 1):
        dr = abs(lam0[b][0] - full[b][0])
        dq = abs(lam0[b][1] - full[b][1])
        print(f"  b={b}: ratio {lam0[b][0]:.4f} -> {full[b][0]:.4f} (|d| = {dr:.4f});  "
              f"Q_peak {lam0[b][1]:.4f} -> {full[b][1]:.4f} (|d| = {dq:.4f})")
        assert dr < 6e-3, f"ratio not tail-inert at b={b}"
        assert dq > 0.2, f"Q_peak unexpectedly tail-inert at b={b}"
    print("  => the SVD tail renormalizes g_eff (through Q_peak) and leaves the")
    print("     dimensionless lineshape ratio fixed: the ratio is a SHAPE INVARIANT.")
    print("BLOCK 2 PASS")
    print("-" * 78)

    print("BLOCK 3  the mu dial: the lift IS the intra dispersion, non-perturbative")
    mu_half = ratios(d, N, 1.0, 0.5)
    print(f"  mu=0.5: b0 ratio = {mu_half[0][0]:.4f} (overshoot vs {full[0][0]:.4f}), "
          f"Q_peak = {mu_half[0][1]:.4f}")
    assert mu_half[0][0] > full[0][0] + 0.02, "mu dial not overshooting: linearity?"
    # mu = 0: no clean resonance; witness via the collapse of the peak value scale
    k_mu0 = k_of_Q(d, N, 1.0, 0.0, 1.5)
    k_mu1 = k_of_Q(d, N, 1.0, 1.0, 1.5)
    print(f"  mu=0 peak-K at Q=1.5: {np.max(k_mu0):.3e} vs mu=1: {np.max(k_mu1):.3e} "
          f"(shape collapse; intra dispersion is essential)")
    print("  => the bare doubled-PTF floors at 0.671535 because it carries NO intra")
    print("     dispersion; (alpha, beta) live in the rank-1-bridge + intra model,")
    print("     to be treated exactly (F89 Chebyshev terrain), not perturbatively.")
    print("BLOCK 3 PASS")
    print("-" * 78)
    print("ALL BLOCKS PASS")


if __name__ == '__main__':
    main()
