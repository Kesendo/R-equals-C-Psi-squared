#!/usr/bin/env python3
"""Edge 5 of the PTF fresh-eyes chain (2026-06-10): carrier-seam retrodiction.

QUESTION. reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md (2026-05-16) reports an
empirical break/hold list for the PTF closure under six bilinear Hamiltonians
(the multi-lens tour, simulations/multi_lens_ptf_carrier.py):
    breaks: IY+YI (Σ ln α = −6.91), YZ+ZY (−5.58), XZ+XZ (−0.52, asymmetric)
    holds:  XY+YX (−0.018), XZ+ZX (+0.018)
Is the whole list predicted by the perturbation's U(1)/Klein content alone, per
the Π-break dichotomy (experiments/PTF_PALINDROME_BREAKING_PERTURBATIONS.md:
closure rides on U(1); breaking U(1) lifts the F4 stationary sector-projectors)?

PROTOCOL FACTS read off the reflection's source (the tour script + captured
results, NOT the fresh-eyes summary):
  - The tour is a REPLACEMENT protocol, not a perturbative one: L_B carries the
    case's bilinear at full strength J = 1.0 INSTEAD of XX+YY (per bond,
    unit coefficient per term, so baseline H_A = Σ_b (XX+YY) at J=1).
  - N = 3 chain, |+−+⟩ product initial state (NOT the bonding mode),
    γ_deph = 0.1, PLUS T1 amplitude damping γ_T1 = 0.01 on every site,
    t ∈ [0, 8], 200 points, eigendecomposition propagation.
  - "Closure" THERE is the PTF α-fit closure: per-site purity
    P_B(i,t) ≈ P_A(i, α_i·t), Σ_i ln α_i ≈ 0, fit bounds α ∈ [0.1, 10].
    It IS the PTF α-closure (same observable family as the canonical PTF),
    but on a different configuration. Note: IY+YI's −6.91 = 3·ln(0.1), i.e.
    ALL THREE α pinned at the LOWER FIT BOUND - a saturation artifact reading
    "much slower / not a rescaling at all", not a meaningful sum.

THIS PROBE therefore runs BOTH readings:
  (R) replacement protocol, reproduced bit-for-bit from the committed tour;
  (C) canonical perturbative protocol (ptf_leak_scaling machinery: N = 5,
      uniform XY chain (J/2)(XX+YY), γ = 0.05, bonding-mode state, RK4,
      perturbation ε·V only in run B, ε = 0.1).

PREDICTOR LENSES (pure algebra): ΔN sector content of V in σ± language,
[V, N̂], N̂-parity [V, (−1)^N̂], Klein bits per term + homogeneity, hidden-U(1)
sublattice-flip scan (does some X/Y-flip frame restore [H, N̂] = 0), and the
F87 spectral-palindrome verdict via fw.classify_pauli_pair (a DIFFERENT axis,
included to show orthogonality).

SMOKING GUN (the F4 mechanism): stationary-mode count of the perturbed/replaced
Liouvillian (|λ| < 1e-10). Baseline N = 5 has exactly N+1 = 6 sector projectors.

RESULT (2026-06-10, banked in the reflection's "Seen again" note and the Π-break
experiment doc): the list is retrodicted by ONE lens, "the closure holds iff the
Liouvillian retains the complete N+1 stationary sector manifold (in any frame)":
  - canonical protocol: ALL five bilinears break (S = +8.3 to +11.3, vs J-defect
    control +0.097); both [V, N̂] ≠ 0 and #stationary < N+1 score 6/6, parity
    only 5/6. XY+YX (ΔN = ±2, parity-conserving) breaks FASTEST: partial
    survival of the parity projectors buys nothing, the fit needs ALL of them.
    The weaker-guardian (parity) hypothesis is REFUTED.
  - tour list: the stationary-manifold count on the replacement Liouvillian
    scores 4/5. XY+YX's tour-hold is a HIDDEN U(1): the sublattice X-flip
    (mask 01010) maps the pure XY+YX chain to a U(1)-conserving chain while
    commuting with Z-dephasing, so the full manifold survives (5/5 at N=4);
    the bare [V, N̂] lens misfires there, the frame-independent count does not.
    The one unretrodicted entry, XZ+ZX's tour-hold, is an |+−+⟩-specific
    trajectory coincidence (breaks hard from the bonding mode, S = +6.90), not
    a closure property. F87 verdicts are 'soft' for all five: the spectral
    palindrome is an orthogonal axis, exactly as the Π-break experiment said.
Caveats: break magnitudes are censored at the α fit bounds (they read "ansatz
dead", not a physical scale); the hidden-frame scan covered X/Y sublattice
flips, but the stationary count is frame-independent, so count = 1 rules out
any hidden full manifold.

Run: python "simulations/ptf_carrier_seam_retrodiction.py"   (~21 s)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw                              # noqa: E402
import ptf_transverse_field_pi_break as ptf         # noqa: E402  committed PTF machinery
import ptf_leak_scaling as pls                      # noqa: E402  committed Edge-2 protocol
import multi_lens_ptf_carrier as mlt               # noqa: E402  the reflection's source

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

EPS = 0.10          # canonical perturbative strength (reflection states none; tour is replacement)
GAMMA_C = 0.05      # canonical PTF γ
N_C = 5             # canonical N (baseline stationary count = N+1 = 6)
STAT_TOL = 1e-10

# The tour's cases (verbatim from multi_lens_ptf_carrier.py) and its captured
# verdicts (simulations/results/multi_lens_ptf_carrier.txt, 2026-05-16).
CASES = [
    ('XY+YX', [('X', 'Y'), ('Y', 'X')]),
    ('IY+YI', [('I', 'Y'), ('Y', 'I')]),
    ('YZ+ZY', [('Y', 'Z'), ('Z', 'Y')]),
    ('XZ+ZX', [('X', 'Z'), ('Z', 'X')]),
    ('XZ+XZ', [('X', 'Z'), ('X', 'Z')]),
]
TOUR_CAPTURED = {   # case -> (Sigma_ln_alpha, max_RMSE, verdict in the reflection)
    'XY+YX': (-0.0184, 0.0434, 'holds'),
    'IY+YI': (-6.9077, 0.0893, 'breaks (α all pinned at 0.1 bound)'),
    'YZ+ZY': (-5.5786, 0.0490, 'breaks'),
    'XZ+ZX': (+0.0176, 0.0428, 'holds'),
    'XZ+XZ': (-0.5213, 0.0949, 'breaks (small, spatially asymmetric)'),
}

LETTER_KLEIN = {'I': (0, 0), 'X': (1, 0), 'Y': (1, 1), 'Z': (0, 1)}


# ---------------------------------------------------------------------------
# Algebra lens helpers
# ---------------------------------------------------------------------------
def number_op_diag(N):
    """Diagonal of N̂ = Σ_l (I − Z_l)/2 in the computational basis = popcount."""
    idx = np.arange(2 ** N, dtype=np.uint64)
    pop = np.zeros(2 ** N, dtype=int)
    for b in range(N):
        pop += ((idx >> np.uint64(b)) & np.uint64(1)).astype(int)
    return pop


def delta_n_weights(V, N):
    """Frobenius weight of V per |ΔN| sector difference. Returns dict |k| -> share."""
    pop = number_op_diag(N)
    dk = pop[:, None] - pop[None, :]
    tot = float(np.sum(np.abs(V) ** 2))
    out = {}
    for k in range(0, N + 1):
        mask = (np.abs(dk) == k)
        w = float(np.sum(np.abs(V[mask]) ** 2))
        if w > 1e-24 * max(tot, 1.0):
            out[k] = w / tot
    return out


def comm_norm_with_diag(V, diag):
    """‖[V, D]‖_F for diagonal D given as a vector."""
    return float(np.linalg.norm(V * (diag[None, :] - diag[:, None])))


def klein_of_term(term):
    a = sum(LETTER_KLEIN[l][0] for l in term) % 2
    b = sum(LETTER_KLEIN[l][1] for l in term) % 2
    ypar = sum(1 for l in term if l == 'Y') % 2
    return a, b, ypar


def flip_scan(H, N, letter='X'):
    """Hidden-U(1) scan: does conjugation by a product of single-site X (or Y)
    flips restore [H', N̂] = 0?  X/Y flips exchange σ+ ↔ σ- locally; conjugation
    by the flip pattern `mask` permutes basis indices i -> i XOR mask (Y adds
    phases (-1)^... which cancel in |H'| sector masks; we conjugate exactly).
    Returns (best_mask, best_comm_norm)."""
    pop = number_op_diag(N)
    d = 2 ** N
    best = (None, np.inf)
    for mask in range(d):
        if letter == 'X':
            perm = np.arange(d) ^ mask
            Hp = H[np.ix_(perm, perm)]
        else:  # Y flips: Y|0>=i|1>, Y|1>=-i|0>  -> phase (-1)^(bits set in i&mask)
            perm = np.arange(d) ^ mask
            ph = np.array([(-1) ** bin(i & mask).count('1') for i in range(d)],
                          dtype=complex)
            Hp = (ph[:, None] * H[np.ix_(perm, perm)] * ph[None, :].conj())
        c = comm_norm_with_diag(Hp, pop.astype(float))
        if c < best[1]:
            best = (mask, c)
        if c < 1e-12:
            return mask, c
    return best


def stationary_count(L, tol=STAT_TOL, n_report=8):
    """Count eigenvalues of L with |λ| < tol; also return smallest |λ|s."""
    ev = np.linalg.eigvals(L)
    mags = np.sort(np.abs(ev))
    return int(np.sum(mags < tol)), mags[:n_report]


def build_V(N, terms, bonds=None):
    """V = Σ_bonds Σ_terms σσ at unit coefficient (the tour's convention)."""
    if bonds is None:
        bonds = [(i, i + 1) for i in range(N - 1)]
    return fw._build_bilinear(N, bonds, [(t[0], t[1], 1.0) for t in terms])


# ---------------------------------------------------------------------------
# Phase 1: algebra lenses
# ---------------------------------------------------------------------------
def phase1():
    print("=" * 96)
    print("PHASE 1  Predictor lenses (pure algebra, V = Σ_bonds term1+term2 at unit coeff, N=5 chain)")
    print("=" * 96)
    pop = number_op_diag(N_C).astype(float)
    par = (-1.0) ** pop
    chain4 = fw.ChainSystem(N=4, gamma_0=GAMMA_C, H_type='xy')

    sigma_reading = {
        'XY+YX': 'XY+YX = −2i(σ⁺σ⁺ − σ⁻σ⁻): pure ΔN = ±2',
        'IY+YI': 'IY+YI = Y_i + Y_j = −i(σ⁺−σ⁻)_i − i(σ⁺−σ⁻)_j: pure ΔN = ±1 (a Y-FIELD, no coupling)',
        'YZ+ZY': 'YZ+ZY = −i(σ⁺−σ⁻)_i Z_j + Z_i(−i)(σ⁺−σ⁻)_j: pure ΔN = ±1',
        'XZ+ZX': 'XZ+ZX = (σ⁺+σ⁻)_i Z_j + Z_i(σ⁺+σ⁻)_j: pure ΔN = ±1',
        'XZ+XZ': 'XZ+XZ = 2·X_i Z_j (directional): pure ΔN = ±1',
    }

    results = {}
    hdr = (f"  {'case':<8s} {'ΔN weights':<22s} {'‖[V,N̂]‖':>10s} {'‖[V,P̂]‖':>10s} "
           f"{'Klein(term1,term2)':<22s} {'homog':>6s} {'y_par':>6s} {'F87':>7s} "
           f"{'X-flip frame':>14s}")
    print(hdr)
    print('  ' + '-' * (len(hdr) - 2))
    for label, terms in CASES:
        V = build_V(N_C, terms)
        w = delta_n_weights(V, N_C)
        cN = comm_norm_with_diag(V, pop)
        cP = comm_norm_with_diag(V, par)   # diag of parity; [V,P] in same diag form
        kl = [klein_of_term(t) for t in terms]
        homog = (kl[0][:2] == kl[1][:2])
        ypars = {k[2] for k in kl}
        try:
            f87 = fw.classify_pauli_pair(chain4, terms)
        except Exception as e:  # 'I' letters etc.
            f87 = f"n/a({type(e).__name__})"
        # hidden-U(1): X-flip scan on the REPLACEMENT chain H = V (tour reading)
        mask, cbest = flip_scan(V, N_C, 'X')
        if cbest > 1e-12:
            masky, cbesty = flip_scan(V, N_C, 'Y')
            if cbesty < cbest:
                mask, cbest = masky, cbesty
        frame = (f"mask={mask:05b}" if cbest < 1e-12 else "none")
        wstr = ', '.join(f"|ΔN|={k}:{v:.3f}" for k, v in sorted(w.items()))
        print(f"  {label:<8s} {wstr:<22s} {cN:>10.3e} {cP:>10.3e} "
              f"{str(kl[0][:2]) + ',' + str(kl[1][:2]):<22s} {str(homog):>6s} "
              f"{str(sorted(ypars)):>6s} {str(f87):>7s} {frame:>14s}")
        results[label] = dict(dN=w, commN=cN, commP=cP, klein=kl, homog=homog,
                              f87=f87, flip=(mask, cbest))
        print(f"           σ±: {sigma_reading[label]}")
    print()
    print("  NOTE the F87 column is the SPECTRAL-palindrome trichotomy (a different axis from")
    print("  the closure); included to test orthogonality. Klein cells: Z=(0,1), X=(1,0), Y=(1,1).")
    print("  'X-flip frame' = sublattice-flip mask restoring [H, N̂]=0 for the REPLACEMENT chain")
    print("  H = V alone (the tour's reading); 'none' = no such frame in the X/Y-flip family.")
    print()
    return results


# ---------------------------------------------------------------------------
# Phase 2: stationary-mode counts (the F4 smoking gun)
# ---------------------------------------------------------------------------
def phase2():
    print("=" * 96)
    print(f"PHASE 2  Stationary modes |λ| < {STAT_TOL:g} of the Liouvillian")
    print("=" * 96)
    out = {'pert': {}, 'repl3_t1': {}, 'repl3': {}, 'repl4': {}}

    # --- (a) canonical perturbative, N = 5, γ = 0.05, no T1 ---
    print(f"  (a) CANONICAL perturbative: N={N_C}, uniform XY (J/2)(XX+YY), γ={GAMMA_C}, "
          f"H_B = H_A + {EPS}·V, no T1")
    H_A = ptf.build_H_xy(N_C, [1.0] * (N_C - 1))
    gl = [GAMMA_C] * N_C
    t0 = time.time()
    n0, m0 = stationary_count(fw.lindbladian_z_dephasing(H_A, gl))
    print(f"      baseline           : {n0} stationary   smallest |λ|: "
          + ' '.join(f"{x:.1e}" for x in m0))
    assert n0 == N_C + 1, f"baseline stationary count {n0} != {N_C + 1}"
    Jd = [1.0] * (N_C - 1)
    Jd[0] += EPS
    nJ, mJ = stationary_count(fw.lindbladian_z_dephasing(ptf.build_H_xy(N_C, Jd), gl))
    print(f"      J-defect control   : {nJ} stationary   smallest |λ|: "
          + ' '.join(f"{x:.1e}" for x in mJ))
    out['pert']['J-defect'] = nJ
    for label, terms in CASES:
        V = build_V(N_C, terms)
        n, m = stationary_count(fw.lindbladian_z_dephasing(H_A + EPS * V, gl))
        out['pert'][label] = n
        print(f"      +{EPS}·({label:<7s})    : {n} stationary   smallest |λ|: "
              + ' '.join(f"{x:.1e}" for x in m))
    print(f"      [wall {time.time() - t0:.0f} s]")
    print()

    # --- (b) replacement chains, N = 4, γ = 0.05, no T1 (size-robust, no tour T1) ---
    print(f"  (b) REPLACEMENT chains H = V at J=1 (tour reading), N=4, γ={GAMMA_C}, no T1 "
          f"(baseline count would be 5)")
    gl4 = [GAMMA_C] * 4
    H_xy4 = build_V(4, [('X', 'X'), ('Y', 'Y')])
    n, m = stationary_count(fw.lindbladian_z_dephasing(H_xy4, gl4))
    print(f"      XX+YY (baseline)   : {n} stationary")
    for label, terms in CASES:
        V = build_V(4, terms)
        n, m = stationary_count(fw.lindbladian_z_dephasing(V, gl4))
        out['repl4'][label] = n
        print(f"      {label:<7s}            : {n} stationary   smallest |λ|: "
              + ' '.join(f"{x:.1e}" for x in m[:6]))
    print()

    # --- (c) the tour's ACTUAL config: N = 3, γ_deph = 0.1 + γ_T1 = 0.01 ---
    print("  (c) the tour's ACTUAL config: N=3, γ_deph=0.1 PLUS γ_T1=0.01 (replacement)")
    for label, terms in [('XX+YY', [('X', 'X'), ('Y', 'Y')])] + CASES:
        L = mlt.build_L(3, [(0, 1), (1, 2)], terms, [0.1] * 3, [0.01] * 3, J=1.0)
        n, m = stationary_count(L)
        out['repl3_t1'][label] = n
        L0 = mlt.build_L(3, [(0, 1), (1, 2)], terms, [0.1] * 3, [0.0] * 3, J=1.0)
        n2, _ = stationary_count(L0)
        out['repl3'][label] = n2
        print(f"      {label:<7s}: with T1: {n} stationary | same config T1=0: {n2} "
              f"(baseline N+1 = 4)")
    print()
    print("  Reading: T1 collapses the stationary manifold to the single steady state for")
    print("  EVERY case - in the tour's own configuration the F4 sector-projector mechanism")
    print("  is structurally absent, yet the tour still found holds/breaks. So whatever the")
    print("  tour's closure rides on, in that config it is not the literal λ=0 count.")
    print()
    return out


# ---------------------------------------------------------------------------
# Phase 3: the closure itself, canonical perturbative protocol (ptf machinery)
# ---------------------------------------------------------------------------
def phase3():
    print("=" * 96)
    print(f"PHASE 3  PTF closure, canonical protocol (ptf_leak_scaling): N={N_C}, γ={GAMMA_C}, "
          f"bonding mode, ε={EPS} only in run B")
    print("=" * 96)
    bloch = ptf.site_bloch_ops(N_C)
    ham = ptf.hamming_matrix(N_C)
    psi0 = fw.bonding_mode_pair_state(N_C, 1)
    rho0 = np.outer(psi0, psi0.conj())
    H_A = ptf.build_H_xy(N_C, [1.0] * (N_C - 1))
    P_A = pls.purity_traj(N_C, H_A, rho0, bloch, ham, GAMMA_C)

    _, _, S_null = ptf.fit_all_sites(pls.T_GRID, P_A, P_A)
    print(f"  A-vs-A null: S = {S_null:+.2e}")

    out = {}
    Jd = [1.0] * (N_C - 1)
    Jd[0] += EPS
    runs = [('J-defect', ptf.build_H_xy(N_C, Jd))]
    runs += [(label, H_A + EPS * build_V(N_C, terms)) for label, terms in CASES]

    print(f"  {'case':<10s} {'α_i':<46s} {'S=Σ ln α':>10s} {'max RMSE':>9s}  bound-pinned?")
    for label, H_B in runs:
        P_B = pls.purity_traj(N_C, H_B, rho0, bloch, ham, GAMMA_C)
        al, rm, S = ptf.fit_all_sites(pls.T_GRID, P_A, P_B)
        pinned = [i for i, a in enumerate(al)
                  if abs(a - ptf.ALPHA_BOUNDS[0]) < 1e-3 or abs(a - ptf.ALPHA_BOUNDS[1]) < 1e-3]
        out[label] = dict(S=S, rmse=float(rm.max()), alphas=al, pinned=pinned)
        astr = ' '.join(f"{a:.3f}" for a in al)
        print(f"  {label:<10s} {astr:<46s} {S:>+10.4f} {rm.max():>9.1e}  "
              f"{('sites ' + str(pinned)) if pinned else 'no'}")
    print()

    # ε-scan on the discovery-sensitive case (ΔN=±2, parity-conserving)
    print(f"  ε-scan, XY+YX (the ΔN=±2 case) and YZ+ZY (a ΔN=±1 case) for contrast:")
    print(f"  {'case':<8s} {'ε':>6s} {'S':>10s} {'maxRMSE':>9s}")
    scan = {}
    for label in ('XY+YX', 'YZ+ZY'):
        terms = dict(CASES)[label]
        V = build_V(N_C, terms)
        for eps in (0.02, 0.05, 0.10, 0.20):
            P_B = pls.purity_traj(N_C, H_A + eps * V, rho0, bloch, ham, GAMMA_C)
            _, rm, S = ptf.fit_all_sites(pls.T_GRID, P_A, P_B)
            scan[(label, eps)] = (S, float(rm.max()))
            print(f"  {label:<8s} {eps:>6.2f} {S:>+10.4f} {rm.max():>9.1e}")
    print()
    out['scan'] = scan
    return out


# ---------------------------------------------------------------------------
# Phase 4: reproduce the reflection's own (replacement) protocol
# ---------------------------------------------------------------------------
def phase4():
    print("=" * 96)
    print("PHASE 4  Replacement protocol reproduced (the tour: N=3, |+−+⟩, γ=0.1, T1=0.01, J=1)")
    print("=" * 96)
    N = 3
    bonds = [(0, 1), (1, 2)]
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho0 = np.outer(psi, psi.conj())
    t_grid = np.linspace(0, 8.0, 200)
    gl, gt1 = [0.1] * N, [0.01] * N

    L_A = mlt.build_L(N, bonds, [('X', 'X'), ('Y', 'Y')], gl, gt1, J=1.0)
    P_A = mlt.per_site_purity_trajectory(L_A, rho0, N, t_grid)

    out = {}
    print(f"  {'case':<8s} {'S here':>9s} {'S captured':>11s} {'match':>6s} {'maxRMSE':>9s}  "
          f"reflection verdict")
    for label, terms in CASES:
        L_B = mlt.build_L(N, bonds, terms, gl, gt1, J=1.0)
        P_B = mlt.per_site_purity_trajectory(L_B, rho0, N, t_grid)
        al = np.zeros(N)
        rm = np.zeros(N)
        for i in range(N):
            a, r, _ = mlt.fit_alpha(t_grid, P_A[i], P_B[i])
            al[i], rm[i] = a, r
        S = float(np.sum(np.log(np.clip(al, 1e-30, None))))
        S_cap, rmse_cap, verdict = TOUR_CAPTURED[label]
        ok = abs(S - S_cap) < 5e-3
        out[label] = dict(S=S, rmse=float(rm.max()), match=ok)
        print(f"  {label:<8s} {S:>+9.4f} {S_cap:>+11.4f} {'OK' if ok else 'DIFF':>6s} "
              f"{rm.max():>9.1e}  {verdict}")
    print()
    return out


# ---------------------------------------------------------------------------
# Phase 4b: is the tour's "hold" pair structural or config-specific?
# ---------------------------------------------------------------------------
def phase4b():
    print("=" * 96)
    print("PHASE 4b  Tour-protocol variants: are the XY+YX / XZ+ZX 'holds' structural?")
    print("=" * 96)
    N = 3
    bonds = [(0, 1), (1, 2)]
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi_pmp = np.kron(plus, np.kron(minus, plus))
    psi_bond = fw.bonding_mode_pair_state(N, 1)
    psi_ppm = np.kron(plus, np.kron(plus, minus))
    t_grid = np.linspace(0, 8.0, 200)

    variants = [
        ('tour (|+−+⟩, T1)',   psi_pmp,  [0.1] * N, [0.01] * N),
        ('|+−+⟩, NO T1',       psi_pmp,  [0.1] * N, [0.0] * N),
        ('bonding mode, T1',   psi_bond, [0.1] * N, [0.01] * N),
        ('|++−⟩, T1',          psi_ppm,  [0.1] * N, [0.01] * N),
    ]
    out = {}
    for vlabel, psi, gl, gt1 in variants:
        rho0 = np.outer(psi, psi.conj())
        L_A = mlt.build_L(N, bonds, [('X', 'X'), ('Y', 'Y')], gl, gt1, J=1.0)
        P_A = mlt.per_site_purity_trajectory(L_A, rho0, N, t_grid)
        print(f"  variant: {vlabel}")
        print(f"    {'case':<8s} {'S':>9s} {'maxRMSE':>9s} {'eff=rms|P_B−P_A|':>17s}")
        for label, terms in CASES:
            L_B = mlt.build_L(N, bonds, terms, gl, gt1, J=1.0)
            P_B = mlt.per_site_purity_trajectory(L_B, rho0, N, t_grid)
            al = np.zeros(N)
            rm = np.zeros(N)
            for i in range(N):
                a, r, _ = mlt.fit_alpha(t_grid, P_A[i], P_B[i])
                al[i], rm[i] = a, r
            S = float(np.sum(np.log(np.clip(al, 1e-30, None))))
            eff = float(np.sqrt(np.mean((P_B - P_A) ** 2)))
            out[(vlabel, label)] = (S, float(rm.max()), eff)
            print(f"    {label:<8s} {S:>+9.4f} {rm.max():>9.1e} {eff:>17.4f}")
        print()
    return out


# ---------------------------------------------------------------------------
# Phase 5: verdict table + lens evaluation
# ---------------------------------------------------------------------------
def phase5(alg, stat, canon, repl):
    print("=" * 96)
    print("PHASE 5  VERDICT TABLE")
    print("=" * 96)
    print(f"  {'case':<8s} {'tour verdict':<34s} {'S_repl':>8s} {'S_canon':>9s} {'RMSE_c':>8s} "
          f"{'#stat(pert)':>11s} {'#stat(repl,N=4)':>15s} {'ΔN':>5s} {'parity':>7s} "
          f"{'Klein':>7s} {'F87':>7s}")
    print('  ' + '-' * 130)
    for label, terms in CASES:
        a = alg[label]
        dn = '/'.join(str(k) for k in sorted(a['dN']))
        parity_ok = a['commP'] < 1e-10
        kl = a['klein'][0][:2] if a['homog'] else 'mixed'
        print(f"  {label:<8s} {TOUR_CAPTURED[label][2]:<34s} {repl[label]['S']:>+8.3f} "
              f"{canon[label]['S']:>+9.4f} {canon[label]['rmse']:>8.1e} "
              f"{stat['pert'][label]:>11d} {stat['repl4'][label]:>15d} {dn:>5s} "
              f"{'cons' if parity_ok else 'broken':>7s} {str(kl):>7s} {str(a['f87']):>7s}")
    print(f"  {'J-def':<8s} {'holds (canonical control)':<34s} {'':>8s} "
          f"{canon['J-defect']['S']:>+9.4f} {canon['J-defect']['rmse']:>8.1e} "
          f"{stat['pert']['J-defect']:>11d}")
    print()

    print("  LENS EVALUATION (canonical perturbative protocol, break := |S| > 0.2 or RMSE > 1e-2):")
    for label in [l for l, _ in CASES]:
        S, rmse = canon[label]['S'], canon[label]['rmse']
        broke = (abs(S) > 0.2) or (rmse > 1e-2)
        a = alg[label]
        lens_a = a['commN'] > 1e-10                      # breaks full U(1)
        lens_b = stat['pert'][label] < N_C + 1           # F4: lost projectors
        lens_c = a['commP'] > 1e-10                      # ΔN-odd content
        print(f"    {label:<8s} measured: {'BREAK' if broke else 'hold ':<5s}  "
              f"(a) [V,N̂]≠0 → {'BREAK' if lens_a else 'hold '} "
              f"{'✓' if lens_a == broke else '✗'}   "
              f"(b) #stat<6 → {'BREAK' if lens_b else 'hold '} "
              f"{'✓' if lens_b == broke else '✗'}   "
              f"(c) parity-odd → {'BREAK' if lens_c else 'hold '} "
              f"{'✓' if lens_c == broke else '✗'}")
    print()
    print("  LENS EVALUATION (TOUR replacement protocol, verdicts as in the reflection):")
    tour_broke = {'XY+YX': False, 'IY+YI': True, 'YZ+ZY': True,
                  'XZ+ZX': False, 'XZ+XZ': True}
    for label in [l for l, _ in CASES]:
        broke = tour_broke[label]
        a = alg[label]
        lens_a = a['commN'] > 1e-10                       # bare [V,N̂] of the replaced H
        lens_b = stat['repl4'][label] < 5                 # full manifold of replacement L (N=4)
        lens_c = a['commP'] > 1e-10                       # parity
        print(f"    {label:<8s} tour: {'BREAK' if broke else 'hold ':<5s}  "
              f"(a) [V,N̂]≠0 → {'BREAK' if lens_a else 'hold '} "
              f"{'✓' if lens_a == broke else '✗'}   "
              f"(b) #stat(repl)<N+1 → {'BREAK' if lens_b else 'hold '} "
              f"{'✓' if lens_b == broke else '✗'}   "
              f"(c) parity-odd → {'BREAK' if lens_c else 'hold '} "
              f"{'✓' if lens_c == broke else '✗'}")
    print()
    print("  Note lens (b) on the replacement reading credits the hidden U(1): XY+YX alone")
    print("  has an X-flip frame (mask 01010) restoring [H, N̂]=0, hence the full N+1")
    print("  stationary manifold, hence 'hold'. XZ+ZX is the residual outlier on the tour list.")


def main():
    t0 = time.time()
    alg = phase1()
    stat = phase2()
    canon = phase3()
    repl = phase4()
    phase4b()
    phase5(alg, stat, canon, repl)
    print(f"\n[total wall time {time.time() - t0:.0f} s]")


if __name__ == "__main__":
    main()
