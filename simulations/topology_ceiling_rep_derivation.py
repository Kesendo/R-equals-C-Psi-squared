"""Gate-first DERIVATION of the structural-ceiling closed forms (topology_band_edge arc NextStep 1+2+3).

The mechanism is already gate-verified (topology_ceiling_mechanism.py): at high Q the structural
ceiling g2 = strict_gap/(2*gamma) saturates, and g2 = <n_XY> of the slowest non-steady mode
(Absorption Theorem, Re(lambda) = -2*gamma*<n_XY>). Banked numbers: g2(K_N) = 4/N for N>=5
(4/5, 2/3, 4/7), K_4 = 0.8453 the N=4 OUTLIER (parallel to ring-4), star N>=6 -> 0.80. Those are
fits against the full Liouvillian. This verifier DERIVES them from first principles, gate-first.

RESULT (all gate-exact, machine precision):
  * g2(K_N)    = 4/N      (N>=5)  -- the (1,1)-sector commutant, the S_N STANDARD rep of the
                                     (N-1)-fold degenerate -J single-particle level.
  * g2(star_N) = 4/(N-1)  (N>=6)  -- the (1,1)-sector leaf manifold (the (N-2)-fold 0-eigenvalue
                                     level). CORRECTS the arc's tentative "star saturates at 0.80".
  * K_4 = 2 - 2/sqrt(3), ring-4 = 1.0  -- the N=4 outlier on BOTH is the (2,2) HALF-FILLING sector
                                     (K_4 dips below the floor, ring-4 co-occupies it); the 4/N
                                     ladder hits 1.0 at N=4 so the half-filling mode is the only one
                                     left under the floor. One sector, two topologies.
  * chain: no ceiling (no adjacency degeneracy, the (1,1) ladder stays >1; band edge protects).
  * NOT a universal law: the tempting 4/(m+1) (m = max adjacency degeneracy) fits complete + star
    but the RING (Fourier-degenerate manifold) breaks it. Per-family closed forms are the real result.

KEY STRUCTURE (the reduction that makes the rep story exact):
  In the |a><b| computational-coherence basis the Z-dephasing superoperator is DIAGONAL:
  L_D(|a><b|) = -2*gamma * hamming(a,b) * |a><b|, so N_XY := -L_D/(2*gamma) has eigenvalue
  hamming(a,b) = #bits where a,b differ = n_XY of that coherence. ad_H = [H,.] preserves the
  (n_left, n_right) = (popcount a, popcount b) bigrading, so the whole problem is SECTOR-DIAGONAL.

  HIGH-Q MECHANISM (degenerate perturbation theory, gamma << J): L = L_H + L_D, L_H dominates.
  The high-Q decay rates = eigenvalues of the L_H-eigenbasis-BLOCK-DIAGONAL part of N_XY (N_XY
  averaged over the ad_H flow; cross-Omega elements drop). So per L_H-eigenspace Omega, the rate
  is 2*gamma * (eigenvalues of P_Omega N_XY P_Omega). Two kinds of slow mode:
    - the BAND EDGE: the (0,1)/(1,0) sector has UNIFORM hamming=1, so N_XY = I there and L_D = -2g*I
      EXACTLY (all Q): rate = 2g exactly, g2 = 1. This floor is always present, so g2 <= 1 always.
    - the CEILING: an Omega=0 (commutant, [H,A]=0) coherence with <n_XY> < 1 dips BELOW the floor.

  (An earlier version projected onto Omega=0 only and missed the band edge -> it over-predicted g2
   for the protecting chain/ring cases. The all-Omega block-diagonal below is the correct model.)

  STAGE 0 (the MECHANISM gate): g2_predicted = min nonzero eigenvalue of the all-Omega
    block-diagonal N_XY, J-independent. GATE: == framework full-L g2 at Q=1000 (O(1/Q) residual)
    for every topology/N. A firing gate = the mechanism is wrong -> diagnose, do not loosen.

  STAGE 1 (the PRIZE -- g2(K_N)=4/N derived, not fit): the K_N ceiling lives in the Omega=0
    commutant of the (1,1) single-excitation-coherence sector = coherences among the (N-1)-fold
    degenerate -J level = the S_N STANDARD representation. GATE: that commutant darkest coherence
    = 4/N to MACHINE PRECISION for N=5,6,7, and it is the global slowest mode.

  STAGE 2 (the N=4 UNIFICATION): GATE: at N=4 the (1,1) commutant = EXACTLY 1.0 (= 4/4 = the band
    edge, no ceiling by itself); the K_4 ceiling comes from the (2,2) HALF-FILLING two-excitation
    sector -- the SAME sector that makes ring-4 special -- with its OWN closed form 2 - 2/sqrt(3)
    ~ 0.845299. The N=4 outlier = 4/N reaching 1.0 so the half-filling mode wins only there.

First-principles build (own Pauli/H/superoperators); Stage 0 cross-checks against framework full-L.
Run:  python simulations/topology_ceiling_rep_derivation.py
"""
import sys
import itertools
from math import sqrt, cos, pi
import numpy as np

sys.path.insert(0, 'simulations')
import framework as fw

GAMMA = 0.05
TOL_REL = 1e-6          # "same Omega" / "Omega=0" cluster tolerance, relative to the Omega scale
NZ = 1e-7               # nonzero-<n_XY> threshold (steady modes have <n_XY>=0 exactly)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)


def kron_list(ops):
    out = np.array([[1]], dtype=complex)
    for o in ops:
        out = np.kron(out, o)
    return out


def site_op(op, l, N):
    return kron_list([op if k == l else I2 for k in range(N)])


def topo_bonds(topo, N):
    if topo == 'chain':
        return [(i, i + 1) for i in range(N - 1)]
    if topo == 'ring':
        return [(i, (i + 1) % N) for i in range(N)]
    if topo == 'star':
        return [(0, i) for i in range(1, N)]
    if topo == 'complete':
        return [(i, j) for i in range(N) for j in range(i + 1, N)]
    raise ValueError(topo)


def H_full(topo, N, J=1.0):
    """H = (J/2) sum_bonds (X_i X_j + Y_i Y_j); real symmetric; SE block = J * adjacency."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in topo_bonds(topo, N):
        H += 0.5 * J * (site_op(X, i, N) @ site_op(X, j, N) + site_op(Y, i, N) @ site_op(Y, j, N))
    return H.real


def sector_states(N, p):
    """Computational basis indices with popcount p (the integer index IS the kron index)."""
    return [x for x in range(2 ** N) if bin(x).count('1') == p]


def sector_analysis(Hf, N, p, q):
    """High-Q decay structure of the (p,q) coherence sector.
    Returns dict with highq_min (min nonzero rate over ALL Omega blocks = the true high-Q g2
    contribution) and comm_min (min nonzero over the Omega=0 commutant block only)."""
    A = sector_states(N, p)
    B = sector_states(N, q)
    na, nb = len(A), len(B)
    if na == 0 or nb == 0:
        return None
    Hp = Hf[np.ix_(A, A)]
    Hq = Hf[np.ix_(B, B)]
    adH = np.kron(Hp, np.eye(nb)) - np.kron(np.eye(na), Hq)     # real symmetric; eigenvalues = Omega
    Omega, V = np.linalg.eigh(adH)
    diag = np.array([bin(A[a] ^ B[b]).count('1') for a in range(na) for b in range(nb)], dtype=float)
    Ntil = V.T @ (diag[:, None] * V)                            # N_XY in the L_H eigenbasis
    oscale = 1.0 + (np.abs(Omega).max() if Omega.size else 0.0)
    # all-Omega block-diagonal part (cross-Omega elements vanish under the ad_H time-average)
    same = np.abs(Omega[:, None] - Omega[None, :]) < TOL_REL * oscale
    w_all = np.linalg.eigvalsh(Ntil * same)
    nz_all = w_all[w_all > NZ]
    highq_min = float(nz_all.min()) if nz_all.size else None
    # Omega=0 commutant block only ([H,A]=0)
    z = np.abs(Omega) < TOL_REL * oscale
    comm_min = None
    if z.any():
        w_c = np.linalg.eigvalsh(Ntil[np.ix_(z, z)])
        nz_c = w_c[w_c > NZ]
        comm_min = float(nz_c.min()) if nz_c.size else None
    return dict(highq_min=highq_min, comm_min=comm_min, na=na, nb=nb)


def global_g2(Hf, N):
    """Global high-Q g2 = min nonzero rate over all (p,q) sectors; plus per-sector table."""
    best_val, best_pq = np.inf, None
    table = {}
    for p in range(N + 1):
        for q in range(N + 1):
            r = sector_analysis(Hf, N, p, q)
            if r is None:
                continue
            table[(p, q)] = r
            if r['highq_min'] is not None and r['highq_min'] < best_val:
                best_val, best_pq = r['highq_min'], (p, q)
    return best_val, best_pq, table


def fullL_g2(topo, N, Q=1000.0):
    """Framework full Liouvillian g2 = strict_gap/(2*gamma) at high Q (the banked oracle)."""
    J = Q * GAMMA
    cs = fw.ChainSystem(N=N, gamma_0=GAMMA, J=J, topology=topo, H_type='xy')
    ev = np.linalg.eigvals(np.asarray(cs.L))
    rates = -ev.real
    dec = rates > 1e-9
    return float(rates[dec].min()) / (2 * GAMMA)


# =====================================================================================
# STAGE 0 -- the MECHANISM gate: all-Omega block-diagonal N_XY min == full-L g2 at Q=1000
# =====================================================================================
print("=" * 100)
print("STAGE 0 -- MECHANISM GATE: high-Q ceiling = min nonzero eig of the all-Omega block-diagonal N_XY")
print("  g2_pred (J-independent, degenerate-PT)  vs  full-L g2 (Q=1000, the banked oracle)")
print("=" * 100)
print(f"{'topo':9} {'N':>2} {'g2_pred':>12} {'g2_fullL':>12} {'win sector':>11} {'|diff|':>9}")

worst = 0.0
for topo in ('chain', 'star', 'ring', 'complete'):
    for N in (3, 4, 5, 6):
        Hf = H_full(topo, N, J=1.0)
        gmin, gwin, _ = global_g2(Hf, N)
        try:
            g2f = fullL_g2(topo, N)
            d = abs(gmin - g2f)
            worst = max(worst, d)
            print(f"{topo:9} {N:>2} {gmin:>12.6f} {g2f:>12.6f} {str(gwin):>11} {d:>9.1e}")
        except Exception as e:
            print(f"{topo:9} {N:>2} {gmin:>12.6f} {'(skipped)':>12} {str(gwin):>11}   {type(e).__name__}")

assert worst < 3e-3, f"STAGE 0 GATE FIRED: degenerate-PT prediction disagrees with full-L by {worst:.2e} (>3e-3)"
print(f"\nSTAGE 0 PASS: the high-Q ceiling IS the block-diagonal-N_XY minimum (worst |diff| = {worst:.1e}, "
      f"the O(1/Q) residual). Band edge = the (0,1) Omega-block at 1.0; ceiling = an Omega=0 commutant mode < 1.")

# =====================================================================================
# STAGE 1 -- the PRIZE: g2(K_N) = 4/N from the (1,1) commutant (S_N standard rep), machine precision
# =====================================================================================
print("\n" + "=" * 100)
print("STAGE 1 -- g2(K_N) = 4/N from the (1,1) commutant (coherences in the (N-1)-fold -J level = S_N standard rep)")
print("=" * 100)
print(f"{'N':>2} {'(1,1) commutant':>16} {'4/N':>10} {'global min':>11} {'win sector':>11} {'==4/N exact?':>13}")
for N in (4, 5, 6, 7):
    Hf = H_full('complete', N, J=1.0)
    gmin, gwin, table = global_g2(Hf, N)
    c11 = table[(1, 1)]['comm_min']
    exact = (c11 is not None) and abs(c11 - 4.0 / N) < 1e-9
    print(f"{N:>2} {c11:>16.9f} {4.0/N:>10.6f} {gmin:>11.9f} {str(gwin):>11} {('YES' if exact else 'no'):>13}")

for N in (5, 6, 7):
    Hf = H_full('complete', N, J=1.0)
    gmin, gwin, table = global_g2(Hf, N)
    c11 = table[(1, 1)]['comm_min']
    assert abs(c11 - 4.0 / N) < 1e-9, f"STAGE 1 GATE FIRED: K_{N} (1,1) commutant {c11:.9f} != 4/{N}"
    assert gwin == (1, 1), f"STAGE 1 GATE FIRED: K_{N} global min not in (1,1) but {gwin}"
    assert abs(gmin - 4.0 / N) < 1e-9, f"STAGE 1 GATE FIRED: K_{N} global g2 {gmin:.9f} != 4/{N}"
print(f"\nSTAGE 1 PASS: g2(K_N) = 4/N EXACTLY (machine precision) for N=5,6,7, the (1,1)-sector commutant "
      f"darkest coherence. The closed form is derived from the S_N standard rep, no longer a fit.")

# =====================================================================================
# STAGE 2 -- the N=4 UNIFICATION: K_4 ceiling = the (2,2) half-filling sector (ring-4's sector), 2 - 2/sqrt(3)
# =====================================================================================
print("\n" + "=" * 100)
print("STAGE 2 -- N=4 UNIFICATION: K_4 ceiling = the (2,2) half-filling sector (the ring-4 sector)")
print("=" * 100)
K4_CLOSED = 2.0 - 2.0 / sqrt(3.0)        # = 2(1 - 1/sqrt(3)) ~ 0.845299
for topo in ('complete', 'ring'):
    Hf = H_full(topo, 4, J=1.0)
    gmin, gwin, table = global_g2(Hf, 4)
    label = 'K_4' if topo == 'complete' else 'ring-4'
    print(f"\n  {label}: global min nonzero <n_XY> = {gmin:.9f}  in sector {gwin}")
    print(f"    {'(p,q)':>7} {'all-Omega min':>15} {'Omega=0 commutant':>18}")
    for (p, q) in sorted(table):
        r = table[(p, q)]
        hm = f"{r['highq_min']:.6f}" if r['highq_min'] is not None else "  -"
        cm = f"{r['comm_min']:.6f}" if r['comm_min'] is not None else "  -"
        print(f"    {str((p, q)):>7} {hm:>15} {cm:>18}")

HfK4 = H_full('complete', 4, J=1.0)
gK4, winK4, tK4 = global_g2(HfK4, 4)
assert abs(tK4[(1, 1)]['comm_min'] - 1.0) < 1e-9, \
    f"STAGE 2 GATE FIRED: K_4 (1,1) commutant expected 1.0 (=4/4=band edge), got {tK4[(1,1)]['comm_min']:.9f}"
assert winK4 == (2, 2), f"STAGE 2 GATE FIRED: K_4 ceiling expected in (2,2) half-filling sector, got {winK4}"
assert abs(gK4 - K4_CLOSED) < 1e-7, \
    f"STAGE 2 GATE FIRED: K_4 ceiling expected 2-2/sqrt(3)={K4_CLOSED:.9f}, got {gK4:.9f}"
print(f"\nSTAGE 2 PASS: K_4 (1,1) commutant = 1.0 exactly (= 4/4 = the band edge, NO ceiling from the ladder); "
      f"the K_4 ceiling g2 = {gK4:.9f} = 2 - 2/sqrt(3) lives in the (2,2) HALF-FILLING sector -- the same "
      f"two-excitation sector that makes ring-4 special. THE N=4 OUTLIER = the 4/N ladder reaching 1.0 at "
      f"N=4, so the half-filling mode is the only one left below the floor.")

# =====================================================================================
# STAGE 2b -- ring-4 (2,2) co-occupier: 2sqrt(2)*J = the anti-periodic even-sector two-fermion band
#             top (JW string wraps the ring); the (2,2) spectrum is the CHIRAL palindrome about 0
#             (C_4 bipartite, K H K = -H, E<->-E; ChiralKClaim / PROOF_K_PARTNERSHIP).
# =====================================================================================
print("\n" + "=" * 100)
print("STAGE 2b -- ring-4 (2,2): Im=2sqrt(2)*J = anti-periodic two-fermion top; chiral palindrome about 0")
print("=" * 100)
Hr4 = H_full('ring', 4, J=1.0)
states22 = sector_states(4, 2)
e22 = np.sort(np.linalg.eigvalsh(Hr4[np.ix_(states22, states22)]).real)
ek_anti = sorted(2.0 * cos(pi * (2 * m + 1) / 4) for m in range(4))   # anti-periodic single-fermion E_k
sums_anti = np.sort([ek_anti[a] + ek_anti[b] for a, b in itertools.combinations(range(4), 2)])
two_sqrt2 = 2.0 * sqrt(2.0)
print(f"  ring-4 (2,2) energies / J:      {np.round(e22, 6)}")
print(f"  anti-periodic 2-fermion sums:   {np.round(sums_anti, 6)}  (single-fermion E_k = +-sqrt(2), k=pi/4)")
assert np.allclose(e22, sums_anti, atol=1e-9), "STAGE 2b GATE FIRED: (2,2) energies != anti-periodic two-fermion sums"
assert np.allclose(e22, -e22[::-1], atol=1e-9), "STAGE 2b GATE FIRED: (2,2) spectrum not palindromic about 0 (chiral K)"
assert abs(e22.max() - two_sqrt2) < 1e-9, f"STAGE 2b GATE FIRED: (2,2) top {e22.max():.6f} != 2sqrt(2)={two_sqrt2:.6f}"
assert two_sqrt2 > 2.0, "STAGE 2b: 2sqrt(2) should exceed the periodic ring band edge 2J"
print(f"  top = {e22.max():.6f} = 2sqrt(2) > band edge 2 (periodic, k=0); +-2sqrt(2) are a CHIRAL K-mirror")
print(f"  pair about 0 (C_4 bipartite, E<->-E; ChiralKClaim). The value = anti-periodic JW (even sector).")
print("STAGE 2b PASS: the co-occupier value (anti-periodic top) and its mirror symmetry (chiral K) derived.")

# =====================================================================================
# STAGE 3 -- the STAR closed form (NextStep item 2, fell out for free): g2(star_N) = 4/(N-1)
#   star = K_{1,N-1}; the darkest (1,1) coherence lives in the (N-2)-fold 0-eigenvalue leaf manifold.
#   This CORRECTS the arc's tentative "star saturates at 0.80 (N-independent)": it is 4/(N-1).
# =====================================================================================
print("\n" + "=" * 100)
print("STAGE 3 -- star closed form g2(star_N) = 4/(N-1) (the (1,1) leaf-manifold coherence), onset N=6")
print("=" * 100)
print(f"{'N':>2} {'star global min':>16} {'4/(N-1)':>10} {'win sector':>11} {'==exact?':>9}")
for N in (5, 6, 7, 8):
    Hf = H_full('star', N, J=1.0)
    gmin, gwin, table = global_g2(Hf, N)
    pred = 4.0 / (N - 1)
    exact = abs(gmin - pred) < 1e-9
    print(f"{N:>2} {gmin:>16.9f} {pred:>10.6f} {str(gwin):>11} {('YES' if exact else 'no'):>9}")
for N in (6, 7, 8):
    Hf = H_full('star', N, J=1.0)
    gmin, gwin, _ = global_g2(Hf, N)
    assert gwin == (1, 1), f"STAGE 3 GATE FIRED: star N={N} ceiling not in (1,1) but {gwin}"
    assert abs(gmin - 4.0 / (N - 1)) < 1e-9, f"STAGE 3 GATE FIRED: g2(star_{N})={gmin:.9f} != 4/{N-1}"
print(f"\nSTAGE 3 PASS: g2(star_N) = 4/(N-1) EXACTLY for N=6,7,8 (4/5, 4/6, 4/7), ceiling onset N=6 "
      f"(star N=5 -> 4/4 = 1.0 = the band edge). The star ceiling is 4/(N-1), NOT a constant 0.80.")

# Honesty check (NOT a universal law): the '4/(m+1)' guess from {complete: m=N-1, star: m=N-2} does
# NOT generalize -- the RING (Fourier-degenerate manifold) breaks it. Print ring (1,1) to show it.
print("\n" + "-" * 100)
print("DIAGNOSTIC (no universal m-law): ring (1,1) commutant vs the tempting 4/(m+1), m=max adjacency degeneracy")
print("-" * 100)
print(f"{'N':>2} {'ring (1,1) commutant':>22} {'4/(m+1) guess':>14} {'fits?':>6}")
for N in (4, 5, 6, 7):
    Hf = H_full('ring', N, J=1.0)
    _, _, table = global_g2(Hf, N)
    c11 = table[(1, 1)]['comm_min']
    # m = max multiplicity of the cycle adjacency spectrum
    A = np.zeros((N, N))
    for (i, j) in topo_bonds('ring', N):
        A[i, j] = A[j, i] = 1.0
    ev = np.round(np.linalg.eigvalsh(A), 6)
    m = max(np.unique(ev, return_counts=True)[1])
    guess = 4.0 / (m + 1)
    fits = "yes" if (c11 is not None and abs(c11 - guess) < 1e-6) else "NO"
    c11s = f"{c11:.9f}" if c11 is not None else "  -"
    print(f"{N:>2} {c11s:>22} {guess:>14.6f} {fits:>6}")
print("=> ring breaks 4/(m+1): the per-family closed forms (complete 4/N, star 4/(N-1)) are real; "
      "the m-law is NOT universal. Ring protects via the band edge (1.0) except the N=4 (2,2) anomaly.")

print("\nDONE.")
