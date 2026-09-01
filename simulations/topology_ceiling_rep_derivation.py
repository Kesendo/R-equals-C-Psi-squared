"""Gate-first DERIVATION of the structural-ceiling closed forms (topology_band_edge arc NextStep 1+2+3).

The mechanism is already gate-verified (topology_ceiling_mechanism.py): at high Q the structural
ceiling g2 = strict_gap/(2*gamma) saturates, and g2 = <n_XY> of the slowest non-steady mode
(Absorption Theorem, Re(lambda) = -2*gamma*<n_XY>). Banked numbers: g2(K_N) = 4/N for N>=5
(4/5, 2/3, 4/7), K_4 = 0.8453 the N=4 OUTLIER (parallel to ring-4), star N>=6 -> 0.80. Those are
fits against the full Liouvillian. This verifier DERIVES them from first principles, gate-first.

MODEL SCOPE, load-bearing and not decoration: the XY network H = (J/2) sum_bonds (XX + YY) under
UNIFORM Z-dephasing; g2 is the ratio to that single gamma. The RING conclusion is exactly where the
XY hypothesis carries the weight. Under isotropic Heisenberg the same 4-cycle DOES ceiling, at the
K_4 value: measured 0.8452992906 at Q=1000 and 0.8452994599 at Q=10^4 against the closed form
2 - 2/sqrt(3) = 0.8452994616, while ring N=5 stays at 1 under both models. See PROOF_STRUCTURAL_CEILING
section 4 and docs/CAUGHT_ERRORS.md 2026-08-22. Every "the ring never ceilings" sentence below is an
XY sentence.

RESULT (all gate-exact, machine precision):
  * g2(K_N)    = 4/N      (N>=5)  -- the (1,1)-sector commutant. N_XY is S_N-equivariant and the
                                     commutant carries V (x) V* MULTIPLICITY-FREE, so by Schur it is
                                     a SCALAR on each isotypic piece: trivial -> 2/N, STANDARD -> 4/N,
                                     the remaining two -> 2, at every N>=4. That fixes the whole
                                     spectrum, not just the minimum, and the minimiser sits in the
                                     standard rep with multiplicity exactly N-1 (gated in STAGE 1).
                                     Careful: the commutant is End(V) (+) span{J}, ONE dimension more
                                     than End(V). Read as "the coherences in the -J level" ALONE it
                                     would give 2/N; the extra dimension is what turns that trivial
                                     mode into the exact zero mode and a 2, leaving 4/N as the
                                     minimum. The top level is load-bearing for the VALUE even though
                                     the minimiser has zero weight on it.
  * g2(star_N) = 4/(N-1)  (N>=4; a CEILING, i.e. < 1, only from N>=6) -- the (1,1)-sector leaf
                                     manifold (the (N-2)-fold 0-eigenvalue level), same Schur
                                     structure, same caveat. At N=3 the form does NOT hold: the value
                                     is 1, not 2. CORRECTS the arc's tentative "star saturates at 0.80".
  * K_4 = 2 - 2/sqrt(3), ring-4 = 1.0  -- the N=4 outlier on BOTH is the (2,2) HALF-FILLING sector
                                     (K_4 dips below the floor, ring-4 lands ON it); the 4/N ladder
                                     hits 1.0 = the floor ITSELF at N=4, so the half-filling mode is
                                     the only one left below it. One sector, two topologies.
  * chain: no ceiling. The (1,1) all-Omega minimum is 2(N-1)/(N+1), which is >= 1 with EQUALITY at
                                     N=3: chain-3 sits ON the band edge, it is not above it. The band
                                     edge protects.
  * NOT a universal law: the tempting 4/(m+1) (m = max adjacency degeneracy) fits complete + star
    but the RING (Fourier-degenerate manifold) breaks it. Per-family closed forms are the real result.

THE MARGINAL CASES ARE NOT FLAT (STAGE 4). Where a closed form equals exactly 1 (star N=5 at 4/4,
ring-4 at 2(N-2)/N), the mode sits ON the floor only in the Q -> infinity LIMIT. At finite Q it is
pulled BELOW the floor, and the pull is a law, measured over Q = 10 .. 3162:
        star-5:  g2 = 1 - 1/Q^2          ring-4:  g2 = 1 - (1/2)/Q^2
The (0,1) band edge itself never moves (STAGE 0a: exactly 1 at every Q). So "the ring has no ceiling"
is a statement about the LIMIT; at any finite Q the ring-4 gap does sit below the floor, by O(1/Q^2),
and the two marginal cases are ONE object, not a ring anomaly beside a star anomaly.

KEY STRUCTURE (the reduction that makes the rep story exact):
  In the |a><b| computational-coherence basis the Z-dephasing superoperator is DIAGONAL:
  L_D(|a><b|) = -2*gamma * hamming(a,b) * |a><b|, so N_XY := -L_D/(2*gamma) has eigenvalue
  hamming(a,b) = #bits where a,b differ = n_XY of that coherence. (Convention: L_l = sqrt(gamma)*Z_l,
  i.e. D[rho] = gamma * sum_l (Z_l rho Z_l - rho); a different jump normalisation rescales the 2.)
  ad_H = [H,.] preserves the (n_left, n_right) = (popcount a, popcount b) bigrading, so the whole
  problem is SECTOR-DIAGONAL.

  THE SECTOR LATTICE HAS EXACT DEGENERACIES, and every gate below is written around that. X^N maps
  the (p,q) sector onto (N-p, N-q) unitarily: [H, X^N] = 0 for XX+YY, and hamming is invariant under
  complementing both indices. The two blocks therefore have IDENTICAL spectra, bit for bit (measured:
  difference exactly 0.0 at complete N=5,6,7). So "which sector wins" is a SET, never an argmin. A
  strict < on tied floats reports rounding, and under a mathematically equivalent implementation it
  flips: star N=6 and N=7 hand the win to (5,5) and (6,6). Every gate below asserts MEMBERSHIP in the
  winning set and asserts the PARTNER is in it too, turning the degeneracy from a hazard into a law.

  HIGH-Q MECHANISM (degenerate perturbation theory, gamma << J): L = L_H + L_D, L_H dominates.
  The high-Q decay rates = eigenvalues of the L_H-eigenbasis-BLOCK-DIAGONAL part of N_XY (N_XY
  averaged over the ad_H flow; cross-Omega elements drop). So per L_H-eigenspace Omega, the rate
  is 2*gamma * (eigenvalues of P_Omega N_XY P_Omega). Two kinds of slow mode:
    - the BAND EDGE: the (0,1)/(1,0) sector has UNIFORM hamming=1, so N_XY = I there and L_D = -2g*I
      EXACTLY (all Q): rate = 2g exactly, g2 = 1. This floor is always present, so g2 <= 1 always.
    - the CEILING: an Omega=0 (commutant, [H,A]=0) coherence with <n_XY> < 1 dips BELOW the floor.

  (An earlier version projected onto Omega=0 only and missed the band edge -> it over-predicted g2
   for the protecting chain/ring cases. The all-Omega block-diagonal below is the correct model.)

  STAGE 0a (the FLOOR is EXACT, not approximate): on (0,1) the hamming vector is identically 1, an
    INTEGER fact, so it is compared exactly rather than to a tolerance. It follows algebraically that
    L restricted there is L_H - 2*gamma*I and every eigenvalue has Re = -2*gamma at every J. The
    eigen-read is printed against its own error model (eps*||L||/(2*gamma)) as a demonstration, not
    as the gate.

  STAGE 0b (the MECHANISM gate): g2_predicted = min nonzero eigenvalue of the all-Omega
    block-diagonal N_XY, J-independent, compared with the full-L oracle. The truncation is SECOND
    order, O(1/Q^2), not first: ad_H is anti-Hermitian, so the first-order correction to the REAL
    part is purely imaginary and cancels (PROOF_STRUCTURAL_CEILING section 1). The gate is therefore
    an error LAW and not a pin: on rows that carry a truncation at all, the measured order
    log(d(Q1)/d(Q2)) / log(Q2/Q1) is 2.00 +- 0.15, which excludes first order (1.00) by six times the
    band. Rows whose prediction is exact carry no truncation and are gated at the eigensolver floor
    instead. The two populations are separated by EIGHT decades, so the split is not a tuned number.

  STAGE 1 (the PRIZE -- g2(K_N)=4/N derived, not fit): the K_N ceiling lives in the Omega=0
    commutant of the (1,1) single-excitation-coherence sector. GATE: the value = 4/N to MACHINE
    PRECISION for N=5,6,7; (1,1) and its X^N partner are both in the winning set; AND the Schur
    statement itself is gated, not narrated: the minimising eigenspace has multiplicity exactly N-1
    and lies in End(V) to machine precision, with zero weight on the top level.

  STAGE 2 (the N=4 UNIFICATION): GATE: at N=4 the (1,1) commutant = EXACTLY 1.0 (= 4/4 = the band
    edge, no ceiling by itself); the K_4 ceiling comes from the (2,2) HALF-FILLING two-excitation
    sector -- the SAME sector that makes ring-4 special -- with its OWN closed form 2 - 2/sqrt(3)
    ~ 0.845299. The N=4 outlier = 4/N reaching 1.0 so the half-filling mode wins only there.
    ((2,2) at N=4 is its own X^N partner, so this winner carries no tie.)

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
    # vec is ROW-major (flat index a*nb + b, the order `diag` below is built in), so the correct
    # vectorisation of [H,.] is Hp (x) I - I (x) Hq^T.  The transpose is not decoration: it is
    # invisible for XX+YY only because that block is real symmetric, and it would silently build a
    # different object for a complex or non-symmetric block (a DM term, a complex field).
    adH = np.kron(Hp, np.eye(nb)) - np.kron(np.eye(na), Hq.T)   # real symmetric; eigenvalues = Omega
    Omega, V = np.linalg.eigh(adH)
    diag = np.array([bin(A[a] ^ B[b]).count('1') for a in range(na) for b in range(nb)], dtype=float)
    Ntil = V.T @ (diag[:, None] * V)                            # N_XY in the L_H eigenbasis
    oscale = 1.0 + (np.abs(Omega).max() if Omega.size else 0.0)
    tol = TOL_REL * oscale
    # all-Omega block-diagonal part (cross-Omega elements vanish under the ad_H time-average).
    # `same` must be an EQUIVALENCE relation.  If the Omega clusters chain (a < b < c with a~b, b~c,
    # a!~c) then Ntil*same is not block diagonal at all, and its eigvalsh is a THIRD object: neither
    # the block spectra nor their transitive closure.  Nothing in the masking enforces that, so gate
    # it.  Omega is sorted, so the cluster labels reproduce `same` exactly iff it is transitive.
    same = np.abs(Omega[:, None] - Omega[None, :]) < tol
    lab = np.cumsum(np.concatenate(([0], (np.diff(Omega) >= tol).astype(int))))
    assert np.array_equal(same, lab[:, None] == lab[None, :]), (
        f"TRANSITIVITY GATE FIRED at (p,q)=({p},{q}): the Omega clusters chain at tol={tol:.2e}, "
        f"so Ntil*same is not block diagonal and its spectrum is not the union of the block spectra")
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


def global_g2(Hf, N, tie=1e-9):
    """Global high-Q g2 = min nonzero rate over all (p,q) sectors; plus per-sector table.

    Returns (value, winners, table) where `winners` is the SET of sectors within `tie` of the
    minimum.  It has to be a set: X^N maps (p,q) to (N-p,N-q) with a bit-identical spectrum, and on
    the protected topologies the whole band-edge family sits at exactly 1.0 (10 tied sectors at
    chain N=5, spread 1e-15).  An argmin over that reports numpy's rounding rather than the physics,
    and it flips under a mathematically equivalent implementation.  `tie` is not a fitted tolerance:
    the ties are EXACT, and the observed spread inside a tied family is ~1e-15, six decades under."""
    table = {}
    for p in range(N + 1):
        for q in range(N + 1):
            r = sector_analysis(Hf, N, p, q)
            if r is not None:
                table[(p, q)] = r
    vals = {k: v['highq_min'] for k, v in table.items() if v['highq_min'] is not None}
    best = min(vals.values())
    winners = frozenset(k for k, v in vals.items() if v - best < tie)
    return best, winners, table


_SWEEP_CACHE = {}


def sweep(topo, N):
    """global_g2 over the whole (p,q) lattice for one (topology, N), computed ONCE.

    The sweep is the expensive object in this file (81 sectors at N=8, the (4,4) block 4900x4900),
    and the stages below each read the same sweep two or three times, across stages as well as
    within them. Memoising it changes no number; it is the reason this verifier is minutes rather
    than tens of minutes."""
    key = (topo, N)
    if key not in _SWEEP_CACHE:
        _SWEEP_CACHE[key] = global_g2(H_full(topo, N, J=1.0), N)
    return _SWEEP_CACHE[key]


def fmt_win(winners, limit=3):
    """Compact rendering of a winning SET, sorted so it never depends on iteration order."""
    w = sorted(winners)
    head = ",".join(str(x).replace(" ", "") for x in w[:limit])
    return head + ("+%d" % (len(w) - limit) if len(w) > limit else "")


def fullL_g2(topo, N, Q=1000.0):
    """Framework full Liouvillian g2 = strict_gap/(2*gamma) at high Q (the banked oracle)."""
    J = Q * GAMMA
    cs = fw.ChainSystem(N=N, gamma_0=GAMMA, J=J, topology=topo, H_type='xy')
    ev = np.linalg.eigvals(np.asarray(cs.L))
    rates = -ev.real
    dec = rates > 1e-9
    return float(rates[dec].min()) / (2 * GAMMA)


def band_edge_read(topo, N, Q):
    """Eigen-read of the (0,1) block of L at J = Q*gamma, plus ||L||_2.  Returned for demonstration
    only: STAGE 0a proves the value is 1 by construction, so what this measures is the eigensolver."""
    J = Q * GAMMA
    Hf = H_full(topo, N, J)
    A = sector_states(N, 0)
    B = sector_states(N, 1)
    na, nb = len(A), len(B)
    adH = np.kron(Hf[np.ix_(A, A)], np.eye(nb)) - np.kron(np.eye(na), Hf[np.ix_(B, B)].T)
    diag = np.array([bin(A[a] ^ B[b]).count('1') for a in range(na) for b in range(nb)], dtype=float)
    L = -1j * adH - 2 * GAMMA * np.diag(diag)
    return float((-np.linalg.eigvals(L).real).min()) / (2 * GAMMA), float(np.linalg.norm(L, 2))


def se_commutant(A):
    """Omega=0 commutant spectrum of N_XY in the (1,1) sector, from the N x N adjacency ALONE.
    The ceiling value is a single-excitation object, so this never touches the 4^N Liouvillian.
    Returns (eigenvalues ascending, the matching modes as N x N matrices).  The Omega=0 criterion is
    the same relative one sector_analysis uses, so the file carries ONE definition of "Omega = 0"."""
    N = A.shape[0]
    adH = np.kron(A, np.eye(N)) - np.kron(np.eye(N), A.T)
    nxy = np.array([0.0 if a == b else 2.0 for a in range(N) for b in range(N)])
    w, V = np.linalg.eigh(adH)
    ker = V[:, np.abs(w) < TOL_REL * (1.0 + np.abs(w).max())]
    ev, U = np.linalg.eigh(ker.T @ (nxy[:, None] * ker))
    return ev, [(ker @ U[:, i]).reshape(N, N) for i in range(len(ev))]


# =====================================================================================
# STAGE 0a -- the FLOOR is EXACT (an integer fact), so it is compared exactly, not gated
# =====================================================================================
print("=" * 100)
print("STAGE 0a -- THE FLOOR IS EXACT: hamming == 1 identically on (0,1), so L there is L_H - 2*gamma*I")
print("=" * 100)
for N in range(2, 9):
    A = sector_states(N, 0)
    B = sector_states(N, 1)
    ham = np.array([bin(a ^ b).count('1') for a in A for b in B])
    assert np.array_equal(ham, np.ones_like(ham)), \
        f"STAGE 0a GATE FIRED: (0,1) hamming is not identically 1 at N={N}: {sorted(set(ham))}"
print("  N = 2..8: the (0,1) hamming vector is identically 1, compared EXACTLY (integers, no tolerance,")
print("  and topology-independent: the sector is (the vacuum) x (the N single-excitation strings)).")
print("  It follows algebraically that N_XY = I there and L|_(0,1) = L_H - 2*gamma*I, so EVERY")
print("  eigenvalue has Re = -2*gamma at every J: g2 <= 1 for every graph and every Q, and the floor")
print("  cannot move. The eigen-read below is a DEMONSTRATION against its own error model, not a gate.")
print(f"  {'topo':9} {'N':>2} {'Q':>7} {'rate/(2g)':>17} {'|dev|':>10} {'eps*||L||/2g':>14} {'ratio':>7}")
for topo in ('chain', 'complete'):
    for N in (4, 6):
        for Q in (1.0, 1000.0):
            r, nrm = band_edge_read(topo, N, Q)
            floor = np.finfo(float).eps * nrm / (2 * GAMMA)
            print(f"  {topo:9} {N:>2} {Q:>7.0f} {r:>17.14f} {abs(r - 1):>10.2e} {floor:>14.2e} "
                  f"{abs(r - 1) / floor:>7.2f}")

# =====================================================================================
# STAGE 0b -- the MECHANISM gate, written as an ERROR LAW: the truncation is O(1/Q^2)
# =====================================================================================
print()
print("=" * 100)
print("STAGE 0b -- MECHANISM GATE: high-Q ceiling = min nonzero eig of the all-Omega block-diagonal N_XY")
print("  g2_pred (J-independent, degenerate-PT)  vs  the full-L oracle at TWO Q, a decade apart.")
print("  order := log(d(Q_lo)/d(Q_hi)) / log(Q_hi/Q_lo).  Second order => 2.00, first order => 1.00.")
print("  Rows whose prediction is EXACT carry no truncation at all and are gated at the eigensolver")
print("  floor instead; the two populations sit many decades apart (measured below), so the split is not")
print("  a tuned threshold but a reading of a bimodal distribution.")
print("=" * 100)
Q_LO, Q_HI = 100.0, 1000.0
EIG_FLOOR = 1e-9            # measured residuals on the exact rows are ~1e-13, four decades under
print(f"{'topo':9} {'N':>2} {'g2_pred':>12} {'g2_fullL':>12} {'winners':>18} "
      f"{'d(100)':>10} {'d(1000)':>10} {'order':>8}")
rows = []
for topo in ('chain', 'star', 'ring', 'complete'):
    for N in (3, 4, 5, 6):
        gmin, gwin, _ = sweep(topo, N)
        # No try/except.  An oracle that cannot be built is a FAILED gate, not a skipped row: the
        # old form updated `worst` only inside the try, so a crash printed "(skipped)" and passed.
        g2_lo = fullL_g2(topo, N, Q=Q_LO)
        g2_hi = fullL_g2(topo, N, Q=Q_HI)
        d_lo, d_hi = abs(gmin - g2_lo), abs(gmin - g2_hi)
        order = np.log(d_lo / d_hi) / np.log(Q_HI / Q_LO) if d_lo > EIG_FLOOR else None
        rows.append((topo, N, gmin, d_lo, d_hi, order))
        print(f"{topo:9} {N:>2} {gmin:>12.6f} {g2_hi:>12.6f} {fmt_win(gwin):>18} "
              f"{d_lo:>10.2e} {d_hi:>10.2e} {(f'{order:.3f}' if order is not None else 'exact'):>8}")

trunc = [r for r in rows if r[5] is not None]
exact = [r for r in rows if r[5] is None]
for (topo, N, gmin, d_lo, d_hi, order) in trunc:
    assert abs(order - 2.0) < 0.15, (
        f"STAGE 0b GATE FIRED: {topo} N={N} truncation order {order:.3f} is not 2.00 +- 0.15 "
        f"(first order would be 1.00); d({Q_LO:.0f})={d_lo:.2e}, d({Q_HI:.0f})={d_hi:.2e}")
for (topo, N, gmin, d_lo, d_hi, order) in exact:
    assert d_hi < EIG_FLOOR, (
        f"STAGE 0b GATE FIRED: {topo} N={N} predicts an EXACT value, but the oracle differs by "
        f"{d_hi:.2e} at Q={Q_HI:.0f}, above the eigensolver floor {EIG_FLOOR:.0e}")
sep = min(r[3] for r in trunc) / max(r[3] for r in exact)
assert sep > 1e4, f"STAGE 0b GATE FIRED: exact and truncating rows are not separated ({sep:.1e})"
print(f"\nSTAGE 0b PASS: the high-Q ceiling IS the block-diagonal-N_XY minimum. The {len(trunc)} rows that")
print(f"carry a truncation all show order 2.00 +- 0.15 (measured "
      f"{min(r[5] for r in trunc):.3f} .. {max(r[5] for r in trunc):.3f}), so the residual is O(1/Q^2)")
print(f"and NOT O(1/Q): ad_H is anti-Hermitian, so the first-order correction to the REAL part is")
print(f"purely imaginary and cancels. The other {len(exact)} rows predict the exact band edge and carry no")
print(f"truncation at all (worst {max(r[4] for r in exact):.1e} at Q={Q_HI:.0f}); the two populations are")
print(f"separated by {sep:.0e}. Ceiling = an Omega=0 commutant mode < 1; floor = STAGE 0a.")

# =====================================================================================
# STAGE 1 -- the PRIZE: g2(K_N) = 4/N from the (1,1) commutant (S_N standard rep), machine precision
# =====================================================================================
print("\n" + "=" * 100)
print("STAGE 1 -- g2(K_N) = 4/N from the (1,1) commutant (coherences in the (N-1)-fold -J level = S_N standard rep)")
print("=" * 100)
print(f"{'N':>2} {'(1,1) commutant':>16} {'4/N':>10} {'global min':>11} {'winners':>18} {'==4/N exact?':>13}")
for N in (4, 5, 6, 7):
    gmin, gwin, table = sweep('complete', N)
    c11 = table[(1, 1)]['comm_min']
    exact = (c11 is not None) and abs(c11 - 4.0 / N) < 1e-9
    print(f"{N:>2} {c11:>16.9f} {4.0/N:>10.6f} {gmin:>11.9f} {fmt_win(gwin):>18} {('YES' if exact else 'no'):>13}")

for N in (5, 6, 7):
    gmin, gwin, table = sweep('complete', N)
    c11 = table[(1, 1)]['comm_min']
    assert abs(c11 - 4.0 / N) < 1e-9, f"STAGE 1 GATE FIRED: K_{N} (1,1) commutant {c11:.9f} != 4/{N}"
    # MEMBERSHIP, not identity: (1,1) and its X^N partner (N-1,N-1) are exactly degenerate, so an
    # argmin between them is decided by rounding.  Assert both are in the winning set: that turns
    # the degeneracy from a hazard into a gated law.
    assert (1, 1) in gwin, f"STAGE 1 GATE FIRED: K_{N} (1,1) not among the winners {sorted(gwin)}"
    assert (N - 1, N - 1) in gwin, \
        f"STAGE 1 GATE FIRED: K_{N} X^N partner ({N-1},{N-1}) not among the winners {sorted(gwin)}"
    assert abs(gmin - 4.0 / N) < 1e-9, f"STAGE 1 GATE FIRED: K_{N} global g2 {gmin:.9f} != 4/{N}"
# The SCHUR statement itself, gated rather than narrated. The value 4/N alone does not test the rep
# claim: a gate on the number would pass whatever the eigenvector did. What the claim asserts is that
# N_XY is a scalar on each isotypic piece of the multiplicity-free V (x) V*, so the minimising
# eigenspace must have multiplicity EXACTLY N-1 = dim of the S_N standard rep and lie inside End(V),
# V = 1^perp the (N-1)-fold -J level, with zero weight on the top level's own coherence.
# N=8 comes free here: the value is a SINGLE-EXCITATION object, so this needs the N x N adjacency and
# never the 4^N Liouvillian. The global-minimum leg stays at N<=7 (the N=8 sweep is a 4900x4900
# half-filling block); the X^N partner (N-1,N-1) carries the same value by conjugation.
print()
print(f"{'N':>2} {'min nonzero':>14} {'4/N':>10} {'mult':>5} {'dim std':>8} "
      f"{'weight in End(V)':>18} {'weight on top':>14}")
for N in (4, 5, 6, 7, 8):
    A = np.ones((N, N)) - np.eye(N)
    ev, modes = se_commutant(A)
    pos = [i for i, x in enumerate(ev) if x > NZ]
    lam = ev[pos[0]]
    mult = sum(1 for i in pos if abs(ev[i] - lam) < 1e-9)
    Ptop = np.ones((N, N)) / N               # the simple + (N-1)J level
    Pstd = np.eye(N) - Ptop                  # V = the (N-1)-fold -J level (the standard rep)
    wV, wT = [], []
    for i in pos[:mult]:
        m = modes[i]
        n2 = float(np.sum(m * m))
        wV.append(float(np.sum((Pstd @ m @ Pstd) ** 2)) / n2)
        wT.append(float(np.sum((Ptop @ m @ Ptop) ** 2)) / n2)
    print(f"{N:>2} {lam:>14.9f} {4.0/N:>10.6f} {mult:>5} {N-1:>8} {min(wV):>18.15f} {max(wT):>14.2e}")
    assert abs(lam - 4.0 / N) < 1e-9, f"STAGE 1 SCHUR GATE FIRED: K_{N} (1,1) commutant {lam:.9f} != 4/{N}"
    assert mult == N - 1, \
        f"STAGE 1 SCHUR GATE FIRED: K_{N} minimiser multiplicity {mult} != dim of the standard rep {N-1}"
    assert min(wV) > 1 - 1e-12, \
        f"STAGE 1 SCHUR GATE FIRED: K_{N} minimiser leaves End(V), weight {min(wV):.15f}"
    # P m P is (sum_ij m_ij / N^2) * J and the minimiser has zero total sum, so this is EXACTLY zero;
    # what is measured is roundoff squared (~1e-31), and 1e-20 is far above it and far below anything real.
    assert max(wT) < 1e-20, f"STAGE 1 SCHUR GATE FIRED: K_{N} minimiser has top-level weight {max(wT):.2e}"
print(f"\nSTAGE 1 PASS: g2(K_N) = 4/N EXACTLY (machine precision) for N=5,6,7 as the global minimum, with")
print(f"the X^N partner tied to the last bit, and at N=8 in the (1,1) commutant. The rep statement is now")
print(f"gated and not narrated: the minimiser has multiplicity exactly N-1 and lies in End(V) to machine")
print(f"precision. Careful with the wording it replaces: End(V) READ ALONE gives 2/N, because the")
print(f"commutant is End(V) (+) span{{J}}, one dimension more, and that extra dimension is what turns the")
print(f"trivial mode into the exact zero and a 2. The closed form is derived from Schur, no longer a fit.")

# =====================================================================================
# STAGE 2 -- the N=4 UNIFICATION: K_4 ceiling = the (2,2) half-filling sector (ring-4's sector), 2 - 2/sqrt(3)
# =====================================================================================
print("\n" + "=" * 100)
print("STAGE 2 -- N=4 UNIFICATION: K_4 ceiling = the (2,2) half-filling sector (the ring-4 sector)")
print("=" * 100)
K4_CLOSED = 2.0 - 2.0 / sqrt(3.0)        # = 2(1 - 1/sqrt(3)) ~ 0.845299
for topo in ('complete', 'ring'):
    gmin, gwin, table = sweep(topo, 4)
    label = 'K_4' if topo == 'complete' else 'ring-4'
    print(f"\n  {label}: global min nonzero <n_XY> = {gmin:.9f}  in sector(s) {sorted(gwin)}")
    print(f"    {'(p,q)':>7} {'all-Omega min':>15} {'Omega=0 commutant':>18}")
    for (p, q) in sorted(table):
        r = table[(p, q)]
        hm = f"{r['highq_min']:.6f}" if r['highq_min'] is not None else "  -"
        cm = f"{r['comm_min']:.6f}" if r['comm_min'] is not None else "  -"
        print(f"    {str((p, q)):>7} {hm:>15} {cm:>18}")

gK4, winK4, tK4 = sweep('complete', 4)
assert abs(tK4[(1, 1)]['comm_min'] - 1.0) < 1e-9, \
    f"STAGE 2 GATE FIRED: K_4 (1,1) commutant expected 1.0 (=4/4=band edge), got {tK4[(1,1)]['comm_min']:.9f}"
# (2,2) at N=4 is its own X^N partner (N-p = 2), so this winner carries no tie; it is the one place
# in the file where a single sector is the whole winning set, and that is asserted rather than assumed.
assert (2, 2) in winK4, f"STAGE 2 GATE FIRED: K_4 ceiling expected in the (2,2) half-filling sector, got {sorted(winK4)}"
assert winK4 == frozenset({(2, 2)}), \
    f"STAGE 2 GATE FIRED: K_4 ceiling should be the UNIQUE winner (it is X^N self-partnered), got {sorted(winK4)}"
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
# The match is the one gate with content here. The two lines that used to follow it were its own
# consequences: sums_anti is ITSELF palindromic about 0 (max|s + s[::-1]| = 4e-16) and its own max IS
# 2sqrt(2), so once the match passes neither could fire. They are kept as readouts of the MEASURED
# spectrum, and the test that actually discriminates is added: the PERIODIC assignment must NOT fit.
assert np.allclose(e22, sums_anti, atol=1e-9), "STAGE 2b GATE FIRED: (2,2) energies != anti-periodic two-fermion sums"
ek_per = sorted(2.0 * cos(2 * pi * m / 4) for m in range(4))          # the even-string (P) sector
sums_per = np.sort([ek_per[a] + ek_per[b] for a, b in itertools.combinations(range(4), 2)])
assert not np.allclose(e22, sums_per, atol=1e-3), (
    "STAGE 2b GATE FIRED: the PERIODIC assignment fits too, so the match does not identify the sector")
assert np.allclose(e22, -e22[::-1], atol=1e-9), "STAGE 2b GATE FIRED: (2,2) spectrum not palindromic about 0 (chiral K)"
assert abs(e22.max() - two_sqrt2) < 1e-9, f"STAGE 2b GATE FIRED: (2,2) top {e22.max():.6f} != 2sqrt(2)={two_sqrt2:.6f}"
# "2sqrt(2) > 2" would be two literals and could never fire. The statement with content compares the
# MEASURED top against the MEASURED periodic band top, both computed above.
assert e22.max() > sums_per.max() + 1e-9, (
    f"STAGE 2b GATE FIRED: the (2,2) top {e22.max():.6f} does not exceed the periodic band top "
    f"{sums_per.max():.6f}, so the wrap does not raise the edge")
print(f"  periodic 2-fermion sums:        {np.round(sums_per, 6)}  (does NOT fit: the discriminator)")
print(f"  top = {e22.max():.6f} = 2sqrt(2) > periodic band top {sums_per.max():.6f}; +-2sqrt(2) are a CHIRAL K-mirror")
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
# The closed form 4/(N-1) is a statement about the (1,1) COMMUTANT and holds from N=4 on; it is a
# CEILING (i.e. the global minimum, below the floor) only from N=6, because 4/(N-1) < 1 needs N >= 6.
# Both legs are shown, and their difference at N=4,5 is the point rather than an inconvenience.
print(f"{'N':>2} {'(1,1) commutant':>16} {'4/(N-1)':>10} {'global min':>12} {'winners':>18} {'ceiling?':>9}")
for N in (3, 4, 5, 6, 7, 8):
    gmin, gwin, table = sweep('star', N)
    c11 = table[(1, 1)]['comm_min']
    pred = 4.0 / (N - 1)
    tag = 'YES' if gmin < 1.0 - 1e-9 else ('form!=' if abs(c11 - pred) > 1e-9 else 'floor')
    print(f"{N:>2} {c11:>16.9f} {pred:>10.6f} {gmin:>12.9f} {fmt_win(gwin):>18} {tag:>9}")
for N in (4, 5, 6, 7, 8):
    _, _, table = sweep('star', N)
    assert abs(table[(1, 1)]['comm_min'] - 4.0 / (N - 1)) < 1e-9, \
        f"STAGE 3 GATE FIRED: star_{N} (1,1) commutant {table[(1,1)]['comm_min']:.9f} != 4/{N-1}"
# NEGATIVE CONTROL for the form's lower edge, so "N>=4" is gated and not fenced in prose: at N=3 the
# leaf manifold has dimension N-2 = 1 and the second branch of lambda2 takes over, so 4/(N-1) fails.
c11_3 = sweep('star', 3)[2][(1, 1)]['comm_min']
assert abs(c11_3 - 4.0 / 2) > 1e-3, \
    f"STAGE 3 NEGATIVE CONTROL FIRED: star_3 (1,1) commutant {c11_3:.9f} should NOT equal 4/(N-1)=2"
for N in (6, 7, 8):
    gmin, gwin, _ = sweep('star', N)
    # MEMBERSHIP, not identity: (1,1) and its X^N partner (N-1,N-1) are exactly degenerate.  This is
    # the assert the old `gwin == (1,1)` got wrong: under a mathematically equivalent implementation
    # the argmin hands the win to (5,5) at N=6 and (6,6) at N=7, and the gate fired on rounding.
    assert (1, 1) in gwin, f"STAGE 3 GATE FIRED: star N={N} (1,1) not among the winners {sorted(gwin)}"
    assert (N - 1, N - 1) in gwin, \
        f"STAGE 3 GATE FIRED: star N={N} X^N partner ({N-1},{N-1}) not among the winners {sorted(gwin)}"
    assert abs(gmin - 4.0 / (N - 1)) < 1e-9, f"STAGE 3 GATE FIRED: g2(star_{N})={gmin:.9f} != 4/{N-1}"
print(f"\nSTAGE 3 PASS: the (1,1) commutant is 4/(N-1) EXACTLY for N=4..8, and it is the global CEILING")
print(f"from N=6 (4/5, 4/6, 4/7), because 4/(N-1) < 1 needs N >= 6. At N=5 it equals 4/4 = 1.0 and lands ON")
print(f"the band edge (STAGE 4 shows it is pulled below at any finite Q); at N=4 it is 4/3 > 1 and the floor")
print(f"protects. At N=3 the form does not hold at all ({c11_3:.6f}, not 2), gated above as a negative")
print(f"control. The star ceiling is 4/(N-1), NOT a constant 0.80.")

# Honesty check (NOT a universal law): the '4/(m+1)' guess from {complete: m=N-1, star: m=N-2} does
# NOT generalize -- the RING (Fourier-degenerate manifold) breaks it. Print ring (1,1) to show it.
print("\n" + "-" * 100)
print("DIAGNOSTIC (no universal m-law): ring (1,1) commutant vs the tempting 4/(m+1), m=max adjacency degeneracy")
print("-" * 100)
print(f"{'N':>2} {'ring (1,1) commutant':>22} {'4/(m+1) guess':>14} {'fits?':>6} {'closed form':>13} {'>= 1?':>6}")
for N in (4, 5, 6, 7):
    _, _, table = sweep('ring', N)
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
    # the ring's own closed form (ring_ceiling_commutant_sweep.py carries it to N=11): the value MOVES
    # at fixed m = 2, which is what a 4/(m+1) law forbids, and it never goes below the band edge.
    closed = 2 * (N - 2) / N if N % 2 == 0 else 2 * (N - 1) / N
    print(f"{N:>2} {c11s:>22} {guess:>14.6f} {fits:>6} {closed:>13.6f} {str(c11 >= 1 - 1e-9):>6}")
    assert abs(c11 - closed) < 1e-9, f"DIAGNOSTIC GATE FIRED: ring N={N} (1,1) commutant {c11} != {closed}"
    assert c11 >= 1 - 1e-9, f"DIAGNOSTIC GATE FIRED: ring N={N} (1,1) commutant dips below the band edge"
    assert m == 2, f"DIAGNOSTIC GATE FIRED: ring N={N} max adjacency degeneracy {m} != 2"
print("=> ring breaks 4/(m+1) in the sharpest way: m = 2 at every N while the value MOVES, 2(N-2)/N (even)")
print("   / 2(N-1)/N (odd), and it never dips below 1. So the per-family closed forms (complete 4/N, star")
print("   4/(N-1)) are real, and ON THE XY NETWORK the ring carries no ceiling: it stays band-edge-protected")
print("   like the chain. Two fences that the flat sentence used to drop, both gated elsewhere in this file:")
print("   (a) MODEL. This is an XY statement. Under isotropic Heisenberg the SAME 4-cycle ceilings, at the")
print("       K_4 value 2 - 2/sqrt(3) (docstring, PROOF_STRUCTURAL_CEILING section 4, CAUGHT_ERRORS")
print("       2026-08-22); ring N=5 stays at 1 under both models.")
print("   (b) LIMIT. 'No ceiling' is a Q -> infinity statement. At N=4 the (2,2) mode lands exactly ON the")
print("       floor in the limit, and at any FINITE Q it is pulled below it, by 1/(2Q^2) (STAGE 4). So it")
print("       does move the gap; what it does not do is survive the limit. star-5 is the same object with")
print("       coefficient 1 instead of 1/2, and the repo has long written the star case correctly.")

# =====================================================================================
# STAGE 4 -- the MARGINAL cases are ONE object: a closed form equal to exactly 1 sits on the floor
#            only in the limit, and is pulled below it at finite Q as 1 - c/Q^2
# =====================================================================================
print("\n" + "=" * 100)
print("STAGE 4 -- the marginal cases (closed form == 1) are pulled BELOW the floor, as 1 - c/Q^2")
print("  The (0,1) band edge itself never moves (STAGE 0a, exact). What moves is the marginal mode that")
print("  LANDS on it: star-5's (1,1) leaf mode at 4/4, ring-4's (2,2) half-filling mode at 2(N-2)/N.")
print("  Gated as a law: c(Q) -> c_closed with the approach itself second order, order 2.00 +- 0.05,")
print("  read on Q = 10 -> 100. The Q = 1000 column is printed but NOT gated: c multiplies the")
print("  eigensolver floor by Q^2, so at Q = 1000 the deviation and the floor are the same size.")
print("=" * 100)
# c = (1 - g2)*Q^2 amplifies the eigensolver floor by Q^2: at Q = 1000 the deviation |c - c_closed| is
# ~1e-6 and so is the amplified floor, so that point is NOISE-LIMITED and the order read across it comes
# out 1.9, not 2.0. The law is therefore gated on the decade Q = 10 -> 100, where both ends sit far above
# the floor; Q = 1000 is printed with its noise floor beside it, and read, not gated.
QS4 = (10.0, 100.0, 1000.0)
Q_LAW = (10.0, 100.0)
print(f"{'case':10} {'sector':>9} {'c_closed':>10} " + "".join(f"{'c(Q=' + f'{q:.0f})':>16}" for q in QS4)
      + f"{'order':>9} {'noise@1e3':>11}")
for label, topo, N, sector, c_closed in (('star-5', 'star', 5, '(1,1)', 1.0),
                                         ('ring-4', 'ring', 4, '(2,2)', 0.5)):
    cs = [(1.0 - fullL_g2(topo, N, Q=Q)) * Q * Q for Q in QS4]
    dev = [abs(c - c_closed) for c in cs]
    i0, i1 = QS4.index(Q_LAW[0]), QS4.index(Q_LAW[1])
    order = np.log(dev[i0] / dev[i1]) / np.log(Q_LAW[1] / Q_LAW[0])
    noise = np.finfo(float).eps * (QS4[-1] ** 2) / GAMMA      # the floor amplified into c at max Q
    print(f"{label:10} {sector:>9} {c_closed:>10.6f} " + "".join(f"{c:>16.9f}" for c in cs)
          + f"{order:>9.3f} {noise:>11.1e}")
    for Q, c in zip(QS4, cs):
        assert c > 0, (
            f"STAGE 4 GATE FIRED: {label} does NOT sit below the floor at Q={Q:.0f} (c={c:.3e}); the "
            f"'co-occupies rather than undercuts' reading would require c <= 0 at finite Q")
    assert abs(order - 2.0) < 0.05, (
        f"STAGE 4 GATE FIRED: {label} coefficient approaches {c_closed} at order {order:.3f} over "
        f"Q = {Q_LAW[0]:.0f} -> {Q_LAW[1]:.0f}, not 2.00 +- 0.05; c = {cs}")
    assert dev[i1] > 10 * noise, (
        f"STAGE 4 GATE FIRED: {label} deviation {dev[i1]:.2e} at Q={Q_LAW[1]:.0f} is not clear of the "
        f"amplified eigensolver floor {noise:.2e}, so the order above would be reading noise")
    assert dev[-1] < 1e-5, \
        f"STAGE 4 GATE FIRED: {label} c(Q=1000) = {cs[-1]:.9f} has not reached {c_closed}"
print("\nSTAGE 4 PASS: star-5 and ring-4 are ONE object, not two anomalies. Wherever a closed form equals")
print("exactly 1 the mode is marginal, and a marginal mode is pulled under the floor at every finite Q,")
print("by c/Q^2 with c = 1 (star-5) and c = 1/2 (ring-4). Only the (0,1) band edge is genuinely immovable,")
print("and that is an exact integer fact (STAGE 0a), not a limit. Consequence for the wording: 'the ring")
print("co-occupies the floor and does not undercut it' is true of the LIMIT and false at every finite Q.")

print("\nDONE.")
