"""Gate-first: the CRACKED BELL, the ring's wrap-bond crack read in the time domain.

Verifier for experiments/THE_CRACKED_BELL.md (2026-08-31). The object: XY ring,
H = (J/2) Sum (XX+YY), uniform Z-dephasing gamma per site, wrap bond cracked to
J' = J*(1-delta). Convention: the hopping book of PROOF_RING_GAP_DOMINANCE.md,
single-excitation block H_se = J*A (A = ring adjacency), J = 1 throughout. This is
NOT the D10/W1Dispersion book (isotropic Heisenberg, Laplacian) and not a seat/rate
detuning (the blind-seat arc's open item is a different object).

Two blocks of one crack:
  The (0,1) block (vacuum coherence |vac><phi|, STAGES A-C): the dissipator is the
  scalar -2*gamma (Absorption Theorem), so phi(t) = e^(-2*gamma*t) e^(-i*H_se*t) phi(0)
  EXACTLY, and the clock is gamma-free in value. Three first-order laws, derived by
  degenerate perturbation theory and gated against exact eigendecomposition:
    (A) SPLITTING: every m <-> N-m pair splits by the SAME Delta_E = 4*delta*J/N
        (|psi_m(j)|^2 = 1/N: a point defect is flat in mode space); the partners are
        the crack-adapted cos/sin standing waves, antinode respectively node ON the cracked bond,
        shifted by -(2*delta*J/N)*(cos(theta_m) +/- 1) in opposite directions.
    (B) REVERSAL: a launched traveling wave of the PERFECT ring fully reverses at
        T_rev = pi/Delta_E = pi*N/(4*delta*J); this is the walk-time step's discarded
        O(delta) reflection (COUPLING_DEFECT_WALK_TIME_STEP: "not read as signal")
        resonantly accumulated by the closed ring until it is the whole signal.
    (C) VISIBILITY: the reversed amplitude carries the exact envelope
        e^(-2*gamma*T_rev) = exp(-pi*N/(2*Q*delta)), Q = J/gamma -- exact in the
        gamma-factorization, exponent accurate to O(delta) relative. Stage C first
        PROVES the two books from below on the full 4^N Lindbladian of the cracked
        ring at N=5 (the (0,1) block pays exactly -2*gamma, entry-wise; the (1,1)
        block exactly -4*gamma off-diagonal, 0 diagonal), then reads the wall end to
        end with the rate taken OFF the generator (C2, delta = 0.2..0.05, a factor 4).
        The exact floor-set observable of PROOF_RING_GAP_DOMINANCE sees delta = 1e-4;
        there the REVERSAL read is exp-blind (10^-1364 at N=4, Q=20) and the early-time
        AMPLITUDE read of the beat deficit, e^(-2 gamma t)|sin(Delta_E t/2)| peaking at
        t = 1/(2 gamma), is 1.8e-4 and LINEAR in delta: 1.7 decades under a 1% bar. The
        modes the floor set actually loses are elsewhere, and C3 asks THEM instead.
  The (1,1) block (single-excitation density, Haken-Strobl, STAGE D -- the block
  MirrorWorld's Cone runs, run mode `warble N [delta]`, pinned by WarbleTests):
  dephasing is NOT scalar here (off-diagonals pay -4*gamma, the diagonal pays
  nothing), and the same crack's clock is gamma-DRESSED: the circulation's zero
  crossing advances under the watching, and the deepest reversal outlives the naive
  scalar model e^(-4*gamma*t)*R_0(t) at its own best point (the dephasing-free
  diagonal feeds the current back). What decides is not which block is read but whether
  the dissipator on it is a SCALAR, and then only for the reads whose times are zeros:
  the (0,1) block's own beat deficit peaks at t = 1/(2*gamma), a time gamma sets (C2c).
  does to the clock; the fast walk-time front never faced this.
  STAGE E: the SPLIT does not have to be read perturbatively. The whole cracked-ring
  spectrum is the zero set of ONE function of k, exactly and at every delta,
  G(k) = (1-u^2) sin(Nk) cos(k) + [(1+u^2) cos(Nk) - 2u] sin(k), u = J'/J = 1-delta;
  u=1 gives cos(Nk)=1 (the perfect ring) and u=0 gives sin((N+1)k)=0 (the OPEN
  N-chain, F2b), so u is a BOUNDARY-CONDITION parameter and the crack is the whole
  road between modulus N and modulus N+1 (not a road between two topologies: for
  every delta < 1 the graph is still a ring). STAGE A is that equation truncated near
  k_m, a quadratic in x = N*q whose two roots are A's committed branch shifts; what
  the linear carry misses is exactly the dispersion's curvature. This covers form 1
  only: T_rev needs the two-level beat on top, the wall needs the (0,1) scalar
  envelope, and gamma does not appear in G at all.
  The mathematics is standard tight-binding scattering (a rank-2 update of the open
  chain, a one-cell transfer matrix), fenced the way THE_BLIND_SITE fences its own
  rank-1 member. What is ours: the repo did not hold the equation for this object; it
  identifies with the amplitude COUPLING_DEFECT_WALK_TIME_STEP.md already carried
  exactly (1/t is G's two coefficients over 2u, so the ring is that scatterer closed,
  Re[e^(-i N k)/t(k)] = 1); and the m-dependence gets a closed form,
  Split_m = (4 delta J/N)[1 + delta(1/2 - 1/(N sin^2 k_m)) + O(delta^2)],
  which makes STAGE D5's measured pin 0.9705 a value of the law, says the flatness is
  O(N*delta) rather than O(delta), and changes sign at N sin^2(k_m) = 2. STAGE A's own
  reason survives untouched and is gated exactly (E5d): the crack's matrix element
  between the partners is 2*delta*J/N for every m, because |psi_m(j)|^2 = 1/N.

Error model (no-rounding convention, case 2): PT truncates at first order (STAGE E
derives that truncation from the exact law instead of assuming it), so
deviations are O(delta^2) on Delta_E (i.e. O(delta) on Delta_E/delta, so the decade
ratio of e(delta) must sit near 0.1, and THAT is gated) and O(delta) on the clocks;
the two hands of the (1,1) clock differ in POWER, not in coefficient (D1b, D1c):
the PEAK read's envelope is ~0.22*delta and the ZERO CROSSING's is ~0.21*delta^2,
each flat over four bins placed across 1.9 decades of delta, so the crossing beats the peak by a factor
~1/delta and that is why it is the hand to read. There is no delta-independent
floor; the ~2e-3 an earlier reading took for one was a single sample of the peak's
O(delta) envelope sitting in a local minimum of its phase oscillation. The scalar-envelope identity in (C) is case 1 and is asserted as an
entry-wise residual on the full Lindbladian, not imported.
STAGE D's dressed-clock and feedback numbers have no closed form; they are pinned
as committed constants (the WalkTime-plateau genre), and the cross-book gate D1
ties the density-matrix beat to STAGE A's spectral splitting through two
independent routes (superoperator expm vs eigh).
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

J = 1.0
FAIL = []


def gate(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def ring_h(n, delta):
    h = np.zeros((n, n))
    for j in range(n - 1):
        h[j, j + 1] = h[j + 1, j] = J
    h[n - 1, 0] = h[0, n - 1] = J * (1.0 - delta)
    return h


def pair_split(n, delta, m):
    """Exact splitting and shifts of the m-pair, matched by proximity to 2cos(theta_m)."""
    ev = np.linalg.eigvalsh(ring_h(n, delta))
    e0 = 2.0 * J * np.cos(2.0 * np.pi * m / n)
    idx = np.argsort(np.abs(ev - e0))[:2]
    pair = np.sort(ev[idx])
    return pair[1] - pair[0], pair - e0


def plane_wave(n, m):
    j = np.arange(n)
    return np.exp(2j * np.pi * m * j / n) / np.sqrt(n)


print("=" * 78)
print("STAGE A: the splitting law Delta_E = 4*delta*J/N, m-independent at O(delta)")
print("=" * 78)
NS = [6, 7, 8, 12, 13, 16, 32]
DELTAS = [1e-5, 1e-4, 1e-3]


def c_closed(n, m):
    """STAGE E's closed form for the leading m-dependence. Defined ONCE, here, because STAGE A
    needs it to name its own exception and STAGE E needs it four more times: when it was written
    out separately in each place, mutating one copy left the other gates testing the old form."""
    return 0.5 - 1.0 / (n * np.sin(2 * np.pi * m / n)**2)


worst_value = 0.0
for n in NS:
    ms = [m for m in range(1, (n + 1) // 2) if 2 * m != n]
    ratios = []
    for m in ms:
        # e(delta) = |Delta_E/delta - 4/N| must shrink LINEARLY: each decade ratio ~ 0.1.
        e = []
        for d in DELTAS:
            split, _ = pair_split(n, d, m)
            e.append(abs(split / d - 4.0 * J / n))
        floor = 1e-9  # eigensolver noise on Delta_E/delta, worst at the SMALLEST delta (1e-5)
        if all(x_ > floor for x_ in e):
            ratios += [e[0] / e[1], e[1] / e[2]]
        split, _ = pair_split(n, 1e-5, m)
        worst_value = max(worst_value, abs(split / (4.0 * 1e-5 * J / n) - 1.0))
    # the LAW gate: linearity means the decade ratio is 0.1; a quadratic (0.01), a
    # square root (0.32) or a constant (1.0) all fail the band. What is gated is the LEADING
    # CORRECTION being linear, which is a statement about c_m (STAGE E) and not about the header
    # claim: e(delta) = (4J/N)*|c_m|*delta + O(delta^2), so wherever c_m passes through zero the
    # correction is quadratic and this band is FALSE while the header claim holds better than
    # anywhere on this grid. That locus is real and near N = 19.03 at m = 1 (E5c finds it); the
    # grid below keeps |c_m| >= 0.073 and the exception is gated on purpose, just after.
    # no vacuous pass: if the floor filter empties the list nothing was measured, so it is a FAIL
    ok = bool(ratios) and all(0.05 < r < 0.2 for r in ratios)
    det = (f"decade ratios in [{min(ratios):.3f}, {max(ratios):.3f}] over {len(ratios)//2} of "
           f"{len(ms)} m values (the rest fell under the noise floor)"
           if ratios else "NOTHING MEASURED: every m fell under the noise floor")
    gate(f"N={n}: e(delta) shrinks one decade per delta decade (linear law)", ok, det)

gate("the first-order VALUE at delta=1e-5 (a report beside the decade law above, not a law: the "
     "bound is 100x the measurement, and a 0.1% wrong 4J/N would still clear it -- what rejects "
     "that is the decade law, not this row)",
     worst_value < 1e-3, f"worst relative deviation = {worst_value:.2e}")
# and the exception, as a PREDICTION rather than a gap in the grid: at the c_m zero the linear
# term is gone, so the decade ratio must LEAVE the band above and fall toward 0.01. N = 19, m = 1
# sits at c_m = +7.9e-4, and it is the one place where the split law is better than the gate for it.
e19 = [abs(pair_split(19, d_, 1)[0] / d_ - 4.0 * J / 19) for d_ in DELTAS]
r19 = [e19[0] / e19[1], e19[1] / e19[2]]
gate("and where c_m vanishes the decade law must FAIL, which is the closed form showing up in "
     "STAGE A: N=19, m=1 leaves the band from below",
     min(r19) < 0.05 and abs(c_closed(19, 1)) < 1e-3,
     f"c_m = {c_closed(19, 1):+.2e}, decade ratios {r19[0]:.3f} and {r19[1]:.3f} against the "
     f"linear band [0.05, 0.2]; the grid above keeps |c_m| >= "
     f"{min(abs(c_closed(n_, 1)) for n_ in NS):.3f}")
print()
print("branch shifts (N=12): even-across-the-crack partner by -(2d/N)(cos+1), odd by -(2d/N)(cos-1)")
n, d = 12, 1e-4
for m in [1, 2, 5]:
    _, shifts = pair_split(n, d, m)
    th = 2.0 * np.pi * m / n
    pred = np.sort(np.array([-(2 * d * J / n) * (np.cos(th) + 1.0),
                             -(2 * d * J / n) * (np.cos(th) - 1.0)]))
    rel = np.max(np.abs(shifts - pred)) / (4 * d * J / n)
    rel2 = np.max(np.abs(pair_split(n, 10 * d, m)[1] - 10 * pred)) / (4 * 10 * d * J / n)
    # the shifts are first order too, so the residual is O(delta) and the LAW is the decade ratio;
    # a bare threshold here would be a number that merely passes (5e-3 against 2.8e-5).
    gate(f"N=12 m={m} branch shifts match PT, and the residual is O(delta): tenfold delta, "
         f"tenfold residual", 5.0 < rel2 / rel < 20.0,
         f"rel={rel:.2e} at delta={d:g} and {rel2:.2e} at {10*d:g}, ratio {rel2/rel:.1f}")
    # the ASSIGNMENT, not just the sorted values: the lower branch is ALWAYS the
    # even-across-the-crack partner (equal signs at the bond's two ends, <B> > 0,
    # loses energy when the bond weakens); the upper is odd (a sign change across
    # the crack, <B> < 0). Which partner feels the bond MORE flips at m = N/4
    # (2cos^2 vs 2sin^2 of theta/2), so parity, not weight, is the invariant.
    ev_, vec_ = np.linalg.eigh(ring_h(n, d))
    e0_ = 2.0 * J * np.cos(th)
    idx_ = np.argsort(np.abs(ev_ - e0_))[:2]
    lo_, hi_ = idx_[np.argsort(ev_[idx_])]
    w_lo = vec_[0, lo_] * vec_[n - 1, lo_]
    w_hi = vec_[0, hi_] * vec_[n - 1, hi_]
    gate(f"N=12 m={m} lower branch is the even-across-the-crack partner",
         w_lo > 0 > w_hi,
         f"bond-end product lower {w_lo:+.4f} vs upper {w_hi:+.4f}")

print()
print("=" * 78)
print("STAGE B: the reversal  T_rev = pi/Delta_E = pi*N/(4*delta*J) (no gamma here at all; that")
print("         it is gamma-free is C1b's scalar dissipator, and C2 reads the wall end to end)")
print("=" * 78)
# The peak read of T_rev, and its error model is the same shape as D1c's: a ONE-SIDED bound
# |dev| <= C*delta is satisfied by every power above one and so certifies nothing about the law it
# names. It is also false here as a bound: the earlier C = 0.5 is violated at 4 of 191 sampled
# delta in [0.001, 0.02] (worst 0.592*delta at delta = 0.0018), and passed only because the three
# sampled points missed them. What IS the law is the flatness of the bin maximum, so that is gated,
# over eight CONTIGUOUS octaves spanning 2.4 decades (a delta^2 law would move the ratio by ~256).
# The rows below are printed as readings, not gated: at the smallest delta the argmax is still the
# nearest grid point rather than the read, and a grid quantum is not a measurement.
def trev_dev(n_, m_, d_, pts=40001):
    """|T_rev*Delta_E/pi - 1| off the pair's exact split, and the peak overlap it is read at."""
    ev_, vec_ = np.linalg.eigh(ring_h(n_, d_))
    c0_ = vec_.conj().T @ plane_wave(n_, m_)
    sp_ = pair_split(n_, d_, m_)[0]
    ts_ = np.linspace(0.0, 1.4 * np.pi / sp_, pts)
    p_ = np.abs(plane_wave(n_, -m_).conj()
                @ (vec_ @ (np.exp(-1j * np.outer(ev_, ts_)) * c0_[:, None]))) ** 2
    return abs(ts_[int(np.argmax(p_))] * sp_ / np.pi - 1.0), float(p_.max())


OCTAVES = [(0.0005 * 2 ** i_, 0.001 * 2 ** i_) for i_ in range(8)]
for (n, m) in [(12, 1), (8, 3), (16, 2)]:
    rows = [(d_, ) + trev_dev(n, m, d_) for d_ in (1e-3, 3e-3, 1e-2)]
    print(f"  reading  N={n} m={m}: " + ", ".join(f"{d_:g}:dev={dv_:.2e}={dv_/d_:.3f}*delta, "
                                                  f"P_-={pk_:.6f}" for d_, dv_, pk_ in rows))
    env = [max(trev_dev(n, m, x_)[0] / x_ for x_ in np.linspace(lo_, hi_, 9)) for lo_, hi_ in OCTAVES]
    # the boolean is invariant under any O(delta) change of reference, so it gates the POWER and
    # not the choice between pi/Delta_E_exact and the first-order pi*N/(4*delta*J): what fixes that
    # choice is E6a (the curve meets eigvalsh) and E5a (the c_m correction), not this row.
    gate(f"B N={n} m={m}: T_rev goes like pi/Delta_E with an O(delta) envelope -- the POWER gated as "
         f"a flat bin maximum over 2.4 decades, the reference fixed elsewhere",
         max(env) / min(env) < 2.5,
         f"|dev|/delta per octave " + ", ".join(f"{lo_:g}:{c_:.3f}" for (lo_, _), c_ in zip(OCTAVES, env))
         + f" -- coefficient {min(env):.2f}..{max(env):.2f}, ratio {max(env)/min(env):.2f} where a "
           f"delta^2 law would give ~256")

# and the completeness, which the document's Fences quote: the peak overlap itself falls with delta
# (admixture into OTHER pairs, which no law about the split speaks to). The displayed quantity is a
# PROBABILITY, so its deficit is the square of an O(delta) amplitude admixture and the power to gate
# is 2, not 1. Measured that way it is a law and not a pin: 1-P over delta^2 is flat across five
# octaves for all three cases, while 1-P over delta moves by 11x to 28x over the same span.
def completeness(n_, m_, delta_):
    ev_, vec_ = np.linalg.eigh(ring_h(n_, delta_))
    c_ = vec_.conj().T @ plane_wave(n_, m_)
    b_ = vec_.conj().T @ plane_wave(n_, -m_)
    sp_ = pair_split(n_, delta_, m_)[0]
    ts_ = np.linspace(0.0, 1.4 * np.pi / sp_, 20001)
    return float((np.abs(np.conj(b_) @ (np.exp(-1j * np.outer(ev_, ts_)) * c_[:, None])) ** 2).max())


COCT = [(0.005 * 2 ** i_, 0.01 * 2 ** i_) for i_ in range(5)]
for (nn_, mm_) in [(12, 1), (8, 3), (16, 2)]:
    sq = [max((1 - completeness(nn_, mm_, x_)) / x_**2 for x_ in np.linspace(lo_, hi_, 9))
          for lo_, hi_ in COCT]
    li = [max((1 - completeness(nn_, mm_, x_)) / x_ for x_ in np.linspace(lo_, hi_, 9))
          for lo_, hi_ in COCT]
    gate(f"B-completeness N={nn_} m={mm_}: the reversal is full only to O(delta^2) in probability, "
         f"the POWER gated as a flat bin maximum over 1.5 decades",
         max(sq) / min(sq) < 2.0 and max(li) / min(li) > 5.0,
         "(1-P)/delta^2 per octave " + ", ".join(f"{lo_:g}:{c_:.3f}" for (lo_, _), c_ in zip(COCT, sq))
         + f" (ratio {max(sq)/min(sq):.2f}), against (1-P)/delta moving by {max(li)/min(li):.1f}x "
           f"over the same span; at delta=0.3 the peak overlap is "
           f"{completeness(nn_, mm_, 0.3):.6f}")

print()
print("=" * 78)
print("STAGE C: the two books from below (full Lindbladian, N=5), then the visibility wall")
print("=" * 78)


def full_liouvillian(n_, delta_, gamma_):
    """The full 4^n Lindblad superoperator of the cracked XY ring, row-major vec."""
    dim = 2 ** n_
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.diag([1.0, -1.0]).astype(complex)

    def site(op, l_):
        mres = np.array([[1.0]], dtype=complex)
        for k_ in range(n_):
            mres = np.kron(mres, op if k_ == l_ else np.eye(2))
        return mres

    h_ = np.zeros((dim, dim), dtype=complex)
    for (a_, b_) in [(a_, a_ + 1) for a_ in range(n_ - 1)] + [(n_ - 1, 0)]:
        w_ = J * (1 - delta_) if (a_, b_) == (n_ - 1, 0) else J
        h_ += (w_ / 2) * (site(sx, a_) @ site(sx, b_) + site(sy, a_) @ site(sy, b_))
    eye_ = np.eye(dim)
    L_ = -1j * (np.kron(h_, eye_) - np.kron(eye_, h_.T))
    D_ = np.zeros((dim * dim, dim * dim), dtype=complex)
    for l_ in range(n_):
        z_ = site(sz, l_)
        D_ += gamma_ * (np.kron(z_, z_.T) - np.kron(eye_, eye_))
    return h_, L_ + D_, D_


# C1: both books, entry-wise. The closed forms rest on the (0,1) block paying exactly
# the scalar -2*gamma (Absorption Theorem) and the Warble's (1,1) block paying exactly
# -4*gamma off-diagonal / 0 diagonal (Haken-Strobl); a wrong rate book (-gamma, -3*gamma)
# passes every spectral gate above, so it is asserted here on the full generator.
nC, dC, gC = 5, 0.15, 0.05
hF, LF, DF = full_liouvillian(nC, dC, gC)
dimC = 2 ** nC
site_idx = [1 << (nC - 1 - j) for j in range(nC)]     # e_j basis index (site 0 leftmost)
hse = np.array([[hF[site_idx[a], site_idx[b]].real for b in range(nC)] for a in range(nC)])
gate("C1a the single-excitation block of the Pauli H IS J*A_cracked, exactly",
     np.max(np.abs(hse - ring_h(nC, dC))) == 0.0,
     f"max|diff| = {np.max(np.abs(hse - ring_h(nC, dC))):.1e}")
# The exact route first: the (0,1) block read entry-wise out of the full Lindbladian IS
# i*h_se^T - 2*gamma*I, with no eigensolver anywhere, so it is compared exactly and not gated.
rows01 = [0 * dimC + site_idx[jj] for jj in range(nC)]
block01 = LF[np.ix_(rows01, rows01)]
target01 = 1j * hse.T - 2.0 * gC * np.eye(nC)
gate("C1b the (0,1) block IS i*h_se^T - 2*gamma*I, entry-wise and exactly",
     np.max(np.abs(block01 - target01)) == 0.0,
     f"max|diff| = {np.max(np.abs(block01 - target01)):.1e}")
evC, vecC = np.linalg.eigh(ring_h(nC, dC))
worst_mode = 0.0
for k in range(nC):
    O = np.zeros((dimC, dimC), dtype=complex)
    for jj in range(nC):
        O[0, site_idx[jj]] = np.conj(vecC[jj, k])     # |vac><phi_k|
    v = O.reshape(-1)
    res = LF @ v - (1j * evC[k] - 2.0 * gC) * v
    worst_mode = max(worst_mode, np.max(np.abs(res)))
gate("C1c and every (0,1) eigenmode pays -2*gamma + i*E_k -- read off the FULL 4^N generator, so this also certifies the block is invariant, which C1b cannot see (eigensolver route, so a residual)",
     worst_mode < 1e-14, f"worst entry-wise residual = {worst_mode:.2e}")
worst_hs = 0.0
for (a, b) in [(0, 1), (1, 4), (2, 3)]:
    O = np.zeros((dimC, dimC), dtype=complex)
    O[site_idx[a], site_idx[b]] = 1.0
    v = O.reshape(-1)
    worst_hs = max(worst_hs, np.max(np.abs(DF @ v + 4.0 * gC * v)))
for a in [0, 3]:
    O = np.zeros((dimC, dimC), dtype=complex)
    O[site_idx[a], site_idx[a]] = 1.0
    worst_hs = max(worst_hs, np.max(np.abs(DF @ O.reshape(-1))))
gate("C1d the (1,1) block pays exactly -4*gamma off-diagonal, 0 diagonal (Haken-Strobl)",
     worst_hs == 0.0, f"worst residual = {worst_hs:.1e}")

# C1e: the circulation is not an observable of the (0,1) block AT ALL. The document's headline
# rests on this and it had no verifier: on |vac><phi| the read sum_a Im rho[a,a+1] touches (1,1)
# entries only, so it is identically zero there, while the same read on |phi><phi| is the ring
# current. An exact route (a structural zero), so an exact comparison.
vacC = np.zeros(2 ** nC, dtype=complex); vacC[0] = 1.0
seC = [np.zeros(2 ** nC, dtype=complex) for _ in range(nC)]
for a_ in range(nC):
    seC[a_][(1 << (nC - 1 - a_))] = 1.0
phiC = sum(plane_wave(nC, 1)[a_] * seC[a_] for a_ in range(nC))
rho01 = np.outer(vacC, phiC.conj())
rho11 = np.outer(phiC, phiC.conj())
circ01 = sum(seC[a_].conj() @ rho01 @ seC[(a_ + 1) % nC] for a_ in range(nC)).imag
circ11 = sum(seC[a_].conj() @ rho11 @ seC[(a_ + 1) % nC] for a_ in range(nC)).imag
gate("C1e the circulation is a (1,1) observable: on the (0,1) coherence it is a STRUCTURAL zero, "
     "while the same read on the single-excitation density is the ring current",
     circ01 == 0.0 and abs(circ11 + np.sin(2 * np.pi / nC)) < 1e-14,
     f"sum_a Im rho[a,a+1] on |vac><phi| is {circ01!r}; on |phi><phi| it is {circ11:.6f}, which is "
     f"-sin(2 pi m/N) = {-np.sin(2 * np.pi / nC):.6f} at N={nC}, m=1")

# C2: the visibility wall. Three versions of this gate failed in three different ways and the third
# is the instructive one, so all three are recorded. (1) It computed |t* - t_first|/t_first and
# called that the visibility exponent: Stage A's split error relabelled, with gamma absent from the
# boolean. (2) The rewrite put gamma on BOTH sides, where it cancelled; and scaling gamma with delta
# to keep the amplitude representable made gamma/delta constant, so a (gamma/delta)^2 law passed
# bit-identically. (3) Fixing gamma fixed that, but Q was then exercised at ONE value, and the
# residual -- normalised by an exponent proportional to gamma -- went red at small gamma where the
# law is equally true, because the gamma-free completeness term sat in its numerator.
# What the wall actually IS, is a composition: the (0,1) block pays exactly -2*gamma (C1b), the beat
# takes T_rev = pi/Delta_E (STAGE B), and Delta_E is 4*delta*J/N to first order (STAGE A). So it is
# gated as a composition, in two halves that fail to different mutations, and NOTHING is typed on
# both sides: the amplitude is PROPAGATED on the full 4^N generator.
def wall_amplitude(n_, delta_, gamma_):
    """Propagate |vac><phi| on the full generator to T_rev and read the counter-wave amplitude."""
    dim_ = 2 ** n_
    vac_ = np.zeros(dim_, dtype=complex); vac_[0] = 1.0
    se_ = [np.zeros(dim_, dtype=complex) for _ in range(n_)]
    for a_ in range(n_):
        se_[a_][1 << (n_ - 1 - a_)] = 1.0
    phi_ = sum(plane_wave(n_, 1)[a_] * se_[a_] for a_ in range(n_))
    chi_ = sum(plane_wave(n_, -1)[a_] * se_[a_] for a_ in range(n_))
    t_star = np.pi / pair_split(n_, delta_, 1)[0]
    rt_ = (expm(full_liouvillian(n_, delta_, gamma_)[1] * t_star)
           @ np.outer(vac_, phi_.conj()).reshape(-1)).reshape(dim_, dim_)
    ev_, vec_ = np.linalg.eigh(ring_h(n_, delta_))
    c_ = vec_.conj().T @ plane_wave(n_, 1)
    ov_ = abs(plane_wave(n_, -1).conj() @ (vec_ @ (np.exp(-1j * ev_ * t_star) * c_)))
    return abs(vac_.conj() @ rt_ @ chi_), np.exp(-2 * gamma_ * t_star) * ov_, t_star


WALL_PTS = [(4, 0.2, 0.02), (4, 0.2, 0.002), (4, 0.1, 0.02), (5, 0.2, 0.02)]
wall_rows = []
for (en_, dd_, gg_) in WALL_PTS:
    amp_, fac_, t_ = wall_amplitude(en_, dd_, gg_)
    wall_rows.append((en_, dd_, gg_, amp_, fac_, abs(amp_ / fac_ - 1),
                      abs(t_ / (np.pi * en_ / (4 * dd_ * J)) - 1)))
gate("C2a the wall's first half, from below: the PROPAGATED (0,1) amplitude at T_rev is exactly "
     "e^(-2*gamma*T_rev) times the unitary overlap, so the envelope really is the scalar (this is "
     "where a wrong rate book dies: 1.5*gamma moves the residual from 1e-15 to 0.26)",
     max(r_[5] for r_ in wall_rows) < 1e-12,
     ", ".join(f"N={n_} d={d_:g} g={g_:g}: |amp| {a_:.6e} vs {f_:.6e}, rel {r_:.1e}"
               for n_, d_, g_, a_, f_, r_, _ in wall_rows))
exp_err = {(n_, d_): e_ for n_, d_, g_, _, _, _, e_ in wall_rows}
# The exponent's residual is |T_rev/T_first - 1|, in which gamma does not appear at all -- so its
# gamma-independence is STRUCTURAL and comparing it across two gamma would be comparing a number to
# itself (an earlier version did exactly that and printed "the two Q agree to 0.0e+00"). What is
# left to gate is the delta law, and that is all this gate claims. The wall's gamma content is C2a's.
gate("C2b and the second half, the exponent: 2*gamma*T_rev against pi*N/(2*Q*delta) differs by "
     "|T_rev/T_first - 1| alone, which gamma does not enter, so what is measurable in it is the "
     "delta law -- it is STAGE A's split error carried into the exponent and is named as that",
     0.4 < exp_err[(4, 0.1)] / exp_err[(4, 0.2)] < 0.65,
     ", ".join(f"N={n_} d={d_:g} g={g_:g}: {e_:.4e}" for n_, d_, g_, _, _, _, e_ in wall_rows)
     + f"; halving delta halves it ({exp_err[(4, 0.1)] / exp_err[(4, 0.2)]:.3f}). The rows at "
       f"gamma={WALL_PTS[0][2]:g} and {WALL_PTS[1][2]:g} are bit-identical BY CONSTRUCTION, not by "
       f"measurement, and are printed to show that")

# C2c: the beat deficit, measured on the TRUE (0,1) evolution rather than on the model formula.
# The earlier version evaluated e^(-2 gamma t)|sin(Delta_E t/2)| and compared it to constants read
# off that same formula, so three of its four conjuncts were arithmetic about a typed expression.
def beat_deficit_true(n_, delta_, gamma_, tmax=400.0, pts=800001):
    """max_t of the true (0,1) counter-wave amplitude: |<psi_-m|phi(t)>| with the scalar envelope."""
    ev_, vec_ = np.linalg.eigh(ring_h(n_, delta_))
    c_ = vec_.conj().T @ plane_wave(n_, 1)
    ts_ = np.linspace(0.0, tmax, pts)
    a_ = np.abs(plane_wave(n_, -1).conj()
                @ (vec_ @ (np.exp(-1j * np.outer(ev_, ts_)) * c_[:, None]))) * np.exp(-2 * gamma_ * ts_)
    return float(a_.max()), float(ts_[int(np.argmax(a_))])


bd, bt = beat_deficit_true(4, 1e-4, 0.05)
bd2, bt2 = beat_deficit_true(4, 1e-4, 0.02)
bd10, _ = beat_deficit_true(4, 1e-3, 0.05)
model = float((np.exp(-2 * 0.05 * np.linspace(0, 400, 800001))
               * np.abs(np.sin((4 * 1e-4 * J / 4) * np.linspace(0, 400, 800001) / 2))).max())
d_beat = brentq(lambda x_: beat_deficit_true(4, x_, 0.05)[0] - 0.01, 1e-4, 0.05)
d_rev = np.pi * 4 / (2.0 * (J / 0.05) * np.log(1 / 0.01))   # N and Q read, not typed
gate("C2c the beat deficit read on the TRUE (0,1) evolution at N=4, Q=20: 1.8e-4 at delta=1e-4, "
     "peaking at t = 1/(2*gamma), linear in delta, and the closed form e^(-2 g t)|sin(dE t/2)| "
     "is right to O(delta)",
     abs(bt - 1 / (2 * 0.05)) < 0.05 and abs(bt2 - 1 / (2 * 0.02)) < 0.05
     and 9.0 < bd10 / bd < 11.0
     and abs(bd / model - 1) < 1e-4 and abs(d_beat - 0.0054) < 5e-4,
     f"true peak {bd:.6e} at t={bt:.3f} (1/(2*gamma) = {1/(2*0.05):.1f}), and at gamma=0.02 the "
     f"peak moves to t={bt2:.3f} against 1/(2*gamma) = {1/(2*0.02):.1f}, so the clock really is "
     f"gamma's; against the closed form's "
     f"{model:.6e}, {abs(bd/model-1):.1e} relative; tenfold delta gives {bd10:.4e}, ratio "
     f"{bd10/bd:.2f}; it clears a 1% bar at delta = {d_beat:.5f} while the reversal read clears "
     f"only at delta = {d_rev:.5f} (that one is the closed form's own root, {np.log10(d_beat/1e-4):.2f} "
     f"and {np.log10(d_rev/1e-4):.2f} decades above the crack the floor set sees)")

print()
print("the visibility wall (reversed amplitude = 10^-x), Q = J/gamma:")
print(f"{'N':>4} {'delta':>8} {'Q':>6} {'x = pi*N/(2*Q*delta)/ln10':>28}")
for (n_, d_, q_) in [(4, 1e-4, 20), (4, 1e-3, 20), (4, 0.068, 20),
                     (12, 1e-2, 20), (12, 1e-2, 200), (32, 1e-3, 200)]:
    x = np.pi * n_ / (2.0 * q_ * d_) / np.log(10.0)
    print(f"{n_:>4} {d_:>8.0e} {q_:>6} {x:>28.1f}")

# C3: the mode that ACTUALLY leaves the floor set. The two reads above both live on the (0,1)
# block, and the floor-set exception the proof's Scope reports is not there: at N=4 the modes that
# leave carry |Im| = 2.8267, while the (0,1) block's |Im| are only {0, 2}. So the question "is the
# floor set's hypersensitivity physical" has to be asked of the departing mode itself, and here it
# is: its rate walks off -2*gamma QUADRATICALLY in delta, three decade ratios of 100, so at the
# delta = 1e-4 crack that empties the set the rate has moved by four parts per million of gamma.
# The claim is the power, gated as decade ratios; the delta = 0 control says the metric is zero
# where it must be. This is a DIFFERENT observable from the beat deficit above, on different modes
# in a different sector, and the two are not compared here beyond noting where each clears 1%.
# The frequency the departing pair sits at is COMPUTED, not typed: it is the anti-periodic
# half-filling sum the proof names, 2*sqrt(2)*J at N=4, and the window is relative so another J
# does not empty it (the earlier literal 2.8266588 silently encoded J = 1 and crashed at J = 2).
DEPART_IM = 2.0 * np.sqrt(2.0) * J


def departing_shift(delta_):
    """|Re + 2*gamma| of the floor mode at |Im| = 2*sqrt(2)*J that the crack pushes off the floor."""
    w_ = np.linalg.eigvals(full_liouvillian(4, delta_, 0.05)[1])
    sel_ = np.abs(np.abs(w_.imag) - DEPART_IM) < 0.02 * DEPART_IM
    if not sel_.any():
        raise AssertionError(f"no mode near |Im| = {DEPART_IM} at delta={delta_}")
    return float(np.abs(w_.real[sel_] + 2 * 0.05).min())


def departing_count(delta_, tol=1e-8):
    """How many modes above the (0,1) block's own |Im| are ON the floor: the proof's exception."""
    w_ = np.linalg.eigvals(full_liouvillian(4, delta_, 0.05)[1])
    return int(((np.abs(w_.real + 2 * 0.05) < tol) & (np.abs(w_.imag) > 0.9 * DEPART_IM)).sum())


shifts = [departing_shift(dd_) for dd_ in (1e-5, 1e-4, 1e-3)]
shift_ratios = [shifts[i_ + 1] / shifts[i_] for i_ in range(2)]
c_quad = shifts[2] / 1e-6
gate("C3 the floor set's exception is hypersensitive but its MODE is not: the departing mode's rate "
     "leaves -2*gamma as O(delta^2), gated as decade ratios of 100",
     all(90.0 < r_ < 110.0 for r_ in shift_ratios) and departing_shift(0.0) < 1e-12
     and all(departing_count(0.0, t_) == 2 for t_ in (1e-9, 1e-8, 1e-7))
     and all(departing_count(1e-4, t_) == 0 for t_ in (1e-9, 1e-8, 1e-7))
     and departing_count(1e-4, 1.5 * shifts[1]) == 2,   # flips back just above the MEASURED shift
     f"|Re+2g| = " + ", ".join(f"{dd_:g}:{s_:.3e}" for dd_, s_ in zip((1e-5, 1e-4, 1e-3), shifts))
     + f", ratios {shift_ratios[0]:.1f} and {shift_ratios[1]:.1f}; delta=0 control "
       f"{departing_shift(0.0):.1e}. The set really is the proof's exception and not a window: "
       f"above the (0,1) block's own |Im| there are exactly {departing_count(0.0)} floor modes at "
       f"delta=0 and {departing_count(1e-4)} at delta=1e-4, at every floor tolerance from 1e-9 to "
       f"1e-7, and flipping back to 2 at {1.5 * shifts[1]:.1e}, one and a half times the measured "
       f"shift, which is where it "
       f"must -- so the threshold is being read rather than carrying the verdict. At the delta=1e-4 crack that empties the set the shift is "
       f"{shifts[1] / 0.05 * 100:.4f}% of gamma, and it reaches 1% of gamma only at "
       f"delta = {np.sqrt(0.01 * 0.05 / c_quad):.4f}")

print()
print("=" * 78)
print("STAGE D: the (1,1) block (Haken-Strobl, the Cone's block): the gamma-DRESSED clock")
print("=" * 78)


def liou11(n, delta, gamma):
    h = ring_h(n, delta)
    eye = np.eye(n)
    lh = -1j * (np.kron(h, eye) - np.kron(eye, h.T))
    dd = np.diag([-4.0 * gamma if a != b else 0.0 for a in range(n) for b in range(n)])
    return lh + dd


# and the seam: every STAGE D reading propagates liou11, an independently written N^2 x N^2
# matrix, while C1d certifies the rate book on the FULL 4^N generator. Nothing tied the two, so a
# wrong rate book in liou11 alone left C1d and D2's law gate green. Tie them entry-wise, once.
def liou11_from_full(n_, delta_, gamma_):
    """The (1,1) block cut out of the full generator, for comparison with the hand-written liou11."""
    dim_ = 2 ** n_
    idx_ = [1 << (n_ - 1 - a_) for a_ in range(n_)]   # single-excitation basis states, popcount 1
    # (the complementary indices dim-1-(1<<k) give the (N-1)-excitation block, which agrees here by
    #  particle-hole symmetry and does NOT under a site-dependent field: a different object.)
    lf_ = full_liouvillian(n_, delta_, gamma_)[1].reshape(dim_, dim_, dim_, dim_)
    return np.array([[lf_[a_, b_, c_, d_] for c_ in idx_ for d_ in idx_]
                     for a_ in idx_ for b_ in idx_])


seam = np.max(np.abs(liou11(5, 0.15, 0.05) - liou11_from_full(5, 0.15, 0.05)))
gate("D0 the propagator STAGE D runs is the (1,1) block of the generator C1 certifies, entry by "
     "entry (an exact route: two constructions of one matrix)",
     seam == 0.0, f"max |liou11 - full_liouvillian's (1,1) block| = {seam!r} at N=5, delta=0.15, "
                  f"gamma=0.05")


def circulation_series(n, m, delta, gamma, tmax, steps):
    psi = plane_wave(n, m)
    rho = np.outer(psi, psi.conj()).reshape(-1)
    P = expm(liou11(n, delta, gamma) * (tmax / steps))
    out = []
    for _ in range(steps + 1):
        r = rho.reshape(n, n)
        out.append(sum(r[a, (a + 1) % n] for a in range(n)).imag)
        rho = P @ rho
    return np.array(out)


def zero_crossing(r, ts):
    for k in range(1, len(ts)):
        if r[k - 1] > 0 >= r[k]:
            return ts[k - 1] + (ts[k] - ts[k - 1]) * r[k - 1] / (r[k - 1] - r[k])
    return np.nan


n, m, d, tmax, steps = 8, 1, 0.15, 58.0, 5800
ts = np.linspace(0, tmax, steps + 1)
r0 = circulation_series(n, m, d, 0.0, tmax, steps)
r0 /= r0[0]

# D1 cross-book: the (1,1) gamma=0 clock against STAGE A's spectral splitting, on both hands.
# The read is shifted by the fast admixture inside an ENVELOPE whose phase oscillates with delta
# (T_zero ~ 1/delta wraps the fast phases), so no single sample and no two-point ratio can see the
# law: each bin is read as a BIN MAXIMUM. Measured, the two hands differ in POWER:
#   the reversal PEAK   |dev| ~ 0.22 * delta
#   the ZERO CROSSING   |dev| ~ 0.21 * delta^2
# each flat over four bins placed across 1.9 decades, so the crossing beats the peak by ~1/delta.
# A one-sided bound of the
# form |dev| <= C*delta cannot tell those apart, which is why both are gated as flat RATIOS.


def zero_dev_pure(n_, m_, delta_):
    """Signed dev of the zero-crossing clock vs pi/(2*DeltaE), pure state, vectorized."""
    h_ = ring_h(n_, delta_)
    ev_, vec_ = np.linalg.eigh(h_)
    psi_ = plane_wave(n_, m_)
    c_ = vec_.conj().T @ psi_
    sp_ = pair_split(n_, delta_, m_)[0]
    ts_ = np.linspace(0.5 * np.pi / (2 * sp_), 1.2 * np.pi / (2 * sp_), 20001)
    F = vec_ @ (np.exp(-1j * np.outer(ev_, ts_)) * c_[:, None])
    cur = np.sum((F * np.roll(F, -1, axis=0).conj()).imag, axis=0)
    r_ = cur / np.sum((psi_ * np.roll(psi_, -1).conj()).imag)
    for k in range(1, len(ts_)):
        if r_[k - 1] > 0 >= r_[k]:
            tz_ = ts_[k - 1] + (ts_[k] - ts_[k - 1]) * r_[k - 1] / (r_[k - 1] - r_[k])
            return tz_ * 2 * sp_ / np.pi - 1.0
    return np.nan


tz_expm = zero_crossing(r0, ts)
sp15 = pair_split(n, d, m)[0]
dev_expm = tz_expm * 2 * sp15 / np.pi - 1.0
dev_pure = zero_dev_pure(n, m, d)
# Two float paths to one real number (a matrix exponential of the superoperator against a
# pure-state reader), so no exact route. What separates them is not float noise (~1e-12) but their
# different time grids and crossing interpolations, ~1e-9, and that is what the threshold is set
# against. Note this runs at gamma = 0, where liou11's dephasing block is identically zero, so it
# ties the two routes to the same UNITARY evolution and says nothing about the dissipator.
gate("D1a route equality: expm superoperator == pure-state reader at delta=0.15 (gamma=0)",
     abs(dev_expm - dev_pure) < 5e-9,
     f"dev_expm={dev_expm:+.9e}, dev_pure={dev_pure:+.9e}, difference "
     f"{abs(dev_expm - dev_pure):.2e}")
def peak_dev_pure(n_, m_, delta_):
    """Relative dev of the PEAK-read clock vs pi/DeltaE; same route as zero_dev_pure."""
    ev_, vec_ = np.linalg.eigh(ring_h(n_, delta_))
    psi_ = plane_wave(n_, m_)
    c_ = vec_.conj().T @ psi_
    sp_ = pair_split(n_, delta_, m_)[0]
    ts_ = np.linspace(0.6 * np.pi / sp_, 1.4 * np.pi / sp_, 20001)
    F_ = vec_ @ (np.exp(-1j * np.outer(ev_, ts_)) * c_[:, None])
    cur_ = np.sum((F_ * np.roll(F_, -1, axis=0).conj()).imag, axis=0)
    r_ = cur_ / np.sum((psi_ * np.roll(psi_, -1).conj()).imag)
    return abs(ts_[int(np.argmin(r_))] * sp_ / np.pi - 1.0)


# Four bins placed across delta = 0.002 to 0.15, i.e. 1.9 decades end to end, of which the bins
# themselves cover 1.1 (three factor-2 windows and a final factor-1.5 one). They are NOT
# contiguous; STAGE B's OCTAVES and COCT are, and say so. The name is historical: each entry is a
# BIN, and what spans decades is the set of them.
DECADES = [(0.002, 0.004), (0.008, 0.016), (0.03, 0.06), (0.10, 0.15)]
cross_env = [max(abs(zero_dev_pure(n, m, dd_)) / dd_**2
                 for dd_ in np.linspace(lo_, hi_, 15)) for lo_, hi_ in DECADES]
peak_env = [max(peak_dev_pure(n, m, dd_) / dd_
                for dd_ in np.linspace(lo_, hi_, 15)) for lo_, hi_ in DECADES]
# What is flat is the POWER. The COEFFICIENT depends on N: measured at m=1, the crossing's
# |dev|/delta^2 runs 0.20-0.22 (N=8), 0.49-0.50 (N=12), 0.82-0.98 (N=16), 1.32-1.56 (N=20), each a
# range of BIN MAXIMA rather than a point, so a band around 0.21 is
# an N=8 statement and is named as one. The law being gated is that each row is flat, i.e. that the
# exponent is 2 and not 1.
cross_other = {en_: [max(abs(zero_dev_pure(en_, 1, dd_)) / dd_**2
                         for dd_ in np.linspace(lo_, hi_, 9)) for lo_, hi_ in DECADES[:3]]
               for en_ in (12, 16, 20)}
gate("D1b (N=8, m=1) the ZERO CROSSING's envelope is O(delta^2): |dev|/delta^2 flat across decades",
     all(0.15 < c < 0.30 for c in cross_env)
     and all(max(v) / min(v) < 2.0 for v in cross_other.values()),
     "|dev|/delta^2 per bin " + ", ".join(f"{lo_:g}-{hi_:g}:{c:.4f}"
                                             for (lo_, hi_), c in zip(DECADES, cross_env))
     + "; the exponent holds away from N=8 with a larger coefficient, "
     + ", ".join(f"N={en_}:{min(v):.2f}-{max(v):.2f}" for en_, v in cross_other.items()))
gate("D1c (N=8, m=1) the PEAK's envelope is O(delta), one power worse: the crossing wins by ~1/delta",
     all(0.15 < c < 0.30 for c in peak_env),
     "|dev|/delta per bin " + ", ".join(f"{lo_:g}-{hi_:g}:{c:.4f}"
                                           for (lo_, hi_), c in zip(DECADES, peak_env))
     + "; the same coefficient in each hand's own power, so at delta=0.15 the peak reads "
       "~7x coarser and the ratio grows as delta falls. No floor: the smallest decade's peak dev is "
     + f"{max(peak_dev_pure(n, m, dd_) for dd_ in np.linspace(0.002, 0.004, 15)):.1e}.")
print(f"  observation: read the zero crossing, not the peak. The single delta=0.15 peak sample, "
      f"{abs(ts[int(np.argmin(r0))] * sp15 / np.pi - 1.0):.2e}, sits in a local minimum of the "
      f"phase oscillation and is NOT a floor: its own decade's envelope is "
      f"{max(peak_dev_pure(n, m, dd_) for dd_ in np.linspace(0.10, 0.15, 15)):.2e}.")

# D2 the dressing: the zero crossing ADVANCES with gamma, monotonically; committed values.
tz = {0.0: zero_crossing(r0, ts)}
rg = {}
for g in [0.01, 0.05]:
    r = circulation_series(n, m, d, g, tmax, steps)
    r /= r[0]
    rg[g] = r
    tz[g] = zero_crossing(r, ts)
gate("D2 the watching dresses this clock: T_zero(0.05) < T_zero(0.01) < T_zero(0)",
     tz[0.05] < tz[0.01] < tz[0.0],
     f"{tz[0.0]:.4f} -> {tz[0.01]:.4f} -> {tz[0.05]:.4f} (committed 20.295 / 19.347 / 16.477)")
gate("D2 committed pins", abs(tz[0.0] - 20.2953) < 0.02 and abs(tz[0.01] - 19.3470) < 0.02
     and abs(tz[0.05] - 16.4765) < 0.02,
     f"{tz[0.0]:.4f} / {tz[0.01]:.4f} / {tz[0.05]:.4f} against 20.2953 / 19.3470 / 16.4765; a pin "
     f"on three committed numbers, not a law (D2 above is the law)")

# D3 the feedback: measured deepest reversal vs the naive scalar model's own best point.
for g, pin in [(0.01, 1.2510), (0.05, 3.8416)]:
    naive = np.exp(-4 * g * ts) * r0
    ratio = rg[g].min() / naive.min()
    gate(f"D3 gamma={g}: deepest R outlives the naive scalar model, ratio committed {pin}",
         abs(ratio - pin) < 0.02 * pin and ratio > 1.0,
         f"measured {rg[g].min():+.6f} vs naive {naive.min():+.6f}, ratio {ratio:.4f}")

# D4 the control: no crack, no warble (plane wave is an eigenstate; expm float only).
rc = circulation_series(n, m, 0.0, 0.0, tmax, 400)
gate("D4 control delta=0: R pinned at 1", np.max(np.abs(rc / rc[0] - 1.0)) < 1e-10,
     f"max|R-1| = {np.max(np.abs(rc / rc[0] - 1.0)):.2e}")

# D2b/D3c the SCOPE of the two STAGE D findings. Both are read at one geometry (N=8, m=1,
# delta=0.15) and the document generalises their DIRECTION while denying their magnitude. That
# generalisation named a set no script contained, so here it is: six geometries, and what is gated
# is the direction only, because the magnitude visibly does not travel (1.8 to 304 in this table).
scope_rows = []
for (en_, em_, dd_, gg_) in [(8, 1, 0.15, 0.05), (8, 1, 0.05, 0.05), (12, 1, 0.10, 0.02),
                             (12, 2, 0.10, 0.05), (16, 3, 0.20, 0.03), (6, 1, 0.30, 0.05)]:
    t_ = 3.0 * np.pi * en_ / (4 * dd_ * J)
    st_ = min(max(2000, int(t_ / 0.01)), 20000)
    ts_ = np.linspace(0.0, t_, st_ + 1)
    r0_ = circulation_series(en_, em_, dd_, 0.0, t_, st_); r0_ /= r0_[0]
    rg_ = circulation_series(en_, em_, dd_, gg_, t_, st_); rg_ /= rg_[0]
    scope_rows.append((en_, em_, dd_, gg_, zero_crossing(r0_, ts_), zero_crossing(rg_, ts_),
                       rg_.min() / (np.exp(-4 * gg_ * ts_) * r0_).min()))
gate("D2b/D3c both STAGE D findings hold in DIRECTION away from their one geometry, and their "
     "magnitudes do not: the clock advances and the feedback ratio exceeds 1 at all six points, "
     "while the ratio itself spans orders (that spread belongs to the six points chosen, not to a "
     "law; what is the law is the two directions)",
     all(a_ > b_ for *_, a_, b_, _ in scope_rows) and all(r_ > 1.0 for *_, r_ in scope_rows)
     and all(r_ > 1.0 for *_, r_ in scope_rows),
     "; ".join(f"N={n_} m={m_} d={d_:g} g={g_:g}: T_zero {a_:.3f}->{b_:.3f}, feedback {r_:.2f}"
               for n_, m_, d_, g_, a_, b_, r_ in scope_rows))

# D3b the mechanism behind D3, which the document asserted with an unverified snapshot: the
# populations keep a standing 2m-harmonic. A single-time FFT sample is worthless here (the printed
# line below computes how far the instantaneous ratio actually travels), so the reading is the
# whole trajectory, and the control is
# the perfect ring, where there is no crack to sustain a harmonic at all.
def harmonic_weights(n_, m_, delta_, gamma_, tmax=58.0, steps=5800, series=False):
    """RMS over the run of the |FFT| of the population profile at q = m and q = 2m."""
    psi_ = plane_wave(n_, m_)
    rho_ = np.outer(psi_, psi_.conj()).reshape(-1)
    P_ = expm(liou11(n_, delta_, gamma_) * (tmax / steps))
    w1_, w2_ = [], []
    for _ in range(steps + 1):
        f_ = np.abs(np.fft.fft(np.real(np.diag(rho_.reshape(n_, n_))) - 1.0 / n_))
        w1_.append(f_[m_]); w2_.append(f_[2 * m_]); rho_ = P_ @ rho_
    w1_, w2_ = np.array(w1_), np.array(w2_)
    if series:
        return w1_, w2_
    return (float(np.sqrt(np.mean(w1_**2))), float(np.sqrt(np.mean(w2_**2))),
            float(np.mean(w2_ > w1_)))


h1, h2, hfrac = harmonic_weights(8, 1, 0.15, 0.05)
hs1, hs2 = harmonic_weights(8, 1, 0.15, 0.05, series=True)
hs_ratio = hs2[hs1 > 1e-10] / hs1[hs1 > 1e-10]
p1, p2, _ = harmonic_weights(8, 1, 0.0, 0.05)
gate("D3b the crack, and only the crack, plants the 2m-harmonic in the populations that D3's "
     "feedback needs: read over the whole run rather than at one time, and machine zero without "
     "the crack",
     h2 > 1.8 * h1 and hfrac > 0.85 and max(p1, p2) < 1e-12,
     f"RMS over the run: q=1 {h1:.4f}, q=2 {h2:.4f}, ratio {h2/h1:.2f}, and q=2 leads at "
     f"{hfrac*100:.1f}% of the samples; the delta=0 control gives {p1:.1e} and {p2:.1e}. A single "
     f"sample says nothing, and these are computed rather than quoted: the instantaneous ratio is "
     f"{hs2[500]/hs1[500]:.1f} at t=5, {hs2[2200]/hs1[2200]:.0f} at t=22, and reaches "
     f"{hs_ratio.max():.0f} over the run")

# D5 the clock's m-dependence at N=12. STAGE E gives the split's exact m-dependence, so the ratio
# of the two clocks is not free: T_zero = pi/(2*Split) makes it Split_1/Split_2, and all the (1,1)
# block adds on top is its own O(delta^2) crossing envelope. That exponent is what is gated, decade
# bin by bin and worst case INSIDE each bin, never a single sample, for D1b's reason on the same
# object: the deviation oscillates in phase as delta moves. The bins are factor-2 windows placed
# across the range, so what the flatness spans is 2.5 decades of delta, not five of them. The committed 0.9705 is then the
# delta = 0.1 reading of the law rather than a constant of its own, and E6b produces it from the
# curve with nothing propagated.
t12, s12 = 60.0, 6000
ts12 = np.linspace(0, t12, s12 + 1)
tzs = []
for mm in [1, 2]:
    r = circulation_series(12, mm, 0.1, 0.0, t12, s12)
    tzs.append(zero_crossing(r / r[0], ts12))
ratio_m = tzs[1] / tzs[0]
law_ratio = pair_split(12, 0.1, 1)[0] / pair_split(12, 0.1, 2)[0]


def ratio_res12(dd):
    """|measured clock ratio / exact-split ratio - 1| at N=12, off the pure reader (gamma = 0)."""
    return abs((1 + zero_dev_pure(12, 2, dd)) / (1 + zero_dev_pure(12, 1, dd)) - 1.0)


DEC12 = [(0.0005, 0.001), (0.002, 0.004), (0.008, 0.016), (0.03, 0.06), (0.10, 0.15)]
res_env = [max(ratio_res12(dd_) / dd_**2 for dd_ in np.linspace(lo_, hi_, 9)) for lo_, hi_ in DEC12]
own_env = {mm: max(abs(zero_dev_pure(12, mm, dd_)) / dd_**2
                   for dd_ in np.linspace(0.10, 0.15, 9)) for mm in (1, 2)}
gate("D5a the clock's m-ratio rides the EXACT split's, the rest being O(delta^2): the residual's "
     "|res|/delta^2 is flat across five bins spanning 2.5 decades",
     max(res_env) / min(res_env) < 2.0,
     "|res|/delta^2 per bin " + ", ".join(f"{lo_:g}-{hi_:g}:{c:.4f}"
                                             for (lo_, hi_), c in zip(DEC12, res_env))
     + "; each clock's own envelope in the top bin, "
     + ", ".join(f"m={mm}:{v:.3f}" for mm, v in own_env.items())
     + f", against the single delta=0.1 samples "
     + ", ".join(f"m={mm}:{abs(zero_dev_pure(12, mm, 0.1)) / 0.01:.3f}" for mm in (1, 2))
     + " -- the m=1 sample reads about half its bin, which is why the bin maximum is what is gated")
gate("D5b and the delta = 0.1 value of that law is the committed pin 0.9705, on both routes",
     abs(ratio_m - 0.9705) < 0.005
     and abs(ratio_m - law_ratio * (1 + zero_dev_pure(12, 2, 0.1))
             / (1 + zero_dev_pure(12, 1, 0.1))) < 1e-8,
     f"propagated {tzs[0]:.3f} / {tzs[1]:.3f} -> {ratio_m:.7f}; the pure reader, a different time "
     f"grid and crossing interpolation, gives "
     f"{law_ratio * (1 + zero_dev_pure(12, 2, 0.1)) / (1 + zero_dev_pure(12, 1, 0.1)):.7f}; the "
     f"exact split alone {law_ratio:.6f} (E6b), the first-order law 1.000000")

print()
print("=" * 78)
print("STAGE E: the crack is exactly solvable, and the exact law is the CHAIN's scatterer closed")
print("=" * 78)
# STAGES A-D read the crack through first-order degenerate PT. The splitting law does not need to
# be read that way: the whole cracked-ring spectrum is the zero set of one function of k, exactly
# and at every delta (E = 2*J*cos(k), u = J'/J = 1 - delta),
#
#   G(k) = (1 - u^2) sin(N k) cos(k) + [(1 + u^2) cos(N k) - 2 u] sin(k)
#        = sin((N+1) k) - u^2 sin((N-1) k) - 2 u sin(k)  = 0,
#
# from the plane-wave ansatz matched across the cracked bond (psi_{-1} = u psi_{N-1},
# psi_N = u psi_0) with the 2x2 determinant set to zero. u = 1 gives cos(N k) = 1, the perfect
# ring; u = 0 gives sin((N+1) k) = 0, the OPEN N-chain. u is a boundary-condition parameter and
# the crack is the whole road between the ring's modulus N and the chain's modulus N+1.
# The mathematics is standard tight-binding scattering (a rank-2 update of the open chain, a
# one-cell transfer matrix); what the repo did not hold is the equation, and what E4 adds is that
# it is the SAME amplitude COUPLING_DEFECT_WALK_TIME_STEP.md already carried exactly, closed into
# a loop. Scope, because the crack has three closed forms and this is one of them: G speaks to the
# SPLIT (STAGE A). T_rev and the visibility wall need the two-level beat and the (0,1) scalar
# envelope on top of it, and gamma does not appear in G at all.

def secular(k, n, u):
    return ((1 - u**2) * np.sin(n * k) * np.cos(k)
            + ((1 + u**2) * np.cos(n * k) - 2 * u) * np.sin(k))

def secular_dk(k, n, u):
    return ((1 - u**2) * (n * np.cos(n * k) * np.cos(k) - np.sin(n * k) * np.sin(k))
            - (1 + u**2) * n * np.sin(n * k) * np.sin(k)
            + ((1 + u**2) * np.cos(n * k) - 2 * u) * np.cos(k))

E_GRID = [(en, ed) for en in [4, 5, 8, 11, 12, 16, 25, 40]
          for ed in [0.001, 0.02, 0.15, 0.3, 0.7, 0.999]]

def interior_k(n, delta):
    """every level as k = arccos(E/2J). For delta > 0 all N of them are strictly inside the band
    (Perron-Frobenius: the weighted row sums are 2J except at the crack, so the spectral radius is
    strictly below 2J), which the caller asserts rather than assumes."""
    ev = np.linalg.eigvalsh(ring_h(n, delta))
    ks = [np.arccos(e / (2 * J)) for e in ev if abs(e) < 2 * J - 1e-12]
    if len(ks) != n:
        raise AssertionError(f"N={n}, delta={delta}: {len(ks)} of {n} levels inside the band")
    return ks

def count_roots(n, u, points=400001):
    """Sign changes of G on (0, pi); the endpoints are excluded because G vanishes there for free.
    A scan can only MISS roots (a pair closer than the spacing), so the count is taken twice, at
    two resolutions, and disagreement raises rather than returning the smaller number quietly."""
    counts = []
    for pts in (points, 2 * points - 1):
        ks = np.linspace(1e-9, np.pi - 1e-9, pts)
        g = secular(ks, n, u)
        counts.append(int(np.count_nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)))
    if counts[0] != counts[1]:
        raise AssertionError(f"N={n}, u={u}: root count is grid-limited ({counts[0]} vs {counts[1]})")
    return counts[0]

def all_roots(n, u, points=400001):
    """every root of G in (0, pi), bracketed on a fine scan and refined. No matrix anywhere."""
    ks = np.linspace(1e-9, np.pi - 1e-9, points)
    g = secular(ks, n, u)
    idx = np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]
    return [brentq(secular, ks[i], ks[i + 1], args=(n, u), xtol=1e-15, rtol=8.9e-16) for i in idx]

# E1 the curve carries every level, and carries ONLY the levels.
# Error model (no-rounding, case 2): no exact route exists, one side being LAPACK's eigenvalue. What
# is gated is the level shift implied by one Newton step along G, at float-noise size across the
# grid. The figure is NOT the eigensolver's backward error: at delta = 0.001 the leading coefficient
# (1 - u^2) = 2e-3 multiplies a nearly cancelling pair, so evaluating G in double precision costs
# about 40x more than the eigenvalue itself is uncertain by (the cancelling pair is inside the OTHER
# coefficient, (1+u^2)cos(Nk) against 2u, and the near-double root leaves G' ~ 2*N*delta to divide
# by). The gate is therefore read as "float noise", not as a physical residual, and the floor grows
# like 1/delta: the same correct law gives 1.2e-12 at delta=1e-5 and 1.9e-11 at delta=1e-6 AT N=12
# (over the whole grid the worst is 4.2e-12 and 5.6e-11, both at N=4), so the
# 1e-13 threshold belongs to this delta range and would have to move with it.
worst_e1, worst_at = 0.0, None
for en, ed in E_GRID:
    eu = 1 - ed
    for k in interior_k(en, ed):
        dE = abs(2 * J * np.sin(k) * secular(k, en, eu) / secular_dk(k, en, eu))
        if dE > worst_e1:
            worst_e1, worst_at = dE, (en, ed)
gate("E1a every eigenvalue sits on the curve: one Newton step moves it by float noise only",
     worst_e1 < 1e-13, f"max |dE| = {worst_e1:.2e} at N={worst_at[0]}, delta={worst_at[1]}")
# and the other direction, which the Newton metric CANNOT see (it is invariant under G -> G*f):
# the curve has no roots besides the spectrum. Exactly N sign changes on (0, pi), at every (N, delta).
counts = [(en, ed, count_roots(en, 1 - ed)) for en, ed in E_GRID]
bad_count = [c for c in counts if c[2] != c[0]]
gate("E1b and ONLY the eigenvalues: exactly N roots on (0, pi), so the zero set IS the spectrum",
     not bad_count, f"checked {len(E_GRID)} (N, delta) points" if not bad_count else f"{bad_count[:3]}")

# E1c what actually closes "and only the eigenvalues". E1b counts sign changes, and a sign count
# is blind to a root of even order, so it cannot exclude one on its own. The determinant identity
# can: det(2J cos k I - H) = J^N * G(k)/sin k makes G's zero set the characteristic polynomial's,
# multiplicities included. It was asserted in the document and gated nowhere.
det_worst = 0.0
for en_ in (4, 5, 7, 8, 12, 17):
    for dd_ in (0.001, 0.1, 0.5, 0.9):
        for kk_ in (0.3, 0.9, 1.5, 2.2, 2.9):
            lhs_ = np.linalg.det(2 * J * np.cos(kk_) * np.eye(en_) - ring_h(en_, dd_))
            rhs_ = J ** en_ * secular(kk_, en_, 1 - dd_) / np.sin(kk_)
            det_worst = max(det_worst, abs(lhs_ - rhs_) / max(abs(rhs_), 1e-12))
gate("E1c and the identity that makes the sign count conclusive: "
     "det(2J cos k I - H) = J^N G(k)/sin k, so G's zeros ARE the characteristic polynomial's",
     det_worst < 1e-12,
     f"max relative |det - J^N G/sin k| = {det_worst:.1e} over N = 4..17, delta = 0.001..0.9, "
     f"five k each")

# E2 the two ends of delta, read off the EQUATION (not merely off the matrix).
ring_on_curve = max(abs(secular(2 * np.pi * m / 12, 12, 1.0)) for m in range(1, 12))
chain_on_curve = max(abs(secular(np.pi * m / 13, 12, 0.0)) for m in range(1, 13))
# E2a is exactly 0.0 rather than small, and for a reason worth naming: at u=1 the (1-u^2) factor
# kills G's first term outright, so this row cannot see it. E2b (u=0) and E1a (every delta) are
# where that term is exercised.
gate("E2a u=1 solves the curve at the ring's k = 2 pi m/N (structurally blind to G's first term, "
     "which (1-u^2) annihilates here)", ring_on_curve < 1e-14,
     f"max |G| = {ring_on_curve:.1e}")
gate("E2b u=0 solves the curve at the OPEN chain's k = pi m/(N+1)", chain_on_curve < 1e-14,
     f"max |G| = {chain_on_curve:.1e}")
e2n = 12
ring_end = np.max(np.abs(np.sort(np.linalg.eigvalsh(ring_h(e2n, 0.0)))
                         - np.sort(2 * J * np.cos(2 * np.pi * np.arange(e2n) / e2n))))
chain_end = np.max(np.abs(np.sort(np.linalg.eigvalsh(ring_h(e2n, 1.0)))
                          - np.sort(2 * J * np.cos(np.pi * np.arange(1, e2n + 1) / (e2n + 1)))))
gate("E2c and the matrix agrees at both ends (the ring's circulant, the chain's F2b spectrum)",
     max(ring_end, chain_end) < 1e-13, f"max |dE| = {max(ring_end, chain_end):.1e}")

# E2d what the ROAD costs, dynamically. The document says the endpoint does not cost the
# circulation and that what dies, dies at delta = 0+ instead. Those four figures had no verifier:
# E2 is purely spectral and evolves nothing. Three claims, from below.
#   (i)  the spectrum is SIMPLE at every delta > 0 (swept, not sampled; since 2026-09-02 also a theorem for every
#        u >= 0 except u = 1, PROOF_CRACKED_RING_EXACT_CURVE Corollary B, G = 2AB; the "measured rather than
#        proved" below is the sweep's own standing, kept as history: this is what makes every
#        eigenstate a standing wave, and it is measured rather than proved).
#   (ii) every eigenstate's circulation is then EXACTLY zero -- an exact route, so == 0.0.
#   (iii) the persistent part dies at delta = 0+ while the oscillating part does not die at all:
#        the time average collapses by three decades between delta = 0 and 0.1, while the RMS
#        falls smoothly across the whole road with nothing happening at the endpoint. max_t |I|
#        is NOT used, because it is attained at t = 0 where I(0) = -sin(2 pi m/N) for EVERY delta:
#        that number is the seed, not a measurement.
# (i) An ABSOLUTE floor would be the wrong shape for the reason E5a states sixty lines down: the
# minimum gap is the split, 4*delta*J/N, wherever the closest pair is a former degenerate pair, so
# a fixed bound fails on the correct law once delta is small. But the identification itself is only
# a SMALL-delta statement: past delta ~ 0.3 the closest pair stops being a split pair and becomes
# the unpaired level against a pair member, and at the endpoint the gap is the OPEN CHAIN's own
# band-edge spacing. Both regimes have closed forms, so both are gated and neither is a threshold.
def gap_ratio(en_, dd_):
    return np.diff(np.linalg.eigvalsh(ring_h(en_, dd_))).min() / (4 * dd_ * J / en_)


gap_tight = [gap_ratio(en_, dd_) for en_ in range(4, 25)
             for dd_ in np.geomspace(1e-6, 1e-4, 8)]
gap_drift = [gap_ratio(en_, dd_) for en_ in range(4, 25)
             for dd_ in np.linspace(0.01, 0.3, 30)]
gap_end = [(en_, np.diff(np.linalg.eigvalsh(ring_h(en_, 1.0))).min(),
            2 * J * np.cos(np.pi / (en_ + 1)) - 2 * J * np.cos(2 * np.pi / (en_ + 1)))
           for en_ in (12, 24, 40, 60)]
# the simplicity floor is RELATIVE to the split, because an absolute one is the shape this file
# rejects for E5a: the smallest gap on any sweep is bounded by the smallest split on it.
SIMPLE_GRID = [(en_, dd_) for en_ in range(4, 25) for dd_ in np.linspace(0.005, 1.0, 50)]
gap_simple = min(np.diff(np.linalg.eigvalsh(ring_h(en_, dd_))).min() / (4 * dd_ * J / en_)
                 for en_, dd_ in SIMPLE_GRID)
gate("E2d(i) the spectrum is simple at every delta > 0 across the sweep, the minimum gap IS the "
     "split while delta is small (1.000 to four digits below 1e-4, drifting by 0.3), and at the "
     "endpoint it is the open chain's own band-edge spacing, exactly",
     gap_simple > 0.25
     and abs(min(gap_tight) - 1.0) < 1e-3 and abs(max(gap_tight) - 1.0) < 1e-3
     and min(gap_drift) < 0.7
     and max(abs(g_ - f_) for _, g_, f_ in gap_end) < 1e-12,
     f"gap/(4*delta*J/N) over N = 4..24: {min(gap_tight):.6f}..{max(gap_tight):.6f} for delta up to "
     f"1e-4, but {min(gap_drift):.4f}..{max(gap_drift):.4f} over 0.01..0.3, so the identification is "
     f"a small-delta statement and is gated as one; "
     + "at delta=1 the gap is the chain's, 2cos(pi/(N+1))-2cos(2pi/(N+1)): "
     + ", ".join(f"N={en_}:{g_:.6f}/{f_:.6f}" for en_, g_, f_ in gap_end)
     + f"; the smallest gap anywhere on the swept road is {gap_simple:.3f} of its own split (gated "
       f"relative, not against a fixed number), against exactly "
       f"{np.diff(np.linalg.eigvalsh(ring_h(12, 0.0))).min():.1e} at delta=0")

# (ii) had to be rebuilt twice. Reading Im(v_a * conj(v_{a+1})) off numpy's eigenvectors is a dtype
# identity; grouping into eigenspaces fixed only the degenerate side, because for a real matrix with
# a simple spectrum b^T CIRC b is a structural zero too. The hypothesis that actually forbids the
# current is REALITY (time reversal), not simplicity, and the discriminating control is a ring with
# a Peierls FLUX: complex H, perfectly simple spectrum, and it carries circulation. Three inputs,
# two of which must come back nonzero, so the gate can fail on both sides.
CIRC = np.zeros((12, 12), dtype=complex)
for a_ in range(12):
    CIRC[a_, (a_ + 1) % 12] += 0.5j          # <psi|CIRC|psi> = sum_a Im(psi_a conj(psi_{a+1}))
    CIRC[(a_ + 1) % 12, a_] += -0.5j


def ring_h_flux(delta_, phi_):
    """The same cracked ring with a Peierls phase on the wrap bond: H complex, time reversal gone."""
    h_ = ring_h(12, delta_).astype(complex)
    h_[11, 0] = J * (1 - delta_) * np.exp(1j * phi_)
    h_[0, 11] = np.conj(h_[11, 0])
    return h_


def eigenspace_circulation(h_, tol=1e-9):
    """max |<CIRC>| over each eigenspace: what can carry current is a space, not a vector."""
    ev_, vec_ = np.linalg.eigh(h_)
    best_, i_ = 0.0, 0
    while i_ < 12:
        j_ = i_
        while j_ + 1 < 12 and ev_[j_ + 1] - ev_[i_] < tol:
            j_ += 1
        b_ = vec_[:, i_:j_ + 1]
        best_ = max(best_, float(np.abs(np.linalg.eigvalsh(b_.conj().T @ CIRC @ b_)).max()))
        i_ = j_ + 1
    return best_


# The flux control breaks REALITY, but it breaks the chain's reflection symmetry too, so on its own
# it cannot say which of the two forbids the current. The separating control is the other way round:
# a ring with RANDOM REAL bond weights has no reflection symmetry left and keeps reality, and it
# carries no circulation either. (On a ring the two are not independent in general -- reflection
# sends the total flux to minus itself and so forces it to 0 or pi -- which is why the flux control
# alone could never have separated them, and why this second one is needed.)
REFL = np.eye(12)[::-1]


def ring_h_profile(seed_):
    """A ring with random REAL bond weights: reflection gone, time reversal kept."""
    rng_ = np.random.default_rng(seed_)
    w_ = 1.0 + 0.6 * rng_.random(12)
    h_ = np.zeros((12, 12))
    for j_ in range(12):
        h_[j_, (j_ + 1) % 12] = h_[(j_ + 1) % 12, j_] = J * w_[j_]
    return h_


circ_cracked = max(eigenspace_circulation(ring_h(12, dd_)) for dd_ in (1e-6, 0.3, 0.95, 1.0))
circ_perfect = eigenspace_circulation(ring_h(12, 0.0))
circ_flux = eigenspace_circulation(ring_h_flux(0.3, 0.7))
flux_gap = float(np.diff(np.linalg.eigvalsh(ring_h_flux(0.3, 0.7))).min())
prof = [(np.linalg.norm(REFL @ ring_h_profile(k_) @ REFL - ring_h_profile(k_)),
         eigenspace_circulation(ring_h_profile(k_))) for k_ in (7, 11, 23)]
flux_refl = float(np.linalg.norm(REFL @ ring_h_flux(0.3, 0.7) @ REFL - ring_h_flux(0.3, 0.7)))
gate("E2d(ii) what forbids the current is REALITY: not simplicity (the flux ring is simple and "
     "carries it) and not the reflection symmetry the crack keeps (a random REAL bond profile has "
     "none and carries none either)",
     circ_cracked == 0.0 and circ_perfect > 0.9 and circ_flux > 0.5
     and all(r_ > 0.5 and c_ == 0.0 for r_, c_ in prof),
     f"max |<C>| over eigenspaces: cracked (real, simple) {circ_cracked!r} at delta = 1e-6, 0.3, "
     f"0.95, 1; perfect (real, DEGENERATE) {circ_perfect:.6f}; flux (COMPLEX, simple, min gap "
     f"{flux_gap:.4f}, and reflection broken too by {flux_refl:.3f}) {circ_flux:.6f}; three random "
     f"real profiles (reflection broken by "
     + ", ".join(f"{r_:.2f}" for r_, _ in prof) + ") carry "
     + ", ".join(f"{c_!r}" for _, c_ in prof)
     + ". So neither simplicity nor reflection is the hypothesis; reality is")

# (iii) rebuilt too. The window was FIXED while the beat period runs like 1/delta, so the boolean
# measured T/T_zero rather than delta; and the raw RMS at delta = 0 is 100% the persistent part the
# same sentence says dies. Now the window is a fixed number of BEATS, the RMS has its mean removed,
# and -- the point the earlier version missed -- the residual average is finite-window truncation,
# so what is gated is that it falls like 1/n_beats rather than that it sits under some bound.
def road_row(delta_, beats):
    t_win = 4000.0 if delta_ == 0.0 else beats * np.pi * 12 / (4 * delta_ * J)
    ev_, vec_ = np.linalg.eigh(ring_h(12, delta_))
    c_ = vec_.conj().T @ plane_wave(12, 1)
    ts_ = np.linspace(0.0, t_win, 200001)
    F_ = vec_ @ (np.exp(-1j * np.outer(ev_, ts_)) * c_[:, None])
    I_ = np.sum((F_ * np.roll(F_, -1, axis=0).conj()).imag, axis=0)
    return float(I_.mean()), float(np.sqrt(((I_ - I_.mean())**2).mean())), float(I_[0])


road = [(dd_,) + road_row(dd_, 20) for dd_ in (0.0, 1e-3, 5e-3, 0.1, 0.3, 0.95, 1.0)]
decay = [max(abs(road_row(dd_, b_)[0]) for dd_ in (0.1, 0.3, 0.95, 1.0)) for b_ in (20, 80, 320)]
decay_ratios = [decay[1] / decay[0], decay[2] / decay[1]]
gate("E2d(iii) the PERSISTENT circulation is not small at delta > 0, it is ZERO: what the finite "
     "window leaves falls like 1/n_beats (four-fold beats, four-fold smaller), while the OSCILLATING "
     "part takes over at the same amplitude, 0.5 DC -> 0.5/sqrt(2) AC, and only then decays",
     abs(road[0][1]) > 0.49 and road[0][2] < 1e-12
     and all(0.15 < r_ < 0.4 for r_ in decay_ratios)
     and abs(road[1][2] - 0.5 / np.sqrt(2)) < 0.01 and road[-1][2] > 0.2
     and all(abs(r_[3] + 0.5) < 1e-12 for r_ in road),
     "at 20 beats: " + ", ".join(f"delta={d_:g} avg {a_:+.5f} acRMS {s_:.4f}" for d_, a_, s_, _ in road)
     + f"; the worst average over delta>0 is {decay[0]:.5f} at 20 beats, {decay[1]:.5f} at 80 and "
       f"{decay[2]:.5f} at 320, ratios {decay_ratios[0]:.3f} and {decay_ratios[1]:.3f} against the "
       f"1/n law's 0.25; I(0) is {road[0][3]:+.4f} at every delta by construction, which is why "
       f"max_t |I| says nothing here"
)

# E3 STAGE A is this equation truncated. Near k_m the law is exactly
#   delta(2-delta) sin(N q) cos(k_m+q) + sin(k_m+q)[2(1-delta)(cos(N q)-1) + delta^2 cos(N q)] = 0,
# and dropping everything past O(delta^2) leaves the QUADRATIC below in x = N q.
def quadratic(x, km, d):
    return np.sin(km) * x**2 - 2 * d * np.cos(km) * x - d**2 * np.sin(km)

# Carrying a root to energy is E = 2J cos(k_m + x/N); STAGE A keeps the linear term. What the
# quadratic's roots therefore miss is exactly the dispersion's own curvature, -J cos(k_m) (x/N)^2,
# and THAT is gated: not a bound on the gap, but the claim that the gap has that exact coefficient,
# the ratio approaching 1 like delta.
worst_e3root, e3curv = 0.0, []
for ed in [0.05, 0.01, 0.002]:
    worst_ratio = 0.0
    for en in [8, 12, 16, 25]:
        for em in range(1, (en + 1) // 2):
            if 2 * em == en:
                continue
            km = 2 * np.pi * em / en
            for s in (+1, -1):
                x = ed * (np.cos(km) + s) / np.sin(km)
                worst_e3root = max(worst_e3root, abs(quadratic(x, km, ed)) / (ed**2 * np.sin(km)))
                lhs = 2 * J * np.cos(km + x / en) - 2 * J * np.cos(km)     # the root, carried to energy
                rhs = -(2 * ed * J / en) * (np.cos(km) + s)                # STAGE A's committed shift
                curv = -J * np.cos(km) * (x / en)**2                       # the dispersion's own term
                if abs(curv) > 1e-18:
                    worst_ratio = max(worst_ratio, abs((lhs - rhs) / curv - 1.0))
    e3curv.append(worst_ratio)
gate("E3a x_pm = delta(cos k_m +- 1)/sin k_m really are the quadratic's roots (it is evaluated)",
     worst_e3root < 1e-12, f"max |quadratic(x_pm)| / delta^2 sin k_m = {worst_e3root:.1e}")
q1, q2 = e3curv[0] / e3curv[1], e3curv[1] / e3curv[2]
gate("E3b what STAGE A's linear shift misses is exactly the dispersion's curvature -J cos(k_m)(x/N)^2",
     3.5 < q1 < 6.5 and 3.5 < q2 < 6.5,
     f"deviation from that coefficient {e3curv[0]:.4f} / {e3curv[1]:.4f} / {e3curv[2]:.4f} as delta "
     f"falls by 5 each step, ratios {q1:.2f}, {q2:.2f} (the O(delta) law, not a bound)")
# the truncation's own error model: halving delta must quarter the level error. The ratio is the claim.
e3n, e3errs = 12, []
for ed in [0.02, 0.01, 0.005]:
    ev, err = np.sort(np.linalg.eigvalsh(ring_h(e3n, ed))), 0.0
    for em in range(1, (e3n + 1) // 2):
        if 2 * em == e3n:
            continue
        km = 2 * np.pi * em / e3n
        xs = sorted(ed * (np.cos(km) + s) / np.sin(km) for s in (+1, -1))
        pred = sorted(2 * J * np.cos(km + np.array(xs) / e3n))
        true = sorted(ev[np.argsort(np.abs(ev - 2 * J * np.cos(km)))[:2]])
        err = max(err, np.max(np.abs(np.array(pred) - np.array(true))))
    e3errs.append(err)
r1, r2 = e3errs[0] / e3errs[1], e3errs[1] / e3errs[2]
gate("E3c the truncation error is O(delta^2): halving delta quarters it (the law, not a number)",
     3.5 < r1 < 4.5 and 3.5 < r2 < 4.5,
     f"errors {e3errs[0]:.2e} / {e3errs[1]:.2e} / {e3errs[2]:.2e}, ratios {r1:.2f}, {r2:.2f}")

# E4 THE JOIN. COUPLING_DEFECT_WALK_TIME_STEP.md carries the exact transmission amplitude of one
# deformed bond on the infinite chain (the leads in which a plane wave is defined),
# t(q) = -2 i u sin(q) / (e^{-i q} - u^2 e^{i q}). Its inverse is exactly G's two coefficients over
# 2u, so the ring's quantization condition is that amplitude with one round trip of phase:
# Re[e^{-i N k}/t(k)] = 1, which is G/(2 u sin k). The chain reads t for a packet that LEAVES; the
# ring closes it, and the upstream reflection that experiment sets aside ("not read as signal") is,
# on the closed ring, the whole of the warble.
def t_chain(k, u):
    return -2j * u * np.sin(k) / (np.exp(-1j * k) - u**2 * np.exp(1j * k))

worst_e4a = 0.0
for k in np.linspace(0.13, np.pi - 0.13, 41):
    for eu in [0.999, 0.85, 0.7, 0.3, 0.05, 1.05, 1.5, 4.0]:
        worst_e4a = max(worst_e4a, abs(1.0 / t_chain(k, eu)
                                       - ((1 + eu**2) / (2 * eu)
                                          + 1j * (1 - eu**2) / np.tan(k) / (2 * eu))))
gate("E4a 1/t is exactly G's two coefficients over 2u (the algebraic identity, at machine eps)",
     worst_e4a < 1e-12, f"max |diff| = {worst_e4a:.1e} over 41 k x 8 u, on both sides of u = 1")
worst_e4b, at_e4b = 0.0, None
for en, ed in E_GRID:
    for k in interior_k(en, ed):
        r = abs(np.real(np.exp(-1j * en * k) / t_chain(k, 1 - ed)) - 1.0)
        if r > worst_e4b:
            worst_e4b, at_e4b = r, (en, ed)
gate("E4b the ring IS that scatterer closed: Re[e^(-i N k)/t(k)] = 1 on every eigenvalue",
     worst_e4b < 1e-8, f"max |dev| = {worst_e4b:.1e} at N={at_e4b[0]}, delta={at_e4b[1]}; the "
                       "residual tracks 1/(2u) as u -> 0, not the band edge (E1a is the sharp form)")

# E5 the m-dependence. STAGE A's flatness is exact where it is stated: the crack's matrix element
# between the pair is |<psi_m|V|psi_-m>| = 2 delta J/N for EVERY m, because |psi_m(j)|^2 = 1/N. (The
# k-split, by contrast, runs like 2 delta/(N sin k_m); that 1/sin is the dispersion's Jacobian
# dk/dE and says nothing about the crack.) What the exact law adds is the NEXT order, which STAGE D
# could only pin as a measured constant:
#
#   Split_m = (4 delta J/N) [ 1 + delta*(1/2 - 1/(N sin^2 k_m)) + O(delta^2) ]
#
# so the flatness is O(N*delta), not O(delta), and the correction CHANGES SIGN at N sin^2(k_m) = 2.
def c_measured(n, delta, m):
    return (pair_split(n, delta, m)[0] / (4 * delta * J / n) - 1.0) / delta


# The error model IS the gate: c_m is the LEADING term, so its residual must fall one decade per
# decade of delta. An absolute bound would be a number, and a wrong one: the residual's coefficient
# grows with N, so any fixed threshold fails on the correct law once N is pushed (at N=300 the
# residual is 5.5e-3 at delta=1e-5, and 5.3e-2 at the delta=1e-4 this grid actually reaches).
# Odd N and N far outside the earlier sweep are included for exactly that reason.
E5_GRID = [(6, 1), (7, 3), (9, 4), (12, 1), (16, 4), (25, 6), (40, 1), (64, 16),
           (100, 49), (101, 50), (300, 1)]
e5ratios = [(en, em, abs(c_measured(en, 1e-3, em) - c_closed(en, em)),
             abs(c_measured(en, 1e-4, em) - c_closed(en, em))) for en, em in E5_GRID]
e5rat = [r3 / r4 for _, _, r3, r4 in e5ratios]
gate("E5a the leading m-dependence has a closed form, c_m = 1/2 - 1/(N sin^2 k_m): its residual "
     "falls a decade per decade of delta",
     all(8.5 < r < 12.0 for r in e5rat),
     f"decade ratios in [{min(e5rat):.2f}, {max(e5rat):.2f}] over N = 6..300, odd and even; worst "
     f"absolute residual at delta=1e-4 is {max(r4 for *_, r4 in e5ratios):.1e}")
# and the closed form's VALUE, reported where the residual is smallest rather than gated by a number:
e5worst, e5at = 0.0, None
for en in [6, 8, 12, 16, 24, 40, 64, 100]:
    for em in range(1, (en + 1) // 2):
        if 2 * em == en:
            continue
        dd = abs(c_measured(en, 1e-5, em) - c_closed(en, em))
        if dd > e5worst:
            e5worst, e5at = dd, (en, em, c_measured(en, 1e-5, em), c_closed(en, em))
gate("E5b and the value itself: c_measured meets the closed form at every m of eight EVEN N from 6 to 100 (the odd N are gated as a law by E5a, not as a value here)",
     e5worst < 1e-2,
     f"max |c_measured - c_closed| = {e5worst:.1e} at delta=1e-5 (worst N={e5at[0]}, m={e5at[1]}: "
     f"{e5at[2]:.5f} vs {e5at[3]:.5f}); the bound is a report, E5a is the gate")
# the flatness is m-flat only at first order, and the closed form says exactly where it breaks:
# c_m > 0 at the band centre, c_m < 0 at the band edge, crossing at N sin^2 k_m = 2 (for m=1 near
# N = 2 pi^2 ~ 19.7). This is why the measured clock deviation flips sign with N at fixed delta.
signs = [(en, np.sign(c_closed(en, 1)),
          np.sign(pair_split(en, 1e-4, 1)[0] - 4 * 1e-4 * J / en))
         for en in [8, 9, 12, 13, 16, 17, 19, 20, 21, 24, 40]]
gate("E5c the correction changes sign at N sin^2(k_m) = 2, and the measured split follows it",
     all(a == b for _, a, b in signs),
     "m=1: " + ", ".join(f"N={n}:{'+' if a > 0 else '-'}" for n, a, _ in signs))
# and the consequence a reader will care about: the (1,1) CLOCK deviation must flip sign with N
# at fixed delta, opposite to c_m. Gated on the pure-state route (D1a ties it to the superoperator).
# NB the reference is the FIRST-ORDER clock pi*N/(8*delta*J), not the exact split: zero_dev_pure
# normalises by the exact DeltaE and so measures the leftover admixture, a different quantity whose
# sign says nothing about c_m. Rescaling by DeltaE_first/DeltaE_exact puts the c_m term back in.
clock_rows = []
for en in range(8, 26):
    km_ = 2 * np.pi / en
    c_closed_ = c_closed(en, 1)
    first_ = 4 * 0.1 * J / en
    dev_ = (1.0 + zero_dev_pure(en, 1, 0.1)) * first_ / pair_split(en, 0.1, 1)[0] - 1.0
    clock_rows.append((en, c_closed_, dev_))
flips = [en for (en, _, a), (_, _, b) in zip(clock_rows, clock_rows[1:]) if a < 0 <= b]
opposite_big = [(en, cc, dv) for en, cc, dv in clock_rows if abs(cc) >= 0.10]
# the locus itself, solved rather than quoted: c_m(N, 1) = 0 is N sin^2(2 pi/N) = 2.
c_zero_N = brentq(lambda x_: c_closed(x_, 1), 10.0, 40.0)
gate("E5e the clock changes sign with N too, but it is NOT c_m's zero that it crosses",
     len(flips) == 1 and flips[0] == 16
     and all(np.sign(dv) == -np.sign(cc) for _, cc, dv in opposite_big)
     and all(a[2] < b[2] for a, b in zip(clock_rows, clock_rows[1:])),
     f"delta=0.1: T_zero's deviation rises monotonically over N = 8..25 and changes sign once, "
     f"between N = {flips[0]} and {flips[0] + 1}. c_m's own zero is at N ~ {c_zero_N:.2f}, and there "
     f"(N = {int(np.floor(c_zero_N)) - 2}, {int(np.floor(c_zero_N)) - 1}, {int(np.floor(c_zero_N))}) "
     f"the two signs AGREE rather than oppose, because near c_m's zero the O(delta^2) admixture "
     f"outweighs delta*c_m. Where |c_m| >= 0.10 the signs are opposite at every N in the sweep.")

# the first-order flatness itself, exactly, so the committed reading is pinned and not merely quoted:
e5mat, e5val = 0.0, 0.0
for en in [11, 12, 16, 21, 24]:
    vals = []
    for em in range(1, (en + 1) // 2):
        if 2 * em == en:
            continue
        psi = plane_wave(en, em)
        V = np.zeros((en, en)); V[en - 1, 0] = V[0, en - 1] = -0.15 * J
        vals.append(abs(np.conj(psi) @ V @ np.conj(psi)))
    e5mat = max(e5mat, max(vals) / min(vals) - 1.0)
    e5val = max(e5val, max(abs(v - 2 * 0.15 * J / en) for v in vals))   # the VALUE, not just flatness
# What this sees is |psi_m(j)|^2 = 1/N and nothing else: it passes for a symmetric rank-two
# perturbation between ANY two distinct sites, and a diagonal one with the matching coefficient
# passes too, so it carries no information about WHICH bond was cracked. That IS STAGE A's stated
# reason (a point defect is flat in mode space), so the name says the reason and not the ring.
gate("E5d STAGE A's REASON is exact -- the plane wave's flat weight gives "
     "|<psi_m|V|psi_-m>| = 2 delta J/N for every m, to machine eps",
     e5mat < 1e-14 and e5val < 1e-16,
     f"max/min - 1 = {e5mat:.1e}; max |value - 2 delta J/N| = {e5val:.1e}")

# E6 the pin, from the law alone. STAGE D5 pins T_zero(m=2)/T_zero(m=1) = 0.9705 at N=12, delta=0.1
# as a MEASURED constant off a propagated (1,1)-block time series. Since T_zero = pi/(2*Split), the
# exact law must produce it with no propagation and no time grid. The roots are taken off a full
# scan of G, never off a window guessed from the first-order answer, and the run is checked against
# eigvalsh at the same point so the "no matrix needed" claim is gated and not merely asserted.
def split_from_law(n, delta, m):
    roots = all_roots(n, 1 - delta)
    if len(roots) != n:
        raise AssertionError(f"N={n}, delta={delta}: found {len(roots)} roots, expected {n}")
    energies = np.array([2 * J * np.cos(k) for k in roots])
    pair = energies[np.argsort(np.abs(energies - 2 * J * np.cos(2 * np.pi * m / n)))[:2]]
    return abs(pair[0] - pair[1])

s1, s2 = split_from_law(12, 0.1, 1), split_from_law(12, 0.1, 2)
route = max(abs(s1 - pair_split(12, 0.1, 1)[0]), abs(s2 - pair_split(12, 0.1, 2)[0]))
tz_law = [np.pi / (2 * s1), np.pi / (2 * s2)]
gate("E6a the curve alone gives the split, with no matrix: it meets eigvalsh at machine eps",
     route < 1e-14, f"max |law - eigvalsh| = {route:.1e}")
gate("E6b and it reproduces D5's measured pin 0.9705 with nothing propagated",
     abs(s1 / s2 - 0.9705) < 0.005,
     f"law {s1/s2:.6f} vs committed 0.9705, and vs the first-order law's 1.000000; T_zero "
     f"{tz_law[0]:.3f} / {tz_law[1]:.3f} vs measured 46.806 / 45.425, the "
     f"{abs(tz_law[0]-46.806)/46.806*100:.2f}% and {abs(tz_law[1]-45.425)/45.425*100:.2f}% left "
     f"over being the (1,1) block's own O(delta^2) crossing deviation AT THIS delta, 0.232 and "
     f"0.104 in that book; the bin worst case is larger for m=1 (0.45), and D5a gates the "
     f"exponent rather than either number")

# E7 the road past u = 1, and only then the controls. The sibling STRENGTHENS its bond, u = 1+delta,
# and the same curve carries that ring too: the in-band levels directly, and the two levels that
# leave the band on the same curve continued to complex k (E = 2J cos k, k = i*kappa or pi + i*kappa).
worst_e7in, worst_e7out, n_out = 0.0, 0.0, 0
for en in [8, 12, 16, 25]:
    for eu in [1.1, 1.5, 2.0, 4.0]:
        for e in np.linalg.eigvalsh(ring_h(en, 1 - eu)):      # ring_h takes delta, so 1-u
            if abs(e) < 2 * J - 1e-12:
                worst_e7in = max(worst_e7in, abs(secular(np.arccos(e / (2 * J)), en, eu)))
            else:
                n_out += 1
                kap = np.arccosh(abs(e) / (2 * J))
                k = 1j * kap if e > 0 else np.pi + 1j * kap
                worst_e7out = max(worst_e7out, abs(secular(k, en, eu)) / abs(np.sinh(en * kap)))
gate("E7a the SAME curve carries the sibling's STRENGTHENED bond, u = 1+delta, in band",
     worst_e7in < 1e-11, f"max |G| = {worst_e7in:.1e} over u = 1.1..4.0")
gate("E7b and the levels that leave the band, on the same curve continued to complex k",
     worst_e7out < 1e-12, f"max scaled |G| = {worst_e7out:.1e} on {n_out} out-of-band levels")
# how MANY leave is a law, and it is not "two at every delta": even N loses the alternating state at
# -2J immediately, odd N has no such level and its second one leaves only at exactly delta = 2/(N-1).
def out_of_band(n, u):
    return sum(1 for e in np.linalg.eigvalsh(ring_h(n, 1 - u)) if abs(e) > 2 * J + 1e-13)

# The residual of that asymptote is the O((N*eps)^2) truncation of sin(N*eps), so a threshold would
# certify nothing: what says the leading form is RIGHT is that the residual falls by 100 when eps
# falls by 10. Both eps are therefore kept and their RATIO is the gate, not their maximum.
asym = {1e-5: 0.0, 1e-6: 0.0}
for en in list(range(4, 60)) + [101, 200]:
    for uu in (1.3, 2.0, 0.7):
        for eps in (1e-5, 1e-6):
            lead = (-eps * (1 + uu) * ((1 - uu) * en + (1 + uu))) if en % 2 else                    (eps * (1 - uu) * ((1 + uu) * en + (1 - uu)))
            asym[eps] = max(asym[eps], abs(secular(np.pi - eps, en, uu) / lead - 1.0))
asym_ratio = asym[1e-5] / asym[1e-6]
# The two factors are NOT the same object, and that difference IS the odd/even split: for odd N the
# vanishing factor is the BRACKET, [(1-u)N + (1+u)] = 0 at u = (N+1)/(N-1); for even N that bracket
# is 2N at u = 1 and never vanishes for u > 0, and what vanishes there is the PREFACTOR (1-u). So
# the even ring sheds its level at u = 1, i.e. at every delta > 0.
gate("E7c-theorem the band bottom: G ~ -eps(1+u)[(1-u)N+(1+u)] for odd N (its BRACKET zero at "
     "u = (N+1)/(N-1)) and +eps(1-u)[(1+u)N+(1-u)] for even (its PREFACTOR zero at u = 1)",
     90.0 < asym_ratio < 110.0,
     f"|G/leading - 1| worst {asym[1e-5]:.1e} at eps=1e-5 and {asym[1e-6]:.1e} at eps=1e-6, ratio "
     f"{asym_ratio:.1f} against the O(eps^2) law's 100, over N = 4..59, 101, 200 and u = 0.7, 1.3, 2.0")

count_ok = []
for en in [5, 7, 8, 9, 10, 11, 12, 15, 25, 27, 31, 35, 41]:
    thr = 1 + 2.0 / (en - 1)
    if en % 2 == 0:
        count_ok.append((en, all(out_of_band(en, uu) == 2 for uu in (1 + 1e-9, 1.001, 1.5, 3.0))))
    else:
        count_ok.append((en, out_of_band(en, 1 + 1e-9) == 1 and out_of_band(en, thr - 1e-5) == 1
                         and out_of_band(en, thr + 1e-5) == 2 and out_of_band(en, thr + 0.5) == 2))
gate("E7c and the count follows that theorem: 2 out at every delta for even N, 1 for odd "
     "N until delta = 2/(N-1)",
     all(ok for _, ok in count_ok),
     f"counted at N = {[en for en, _ in count_ok]}, sharp to 1e-5 in u; the theorem above covers every other N")
# the controls, measured through E1a's own metric so the control runs through the gate it controls.
def newton_worst(n, d, u_used):
    return max(abs(2 * J * np.sin(k) * secular(k, n, u_used) / secular_dk(k, n, u_used))
               for k in interior_k(n, d))
mut_u = newton_worst(12, 0.15, 1 + 0.15)          # the ring's u is 1-delta; feed the curve 1+delta
def newton_worst_len(n, d, n_used):
    """E1a's own metric with the round trip lengthened: the curve for n_used read on an n-ring."""
    u_ = 1 - d
    return max(abs(2 * J * np.sin(k) * secular(k, n_used, u_) / secular_dk(k, n_used, u_))
               for k in interior_k(n, d))
mut_len_newton = newton_worst_len(12, 0.15, 14)
gate("E7d control: the curve must be fed the ring's OWN u; 1+delta on a 1-delta ring leaves it",
     mut_u > 1e-3, f"Newton step {mut_u:.2e} against E1a's 1e-13")
gate("E7e control: the round trip is N; lengthening it to N+2 leaves the curve (same metric)",
     mut_len_newton > 1e-3, f"Newton step {mut_len_newton:.2e} against E1a's 1e-13")

print()
if FAIL:
    print(f"{len(FAIL)} GATE(S) FAILED:", *FAIL, sep="\n  ")
    raise SystemExit(1)
print("ALL GATES PASS. One crack, three closed forms (4*delta*J/N spectral; pi*N/(4*delta*J)")
print("and exp(-pi*N/(2*Q*delta)) on the (0,1) block), and the (1,1) block's gamma-dressed")
print("clock (advanced zero crossing, diagonal feedback over the naive envelope). The exact")
print("floor-set observable sees the delta = 1e-4 crack; there the reversal read is")
print("exp-blind and even the amplitude read of the early beat deficit is 1.8e-4, under a")
print("1% bar. What decides is whether the dissipator on the block read is a SCALAR, and then")
print("only for the reads whose times are zeros; on the")
print("(1,1) block the two hands differ in POWER: peak ~0.22*delta, crossing ~0.21*delta^2.")
print()
print("And the crack is exactly solvable: one quantization curve carries the whole spectrum")
print("at every delta, the SPLIT is its truncation, delta=1 is the open chain, the curve is")
print("the walk-time chain's own scatterer closed into a loop, and the m-dependence STAGE D")
print("could only pin as a number has the closed form c_m = 1/2 - 1/(N sin^2 k_m).")
