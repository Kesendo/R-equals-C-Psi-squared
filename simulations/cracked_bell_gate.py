"""Gate-first: the CRACKED BELL, the ring's wrap-bond crack read in the time domain.

Verifier for experiments/THE_CRACKED_BELL.md (2026-08-31). The object: XY ring,
H = (J/2) Sum (XX+YY), uniform Z-dephasing gamma per site, wrap bond cracked to
J' = J*(1-delta). Convention: the hopping book of PROOF_RING_GAP_DOMINANCE.md,
single-excitation block H_se = J*A (A = ring adjacency), J = 1 throughout. This is
NOT the D10/W1Dispersion book (isotropic Heisenberg, Laplacian) and not a seat/rate
detuning (the blind-seat arc's open item is a different object).

Two pages of one crack:
  The (0,1) page (vacuum coherence |vac><phi|, STAGES A-C): the dissipator is the
  scalar -2*gamma (Absorption Theorem), so phi(t) = e^(-2*gamma*t) e^(-i*H_se*t) phi(0)
  EXACTLY, and the clock is gamma-free in value. Three first-order laws, derived by
  degenerate perturbation theory and gated against exact eigendecomposition:
    (A) SPLITTING: every m <-> N-m pair splits by the SAME Delta_E = 4*delta*J/N
        (|psi_m(j)|^2 = 1/N: a point defect is flat in mode space); the partners are
        the cut-adapted cos/sin standing waves, node/antinode ON the cracked bond,
        shifted by -(2*delta*J/N)*(cos(theta_m) +/- 1) in opposite directions.
    (B) REVERSAL: a launched traveling wave of the PERFECT ring fully reverses at
        T_rev = pi/Delta_E = pi*N/(4*delta*J); this is the walk-time step's discarded
        O(delta) reflection (COUPLING_DEFECT_WALK_TIME_STEP: "not read as signal")
        resonantly accumulated by the closed ring until it is the whole signal.
    (C) VISIBILITY: the reversed amplitude carries the exact envelope
        e^(-2*gamma*T_rev) = exp(-pi*N/(2*Q*delta)), Q = J/gamma -- exact in the
        gamma-factorization, exponent accurate to O(delta) relative. Stage C first
        PROVES the two books from below on the full 4^N Lindbladian of the cracked
        ring at N=5 (the (0,1) page pays exactly -2*gamma, entry-wise; the (1,1)
        page exactly -4*gamma off-diagonal, 0 diagonal), then gates the exponent's
        O(delta) law across a delta decade. The exact floor-set observable of
        PROOF_RING_GAP_DOMINANCE sees delta = 1e-4; there the REVERSAL read is
        exp-blind (10^-1364 at N=4, Q=20) and even the most sensitive early-time
        read of the beat deficit, ~0.068*(4*delta*Q/N)^2 at t = 1/gamma, is 2.7e-7:
        polynomially blind, so everything physical is blind at that crack.
  The (1,1) page (single-excitation density, Haken-Strobl, STAGE D -- the block
  MirrorWorld's Cone runs, run mode `warble N [delta]`, pinned by WarbleTests):
  dephasing is NOT scalar here (off-diagonals pay -4*gamma, the diagonal pays
  nothing), and the same crack's clock is gamma-DRESSED: the circulation's zero
  crossing advances under the watching, and the deepest reversal outlives the naive
  scalar model e^(-4*gamma*t)*R_0(t) at its own best point (the dephasing-free
  diagonal feeds the current back). Which page is read decides what the watching
  does to the clock; the fast walk-time front never faced this.

Error model (no-rounding convention, case 2): PT truncates at first order, so
deviations are O(delta^2) on Delta_E (i.e. O(delta) on Delta_E/delta, so the decade
ratio of e(delta) must sit near 0.1, and THAT is gated) and O(delta) on the clocks;
peak reads additionally carry the delta-independent flat-hand floor (~2e-3), stated
where gated. The scalar-envelope identity in (C) is case 1 and is asserted as an
entry-wise residual on the full Lindbladian, not imported.
STAGE D's dressed-clock and feedback numbers have no closed form; they are pinned
as committed constants (the WalkTime-plateau genre), and the cross-book gate D1
ties the density-matrix beat to STAGE A's spectral splitting through two
independent routes (superoperator expm vs eigh).
"""

import numpy as np
from scipy.linalg import expm

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
for n in NS:
    ms = [m for m in range(1, (n + 1) // 2) if 2 * m != n]
    ratios = []
    for m in ms:
        # e(delta) = |Delta_E/delta - 4/N| must shrink LINEARLY: each decade ratio ~ 0.1.
        e = []
        for d in DELTAS:
            split, _ = pair_split(n, d, m)
            e.append(abs(split / d - 4.0 * J / n))
        floor = 1e-9  # eigensolver noise on Delta_E/delta at delta=1e-5
        if e[1] > floor and e[2] > floor:
            ratios += [e[0] / e[1], e[1] / e[2]]
        split, _ = pair_split(n, 1e-5, m)
        rel = abs(split / (4.0 * 1e-5 * J / n) - 1.0)
        if rel >= 1e-3:
            gate(f"N={n} m={m} first-order value", False, f"rel={rel:.2e}")
    # the LAW gate: linearity means the decade ratio is 0.1; a quadratic (0.01), a
    # square root (0.32) or a constant (1.0) all fail the band.
    ok = all(0.05 < r < 0.2 for r in ratios) if ratios else True
    det = (f"decade ratios in [{min(ratios):.3f}, {max(ratios):.3f}], all pairs"
           if ratios else "all below noise floor")
    gate(f"N={n}: e(delta) shrinks one decade per delta decade (linear law)", ok, det)

print()
print("branch shifts (N=12): even-across-the-crack partner by -(2d/N)(cos+1), odd by -(2d/N)(cos-1)")
n, d = 12, 1e-4
for m in [1, 2, 5]:
    _, shifts = pair_split(n, d, m)
    th = 2.0 * np.pi * m / n
    pred = np.sort(np.array([-(2 * d * J / n) * (np.cos(th) + 1.0),
                             -(2 * d * J / n) * (np.cos(th) - 1.0)]))
    rel = np.max(np.abs(shifts - pred)) / (4 * d * J / n)
    gate(f"N=12 m={m} branch shifts match PT", rel < 5e-3, f"rel={rel:.2e}")
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
print("STAGE B: the reversal  T_rev = pi/Delta_E = pi*N/(4*delta*J), gamma-free (0,1) page")
print("=" * 78)
for (n, m) in [(12, 1), (8, 3), (16, 2)]:
    for d in [1e-3, 3e-3, 1e-2]:
        h = ring_h(n, d)
        ev, vec = np.linalg.eigh(h)
        phi0 = plane_wave(n, m)
        psi_minus = plane_wave(n, -m)
        c0 = vec.conj().T @ phi0
        split, _ = pair_split(n, d, m)
        t_pred = np.pi / split
        ts = np.linspace(0.0, 1.4 * t_pred, 4001)
        amp = (psi_minus.conj() @ (vec @ (np.exp(-1j * np.outer(ev, ts)) * c0[:, None])))
        p_minus = np.abs(amp) ** 2
        t_rev = ts[int(np.argmax(p_minus))]
        dev = abs(t_rev * split / np.pi - 1.0)
        # a PEAK read: the flat-hand floor (~2e-3, delta-independent) plus O(delta) leakage.
        gate(f"N={n} m={m} d={d:.0e}: T_rev*DeltaE/pi = 1 (peak read: floor + O(delta))",
             dev < 2e-3 + 0.5 * d,
             f"dev={dev:.2e}, peak P_-={p_minus.max():.6f}")

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


# C1: both books, entry-wise. The closed forms rest on the (0,1) page paying exactly
# the scalar -2*gamma (Absorption Theorem) and the Warble's (1,1) page paying exactly
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
evC, vecC = np.linalg.eigh(ring_h(nC, dC))
worst_mode = 0.0
for k in range(nC):
    O = np.zeros((dimC, dimC), dtype=complex)
    for jj in range(nC):
        O[0, site_idx[jj]] = np.conj(vecC[jj, k])     # |vac><phi_k|
    v = O.reshape(-1)
    res = LF @ v - (1j * evC[k] - 2.0 * gC) * v
    worst_mode = max(worst_mode, np.max(np.abs(res)))
gate("C1b every (0,1) mode pays exactly -2*gamma + i*E_k on the full Lindbladian",
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
gate("C1c the (1,1) page pays exactly -4*gamma off-diagonal, 0 diagonal (Haken-Strobl)",
     worst_hs == 0.0, f"worst residual = {worst_hs:.1e}")

# C2: the visibility wall's exponent law. With the envelope proven in C1, the reversed
# amplitude at t* = pi/Delta_E_exact is amp0 * e^(-2*gamma*t*); against the closed form
# exp(-pi*N/(2*Q*delta)) the whole deviation is the PT error of Delta_E inside the
# EXPONENT: |ln(measured/closed)| / (2*gamma*t_first) = |t* - t_first|/t_first = O(delta),
# gated across a delta decade (the ratio itself is delta-independent in the amplitude,
# which is why an amplitude-ratio tolerance would be the wrong model).
gamma = 0.05
n, m = 12, 1
exp_rels = []
for d in [1e-2, 1e-3]:
    split, _ = pair_split(n, d, m)
    t_star = np.pi / split
    t_first = np.pi * n / (4.0 * d * J)
    exp_rels.append(abs(t_star - t_first) / t_first)
gate(f"C2 visibility exponent is O(delta) relative: dev(1e-2), dev(1e-3), decade ratio",
     exp_rels[0] < 0.5e-2 * 5 and exp_rels[1] < 0.5e-3 * 5
     and 5.0 < exp_rels[0] / exp_rels[1] < 20.0,
     f"{exp_rels[0]:.2e}, {exp_rels[1]:.2e}, ratio {exp_rels[0] / exp_rels[1]:.1f}")
d = 1e-2
split, _ = pair_split(n, d, m)
wall_form = np.exp(-np.pi * gamma * n / (2.0 * d * J))
amp_ratio = np.exp(-2.0 * gamma * np.pi / split) / wall_form
print(f"  at N={n}, d={d}, Q={J / gamma:.0f}: envelope at t* over closed form = "
      f"{amp_ratio:.6f} (the O(delta) exponent error amplified by 2*gamma*t_first = "
      f"{2 * gamma * np.pi * n / (4 * d * J):.1f}), amplitude = "
      f"{np.exp(-2.0 * gamma * np.pi / split):.3e}")

print()
print("the visibility wall (reversed amplitude = 10^-x), Q = J/gamma:")
print(f"{'N':>4} {'delta':>8} {'Q':>6} {'x = pi*N/(2*Q*delta)/ln10':>28}")
for (n_, d_, q_) in [(4, 1e-4, 20), (4, 1e-3, 20), (4, 0.068, 20),
                     (12, 1e-2, 20), (12, 1e-2, 200), (32, 1e-3, 200)]:
    x = np.pi * n_ / (2.0 * q_ * d_) / np.log(10.0)
    print(f"{n_:>4} {d_:>8.0e} {q_:>6} {x:>28.1f}")

print()
print("=" * 78)
print("STAGE D: the (1,1) page (Haken-Strobl, the Cone's block): the gamma-DRESSED clock")
print("=" * 78)


def liou11(n, delta, gamma):
    h = ring_h(n, delta)
    eye = np.eye(n)
    lh = -1j * (np.kron(h, eye) - np.kron(eye, h.T))
    dd = np.diag([-4.0 * gamma if a != b else 0.0 for a in range(n) for b in range(n)])
    return lh + dd


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

# D1 cross-book: the (1,1) gamma=0 clock against STAGE A's spectral splitting.
# The clock is read at the ZERO CROSSING, the steep hand. The O(delta) fast admixture
# shifts the read within an O(delta) ENVELOPE whose phase oscillates with delta
# (T_zero ~ 1/delta wraps the fast phases), so a two-point ratio cannot see the law;
# the gate is the envelope |dev| <= C*delta across a delta sweep (measured C ~ 0.021),
# after a route-equality gate ties the fast pure-state reader to the superoperator.
# (The reversal PEAK is the flat hand: the same admixture on a curvature falling as
# delta^2 leaves a delta-INDEPENDENT relative floor ~2e-3; printed as an observation.)


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
gate("D1a route equality: expm superoperator == pure-state reader at delta=0.15",
     abs(dev_expm - dev_pure) < 5e-4,
     f"dev_expm={dev_expm:+.2e}, dev_pure={dev_pure:+.2e}")
sweep = [0.01, 0.02, 0.0375, 0.075, 0.11, 0.15, 0.2]
devs = [zero_dev_pure(n, m, dd_) for dd_ in sweep]
worst = max(abs(v) / dd_ for v, dd_ in zip(devs, sweep))
gate("D1b cross-book envelope: |T_zero*2*DeltaE/pi - 1| <= 0.03*delta across the sweep",
     worst < 0.03,
     "signed devs " + ", ".join(f"{dd_:g}:{v:+.1e}" for dd_, v in zip(sweep, devs))
     + f"; worst |dev|/delta = {worst:.4f}")
t_peak_dev = abs(ts[int(np.argmin(r0))] * sp15 / np.pi - 1.0)
print(f"  observation: the reversal PEAK is the flat hand -- its relative dev "
      f"({t_peak_dev:.2e} at delta=0.15) sits on the delta-independent "
      f"wiggle-to-curvature floor ~2e-3. Read the zero crossing, not the peak.")

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
     and abs(tz[0.05] - 16.4765) < 0.02)

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

# D5 m-flatness in time (N=12, delta=0.1): the same clock for every m, to O(delta).
t12, s12 = 60.0, 6000
ts12 = np.linspace(0, t12, s12 + 1)
tzs = []
for mm in [1, 2]:
    r = circulation_series(12, mm, 0.1, 0.0, t12, s12)
    tzs.append(zero_crossing(r / r[0], ts12))
ratio_m = tzs[1] / tzs[0]
gate("D5 m-flat clock, the law: |T_zero(m=2)/T_zero(m=1) - 1| = O(delta)",
     abs(ratio_m - 1.0) < 0.6 * 0.1,
     f"{tzs[0]:.3f} vs {tzs[1]:.3f}, ratio {ratio_m:.4f}")
gate("D5 m-flat clock, the pin: committed 0.9705 (the genuine second-order m-dependence)",
     abs(ratio_m - 0.9705) < 0.005)

print()
if FAIL:
    print(f"{len(FAIL)} GATE(S) FAILED:", *FAIL, sep="\n  ")
    raise SystemExit(1)
print("ALL GATES PASS. One crack, three closed forms (4*delta*J/N spectral; pi*N/(4*delta*J)")
print("and exp(-pi*N/(2*Q*delta)) on the (0,1) page), and the (1,1) page's gamma-dressed")
print("clock (advanced zero crossing, diagonal feedback over the naive envelope). The exact")
print("floor-set observable sees the delta = 1e-4 crack; there the reversal read is")
print("exp-blind and even the early-time deficit read (~0.068*(4*delta*Q/N)^2) is 2.7e-7.")
print("Which page is read decides what the watching does to the clock.")
