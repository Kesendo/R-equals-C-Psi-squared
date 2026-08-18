"""The corner-beat committed gate (G1/G3/G4/G5 + G2 demos + --frozen).

Companion to experiments/CORNER_BEAT_HARDWARE_PREDICTION.md (v7.2+).
Gate v2.1: two promotion review rounds folded (hop bias channel,
two-stage selection on the verdict metric, per-M retention statistics,
s2 realization error, difference-quantity G4, nominal-basis G5, readout
confusion, sqrt-weighting, the GAUGE PIN in floquet_modes, held-out H0
splits, and the --frozen mode that produces the pre-registration's
governing frozen-configuration numbers from committed code).

What it computes, per the pre-registration's §8:

G1  --g1   Joint working-point optimization (point x grid density x
           shots x gates/block) at the recorded p2 = 0.5%/2q, counts-
           level MC through the committed eigenchannel estimator with
           the REAL finite-M randomized-RZ channel (per-binding phase
           vectors; realization error included by construction).
           Two-stage: screen on the worst f_leak end, confirm the
           shortlist on held-out seeds with the verdict metric
           P(d > theta_D), plus C-prime/W lines and the p2 row.
G3  --g3   Dose certificates on the frozen phase tables: per-site and
           TWO-SITE realized retention over the M bindings (lit/lit,
           lit/dark, dark/dark), the concentrator pass criterion
           (every step within 0.02), and the s2 realization error
           across independent channel realizations vs M.
G4  --g4   Scattered-background budget: site-scattered T1-equivalent
           profiles on BOTH corner arms and U (spurious split), 95th
           percentiles for the band widening (in s2 units).
G5  --g5   J-disorder and J-calibration: s2 under bond-J scatter
           sigma_J (frozen-draw systematic), and d(s2)/d(J*dt).
G2  --g2   Demos: sorted-vs-invariant bias, the (1,2)/(2,3)
           difference-null readouts under a generic background.

Model: one-magnon sector (number conservation exact for the flown
layers), Strang step, per-binding pure-state trajectories with
element-wise RZ phases (coherence (l,m) picks up phi_l - phi_m; dose
sigma_l = sqrt(4 gamma_l dt) gives e^{-2(gamma_l+gamma_m) dt}), gate
error as the pre-registration's bracket: leakage f_leak (post-selected
away, costs counts), within-sector residue half Z-like dephasing at the
Strang gate-count profile (2,3,3,3,3,2), half depolarization toward
Id/6 (applied at the counts layer as a mixing fraction).

Run: python simulations/corner_beat_gate.py [--frozen|--g1|--g2|--g3|--g4|--g5|--all]
"""

import argparse
import numpy as np

N = 6
D6 = 1 << N
DY = [(1, 2), (2, 3), (3, 5)]
GCOUNT = np.array([2, 3, 3, 3, 3, 2], dtype=float)
UNIF = np.ones(N)
CORN = np.zeros(N); CORN[[0, 3, 4]] = 2.0     # maximizing transversal
CORNP = np.zeros(N); CORNP[[0, 1, 2]] = 2.0   # non-maximizing
P2_RECORDED = 0.005          # the recorded ~0.5%/2q
F_LEAK_BRACKET = (8 / 15, 0.9)
SEED = 20260816


def bond_H_1m(bonds):
    H = np.zeros((D6, D6))
    for a in range(D6):
        for l in bonds:
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a, a] += za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2
    cfg = [1 << l for l in range(N)]
    return np.array([[H[a, b] for b in cfg] for a in cfg])


Ho = bond_H_1m([0, 2, 4])
He = bond_H_1m([1, 3])
E0, V0 = np.linalg.eigh(Ho + He)


# --------------------------------------------------------- fit axis
# The axis the channel fits run on, as a SWITCH rather than an
# expression repeated at every fit site. Section 5 of the
# pre-registration registers the refit on the REALIZED-DOSE axis
# t_eff = max(n-1, 0)*dt "in gate and runner together": the final
# injected RZ layer is a Z-basis no-op, so a depth-n cell carries n-1
# injection layers, and fitting on the nominal n*dt leaves depth 0 (the
# fit's highest-weight point) sitting ~5% below the model, biasing
# every rate low and giving away ~16.6% of s2 noiselessly. Depth 0
# carries no injection layer at all and stays at 0 (Amendment 1.4).
#
# This function is the gate half of the machinery item
# `realized_dose_time_axis_refit`. Until it existed the gate built
# `[i*k*dt for i in range(steps//k+1)]` inline at four fit sites, none
# of them switchable, and --certify could not see it: certify hands ONE
# axis array to both fit implementations, so it certified that the two
# FITS agree on a given axis and never that the gate CONSTRUCTS the
# same axis as the runner. Certify now compares the axes themselves.
#
# The default stays "nominal" so this refactor changes no committed
# number; the refreeze throws the switch on both sides at once, and the
# runner's own `time_axis` manifest entry is what forces it.
TIME_AXIS = "nominal"


def fit_time_axis(steps, k, dt, axis=None):
    """Grid times for the channel fits. steps//k + 1 points at spacing
    k Trotter steps."""
    axis = TIME_AXIS if axis is None else axis
    n = [i * k for i in range(steps // k + 1)]
    if axis == "nominal":
        return np.array([s * dt for s in n])
    if axis == "realized_dose":
        return np.array([max(s - 1, 0) * dt for s in n])
    raise RuntimeError(f"unknown time axis: {axis!r}")


def strang(dt):
    Eo, Uo = np.linalg.eigh(Ho)
    Ee, Ue = np.linalg.eigh(He)
    Ah = Uo @ np.diag(np.exp(-1j * Eo * dt / 2)) @ Uo.T
    B = Ue @ np.diag(np.exp(-1j * Ee * dt)) @ Ue.T
    return Ah @ B @ Ah


def floquet_modes(U, dt):
    lam, W = np.linalg.eig(U)
    eps = -np.angle(lam) / dt
    order = []
    for kk in range(N):
        ov = np.abs(V0[:, kk] @ W.conj())
        for cand in np.argsort(-ov):
            if cand not in order:
                order.append(cand); break
    Wm = W[:, order]
    for kk in range(N):
        ph = Wm[np.argmax(np.abs(Wm[:, kk])), kk]
        Wm[:, kk] = Wm[:, kk] * np.exp(-1j * np.angle(ph))
    Wr = np.real(Wm)
    # GAUGE PIN (v2.1): the sign of each Floquet mode is fixed by its
    # overlap with the continuous-H mode, NOT by the largest component
    # (which flips with J*dt and swaps the +- channel labels).
    for kk in range(N):
        if V0[:, kk] @ Wr[:, kk] < 0:
            Wr[:, kk] = -Wr[:, kk]
    return eps[order], Wr


# ------------------------------------------------------- flown channel
def phase_tables(gamma, steps, m_bind, rng):
    """Per-binding i.i.d. phase vectors, sigma_l = sqrt(4 gamma_l dt)
    per step (dt folded into gamma*dt by the caller)."""
    sig = np.sqrt(np.clip(4 * gamma, 0, None))     # gamma already * dt
    return rng.normal(0.0, 1.0, size=(m_bind, steps, N)) * sig[None, None, :]


def run_bindings(U, dt, err_deph_dt, p_hop, steps, k, prep_vec, phases,
                 rng):
    """Vectorized per-binding pure-state trajectories. phases: the frozen
    engineered tables [M, steps, N]. Gate-error within-sector residue,
    BOTH faces (the v2 repair of a bracket that was exactly s2-neutral):
    - Z-like half: extra random phases, variance 4*err_deph_dt per site
      (stochastic unraveling of dephasing at the gate-count profile);
    - hop half: with prob p_hop per step per binding, an XX/YY-type
      two-site error SWAPS the amplitudes on a random bond (a magnon
      hop, arm-dependent, the bias channel the review demanded).
    Returns occupations per grid point per binding: [T, M, N]."""
    M = phases.shape[0]
    Psi = np.tile(prep_vec[None, :], (M, 1)).astype(complex)
    sig_err = np.sqrt(np.clip(4 * err_deph_dt, 0, None))
    out = []
    for s in range(steps + 1):
        if s % k == 0:
            out.append(np.abs(Psi) ** 2)
        Psi = Psi @ U.T
        ph = phases[:, s, :] if s < phases.shape[1] else np.zeros((M, N))
        ph_err = rng.normal(0.0, 1.0, size=(M, N)) * sig_err[None, :]
        Psi = Psi * np.exp(1j * (ph + ph_err))
        hop = rng.random(M) < p_hop
        if hop.any():
            bonds = rng.integers(0, N - 1, size=int(hop.sum()))
            idx = np.where(hop)[0]
            for b in range(N - 1):
                rows = idx[bonds == b]
                if rows.size:
                    Psi[np.ix_(rows, [b, b + 1])] = (
                        Psi[np.ix_(rows, [b + 1, b])])
    return np.array(out)


P01 = 0.015   # readout confusion, symmetric per-qubit (Class-1 bound)
P10 = 0.015


def cell_counts(occ_t, surv, shots, m_bind, rng):
    """Counts for one (arm, depth-grid) cell: occ_t [T, M, N] per-binding
    occupations. Applied at the counts layer: leakage survival,
    readout confusion (a one-magnon shot survives post-selection with
    prob (1-P10)(1-P01)^(N-1) reading the true site; with prob
    P10*P01/... it is rerouted to a random wrong site, else lost).
    Returns post-selected occupations [T, N] and kept weights [T]."""
    T = occ_t.shape[0]
    shots_bind = max(1, shots // m_bind)
    keep_ro = (1 - P10) * (1 - P01) ** (N - 1)
    wrong_ro = P10 * (N - 1) * P01 * (1 - P01) ** (N - 2)
    Yn, Wn = [], []
    for t in range(T):
        p = occ_t[t]                                  # [M, N]
        p = np.clip(p, 0, None)
        p = p / p.sum(axis=1, keepdims=True)
        cum = np.cumsum(p, axis=1)                    # [M, N]
        u = rng.random((p.shape[0], shots_bind))
        idx = (u[:, :, None] > cum[:, None, :]).sum(axis=2)   # [M, shots]
        keep = rng.random(idx.shape) < min(surv[t], 1.0) * keep_ro
        sites = idx[keep]
        wr = rng.random(idx.shape) < min(surv[t], 1.0) * wrong_ro
        if wr.any():
            sites = np.concatenate(
                [sites, rng.integers(0, N, size=int(wr.sum()))])
        if sites.size < 50:
            Yn.append(np.full(N, np.nan)); Wn.append(0.0)
        else:
            cnt = np.bincount(sites, minlength=N).astype(float)
            Yn.append(cnt / sites.size); Wn.append(float(sites.size))
    return np.array(Yn), np.array(Wn)


# ---------------------------------------------------------- estimator
def fit_rate(traces, ts, om, r_max):
    def sse(r, dw):
        tot = 0.0
        for (y, w) in traces:
            ok = ~np.isnan(y) & (w > 0)
            if ok.sum() < 5:
                continue
            t = ts[ok]; yy = y[ok]; ww = w[ok] / max(w[ok].max(), 1.0)
            sw = np.sqrt(ww)
            env = np.exp(-r * t)
            X = np.column_stack([np.ones_like(t),
                                 env * np.cos((om + dw) * t),
                                 env * np.sin((om + dw) * t)])
            try:
                coef, *_ = np.linalg.lstsq(X * sw[:, None], yy * sw,
                                           rcond=None)
            except np.linalg.LinAlgError:
                return np.inf
            tot += float(((yy - X @ coef) ** 2 * ww).sum())
        return tot
    rs = np.linspace(0.0, r_max, 25)
    dws = np.linspace(-0.06, 0.06, 7)
    _, r0, dw0 = min(((sse(r, dw), r, dw) for r in rs for dw in dws))
    for span in (r_max / 12, r_max / 60):
        rs = np.linspace(max(0, r0 - span), r0 + span, 9)
        dws = np.linspace(dw0 - 0.02, dw0 + 0.02, 5)
        _, r0, dw0 = min(((sse(r, dw), r, dw) for r in rs for dw in dws))
    return r0


def s2_of(rates):
    r = np.asarray(rates, dtype=float)
    return float(((r[0] - r[1]) ** 2 + (r[0] - r[2]) ** 2
                  + (r[1] - r[2]) ** 2) / 6)


def arm_s2(Yp, Wp, A, dyf, ts, arm, r_max):
    yd = [6 * Yp[p] @ A for p in range(3)]
    om = float(np.mean(dyf))
    if arm == "U":
        rs = [fit_rate([(yd[p][:, p], Wp[p])], ts, dyf[p], r_max)
              for p in range(3)]
        return s2_of(rs)
    if arm == "C":
        r0 = fit_rate([(yd[0][:, 0], Wp[0])], ts, dyf[0], r_max)
        rp = fit_rate([(yd[1][:, 1] + yd[1][:, 2], Wp[1]),
                       (yd[2][:, 1] + yd[2][:, 2], Wp[2])], ts, om, r_max)
        rm = fit_rate([(yd[1][:, 1] - yd[1][:, 2], Wp[1]),
                       (yd[2][:, 1] - yd[2][:, 2], Wp[2])], ts, om, r_max)
        return s2_of([r0, rp, rm])
    if arm == "Cp":
        r0 = fit_rate([(yd[1][:, 1], Wp[1])], ts, dyf[1], r_max)
        omp = float(np.mean([dyf[0], dyf[2]]))
        rp = fit_rate([(yd[0][:, 0] + yd[0][:, 2], Wp[0]),
                       (yd[2][:, 0] + yd[2][:, 2], Wp[2])], ts, omp, r_max)
        rm = fit_rate([(yd[0][:, 0] - yd[0][:, 2], Wp[0]),
                       (yd[2][:, 0] - yd[2][:, 2], Wp[2])], ts, omp, r_max)
        return s2_of([r0, rp, rm])


# ---------------------------------------------------------------- G1
def config_mc(Q, jdt, steps, k, p2, f_leak, gpb, shots, m_bind, reps, rng,
              null_c_as_u=False):
    gbar = 1.0 / Q
    dt = jdt
    U = strang(dt)
    em, Wm = floquet_modes(U, dt)
    dyf = [em[i] - em[j] for (i, j) in DY]
    A = np.array([[Wm[l, i] * Wm[l, j] for (i, j) in DY] for l in range(N)])
    eps_step = 8 * gpb * p2
    Tn = steps // k + 1
    surv = np.exp(-eps_step * f_leak * k * np.arange(Tn))
    eps_in = eps_step * (1 - f_leak)
    # within-sector residue: half Z-like (phase variance 4*e = 2*eps_in/2
    # per mean site => per-step pair-coherence damping e^{-eps_in/2},
    # the PINNED convention), half XX/YY hop errors (bias-capable)
    err_deph_dt = GCOUNT / GCOUNT.mean() * (eps_in / 2) / 4
    p_hop = eps_in / 2
    ts = fit_time_axis(steps, k, dt)
    r_max = 8 * gbar + eps_in / (2 * dt)
    arms = {"U": UNIF * gbar, "C": CORN * gbar, "Cp": CORNP * gbar}
    if null_c_as_u:
        arms = {"U": UNIF * gbar, "C": UNIF * gbar, "Cp": CORNP * gbar}
    stats = {a: [] for a in arms}
    for _ in range(reps):
        for a, g in arms.items():
            Yp, Wp = [], []
            for d in DY:
                i, j = d
                prep = (Wm[:, i] + Wm[:, j]) / np.sqrt(2)
                ph = phase_tables(g * dt, steps, m_bind, rng)
                occ = run_bindings(U, dt, err_deph_dt, p_hop, steps, k,
                                   prep, ph, rng)
                Y, W = cell_counts(occ, surv, shots, m_bind, rng)
                Yp.append(Y); Wp.append(W)
            # under H0 the C arm is FLOWN with the U profile but ANALYZED
            # with the C-arm estimator (the flown pipeline); the arm label
            # for the estimator is therefore always the nominal one:
            stats[a].append(arm_s2(Yp, Wp, A, dyf, ts, a, r_max))
    return {a: np.array(v) for a, v in stats.items()}


def cmd_g1(quick=False):
    """Two-stage selection (the v2 repair of a winner's-curse pick):
    stage 1 screens all configs at 40 reps, BOTH f_leak ends, on the
    verdict-rule metric P(d > theta_D) with theta_D from each config's
    own H0; stage 2 confirms the top-2 at 200 reps on a HELD-OUT seed.
    Also prints the p2 sensitivity row (Class-1 guard input)."""
    reps1 = 12 if quick else 40
    reps2 = 40 if quick else 200
    print("G1 v2: stage-1 screen (both f_leak ends), metric = worst-end "
          "power and P(d > theta_D)\n")
    configs = []
    for (Q, jdt) in [(10, 0.15), (12, 0.10)]:
        for (steps, k) in [(60, 5), (60, 3)]:
            for shots in [8192, 16384]:
                for gpb in [3, 2]:
                    configs.append((Q, jdt, steps, k, shots, gpb))
    table = []
    for cfg in configs:
        Q, jdt, steps, k, shots, gpb = cfg
        worst_pw, worst_mu, worst_sd = np.inf, None, None
        for fl in F_LEAK_BRACKET:
            rng = np.random.default_rng(SEED + hash(cfg) % 1000)
            st = config_mc(Q, jdt, steps, k, P2_RECORDED, fl, gpb,
                           shots, 1024, reps1, rng)
            ds = st["C"] - st["U"]
            mu, sd = np.nanmean(ds), np.nanstd(ds, ddof=1)
            pw = mu / sd if sd > 0 else np.nan
            if pw < worst_pw:
                worst_pw, worst_mu, worst_sd = pw, mu, sd
        Tn = steps // k + 1
        mins = 4 * Tn * 3 * shots * 0.00032 / 60
        table.append((worst_pw, mins, cfg, worst_mu, worst_sd))
        print(f"  Q={Q:2} dt={jdt:.2f} {Tn:2}x{k} {shots:5} gpb{gpb} | "
              f"worst-end power {worst_pw:5.2f} | {mins:5.1f} min")
    # shortlist: power-per-minute efficient, worst-end power >= 2.5
    short = sorted([r for r in table if r[0] >= 2.5],
                   key=lambda r: r[1])[:2]
    if not short:
        print("\nNO CONFIG reaches worst-end power 2.5: PARK signal.")
        return
    print("\nstage-2 confirmation (held-out seed, %d reps, worst f_leak "
          "per config, verdict metric):" % reps2)
    for (pw1, mins, cfg, _, _) in short:
        Q, jdt, steps, k, shots, gpb = cfg
        # H0 for this config -> theta_D
        rngH = np.random.default_rng(97001 + hash(cfg) % 997)
        best_fl = None; best_pw = np.inf
        for fl in F_LEAK_BRACKET:
            rng = np.random.default_rng(53000 + hash((cfg, fl)) % 997)
            st = config_mc(Q, jdt, steps, k, P2_RECORDED, fl, gpb,
                           shots, 1024, reps2 // 2, rng)
            ds = st["C"] - st["U"]
            pw = np.nanmean(ds) / np.nanstd(ds, ddof=1)
            if pw < best_pw:
                best_pw, best_fl = pw, fl
        st0 = config_mc(Q, jdt, steps, k, P2_RECORDED, best_fl, gpb,
                        shots, 1024, reps2 // 2, rngH, null_c_as_u=True)
        d0 = st0["C"] - st0["U"]
        # ONE threshold from the first H0 half; the second half is the
        # held-out false-rate sample (the frozen-mode template)
        th = 3 * np.nanstd(d0[:len(d0) // 2], ddof=1)
        rng = np.random.default_rng(64000 + hash(cfg) % 997)
        st = config_mc(Q, jdt, steps, k, P2_RECORDED, best_fl, gpb,
                       shots, 1024, reps2, rng)
        ds = st["C"] - st["U"]
        mu, sd = np.nanmean(ds), np.nanstd(ds, ddof=1)
        pdet = float(np.mean(ds > th))
        # held-out H0 split (v2.1): theta_D from the first half, the
        # false rate measured on the second half
        h1, h2 = d0[:len(d0) // 2], d0[len(d0) // 2:]
        th_h = 3 * np.nanstd(h1, ddof=1)
        Tn = steps // k + 1
        print(f"  Q={Q:2} dt={jdt:.2f} {Tn:2}x{k} {shots:5} gpb{gpb} "
              f"fl={best_fl:.2f}: d = {mu:+.5f} +- {sd:.5f} "
              f"(power {mu/sd:.2f}), theta_D {th:.5f}, "
              f"P(d>theta_D) = {pdet:.3f} (n={len(ds)}), "
              f"H0 false rate {np.mean(d0 > th):.3f} (in-sample), "
              f"{np.mean(h2 > th_h):.3f} (held-out)")
        dW = st["C"] - st["Cp"]
        dCp = st["Cp"] - st["U"]
        print(f"      C' lines: s2Cp = {np.nanmean(st['Cp']):+.5f} +- "
              f"{np.nanstd(st['Cp'], ddof=1):.5f}; W exceedance "
              f"s2C - s2Cp = {np.nanmean(dW):+.5f} +- "
              f"{np.nanstd(dW, ddof=1):.5f} (power "
              f"{np.nanmean(dW)/np.nanstd(dW, ddof=1):.2f}); "
              f"s2Cp - s2U = {np.nanmean(dCp):+.5f} +- "
              f"{np.nanstd(dCp, ddof=1):.5f}")
    # p2 sensitivity at the frozen grid (16384 shots = the frozen
    # budget; the 8192 run was measured conservative)
    Q, jdt, steps, k, shots, gpb = (10, 0.15, 60, 3, 16384, 2)
    print("\np2 sensitivity (worst f_leak, %d reps): the Class-1 guard "
          "input" % reps1)
    for p2 in [0.005, 0.0075, 0.01]:
        rng = np.random.default_rng(71000 + int(p2 * 1e4))
        st = config_mc(Q, jdt, steps, k, p2, F_LEAK_BRACKET[0], gpb,
                       shots, 1024, reps1, rng)
        ds = st["C"] - st["U"]
        print(f"  p2={p2:.4f}: d = {np.nanmean(ds):+.5f} +- "
              f"{np.nanstd(ds, ddof=1):.5f} (power "
              f"{np.nanmean(ds)/np.nanstd(ds, ddof=1):.2f})")


def cmd_frozen():
    """The governing frozen-configuration numbers (committed source for
    the pre-registration's §8a/§4/status): Q = 10, J*dt = 0.15, 21 x 3
    grid, 16384 shots/(arm, depth, prep), fractional (2 gates/block),
    M = 1024; both f_leak ends, 200 reps, ONE held-out-half theta_D."""
    Q, jdt, steps, k, shots, gpb = (10, 0.15, 60, 3, 16384, 2)
    print("FROZEN CONFIG: Q=10 dt=0.15 21x3 16384 gpb2 M=1024; 200 reps,"
          " both f_leak ends, held-out-half theta_D")
    for fl in F_LEAK_BRACKET:
        rngH = np.random.default_rng(88001 + int(fl * 100))
        st0 = config_mc(Q, jdt, steps, k, P2_RECORDED, fl, gpb, shots,
                        1024, 100, rngH, null_c_as_u=True)
        d0 = st0["C"] - st0["U"]
        h1, h2 = d0[:50], d0[50:]
        th = 3 * np.nanstd(h1, ddof=1)
        rng = np.random.default_rng(99001 + int(fl * 100))
        st = config_mc(Q, jdt, steps, k, P2_RECORDED, fl, gpb, shots,
                       1024, 200, rng)
        ds = st["C"] - st["U"]
        dW = st["C"] - st["Cp"]
        print(f"fl={fl:.2f}: d = {np.nanmean(ds):+.5f} +- "
              f"{np.nanstd(ds, ddof=1):.5f} "
              f"(power {np.nanmean(ds)/np.nanstd(ds, ddof=1):.2f}), "
              f"theta_D {th:.5f}, P(d>theta_D) = {np.mean(ds > th):.3f},"
              f" H0 false (held-out half) {np.mean(h2 > th):.3f}")
        print(f"        s2Cp = {np.nanmean(st['Cp']):+.5f} +- "
              f"{np.nanstd(st['Cp'], ddof=1):.5f}; W exceedance "
              f"{np.nanmean(dW):+.5f} +- {np.nanstd(dW, ddof=1):.5f} "
              f"(power {np.nanmean(dW)/np.nanstd(dW, ddof=1):.2f})")
    # noiseless dressed references under the pinned estimator (the f_C /
    # f_Cp scaffolding numbers)
    gbar = 1.0 / Q
    U = strang(jdt)
    em, Wm = floquet_modes(U, jdt)
    dyf = [em[i] - em[j] for (i, j) in DY]
    A = np.array([[Wm[l, i] * Wm[l, j] for (i, j) in DY]
                  for l in range(N)])
    ts = fit_time_axis(steps, k, jdt)
    r_max = 8 * gbar
    rng0 = np.random.default_rng(1)
    for arm, prof in [("C", CORN), ("Cp", CORNP), ("U", UNIF)]:
        Yp, Wp = [], []
        for d in DY:
            i, j = d
            prep = (Wm[:, i] + Wm[:, j]) / np.sqrt(2)
            ph = np.zeros((1, steps, N))
            mask_g = prof * gbar * jdt
            # deterministic ensemble evolution: use the exact mask
            rho_mask = np.array([[np.exp(-2 * jdt * (prof[l] + prof[m])
                                          * gbar) if l != m else 1.0
                                  for m in range(N)] for l in range(N)])
            rho = np.outer(prep, prep).astype(complex)
            Y = []
            for s in range(steps + 1):
                if s % k == 0:
                    Y.append(np.real(np.diag(rho)).copy())
                rho = U @ rho @ U.conj().T
                rho = rho * rho_mask
            Yp.append(np.array(Y)); Wp.append(np.ones(len(ts)))
        s2 = arm_s2(Yp, Wp, A, dyf, ts, arm, r_max)
        ideal = {"C": (2 / 3 * gbar) ** 2, "Cp": (2 / 3 * gbar) ** 2 / 3,
                 "U": np.nan}[arm]
        print(f"dressed ref {arm:2}: s2 = {s2:.5f}"
              + (f"  ({s2 / ideal:.2f}x bare ideal)"
                 if np.isfinite(ideal) else ""))


# ---------------------------------------------------------------- G3
def cmd_g3(quick=False):
    """v2: (a) retention PASS probabilities over many table draws (the
    single-draw maxima of v1 are not freezable constants); (b) the
    quantity §8 actually freezes M from: the s2 realization error across
    independent frozen channel realizations, measured on the estimator,
    noiseless counts (channel systematic only)."""
    Q, jdt, steps, k = 10, 0.15, 60, 5
    gbar = 1 / Q
    draws = 10 if quick else 40
    print("G3 v2 (a): retention certificate pass probability, lit/lit "
          "pair, criterion 0.02 (note: transplanted from the concentrator"
          " at a DIFFERENT dose; a per-M criterion is the honest freeze):")
    for m_bind in [256, 512, 1024, 4096]:
        rng = np.random.default_rng(SEED + m_bind)
        g = CORN * gbar * jdt
        devs = []
        for _ in range(draws):
            ph = phase_tables(g, steps, m_bind, rng)
            ret = np.abs(np.mean(np.exp(1j * (ph[:, :, 0] - ph[:, :, 3])),
                                 axis=0))
            devs.append(np.max(np.abs(ret - np.exp(-2 * (g[0] + g[3])))))
        devs = np.array(devs)
        print(f"  M={m_bind:4}: median max|dev| {np.median(devs):.4f} "
              f"range [{devs.min():.4f},{devs.max():.4f}] "
              f"P(PASS@0.02) = {np.mean(devs < 0.02):.2f}")
    print("\nG3 v2 (b): s2(C) realization spread across frozen channels "
          "(estimator, noiseless counts):")
    U = strang(jdt)
    em, Wm = floquet_modes(U, jdt)
    dyf = [em[i] - em[j] for (i, j) in DY]
    A = np.array([[Wm[l, i] * Wm[l, j] for (i, j) in DY] for l in range(N)])
    ts = fit_time_axis(steps, k, jdt)
    r_max = 8 * gbar
    for m_bind in [256, 1024, 4096]:
        rng = np.random.default_rng(SEED + 10 * m_bind)
        vals = []
        for _ in range(max(8, draws // 2)):
            Yp, Wp = [], []
            for d in DY:
                i, j = d
                prep = (Wm[:, i] + Wm[:, j]) / np.sqrt(2)
                ph = phase_tables(CORN * gbar * jdt, steps, m_bind, rng)
                occ = run_bindings(U, jdt, np.zeros(N), 0.0, steps, k,
                                   prep, ph, rng)
                Yp.append(occ.mean(axis=1))
                Wp.append(np.ones(len(ts)))
            vals.append(arm_s2(Yp, Wp, A, dyf, ts, "C", r_max))
        vals = np.array(vals)
        print(f"  M={m_bind:4}: s2(C) {np.mean(vals):.5f} +- "
              f"{np.std(vals, ddof=1):.5f}  (realization sd = "
              f"{100*np.std(vals, ddof=1)/np.mean(vals):.1f}% of centre)")


# ---------------------------------------------------------------- G4
def cmd_g4():
    """v2: report the VERDICT quantity (the difference s2C - s2U under a
    COMMON scattered background, drawn once per trial) at three scatter
    scales, plus the per-arm shifts; propose the T1-CLEAN bound as the
    largest scale keeping the 95% band under half the base difference."""
    rng = np.random.default_rng(SEED)
    Q = 10
    gbar = 1 / Q
    occ = np.einsum('li,lj->lij', V0, V0)

    def comp(gamma):
        M = np.zeros((3, 3))
        for x, (i, j) in enumerate(DY):
            for y, (i2, j2) in enumerate(DY):
                M[x, y] = 2 * sum(
                    gamma[l] * (occ[l, i, i2] * ((j == j2) - occ[l, j, j2])
                                + ((i == i2) - occ[l, i, i2]) * occ[l, j, j2])
                    for l in range(N))
        return M

    def s2c(gamma):
        return s2_of(np.sort(np.linalg.eigvalsh(comp(gamma))))
    baseD = s2c(CORN * gbar) - s2c(UNIF * gbar)
    print("G4 v2: common scattered background (T1 as Gam_l/4-equivalent "
          "dephasing; 'scale' = Gam_mean/4 in units gbar), 4000 draws")
    print(f"  base difference s2C - s2U = {baseD:.4e}\n")
    clean = None
    for scale in [0.125, 0.25, 0.5, 1.0]:
        dd = []
        for _ in range(4000):
            bg = rng.uniform(0, 2 * scale * gbar, N)
            dd.append(s2c(CORN * gbar + bg) - s2c(UNIF * gbar + bg) - baseD)
        dd = np.array(dd)
        lo, hi = np.percentile(dd, [2.5, 97.5])
        ok = max(abs(lo), abs(hi)) < baseD / 2
        frac = 100 * max(abs(lo), abs(hi)) / baseD
        print(f"  scale {scale:5.3f}*gbar: 95% band on Delta(d) "
              f"[{lo:+.2e},{hi:+.2e}]  ({frac:.0f}% of base)  "
              f"{'OK' if ok else 'TOO WIDE'}")
        if ok:
            clean = scale
    print(f"\nproposed T1-CLEAN bound: scattered Gam_l/4 profile scale "
          f"<= {clean}*gbar (95% band under half the base difference)"
          if clean else "\nNO scale passes: CLEAN bound must be tighter"
          " than 0.125*gbar or the band widens")


# ---------------------------------------------------------------- G5
def cmd_g5(quick=False):
    """v2: sigma_J systematic measured at the ESTIMATOR level with the
    preparation and demodulation frozen at the NOMINAL J (the flight's
    situation; v1's compressed shortcut let the basis follow the
    disorder and its base disagreed with §4's authority table)."""
    print("G5 v2: bond-J scatter at the estimator level, nominal-J prep "
          "and demod frequencies, noiseless counts")
    Q, jdt, steps, k = 10, 0.15, 60, 5
    gbar = 1 / Q
    U0 = strang(jdt)
    em, Wm = floquet_modes(U0, jdt)      # NOMINAL basis, frozen
    dyf = [em[i] - em[j] for (i, j) in DY]
    A = np.array([[Wm[l, i] * Wm[l, j] for (i, j) in DY] for l in range(N)])
    ts = fit_time_axis(steps, k, jdt)
    r_max = 8 * gbar
    mask = np.array([[np.exp(-2 * jdt * (CORN[l] + CORN[m]) * gbar)
                      if l != m else 1.0 for m in range(N)]
                     for l in range(N)])

    def s2_disordered(Jb):
        Hs = [bond_H_1m([b]) * Jb[b] for b in range(N - 1)]
        Ho_ = Hs[0] + Hs[2] + Hs[4]
        He_ = Hs[1] + Hs[3]
        Eo, Uo = np.linalg.eigh(Ho_)
        Ee, Ue = np.linalg.eigh(He_)
        Ah = Uo @ np.diag(np.exp(-1j * Eo * jdt / 2)) @ Uo.T
        B = Ue @ np.diag(np.exp(-1j * Ee * jdt)) @ Ue.T
        Ud = Ah @ B @ Ah
        Yp, Wp = [], []
        for (i, j) in DY:
            prep = (Wm[:, i] + Wm[:, j]) / np.sqrt(2)   # NOMINAL orbitals
            rho = np.outer(prep, prep).astype(complex)
            Y = []
            for s in range(steps + 1):
                if s % k == 0:
                    Y.append(np.real(np.diag(rho)).copy())
                rho = Ud @ rho @ Ud.conj().T
                rho = rho * mask
            Yp.append(np.array(Y)); Wp.append(np.ones(len(ts)))
        return arm_s2(Yp, Wp, A, dyf, ts, "C", r_max)

    base = s2_disordered(np.ones(5))
    print(f"  base (nominal J, estimator-dressed): s2(C) = {base:.5f}")
    rng = np.random.default_rng(SEED)
    ndr = 20 if quick else 100
    for sig in [0.01, 0.02]:
        vals = [s2_disordered(1 + rng.normal(0, sig, 5))
                for _ in range(ndr)]
        d = np.array(vals) - base
        print(f"  sigma_J={sig:.2f}: shift median {np.median(d):+.5f} 95% "
              f"[{np.percentile(d, 2.5):+.5f},{np.percentile(d, 97.5):+.5f}]"
              f" ({100*np.percentile(np.abs(d), 95)/base:.0f}% of base)")
    for eps in [0.05]:
        up = s2_disordered(np.ones(5) * (1 + eps))
        dn = s2_disordered(np.ones(5) * (1 - eps))
        print(f"  d(s2)/s2 per +-5% J (all bonds): {up/base-1:+.3f} / "
              f"{dn/base-1:+.3f}")


# ---------------------------------------------------------------- G2
def cmd_g2():
    rng = np.random.default_rng(SEED)
    print("G2 demos: difference-nulls under a generic background")
    occ = np.einsum('li,lj->lij', V0, V0)

    def comp(gamma):
        M = np.zeros((3, 3))
        for x, (i, j) in enumerate(DY):
            for y, (i2, j2) in enumerate(DY):
                M[x, y] = 2 * sum(
                    gamma[l] * (occ[l, i, i2] * ((j == j2) - occ[l, j, j2])
                                + ((i == i2) - occ[l, i, i2]) * occ[l, j, j2])
                    for l in range(N))
        return M
    gbar = 0.1
    for _ in range(3):
        bg = rng.uniform(0, 0.05, N)
        GU = comp(UNIF * gbar + bg)
        GC = comp(CORN * gbar + bg)
        GCp = comp(CORNP * gbar + bg)
        d1 = GC - GU
        d2 = GCp - GU
        print(f"  bg draw: (1,2) row of G(C)-G(U): {np.round(d1[0], 8)} "
              f"(exact 0);  (2,3) row of G(Cp)-G(U): {np.round(d2[1], 8)} "
              f"(exact 0)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--g1", action="store_true")
    ap.add_argument("--frozen", action="store_true")
    ap.add_argument("--g2", action="store_true")
    ap.add_argument("--g3", action="store_true")
    ap.add_argument("--g4", action="store_true")
    ap.add_argument("--g5", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.g1 or a.all:
        cmd_g1(quick=a.quick)
    if a.frozen or a.all:
        cmd_frozen()
    if a.g3 or a.all:
        cmd_g3(quick=a.quick)
    if a.g4 or a.all:
        cmd_g4()
    if a.g5 or a.all:
        cmd_g5(quick=a.quick)
    if a.g2 or a.all:
        cmd_g2()
    if not any([a.g1, a.g2, a.g3, a.g4, a.g5, a.frozen, a.all]):
        print("pick a mode: --frozen --g1 --g2 --g3 --g4 --g5 --all"
              " [--quick]")
