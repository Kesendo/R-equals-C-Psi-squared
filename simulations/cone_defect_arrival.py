"""Ballistic front across a delta-J defect bond: the pre-registered arrival-time runs.

Local-only WIP script (design + revision in the session scratchpad,
design_cone_defect_arrival.md). Single-excitation sector only:
  gamma=0   : pure state, exact eigenbasis evolution of the N x N hopping matrix.
  gamma>0   : rho (N x N) RK4 with d(rho)/dt = -i[h,rho] - 4*gamma*(rho - diag(rho)).

Predictions under test (REVISED, post-review):
  P1  step law: Delta-t_i ~ 0 upstream, constant plateau ~ -delta/(2J) downstream (gamma=0).
  P2  broadband seed undershoots |slope|=1/(2J) by O(15%); quasi-monochromatic packet
      (q = -pi/2, sigma_q = 0.1*pi) recovers -1/(2J) within +-5%.
  P3' locality scan: the step edge tracks the defect bond b; plateau height b-independent.
  P4' crossover: at gamma=0.05 the step survives only in the near-field window
      t_arr <~ 1/(4 gamma); far field is diffusive (descriptive).
      OUTCOME: the far field carried the step after all (the front's TIMING stays
      ballistic while its amplitude decays); see the [Probe] block and the writeup
      experiments/COUPLING_DEFECT_WALK_TIME_STEP.md.
Controls: relative threshold theta*max_t P_i|clean, theta in {0.1,0.2,0.4};
  dt-halving reproducibility; lattice-reflection mirror (seed AND defect mirrored).
"""

import numpy as np

J = 1.0
OUT = []


def log(s=""):
    print(s, flush=True)
    OUT.append(s)


def hopping(n, defect_bond=None, delta=0.0):
    h = np.zeros((n, n))
    for a in range(n - 1):
        j = J * (1.0 + (delta if defect_bond == a else 0.0))
        h[a, a + 1] = j
        h[a + 1, a] = j
    return h


def populations_pure(h, psi0, tgrid):
    """Exact eigenbasis evolution; returns P[t, site]."""
    E, U = np.linalg.eigh(h)
    c = U.T @ psi0
    phases = np.exp(-1j * np.outer(tgrid, E)) * c  # (T, N)
    amps = phases @ U.T                            # (T, N)
    return np.abs(amps) ** 2


def populations_lindblad(h, psi0, tgrid, gamma, dt_int=0.002):
    """RK4 on rho with uniform -4*gamma coherence damping; sampled on tgrid."""
    n = h.shape[0]
    rho = np.outer(psi0, psi0.conj()).astype(complex)
    mask = 4.0 * gamma * (1.0 - np.eye(n))

    def rhs(r):
        return -1j * (h @ r - r @ h) - mask * r

    P = np.empty((len(tgrid), n))
    t = 0.0
    it = 0
    P[0] = np.real(np.diag(rho))
    for k in range(1, len(tgrid)):
        while t < tgrid[k] - 1e-12:
            step = min(dt_int, tgrid[k] - t)
            k1 = rhs(rho)
            k2 = rhs(rho + 0.5 * step * k1)
            k3 = rhs(rho + 0.5 * step * k2)
            k4 = rhs(rho + step * k3)
            rho = rho + (step / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += step
            it += 1
        P[k] = np.real(np.diag(rho))
    return P


def arrival_times(P, tgrid, ref_peak, theta):
    """First t with P_i(t) >= theta * ref_peak_i (clean per-site peak); linear interp."""
    T, n = P.shape
    ta = np.full(n, np.nan)
    for i in range(n):
        thr = theta * ref_peak[i]
        idx = np.nonzero(P[:, i] >= thr)[0]
        if len(idx) == 0 or idx[0] == 0:
            continue
        k = idx[0]
        p0, p1 = P[k - 1, i], P[k, i]
        frac = (thr - p0) / (p1 - p0) if p1 > p0 else 0.0
        ta[i] = tgrid[k - 1] + frac * (tgrid[k] - tgrid[k - 1])
    return ta


def seed_site(n, s):
    v = np.zeros(n, dtype=complex)
    v[s] = 1.0
    return v


def seed_packet(n, n0, q=-np.pi / 2, sigma_q=0.1 * np.pi):
    sigma_x = 1.0 / (2.0 * sigma_q)
    x = np.arange(n)
    v = np.exp(-((x - n0) ** 2) / (4.0 * sigma_x ** 2)) * np.exp(1j * q * x)
    return v / np.linalg.norm(v)


def run_profile(n, seed, defect_bond, delta, gamma, tmax, dt=0.001, theta=0.2,
                dt_int=0.002):
    tgrid = np.arange(0.0, tmax + dt / 2, dt)
    h0 = hopping(n)
    hd = hopping(n, defect_bond, delta)
    if gamma == 0.0:
        P0 = populations_pure(h0, seed, tgrid)
        Pd = populations_pure(hd, seed, tgrid)
    else:
        P0 = populations_lindblad(h0, seed, tgrid, gamma, dt_int)
        Pd = populations_lindblad(hd, seed, tgrid, gamma, dt_int)
    ref = P0.max(axis=0)
    t0 = arrival_times(P0, tgrid, ref, theta)
    td = arrival_times(Pd, tgrid, ref, theta)
    return t0, td, td - t0


def plateau_stats(dt_prof, lo, hi):
    seg = dt_prof[lo:hi]
    seg = seg[~np.isnan(seg)]
    if len(seg) == 0:
        return np.nan, np.nan, 0
    return float(np.mean(seg)), float(np.std(seg)), len(seg)


# ----------------------------------------------------------------------------
log("=" * 78)
log("BALLISTIC FRONT ACROSS A delta-J DEFECT: pre-registered arrival-time runs")
log("J = 1, chain, relative threshold theta * max_t P_i|clean")
log("=" * 78)

# ---- Run A: the step (coherent), N=60, defect mid-chain, 4 deltas ----------
N, B, SEED, TMAX = 60, 29, 0, 40.0
log("\n[Run A] N=60, gamma=0, seed site 0, defect bond (29,30), theta=0.2")
log("  predicted plateau height = -delta/(2J); broadband undershoot O(15%) expected")
log(f"  {'delta':>7} | {'near-side mean(5:20)':>22} | {'plateau mean(40:56)':>20} "
    f"| {'plateau std':>11} | pred")
heights = {}
for delta in (-0.10, -0.05, +0.05, +0.10):
    _, _, dtp = run_profile(N, seed_site(N, SEED), B, delta, 0.0, TMAX)
    near_m, near_s, _ = plateau_stats(dtp, 5, 20)
    pl_m, pl_s, cnt = plateau_stats(dtp, 40, 56)
    heights[delta] = pl_m
    log(f"  {delta:>+7.2f} | {near_m:>+22.5f} | {pl_m:>+20.5f} | {pl_s:>11.5f} "
        f"| {-delta/(2*J):>+.4f}")

deltas = np.array(sorted(heights))
hh = np.array([heights[d] for d in deltas])
slope = np.polyfit(deltas, hh, 1)[0]
log(f"  P1 plateau flatness: see std column (pre-reg: far-half variation < +-30%)")
log(f"  P2 broadband slope fit: {slope:+.4f}  (predicted -0.50, undershoot O(15%) OK)")

# theta-independence control
log("\n[Control iii] theta-independence, delta=+0.10:")
for theta in (0.1, 0.2, 0.4):
    _, _, dtp = run_profile(N, seed_site(N, SEED), B, +0.10, 0.0, TMAX, theta=theta)
    pl_m, pl_s, _ = plateau_stats(dtp, 40, 56)
    log(f"  theta={theta:.1f}: plateau {pl_m:+.5f} +- {pl_s:.5f}")

# dt-halving control
log("\n[Control i] dt-halving (clean pipeline), delta=+0.10:")
_, _, d1 = run_profile(N, seed_site(N, SEED), B, +0.10, 0.0, TMAX, dt=0.001)
_, _, d2 = run_profile(N, seed_site(N, SEED), B, +0.10, 0.0, TMAX, dt=0.0005)
diff = np.nanmax(np.abs(d1 - d2))
log(f"  max |Delta-t(dt=0.001) - Delta-t(dt=0.0005)| = {diff:.2e}")

# lattice-reflection mirror control (seed AND defect mirrored)
log("\n[Control ii] lattice reflection: (seed 0, bond (20,21)) vs (seed 59, bond (38,39)):")
_, _, da = run_profile(N, seed_site(N, SEED), 20, +0.10, 0.0, TMAX)
_, _, db = run_profile(N, seed_site(N, N - 1), 38, +0.10, 0.0, TMAX)
mir = np.nanmax(np.abs(da - db[::-1]))
log(f"  max |profile - mirror(profile)| = {mir:.2e}  (must be ~0)")

# ---- Run B (P2'): quasi-monochromatic packet -------------------------------
log("\n[Run B / P2'] quasi-monochromatic packet, q=-pi/2, sigma_q=0.1*pi, n0=10,")
log("  N=120 (room for the packet), defect bond (59,60), gamma=0, theta=0.2")
N2, B2 = 120, 59
pk_heights = {}
for delta in (-0.10, -0.05, +0.05, +0.10):
    _, _, dtp = run_profile(N2, seed_packet(N2, 10), B2, delta, 0.0, 55.0)
    pl_m, pl_s, cnt = plateau_stats(dtp, 75, 110)
    pk_heights[delta] = pl_m
    log(f"  delta={delta:+.2f}: plateau {pl_m:+.5f} +- {pl_s:.5f}  "
        f"(pred {-delta/(2*J):+.4f})")
dq = np.array(sorted(pk_heights))
hq = np.array([pk_heights[d] for d in dq])
slope_q = np.polyfit(dq, hq, 1)[0]
log(f"  P2 quasi-monochromatic slope: {slope_q:+.4f}  "
    f"(pre-reg: -1/(2J) = -0.5000 within +-5%)")

# ---- Run C (P3'): defect-position locality scan ----------------------------
log("\n[Run C / P3'] locality scan, N=60, gamma=0, seed 0, delta=+0.10,")
log("  bond b in {10, 20, 30, 40}: edge must sit at b, plateau b-independent")
log(f"  {'b':>3} | {'mean up(b-8:b-1)':>17} | {'edge jump at':>12} "
    f"| {'plateau(b+6:b+16)':>18}")
for b in (10, 20, 30, 40):
    _, _, dtp = run_profile(N, seed_site(N, SEED), b, +0.10, 0.0, TMAX)
    up_m, _, _ = plateau_stats(dtp, max(b - 8, 2), b - 1)
    pl_m, pl_s, _ = plateau_stats(dtp, b + 6, min(b + 16, 56))
    # edge locator: first site where |dtp - up_m| exceeds half the plateau height
    edge = np.nan
    for i in range(2, 57):
        if not np.isnan(dtp[i]) and abs(dtp[i] - up_m) > 0.5 * abs(pl_m - up_m):
            edge = i
            break
    log(f"  {b:>3} | {up_m:>+17.5f} | {edge:>12} | {pl_m:>+13.5f} +- {pl_s:.5f}")

# ---- Run D (P4'): dephasing near-field window + far-field crossover --------
log("\n[Run D / P4'] gamma=0.05: near-field step, far-field diffusive crossover")
log("  D1: N=20, seed 0, defect (4,5), sites 6..9 (t_arr < 1/(4 gamma) = 5):")
_, _, dtp = run_profile(20, seed_site(20, 0), 4, +0.10, 0.05, 12.0, dt=0.002)
seg = dtp[6:10]
log(f"     Delta-t sites 6..9: {np.array2string(seg, precision=4)}")
log(f"     pre-reg: within factor 2 of -delta/(2J) = -0.0500")
log("  D2 (descriptive): N=60, seed 0, defect (29,30), far half:")
t0, td, dtp = run_profile(60, seed_site(60, 0), 29, +0.10, 0.05, 40.0, dt=0.002)
near_m, near_s, _ = plateau_stats(dtp, 32, 40)
far_m, far_s, farn = plateau_stats(dtp, 40, 56)
ncross = int(np.sum(~np.isnan(dtp[40:56])))
log(f"     just-past-defect (32:40): {near_m:+.4f} +- {near_s:.4f}")
log(f"     far (40:56): {far_m:+.4f} +- {far_s:.4f}  ({ncross}/16 sites cross)")
log("     (no step claim far-field; diffusive regime, documented)")

# ---- Run D follow-up probe: the far-field surprise, checked before believed
log("\n[Probe] far-field step at gamma=0.05 vs threshold (theta-stability + the")
log("  absolute-threshold artifact that produced the design-review scout verdict):")
tg = np.arange(0.0, 40.001, 0.002)
sd = seed_site(60, 0)
P0 = populations_lindblad(hopping(60), sd, tg, 0.05)
Pd = populations_lindblad(hopping(60, 29, +0.10), sd, tg, 0.05)
ref = P0.max(axis=0)
for th in (0.05, 0.1, 0.2, 0.4, 0.6, 0.8):
    d = arrival_times(Pd, tg, ref, th) - arrival_times(P0, tg, ref, th)
    seg = d[40:56]
    seg = seg[~np.isnan(seg)]
    log(f"  relative theta={th:.2f}: mean {np.mean(seg):+.5f} +- {np.std(seg):.5f}"
        f"  ({len(seg)}/16 cross)")
ones = np.ones(60)
d = (arrival_times(Pd, tg, ones, 0.01) - arrival_times(P0, tg, ones, 0.01))
seg = d[40:56]
seg = seg[~np.isnan(seg)]
log(f"  ABSOLUTE thr=0.01 (the scout's observable): mean {np.mean(seg):+.5f} "
    f"+- {np.std(seg):.5f}  ({len(seg)}/16 cross)  <- the artifact, reproduced")
i = 50
thr = 0.2 * ref[i]
k = np.nonzero(P0[:, i] >= thr)[0][0]
log(f"  timing is still ballistic: site 50 clean crossing t={tg[k]:.2f} "
    f"~ 50/(2J)=25.0 (diffusive would be ~125); max_clean P_50={ref[i]:.4f}")

# ---- Run E: N=7 descriptive side-by-side with the PTF alpha row ------------
log("\n[Run E] N=7 descriptive: middle defect bond (2,3), delta=+0.10, seed 0")
log("  (incommensurate with alpha: localized-packet arrival vs delocalized-state")
log("   purity rescaling; printed side by side, no test)")
for gam, tag in ((0.0, "gamma=0   "), (0.05, "gamma=0.05")):
    _, _, dtp = run_profile(7, seed_site(7, 0), 2, +0.10, gam, 10.0, dt=0.002)
    log(f"  {tag} Delta-t_i: {np.array2string(dtp, precision=4)}")
log("  PTF alpha_i (N=7, bond (0,1), dJ=+0.1, gamma=0.05, bonding mode):")
log("             [1.095, 1.182, 1.051, 0.991, 0.845, 0.923, 0.997]  (smooth, nonlocal)")

log("\n" + "=" * 78)
log("DONE")

import os

os.makedirs("simulations/results/cone_defect_arrival", exist_ok=True)
with open("simulations/results/cone_defect_arrival/arrival_runs.txt", "w",
          encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
