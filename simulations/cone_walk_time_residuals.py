"""Walk-time step, the two named-open items resolved (follow-up to
cone_defect_arrival.py; see experiments/COUPLING_DEFECT_WALK_TIME_STEP.md).

The original experiment landed the first-order step law -delta/(2J) and named two
residuals open. This script pins both, reusing the same single-excitation engine
(the primitives below are the same tight-binding hopping / Haken-Strobl RK4 as
cone_defect_arrival.py).

  ITEM 1 (the +- second-order packet residual). The exact single-q front delay is
     Delta_t_sq(delta) = -(u^2 - 1) / (2J (u^2 + 1)),  u = J'/J = 1 + delta,
     with Taylor -delta/(2J) + delta^2/(4J) - 0*delta^3, i.e. an intrinsic second
     order c2 = +1/(4J) = +0.25 and NO delta^3 asymmetry. Section 1 verifies this
     closed form against the numeric Wigner delay, then shows the large measured
     residual (clean-ref c2 ~ +1.18, a strong +- asymmetry) is an artifact of the
     relative-threshold arrival method: the transmitted front is dimmed by the
     reflected fraction |r|^2 ~ delta^2, so it crosses theta * max_t P|clean LATER
     (an even O(delta^2) delay). Referencing the threshold to the defect run's own
     peak removes it and collapses c2 to single-q order.

  ITEM 2 (front amplitude survival vs gamma). Section 2 checks the coherent (gamma=0)
     front peak against the Airy/Bessel caustic scaling max_t P_n ~ C n^(-2/3), then
     measures the dephasing survival g = max_t P_n(gamma) / max_t P_n(0) in the
     first-arrival window. It falls off exponentially in the dose K = gamma * t_arr =
     gamma * n / (2J) as g ~ exp(-A K), A ~ 2.8 ~ 0.7 * gamma_phi (gamma_phi = 4 gamma),
     collapsing onto K only approximately (a ~15-20% systematic residual). The naive
     full-rate corollary g = exp(-4 gamma t_arr) = exp(-4K), A = 4, is falsified: the
     front survives several times better. Section 2c pins the mechanism: the COHERENT
     front (the noise-averaged amplitude |<a_n>|^2 = e^(-4 gamma t) P_n^coh) pays the
     full 4 gamma rate; the survival boost is incoherent population that dephasing
     redistributes toward the front, not the caustic decaying more slowly.

     The measured A ~ 2.8 is a pre-asymptotic window value, not a constant of the model:
     the sequel verifier cone_front_survival_asymptote.py (sections [7]/[8]) closes the
     n^(-1/3) coefficient in Airy form and derives the true fixed-gamma ceiling
     A_inf(gamma) = 4 - phi(2J)/gamma < 4 (the experiment doc's third follow-up).

Run:  python simulations/cone_walk_time_residuals.py
Writes simulations/results/cone_defect_arrival/walk_time_residuals.txt
"""

import os
import numpy as np
from scipy.special import jn

J = 1.0
OUT = []


def log(s=""):
    print(s, flush=True)
    OUT.append(s)


# ---- single-excitation engine (identical to cone_defect_arrival.py) --------
def hopping(n, defect_bond=None, delta=0.0):
    h = np.zeros((n, n))
    for a in range(n - 1):
        j = J * (1.0 + (delta if defect_bond == a else 0.0))
        h[a, a + 1] = j
        h[a + 1, a] = j
    return h


def populations_pure(h, psi0, tgrid):
    E, U = np.linalg.eigh(h)
    c = U.T @ psi0
    phases = np.exp(-1j * np.outer(tgrid, E)) * c
    amps = phases @ U.T
    return np.abs(amps) ** 2


def populations_lindblad(h, psi0, tgrid, gamma, dt_int=0.004):
    n = h.shape[0]
    rho = np.outer(psi0, psi0.conj()).astype(complex)
    mask = 4.0 * gamma * (1.0 - np.eye(n))

    def rhs(r):
        return -1j * (h @ r - r @ h) - mask * r

    P = np.empty((len(tgrid), n))
    t = 0.0
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
        P[k] = np.real(np.diag(rho))
    return P


def seed_site(n, s):
    v = np.zeros(n, dtype=complex)
    v[s] = 1.0
    return v


def seed_packet(n, n0, q=-np.pi / 2, sigma_q=0.1 * np.pi):
    sigma_x = 1.0 / (2.0 * sigma_q)
    x = np.arange(n)
    v = np.exp(-((x - n0) ** 2) / (4.0 * sigma_x ** 2)) * np.exp(1j * q * x)
    return v / np.linalg.norm(v)


def arrival_times(P, tgrid, ref_peak, theta):
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


def even_odd_coeffs(deltas, getval):
    """Split the arrival excess over first order into c2*delta^2 (even) and
    c3*delta^3 (odd); returns (c2, c3)."""
    ev, od = [], []
    for delta in deltas:
        ex_p = getval(+delta) - (-delta / 2)
        ex_m = getval(-delta) - (+delta / 2)
        ev.append(0.5 * (ex_p + ex_m))
        od.append(0.5 * (ex_p - ex_m))
    dd = np.array(deltas); ev = np.array(ev); od = np.array(od)
    c2 = np.sum(dd ** 2 * ev) / np.sum(dd ** 4)
    c3 = np.sum(dd ** 3 * od) / np.sum(dd ** 6)
    return c2, c3


# ============================================================================
log("=" * 78)
log("WALK-TIME STEP: the two named-open residuals, pinned")
log("J = 1, single-excitation sector")
log("=" * 78)

# ---------------------------------------------------------------------------
# SECTION 1 : ITEM 1, the second-order packet residual
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[1] ITEM 1 : the +- second-order arrival residual")
log("-" * 78)

# 1a. the exact single-q front delay vs the closed form (independent of the packet)
log("\n[1a] exact single-q front delay: numeric Wigner delay vs the closed form")
log("     Delta_t_sq = -(u^2-1)/(2J(u^2+1)),  u = 1+delta   [Taylor c2=+0.25, c3=0]")
q0 = -np.pi / 2


def phi_tau(q, u):
    tau = -2j * u * np.sin(q) / (np.exp(-1j * q) - u ** 2 * np.exp(1j * q))
    return np.angle(tau)


def wigner_dt(u):
    h = 1e-5
    dphidq = (phi_tau(q0 + h, u) - phi_tau(q0 - h, u)) / (2 * h)
    vg = 2 * J * abs(np.sin(q0))
    return -(-dphidq) / vg


log(f"     {'delta':>6} {'numeric':>10} {'closed form':>12} {'excess/delta^2':>15}")
for d in (0.02, 0.05, 0.10, 0.15, -0.05, -0.10, -0.15):
    u = 1 + d
    num = wigner_dt(u)
    closed = -(u ** 2 - 1) / (2 * J * (u ** 2 + 1))
    log(f"     {d:>+6.2f} {num:>+10.6f} {closed:>+12.6f} {(closed + d / 2) / d ** 2:>+15.4f}")

# 1b. the packet measurement: is the large residual a relative-threshold artifact?
log("\n[1b] packet measurement (q=-pi/2, sigma_q=0.1*pi, N=120, defect bond (59,60)):")
log("     compare arrival with threshold relative to the CLEAN peak (the paper's")
log("     method) vs the DEFECT run's own peak (removes the amplitude-dimming delay)")
N1, B1, N0, SIG = 120, 59, 10, 0.10 * np.pi
LO, HI = 75, 110
DELTAS = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15]
tg = np.arange(0.0, 55.0 + 1e-9, 0.001)
pk = seed_packet(N1, N0, sigma_q=SIG)
Pc = populations_pure(hopping(N1), pk, tg)
refc = Pc.max(axis=0)
t0 = arrival_times(Pc, tg, refc, 0.2)


def plateau(vec):
    seg = vec[LO:HI]
    return float(np.mean(seg[~np.isnan(seg)]))


d_cref, d_dref, dimming = {}, {}, []
log(f"     {'delta':>6} {'peakRatio':>10} {'dt(clean-ref)':>14} {'dt(defect-ref)':>15}")
for delta in DELTAS + [-d for d in DELTAS]:
    Pd = populations_pure(hopping(N1, B1, delta), pk, tg)
    refd = Pd.max(axis=0)
    ratio = float(np.mean(refd[LO:HI] / refc[LO:HI]))
    d_cref[delta] = plateau(arrival_times(Pd, tg, refc, 0.2) - t0)
    d_dref[delta] = plateau(arrival_times(Pd, tg, refd, 0.2) - t0)
    log(f"     {delta:>+6.2f} {ratio:>10.5f} {d_cref[delta]:>+14.5f} {d_dref[delta]:>+15.5f}")
    if delta > 0:
        dimming.append((1 - 0.5 * (ratio + float(np.mean(
            populations_pure(hopping(N1, B1, -delta), pk, tg).max(axis=0)[LO:HI]
            / refc[LO:HI])))) / delta ** 2)

c2c, c3c = even_odd_coeffs(DELTAS, lambda d: d_cref[d])
c2d, c3d = even_odd_coeffs(DELTAS, lambda d: d_dref[d])
log(f"\n     peak dimming 1-R ~ a*delta^2, a = {np.mean(dimming):.3f}  "
    f"(reflected |r|^2 = delta^2 predicts a ~ 1)")
log(f"     clean-ref : c2 = {c2c:+.4f}   c3 = {c3c:+.4f}")
log(f"     defect-ref: c2 = {c2d:+.4f}   c3 = {c3d:+.4f}   [single-q: c2=+0.25, c3=0]")
log(f"     => the amplitude artifact accounts for {c2c - c2d:+.4f} of the even residual")
log("     => the residual mechanism 'left unpinned' is the relative-threshold method")

# ---------------------------------------------------------------------------
# SECTION 2 : ITEM 2, front amplitude survival vs gamma
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[2] ITEM 2 : front amplitude survival vs gamma")
log("-" * 78)

N = 130
SITES = [15, 20, 25, 30, 35, 40, 45, 50, 55]
GAMMAS = [0.01, 0.02, 0.03, 0.05, 0.08]
TMAX = 40.0
seed = seed_site(N, 0)
tg2 = np.arange(0.0, TMAX + 1e-9, 0.02)


def windowed_peak(P, n):
    """Front peak inside the first-arrival window [t_front-2, t_front+4]; also
    flags when the max sits at the upper edge (the diffusive-onset peak has
    drifted past the window, so the value is not a genuine interior maximum)."""
    tfront = n / (2 * J)
    lo = max(0, int((tfront - 2.0) / 0.02))
    hi = min(len(tg2), int((tfront + 4.0) / 0.02))
    seg = P[lo:hi, n]
    kmax = int(np.argmax(seg))
    return seg.max(), (kmax >= len(seg) - 2)


# 2a. coherent Airy scaling
log("\n[2a] coherent (gamma=0) front peak vs Airy/Bessel caustic max_t J_n(2Jt)^2")
P0 = populations_pure(hopping(N), seed, tg2)
mp0 = {n: windowed_peak(P0, n)[0] for n in SITES}
log(f"     {'n':>3} {'max_t P_n (boundary seed)':>26} {'max_t J_n^2 (interior)':>22} "
    f"{'P_n * n^(2/3)':>14}")
for n in SITES:
    jref = (jn(n, 2 * J * tg2) ** 2).max()
    log(f"     {n:>3} {mp0[n]:>26.6f} {jref:>22.6f} {mp0[n] * n ** (2 / 3):>14.4f}")
log("     (boundary seed = the image-sum-enhanced front, coefficient ~1.8;")
log("      interior Bessel coefficient ~0.46; both scale as n^(-2/3) asymptotically)")

# 2b. dephasing survival collapse onto the dose K = gamma*n/(2J)
log("\n[2b] survival g = max_t P_n(gamma)/max_t P_n(0) vs dose K = gamma*n/(2J)")
log("     naive corollary g = exp(-4 gamma t_arr) = exp(-4K) would give A = 4")
log("     window-edge hits (peak drifted past the window) are flagged and excluded")
Ks, mlng = [], []
log(f"     {'n':>3} {'gamma':>6} {'K':>7} {'g':>8} {'-ln g':>8} {'(-ln g)/K':>10} {'edge':>6}")
for gam in GAMMAS:
    Pg = populations_lindblad(hopping(N), seed, tg2, gam)
    for n in SITES:
        pk, edge = windowed_peak(Pg, n)
        g = pk / mp0[n]
        K = gam * n / (2 * J)
        if g > 0 and not edge:
            Ks.append(K); mlng.append(-np.log(g))
        log(f"     {n:>3} {gam:>6.2f} {K:>7.3f} {g:>8.4f} {-np.log(g):>8.4f} "
            f"{-np.log(g) / K:>10.4f} {'EDGE' if edge else '':>6}")
Ks = np.array(Ks); mlng = np.array(mlng)
p, la = np.polyfit(np.log(Ks), np.log(mlng), 1)
Alin = np.sum(Ks * mlng) / np.sum(Ks * Ks)
small = Ks < 0.4
A_s = np.sum(Ks[small] * mlng[small]) / np.sum(Ks[small] ** 2)
spread = mlng / Ks
log(f"\n     (edge-hit points excluded) power fit  -ln g = {np.exp(la):.3f} * K^{p:.3f}")
log(f"     linear fit -ln g = {Alin:.3f} * K   (small-dose K<0.4 slope {A_s:.3f})")
log(f"     collapse onto K is APPROXIMATE: (-ln g)/K spreads {spread.min():.2f}-{spread.max():.2f}")
log(f"     => A ~ {Alin:.2f} ({Alin/4:.2f}*gamma_phi) to {A_s:.2f} small-dose ({A_s/4:.2f}*gamma_phi),")
log(f"        ~0.7 (between 2/3 and 3/4) of the coherence dose; naive full-rate A=4 falsified.")

# 2c. the mechanism: the COHERENT front pays the full rate; the boost is incoherent
log("\n[2c] mechanism (gamma=0.05): the noise-averaged amplitude <a_n> damps uniformly")
log("     at Gamma/2 = 2 gamma, so |<a_n>|^2 = e^(-4 gamma t) * P_n^coh (the coherent")
log("     front). Compare its survival g_coh to the full survival g_full:")
gam = 0.05
Pfull = populations_lindblad(hopping(N), seed, tg2, gam)
Pcohd = P0 * np.exp(-4 * gam * tg2)[:, None]
log(f"     {'n':>3} {'K':>6} {'naive e^-4K':>12} {'g_coh':>9} {'g_full':>9} {'full/coh':>9}")
for n in [15, 30, 50, 55]:
    gcoh = windowed_peak(Pcohd, n)[0] / mp0[n]
    gfull = windowed_peak(Pfull, n)[0] / mp0[n]
    K = gam * n / (2 * J)
    log(f"     {n:>3} {K:>6.3f} {np.exp(-4 * K):>12.5f} {gcoh:>9.5f} {gfull:>9.5f} "
        f"{gfull / gcoh:>9.2f}")
log("     => g_coh ~ naive (the coherent front pays the FULL 4 gamma rate, even a hair")
log("        more). The front survives several times better because incoherent population,")
log("        redistributed by dephasing, refills near the front; the reduced effective A")
log("        is a property of the total arrival-window population, not a slower caustic.")

log("\n" + "=" * 78)
log("DONE")

os.makedirs("simulations/results/cone_defect_arrival", exist_ok=True)
with open("simulations/results/cone_defect_arrival/walk_time_residuals.txt", "w",
          encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
