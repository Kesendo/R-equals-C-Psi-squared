"""Front-survival exponent A: does it settle on a clean fraction of the coherence
dose gamma_phi = 4 gamma, or is it a slowly-varying function of (n, gamma)?

This verifier pins the answer to the last open item of
experiments/COUPLING_DEFECT_WALK_TIME_STEP.md (the front-survival A). The prior
script simulations/cone_walk_time_residuals.py measured A ~ 2.8 ~ 0.7 gamma_phi in
the n ~ 15-55 window and left it as an apparent constant. It is not a constant.

  QUESTION. For the first-arrival front peak survival
     g(n, gamma) = max_t P_n(gamma) / max_t P_n(0) ~ exp(-A K),  K = gamma n / (2J),
  does A settle on a clean sub-gamma_phi constant (8/3? 2.8?) or is it structural?

  ANSWER. No sub-gamma_phi constant. A_eff(n, gamma) is the caustic-corrected
  exponent 4 - O(n^{-1/3}) (up in n) - O(gamma) (down in gamma); the only closed
  form is the asymptote A_inf = gamma_phi = 4, approached from below like n^{-1/3}
  (the Airy-caustic exponent). The measured 2.55-3.05 is the pre-asymptotic surface.
  Crucially, most of the apparent gamma-drift is a MEASUREMENT ARTIFACT: the
  standard ratio-of-windowed-maxima method mixes in the front peak-time shift; the
  physical same-arrival-instant survival is nearly gamma-independent.

The exact single-excitation dephasing model:
     rho_dot = -i[h, rho] - Gamma (rho - diag rho),  Gamma = 4 gamma = gamma_phi,
  h nearest-neighbour hopping J = 1 (band 2 cos q, front speed 2J). The clean
  amplitude propagator on the infinite chain is G_{nm}(tau) = (-i)^{|n-m|}
  J_{|n-m|}(2 J tau), so |G|^2 = J_{|n-m|}(2 J tau)^2.

The exact renewal / ladder representation (proven, validated below to grid tolerance):
     P_n(t) = e^{-Gamma t} S_n(t),
     S_n(t) = |G_{n0}(t)|^2 + Gamma int_0^t ds sum_m |G_{nm}(t-s)|^2 S_m(s).
  Momentum-Laplace closed form Shat(p, z) = 1 / (sqrt(z^2 + a^2) - Gamma),
  a(p) = 4 J sin(p/2). j = 0 is the coherent front |<a_n>|^2 = e^{-Gamma t} J_n(2Jt)^2;
  j >= 1 is the incoherent halo that refills near the front.

This script verifies the 2026-07-13 derivation (walk_time_A_derivation.md), which was
independently reviewed by a mathematical referee (the I_1 caustic power-counting and
the renewal exactness) and a physics referee (the same-arrival-instant artifact split).

Run:  python simulations/cone_front_survival_asymptote.py
Writes simulations/results/cone_defect_arrival/front_survival_asymptote.txt
"""

import os
import numpy as np
from scipy.special import jn, jv, airy
from scipy.integrate import quad

J = 1.0
OUT = []


def log(s=""):
    print(s, flush=True)
    OUT.append(s)


# ---- single-excitation engine (same tight-binding / Haken-Strobl RK4 as -------
# ---- cone_walk_time_residuals.py) ---------------------------------------------
def hopping(n):
    h = np.zeros((n, n))
    for a in range(n - 1):
        h[a, a + 1] = h[a + 1, a] = J
    return h


def lindblad_P(n, seed, tgrid, gamma, dt_int=0.02):
    """Direct RK4 Lindblad on the N x N single-excitation density matrix; returns
    P_n(t) = rho_nn(t) on tgrid. This is the truth the renewal solution is tested
    against."""
    h = hopping(n)
    rho = np.zeros((n, n), complex)
    rho[seed, seed] = 1.0
    mask = 4.0 * gamma * (1.0 - np.eye(n))

    def rhs(r):
        return -1j * (h @ r - r @ h) - mask * r

    P = np.empty((len(tgrid), n))
    P[0] = np.real(np.diag(rho))
    t = 0.0
    for k in range(1, len(tgrid)):
        while t < tgrid[k] - 1e-12:
            s = min(dt_int, tgrid[k] - t)
            k1 = rhs(rho); k2 = rhs(rho + 0.5 * s * k1)
            k3 = rhs(rho + 0.5 * s * k2); k4 = rhs(rho + s * k3)
            rho = rho + (s / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += s
        P[k] = np.real(np.diag(rho))
    return P


def renewal_P(dmax, tgrid, gamma, Np=1024):
    """Solve the renewal equation by per-momentum scalar Volterra iteration and
    inverse-transform to displacement d = 0..dmax from an interior seed on the
    infinite chain. Returns P_d(t) = e^{-Gamma t} S_d(t). Uses the Graf identity
    Khat(p, tau) = J_0(4 J tau sin(p/2))."""
    Gamma = 4.0 * gamma
    h = tgrid[1] - tgrid[0]
    p = np.arange(Np) * 2 * np.pi / Np
    a = 4.0 * J * np.abs(np.sin(p / 2))
    nt = len(tgrid)
    K = jn(0, np.outer(tgrid, a))            # Khat(p, tau), shape (nt, Np)
    Shat = np.zeros((nt, Np))
    Shat[0] = K[0]                           # = 1
    for k in range(1, nt):
        conv = 0.5 * K[k] * Shat[0]
        if k > 1:
            conv = conv + np.sum(K[k - 1:0:-1] * Shat[1:k], axis=0)
        Shat[k] = (K[k] + Gamma * h * conv) / (1.0 - 0.5 * Gamma * h)
    d = np.arange(dmax + 1)
    phase = np.exp(1j * np.outer(d, p)) / Np
    S = np.real(Shat @ phase.T)
    return S * np.exp(-Gamma * tgrid)[:, None], Shat


def parabolic_peak(P_col, tgrid, tf, half_lo=2.0, half_hi=4.0):
    """3-point parabolic-refined peak (value and time) of a column P_col in the
    first-arrival window [tf - half_lo, tf + half_hi]."""
    dt = tgrid[1] - tgrid[0]
    lo = max(1, int((tf - half_lo) / dt))
    hi = min(len(tgrid) - 1, int((tf + half_hi) / dt))
    seg = P_col[lo:hi]
    i = lo + int(np.argmax(seg))
    y0, y1, y2 = P_col[i - 1], P_col[i], P_col[i + 1]
    denom = y0 - 2 * y1 + y2
    if denom < 0:
        frac = 0.5 * (y0 - y2) / denom
        return y1 - 0.25 * (y0 - y2) * frac, tgrid[i] + frac * dt
    return y1, tgrid[i]


# ============================================================================
log("=" * 78)
log("FRONT-SURVIVAL EXPONENT A: does it settle on a clean fraction of gamma_phi?")
log("Exact single-excitation dephasing, J = 1, Gamma = 4 gamma = gamma_phi")
log("Answer: no sub-gamma_phi constant; A_inf = gamma_phi = 4, approached ~ n^{-1/3}")
log("=" * 78)

# ---------------------------------------------------------------------------
# SECTION 1 : RENEWAL EXACTNESS (renewal Volterra vs direct RK4 Lindblad)
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[1] RENEWAL EXACTNESS : per-momentum Volterra vs direct RK4 Lindblad")
log("-" * 78)
log("    P_n = e^{-Gamma t}[|G_{n0}|^2 + Gamma int_0^t ds sum_m |G_{nm}(t-s)|^2 S_m(s)]")
log("    direct RK4 on the N=121 density matrix (interior seed) is the truth")

GAMMA1 = 0.05
NCHAIN = 121
SEED1 = 60
tg1 = np.arange(0.0, 22.0 + 1e-9, 0.02)
Pdir = lindblad_P(NCHAIN, SEED1, tg1, GAMMA1)
Pren, Shat1 = renewal_P(40, tg1, GAMMA1)
mask_t = tg1 <= 18.0             # pre-boundary-reflection window

log(f"\n    gamma = {GAMMA1}, compare over t <= 18 (before boundary reflections):")
log(f"    {'d':>4} {'peak_dir':>10} {'peak_ren':>10} {'max|dP|':>10}")
maxdev1 = 0.0
for d in [5, 10, 15, 20, 25, 30]:
    pd = Pdir[:, SEED1 + d][mask_t]
    pr = Pren[:, d][mask_t]
    dev = float(np.max(np.abs(pd - pr)))
    maxdev1 = max(maxdev1, dev)
    log(f"    {d:>4} {pd.max():>10.5f} {pr.max():>10.5f} {dev:>10.2e}")
log(f"\n    max deviation over all d, t<=18 : {maxdev1:.2e}  (grid-limited, dt=0.02)")
TOL1 = 5e-3
assert maxdev1 < TOL1, f"renewal-vs-direct deviation {maxdev1:.2e} exceeds {TOL1:.0e}"
log(f"    PASS: renewal representation reproduces the direct Lindblad below {TOL1:.0e}")

# ---------------------------------------------------------------------------
# SECTION 2 : GREEN'S FUNCTION SANITY (conservation + gamma=0 limit)
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[2] GREEN'S FUNCTION SANITY : probability conservation + gamma=0 = J_n(2Jt)^2")
log("-" * 78)

# 2a. probability conservation from the trace-preserving direct RK4
Ptot = Pdir.sum(axis=1)
cons_dev = float(np.max(np.abs(Ptot - 1.0)))
log(f"\n[2a] sum_n P_n(t) from the direct RK4 (interior seed, pre-reflection):")
log(f"     max_t | sum_n P_n - 1 | = {cons_dev:.2e}")
TOL_CONS = 1e-12
assert cons_dev < TOL_CONS, f"probability not conserved: {cons_dev:.2e}"
log(f"     PASS: trace preserved below {TOL_CONS:.0e}")

# 2b. gamma = 0 renewal reproduces the clean front J_n(2Jt)^2
tg2 = np.arange(0.0, 22.0 + 1e-9, 0.02)
Pclean, _ = renewal_P(40, tg2, 0.0)
log("\n[2b] gamma = 0 renewal S_n(t) vs the clean front J_n(2Jt)^2:")
log(f"     {'d':>4} {'max|S_ren - J^2|':>18}")
maxdev2 = 0.0
for d in [5, 10, 20, 30]:
    jref = jn(d, 2 * J * tg2) ** 2
    dev = float(np.max(np.abs(Pclean[:, d] - jref)))
    maxdev2 = max(maxdev2, dev)
    log(f"     {d:>4} {dev:>18.2e}")
TOL2 = 5e-3
assert maxdev2 < TOL2, f"gamma=0 limit off by {maxdev2:.2e}"
log(f"     max deviation {maxdev2:.2e}; PASS below {TOL2:.0e} (p-grid Np=1024 limited)")

# ---------------------------------------------------------------------------
# SECTION 3 : COHERENT PIECE (j=0 = |<a_n>|^2, A_coh -> 4 from above)
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[3] COHERENT PIECE : j=0 term = e^{-Gamma t} J_n(2Jt)^2 = |<a_n>|^2, A_coh -> 4+")
log("-" * 78)

# 3a. the noise-averaged amplitude damps uniformly at Gamma/2 = 2 gamma, so
#     |<a_n>(t)|^2 = e^{-Gamma t} J_n(2Jt)^2. Verify by evolving the amplitude ODE
#     a_dot = -i h a - (Gamma/2) a directly.
def amplitude_pop(n, seed, tgrid, gamma):
    h = hopping(n)
    a = np.zeros(n, complex); a[seed] = 1.0
    Ph = np.empty((len(tgrid), n)); Ph[0] = np.abs(a) ** 2
    dt_int = 0.02; t = 0.0
    half = 2.0 * gamma  # Gamma/2

    def rhs(v):
        return -1j * (h @ v) - half * v

    for k in range(1, len(tgrid)):
        while t < tgrid[k] - 1e-12:
            s = min(dt_int, tgrid[k] - t)
            k1 = rhs(a); k2 = rhs(a + 0.5 * s * k1)
            k3 = rhs(a + 0.5 * s * k2); k4 = rhs(a + s * k3)
            a = a + (s / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += s
        Ph[k] = np.abs(a) ** 2
    return Ph

tg3 = np.arange(0.0, 22.0 + 1e-9, 0.02)
Pamp = amplitude_pop(NCHAIN, SEED1, tg3, GAMMA1)
coh_dev = 0.0
for d in [10, 20, 30]:
    ref = np.exp(-4 * GAMMA1 * tg3) * jn(d, 2 * J * tg3) ** 2
    coh_dev = max(coh_dev, float(np.max(np.abs(Pamp[:, SEED1 + d] - ref))))
log(f"\n[3a] |<a_n>(t)|^2 (amplitude ODE with uniform decay Gamma/2 = 2 gamma) vs")
log(f"     e^{{-Gamma t}} J_n(2Jt)^2 : max deviation {coh_dev:.2e}")
TOL3 = 5e-4
assert coh_dev < TOL3, f"coherent-front identity off by {coh_dev:.2e}"
log(f"     PASS below {TOL3:.0e}: the j=0 ladder term IS the coherent front |<a_n>|^2")

# 3b. A_coh(n) = 8 J t*/n with 2 J t* = n + 0.8086 n^{1/3} (first Airy maximum)
log("\n[3b] A_coh(n) = 8 J t*/n, first Airy max 2 J t* = n + 0.8086 n^{1/3}, -> 4+")
AIRY1 = 0.8086  # largest maximum of J_n at x = n + 0.8086 n^{1/3}
log(f"     {'n':>4} {'t* numeric':>11} {'t* formula':>11} {'A_coh num':>10} {'A_coh form':>11}")
Acoh_vals = []
for n in [30, 50, 100]:
    x = np.linspace(n - 2.0, n + 6.0 * n ** (1 / 3) + 6.0, 400000)
    v = jv(n, x) ** 2
    xstar = x[int(np.argmax(v))]
    tstar_num = xstar / (2 * J)
    tstar_form = (n + AIRY1 * n ** (1 / 3)) / (2 * J)
    Acoh_num = 8 * J * tstar_num / n
    Acoh_form = 8 * J * tstar_form / n
    Acoh_vals.append(Acoh_num)
    log(f"     {n:>4} {tstar_num:>11.4f} {tstar_form:>11.4f} {Acoh_num:>10.4f} {Acoh_form:>11.4f}")
assert Acoh_vals[0] > Acoh_vals[1] > Acoh_vals[2] > 4.0, "A_coh not monotone -> 4+"
log(f"     A_coh decreases {Acoh_vals[0]:.3f} -> {Acoh_vals[1]:.3f} -> {Acoh_vals[2]:.3f},")
log("     monotone toward 4 from ABOVE (the caustic peak sits inside the cone)")

# ---------------------------------------------------------------------------
# SECTION 4 : I_1 SATURATION AND BULK DOMINANCE
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[4] I_1 SATURATION : the first-order refill I_1(n, t*) -> const; bulk dominates")
log("-" * 78)
log("    I_1(n,t*) = int_0^{t*} ds sum_m |G_{nm}(t*-s)|^2 |G_{m0}(s)|^2, seed at 0")
log("    exact 2D Bessel sum saturates to an n-independent constant ~ 0.277")


def peak_time_val(n):
    x = np.linspace(n - 2.0, n + 8.0 * n ** (1 / 3) + 8.0, 200000)
    v = jv(n, x) ** 2
    i = int(np.argmax(v))
    return x[i] / (2 * J), v[i]     # t*, |G_n|^2 at the peak


def I1_exact(n):
    tstar, Gn2 = peak_time_val(n)
    Ns = 1600
    s = (np.arange(Ns) + 0.5) * tstar / Ns
    ds = tstar / Ns
    m = np.arange(-8, n + 9)
    total = 0.0
    for sv in s:
        leg2 = jv(m, 2 * J * sv) ** 2               # |G_{m0}(s)|^2
        leg1 = jv(n - m, 2 * J * (tstar - sv)) ** 2  # |G_{nm}(t*-s)|^2
        total += np.dot(leg1, leg2)
    return total * ds, Gn2, tstar

log(f"\n    {'n':>5} {'I_1 exact':>11} {'|G_n|^2':>10} {'|G_n|^2 n^{2/3}':>15} "
    f"{'coeff = 8 I_1/(|G_n|^2 n^{2/3})':>32}")
I1_last = None
Gn2coeff_last = None
for n in [15, 30, 55, 90, 120, 200]:
    i1, gn2, tstar = I1_exact(n)
    gcoeff = gn2 * n ** (2 / 3)
    coeff = 8 * i1 / gcoeff
    I1_last = i1
    Gn2coeff_last = gcoeff
    log(f"    {n:>5} {i1:>11.5f} {gn2:>10.6f} {gcoeff:>15.5f} {coeff:>32.5f}")
assert 0.27 < I1_last < 0.28, f"I_1 did not saturate near 0.277 (got {I1_last:.4f})"
log(f"\n    I_1 saturates to {I1_last:.5f} (n-independent). |G_n|^2 n^{{2/3}} -> 0.4553")
log(f"    => coefficient of n^{{-1/3}} in (ln B)/K is 8 I_1 / 0.4553 = "
    f"{8 * I1_last / 0.4553:.4f}")

# 4b. the joint-caustic-sliver-only version (uniform Airy, both legs on the cone)
log("\n[4b] joint-caustic-sliver-only value (reduced uniform-Airy double integral):")
log("     C = 2^{1/3} int_0^1 dsigma sigma^{-1/3}(1-sigma)^{-2/3} Phi(sigma),")
log("     Phi(sigma) = int dv Ai(-2^{1/3} v)^2 Ai(2^{1/3} v (sigma/(1-sigma))^{1/3})^2")


def Ai(x):
    return airy(x)[0]


def Phi(sigma):
    r = (sigma / (1 - sigma)) ** (1 / 3)
    f = lambda v: Ai(-2 ** (1 / 3) * v) ** 2 * Ai(2 ** (1 / 3) * v * r) ** 2
    val, _ = quad(f, -80, 80, limit=1200)
    return val

Cint, _ = quad(lambda s: s ** (-1 / 3) * (1 - s) ** (-2 / 3) * Phi(s),
               1e-7, 1 - 1e-7, limit=600)
C_sliver = 2 ** (1 / 3) * Cint
ratio = C_sliver / I1_last
log(f"     caustic-sliver C = {C_sliver:.5f}   vs   exact I_1 = {I1_last:.5f}")
log(f"     ratio C / I_1 = {ratio:.3f}  (~30%)")
assert 0.07 < C_sliver < 0.10, f"caustic sliver off (got {C_sliver:.4f})"
assert 0.20 < ratio < 0.40, f"sliver/exact ratio off (got {ratio:.3f})"
log("     => the joint caustic fixes only the EXPONENT n^{-1/3}; the bulk cone")
log("        interior supplies ~70% of the constant. Power-counting, not the value.")

# ---------------------------------------------------------------------------
# SECTION 5 : THE ARTIFACT SPLIT (same-arrival-instant vs ratio-of-maxima)
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[5] THE ARTIFACT SPLIT : same-arrival-instant A vs windowed ratio-of-maxima A")
log("-" * 78)
log("    interior seed n=40, gamma sweep. (i) same-instant: numerator and")
log("    denominator BOTH at the CLEAN caustic time t*_0. (ii) max-ratio: the")
log("    doc method, ratio of front-windowed maxima (front window declared).")
log("    The gamma-drift of (ii) is mostly the front peak-time shift, not physics.")

NI = 161
SI = 80
N5 = 40
tg5 = np.arange(0.0, N5 / (2 * J) + 10.0 + 1e-9, 0.01)
site5 = SI + N5
tf5 = N5 / (2 * J)
# clean run: clean caustic peak value and time
Pc5 = lindblad_P(NI, SI, tg5, 0.0)
clean_peak5, tstar0 = parabolic_peak(Pc5[:, site5], tg5, tf5)
# front window: [0.5 t*_0, t*_0 + 4 n^{1/3}/(2J)] (excludes the late diffusive bump)
win_lo = 0.5 * tstar0
win_hi = tstar0 + 4.0 * (N5 ** (1 / 3)) / (2 * J)
log(f"\n    clean caustic: t*_0 = {tstar0:.4f}, peak = {clean_peak5:.6f}")
log(f"    front window declared: [{win_lo:.3f}, {win_hi:.3f}]")
log("    'edge' flags when the windowed max sits at the window upper edge: the")
log("    front peak has drifted past and the diffusive bump is entering (breakdown).")
log(f"    {'gamma':>7} {'K':>7} {'A_sameinstant':>14} {'A_maxratio':>11} "
    f"{'front dt':>9} {'edge':>5}")
same_col = []
maxr_col = []
edge_col = []
GAMMAS5 = [0.002, 0.01, 0.03, 0.05, 0.08]
for gam in GAMMAS5:
    Pg5 = lindblad_P(NI, SI, tg5, gam)
    col = Pg5[:, site5]
    K = gam * N5 / (2 * J)
    # (i) same arrival instant: interpolate P(gamma) at the CLEAN t*_0
    idx = tstar0 / (tg5[1] - tg5[0])
    i0 = int(np.floor(idx)); fr = idx - i0
    P_at = col[i0] * (1 - fr) + col[i0 + 1] * fr
    A_same = -np.log(P_at / clean_peak5) / K
    # (ii) windowed ratio of maxima
    dt5 = tg5[1] - tg5[0]
    lo = max(1, int(win_lo / dt5)); hi = min(len(tg5) - 1, int(win_hi / dt5))
    seg = col[lo:hi]
    jj = int(np.argmax(seg))
    j = lo + jj
    edge = jj >= len(seg) - 2   # windowed max at the upper edge = breakdown
    y0, y1, y2 = col[j - 1], col[j], col[j + 1]
    den = y0 - 2 * y1 + y2
    if den < 0:
        frac = 0.5 * (y0 - y2) / den
        pw = y1 - 0.25 * (y0 - y2) * frac
        tw = tg5[j] + frac * dt5
    else:
        pw = y1; tw = tg5[j]
    A_max = -np.log(pw / clean_peak5) / K
    same_col.append(A_same); maxr_col.append(A_max); edge_col.append(edge)
    log(f"    {gam:>7.3f} {K:>7.3f} {A_same:>14.4f} {A_max:>11.4f} {tw - tstar0:>+9.4f} "
        f"{'EDGE' if edge else '':>5}")

# same-instant is robust across the whole range; report its full drift.
drift_same = same_col[-1] - same_col[0]
# max-ratio: report the drift to gamma=0.05 (index 3, the derivation's figure), where
# the windowed peak is still a genuine front feature; beyond that the method breaks.
i05 = GAMMAS5.index(0.05)
drift_max_05 = maxr_col[i05] - maxr_col[0]
drift_max_full = maxr_col[-1] - maxr_col[0]
log(f"\n    same-instant A drift, full range gamma 0.002 -> 0.08 : {drift_same:+.4f}")
log(f"    max-ratio   A drift, gamma 0.002 -> 0.05 (robust)    : {drift_max_05:+.4f}")
log(f"    max-ratio   A drift, full range 0.002 -> 0.08        : {drift_max_full:+.4f}")
log("    (the front peak-time shift front dt jumps once dephasing kills the coherent")
log("     front; at gamma >= 0.05 the windowed max is the diffusive bump, not the front)")
assert abs(drift_same) < 0.30, f"same-instant drift too large: {drift_same:+.4f}"
assert abs(drift_max_05) > 2.5 * abs(drift_same), \
    f"max-ratio drift not dominated by the artifact ({drift_max_05:+.4f} vs {drift_same:+.4f})"
log("    => the PHYSICAL (same-instant) A is nearly gamma-independent (~0.14); the")
log("       measured ~0.55 gamma-drift is mostly the front peak-time-shift artifact.")

# ---------------------------------------------------------------------------
# SECTION 6 : THE n-CLIMB (A_eff climbs toward 4 like n^{-1/3})
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[6] THE n-CLIMB : A_eff(n) at small gamma climbs toward A_inf = 4 like n^{-1/3}")
log("-" * 78)
log("    exact renewal solution (interior seed, infinite chain), gamma = 0.002")
log("    (first-order regime). 4 - A_eff should fall as n^{-p}, p ~ 1/3 with a")
log("    positive n^{-2/3} sub-leading term shallowing the apparent slope in-window.")

tg6 = np.arange(0.0, 42.0 + 1e-9, 0.03)
NS6 = [15, 30, 55]
P0_6, _ = renewal_P(60, tg6, 0.0)
GAM6 = 0.002
Pg6, _ = renewal_P(60, tg6, GAM6)
log(f"\n    {'n':>4} {'K':>8} {'g':>9} {'-ln g':>9} {'A_eff':>9} {'4 - A_eff':>10}")
A6 = []
for n in NS6:
    mp0, _ = parabolic_peak(P0_6[:, n], tg6, n / (2 * J))
    pk, _ = parabolic_peak(Pg6[:, n], tg6, n / (2 * J))
    g = pk / mp0
    K = GAM6 * n / (2 * J)
    A = -np.log(g) / K
    A6.append(A)
    log(f"    {n:>4} {K:>8.4f} {g:>9.5f} {-np.log(g):>9.5f} {A:>9.4f} {4 - A:>10.4f}")
assert A6[0] < A6[1] < A6[2] < 4.0, "A_eff not a monotone climb toward 4"
# fit apparent exponent p on 4 - A_eff = c n^{-p}
lp = np.polyfit(np.log(NS6), np.log([4 - a for a in A6]), 1)
p_fit = -lp[0]
c_fit = np.exp(lp[1])
log(f"\n    A_eff climbs {A6[0]:.4f} -> {A6[1]:.4f} -> {A6[2]:.4f} (toward 4)")
log(f"    fit 4 - A_eff = {c_fit:.3f} n^{{-{p_fit:.3f}}}  (apparent p ~ 0.26 in-window,")
log("    the pure Airy exponent 1/3 shallowed by the positive n^{-2/3} sub-leading)")
log(f"    asymptotic n^{{-1/3}} coefficient (from I_1): 8 I_1 / 0.4553 = "
    f"{8 * I1_last / 0.4553:.4f}")
assert 0.15 < p_fit < 0.45, f"apparent exponent out of range: {p_fit:.3f}"

log("\n" + "=" * 78)
log("SUMMARY: no sub-gamma_phi constant. A_inf = gamma_phi = 4, approached ~ n^{-1/3}.")
log("The single clean number is the asymptote; 2.55-3.05 is the pre-asymptotic surface,")
log("and most of the measured gamma-drift is the max-ratio peak-shift artifact.")
log("=" * 78)
log("DONE")

os.makedirs("simulations/results/cone_defect_arrival", exist_ok=True)
with open("simulations/results/cone_defect_arrival/front_survival_asymptote.txt", "w",
          encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
