"""Front-survival exponent A: does it settle on a clean fraction of the coherence
dose gamma_phi = 4 gamma, or is it a slowly-varying function of (n, gamma)?

This verifier pins the answer to the last open item of
experiments/COUPLING_DEFECT_WALK_TIME_STEP.md (the front-survival A). The prior
script simulations/cone_walk_time_residuals.py measured A ~ 2.8 ~ 0.7 gamma_phi in
the n ~ 15-55 window and left it as an apparent constant. It is not a constant.

  QUESTION. For the first-arrival front peak survival
     g(n, gamma) = max_t P_n(gamma) / max_t P_n(0) ~ exp(-A K),  K = gamma n / (2J),
  does A settle on a clean sub-gamma_phi constant (8/3? 2.8?) or is it structural?

  ANSWER (corrected 2026-07-13, refuting the prior "A_inf = gamma_phi = 4" reading).
  No sub-gamma_phi constant, and NO single clean asymptote 4 either. Two regimes:

    * Pre-asymptotic window n << n* ~ 6 (gamma/J)^{-3/2}. A_eff climbs toward 4
      like A_eff ~ 4 - 4.864 n^{-1/3} (the Airy-caustic exponent). The n^{-1/3}
      coefficient is 8 I_1 / (2^{2/3} Ai(-alpha)^2) = 4.864 with the closed-form
      single-refill constant I_1 = 1/12 + (1/4) int_0^{2c} Ai(-w) dw = 0.27694424,
      2c = 2^{2/3} alpha, alpha = |first zero of Ai'| = 1.0187929716.

    * The TRUE fixed-gamma ceiling is strictly below 4:
        A_inf(gamma) = 4 - phi(2J)/gamma,
        phi(2J) = sqrt(Gamma (Gamma + 4J)) - 4J arcsinh sqrt(Gamma/(4J)),  Gamma = 4 gamma.
      Small gamma: A_inf = 4 - (8/3) sqrt(gamma/J). Values 3.881 (0.002), 3.412
      (0.05), 3.181 (0.10), 2.872 (0.20). The single-refill argument that gave
      "A -> 4" resummed the wrong way: ln S_n(t*) grows LINEARLY in n at rate
      phi(2J)/2 (a large-deviation rate, the dominant pole of Shat), not O(ln n).
      The incoherent halo is an exp-in-n boost e^{(phi/2) n} of the front, lowering
      A by the constant phi(2J)/gamma. A_eff -> 4 holds only as gamma -> 0.

  The gamma-drift of the K-collapse measured in cone_walk_time_residuals.py remains
  mostly a MEASUREMENT ARTIFACT of the ratio-of-windowed-maxima method (the front
  peak-time shift); the same-arrival-instant survival is nearly gamma-independent.

  Both the closed form for I_1 and the fixed-gamma refutation passed independent
  math and physics referee rounds on 2026-07-13 (both recomputed from scratch).

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

This script verifies the 2026-07-13 renewal derivation (session work), which was
independently reviewed by a mathematical referee (the I_1 caustic power-counting and
the renewal exactness) and a physics referee (the same-arrival-instant artifact split).

The gamma ~ J fence is now docked (section [10]): at experimental distances the edge is the
dose K_deg ~ 0.6 (theta = 0.2), gamma_deg ~ 1.2 J/n; gamma ~ J is only the n -> 1 EP corner.
The o(n) is closed by section [9] (referee-confirmed 2026-07-13: independent
SymPy + steepest-descent + direct-Lindblad + Haken-Strobl cross-checks). Section [9] pins
the prefactor S_n(t*_0) = C(gamma) n^{-1/2} e^{(phi/2J)n} with the closed constant
C(gamma) = (2 pi)^{-1/2}(gamma/(gamma+J))^{1/4}, the approach law A_eff = A_inf +
(2J/(gamma n))[-(1/6) ln n + ln(c_G/C)], and the third reading A_max = 0 (the Haken-Strobl
diffusive peak-value exponent), so all three readings of "front survival" are settled.

Run:  python simulations/cone_front_survival_asymptote.py
Writes simulations/results/cone_defect_arrival/front_survival_asymptote.txt
"""

import os
import numpy as np
from scipy.special import jn, jv, ive, airy, ai_zeros
from scipy.integrate import quad

J = 1.0
OUT = []

# ---- shared closed-form constants (the caustic peak and the single-refill I_1) ---
ALPHA = float(abs(ai_zeros(1)[1][0]))           # |first zero of Ai'| = 1.0187929716
TWO_C = 2.0 ** (2.0 / 3.0) * ALPHA              # 2c = 2^{2/3} alpha = 1.6172330
CAUSTIC_PEAK = 2.0 ** (2.0 / 3.0) * airy(-ALPHA)[0] ** 2   # |G_n|^2 n^{2/3} -> 0.45547
_ANTIDER, _ = quad(lambda w: airy(-w)[0], 0.0, TWO_C, limit=400)
I1_CLOSED = 1.0 / 12.0 + 0.25 * _ANTIDER        # single-refill constant = 0.27694424
COEFF_N13 = 8.0 * I1_CLOSED / CAUSTIC_PEAK       # n^{-1/3} coefficient of (ln B)/K = 4.864


def phi_of_gamma(gamma):
    """Large-deviation rate phi(2J) at the front v = 2J (Legendre of the tilted
    pole mu(theta) = sqrt(Gamma^2 + 16 J^2 sinh^2(theta/2)), Gamma = 4 gamma).
    The fixed-gamma ceiling is A_inf(gamma) = 4 - phi(2J)/gamma."""
    G = 4.0 * gamma
    return np.sqrt(G * (G + 4.0 * J)) - 4.0 * J * np.arcsinh(np.sqrt(G / (4.0 * J)))


def C_of_gamma(gamma):
    """Closed prefactor constant of the halo-dominated front (section [9]):
    S_n(t*_0) ~ C(gamma) n^{-1/2} e^{(phi/2J) n}. The single interior Gaussian saddle
    gives C(gamma) = (Gamma/mu*) sqrt(J/(pi mu''(theta*))) = (2 pi)^{-1/2}
    (gamma/(gamma+J))^{1/4}, with mu''(theta*) = 2J sqrt(Gamma/(Gamma+4J))."""
    return (2.0 * np.pi) ** -0.5 * (gamma / (gamma + J)) ** 0.25


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
log("Answer: NO. Pre-asymptotic climb 4 - 4.864 n^{-1/3}; true fixed-gamma ceiling")
log("A_inf(gamma) = 4 - phi(2J)/gamma < 4 (= 4 - (8/3)sqrt(gamma/J) for small gamma).")
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
log(f"\n    I_1 saturates to {I1_last:.5f} (n-independent). |G_n|^2 n^{{2/3}} -> "
    f"{CAUSTIC_PEAK:.5f} = 2^{{2/3}} Ai(-alpha)^2 (sharpened caustic peak)")
log(f"    => coefficient of n^{{-1/3}} in (ln B)/K is 8 I_1 / {CAUSTIC_PEAK:.5f} = "
    f"{COEFF_N13:.4f} (using the closed-form I_1 = {I1_CLOSED:.8f}, section [7])")

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
log("[6] THE n-CLIMB : A_eff(n) at small gamma climbs like n^{-1/3} toward the")
log("    gamma=0.002 ceiling A_inf(0.002) = 3.881 (not resolved from 4 in this window)")
log("-" * 78)
log("    exact renewal solution (interior seed, infinite chain), gamma = 0.002")
log("    (first-order regime). A_inf(0.002) - A_eff should fall as n^{-p}, p ~ 1/3")
log("    with a positive n^{-2/3} sub-leading shallowing the apparent slope in-window.")
log("    NOTE: at these n (<= 55) the ceiling 3.881 is not resolved from 4; the climb")
log("    heads to 3.881, not to 4. The crossover n* ~ 67000 (section [8c]) is far above.")

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
A_INF_6 = 4.0 - phi_of_gamma(GAM6) / GAM6      # fixed-gamma ceiling at gamma=0.002
assert A6[0] < A6[1] < A6[2] < 4.0, "A_eff not a monotone climb toward the ceiling"
# fit apparent exponent p on 4 - A_eff = c n^{-p} (4 stands in for 3.881 in-window)
lp = np.polyfit(np.log(NS6), np.log([4 - a for a in A6]), 1)
p_fit = -lp[0]
c_fit = np.exp(lp[1])
log(f"\n    A_eff climbs {A6[0]:.4f} -> {A6[1]:.4f} -> {A6[2]:.4f} (toward the")
log(f"    gamma=0.002 ceiling A_inf = {A_INF_6:.4f}, indistinguishable from 4 here)")
log(f"    fit 4 - A_eff = {c_fit:.3f} n^{{-{p_fit:.3f}}}  (apparent p ~ 0.26 in-window,")
log("    the pure Airy exponent 1/3 shallowed by the positive n^{-2/3} sub-leading)")
log(f"    asymptotic n^{{-1/3}} coefficient (from closed-form I_1): 8 I_1 / "
    f"{CAUSTIC_PEAK:.5f} = {COEFF_N13:.4f}")
assert 0.15 < p_fit < 0.45, f"apparent exponent out of range: {p_fit:.3f}"

# ---------------------------------------------------------------------------
# SECTION 7 : THE CLOSED FORM OF THE SINGLE-REFILL CONSTANT I_1
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[7] CLOSED FORM OF I_1 : I_1 = 1/12 + (1/4) int_0^{2c} Ai(-w) dw = 0.27694424")
log("-" * 78)
log("    2c = 2^{2/3} alpha, alpha = |first zero of Ai'|. The n^{-1/3} coefficient of")
log("    (ln B)/K is 8 I_1 / (2^{2/3} Ai(-alpha)^2). Verified vs the exact 1D momentum")
log("    representation I_1 = (1/2pi) int_0^{pi/2} cos(2 n th) sin(4 J t* sin th)/sin th dth.")

log(f"\n    alpha = |first zero of Ai'|   = {ALPHA:.10f}")
log(f"    2c    = 2^{{2/3}} alpha          = {TWO_C:.10f}")
log(f"    int_0^{{2c}} Ai(-w) dw         = {_ANTIDER:.10f}")
log(f"    I_1 (closed) = 1/12 + (1/4) int = {I1_CLOSED:.8f}   (target 0.27694424)")
assert abs(I1_CLOSED - 0.27694424) < 1e-6, f"I_1 closed form off: {I1_CLOSED:.8f}"

# endpoint / caustic split (E = 1/8 exact, caustic = I_1 - 1/8)
E_END = 1.0 / 8.0
C_CAUSTIC = I1_CLOSED - E_END
log(f"\n    endpoint/caustic split: E = 1/8 = {E_END:.5f} (exact, no stationary point),")
log(f"    caustic = I_1 - 1/8 = {C_CAUSTIC:.5f}  (referee numeric C_term quad gave 0.15189)")


# exact 1D momentum representation and the numeric peak time t*(n)
def I1_1D(n, t):
    """Exact 1D momentum-Laplace form of I_1(n, t) (Graf identity + the closed
    convolution int_0^t J0(a(t-s)) J0(as) ds = sin(at)/a)."""
    f = lambda th: np.cos(2 * n * th) * np.sin(4 * J * t * np.sin(th)) / np.sin(th)
    val, _ = quad(f, 1e-12, np.pi / 2, limit=8000)
    return val / (2 * np.pi)


def peak_time_1d(n):
    x = np.linspace(n - 2.0, n + 8.0 * n ** (1 / 3) + 8.0, 400000)
    v = jv(n, x) ** 2
    return x[int(np.argmax(v))] / (2 * J)


# cross-check the 1D representation against the raw 2D Bessel double sum at n=15
i1_2d_15, _, t15_2d = I1_exact(15)
i1_1d_15 = I1_1D(15, peak_time_1d(15))
dev_1d_2d = abs(i1_1d_15 - i1_2d_15)
log(f"\n    cross-check at n=15: 1D repr = {i1_1d_15:.6f}, 2D Bessel sum = {i1_2d_15:.6f},")
log(f"    deviation = {dev_1d_2d:.2e}  (the two exact forms of the same object agree)")
TOL_1D = 1e-3
assert dev_1d_2d < TOL_1D, f"1D-vs-2D I_1 mismatch {dev_1d_2d:.2e}"

# extrapolate the exact 1D I_1(n) to n->inf and compare to the closed form
log(f"\n    exact 1D I_1(n) at the front peak, with an n^{{-2/3}} + n^{{-4/3}} extrapolation:")
log(f"    {'n':>6} {'t*':>10} {'I_1 (1D exact)':>15}")
NS7 = [400, 1600, 6400]
vals7 = []
for n in NS7:
    ts = peak_time_1d(n)
    v = I1_1D(n, ts)
    vals7.append(v)
    log(f"    {n:>6} {ts:>10.4f} {v:>15.7f}")
# Richardson-style extrapolation on the three points: I_1 = c0 + c1 n^{-2/3} + c2 n^{-4/3}
M = np.array([[1.0, n ** (-2 / 3), n ** (-4 / 3)] for n in NS7])
c_extrap = np.linalg.solve(M, np.array(vals7))
I1_extrap = c_extrap[0]
dev_close = abs(I1_CLOSED - I1_extrap)
log(f"\n    extrapolated I_1_inf (n^{{-2/3}} + n^{{-4/3}}) = {I1_extrap:.8f}")
log(f"    closed form                                 = {I1_CLOSED:.8f}")
log(f"    | closed - extrapolated |                   = {dev_close:.2e}")
TOL7 = 1e-4
assert dev_close < TOL7, f"closed I_1 vs extrapolated off by {dev_close:.2e} (> {TOL7:.0e})"
log(f"    PASS below {TOL7:.0e}: the closed form is the n->inf limit of the exact I_1")

# ---------------------------------------------------------------------------
# SECTION 8 : THE FIXED-GAMMA CEILING (the refutation of A_inf = 4)
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[8] THE FIXED-GAMMA CEILING : A_inf(gamma) = 4 - phi(2J)/gamma < 4, NOT 4")
log("-" * 78)

# 8a. closed-form ceiling table
log("\n[8a] A_inf(gamma) = 4 - phi(2J)/gamma, phi(2J) = sqrt(Gamma(Gamma+4J))")
log("     - 4J arcsinh sqrt(Gamma/4J), Gamma = 4 gamma; small-gamma form 4 - (8/3)sqrt(gamma):")
log(f"     {'gamma':>7} {'phi(2J)':>10} {'A_inf (closed)':>15} {'4-(8/3)sqrt(g)':>15} "
    f"{'expected':>9} {'dev':>8}")
GAMMAS8 = [0.002, 0.01, 0.05, 0.10, 0.20]
EXPECT8 = [3.881, 3.734, 3.412, 3.181, 2.872]
maxdev8 = 0.0
for gam, exp in zip(GAMMAS8, EXPECT8):
    phi = phi_of_gamma(gam)
    A_inf = 4.0 - phi / gam
    A_small = 4.0 - (8.0 / 3.0) * np.sqrt(gam / J)
    dev = abs(A_inf - exp)
    maxdev8 = max(maxdev8, dev)
    log(f"     {gam:>7.3f} {phi:>10.6f} {A_inf:>15.4f} {A_small:>15.4f} "
        f"{exp:>9.3f} {dev:>8.4f}")
assert maxdev8 < 0.03, f"closed-form A_inf drifts from expected by {maxdev8:.4f}"
log(f"     max deviation from the stated expectations: {maxdev8:.4f} (rounding of the")
log("     expected list; the closed-form values are the correct ceiling)")

# 8b. the DISCRIMINATOR: exact renewal ladder slope d(ln S)/dn exceeds the
#     single-refill ceiling and climbs toward phi/2 (the sub-linear law falsified)
log("\n[8b] DISCRIMINATOR : exact renewal ladder S_n(t*_0), t*_0 = n/(2J), gamma = 0.10.")
log("     d(ln S)/dn climbs (toward phi/2), EXCEEDING the single-refill ceiling; the")
log("     single-refill law (slope ~ 0.162 n^{-1/3} - (2/3)/n, DECREASING) is falsified.")


def renewal_ladder(dmax, tgrid, gamma, Np=2048):
    """Per-momentum Volterra for the ladder sum S_d(t) (NO e^{-Gamma t} factor).
    Same solver as section [1]/[6], returns S[t_index, d]."""
    Gamma = 4.0 * gamma
    h = tgrid[1] - tgrid[0]
    p = np.arange(Np) * 2 * np.pi / Np
    a = 4.0 * J * np.abs(np.sin(p / 2))
    nt = len(tgrid)
    K = jn(0, np.outer(tgrid, a))
    Shat = np.zeros((nt, Np))
    Shat[0] = K[0]
    for k in range(1, nt):
        conv = 0.5 * K[k] * Shat[0]
        if k > 1:
            conv = conv + np.sum(K[k - 1:0:-1] * Shat[1:k], axis=0)
        Shat[k] = (K[k] + Gamma * h * conv) / (1.0 - 0.5 * Gamma * h)
    d = np.arange(dmax + 1)
    phase = np.exp(1j * np.outer(d, p)) / Np
    return np.real(Shat @ phase.T)


GAM8 = 0.10
NS8 = np.arange(40, 131, 10)
DT8 = 0.05
tg8 = np.arange(0.0, NS8[-1] / (2 * J) + 2.0 + 1e-9, DT8)
S8 = renewal_ladder(int(NS8[-1]), tg8, GAM8)
lnS8 = np.array([np.log(S8[int(round((n / (2 * J)) / DT8)), n]) for n in NS8])
# central-difference local slope at the interior nodes n = 50..120
mid8 = NS8[1:-1]
slope8 = (lnS8[2:] - lnS8[:-2]) / (NS8[2:] - NS8[:-2])
sr8 = 0.162 * mid8.astype(float) ** (-1 / 3) - (2.0 / 3.0) / mid8   # single-refill slope
sr_ceiling = float(sr8.max())
phi8 = phi_of_gamma(GAM8)
log(f"\n     phi(2J) = {phi8:.5f}, phi/2 = {phi8 / 2:.5f}, A_inf = {4 - phi8 / GAM8:.4f}")
log(f"     {'n':>5} {'ln S_n':>10} {'d(lnS)/dn':>11} {'single-refill slope':>21}")
for i, n in enumerate(mid8):
    log(f"     {n:>5} {lnS8[i + 1]:>10.4f} {slope8[i]:>11.5f} {sr8[i]:>21.5f}")
slope_at_120 = float(slope8[-1])
log(f"\n     measured slope at n=120 = {slope_at_120:.5f} (RISING toward phi/2 = {phi8/2:.5f})")
log(f"     single-refill ceiling (max over window) = {sr_ceiling:.5f} (and DECREASING)")
assert slope8[-1] > slope8[0], "renewal ladder slope not climbing (LD signature missing)"
assert slope_at_120 > sr_ceiling, (
    f"slope at n=120 ({slope_at_120:.5f}) does not exceed the single-refill "
    f"ceiling ({sr_ceiling:.5f}) - sub-linear law not falsified")
log("     PASS: the measured slope exceeds the single-refill ceiling and climbs, while")
log("     the single-refill prediction falls => ln S_n grows LINEARLY, A_inf < 4.")

# 8c. the crossover scale n* : all section [5]/[6] windows sit far below it
log("\n[8c] CROSSOVER n* ~ 6 (gamma/J)^{-3/2} (single-refill n^{2/3} = linear phi/2 n):")
log(f"     {'gamma':>7} {'n* = 6 (g/J)^-3/2':>18}")
for gam in [0.002, 0.05, 0.10]:
    n6 = 6.0 * (gam / J) ** (-1.5)
    log(f"     {gam:>7.3f} {n6:>18.0f}")
log("     (all section [5]/[6] windows use n <= 55, far below even the smallest n* ~ 190;")
log("     the climb there is pre-asymptotic, not a plateau at 4.)")

# ---------------------------------------------------------------------------
# SECTION 9 : THE o(n) : prefactor, approach law, and the third reading
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[9] THE o(n) : prefactor (power -1/2, closed constant), approach law, A_max = 0")
log("-" * 78)
log("    ln S_n(t*_0) = (phi/2J) n - (1/2) ln n + ln C(gamma) + O(1/n), t*_0 = n/(2J),")
log("    C(gamma) = (2 pi)^{-1/2} (gamma/(gamma+J))^{1/4} = (Gamma/mu*) sqrt(J/(pi mu''(th*))).")
log("    Referee-confirmed 2026-07-13 (SymPy + steepest-descent + direct-Lindblad + Haken-Strobl).")


# 9a. the saddle constants: mu''(theta*) closed form vs 5-point finite difference
log("\n[9a] SADDLE CONSTANTS : mu(theta) = sqrt(Gamma^2 + 16 J^2 sinh^2(theta/2)),")
log("     theta* = 2 arcsinh sqrt(Gamma/4J); mu''(theta*) = 2J sqrt(Gamma/(Gamma+4J)) checked")
log("     vs a 5-point finite-difference mu'' to <= 1e-8; C(gamma) = (2pi)^-1/2 (g/(g+J))^1/4.")


def _mu(theta, gamma):
    G = 4.0 * gamma
    return np.sqrt(G * G + 16.0 * J * J * np.sinh(theta / 2.0) ** 2)


log(f"\n     {'gamma':>7} {'theta*':>9} {'mu\'(FD)':>9} {'mu\"(FD)':>12} "
    f"{'mu\"(closed)':>12} {'dev':>10} {'C(gamma)':>10}")
maxdev9a = 0.0
for gam in [0.05, 0.10, 0.20]:
    G = 4.0 * gam
    th = 2.0 * np.arcsinh(np.sqrt(G / (4.0 * J)))
    hh = 1e-3
    m_m2, m_m1, m_0, m_p1, m_p2 = (_mu(th - 2 * hh, gam), _mu(th - hh, gam),
                                   _mu(th, gam), _mu(th + hh, gam), _mu(th + 2 * hh, gam))
    mu1_fd = (m_m2 - 8 * m_m1 + 8 * m_p1 - m_p2) / (12.0 * hh)        # -> 2J
    mu2_fd = (-m_m2 + 16 * m_m1 - 30 * m_0 + 16 * m_p1 - m_p2) / (12.0 * hh * hh)
    mu2_cl = 2.0 * J * np.sqrt(G / (G + 4.0 * J))
    dev = abs(mu2_fd - mu2_cl)
    maxdev9a = max(maxdev9a, dev)
    log(f"     {gam:>7.2f} {th:>9.5f} {mu1_fd:>9.5f} {mu2_fd:>12.7f} "
        f"{mu2_cl:>12.7f} {dev:>10.2e} {C_of_gamma(gam):>10.6f}")
TOL9A = 1e-8
assert maxdev9a < TOL9A, f"mu'' finite-difference vs closed off by {maxdev9a:.2e}"
log(f"     max |mu''(FD) - mu''(closed)| = {maxdev9a:.2e}; PASS below {TOL9A:.0e} (mu'(FD) = 2J)")

# 9b. the -1/2 power and the constant, read directly off the exact renewal ladder
log("\n[9b] POWER -1/2 + CONSTANT : exact renewal ladder (the section [8b] solver), gamma=0.10,")
log("     lnC_est = ln S_n(t*_0) - (phi/2J) n + (1/2) ln n converges to ln C(0.1) = -1.5184")
log("     iff the power is exactly -1/2 and the constant is C(0.1).")
GAM9 = 0.10
phi9 = phi_of_gamma(GAM9)
lnC9 = np.log(C_of_gamma(GAM9))
DT9 = 0.03
NP9 = 4096
NS9B = [80, 100, 120, 140, 160]
tg9 = np.arange(0.0, NS9B[-1] / (2 * J) + 1.0 + 1e-9, DT9)
S9 = renewal_ladder(NS9B[-1], tg9, GAM9, Np=NP9)
log(f"\n     phi/2J = {phi9 / (2 * J):.6f}, ln C(0.1) = {lnC9:.4f}, A_inf(0.1) = {4 - phi9 / GAM9:.4f}")
log(f"     {'n':>5} {'S_n(t*_0)':>13} {'lnC_est':>10} {'lnC_est - lnC':>14}")
lnCest_160 = None
for n in NS9B:
    ti = int(round((n / (2 * J)) / DT9))
    Sn = S9[ti, n]
    lnCest = np.log(Sn) - (phi9 / (2 * J)) * n + 0.5 * np.log(n)
    if n == 160:
        lnCest_160 = lnCest
    log(f"     {n:>5} {Sn:>13.5e} {lnCest:>10.5f} {lnCest - lnC9:>14.5f}")
dev9b = abs(lnCest_160 - lnC9)
assert dev9b <= 0.05, f"lnC_est(n=160) off from ln C by {dev9b:.4f} (> 0.05)"
log(f"\n     |lnC_est(n=160) - ln C(0.1)| = {dev9b:.4f}; PASS below 0.05 => the power is exactly")
log(f"     -1/2, the constant is C(0.1) = {C_of_gamma(GAM9):.6f} (a single interior Gaussian saddle).")

# 9c. the approach law A_eff(n) = A_inf + (2J/(gamma n))[-(1/6) ln n + ln(c_G/C)]
log("\n[9c] APPROACH LAW : A_eff(n) = A_inf + (2J/(gamma n))[-(1/6) ln n + ln(c_G/C(gamma))],")
log("     c_G = 2^{2/3} Ai(-alpha)^2 = 0.45547 (the coherent-peak denominator). vs renewal A_eff.")
A_inf9 = 4.0 - phi9 / GAM9
cG_over_C = CAUSTIC_PEAK / C_of_gamma(GAM9)
log(f"\n     A_inf(0.1) = {A_inf9:.4f}, c_G = {CAUSTIC_PEAK:.5f}, c_G/C = {cG_over_C:.4f}")
log(f"     {'n':>5} {'A_eff(renewal)':>15} {'A_eff(law)':>11} {'|diff|':>9}")
maxdiff9c = 0.0
for n in [60, 100, 160]:
    ti = int(round((n / (2 * J)) / DT9))
    Sn = S9[ti, n]
    K = GAM9 * n / (2 * J)
    g = (np.exp(-4 * GAM9 * (n / (2 * J))) * Sn) / (CAUSTIC_PEAK * n ** (-2 / 3))
    A_num = -np.log(g) / K
    A_law = A_inf9 + (2 * J / (GAM9 * n)) * (-(1.0 / 6.0) * np.log(n) + np.log(cG_over_C))
    diff = abs(A_num - A_law)
    maxdiff9c = max(maxdiff9c, diff)
    log(f"     {n:>5} {A_num:>15.4f} {A_law:>11.4f} {diff:>9.4f}")
assert maxdiff9c <= 0.02, f"approach-law vs renewal A_eff off by {maxdiff9c:.4f} (> 0.02)"
log(f"     max |A_eff(renewal) - A_eff(law)| = {maxdiff9c:.4f}; PASS below 0.02")

# the non-monotone crossing of A_inf (the coherent-peak denominator reading)
n_x = float(np.exp(6.0 * np.log(cG_over_C)))
sign_lo = sign_hi = None
prev = None
for n in range(50, 141):
    ti = int(round((n / (2 * J)) / DT9))
    Sn = S9[ti, n]
    K = GAM9 * n / (2 * J)
    g = (np.exp(-4 * GAM9 * (n / (2 * J))) * Sn) / (CAUSTIC_PEAK * n ** (-2 / 3))
    diff = -np.log(g) / K - A_inf9
    if prev is not None and prev > 0.0 >= diff:
        sign_lo, sign_hi = n - 1, n
        break
    prev = diff
log(f"\n     predicted crossing n_x = exp(6 ln(c_G/C)) = {n_x:.1f}")
log(f"     measured A_eff - A_inf sign change (+ -> -) between n = {sign_lo} and n = {sign_hi}")
assert sign_lo is not None and 70 <= sign_lo <= 110, \
    f"A_eff - A_inf crossing outside the expected n = 80..100 window (got {sign_lo})"
log("     convention flag: the crossing exists for the coherent-PEAK denominator reading;")
log("     a same-instant denominator removes it (A_eff approaches A_inf monotonically from")
log("     below). A_inf and the -(1/6) ln n term are convention-robust; only the additive")
log("     O(ln n / n) constant (hence the crossing location) is not.")

# 9d. the third reading, A_max = 0 : the exact Haken-Strobl diffusive peak height
log("\n[9d] THIRD READING A_max = 0 : exact Haken-Strobl diffusive walk")
log("     P_n(t) = e^{-2Dt} I_n(2Dt) = ive(n, 2Dt), D = 2 J^2 / Gamma. The GLOBAL max over t")
log("     is the diffusive plateau: t_peak = n^2/(2D), peak height ~ 1/n (algebraic).")
GAMD = 0.10
DD = 2.0 * J * J / (4.0 * GAMD)                 # D = 2 J^2 / Gamma
log(f"\n     gamma = {GAMD}, Gamma = {4 * GAMD}, D = 2 J^2 / Gamma = {DD:.4f}")
log(f"     {'n':>5} {'t_peak(num)':>12} {'n^2/(2D)':>11} {'dev %':>8} {'Pmax*n':>10}")
maxtdev = 0.0
pmaxn_320 = None
for n in [80, 160, 320]:
    x = np.linspace(0.6 * n * n, 1.4 * n * n, 60001)
    dx = x[1] - x[0]
    y = ive(n, x)
    i = int(np.argmax(y))
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    den = y0 - 2.0 * y1 + y2
    xpk = x[i] + 0.5 * (y0 - y2) / den * dx if den < 0 else x[i]
    t_peak = xpk / (2.0 * DD)
    t_hs = n * n / (2.0 * DD)
    devpc = abs(t_peak - t_hs) / t_hs * 100.0
    maxtdev = max(maxtdev, devpc)
    Pmax = float(ive(n, xpk))
    if n == 320:
        pmaxn_320 = Pmax * n
    log(f"     {n:>5} {t_peak:>12.2f} {t_hs:>11.2f} {devpc:>8.4f} {Pmax * n:>10.5f}")
assert maxtdev <= 1.0, f"t_peak off from n^2/(2D) by {maxtdev:.3f}% (> 1%)"
HS_LIMIT = np.exp(-0.5) / np.sqrt(2.0 * np.pi)   # 0.24197
dev_pmaxn = abs(pmaxn_320 - HS_LIMIT) / HS_LIMIT * 100.0
assert dev_pmaxn <= 1.0, f"Pmax*n(320) off from e^-1/2/sqrt(2pi) by {dev_pmaxn:.3f}% (> 1%)"
log(f"\n     t_peak = n^2/(2D) to {maxtdev:.4f}% (<= 1%); Pmax*n(320) = {pmaxn_320:.5f} vs")
log(f"     e^{{-1/2}}/sqrt(2pi) = {HS_LIMIT:.5f} ({dev_pmaxn:.4f}% <= 1%). Peak height ~ 1/n:")
log("     the global-max peak-value survival exponent is 0 (algebraic decay, psi = Gamma).")
log("\n     TRICHOTOMY: the fixed-time ceiling A_inf = 4 - phi/gamma is the arrival-instant")
log("     amplitude; a narrow front window interpolates between A_inf and 0; the global-max")
log("     exponent is 0. The section [5] EDGE exclusions are exactly this diffusive plateau")
log("     entering the front window (the windowed max jumping off the front onto the plateau).")

# ---------------------------------------------------------------------------
# SECTION 10 : THE FENCE (where the ballistic arrival reading dies: a dose)
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[10] THE FENCE : where the ballistic arrival reading dies (a dose, not a ratio)")
log("-" * 78)
log("     interior-seed renewal front-window survival g(n,gamma) = max over [0.5,1.3] n/(2J)")
log("     of P_n(t), divided by the clean (gamma=0) front peak. The arrival reading dies when")
log("     g drops below theta. The edge is a fixed dephasing DOSE K_deg = gamma_deg n /(2J) ~ 0.6,")
log("     NOT a fixed gamma/J. Docks the gamma ~ J fence of COUPLING_DEFECT_WALK_TIME_STEP.md")
log("     onto the owned coherence-horizon dispersion (EP gamma*(q) = J q / 2).")


def _front_g(n, gamma, tg, Pc_peak):
    """Front-window survival g(n,gamma): max over [0.5,1.3] n/(2J) of the interior
    renewal P_n(t), divided by the clean front peak Pc_peak. Monotone decreasing in gamma."""
    Pg, _ = renewal_P(n, tg, gamma)
    col = Pg[:, n]
    dt = tg[1] - tg[0]
    tball = n / (2.0 * J)
    lo = int(0.5 * tball / dt)
    hi = int(1.3 * tball / dt)
    return float(np.max(col[lo:hi])) / Pc_peak


def _gamma_deg(n, theta, tg, Pc_peak, glo=0.02, ghi=0.16):
    """Bisect gamma_deg where g(n,gamma) crosses theta (g monotone decreasing in gamma)."""
    for _ in range(16):
        gm = 0.5 * (glo + ghi)
        if _front_g(n, gm, tg, Pc_peak) > theta:
            glo = gm
        else:
            ghi = gm
        if ghi - glo < 1e-4:
            break
    return 0.5 * (glo + ghi)


# 10a. n = 20, 30, 40 : the degeneration dose K_deg at theta = 0.2
log("\n[10a] DEGENERATION DOSE : theta = 0.2, interior seed n = 20, 30, 40 (exact renewal).")
log("      expect K_deg ~ 0.69 / 0.64 / 0.61 (weakly decreasing: the pre-asymptotic A-climb).")
log(f"      {'n':>4} {'gamma_deg':>10} {'Q_deg=J/g':>10} {'K_deg':>8}")
THETA_A = 0.2
Kdeg_a = []
gdeg_a = []
for n in [20, 30, 40]:
    tball = n / (2.0 * J)
    tgA = np.arange(0.0, 1.5 * tball + 1e-9, 0.04)
    Pc0, _ = renewal_P(n, tgA, 0.0)
    Pc_peak = float(np.max(Pc0[:, n]))
    gdeg = _gamma_deg(n, THETA_A, tgA, Pc_peak)
    Qdeg = J / gdeg
    Kdeg = gdeg * n / (2.0 * J)
    Kdeg_a.append(Kdeg)
    gdeg_a.append(gdeg)
    log(f"      {n:>4} {gdeg:>10.5f} {Qdeg:>10.3f} {Kdeg:>8.4f}")
for n, K in zip([20, 30, 40], Kdeg_a):
    assert 0.5 <= K <= 0.8, f"K_deg(n={n}) = {K:.4f} outside [0.5, 0.8]"
log("\n      K_deg ~ 0.6 across n = 20..40 (weakly decreasing): a fixed DOSE, not a fixed")
log("      gamma/J. gamma_deg ~ 1.2 J/n (from K_deg ~ 0.6): the degeneration gamma falls with n.")

# 10b. theta-scan at n = 20 : K_deg moves with ln(1/theta), NOT a single-A identity
log("\n[10b] THETA-SCAN at n = 20 : K_deg vs theta. The implied A = ln(1/theta)/K_deg is NOT")
log("      constant (ln g is concave in K: the halo boost softens the slope at larger dose),")
log("      so K_deg tracks ln(1/theta) in DIRECTION only, not as a single-A identity.")
n_b = 20
tball_b = n_b / (2.0 * J)
tgB = np.arange(0.0, 1.5 * tball_b + 1e-9, 0.04)
Pc0b, _ = renewal_P(n_b, tgB, 0.0)
Pc_peak_b = float(np.max(Pc0b[:, n_b]))
log(f"\n      {'theta':>6} {'gamma_deg':>10} {'K_deg':>8} {'A=ln(1/theta)/K_deg':>20}")
Kdeg_b = []
for theta in [0.1, 0.2, 0.4]:
    gdeg = _gamma_deg(n_b, theta, tgB, Pc_peak_b)
    Kdeg = gdeg * n_b / (2.0 * J)
    A_impl = np.log(1.0 / theta) / Kdeg
    Kdeg_b.append(Kdeg)
    log(f"      {theta:>6.1f} {gdeg:>10.5f} {Kdeg:>8.4f} {A_impl:>20.4f}")
assert Kdeg_b[0] > Kdeg_b[1] > Kdeg_b[2], f"K_deg not decreasing with theta: {Kdeg_b}"
log(f"\n      K_deg falls {Kdeg_b[0]:.3f} -> {Kdeg_b[1]:.3f} -> {Kdeg_b[2]:.3f} as theta rises")
log("      0.1 -> 0.2 -> 0.4; the implied A climbs ~2.1 -> 2.5 (a ~20% drift), so ln(1/theta)")
log("      is the right direction but K_deg = ln(1/theta)/A is a ~25% approximation, not exact.")

# 10c. the docking summary : no n-free front EP; the n-free content is the small-N corner
log("\n[10c] DOCKING SUMMARY (onto the owned coherence-horizon dispersion gamma*(q) = J q / 2):")
log("      * NO n-free front EP. The front's coherence content is q -> 0 in the owned")
log("        dispersion (group velocity domega/dq = 2J, the ballistic speed), where the EP")
log("        gamma*(q) = J q / 2 -> 0. The naive band-momentum substitution gamma* = pi J / 4")
log("        (q = pi/2) is a wrong-variable reading: q is the coherence centre-of-mass")
log("        wavevector, not a band momentum, and pi/2 is outside the q -> 0 validity.")
log("      * The n-free content is the small-N EP corner Q*(2) = 1 (gamma = J, N=2 critical")
log("        damping) and Q*(3) = sqrt(2) = 1.414 (the ibm_kingston brackets Q = 1 and Q ~ 1.5),")
log("        and the ENAQT front-dominant window L* = Q/2 (ballistic edge = diffusive spread).")
Qdeg_over_n = [(J / g) / n for g, n in zip(gdeg_a, [20, 30, 40])]
slope_Q = float(np.mean(Qdeg_over_n))
log(f"      * measured Q_deg / n = {Qdeg_over_n[0]:.3f}, {Qdeg_over_n[1]:.3f}, {Qdeg_over_n[2]:.3f} "
    f"(mean ~ {slope_Q:.2f}); extrapolating Q_deg ~ {slope_Q:.2f} n to n = 1, 2 gives")
log(f"        Q ~ {slope_Q * 1:.2f} (n=1) and {slope_Q * 2:.2f} (n=2), landing on the Q*(2) = 1 /")
log("        Q*(3) = 1.41 bracket (the ibm_kingston Q = 1 and Q ~ 1.5).")
log("      => the fence is a DOSE K_deg ~ 0.6, gamma_deg ~ 1.2 J/n; gamma ~ J is only the")
log("         n -> 1 EP corner of the same coherence-horizon ladder (large-N slope 2/pi).")

log("\n" + "=" * 78)
log("SUMMARY: A_inf = gamma_phi = 4 is REFUTED. Pre-asymptotically A_eff climbs like")
log("4 - 4.864 n^{-1/3} (closed-form I_1 = 0.27694424); the true fixed-gamma ceiling is")
log("A_inf(gamma) = 4 - phi(2J)/gamma < 4 (small gamma: 4 - (8/3)sqrt(gamma/J)). ln S_n")
log("grows LINEARLY at rate phi/2, not O(ln n); the halo is an exp-in-n front boost. The")
log("measured gamma-drift of the K-collapse remains mostly the max-ratio peak-shift artifact.")
log("The o(n) is CLOSED (section [9]): S_n(t*_0) = C(gamma) n^{-1/2} e^{(phi/2J)n} with the")
log("closed constant C(gamma) = (2 pi)^{-1/2}(gamma/(gamma+J))^{1/4}; A_eff approaches A_inf as")
log("-(J/(3 gamma n)) ln n (alpha_ln = -1/6), and the global-max peak-value exponent is 0")
log("(Haken-Strobl). The gamma ~ J fence is docked (section [10]): the edge is a dephasing dose")
log("K_deg ~ 0.6 (theta = 0.2, gamma_deg ~ 1.2 J/n); gamma ~ J is only the n -> 1 EP corner.")
log("=" * 78)
log("DONE")

os.makedirs("simulations/results/cone_defect_arrival", exist_ok=True)
with open("simulations/results/cone_defect_arrival/front_survival_asymptote.txt", "w",
          encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
