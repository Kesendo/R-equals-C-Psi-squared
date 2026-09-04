"""Converged re-measurement of experiments/STAR_TOPOLOGY_OBSERVERS.md.

The March 2026 readings in that document were taken with `run_star_topology` at
its default `sample_every=20`: the observable is recorded every 0.1 time units
while the oscillation period is about 0.66, so peaks are systematically missed.
Because the peak of CPsi is a supremum over continuous time, an undersampled
peak is always too LOW and a threshold read from it always too HIGH, which is
why every quantity re-measured here moved in one direction.

This script computes the converged values by a route that is not the document's
own:

  * an EXACT propagator, expm(L*dt) applied to the vectorised density matrix,
    so the only discretisation left is where the peak is sampled and not how the
    state got there;
  * the repo's RK4 (`star_topology_v3.py`) run beside it as an independent
    second route, at dt=0.001 recording every step;
  * the peak over t refined by a parabola through the sampled maximum;
  * window edges located by linear interpolation, not by counting grid points.

Every number that document labels "converged" is printed by one of the blocks
below. Run one block by name, or `all` for the lot (about ten minutes).

    python simulations/star_topology_converged.py threshold
    python simulations/star_topology_converged.py all

Conventions match `star_topology_v3.py`: H = J_SA(sigma_S.sigma_A) +
J_SB(sigma_S.sigma_B), Lindblad operators sqrt(gamma_q)*sigma_z, initial state
Bell_SA (x) |+>_B, and the observable is CPsi = concurrence * Psi on a reduced
pair (NOT R = C*Psi^2, which the same module also computes).
"""

import sys
import time
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, "simulations")
import star_topology_v3 as v3


# --------------------------------------------------------------- the exact route

def liouvillian(H, gammas, n_qubits):
    """Vectorised generator, column-stacking convention vec(ABC) = (C^T kron A) vec(B)."""
    d = 2 ** n_qubits
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for q, g in enumerate(gammas):
        if g <= 0:
            continue
        Z = v3.op_on_qubit(v3.sz, q, n_qubits)
        L += g * (np.kron(Z.T, Z) - np.kron(Id, Id))
    return L


def cpsi_series(n_obs, J_SA, J_SB, gammas, J_AB=0.0, dt=0.002, t_max=5.0, pair=(1, 2)):
    """CPsi of one reduced pair over time, exact in dt."""
    n_q = 1 + n_obs
    d = 2 ** n_q
    H = v3.star_hamiltonian_n(n_obs, J_SA=J_SA, J_SB=J_SB, J_AB=J_AB)
    P = expm(liouvillian(H, gammas, n_q) * dt)
    vec = v3.bell_sa_plus_rest(n_obs).reshape(-1, order="F").astype(complex)
    steps = int(round(t_max / dt))
    out = np.empty(steps + 1)
    for s in range(steps + 1):
        r = vec.reshape(d, d, order="F")
        r = 0.5 * (r + r.conj().T)
        out[s] = v3.pair_metrics(v3.partial_trace_keep(r, keep=list(pair), n_qubits=n_q))["cpsi"]
        if s < steps:
            vec = P @ vec
    return np.arange(steps + 1) * dt, out


def refined_peak(ts, vals):
    """Peak value and location, parabola through the sampled maximum."""
    i = int(np.argmax(vals))
    if 0 < i < len(vals) - 1:
        y0, y1, y2 = vals[i - 1], vals[i], vals[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            s = 0.5 * (y0 - y2) / den
            return y1 - 0.25 * (y0 - y2) * s, ts[i] + s * (ts[1] - ts[0])
    return vals[i], ts[i]


def peak(J_SA, J_SB, gammas, J_AB=0.0, dt=0.002, t_max=5.0, n_obs=2):
    return refined_peak(*cpsi_series(n_obs, J_SA, J_SB, gammas, J_AB, dt, t_max))[0]


def windows(J_SA, J_SB, gammas, dt=0.00025, t_max=40.0, level=0.25):
    """Lengths of every excursion above `level`, with the edges interpolated."""
    ts, ab = cpsi_series(2, J_SA, J_SB, gammas, dt=dt, t_max=t_max)
    a = ab >= level
    spans = []
    i, n = 0, len(a)
    while i < n:
        if a[i]:
            j = i
            while j + 1 < n and a[j + 1]:
                j += 1
            lo, hi = ts[i], ts[j]
            if i > 0:
                lo = ts[i - 1] + (level - ab[i - 1]) / (ab[i] - ab[i - 1]) * dt
            if j + 1 < n:
                hi = ts[j] + (ab[j] - level) / (ab[j] - ab[j + 1]) * dt
            spans.append(hi - lo)
            i = j + 1
        else:
            i += 1
    return spans


def bisect(f, lo, hi, target=0.25, tol=2e-6):
    """Smallest x in [lo, hi] with f(x) >= target. f must be monotone on the bracket."""
    if f(lo) >= target or f(hi) < target:
        return None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if f(mid) >= target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def r2(y, yhat):
    return 1 - np.sum((y - yhat) ** 2) / np.sum((y - np.mean(y)) ** 2)


# ------------------------------------------------------------------------ blocks

def threshold():
    """Section 4.10: the threshold, both routes, and the grid it needs."""
    g = [0.05] * 3

    def march(J):
        rec = v3.run_star_topology(2, J_SA=1.0, J_SB=J, gammas=g, dt=0.005,
                                   t_max=5.0, sample_every=20)
        return max(rec.pairs["AB"].cpsi)

    print("March's own grid (dt=0.005, sample_every=20) reproduces to the last digit:")
    for J in (1.4650, 1.4655):
        print("   J_SB=%.4f: %.6f" % (J, march(J)))
    lo, hi = 1.4600, 1.4700
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if march(mid) >= 0.25:
            hi = mid
        else:
            lo = mid
    print("   so that grid puts the threshold at %.6f, not at the 1.466 March states"
          % (0.5 * (lo + hi)))

    print("converged, exact propagator, step halved:")
    for dt in (0.008, 0.004, 0.002, 0.001):
        print("   dt=%.3f: threshold J_SB = %.6f"
              % (dt, bisect(lambda J: peak(1.0, J, g, dt=dt), 1.30, 1.60)))

    def rk4(J):
        rec = v3.run_star_topology(2, J_SA=1.0, J_SB=J, gammas=g, dt=0.001,
                                   t_max=5.0, sample_every=1)
        return refined_peak(np.array(rec.t), np.array(rec.pairs["AB"].cpsi))[0]

    print("converged, the repo's own RK4 recording every step:")
    print("   dt=0.001: threshold J_SB = %.6f" % bisect(rk4, 1.30, 1.60))


def scale():
    """Section 4.10: Q = J_SA/gamma is the knob, and absolute coupling is not a second one."""
    print("peak CPsi_AB at ratio 2 with Q = 20 held, everything scaled by s:")
    for gam in (0.0125, 0.025, 0.05, 0.1, 0.5, 1.0, 2.0):
        s = gam / 0.05
        print("   gamma=%-7.4f J=(%.4f, %.4f): peak = %.9f"
              % (gam, 1.0 * s, 2.0 * s,
                 peak(1.0 * s, 2.0 * s, [gam] * 3, dt=0.002 / s, t_max=5.0 / s)))
    print("the same check on a FIXED grid, which is the one that can fail:")
    print("   (scaling the grid with s makes expm(L(sJ,sg)*dt/s) the SAME matrix as")
    print("    expm(L(J,g)*dt), so the runs above are arithmetically identical)")
    for gam in (0.0125, 0.025, 0.05, 0.1, 0.5):
        s_ = gam / 0.05
        print("   gamma=%-7.4f J=(%.4f, %.4f) at dt=0.001, t_max=8: peak = %.9f"
              % (gam, 1.0 * s_, 2.0 * s_, peak(1.0 * s_, 2.0 * s_, [gam] * 3, dt=0.001, t_max=8.0)))
    print("March's two 'absolute coupling' rows, and the same two Q reached by moving gamma:")
    for J_SA, J_SB, gam in ((0.50, 0.75, 0.05), (1.00, 1.50, 0.05), (1.00, 1.50, 0.10)):
        p = peak(J_SA, J_SB, [gam] * 3, dt=0.002)
        print("   J=(%.2f, %.2f) gamma=%.2f  Q=%3.0f  peak = %.6f  %s"
              % (J_SA, J_SB, gam, J_SA / gam, p, "crosses" if p >= 0.25 else "does NOT cross"))


def noise():
    """Section 4.9: the two gamma scans, the fixed-product test, and the window column."""
    print("the two scans (gamma_S=0.05, J=(1,2)); Window = longest single excursion")
    for lbl, gA, gB in (("gamma_A=0.1, gamma_B=0.001", 0.1, 0.001),
                        ("gamma_A=0.1, gamma_B=0.1", 0.1, 0.1),
                        ("gamma_A=0.1, gamma_B=0.2", 0.1, 0.2),
                        ("gamma_A=0.001, gamma_B=0.1", 0.001, 0.1),
                        ("gamma_A=0.2, gamma_B=0.1", 0.2, 0.1)):
        p = peak(1.0, 2.0, [0.05, gA, gB], dt=0.002, t_max=40.0)
        w = windows(1.0, 2.0, [0.05, gA, gB])
        print("   %-28s peak=%.4f  window=%.5f  lobes=%d"
              % (lbl, p, max(w) if w else 0.0, len(w)))
    print("fixed-product test, gamma_A * gamma_B = 0.0025:")
    for gA, gB in ((0.005, 0.5), (0.05, 0.05), (0.5, 0.005)):
        p = peak(1.0, 2.0, [0.05, gA, gB], dt=0.002, t_max=40.0)
        print("   gamma_A=%-6s gamma_B=%-6s peak=%.4f  %s"
              % (gA, gB, p, "crosses" if p >= 0.25 else "does NOT cross"))
    # dt=0.001 rather than 0.004: the longest window and the lobe count are
    # step-stable at 0.004 but the TOTAL is not, and it is the total the
    # document's table prints (9.4237 at 0.004 against 9.4259 at 0.001 and
    # 9.4260 at 0.0005, for gamma=0.001).
    print("the window column, all three rates equal, converged in run length (t_max=160):")
    for gam in (0.001, 0.01, 0.05, 0.1, 0.2):
        w = windows(1.0, 2.0, [gam] * 3, dt=0.001, t_max=160.0)
        print("   gamma=%-7s longest=%.3f  total=%.3f  lobes=%d"
              % (gam, max(w) if w else 0.0, sum(w), len(w)))
    print("   the same total on shorter runs at gamma=0.001, which is what March had:")
    for t_max in (5.0, 20.0, 40.0, 80.0):
        w = windows(1.0, 2.0, [0.001] * 3, dt=0.001, t_max=t_max)
        print("      t_max=%-5s total=%.3f  lobes=%d" % (t_max, sum(w), len(w)))


def curve():
    """Sections 4.10 and 8.5: the threshold as a function of Q, and the fit."""
    doc = {0.001: 1.183, 0.010: 1.247, 0.020: 1.296, 0.050: 1.466, 0.070: 1.634,
           0.100: 1.820, 0.120: 1.929, 0.150: 2.146, 0.170: 2.253, 0.200: 2.460}
    print("%7s %7s %10s %7s" % ("gamma", "Q", "converged", "March"))
    gam, th = [], []
    for g, dj in doc.items():
        t = bisect(lambda J: peak(1.0, J, [g] * 3, dt=0.002,
                                  t_max=14.0 if g < 0.01 else 6.0), 1.05, 3.4)
        gam.append(g)
        th.append(t)
        print("%7.3f %7.1f %10.5f %7.3f" % (g, 1.0 / g, t, dj))
    gam, th = np.array(gam), np.array(th)

    best = None
    for c in np.arange(0.5, 2.001, 0.001):
        A = np.vstack([gam ** c, np.ones_like(gam)]).T
        coef = np.linalg.lstsq(A, th, rcond=None)[0]
        s = r2(th, A @ coef)
        if best is None or s > best[0]:
            best = (s, c, coef)
    s, c, coef = best
    res = np.max(np.abs(th - (coef[0] * gam ** c + coef[1])))
    print("  power fit: %.4f*gamma^%.4f + %.4f   R2=%.6f  maxres=%.5f"
          % (coef[0], c, coef[1], s, res))
    A = np.vstack([gam, np.ones_like(gam)]).T
    coef = np.linalg.lstsq(A, th, rcond=None)[0]
    print("  linear   : %.4f*gamma + %.4f   R2=%.6f  maxres=%.5f"
          % (coef[0], coef[1], r2(th, A @ coef), np.max(np.abs(th - A @ coef))))
    d = np.array(list(doc.values()))
    print("  March's published fit 7.35*gamma^1.08 + 1.18 on March's own ten rows:"
          " maxres=%.5f, so it does not come from that table"
          % np.max(np.abs(d - (7.35 * gam ** 1.08 + 1.18))))

    print("the quiet end: no plateau, and where the reading stops being about the system")
    for g, t_max in ((0.001, 20.0), (0.0003, 30.0), (0.0001, 40.0)):
        print("   gamma=%-8s Q=%7.0f: %.5f"
              % (g, 1.0 / g,
                 bisect(lambda J: peak(1.0, J, [g] * 3, dt=0.002, t_max=t_max), 1.05, 1.45)))
    for g in (0.05, 0.01, 0.001, 0.0001, 0.0):
        row = [bisect(lambda J: peak(1.0, J, [g] * 3, dt=0.002, t_max=t), 1.02, 1.60, tol=1e-5)
               for t in (5.0, 10.0, 20.0, 40.0, 80.0)]
        print("   gamma=%-8s threshold at t_max 5/10/20/40/80: %s"
              % (g, " ".join("%.5f" % x for x in row)))
    print("   and at gamma=0 with SYMMETRIC coupling, which never crosses however long the run:")
    for t_max in (5.0, 20.0, 80.0, 200.0):
        print("      t_max=%6s: peak = %.5f"
              % (t_max, peak(1.0, 1.0, [0.0] * 3, dt=0.002, t_max=t_max)))


def coarse():
    """Sections 4.9 and 4.10: the coarse ratio scan, and the boundaries Section 7 quotes."""
    print("the coarse ratio scan of Section 4.10 (gamma=0.05, J_SA=1), converged:")
    for J, march in ((1.30, 0.226), (1.40, 0.242), (1.47, 0.251), (2.00, 0.329), (3.00, 0.406)):
        p = peak(1.0, J, [0.05] * 3, dt=0.002)
        print("   J_SB=%.2f: converged %.4f   March %.3f   %s"
              % (J, p, march, "crosses" if p >= 0.25 else "does not"))
    print("the two noise boundaries at matched partner rates (J=(1,2), gamma_S=0.05):")
    print("   %-14s %-12s %-12s %s" % ("partner rate", "gamma_A", "gamma_B", "sender/receiver"))
    for pr in (0.005, 0.05, 0.1, 0.2):
        a = bisect(lambda g: -peak(1.0, 2.0, [0.05, g, pr], dt=0.002, t_max=40.0),
                   0.02, 0.9, target=-0.25)
        b = bisect(lambda g: -peak(1.0, 2.0, [0.05, pr, g], dt=0.002, t_max=40.0),
                   0.02, 0.9, target=-0.25)
        print("   %-14s %-12.4f %-12.4f %.2f%s"
              % (pr, a, b, b / a, "   INVERTED" if b < a else ""))
    print("where the two boundaries cross, which is a solvable point and not a sampled row.")
    print("They can only meet where both rates are equal, so it is the gamma solving")
    print("peak CPsi_AB(gamma_A = gamma_B = gamma, gamma_S = 0.05) = 1/4:")
    sym = bisect(lambda x: -peak(1.0, 2.0, [0.05, x, x], dt=0.002, t_max=40.0),
                 0.05, 0.5, target=-0.25, tol=1e-6)
    a = bisect(lambda x: -peak(1.0, 2.0, [0.05, x, sym], dt=0.002, t_max=40.0),
               0.02, 0.9, target=-0.25, tol=1e-6)
    b = bisect(lambda x: -peak(1.0, 2.0, [0.05, sym, x], dt=0.002, t_max=40.0),
               0.02, 0.9, target=-0.25, tol=1e-6)
    print("   gamma = %.5f, and read as a partner rate it gives gamma_A = %.5f,"
          " gamma_B = %.5f, ratio %.4f" % (sym, a, b, b / a))

    print("the two large-Q thresholds read at the Section 8.5 table's own t_max, for comparison:")
    for g in (0.0003, 0.0001):
        print("   gamma=%-8s t_max=6: %.5f"
              % (g, bisect(lambda J: peak(1.0, J, [g] * 3, dt=0.002, t_max=6.0), 1.05, 1.45)))


def observers():
    """Section 8.1: the N-observer thresholds."""
    print("N=2 and N=3 thresholds, and the grid spread at N=3:")
    print("   N=2: %.6f" % bisect(lambda J: peak(1.0, J, [0.05] * 3, dt=0.002), 1.30, 1.60))
    for dt in (0.004, 0.002, 0.001):
        print("   N=3 at dt=%s: %.6f"
              % (dt, bisect(lambda J: peak(1.0, J, [0.05] * 4, dt=dt, t_max=6.0, n_obs=3),
                            3.5, 3.8)))
    print("   N=3 monotonicity, checked rather than assumed:")
    prev = None
    for J in np.arange(2.0, 5.01, 0.125):
        p = peak(1.0, J, [0.05] * 4, dt=0.004, t_max=6.0, n_obs=3)
        assert prev is None or p > prev, "not monotone at J_SB=%s" % J
        if prev is None:
            first = p
        prev = p
    print("      rises without a dip from %.5f at J_SB=2.0 to %.5f at 5.0" % (first, prev))
    print("   N=4 at March's two probes:")
    for J in (2.0, 4.5):
        print("      J_SB=%s: %.4f" % (J, peak(1.0, J, [0.05] * 5, dt=0.006, t_max=6.0, n_obs=4)))


def direct():
    """Section 8.3: direct A-B coupling, and why the shadow column reads exactly zero."""
    print("   threshold at J_AB=0.0: %.5f"
          % bisect(lambda J: peak(1.0, J, [0.05] * 3, J_AB=0.0, dt=0.002), 1.0, 1.6, tol=1e-5))
    print("   threshold at J_AB=0.5: %.5f"
          % bisect(lambda J: peak(1.0, J, [0.05] * 3, J_AB=0.5, dt=0.002), 1.0, 1.6, tol=1e-5))
    print("   dominance crossover, J_SB=0: J_AB = %.4f"
          % bisect(lambda x: peak(1.0, 0.0, [0.05] * 3, J_AB=x, dt=0.002), 0.4, 1.2, tol=1e-4))
    print("   the shadow sweep as the committed code runs it, with WHERE the peak falls,")
    print("   which is the load-bearing step: an exact zero means the peak sits at t_measure.")
    for J_AB in (0.0, 0.1, 0.3, 0.5, 1.0):
        r = v3.sweep_direct_coupling([J_AB])[0]
        arg = []
        for meas in (None, 1.0):
            rec = v3.run_star_topology(2, J_SA=1.0, J_SB=1.466, J_AB=J_AB, gamma=0.05,
                                       dt=0.005, t_max=5.0, measure_a_at=meas)
            ts = [t for t in rec.t if t >= 1.0]
            vals = [v for t, v in zip(rec.t, rec.pairs["SB"].R) if t >= 1.0]
            arg.append(ts[int(np.argmax(vals))])
        print("      J_AB=%-5s suppression %-12s argmax of R_SB over t>=1: %.3f without, %.3f with"
              % (J_AB, str(r["shadow_relative_change_percent"]) + "%", arg[0], arg[1]))


def frequency():
    """Section 4.12: the closed form and the gamma shift."""
    H = v3.star_hamiltonian_n(2, J_SA=1.0, J_SB=2.0)
    closed = (1.0 + 2.0 + np.sqrt(1.0 - 2.0 + 4.0)) / np.pi
    ev = np.linalg.eigvalsh(H)
    print("   closed form at J=(1,2): %.6f;  max Bohr frequency / 2pi: %.6f"
          % (closed, (ev.max() - ev.min()) / (2 * np.pi)))
    def bohr(gam, dp=10):
        ev = np.linalg.eigvals(liouvillian(H, [gam] * 3, 3))
        w = np.sort(np.unique(np.round(np.abs(ev.imag), dp)))
        return w[w > 1e-9]

    a, b = bohr(0.0), bohr(0.005)
    print("   at gamma=0.005 the generator carries the dominant frequency both ways: %.8f and %.8f"
          % (a[-1], b[b < a[-1] - 1e-9].max()))
    print("      shift of the dominant pair: %.2e in f units"
          % ((b[b < a[-1] - 1e-9].max() - a[-1]) / (2 * np.pi)))
    print("      shift of the smallest:      %.2e in f units" % ((b[0] - a[0]) / (2 * np.pi)))
    print("   second order in gamma; the coefficient is per mode, so name the mode:")
    print("      %-9s %-14s %-14s" % ("gamma", "dominant", "smallest"))
    for g in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
        c = bohr(g)
        dom = (a[-1] - c[c < a[-1] - 1e-9].max()) / g ** 2
        print("      %-9s %-14.4f %-14.4f" % (g, dom, abs(c[0] - a[0]) / g ** 2))


def others():
    """Sections 4.2, 4.3, 4.5 and 4.11: the numbers that are not thresholds."""
    from star_topology_v2 import (star_hamiltonian, dephasing_ops, make_state,
                                  rk4_step, ptrace, concurrence, psi_norm)

    def run(state, J_SB, gam, t_max=5.0, dt=0.001):
        H = star_hamiltonian(1.0, J_SB)
        L = dephasing_ops(gam, gam, gam)
        rho = make_state(state)
        rsum, ab = [], []
        steps = int(t_max / dt)
        for s in range(steps + 1):
            sa, sb, p = ptrace(rho, 2), ptrace(rho, 1), ptrace(rho, 0)
            rsum.append(concurrence(sa) * psi_norm(sa) ** 2 + concurrence(sb) * psi_norm(sb) ** 2)
            ab.append(concurrence(p) * psi_norm(p))
            if s < steps:
                rho = rk4_step(rho, H, L, dt)
        return np.array(rsum), np.array(ab)

    r, _ = run("Bell_SA+B", 1.0, 0.05)
    print("   4.3  R_SA+R_SB: init %.4f, peak %.4f, ratio %.2fx" % (r[0], r.max(), r.max() / r[0]))
    r, _ = run("0++", 1.0, 0.05)
    print("   4.5  |0++> R_SA+R_SB: init %.4f, peak %.4f" % (r[0], r.max()))
    r, _ = run("W", 1.0, 0.05)
    print("   4.2  W state R_SA+R_SB: init %.4f, at t=5 %.4f" % (r[0], r[-1]))
    for t_max in (2.0, 5.0, 8.0):
        _, ab = run("0++", 10.0, 0.001, t_max=t_max)
        print("   4.11 |0++> at J_SB=10, gamma=0.001, t_max=%s: peak CPsi_AB = %.5f  %s"
              % (t_max, ab.max(), "CROSSES" if ab.max() >= 0.25 else "no"))

    # Section 5.10, the echo: AB alive while both SA and SB read zero.
    ts, ab = cpsi_series(2, 1.0, 2.0, [0.05] * 3, dt=0.002, t_max=5.0)
    _, sa = cpsi_series(2, 1.0, 2.0, [0.05] * 3, dt=0.002, t_max=5.0, pair=(0, 1))
    _, sb = cpsi_series(2, 1.0, 2.0, [0.05] * 3, dt=0.002, t_max=5.0, pair=(0, 2))
    i = int(round(2.2 / 0.002))
    print("   5.10 at t=2.2: CPsi_SA = %.6f, CPsi_SB = %.6f, CPsi_AB = %.6f"
          % (sa[i], sb[i], ab[i]))


BLOCKS = {"threshold": threshold, "scale": scale, "noise": noise, "curve": curve,
          "coarse": coarse, "observers": observers, "direct": direct, "frequency": frequency, "others": others}

if __name__ == "__main__":
    want = sys.argv[1] if len(sys.argv) > 1 else "all"
    if want != "all" and want not in BLOCKS:
        print("blocks: " + ", ".join(BLOCKS) + ", or all")
        sys.exit(1)
    for name in (list(BLOCKS) if want == "all" else [want]):
        print("\n=== %s: %s" % (name, BLOCKS[name].__doc__.splitlines()[0]))
        t0 = time.time()
        BLOCKS[name]()
        print("    [%.0f s]" % (time.time() - t0))
