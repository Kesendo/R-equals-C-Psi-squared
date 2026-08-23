"""Reading 3 of experiments/MEDIATOR_NOISE_GATE_LEVEL_THREE.md, plus the three
convergence checks the page's Reading 1 and Reading 2 rely on.

Reading 3 tests, and refutes, a mechanism.  F64 (docs/ANALYTICAL_FORMULAS.md)
says a single dephased site enters weighted by its own amplitude, |a_B|^2, and
is "not diminished by intervening sites".  The mediator on these chains is the
exact centre of a uniform open chain, so reflection symmetry about the centre
forces every antisymmetric single-excitation eigenvector to vanish there
exactly: the centre is the one site with nodes.  That is the only structural
difference between it and any other site, since the mean weight is exactly 1/N
everywhere by completeness.

The script computes those weights from the single-excitation block of the SAME
Heisenberg chain (the ZZ term shifts the weights but leaves the nodes exact),
and then puts the noise on each site in turn.  The result is a null: the site
with nodes is not anomalous, and the static weight is flat across the chain
while the measured response spans a factor of 8.6.

Also here, because the page states them and nothing else checks them:
  - integrator convergence of the 15.6 % span (dt refined 0.05 -> 0.00625)
  - window invariance of the same span (t_max = 20, 40, 80)
  - the dimensionless response at the operating point (log-log derivatives)
  - the unitary-window divergence that makes a quiet-background comparison
    impossible with a peak-over-window observable
"""

import sys
import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0])

from bridge_sector import (sector_basis, SectorPropagator, bell_on_vacuum,
                           mutual_information, mediator_bridge, build_h)

A3, B3 = [0, 1, 2, 3, 4], [6, 7, 8, 9, 10]
S3, I3 = sector_basis(11)
RHO3 = bell_on_vacuum(11, 0, 1, S3, I3)
BONDS3 = mediator_bridge(3)


def peak(j_meta=1.0, gamma_m=0.05, t_max=20.0, dt=0.05,
         keep_a=None, keep_b=None):
    keep_a = keep_a or A3
    keep_b = keep_b or B3
    g = [0.05] * 11
    g[5] = gamma_m
    p = SectorPropagator(11, mediator_bridge(3, j_meta=j_meta), g, S3, I3)
    best = [0.0, 0.0]

    def cb(t, rho):
        mi = mutual_information(rho, 11, S3, keep_a, keep_b)
        if mi > best[0]:
            best[0], best[1] = mi, t

    p.propagate_every_step(RHO3, t_max, dt, cb)
    return best[0], best[1]


def span(t_max=20.0, dt=0.05):
    hi = peak(gamma_m=0.0, t_max=t_max, dt=dt)[0]
    lo = peak(gamma_m=0.5, t_max=t_max, dt=dt)[0]
    return (hi - lo) / hi * 100


# ------------------------------------------------------------- Reading 3
def single_excitation_block(n):
    """The one-excitation block of the uniform Heisenberg chain, built through
    the same primitive the propagator uses."""
    states = [1 << (n - 1 - j) for j in range(n)]
    index = {s: i for i, s in enumerate(states)}
    bonds = [(i, i + 1, 1.0) for i in range(n - 1)]
    return build_h(n, bonds, states, index)


print("=" * 70)
print("READING 3: the mediator's amplitude weight (F64)")
print("=" * 70)
weights = {}
for n in (5, 11):
    h = single_excitation_block(n)
    w, v = np.linalg.eigh(h)
    med = (n - 1) // 2
    a2 = np.abs(v[med, :]) ** 2
    nodes = int((a2 < 1e-12).sum())
    nz = np.sort(np.unique(np.round(a2[a2 > 1e-12], 9)))[::-1]
    weights[n] = nz[0]
    print(f"\n  N={n}, mediator = site {med} (the exact centre)")
    print("    |a_M|^2 per mode: " + " ".join(f"{x:.4f}" for x in a2))
    print(f"    exact nodes at the mediator: {nodes} of {n}"
          f"   (floor(N/2) = {n // 2})")
    print(f"    distinct nonzero weights: "
          + ", ".join(f"{x:.6f}" for x in nz))
    print(f"    mean over all modes: {a2.mean():.6f}  (= 1/N = {1 / n:.6f})")

print(f"\n  NOTE: the two distributions differ by the single factor "
      f"{weights[11] / weights[5]:.4f} = 5/11,")
print("  and the mean weight is exactly 1/N at EVERY site of every chain")
print("  (completeness / N), so neither the mean nor the dominant weight")
print("  distinguishes the mediator from any other site.  Only the nodes do,")
print("  and the test below asks whether that structural difference matters.")

# ------------------- the test: put the noise on each site in turn
print("\n" + "=" * 70)
print("THE TEST: response to gamma on each site in turn (N=11)")
print("=" * 70)
_states = [1 << (11 - 1 - j) for j in range(11)]
_index = {s_: i for i, s_ in enumerate(_states)}
_h = build_h(11, [(i, i + 1, 1.0) for i in range(10)], _states, _index)
_, _v = np.linalg.eigh(_h)


def peak_site(site, g_site):
    g = [0.05] * 11
    g[site] = g_site
    p = SectorPropagator(11, BONDS3, g, S3, I3)
    best = [0.0]

    def cb(t, rho, b=best):
        mi = mutual_information(rho, 11, S3, A3, B3)
        if mi > b[0]:
            b[0] = mi

    p.propagate_every_step(RHO3, 20.0, 0.05, cb)
    return best[0]


print(f"  {'site':>5}  {'nodes':>6}  {'max|a|^2':>9}  {'span %':>8}")
print("  " + "-" * 34)
for j in range(11):
    a2 = np.abs(_v[j, :]) ** 2
    hi, lo = peak_site(j, 0.0), peak_site(j, 0.5)
    mark = "   <- the mediator" if j == 5 else ""
    print(f"  {j:>5}  {int((a2 < 1e-12).sum()):>6}  {a2.max():>9.4f}"
          f"  {(hi - lo) / hi * 100:>8.2f}{mark}")
print("\n  The one site with nodes is not anomalous: it sits between its")
print("  neighbours and above site 4.  Static weight is flat to 2 % across")
print("  the chain while the response spans a factor of 8.6, so the amplitude")
print("  reading cannot carry the response.  The mediator is an ordinary site.")

# ------------------------------------------------- convergence of Reading 1
print("\n" + "=" * 70)
print("READING 1 CONVERGENCE: integrator, and window")
print("=" * 70)
print(f"  {'RK4 dt':>10}  {'span %':>9}")
print("  " + "-" * 22)
for dt in (0.05, 0.025, 0.0125, 0.00625):
    print(f"  {dt:10.5f}  {span(dt=dt):9.3f}")
print(f"\n  {'t_max':>10}  {'span %':>9}")
print("  " + "-" * 22)
for tm in (20.0, 40.0, 80.0):
    print(f"  {tm:10.1f}  {span(t_max=tm):9.3f}")

# ------------------------------------- dimensionless response, Reading 1
print("\n" + "=" * 70)
print("READING 1: dimensionless response at the operating point")
print("=" * 70)
h = 0.05
s_j = (np.log(peak(j_meta=np.exp(h))[0])
       - np.log(peak(j_meta=np.exp(-h))[0])) / (2 * h)
s_g = (np.log(peak(gamma_m=0.05 * np.exp(h))[0])
       - np.log(peak(gamma_m=0.05 * np.exp(-h))[0])) / (2 * h)
print(f"  central difference in log space, h = {h}")
print(f"  d lnMI / d lnJ_meta   = {s_j:+.4f}")
print(f"  d lnMI / d lngamma_M  = {s_g:+.4f}")
print(f"  ratio = {abs(s_j / s_g):.1f} : 1")

# --------------------------- why the quiet-background comparison cannot be made
print("\n" + "=" * 70)
print("WHY A QUIET-BACKGROUND COMPARISON FAILS WITH THIS OBSERVABLE")
print("=" * 70)
print("  With every gamma zero the dynamics is unitary: I(A:B) never settles,")
print("  so 'peak over the window' is only 'largest sample so far'.")
print(f"  {'t_max':>8}  {'peak':>10}  {'t*':>8}")
print("  " + "-" * 30)
for tm in (20.0, 40.0, 80.0):
    p = SectorPropagator(11, BONDS3, [0.0] * 11, S3, I3)
    best = [0.0, 0.0]

    def cb(t, rho, b=best):
        mi = mutual_information(rho, 11, S3, A3, B3)
        if mi > b[0]:
            b[0], b[1] = mi, t

    p.propagate_every_step(RHO3, tm, 0.05, cb)
    print(f"  {tm:8.1f}  {best[0]:10.6f}  {best[1]:8.2f}")
print("  The reference arm is still climbing at the window edge.")

# ------------- Reading 2: the far sweep is a property of the window
print("\n" + "=" * 70)
print("READING 2: the large-gamma_M rows move with the window")
print("=" * 70)


def peak_wide(gamma_m, t_max):
    g = [0.05] * 11
    g[5] = gamma_m
    p = SectorPropagator(11, BONDS3, g, S3, I3)
    dt = p.stable_dt(0.05)
    best = [0.0, 0.0]

    def cb(t, rho, b=best):
        mi = mutual_information(rho, 11, S3, A3, B3)
        if mi > b[0]:
            b[0], b[1] = mi, t

    p.propagate_every_step(RHO3, t_max, dt, cb)
    return best[0], best[1], dt


ref = peak_wide(0.0, 80.0)[0]
print(f"  gamma_M = 0 is window-stable at {ref:.6f}; half of it is {ref / 2:.4f}")
print(f"  {'gamma_M':>8}  {'t_max=20':>9}  {'t_max=80':>9}  {'t* (80)':>8}  {'dt':>9}")
print("  " + "-" * 52)
for gm in (5, 10, 20, 50, 100, 200):
    a = peak_wide(gm, 20.0)
    b = peak_wide(gm, 80.0)
    print(f"  {gm:>8}  {a[0]:9.4f}  {b[0]:9.4f}  {b[1]:8.2f}  {b[2]:9.5f}")
print("  Every t* in the t_max=80 column is at the window edge, so even that")
print("  column is not converged.  No closure threshold follows from this")
print("  observable; only the ordering at a fixed window does.")
