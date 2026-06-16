"""Physics-first review verifier for the TOPOLOGY BAND EDGE spec
(docs/superpowers/specs/2026-06-16-topology-band-edge-design.md).

Independent ground truth for:
  (LAW)  band edge = J * adjacency spectral radius  (chain 2cos(pi/(N+1)), star sqrt(N-1), ring 2)
  (MAP)  star N>=5 "NEVER the gap mode at ANY Q"; ring N=4 "NEVER" -- the load-bearing novel claims.

The attack: the |vac><psi_k| single-excitation coherence sits at Re = -2g EXACTLY for ANY
topology (L_D = -2g*I on that whole space, n_XY=1). So an oscillating band-edge mode ALWAYS
exists at the 2g floor. The question is only whether a DIFFERENT mode is a hair SLOWER than 2g,
stealing the strict-gap label so the clock reads omega_mem=0. We therefore print, per (topo,N,Q):
  - the strict gap and clock omega_mem (tol 1e-6, the witness convention)
  - whether a band-edge oscillation exists AT the 2g floor (|rate-2g| < 1e-4)
  - the few slowest modes (rate/2g, |Im|/J) so we can SEE what holds the floor
  - the Q-dependence (does gap/2g -> 1 as Q grows? => horizon, contradicts "structural")
"""
import sys
from math import cos, pi, sqrt
import numpy as np

sys.path.insert(0, 'simulations')
import framework as fw

GAMMA = 0.05
TOLZ = 1e-7
CLOCK_TOL = 1e-6          # the clock witness's rate-==-gap tolerance (Symphony.cs:165)


def adjacency_radius(N, bonds):
    A = np.zeros((N, N))
    for (i, j) in bonds:
        A[i, j] = 1.0
        A[j, i] = 1.0
    return max(abs(np.linalg.eigvalsh(A)))


def closed_form_rho(topo, N):
    if topo == 'chain':
        return 2 * cos(pi / (N + 1))
    if topo == 'star':
        return sqrt(N - 1)
    if topo == 'ring':
        return 2.0
    raise ValueError(topo)


def se_coherence_eigs(chain):
    """Eigenvalues of the |vac><single-excitation| coherence block, built directly:
    lambda = -2*gamma + i*E_k, E_k = single-particle (J*adjacency) spectrum.
    Confirms these modes sit at Re = -2g EXACTLY, any topology."""
    N, J = chain.N, chain.J
    A = np.zeros((N, N))
    for (i, j) in chain.bonds:
        A[i, j] = A[j, i] = 1.0
    Ek = np.linalg.eigvalsh(J * A)
    return -2 * GAMMA + 1j * Ek      # the vac<->SE coherence line


def analyze(topo, N, Q):
    J = Q * GAMMA
    chain = fw.ChainSystem(N=N, gamma_0=GAMMA, J=J, topology=topo, H_type='xy')
    lam = np.linalg.eigvals(chain.L)
    re, im = lam.real, lam.imag

    dec = re < -TOLZ
    rate = -re                                   # decay rate (positive)
    gap = rate[dec].min()                        # strict slowest
    g2 = gap / (2 * GAMMA)

    # clock omega_mem: max|Im| among modes at the STRICT gap (witness convention)
    at_gap = dec & (np.abs(rate - gap) < CLOCK_TOL)
    omega_strict = np.abs(im[at_gap]).max() if at_gap.any() else 0.0

    # is there a band-edge oscillation AT the 2g floor (even if not the strict gap)?
    at_floor = np.abs(rate - 2 * GAMMA) < 1e-4 * (2 * GAMMA)
    omega_at_floor = np.abs(im[at_floor]).max() if at_floor.any() else 0.0

    band = J * closed_form_rho(topo, N)

    # the 6 slowest decaying modes
    order = np.argsort(rate[dec])
    rr = rate[dec][order][:6]
    ii = np.abs(im[dec][order][:6])
    slowest = list(zip(rr / (2 * GAMMA), ii / J))

    return dict(J=J, band=band, g2=g2,
                omega_strict=omega_strict, omega_strict_over_band=omega_strict / band if band else 0,
                omega_at_floor=omega_at_floor, omega_floor_over_band=omega_at_floor / band if band else 0,
                slowest=slowest)


# ---------- LAW ----------
print("=" * 100)
print("LAW: J*rho(adjacency) == closed form == max|Im| of vac<->SE coherence line")
print("=" * 100)
for topo in ('chain', 'star', 'ring'):
    for N in (3, 4, 5, 6):
        ch = fw.ChainSystem(N=N, gamma_0=GAMMA, J=1.0, topology=topo, H_type='xy')
        rho = adjacency_radius(N, ch.bonds)
        cf = closed_form_rho(topo, N)
        se = se_coherence_eigs(ch)
        se_maxim = max(abs(se.imag))
        print(f"  {topo:5s} N={N}: rho={rho:.6f}  closed={cf:.6f}  d={abs(rho-cf):.1e}   "
              f"max|Im(SE-line)|={se_maxim:.6f}  Re(SE-line)={se.real[0]:+.4f} (==-2g={-2*GAMMA})")

# ---------- MAP ----------
print("\n" + "=" * 100)
print("MAP: is the band edge the GAP MODE?  (g2=gap/2g; w_strict=clock omega_mem at strict gap;")
print("     w_floor=max|Im| of any mode sitting AT the 2g floor; both / band edge)")
print("=" * 100)
for topo in ('chain', 'star', 'ring'):
    print(f"\n########## {topo} ##########")
    for N in (3, 4, 5):
        print(f"  --- N={N}  (band edge J*rho, closed rho={closed_form_rho(topo,N):.4f}) ---")
        print("    Q       g2       w_strict/band   w_floor/band   slowest 6 modes (rate/2g, |Im|/J)")
        for Q in (1.0, 5.0, 20.0, 100.0, 1000.0):
            r = analyze(topo, N, Q)
            sl = "  ".join(f"({a:.4f},{b:.3f})" for a, b in r['slowest'])
            print(f"    {Q:<7.0f} {r['g2']:.5f}   {r['omega_strict_over_band']:.4f}          "
                  f"{r['omega_floor_over_band']:.4f}         {sl}")

# Targeted N=6 spot-checks (4096^2 dense; star ceiling + ring even-fill, high Q only)
print("\n########## N=6 spot-checks (the 'NEVER at any Q' stress) ##########")
print("    topo   Q       g2       w_strict/band   w_floor/band   slowest 6 (rate/2g, |Im|/J)")
for topo in ('star', 'ring'):
    for Q in (20.0, 1000.0):
        r = analyze(topo, 6, Q)
        sl = "  ".join(f"({a:.4f},{b:.3f})" for a, b in r['slowest'])
        print(f"    {topo:5s}  {Q:<7.0f} {r['g2']:.5f}   {r['omega_strict_over_band']:.4f}          "
              f"{r['omega_floor_over_band']:.4f}         {sl}")

print("\nDONE.")
