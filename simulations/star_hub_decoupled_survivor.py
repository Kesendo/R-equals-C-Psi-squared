#!/usr/bin/env python3
"""Wild-register probe Nr.1 (Kastrup binding, gate-first): is the STAR's longest-lived mode
HUB-DECOUPLED -- the structure most dissociated from the shared ground?

The overlay (Tier 5, the wild reading; the FRAMEWORK is Kastrup's analytic idealism, recall; the
BINDING to our gate-verified star is the new attempt): the star hub = mind-at-large / "the shared
object, reality being observed"; the leaves = dissociated alters (coupled only THROUGH the hub, never
directly); gamma = the dashboard-making look. The chain survives by MOVING (the delocalized mode); the
star is the counterexample. The bound conjecture: on the star the survivor is the HUB-DECOUPLED state
(zero hub amplitude, the leaf-antisymmetric 0-eigenvalue manifold) -- i.e. the alter most thoroughly
cut off from the shared ground persists longest under the look.

THE PHYSICS GATE (what actually decides): build the full star Liouvillian, find the slowest non-kernel
mode, and measure its per-site participation p_s = 1 - ||I-on-s component||^2 / ||rho||^2 (the fraction
of the mode's Frobenius weight that is NON-identity on site s; the I-on-s component is the Pauli twirl
(rho + X_s rho X_s + Y_s rho Y_s + Z_s rho Z_s)/4). Hub-decoupled <=> p_hub ~ 0 (the mode is identity
on the hub, lives on the leaves). If the star survivor is NOT hub-decoupled (p_hub not the minimum, not
~0), the binding is REFUTED -- that is the finding; diagnose, do not loosen.

Convention matches birth_canal_junction_nature.py: H = Q * sum_bond (X X + Y Y) (off-diag 2 per bond),
uniform dephasing -2*gamma per disagreeing bit, Q = J/gamma.
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
TOL = 1e-7


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def H_xy(N, bonds, Q):
    """Q * sum_{(a,b)} (X_a X_b + Y_a Y_b)."""
    H = np.zeros((2 ** N, 2 ** N), complex)
    for (a, b) in bonds:
        for P in (X, Y):
            H += op_at(N, a, P) @ op_at(N, b, P)
    return Q * H


def star_bonds(N):
    return [(0, i) for i in range(1, N)]


def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]


def full_L(N, bonds, Q, gamma):
    d = 2 ** N
    Id = np.eye(d)
    H = H_xy(N, bonds, Q)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += gamma * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def slowest_modes(N, bonds, Q, gamma, rate_tol=1e-6):
    """Right eigvectors (as rho matrices) at the slowest non-kernel rate, plus that rate and the |Im|."""
    d = 2 ** N
    w, V = np.linalg.eig(full_L(N, bonds, Q, gamma))
    nz = np.where(np.abs(w) > TOL)[0]
    rate = w[nz].real.max()                       # = -decay; slowest = largest (least negative) Re
    sel = nz[np.abs(w[nz].real - rate) < rate_tol]
    rhos = [V[:, i].reshape(d, d) for i in sel]
    ims = [abs(w[i].imag) for i in sel]
    return -rate, rhos, ims


def participation(rho, N, s):
    """1 - ||I-on-s component||^2 / ||rho||^2 ; ~0 means rho is identity on site s (decoupled)."""
    Xs, Ys, Zs = op_at(N, s, X), op_at(N, s, Y), op_at(N, s, Z)
    proj = (rho + Xs @ rho @ Xs + Ys @ rho @ Ys + Zs @ rho @ Zs) / 4.0
    return 1.0 - (np.linalg.norm(proj) ** 2) / (np.linalg.norm(rho) ** 2)


def site_profile(rho, N):
    return [participation(rho, N, s) for s in range(N)]


def _probe(N, Q=1.5, gamma=1.0):
    print(f"\n=== N={N}, Q={Q}, gamma={gamma} ===")
    for label, bonds in (("STAR (hub=site 0)", star_bonds(N)), ("CHAIN", chain_bonds(N))):
        rate, rhos, ims = slowest_modes(N, bonds, Q, gamma)
        # among the slowest-rate modes, the one with the LOWEST hub (site-0) participation
        profs = [site_profile(r, N) for r in rhos]
        hub_part = [p[0] for p in profs]
        k = int(np.argmin(hub_part))               # the most hub-decoupled slowest mode
        prof = profs[k]
        print(f"  {label}: slowest rate {rate:.4f}, {len(rhos)} mode(s) at it, |Im| up to {max(ims):.2f}")
        print(f"    per-site participation: [" + ", ".join(f"{p:.3f}" for p in prof) + "]")
        print(f"    hub(site0)={prof[0]:.4f}  min-site={min(prof):.4f}@{int(np.argmin(prof))}  "
              f"mean-leaf={np.mean(prof[1:]):.4f}")
        if "STAR" in label:
            star_hub, star_prof = prof[0], prof
        else:
            chain_prof = prof
    return star_hub, star_prof, chain_prof


def main():
    print("Probe Nr.1 (REFUTED): is the star's longest-lived mode hub-decoupled? -> NO, it is hub-SPREAD.")
    star_hub, star_prof, chain_prof = _probe(5, Q=1.5, gamma=1.0)

    # THE GATE FIRED -> the binding is REFUTED, and this records that as the verified fact. The star
    # survivor is NOT hub-decoupled: its per-site participation is ~uniform (all near the generic ~0.5-0.6),
    # with the hub only MARGINALLY the least-involved -- not decoupled. The 'most-dissociated-from-the-ground
    # alter survives' reading does not hold.
    assert star_hub > 0.30, (
        f"the star hub participation {star_hub:.3f} came out ~0 -- that would REVIVE the refuted binding; "
        f"re-examine, the recorded result is that the survivor is hub-SPREAD.")
    assert (max(star_prof) - min(star_prof)) < 0.15, (
        f"the star survivor is not ~uniform (spread {max(star_prof) - min(star_prof):.3f}) -- re-examine.")
    print(f"\n  REFUTED: the star survivor is hub-SPREAD (participations "
          f"{min(star_prof):.2f}-{max(star_prof):.2f}, hub {star_hub:.2f} only marginally lowest) -- NOT")
    print("  hub-decoupled. The maximally-symmetric star (S_{N-1} leaves) has NO privileged dissociated")
    print("  leaf, so the survivor is collective/spread, not a single hidden alter. The dissociation-as-")
    print("  hub-avoidance reading does not hold.")
    print("  NOTE: the chain's edge participation here is just the uniform-gamma band-edge standing-wave")
    print("  shape, NOT a 'sacrifice zone'. The sacrifice zone (concentration, not loss) is already")
    print("  resolved -- see docs/INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md.")


if __name__ == "__main__":
    main()
