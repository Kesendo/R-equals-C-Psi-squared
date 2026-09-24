"""
The February pair tables in their own book: the pairwise correlation bridge
===========================================================================
The February 2026 pair tables of SUBSYSTEM_CROSSING, DYNAMIC_ENTANGLEMENT §5-§6 and
SIMULATION_EVIDENCE §7.2 are QuTiP `mesolve` runs, and their pair readings use a pairwise
correlation bridge,

    C_corr = clip((P_AB - P_A*P_B) / (1 - P_A*P_B), 0, 1),

times Psi = l1/3 of the pair's reduced density matrix. The retired delta_calc MCP tool
carries the same formula in its subsystem-crossing routine (its source is kept outside the
repo and is transcribed here, not imported), and the tool's own runs of the |0+0+> trajectory
are DYNAMIC_ENTANGLEMENT §9 and §11.2. This script recomputes both, with two propagators:

  exact propagation (the matrix exponential), for the QuTiP tables:
    SUBSYSTEM_CROSSING   §3.2 W4 pair at t = 0
                         §3.3 Bell+ x Bell+ pair (0,1): downward crossing at t = 0.073
    SIMULATION_EVIDENCE  §7.2 |0+0+> pair (0,2) at gamma = 0.05: the fifteen-entry table
    DYNAMIC_ENTANGLEMENT §5.1 unitary and §5.2 gamma = 0.05 tables, §6 pair diagonal
  the tool's Euler loop (dt 0.01, Hermitian part, negative eigenvalues clipped, trace
  renormalized, t_max as in each call), for DYNAMIC_ENTANGLEMENT §9.2 (the tool's Ising,
  an open chain), §9.1's N = 4 column, the §9.3 gamma sweep and §11.2's Bell-pair line,
  and without the clipping, as the routine's first version ran, for that version's
  Bell-pair crossing and runaway maximum;
and it checks the claims the pages make about these runs under exact propagation: the
unitary period, the §6 population bound, what §9's numbers become without the Euler step,
this script's own sigma_x run, which tests §5.3's prediction that sigma_x dephasing moves
the surviving pair to (1,3), and the Bell+ x Bell+ numbers of
hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md (tool run and exact).

Conventions (the tool's and the tables'): H = J * sum over ring bonds of (XX + YY + ZZ) with
J = 1 in the Pauli book, N = 4, qubit 0 the most significant bit; local jumps Z_k at rate
gamma, D[rho] = gamma * sum_k (Z_k rho Z_k - rho) (the Lindblad book); the sigma_x run swaps
Z_k for X_k. The canonical pair reading of the repo is a different book (Wootters
concurrence * l1/3, simulations/subsystem_crossing_pairs.py).

Import-inert; prints only.
"""

from itertools import combinations
import sys

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar

N = 4
J = 1.0
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
ZERO = np.array([1, 0], dtype=complex)
PLUS = np.array([1, 1], dtype=complex) / np.sqrt(2)
PAIRS = list(combinations(range(N), 2))


def site(P, k):
    m = np.array([[1.0]], dtype=complex)
    for i in range(N):
        m = np.kron(m, P if i == k else I2)
    return m


def heisenberg_ring():
    H = np.zeros((2**N, 2**N), dtype=complex)
    for i in range(N):
        j = (i + 1) % N
        for P in (X, Y, Z):
            H += J * site(P, i) @ site(P, j)
    return H


def ising_chain():
    """The tool's 'ising': J * sum over the open chain's bonds of ZZ."""
    return sum(J * site(Z, i) @ site(Z, i + 1) for i in range(N - 1))


def liouvillian(gamma, jump=Z, H=None):
    d = 2**N
    Id = np.eye(d)
    H = heisenberg_ring() if H is None else H
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Lk = site(jump, k)
        L += gamma * (np.kron(Lk, Lk.conj()) - np.kron(Id, Id))
    return L


def reduced(rho, keep):
    t = rho.reshape([2] * (2 * N))
    for i in sorted((q for q in range(N) if q not in keep), reverse=True):
        t = np.trace(t, axis1=i, axis2=i + t.ndim // 2)
    k = len(keep)
    return t.reshape(2**k, 2**k)


def purity(r):
    return float(np.real(np.trace(r @ r)))


def pair_reading(rho, i, j):
    """(C_corr, Psi, CPsi, diagonal) of the pair (i, j) in the pairwise bridge."""
    rab = reduced(rho, [i, j])
    pa, pb, pab = purity(reduced(rho, [i])), purity(reduced(rho, [j])), purity(rab)
    den = 1.0 - pa * pb
    c = min(1.0, max(0.0, (pab - pa * pb) / den)) if den > 1e-12 else 0.0
    l1 = float(np.abs(rab).sum() - np.abs(np.diag(rab)).sum())
    return c, l1 / 3.0, c * l1 / 3.0, np.real(np.diag(rab))


def product(*kets):
    psi = kets[0]
    for k in kets[1:]:
        psi = np.kron(psi, k)
    return np.outer(psi, psi.conj())


def evolve(L, rho0, t):
    return (expm(L * t) @ rho0.ravel()).reshape(2**N, 2**N)


def crossing(L, rho0, pair, lo, hi):
    """Bisect the time where the pair's CPsi passes 1/4 between lo and hi."""
    f = lambda t: pair_reading(evolve(L, rho0, t), *pair)[2] - 0.25
    flo = f(lo)
    for _ in range(60):
        mid = (lo + hi) / 2
        if (f(mid) > 0) == (flo > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def passages(series, step):
    """The tool's detector on each pair: first upward and first downward passage of 1/4,
    linearly interpolated between grid points, the number of upward passages, the maximum
    and the time of the maximum."""
    out = {}
    for n, p in enumerate(PAIRS):
        c = series[:, n]
        ups = np.where((c[:-1] < 0.25) & (c[1:] >= 0.25))[0]
        downs = np.where((c[:-1] >= 0.25) & (c[1:] < 0.25))[0]
        up = (ups[0] + (0.25 - c[ups[0]]) / (c[ups[0] + 1] - c[ups[0]])) * step if len(ups) else None
        down = (downs[0] + (c[downs[0]] - 0.25) / (c[downs[0]] - c[downs[0] + 1])) * step if len(downs) else None
        k = int(np.argmax(c))
        out[p] = {"up": up, "down": down, "n_up": len(ups), "max": float(c[k]), "t_max": k * step}
    return out


def first_peak(L, rho0, pair, lo, hi):
    """Maximum of the pair's CPsi on (lo, hi), located by a bounded scalar search."""
    f = lambda t: -pair_reading(evolve(L, rho0, t), *pair)[2]
    return -minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-10}).fun


def pairs_crossing(res):
    return [p for p, r in res.items() if r["up"] is not None or r["down"] is not None]


def exact_series(L, rho0, t_max, step=0.001):
    U = expm(L * step)
    v = rho0.ravel().copy()
    n = int(round(t_max / step))
    series = np.zeros((n + 1, len(PAIRS)))
    for s in range(n + 1):
        rho = v.reshape(2**N, 2**N)
        series[s] = [pair_reading(rho, *p)[2] for p in PAIRS]
        v = U @ v
    return series


def euler_series(rho0, gamma, t_max, H=None, dt=0.01, clip=True):
    """The tool's loop, transcribed: an Euler step of the Lindblad equation, then the
    Hermitian part, negative eigenvalues set to zero, and the trace renormalized. With
    clip=False, the plain Euler step of the routine's first version."""
    H = heisenberg_ring() if H is None else H
    jumps = [site(Z, k) for k in range(N)]
    n = int(t_max / dt)
    rho = rho0.copy()
    series = np.zeros((n + 1, len(PAIRS)))
    for s in range(n + 1):
        series[s] = [pair_reading(rho, *p)[2] for p in PAIRS]
        if s < n:
            drho = -1j * (H @ rho - rho @ H)
            for Lk in jumps:
                LdL = Lk.conj().T @ Lk
                drho += gamma * (Lk @ rho @ Lk.conj().T - 0.5 * (LdL @ rho + rho @ LdL))
            rho = rho + dt * drho
            if clip:
                rho = 0.5 * (rho + rho.conj().T)
                w, v = np.linalg.eigh(rho)
                rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
                rho /= np.real(np.trace(rho))
    return series


VALUES = []  # printed page values, compared at the page's precision
CLAIMS = []  # the pages' statements about these runs


def check(label, page, value, digits=3):
    ok = round(value, digits) == round(page, digits)
    print(f"  {label:<48s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")
    VALUES.append((label, ok))


def check_count(label, page, value):
    ok = value == page
    print(f"  {label:<48s} page {page}  computed {value}  {'match' if ok else 'DIFFERS'}")
    VALUES.append((label, ok))


def claim(label, ok, detail):
    print(f"  {label:<48s} {detail}  {'holds' if ok else 'FAILS'}")
    CLAIMS.append((label, ok))


def main():
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    rho_a = product(ZERO, PLUS, ZERO, PLUS)
    L05 = liouvillian(0.05)

    print("SUBSYSTEM_CROSSING §3.2, W4 pair at t = 0 (no propagator enters)")
    w = np.zeros(2**N, dtype=complex)
    for k in range(N):
        w[1 << k] = 0.5
    c, psi, cpsi, _ = pair_reading(np.outer(w, w.conj()), 0, 1)
    check("C_corr(0), the February table's 0.180", 0.180, c)
    check("CPsi(0)", 0.030, cpsi)
    print(f"  (C_corr(0) = 7/39 = {7 / 39:.6f} exactly: P_AB = 1/2, P_A = P_B = 5/8)")

    print("SUBSYSTEM_CROSSING §3.3, Bell+ x Bell+ pair (0,1), gamma = 0.05")
    bell = np.zeros(4, dtype=complex)
    bell[0] = bell[3] = 1 / np.sqrt(2)
    rho_bb = np.outer(np.kron(bell, bell), np.kron(bell, bell).conj())
    check("downward crossing t", 0.073, crossing(L05, rho_bb, (0, 1), 0.0, 0.3))

    print("SIMULATION_EVIDENCE §7.2, |0+0+> pair (0,2), gamma = 0.05")
    table = {0.1: (0.049, 0.301, 0.015), 0.2: (0.187, 0.595, 0.111), 0.286: (0.327, 0.768, 0.251),
             0.35: (0.380, 0.834, 0.316), 0.5: (0.240, 0.787, 0.189)}
    for t, (pc, pp, pcp) in table.items():
        c, psi, cpsi, _ = pair_reading(evolve(L05, rho_a, t), 0, 2)
        check(f"t = {t}: C_corr", pc, c)
        check(f"t = {t}: Psi", pp, psi)
        check(f"t = {t}: CPsi", pcp, cpsi)

    print("DYNAMIC_ENTANGLEMENT §5.1 (unitary, t_max 3) and §5.2 (gamma = 0.05, t_max 3), exact, grid 0.001")
    L0 = liouvillian(0.0)
    r0 = passages(exact_series(L0, rho_a, 3.0), 0.001)
    for pair, (pt, pm) in {(0, 1): (0.073, 0.285), (0, 2): (0.260, 0.435), (1, 3): (0.679, 0.341)}.items():
        up = r0[pair]["up"]
        check(f"unitary {pair}: first crossing", pt, crossing(L0, rho_a, pair, up - 0.001, up + 0.001))
        check(f"unitary {pair}: max CPsi", pm, r0[pair]["max"])
    check_count("unitary: pairs crossing", 6, len(pairs_crossing(r0)))
    r5 = passages(exact_series(L05, rho_a, 3.0), 0.001)
    for pair, pm in {(0, 1): 0.247, (0, 2): 0.320, (1, 3): 0.224}.items():
        check(f"gamma 0.05 {pair}: max CPsi", pm, r5[pair]["max"])
    check_count("gamma 0.05: pairs crossing", 1, len(pairs_crossing(r5)))
    t_cross = crossing(L05, rho_a, (0, 2), 0.2, 0.3)
    check("gamma 0.05 (0,2): crossing t", 0.285, t_cross)
    energies = sorted(set(np.round(np.linalg.eigvalsh(heisenberg_ring()), 9)))
    claim("unitary period pi/2: energies differ by 4s", all(abs(e / 4 - round(e / 4)) < 1e-9 for e in energies),
          f"energies {[float(e) for e in energies]}")
    rp = passages(exact_series(L0, rho_a, np.pi / 2 - 0.001), 0.001)
    per_period = [rp[p]["n_up"] for p in ((0, 1), (0, 2), (1, 3))]
    claim("upward passages per period: (0,1) 2, (0,2) 2, (1,3) 1", per_period == [2, 2, 1],
          f"computed {per_period}")
    print(f"  (the (0,2) maximum at gamma 0.05 sits at t = {r5[(0, 2)]['t_max']:.3f})")

    print("DYNAMIC_ENTANGLEMENT §6, pair (0,2) diagonal at t = 0.286, gamma = 0.05")
    diag = pair_reading(evolve(L05, rho_a, 0.286), 0, 2)[3]
    for label, page, value in zip(("|00>", "|01>", "|10>", "|11>"), (0.425, 0.257, 0.257, 0.061), diag):
        check(label, page, value)
    claim("t = 0.286 lies past the crossing", 0.286 > t_cross, f"crossing at t = {t_cross:.6f}")
    # Heisenberg exchange and Z jumps conserve the excitation number; |0+0+> weighs 1/4, 1/2, 1/4
    # on the sectors n = 0, 1, 2. |00> on (0,2) contains the whole n = 0 sector, |0000>, and |11>
    # on (0,2) lies inside n = 2, so P(|00>) >= 1/4 >= P(|11>) at every t.
    U = expm(L05 * 0.01)
    v = rho_a.ravel().copy()
    p00, p11 = 1.0, 0.0
    for _ in range(2001):
        d = pair_reading(v.reshape(2**N, 2**N), 0, 2)[3]
        p00, p11 = min(p00, d[0]), max(p11, d[3])
        v = U @ v
    claim("P(|00>) >= 1/4 >= P(|11>) over t <= 20", p00 >= 0.25 >= p11,
          f"min P(|00>) = {p00:.6f}, max P(|11>) = {p11:.6f}")

    print("DYNAMIC_ENTANGLEMENT §5.3's prediction, tested: sigma_x dephasing, gamma = 0.05, exact, t_max 5")
    rx = passages(exact_series(liouvillian(0.05, X), rho_a, 5.0), 0.001)
    claim("sigma_x: only (0,2) crosses, not (1,3)", pairs_crossing(rx) == [(0, 2)],
          f"pairs crossing {pairs_crossing(rx)}")
    check("sigma_x (0,2): max CPsi", 0.335, rx[(0, 2)]["max"])
    check("sigma_x (1,3): max CPsi", 0.240, rx[(1, 3)]["max"])
    check("sigma_x ring neighbour (0,1): max CPsi", 0.249, rx[(0, 1)]["max"])
    rz = passages(exact_series(L05, rho_a, 5.0), 0.001)
    claim("sigma_x raises every pair's maximum", all(rx[p]["max"] > rz[p]["max"] for p in PAIRS),
          "sigma_x over sigma_z: " + ", ".join(f"{rx[p]['max'] - rz[p]['max']:+.4f}" for p in PAIRS))
    rz04 = passages(exact_series(liouvillian(0.04), rho_a, 8.0), 0.001)
    rx04 = passages(exact_series(liouvillian(0.04, X), rho_a, 8.0), 0.001)
    claim("gamma 0.04: sigma_z five pairs, sigma_x adds (1,3)",
          pairs_crossing(rz04) == [p for p in PAIRS if p != (1, 3)] and pairs_crossing(rx04) == PAIRS,
          f"sigma_z {len(pairs_crossing(rz04))} without (1,3): {(1, 3) not in pairs_crossing(rz04)},"
          f" sigma_x {len(pairs_crossing(rx04))}")

    print("DYNAMIC_ENTANGLEMENT §9, the tool's Euler loop, dt 0.01")
    sweep = {0.01: (10.0, 6, 0.260, 0.408, 13), 0.05: (5.0, 5, 0.275, 0.339, 1),
             0.10: (5.0, 1, 0.304, 0.276, 1), 0.20: (5.0, 0, None, 0.193, 0)}
    for gamma, (t_max, n_page, up_page, max_page, osc_page) in sweep.items():
        re = passages(euler_series(rho_a, gamma, t_max), 0.01)
        check_count(f"gamma {gamma:.2f}: pairs crossing (of 6)", n_page, len(pairs_crossing(re)))
        if up_page is not None:
            check(f"gamma {gamma:.2f} (0,2): t_cross_up", up_page, re[(0, 2)]["up"])
        check(f"gamma {gamma:.2f} (0,2): max CPsi", max_page, re[(0, 2)]["max"])
        check_count(f"gamma {gamma:.2f} (0,2): upward passages", osc_page, re[(0, 2)]["n_up"])
        if gamma == 0.05:
            check("gamma 0.05 ring neighbour (0,1): max CPsi", 0.251, re[(0, 1)]["max"])
            check("gamma 0.05 (1,3): max CPsi", 0.234, re[(1, 3)]["max"])
            check("gamma 0.05 (0,2): t_cross_down (§11.2)", 0.470, re[(0, 2)]["down"])
    ri = passages(euler_series(rho_a, 0.05, 5.0, H=ising_chain()), 0.01)
    check("Ising (open chain) (1,3): max CPsi (§9.2)", 0.068, ri[(1, 3)]["max"])
    rb = passages(euler_series(rho_bb, 0.05, 5.0), 0.01)
    check("Bell+ x Bell+ (0,1): t_cross_down (§11.2)", 0.072, rb[(0, 1)]["down"])
    check("MEDIATOR §1.3, tool run: neighbours (0,3) max", 0.13, rb[(0, 3)]["max"], digits=2)
    check("MEDIATOR §1.3, tool run: diagonals (0,2) max", 0.11, rb[(0, 2)]["max"], digits=2)
    ru = passages(euler_series(rho_bb, 0.05, 5.0, clip=False), 0.01)
    check("first, unclipped version: (0,1) t_cross_down", 0.077, ru[(0, 1)]["down"])
    check("first, unclipped version: max CPsi over pairs", 2.016, max(r["max"] for r in ru.values()))

    print("MEDIATOR_AS_QUANTUM_TRANSISTOR §1.3 and §4.2, Bell+ x Bell+ under exact propagation, gamma = 0.05")
    sb = exact_series(L05, rho_bb, 5.0)
    xb = passages(sb, 0.001)
    check("neighbours (0,3): max CPsi, exact", 0.132, xb[(0, 3)]["max"])
    check("diagonals (0,2): max CPsi, exact", 0.101, xb[(0, 2)]["max"])
    for t, page in ((0.8, 0.034), (1.6, 0.136), (2.4, 0.013)):
        check(f"(0,1) at t = {t}", page, sb[int(round(t / 0.001)), PAIRS.index((0, 1))])
    c01 = sb[:, PAIRS.index((0, 1))]
    k = int(np.ceil(xb[(0, 1)]["down"] / 0.001))
    while k + 1 < len(c01) and c01[k + 1] <= c01[k]:
        k += 1  # the first minimum after the crossing; the echoes come after it
    after = c01[k:]
    check("(0,1): largest echo after the crossing", 0.143, float(after.max()))
    check("(0,1): time of the largest echo", 1.57, (k + int(np.argmax(after))) * 0.001, digits=2)
    claim("(0,1) stays below 1/4 after its crossing", float(after.max()) < 0.25,
          f"first minimum at t = {k * 0.001:.3f}, largest echo {float(after.max()):.4f}"
          f" at t = {(k + int(np.argmax(after))) * 0.001:.3f}")

    print("The same §9 cells under exact propagation")
    x01 = passages(exact_series(liouvillian(0.01), rho_a, 10.0), 0.001)
    claim("gamma 0.01 (0,2): 5 upward passages in t <= 10", x01[(0, 2)]["n_up"] == 5,
          f"computed {x01[(0, 2)]['n_up']}")
    check("gamma 0.01 (0,2): t_cross_up, exact", 0.264, x01[(0, 2)]["up"])
    check("gamma 0.01 (0,2): max CPsi, exact", 0.407, x01[(0, 2)]["max"])
    x00 = passages(exact_series(L0, rho_a, 10.0), 0.001)
    claim("unitary (0,2): 13 upward passages in t <= 10", x00[(0, 2)]["n_up"] == 13,
          f"computed {x00[(0, 2)]['n_up']}")
    x09 = passages(exact_series(liouvillian(0.09), rho_a, 5.0), 0.001)
    x10 = passages(exact_series(liouvillian(0.10), rho_a, 5.0), 0.001)
    claim("gamma_c of (0,2) between 0.09 and 0.10", x09[(0, 2)]["max"] >= 0.25 > x10[(0, 2)]["max"],
          f"max (0,2) {x09[(0, 2)]['max']:.4f} at 0.09, {x10[(0, 2)]['max']:.4f} at 0.10")
    check("gamma 0.09 (0,2): max CPsi, exact", 0.258, x09[(0, 2)]["max"])
    gamma_c = brentq(lambda g: first_peak(liouvillian(g), rho_a, (0, 2), 0.2, 0.6) - 0.25, 0.09, 0.10, xtol=1e-9)
    check("gamma_c of (0,2), exact", 0.096, gamma_c)
    claim("gamma 0.10: no pair crosses", pairs_crossing(x10) == [], f"pairs crossing {pairs_crossing(x10)}")
    check("gamma 0.10 (0,2): max CPsi, exact", 0.245, x10[(0, 2)]["max"])
    x20 = passages(exact_series(liouvillian(0.20), rho_a, 5.0), 0.001)
    check("gamma 0.20 (0,2): max CPsi, exact", 0.157, x20[(0, 2)]["max"])
    # The |0> qubits 0 and 2 are fixed by the diagonal H and the Z jumps, so each ZZ bond acts as a
    # local field on qubit 1 or 3 and the state stays a product: every pair state equals the product
    # of its one-site states, and C_corr = 0. The float propagation leaves rounding in that equality,
    # a few eps per step; the gate asks for no more than one eps (2.2e-16) per step, 500 steps.
    Li = liouvillian(0.05, H=ising_chain())
    U = expm(Li * 0.01)
    v = rho_a.ravel().copy()
    residual = 0.0
    for _ in range(501):
        rho = v.reshape(2**N, 2**N)
        for i, j in PAIRS:
            prod = np.kron(reduced(rho, [i]), reduced(rho, [j]))
            residual = max(residual, float(np.abs(reduced(rho, [i, j]) - prod).max()))
        v = U @ v
    claim("Ising (open chain): every pair state a product", residual <= 500 * 2.2e-16,
          f"max |rho_AB - rho_A x rho_B| = {residual:.1e}")

    matched = sum(ok for _, ok in VALUES)
    held = sum(ok for _, ok in CLAIMS)
    print(f"\n{matched} of {len(VALUES)} printed values match at the page's precision;"
          f" {held} of {len(CLAIMS)} claims hold.")
    for label, ok in VALUES + CLAIMS:
        if not ok:
            print(f"  not reproduced: {label}")


if __name__ == "__main__":
    main()
