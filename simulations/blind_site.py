"""Blind sites: which single-excitation states a Z-dephasing channel cannot touch.

"Blind" and not "dark": in this repo "dark" already names a mode of zero light
content <n_XY> (the AT entry of docs/ANALYTICAL_FORMULAS.md) and an F135/F136 record
class, and neither is this.

WHICH OBJECT.  Two, and the arc `site_resolved_vacuum_block` says to name them
apart every time (OpenArcsRegistry.cs:4846-4849).  The PROPAGATION runs on rho in
the (1,1) joint-popcount block, the N x N Haken-Strobl density block.  The
COUNTING runs on h_SE, the N-dimensional single-excitation Hamiltonian, and the
blind dimension reported here is a dimension of that Hilbert space, i.e. a count
of STATES, not of operators.  The blind operator space inside the (1,1) block is
larger and is not computed here.  A THIRD object appears in run_coherence: the
(0,1) vacuum-coherence block of F152, where the same count holds.  None of the
three is the Pauli w=1 sector, which is not L-invariant.

The Hamiltonian is HEISENBERG, sum_bonds J (XX + YY + ZZ), carrying the ZZ degree
diagonal: h_SE = (N-1) J Id - 2 J L with L the path Laplacian, so the modes are
the Neumann half-integer cosines at modulus N (PROOF_UNIFORM_LAW.md:182-186,
D10:134/180, F2:111).  It is NOT the XY chain of F2b at modulus N+1, which is
what Cone.cs and cone_defect_arrival.py build.

Gamma book, and it depends on the SUPPORT.  Two single-excitation configurations
differ in exactly two bits, so a coherence between them decays at
2(gamma_a + gamma_b).  With gamma on EVERY site that is 4*gamma, the Gamma = 4*gamma
of F126.  With gamma on ONE site f, which is this page's case, the largest rate
present is 2*gamma_f and never 4*gamma_f, so F126's book does not describe these
runs even though the engine convention is the same.

Conventions inherited verbatim from simulations/bridge_sector.py, which
reproduces the C# engine (RK4, dephasing mask, MSB-first qubit order).

Runs:  prep | support | dimension | identity | branch | coherence | scope | all
Usage: python simulations/blind_site.py [run]
"""

from __future__ import annotations

import itertools
import sys
from math import gcd

import numpy as np

sys.path.insert(0, __file__.replace("\\", "/").rsplit("/", 1)[0])
from bridge_sector import SectorPropagator, bit, build_h, mediator_bridge

PRIMES = ((1 << 31) - 1, (1 << 61) - 1)     # two reductions; both can understate rational rank


# ----------------------------------------------------------------------
# sectors
# ----------------------------------------------------------------------

def se_basis(n: int):
    """Single-excitation bitstrings, qubit 0 = MSB (the C# convention)."""
    states = sorted(1 << (n - 1 - a) for a in range(n))
    return states, {s: i for i, s in enumerate(states)}


def site_rows(n: int, index):
    return [index[1 << (n - 1 - a)] for a in range(n)]


def sector(n: int, popcount: int):
    states = sorted(sum(1 << (n - 1 - a) for a in c)
                    for c in itertools.combinations(range(n), popcount))
    return states, {s: i for i, s in enumerate(states)}


def chain(n: int, j: float = 1.0):
    return [(i, i + 1, j) for i in range(n - 1)]


def ket(n: int, site: int, states, index):
    v = np.zeros(len(states), dtype=complex)
    v[index[1 << (n - 1 - site)]] = 1.0
    return v


# ----------------------------------------------------------------------
# run 1: the preparations
# ----------------------------------------------------------------------

def preparations(n: int, states, index):
    """One-sided, the two coherent two-sided states, and the incoherent control.

    MIX carries the same two wavepackets, the same energy and the same Sum-gamma
    as SYM and ANTI and has no phase relation, so it cannot interfere: whatever
    SYM and ANTI do that MIX does not is interference and nothing else.  By
    linearity MIX = (SYM + ANTI)/2 exactly, which the run checks rather than
    assumes.
    """
    a, b = ket(n, 0, states, index), ket(n, n - 1, states, index)
    sym, anti = (a + b) / np.sqrt(2.0), (a - b) / np.sqrt(2.0)
    return {
        "ONE":  np.outer(a, a.conj()),
        "SYM":  np.outer(sym, sym.conj()),
        "ANTI": np.outer(anti, anti.conj()),
        "MIX":  0.5 * (np.outer(a, a.conj()) + np.outer(b, b.conj())),
    }


def run_prep(n=11, gammas=(0.0, 0.05, 0.2), t_max=20.0, dt=0.01, n_samples=201):
    states, index = se_basis(n)
    rows = site_rows(n, index)
    f = (n - 1) // 2
    for gamma in gammas:
        prop = SectorPropagator(n, chain(n), [gamma] * n, states, index)
        times = list(np.linspace(0.0, t_max, n_samples))
        occ = {}
        for name, rho0 in preparations(n, states, index).items():
            acc = []
            prop.propagate(rho0, t_max, dt, times,
                           lambda t, rho, acc=acc: acc.append(
                               [rho[r, r].real for r in rows]))
            occ[name] = np.array(acc)
        t = np.array(times[:len(occ["ONE"])])

        resid = np.abs(occ["MIX"] - 0.5 * (occ["SYM"] + occ["ANTI"])).max()
        print()
        print(f"=== gamma = {gamma}  N = {n}  J = 1  uniform chain  t_max = {t_max} ===")
        print(f"  linearity  max|MIX - (SYM+ANTI)/2| = {resid:.3e}")
        print()
        print(f"  occupancy of the reversal-fixed site {f}, every 10th sample")
        print("      t      ONE       MIX       SYM       ANTI    SYM/MIX")
        for k in range(0, len(t), 10):
            denom = occ["MIX"][k, f]
            r = occ["SYM"][k, f] / denom if denom > 1e-12 else float("nan")
            print(f"  {t[k]:6.2f}  {occ['ONE'][k, f]:.6f}  {occ['MIX'][k, f]:.6f}"
                  f"  {occ['SYM'][k, f]:.6f}  {occ['ANTI'][k, f]:.6f}  {r:9.6f}")
        sw = occ["MIX"][:, f]
        print(f"  MIX at site {f} ranges over [{sw.min():.6f}, {sw.max():.6f}]:"
              f" the ratio column above holds across that whole range, not at a")
        print("  plateau.  A max/min RATIO is not quoted because the occupancy")
        print("  passes through zero, which makes the ratio arbitrarily large.")
        print()
        print(f"  time-averaged occupancy profile over t in [0, {t_max:.0f}]")
        print("   site " + " ".join(f"{s:8d}" for s in range(n)))
        for name in ("ONE", "MIX", "SYM", "ANTI"):
            print(f"  {name:5s} " + " ".join(f"{v:8.5f}" for v in occ[name].mean(axis=0)))
        d = occ["SYM"].mean(axis=0) - occ["MIX"].mean(axis=0)
        print("  S-M   " + " ".join(f"{v:+8.5f}" for v in d))
        print(f"  the ten off-centre deficits sum to {-(d.sum() - d[f]):.5f} against the"
              f" centre's {d[f]:+.5f}: trace conservation, not an independent fact")


# ----------------------------------------------------------------------
# run 2: which dephasing support the blind state survives
# ----------------------------------------------------------------------

def blind_trajectory(n, bonds, gammas, t_max=20.0, dt=0.01):
    """Propagate the reversal-odd end state; return (max occupancy on the
    reversal-fixed site, min purity, dt actually used)."""
    states, index = se_basis(n)
    rows = site_rows(n, index)
    prop = SectorPropagator(n, bonds, gammas, states, index)
    dt = prop.stable_dt(dt)
    f = (n - 1) // 2
    psi = (ket(n, 0, states, index) - ket(n, n - 1, states, index)) / np.sqrt(2.0)
    occ, pur = [], []

    def cb(t, rho):
        occ.append(rho[rows[f], rows[f]].real)
        pur.append(np.trace(rho @ rho).real)

    prop.propagate_every_step(np.outer(psi, psi.conj()), t_max, dt, cb)
    return max(occ), min(pur), dt


def deviation(n, bonds, gammas, rho0, t_max=20.0, dt=0.01):
    """max over the trajectory of |rho_gamma(t) - rho_0(t)|, the SAME rho0 run
    once with the given gammas and once at gamma = 0.

    This, not the purity column, is what decides whether a dissipator is a
    no-op.  The RK4 truncation floor on purity here is about 1.8e-7, so two
    identical nine-digit purities bound the effect only at that level; this
    difference resolves it at machine level.
    """
    states, index = se_basis(n)
    prop_g = SectorPropagator(n, bonds, gammas, states, index)
    prop_0 = SectorPropagator(n, bonds, [0.0] * n, states, index)
    step = prop_g.stable_dt(dt)
    a, b = [], []
    prop_g.propagate_every_step(rho0.copy(), t_max, step, lambda t, r: a.append(r.copy()))
    prop_0.propagate_every_step(rho0.copy(), t_max, step, lambda t, r: b.append(r.copy()))
    return max(np.abs(x - y).max() for x, y in zip(a, b))


def run_support(n=11):
    f = (n - 1) // 2
    print(f"N = {n}, uniform chain, J = 1, reversal-fixed site f = {f}")
    print(f"preparation (|0> - |{n - 1}>)/sqrt(2), gamma = 0.05 on the listed support")
    print("dt is printed because the identical purities are evidence only if the")
    print("integrator step is the same in every row.")
    print()
    print(f"{'dephasing support':<36} {'max occ at f':>14} {'min purity':>12} {'dt':>8}")
    for label, support in [
        ("none (gamma = 0)",                       []),
        (f"the reversal-fixed site {{{f}}} only",   [f]),
        ("the left end {0} only",                  [0]),
        (f"the symmetric pair {{0, {n-1}}}",        [0, n - 1]),
        (f"the asymmetric pair {{1, {n-2}}}",       [1, n - 2]),
        ("every site (uniform)",                   list(range(n))),
    ]:
        g = [0.0] * n
        for s in support:
            g[s] = 0.05
        occ, pur, used = blind_trajectory(n, chain(n), g)
        print(f"{label:<36} {occ:14.3e} {pur:12.9f} {used:8.5f}")

    sym = mediator_bridge(3)
    broken = [(a, b, 1.6 if (a, b) == (1, 2) else j) for (a, b, j) in sym]
    matched = [(a, b, 1.6 if (a, b) in ((1, 2), (8, 9)) else j) for (a, b, j) in sym]
    print()
    print(f"MediatorBridge level 3, N = {n}: site {f} is the meta-mediator, and at")
    print("uniform J = 1 that topology IS the uniform chain.")
    print()
    print("Three arms.  Detuning bond (1,2) alone breaks the mirror, but it also")
    print("changes the couplings, so on its own it cannot say which of the two did")
    print("the work.  The third arm detunes (1,2) AND (8,9), a LARGER perturbation")
    print("(sum J 11.2 against 10.6) that PRESERVES the mirror.  If the mirror is")
    print("what matters, the third arm stays blind and the second does not.")
    print()
    print(f"{'placement':<44} {'gamma':>7} {'max occ at f':>14} {'min purity':>13} {'dt':>8}")
    for gam in (0.05, 0.5, 5.0, 50.0):
        for label, site, bonds in (
            (f"mediator (site {f}), mirror intact", f, sym),
            (f"one site off (site {f - 1}), mirror intact", f - 1, sym),
            (f"mediator (site {f}), mirror BROKEN", f, broken),
            (f"mediator (site {f}), mirror KEPT, bigger detune", f, matched),
        ):
            g = [0.0] * n
            g[site] = gam
            occ, pur, used = blind_trajectory(n, bonds, g)
            print(f"{label:<44} {gam:7.2f} {occ:14.3e} {pur:13.9f} {used:8.5f}")
        print()
    print("The purity column CANNOT establish a no-op: the RK4 truncation floor")
    print("here is about 1.8e-7, so identical nine-digit purities bound the effect")
    print("only at that level.  The observable that decides it is the difference")
    print("against the same state run at gamma = 0, over the whole trajectory:")
    print()
    states, index = se_basis(n)
    psi = (ket(n, 0, states, index) - ket(n, n - 1, states, index)) / np.sqrt(2.0)
    anti = np.outer(psi, psi.conj())
    print(f"  {'placement':<44} {'gamma':>7} {'max|rho_g - rho_0|':>20}")
    for gam in (0.05, 0.5, 5.0, 50.0):
        for label, site, bonds in (
            (f"mediator (site {f})", f, sym),
            (f"one site off (site {f - 1})", f - 1, sym),
        ):
            g = [0.0] * n
            g[site] = gam
            print(f"  {label:<44} {gam:7.2f} {deviation(n, bonds, g, anti):20.3e}")
    print()
    print("Note the 'one site off' purity ladder is NOT monotone in gamma")
    print("(0.574, 0.133, 0.144, 0.632): strong local dephasing begins to freeze")
    print("the state it is destroying.  Reported, not smoothed.")


# ----------------------------------------------------------------------
# run 3: modular blind counts against the exact uniform-chain law
# ----------------------------------------------------------------------

def build_h_int(n, bonds, states, index):
    """Integer H on the sector; requires integer J. Python ints, no overflow."""
    m = len(states)
    h = [[0] * m for _ in range(m)]
    for col, s in enumerate(states):
        for (a, b, j) in bonds:
            ja = int(j)
            assert ja == j, "integer J required for the modular rank"
            ba, bb = bit(s, a, n), bit(s, b, n)
            h[col][col] += ja * (1 - 2 * ba) * (1 - 2 * bb)
            if ba != bb:
                flipped = s ^ (1 << (n - 1 - a)) ^ (1 << (n - 1 - b))
                h[index[flipped]][col] += 2 * ja
    return h


def rref_modp(rows, m, p):
    rows = [r[:] for r in rows]
    rank, col = 0, 0
    while col < m and rank < len(rows):
        piv = next((r for r in range(rank, len(rows)) if rows[r][col] % p), None)
        if piv is None:
            col += 1
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][col], p - 2, p)
        rows[rank] = [(x * inv) % p for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col] % p:
                fq = rows[r][col]
                rows[r] = [(x - fq * y) % p for x, y in zip(rows[r], rows[rank])]
        rank += 1
        col += 1
    return rows[:rank]


def blind_dim(n, bonds, popcount, sites, p=PRIMES[0]):
    """dim of the largest H-invariant subspace inside the intersection of
    ker(n_k) over k in `sites`.

    Z_k = I - 2 n_k has eigenvalue +1 on the configurations avoiding k and -1 on
    those occupying k, and the dissipator is silent on any rho with no coherence
    ACROSS those eigenspaces.  A PURE state, hence a subspace, is silent iff it
    lies wholly in one eigenspace.  In the single-excitation sector the -1
    eigenspace is the single ray |k>, which is not H-invariant whenever site k
    carries an incident bond with J != 0 (true for every topology used here, but
    NOT automatic: an isolated site would make that branch one-dimensional).  So
    only the +1 branch contributes, and the TRAJECTORY is untouched iff the state
    is orthogonal to the block Krylov space generated by the occupied-at-k
    configurations.  (H Hermitian is needed and
    holds: <c|H^m|v> = <H^m c|v>, so the complement is the intersection of
    ker(n_k H^m).)

    Taken as a rank over GF(p) on the integer sector matrix at J = 1.  A
    float rank miscounts here: the first version of this routine used an SVD and
    returned blind dimensions that were not mirror-symmetric on a mirror-
    symmetric chain, which is impossible.  The iteration runs to a fixed point of
    the RANK, not of the row count; testing the row count exits early and
    understates the Krylov space.

    One prime gives rank_p <= rank_Q, so the returned dimension is an UPPER bound
    on rational blindness.  `run_dimension` repeats values at a second prime
    to catch some bad reductions; agreement alone is not an exact certificate.
    """
    states, index = sector(n, popcount)
    h = build_h_int(n, bonds, states, index)
    m = len(states)
    frontier = [[1 if i == s else 0 for i in range(m)]
                for s, st in enumerate(states) if any(bit(st, k, n) for k in sites)]
    if not frontier:
        return m
    basis = []
    for _ in range(m + 2):
        basis = rref_modp(basis + frontier, m, p)
        if len(basis) >= m:
            break
        frontier = [[sum(h[i][j] * v[j] for j in range(m)) % p for i in range(m)]
                    for v in basis]
        if len(rref_modp(basis + frontier, m, p)) == len(basis):
            break
    return m - len(basis)


def mode_set(n, j):
    """Indices of the Neumann modes with an exact node at site j.

    u_m(j) proportional to cos(pi m (2j+1)/(2N)) vanishes iff (2j+1)m = N mod 2N.
    """
    return {m for m in range(n) if ((2 * j + 1) * m - n) % (2 * n) == 0}


def run_dimension(n_max=21):
    print("Modular blind-dimension reads per site, uniform chain, single excitation.")
    print("Two-prime agreement checks reductions; the analytic node law supplies")
    print("the exact uniform-chain value, since either prime can overstate it.")
    print()
    print("closed form:  dim = (gcd(2j+1, N) - 1) / 2")
    print()
    bad_form, bad_prime = [], []
    for n in range(3, n_max + 1):
        bonds = chain(n, 1)
        first = [blind_dim(n, bonds, 1, [j], PRIMES[0]) for j in range(n)]
        second = [blind_dim(n, bonds, 1, [j], PRIMES[1]) for j in range(n)]
        pred = [(gcd(2 * j + 1, n) - 1) // 2 for j in range(n)]
        if first != pred:
            bad_form.append(n)
        if first != second:
            bad_prime.append(n)
        mirror = all(first[j] == first[n - 1 - j] for j in range(n))
        print(f"  N = {n:>2}  {'ok ' if first == pred else 'BAD'}  "
              f"mirror {str(mirror):5s}  {first}")
    print()
    print(f"  modular rows match the closed form, N = 3..{n_max}: {not bad_form}")
    print(f"  the two primes agree everywhere:                 {not bad_prime}")

    print()
    print("The intersection law: a blind subspace is a SPAN OF MODES, so for a")
    print("support S the dimension is the number of modes with a node at every")
    print("site of S.  Checked against the rank for all 1-, 2- and 3-site supports.")
    bad_int = []
    for n in (9, 11, 12, 15):
        for r in (1, 2, 3):
            for s in itertools.combinations(range(n), r):
                pred = len(set.intersection(*[mode_set(n, j) for j in s]))
                if pred != blind_dim(n, chain(n, 1), 1, list(s)):
                    bad_int.append((n, s))
    print(f"  intersection law holds at N = 9, 11, 12, 15: {not bad_int}")
    print("  worked cases at N = 15:")
    for s in ((7,), (2,), (12,), (2, 12), (7, 2), (1, 4, 13)):
        print(f"    S = {str(list(s)):<12} modes {sorted(set.intersection(*[mode_set(15, j) for j in s]))}"
              f"  dim {blind_dim(15, chain(15, 1), 1, list(s))}")

    print()
    print("The same rank in the popcount-2 sector (the Bell-on-vacuum sector), N = 11.")
    d2 = [blind_dim(11, chain(11, 1), 2, [k]) for k in range(11)]
    print(f"  blind dimension per site: {d2}")
    print("  (a full-rank Krylov space certifies these zeros at one prime already,")
    print("   since rank_p <= rank_Q and rank_p is already maximal.)")


# ----------------------------------------------------------------------
# run 4: the two identities the prose leans on
# ----------------------------------------------------------------------

def run_identity():
    print("Identity 1: h_SE = (N-1) J Id - 2 J L, the path Laplacian")
    print("(PROOF_UNIFORM_LAW.md:182-186 states this and gates it entry-exactly;")
    print("recomputed here only so this page's numbers rest on it visibly)")
    print()
    print(f"  {'N':>3} {'max off-diagonal':>18} {'diagonal spread':>17} {'constant':>10}")
    for n in (5, 8, 11, 12, 15):
        states, index = se_basis(n)
        rows = site_rows(n, index)
        h = build_h(n, chain(n), states, index).real[np.ix_(rows, rows)]
        a = np.zeros((n, n))
        for i in range(n - 1):
            a[i, i + 1] = a[i + 1, i] = 1.0
        lap = np.diag(a.sum(1)) - a
        c = h + 2.0 * lap
        off = np.abs(c - np.diag(np.diag(c))).max()
        print(f"  {n:>3} {off:18.2e} "
              f"{np.diag(c).max() - np.diag(c).min():17.2e} {np.diag(c)[0]:10.1f}")

    print()
    print("Identity 2: the spectrum against Lemma 5, lambda_k = 4cos(k pi/N) + N - 5")
    print("(PROOF_R90_FROZEN_DIVISOR.md:178; an eigensolver is used HERE and only")
    print("here, to compare two published closed forms, never to count anything)")
    print()
    for n in (5, 8, 11, 12, 15):
        states, index = se_basis(n)
        rows = site_rows(n, index)
        h = build_h(n, chain(n), states, index).real[np.ix_(rows, rows)]
        dev = np.abs(np.sort(np.linalg.eigvalsh(h))
                     - np.sort([4 * np.cos(k * np.pi / n) + n - 5
                                for k in range(n)])).max()
        print(f"  N = {n:>3}  max deviation {dev:.3e}")

    print()
    print("Identity 3: at odd N the centre's blind subspace IS the reversal-odd")
    print("site space, as a SET, not merely the same integer.  Exact integer")
    print("arithmetic: every vector e_j - e_(N-1-j) must have (H^m v)[f] == 0.")
    print()
    for n in (5, 7, 9, 11, 13, 15, 17, 19, 21):
        f = (n - 1) // 2
        states, index = sector(n, 1)
        hs = build_h_int(n, chain(n, 1), states, index)
        row = [index[1 << (n - 1 - a)] for a in range(n)]
        h = [[hs[row[a]][row[b]] for b in range(n)] for a in range(n)]
        ok = True
        for j in range(n // 2):
            v = [0] * n
            v[j], v[n - 1 - j] = 1, -1
            for _ in range(n):
                if v[f] != 0:
                    ok = False
                    break
                v = [sum(h[i][k] * v[k] for k in range(n)) for i in range(n)]
            if not ok:
                break
        d = blind_dim(n, chain(n, 1), 1, [f])
        print(f"  N = {n:>3}  every reversal-odd vector blind: {str(ok):5s}"
              f"   dim reversal-odd {n // 2:>2}  blind dim {d:>2}"
              f"   same set: {ok and d == n // 2}")


# ----------------------------------------------------------------------
# run 5: the branch this page's first draft got wrong
# ----------------------------------------------------------------------

def run_branch(n=11):
    """Z_k rho Z_k = rho does NOT require lying in ker(n_k).

    Z_k has eigenvalues +1 (on the configurations avoiding k) and -1 (on those
    occupying k), and the dissipator is silent on any rho with no coherence
    ACROSS those two eigenspaces.  ker(n_k) is the +1 branch only.

    For a SUBSPACE OF PURE STATES the count is unaffected, because a pure state
    blind at k must lie wholly inside one eigenspace, and in the
    single-excitation sector the -1 eigenspace is the single ray |k>, which is
    not H-invariant (H hops the excitation off k), so the largest H-invariant
    subspace inside it is 0.  For MIXED states the blind set is strictly larger,
    and the identity is the standard example.  This run shows both, because the
    first draft of the page stated the wrong iff and then reconciled the
    all-sites theorem with it.
    """
    states, index = se_basis(n)
    rows = site_rows(n, index)
    f = (n - 1) // 2

    print("A. the identity is blind under the FULL support, non-trivially")
    ident = np.eye(n, dtype=complex) / n
    print(f"   uniform gamma = 0.5 on all {n} sites:"
          f"  max|rho(t) - I/N| = {deviation(n, chain(n), [0.5] * n, ident):.3e}")
    print("   So the blind set at full support is NOT empty: it is the sector")
    print("   steady state, which is exactly the one state per sector that")
    print("   SYMMETRY_CENSUS.md:101 and PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md")
    print("   already report.  What is empty at full support is the blind space")
    print("   of NON-STATIONARY trajectories.")

    print()
    print("B. blindness is about SUPPORT, not about purity")
    a = (ket(n, 0, states, index) - ket(n, n - 1, states, index)) / np.sqrt(2.0)
    b = (ket(n, 1, states, index) - ket(n, n - 2, states, index)) / np.sqrt(2.0)
    mixed = 0.6 * np.outer(a, a.conj()) + 0.4 * np.outer(b, b.conj())
    print(f"   a mixture of two reversal-odd states, purity"
          f" {np.trace(mixed @ mixed).real:.4f}")
    for gam in (0.05, 50.0):
        g = [0.0] * n
        g[f] = gam
        print(f"   gamma_{f} = {gam:>6}:  max|rho_g - rho_0| ="
              f" {deviation(n, chain(n), g, mixed):.3e}")
    print("   A state at purity 0.52 is exactly as blind as the pure one, so")
    print("   'purity protects' is the wrong label; the support does.")

    print()
    print("C. sparse-support control for PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md")
    print("   The theorem requires gamma_k > 0 at EVERY site.  This run lies")
    print("   outside that scope: the blind initial state does not converge")
    print("   to P_w / d_w under centre-only dephasing:")
    anti = np.outer(a, a.conj())
    g = [0.0] * n
    g[f] = 0.5
    prop = SectorPropagator(n, chain(n), g, states, index)
    out = []
    prop.propagate_every_step(anti.copy(), 200.0, prop.stable_dt(0.01),
                              lambda t, r: out.append(r))
    end = out[-1]
    print(f"   gamma = 0.5 on site {f} only, t = 200:"
          f"  |rho - I/N| = {np.abs(end - ident).max():.3e},"
          f"  purity = {np.trace(end @ end).real:.6f}")
    print("   The trajectory stays pure and never reaches the maximally mixed")
    print("   sector state.  The proof's hypothesis, gamma_k > 0 at every site,")
    print("   is what its Step 2(a) spends; this support lies outside it.")

    print()
    print("D. the Bell-on-vacuum family, dynamically rather than by rank")
    print("   The popcount-2 rank says 0 at every site, but that rank covers the")
    print("   (2,2) block only, while this state also occupies (0,0) and (0,2).")
    print("   The direct check, with the state propagated whole:")
    from bridge_sector import sector_basis, bell_on_vacuum
    st2, ix2 = sector_basis(n)
    rho2 = bell_on_vacuum(n, 0, 1, st2, ix2)
    for site in (f, 0):
        g = [0.0] * n
        g[site] = 0.5
        pg = SectorPropagator(n, chain(n), g, st2, ix2)
        p0 = SectorPropagator(n, chain(n), [0.0] * n, st2, ix2)
        step = pg.stable_dt(0.01)
        aa, bb = [], []
        pg.propagate_every_step(rho2.copy(), 20.0, step, lambda t, r: aa.append(r.copy()))
        p0.propagate_every_step(rho2.copy(), 20.0, step, lambda t, r: bb.append(r.copy()))
        d = max(np.abs(x - y).max() for x, y in zip(aa, bb))
        print(f"   gamma = 0.5 on site {site:>2}:  max|rho_g - rho_0| = {d:.3e}")
    print("   Nonzero at both, so the Bell-on-vacuum family has no blind state")
    print("   anywhere, which is why the mediator looked ordinary to it.")


# ----------------------------------------------------------------------
# run 6: the same count inside F64's own block
# ----------------------------------------------------------------------

def run_coherence():
    """The count transfers to the (0,1) vacuum-coherence block, exactly.

    EQ-015's identity is -Re(lambda_k) = 2 gamma_B |v_k(B)|^2 exactly, mode by
    mode, at every gamma_B.  Its generator L_coh = i H_1 - 2 gamma_B |B><B|
    (EMERGING_QUESTIONS.md:508) is derived for an XY Hamiltonian.  For HEISENBERG
    the vacuum bra carries E_vac = (N-1) J, so the block generator is
    M = -i (h_SE - E_vac Id) - 2 gamma_B |B><B|, which is exactly F152's
    M = -2i L_J - 2 diag(gamma) (ANALYTICAL_FORMULAS.md:7058), since
    h_SE = (N-1) J Id - 2 J L.  This run uses F152's form.  The XY form differs
    by a sign on the Hamiltonian part and an imaginary shift, neither of which
    touches Re(lambda), so the count is the same either way; the correct
    generator is used anyway because the whole point of this page is that the ZZ
    diagonal is what makes the chain Heisenberg.

    The bijection is two lines and needs no eigensolver.  If v is an H_1
    eigenvector with v(B) = 0 then L_coh v = i H_1 v, so v is an L_coh
    eigenvector with Re(lambda) = 0.  Conversely Re(lambda) = 0 forces
    |v(B)|^2 = 0 by the identity, and then L_coh v = i H_1 v = lambda v makes v
    an H_1 eigenvector with a node at B.  So the undamped modes of the
    coherence block are exactly the node modes counted in the density block,
    and the closed form applies there unchanged.

    The run below is a CHECK of that bijection, not its proof; it uses a
    general eigensolver on a non-normal matrix, which is why the argument above
    is what carries the claim.
    """
    print("The same count in the (0,1) coherence block, F64's own object.")
    print("blind there = Re(lambda) = 0 for F152 M = -2i L_J - 2 gamma_B |B><B|")
    print()
    print("The 1e-9 cut is reported with the gap it sits in: the largest |Re|")
    print("counted as zero beside the smallest |Re| excluded.  A threshold that")
    print("merely passes is worthless; these are decades apart.")
    print()
    print(f"  {'N':>3} {'site':>5} {'closed form':>12} {'Re=0 at 0.01':>13}"
          f" {'Re=0 at 0.7':>12} {'largest counted':>16} {'smallest excluded':>18}")
    bad = 0
    for n in (9, 11, 12, 15):
        states, index = se_basis(n)
        rows = site_rows(n, index)
        h1 = build_h(n, chain(n), states, index).real[np.ix_(rows, rows)]
        for j in range(n):
            pred = (gcd(2 * j + 1, n) - 1) // 2
            proj = np.zeros((n, n))
            proj[j, j] = 1.0
            counts, worst, nearest = [], 0.0, float("inf")
            for gam in (0.01, 0.7):
                lap = (np.eye(n) * (n - 1) - h1) / 2.0        # h_SE = (N-1)Id - 2L
                m = -2j * lap - 2 * gam * proj                # F152
                re = np.abs(np.linalg.eigvals(m).real)
                counts.append(int((re < 1e-9).sum()))
                if counts[-1]:
                    worst = max(worst, re[re < 1e-9].max())
                if (re >= 1e-9).any():
                    nearest = min(nearest, re[re >= 1e-9].min())
            if counts[0] != pred or counts[1] != pred:
                bad += 1
            if pred > 0:
                print(f"  {n:>3} {j:>5} {pred:>12} {counts[0]:>13} {counts[1]:>12}"
                      f" {worst:>16.2e} {nearest:>18.2e}")
    print()
    print(f"  mismatches over EVERY site at N = 9, 11, 12, 15, both rates: {bad}")


# ----------------------------------------------------------------------
# run 7: what a blind seat does to the asymptotic sector projection
# ----------------------------------------------------------------------

def _row_krylov_rank_modp(ad, seeds, p=PRIMES[0]):
    """Rank over GF(p) of the smallest row space that contains the unit rows
    `seeds` and is closed under right multiplication by the integer matrix
    `ad`.  Its complement is the largest ad-invariant subspace on which all
    those rows vanish.  A reduction mod p can only LOWER a rank, so m minus
    this rank is an UPPER bound on the rational dimension of that subspace.
    Overflow: entries stay below p = 2^31 - 1, so a product of two fits int64,
    and ad has small integer entries, so a row times ad fits as well."""
    m = ad.shape[0]
    basis = {}                                  # pivot -> fully reduced row
    queue = []
    for s in seeds:
        v = np.zeros(m, dtype=np.int64)
        v[s] = 1
        queue.append(v)
    while queue:
        v = queue.pop() % p
        for c, b in basis.items():
            if v[c]:
                v = (v - v[c] * b) % p
        nz = np.flatnonzero(v)
        if nz.size == 0:
            continue
        c = int(nz[0])
        v = (v * pow(int(v[c]), p - 2, p)) % p
        for c2 in basis:
            if basis[c2][c]:
                basis[c2] = (basis[c2] - basis[c2][c] * v) % p
        basis[c] = v
        queue.append((v @ ad) % p)
    return len(basis)


def _undamped_bound(blocks, support, n):
    """Upper bound mod p on the dimension of the UNDAMPED space of L, the span
    of every eigenvector with Re(lambda) = 0.  The jumps Z_k are Hermitian, so
    D is self-adjoint and non-positive; Re(lambda) = 0 then forces D v = 0, and
    the undamped space is the largest ad_H-invariant subspace inside ker D.  In
    a basis of cells |a><b| in which D is diagonal, ker D is spanned by the cells
    that agree on the support, so the bound is (cells) minus the row-Krylov rank
    of the disagreeing cells.  `blocks` yields (reps_a, h_a, reps_b, h_b) with
    h the integer matrix of H on the chosen basis of each side."""
    total = 0
    for ra, ha, rb, hb in blocks:
        ad = (np.kron(ha, np.eye(len(rb), dtype=np.int64))
              - np.kron(np.eye(len(ra), dtype=np.int64), hb.T))
        seeds = [ia * len(rb) + ib for ia, a in enumerate(ra) for ib, b in enumerate(rb)
                 if any(bit(a ^ b, k, n) for k in support)]
        total += len(ra) * len(rb) - (_row_krylov_rank_modp(ad, seeds) if seeds else 0)
    return total


def _popcount_blocks(n, bonds):
    """Every joint-popcount block (p, q) in the computational basis."""
    secs = []
    for w in range(n + 1):
        st, ix = sector(n, w)
        secs.append((st, np.array(build_h_int(n, bonds, st, ix), dtype=np.int64)))
    for sa, ha in secs:
        for sb, hb in secs:
            yield sa, ha, sb, hb


def _parity_sectors(n, bonds, mirror_of, sign):
    """Per popcount w, the mirror-parity part of sign +1 or -1, in the basis
    e_i = |i> + sign |M i> (for sign +1 also |i> when M fixes i), one e_i per
    mirror orbit.  H commutes with M, so H e_i = sum_k alpha[k, i] e_k with
    integer alpha, and a Z_k at a mirror-fixed seat stays diagonal here."""
    secs = []
    for w in range(n + 1):
        st, ix = sector(n, w)
        h = build_h_int(n, bonds, st, ix)
        reps = [s for s in st if s < mirror_of[s] or (sign > 0 and s == mirror_of[s])]
        al = np.zeros((len(reps), len(reps)), dtype=np.int64)
        for ci, i in enumerate(reps):
            for ck, k in enumerate(reps):
                val = h[ix[k]][ix[i]]
                if mirror_of[i] != i:
                    val += sign * h[ix[k]][ix[mirror_of[i]]]
                al[ck, ci] = val
        secs.append((reps, al))
    return secs


def _even_blocks(n, bonds, mirror_of):
    """Every block restricted to the mirror-EVEN space on both sides."""
    secs = _parity_sectors(n, bonds, mirror_of, +1)
    for ra, ha in secs:
        for rb, hb in secs:
            yield ra, ha, rb, hb.T          # _undamped_bound transposes the bra side back


def run_scope(n=5, times=(200.0, 400.0), dt=0.02):
    """PROOF_ASYMPTOTIC_SECTOR_PROJECTION predicts rho(inf) = sum_w p_w P_w/d_w.

    Its argument spends gamma_k > 0 at EVERY site (Step 2 part (a), and the
    separate Step 2b for the off-diagonal blocks) and a CONNECTED graph (Step 2
    part (b)).  That hypothesis is sufficient but not necessary, and this run
    decides, support by support, which of its supports keep the conclusion.

    Every support is run at two times, and the two residuals are only a finite
    read.  What decides the limit is the certificate after the table: the
    undamped space of L, bounded from above by a GF(p) rank and from below by
    the N+1 sector projectors P_w, which are stationary for every support.
    Where the two bounds meet, the undamped space IS the span of the P_w, and
    every initial state converges to the per-sector mixed state.  Where they do
    not (the centre seat), a sharper lower bound is assembled from three
    families in different parity blocks, the even side, the odd parts and the
    parity-crossing coherences of the vacuum and the all-excited state, the
    last checked over Q, and compared with the same upper bound.

    Full Hilbert space, not a sector: the predicted limit mixes every popcount,
    so the state is propagated whole.  |+>^N is used because it puts weight in
    all N+1 sectors and carries every inter-sector coherence.
    """
    states = list(range(1 << n))
    index = {s_: i for i, s_ in enumerate(states)}
    bonds = chain(n)
    psi = np.ones(1 << n, dtype=complex)
    psi /= np.linalg.norm(psi)
    rho0 = np.outer(psi, psi.conj())

    d = 1 << n
    sector_idx = {w: [i for i in range(d) if bin(i).count("1") == w]
                  for w in range(n + 1)}
    pred = np.zeros((d, d), dtype=complex)
    for w, idx in sector_idx.items():
        pw = sum(rho0[i, i].real for i in idx)
        for i in idx:
            pred[i, i] += pw / len(idx)

    def run_to(g, t_max):
        prop = SectorPropagator(n, bonds, g, states, index)
        last = [None]
        prop.propagate_every_step(rho0.copy(), t_max, prop.stable_dt(dt),
                                  lambda t, r, l=last: l.__setitem__(0, r))
        return last[0]

    print(f"N = {n} chain, |+>^N: does rho reach sum_w p_w P_w/d_w?")
    print(f"blind single-excitation dimensions per seat, (gcd(2j+1, N) - 1)/2: "
          f"{[(gcd(2 * j + 1, n) - 1) // 2 for j in range(n)]}")
    print()
    print("The blind column is the INTERSECTION over the support, not a sum: a")
    print("seat with a blind subspace loses it as soon as a second seat is")
    print("dephased whose modes do not also vanish there.")
    print()
    t0, t1 = times
    supports = ([0, 1, 2, 3, 4], [0, 1, 2, 3], [1, 2, 3], [0, 2, 4],
                [0, 4], [1, 3], [0], [2])
    print(f"  {'dephasing support':<24} {'blind dim':>9} {'dev at t=' + str(int(t0)):>16}"
          f" {'dev at t=' + str(int(t1)):>16} {'ratio':>9} {'finite read':>11}")
    late = {}
    for sup in supports:
        g = [0.0] * n
        for site in sup:
            g[site] = 0.5
        a = np.abs(run_to(g, t0) - pred).max()
        rho1 = run_to(g, t1)
        late[tuple(sup)] = rho1
        b = np.abs(rho1 - pred).max()
        ratio = a / b if b > 0 else float("inf")
        finite_read = "at floor" if b < 1e-14 else ("falling" if ratio > 10.0 else "plateau")
        print(f"  {str(sup):<24} {len(set.intersection(*[mode_set(n, j) for j in sup])):>9}"
              f" {a:>16.3e} {b:>16.3e} {ratio:>9.1e} {finite_read:>11}")
    print()
    print("  These labels describe only the two sampled times: a residual falls")
    print("  by >10, sits below the 1e-14 grid floor, or stays on a plateau.")
    print("  The all-sites row's 1e-176 is a tiny exponential tail, reported")
    print("  as printed and not as a precision claim.  What decides the limit is")
    print("  the certificate below, not these two snapshots.")

    print()
    print("  Certificate.  The jumps Z_k are Hermitian, so D is self-adjoint and")
    print("  non-positive, and an eigenvalue on the imaginary axis forces D v = 0:")
    print("  the UNDAMPED space of L is the largest ad_H-invariant subspace inside")
    print("  ker D.  A GF(p) rank bounds its dimension from ABOVE; the N+1 sector")
    print("  projectors P_w, stationary for every support, bound it from BELOW.")
    print("  Where the two meet, the undamped space is the span of the P_w and")
    print("  EVERY initial state converges to sum_w p_w P_w/d_w.")
    print()
    print(f"  {'dephasing support':<24} {'undamped <= (mod p)':>20} {'>= (the P_w)':>13}   verdict")
    beyond = []
    for sup in supports:
        ub = _undamped_bound(_popcount_blocks(n, bonds), sup, n)
        if ub == n + 1:
            verdict = "every state converges"
        else:
            verdict = "the bounds do not meet"
            beyond.append(sup)
        print(f"  {str(sup):<24} {ub:>20} {n + 1:>13}   {verdict}")

    if not beyond:
        return
    mir = [int(format(bb, "0%db" % n)[::-1], 2) for bb in range(d)]
    print()
    print("  Where the bounds do not meet, one commutator names an extra invariant.")
    print()
    perm = np.zeros((d, d))
    for bb in range(d):
        perm[mir[bb], bb] = 1.0
    for k in range(n):
        zk = np.diag([1.0 - 2 * ((bb >> (n - 1 - k)) & 1) for bb in range(d)])
        c = np.abs(zk @ perm - perm @ zk).max()
        print(f"    seat {k}: max|[Z_k, M]| = {c:.1f}"
              f"   {'preserves mirror parity' if c < 1e-12 else 'mixes it'}")
    print()
    print("  M is the site reversal.  Z_k commutes with it only at the mirror-fixed")
    print("  seat, so only a support inside the fixed set conserves the two parity")
    print("  projectors, constants of motion beyond the N+1 sector populations.")
    print("  What they do to the limit depends on where the state sits.  On the")
    print("  mirror-EVEN side the same certificate decides it, run on operators")
    print("  that live on the even space on both sides (the basis |i> + |M i>,")
    print("  in which a jump at a fixed seat stays diagonal):")
    print()
    for sup in beyond:
        assert all(k == n - 1 - k for k in sup), "even-side certificate needs a fixed-seat support"
        n_even = sum(1 for w in range(n + 1)
                     if any(s <= mir[s] for s in sector(n, w)[0]))
        ub_even = _undamped_bound(_even_blocks(n, bonds, mir), sup, n)
        meet = "the bounds meet" if ub_even == n_even else "the bounds do NOT meet"
        print(f"    support {sup}: undamped on the even side <= {ub_even} (mod p),"
              f" >= {n_even} (the even part of each P_w): {meet}")
    print()
    print("  So a state with no weight in the mirror-odd space converges to the")
    print("  maximally mixed state on the even part of each sector.  The odd side")
    print("  need not mix at all.  At popcount 1 the mirror-odd space has the basis")
    print("  e_j - e_(N-1-j), j < (N-1)/2, with no amplitude at the centre, so the")
    print("  centre jump acts on it as the identity and only H moves it.  This is")
    print("  the centre's blind subspace (run_identity, identity 3).  H there, as")
    print("  exact integers:")
    st1, ix1 = sector(n, 1)
    h1 = build_h_int(n, bonds, st1, ix1)
    half = (n - 1) // 2

    def at(site):
        return ix1[1 << (n - 1 - site)]

    odd_basis = []
    for j in range(half):
        v = [0] * n
        v[at(j)], v[at(n - 1 - j)] = 1, -1
        odd_basis.append(v)
    images = [[sum(h1[r][c] * v[c] for c in range(n)) for r in range(n)] for v in odd_basis]
    closed = all(img[at(half)] == 0 and all(img[at(k)] == -img[at(n - 1 - k)] for k in range(half))
                 for img in images)
    h_odd = [[images[jj][at(k)] for jj in range(half)] for k in range(half)]
    print(f"    H_odd = {h_odd}   (image stays in the odd span: {closed})")
    print("  A pure state there that is not an eigenvector of H_odd moves")
    print("  unitarily and forever, at every gamma.")
    print()
    print("  Between popcounts 1 and N-1 the odd side does relax.  The same certificate")
    print("  on each popcount's odd part (odd on both sides), with the part's own")
    print("  identity as the lower bound, and the whole of End(odd part) at popcounts")
    print("  1 and N-1, where the centre jump acts on the odd part as +1 or -1:")
    odd_parts = _parity_sectors(n, bonds, mir, -1)
    for sup in beyond:
        for w, (reps, al) in enumerate(odd_parts):
            if not reps:
                continue
            ub = _undamped_bound([(reps, al, reps, al.T)], sup, n)
            lb = len(reps) ** 2 if w in (1, n - 1) else 1
            print(f"    support {sup}, popcount {w}: odd part dim {len(reps)},"
                  f" undamped <= {ub} (mod p), >= {lb}:"
                  f" {'the bounds meet' if ub == lb else 'the bounds do NOT meet'}")

    from fractions import Fraction
    print()
    print("  The parity-diagonal parts account for the even side and the odd parts,")
    print("  and they leave the parity-crossing operators, even on one side and odd")
    print("  on the other.  There the vacuum and the all-excited state, both H-eigenstates")
    print("  and both mirror-even, pair with the odd parts at popcounts 1 and N-1:")
    print("  the odd part at popcount 1 has no amplitude at the centre and the one at")
    print("  N-1 is filled there, like the vacuum and the all-excited state respectively,")
    print("  so each |0..0><o|, |o><0..0|, |1..1><o'|, |o'><1..1| agrees on the centre and")
    print("  D kills none of them.  Built here as integer matrices and checked over")
    print("  Q, no prime: every entry agrees on the support, and [H, X] stays in their")
    print("  span, so the span is ad_H-invariant inside ker D and hence undamped.")
    h_full = np.array(build_h_int(n, bonds, list(range(d)), {s_: s_ for s_ in range(d)}),
                      dtype=np.int64)
    top = d - 1
    cross_ops = []
    for w_side, anchor in ((1, 0), (n - 1, top)):
        for reps in (s_ for s_ in sector(n, w_side)[0] if s_ < mir[s_]):
            o = np.zeros(d, dtype=np.int64)
            o[reps], o[mir[reps]] = 1, -1
            e = np.zeros(d, dtype=np.int64)
            e[anchor] = 1
            cross_ops += [np.outer(e, o), np.outer(o, e)]

    def rank_q(mats):
        rows = []
        for mat in mats:
            nz = {int(k): Fraction(int(v)) for k, v in enumerate(mat.ravel()) if v}
            for piv, prow in rows:
                if piv in nz:
                    f = nz[piv] / prow[piv]
                    for k, v in prow.items():
                        nz[k] = nz.get(k, 0) - f * v
                        if nz[k] == 0:
                            del nz[k]
            if nz:
                rows.append((min(nz), nz))
        return len(rows)

    for sup in beyond:
        in_ker_d = all(all(bit(a_, k, n) == bit(b_, k, n) for k in sup)
                       for mat in cross_ops for a_, b_ in zip(*np.nonzero(mat)))
        comms = [h_full @ mat - mat @ h_full for mat in cross_ops]
        r_ops, r_both = rank_q(cross_ops), rank_q(cross_ops + comms)
        # the crossing family counts only if both checks pass; the even side and
        # the odd parts count only where their own certificates above met
        crossing = r_ops if (in_ker_d and r_both == r_ops) else 0
        lb_even = sum(1 for w in range(n + 1) if any(s <= mir[s] for s in sector(n, w)[0]))
        ub_even = _undamped_bound(_even_blocks(n, bonds, mir), sup, n)
        lb_odd = 0
        for w, (reps, al) in enumerate(odd_parts):
            if not reps:
                continue
            lb_w = len(reps) ** 2 if w in (1, n - 1) else 1
            if _undamped_bound([(reps, al, reps, al.T)], sup, n) == lb_w:
                lb_odd += lb_w
        lb_total = (lb_even if ub_even == lb_even else 0) + lb_odd + crossing
        ub = _undamped_bound(_popcount_blocks(n, bonds), sup, n)
        print(f"    support {sup}: {len(cross_ops)} crossing operators, rank over Q {r_ops},"
              f" in ker D: {in_ker_d}, rank with their commutators {r_both}"
              f" ({'invariant' if r_both == r_ops else 'NOT invariant'})")
        print(f"    support {sup}: undamped >= {lb_total} (even side + odd parts + crossing),"
              f" <= {ub} (mod p): {'the bounds meet' if ub == lb_total else 'the bounds do NOT meet'}")
    print("  The three families sit in different parity blocks, so their dimensions")
    print("  add.  Where the bounds meet, the undamped space at the centre is exactly")
    print("  the even part of each sector projector, all of End(odd part) at popcounts")
    print("  1 and N-1, the identity of the odd part at the popcounts between, and the")
    print("  vacuum and all-excited coherences with the odd parts.  Whatever part of a")
    print("  state lies outside that space decays, and nothing else survives.")

    from math import comb
    proj = {}
    for w, idx in sector_idx.items():
        seen, vecs = set(), {+1: [], -1: []}
        for i in idx:
            j = mir[i]
            if i in seen:
                continue
            seen.add(i)
            seen.add(j)
            e = np.zeros(d)
            e[i] = 1.0
            f = np.zeros(d)
            f[j] = 1.0
            if i == j:
                vecs[+1].append(e)
            else:
                vecs[+1].append((e + f) / np.sqrt(2.0))
                vecs[-1].append((e - f) / np.sqrt(2.0))
        for sgn in (+1, -1):
            if vecs[sgn]:
                V = np.array(vecs[sgn]).T
                proj[(w, sgn)] = V @ V.T
    refined = np.zeros((d, d), dtype=complex)
    for key, pr_ in proj.items():
        refined += np.trace(pr_ @ rho0).real * pr_ / np.trace(pr_).real
    odd_initial = sum(np.trace(pr_ @ rho0).real
                      for (w, sgn), pr_ in proj.items() if sgn == -1)
    print()
    print(f"  |+>^N is mirror-even: initial mirror-odd weight {odd_initial:.3e}.")
    print("  Its limit is therefore the per-(sector, even part) mixed state, which")
    print("  differs from the per-sector one EXACTLY by, sector by sector:")
    worst = Fraction(0)
    for w in range(n + 1):
        cw = comb(n, w)
        fixed = sum(1 for s in sector(n, w)[0] if mir[s] == s)
        dp = (cw + fixed) // 2
        pw = Fraction(cw, 2 ** n)
        on_fixed = pw / dp - pw / cw
        worst = max(worst, abs(on_fixed))
        if cw == fixed:
            print(f"    w = {w}: no difference, every configuration is mirror-symmetric")
            continue
        line = f"    w = {w}: {on_fixed} on each of {fixed} mirror-symmetric configurations"
        if cw > fixed:
            on_pair = pw / (2 * dp) - pw / cw
            coh = pw / (2 * dp)
            worst = max(worst, abs(on_pair), coh)
            line += f"; {on_pair} on each of {cw - fixed} with a partner; {coh} on each partner coherence"
        print(line)
    print(f"  largest exact entry of the difference: {worst} = {float(worst):.3e}")
    for sup in beyond:
        rho1 = late[tuple(sup)]
        print(f"    support {sup}, the propagated state at t = {t1:g}:")
        print(f"      max|state - per-sector mixed state|         = {np.abs(rho1 - pred).max():.3e}")
        print(f"      max|state - per-(sector, even part) limit|  = {np.abs(rho1 - refined).max():.3e}")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("prep", "all"):
        run_prep()
    if which in ("support", "all"):
        print()
        run_support()
    if which in ("dimension", "all"):
        print()
        run_dimension()
    if which in ("identity", "all"):
        print()
        run_identity()
    if which in ("branch", "all"):
        print()
        run_branch()
    if which in ("coherence", "all"):
        print()
        run_coherence()
    if which in ("scope", "all"):
        print()
        run_scope()


if __name__ == "__main__":
    main()
