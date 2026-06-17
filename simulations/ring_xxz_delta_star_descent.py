"""Ring XXZ handover Delta*(N): the periodic twin of the open-chain Delta* descent.

Companion to simulations/xxz_delta_star_descent.py (the open chain). Same machinery, one
bond added (the wrap N-1 <-> 0). The headline: the ring Delta*(N) is NON-MONOTONE -- both
parities rise to a peak ~1.31-1.33 near N=9-10 then DESCEND through N=14; the ring crosses
ABOVE the chain near N=7-8; N=4 has NO handover (a full-block Delta=0 tangency); and the
N->inf limit is OPEN at N<=14. Qualitatively unlike the chain's monotone descent to Delta=1.

The frame (reflections/ON_THE_ONE_DIAGONAL.md): the band-edge floor 2*gamma is the first rung
of the one diagonal popcount(i^j) that the light ever touches -- universal, topology-free (the
light's question). Delta*(N) is the Hamiltonian's argument about that fixed floor -- hence
topology-dependent. The diagonal is one; the climb is many.

gamma->0 reduction: Delta* <=> gap(R)=2, R the Z-coupled classical rate matrix among the
half-filling XXZ eigenstates (a Pauli/Fermi-golden-rule relaxation). Convention matches
xxz_delta_star.py: H = J*sum_<ij>(X_i X_j + Y_i Y_j) + Delta*sum_<ij> Z_i Z_j; hopping 2J;
uniform Z-dephasing gamma; coherence |a><b| decays at 2*gamma*popcount(a^b).

Self-contained (no framework import), folding the Round-1/Round-2 review probes.
"""
import numpy as np
from itertools import combinations

# ring Delta* (gamma->0 reduction). N=4 has NO handover (full-block Delta=0 tangency, see
# check_n4_no_handover). Computed once via descent_sequence(14); regenerate with that call.
RING_DSTAR_SEQUENCE = {5: 1.25633, 6: 1.18614, 7: 1.32111, 8: 1.29481, 9: 1.33095,
                       10: 1.30791, 11: 1.31596, 12: 1.29608, 13: 1.29499, 14: 1.27816}
# N values cheap enough to recompute live every run (high-N grid scans are too slow for main).
FAST_REVAL = (5, 6, 7, 8)
# chain Delta* (from xxz_delta_star_descent.DSTAR_SEQUENCE, the landed open chain) for the
# ring-vs-chain crossing assert: ring < chain for N<=7, ring > chain for N>=8.
CHAIN_DSTAR = {4: 1.619612, 5: 1.527984, 6: 1.384629, 7: 1.330070, 8: 1.272426,
               9: 1.247380, 10: 1.215778, 11: 1.199583, 12: 1.179327, 13: 1.168273, 14: 1.153892}


def bonds(N, ring):
    b = [(i, i + 1) for i in range(N - 1)]
    if ring and N > 2:
        b.append((N - 1, 0))  # the wrap bond
    return b


def Hp(N, p, Delta, ring, J=1.0):
    """p-excitation block of H, real symmetric."""
    states = [sum(1 << i for i in c) for c in combinations(range(N), p)]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for a, b in bonds(N, ring):
        for s in states:
            za = 1 - 2 * ((s >> a) & 1)
            zb = 1 - 2 * ((s >> b) & 1)
            H[idx[s], idx[s]] += Delta * za * zb
            if (s >> a) & 1 and not (s >> b) & 1:
                H[idx[(s & ~(1 << a)) | (1 << b)], idx[s]] += 2 * J
            if (s >> b) & 1 and not (s >> a) & 1:
                H[idx[(s & ~(1 << b)) | (1 << a)], idx[s]] += 2 * J
    return H, states


def reduction_R(N, Delta, ring):
    """gamma-free classical rate matrix R among half-filling XXZ eigenstates.
    R_ab = sum_k |<E_a|Z_k|E_b>|^2 (a!=b, gain); R_aa = -sum_{b!=a} R_ba (= -4 sum_k Var_a(n_k)).
    The Z<->n factor of 4 is load-bearing (Var(Z_k)=4 Var(n_k)); harmonizing the gain to n_k
    silently breaks the zero column sums by 4x (check_generator catches it)."""
    p = (N + 1) // 2
    H, states = Hp(N, p, Delta, ring)
    w, V = np.linalg.eigh(H)
    dim = len(states)
    R = np.zeros((dim, dim))
    for k in range(N):
        z = np.array([1 - 2 * ((s >> k) & 1) for s in states], float)
        Zke = V.T @ (z[:, None] * V)
        R += Zke ** 2
    np.fill_diagonal(R, 0.0)
    R -= np.diag(R.sum(axis=0))
    return R, w, states, V


def slow_rate(R, tol=1e-9):
    ev = np.linalg.eigvalsh(R)
    nz = ev[ev < -tol]
    return -nz.max() if len(nz) else 0.0


def gapR(N, Delta, ring):
    R, _, _, _ = reduction_R(N, Delta, ring)
    return slow_rate(R)


def _bisect_cross(N, ring, lo, hi, tol=1e-7):
    f = lambda D: gapR(N, D, ring) - 2
    flo = f(lo)
    while hi - lo > tol:
        m = 0.5 * (lo + hi)
        fm = f(m)
        if (fm > 0) == (flo > 0):
            lo, flo = m, fm
        else:
            hi = m
    return 0.5 * (lo + hi)


def crossings_of_2(N, ring, lo=-1.0, hi=3.0, step=0.01):
    """Scan gap(R) over [lo,hi]; return ALL downward & upward crossings of the floor 2, and the
    grid max. Non-monotone-aware: does NOT assume a single Neel-side crossing (Round-1's bug)."""
    grid = np.arange(lo, hi + 1e-9, step)
    g = np.array([gapR(N, D, ring) for D in grid])
    downs, ups = [], []
    for i in range(len(grid) - 1):
        a, b = g[i] - 2, g[i + 1] - 2
        if a > 0 >= b:
            downs.append(_bisect_cross(N, ring, grid[i], grid[i + 1]))
        elif a <= 0 < b:
            ups.append(_bisect_cross(N, ring, grid[i], grid[i + 1]))
    return downs, ups, float(g.max())


def delta_star_ring(N, lo=1.1, hi=1.5, step=0.01):
    """Ring Delta* = the Neel-side downward crossing of the floor. None if no crossing."""
    downs, _, _ = crossings_of_2(N, ring=True, lo=lo, hi=hi, step=step)
    return downs[0] if downs else None


def descent_sequence(N_max, N_min=5):
    """Regenerate RING_DSTAR_SEQUENCE (live, expensive at high N). N=4 excluded (no handover)."""
    return {N: delta_star_ring(N) for N in range(N_min, N_max + 1)}


def leb_full_block(N, Delta, ring, gamma):
    """Finite-gamma full (p,p) half-filling coherence-block Lebensader rate (the cross-check truth)."""
    p = (N + 1) // 2
    H, states = Hp(N, p, Delta, ring)
    n = len(states)
    L = -1j * (np.kron(H, np.eye(n)) - np.kron(np.eye(n), H.T))
    deph = np.array([-2.0 * gamma * bin(states[a] ^ states[b]).count("1")
                     for a in range(n) for b in range(n)])
    ev = np.linalg.eigvals(L + np.diag(deph))
    nz = ev[ev.real < -1e-12]
    return -nz.real.max() if len(nz) else 0.0


def full_slowest_rate(N, Delta, ring, gamma):
    """Finite-gamma slowest rate of the FULL 4^N Liouvillian (no sector reduction)."""
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1]).astype(complex)

    def site(op, l):
        m = np.array([[1.0 + 0j]])
        for k in range(N):
            m = np.kron(m, op if k == l else I2)
        return m

    d = 2 ** N
    H = np.zeros((d, d), complex)
    for i, j in bonds(N, ring):
        H += site(X, i) @ site(X, j) + site(Y, i) @ site(Y, j) + Delta * (site(Z, i) @ site(Z, j))
    L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    deph = np.array([-2.0 * gamma * bin(a ^ b).count("1") for a in range(d) for b in range(d)])
    ev = np.linalg.eigvals(L + np.diag(deph))
    rate = -ev.real
    nz = rate[rate > 1e-9]
    return nz.min()


def magnon_energies(N, ring, J=1.0):
    H, _ = Hp(N, 1, 0.0, ring, J)  # Delta=0 single magnon = tight-binding
    return np.sort(np.linalg.eigvalsh(H))


# ----------------------------- self-validation -----------------------------

def check_generator():
    """[3] ring R is a valid generator: real-symmetric, zero column sums, single zero mode; and
    the Z<->n factor of 4 (-R_aa == 4 sum_k Var_a(n_k))."""
    for N in (5, 6):
        R, _, states, V = reduction_R(N, 1.3, ring=True)
        assert np.max(np.abs(R - R.T)) < 1e-12, f"N={N}: R not symmetric"
        assert np.max(np.abs(R.sum(axis=0))) < 1e-10, f"N={N}: R columns not zero-sum"
        assert int(np.sum(np.abs(np.linalg.eigvalsh(R)) < 1e-10)) == 1, f"N={N}: not one zero mode"
        var4 = np.zeros(len(states))
        for k in range(N):
            nk = np.array([(s >> k) & 1 for s in states], float)
            ne = (V.T @ (nk[:, None] * V)).diagonal()
            var4 += ne * (1 - ne)
        assert np.max(np.abs(-np.diag(R) - 4 * var4)) < 1e-9, f"N={N}: factor-4 broken"
    print("[3] ring R is a valid generator (sym, zero col-sum, one zero mode, Z<->n factor 4).  OK")


def check_momenta():
    """[6a] ring single-magnon spectrum = 2*(2J)*cos(2pi k/N) (the 2cos(2pi/N) momentum family),
    NOT the chain's 2cos(pi/(N+1))."""
    for N in (4, 5, 6, 7):
        got = magnon_energies(N, ring=True)
        pred = np.sort([4.0 * np.cos(2 * np.pi * k / N) for k in range(N)])
        assert np.max(np.abs(got - pred)) < 1e-10, f"ring N={N} momenta off the 2cos(2pi k/N) family"
    print("[6a] ring single-magnon momenta = 2cos(2pi k/N) family (bit-exact).  OK")


def check_floor_universal():
    """[1] the band-edge floor sits at Re=-2*gamma EXACTLY for every Delta (the one diagonal's
    first rung; topology- and Delta-independent). Full 4^N ring L at N=4."""
    gamma = 0.05
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1]).astype(complex)

    def site(op, l):
        m = np.array([[1.0 + 0j]])
        for k in range(4):
            m = np.kron(m, op if k == l else I2)
        return m

    d = 2 ** 4
    for Delta in (0.0, 0.5, 1.0, 1.5, 2.0):
        H = np.zeros((d, d), complex)
        for i, j in bonds(4, ring=True):
            H += site(X, i) @ site(X, j) + site(Y, i) @ site(Y, j) + Delta * (site(Z, i) @ site(Z, j))
        L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
        deph = np.array([-2.0 * gamma * bin(a ^ b).count("1") for a in range(d) for b in range(d)])
        ev = np.linalg.eigvals(L + np.diag(deph))
        assert np.min(np.abs(ev.real + 2 * gamma)) < 1e-9, f"Delta={Delta}: no mode exactly at -2g"
    print("[1] band-edge floor at Re=-2g exactly for all Delta (the one diagonal, topology-free).  OK")


def check_sector_hosts_lebensader():
    """[2] above Delta*, the ring half-filling (p,p) block's slowest rate IS the full 4^N ring
    Liouvillian's global slowest (bit-exact) -- i.e. the half-filling sector hosts the ring
    Lebensader (the load-bearing premise the Q-axis (2,2)-vs-(1,1) caution put in doubt)."""
    gamma = 0.05
    for N, Delta in [(4, 2.2), (5, 2.0)]:  # above Delta* (deep Neel); N=4 has no handover but the
        sec = leb_full_block(N, Delta, ring=True, gamma=gamma)   # half-filling block is still the
        full = full_slowest_rate(N, Delta, ring=True, gamma=gamma)  # global slowest at large Delta.
        assert abs(sec - full) < 1e-9, f"N={N}: half-filling sector {sec} != full 4^N global slowest {full}"
    print("[2] ring half-filling sector hosts the Lebensader (sector == full 4^N, bit-exact, above Delta*).  OK")


def check_reduction_certified():
    """[1]+[4] the gamma->0 reduction is certified against the full block where the half-filling
    sector hosts the Lebensader: ratio gamma*gap(R)/full -> 1 at the GENERIC N (5,6,7,8), and is
    STUCK ~0.985 at N=4 (the K2,2 accidental degeneracy beyond the generic translation pairs --
    why N=4 must be read off the full block, check_n4_no_handover, not the reduction). The generic
    +-k momentum degeneracy is benign (ratio -> 1); only accidental degeneracy breaks it."""
    for N, dstar in [(5, 1.25633), (6, 1.18614), (7, 1.32111), (8, 1.29481)]:
        g0 = gapR(N, dstar, ring=True)
        ratios = [leb_full_block(N, dstar, ring=True, gamma=g) / (g * g0) for g in (0.05, 0.025, 0.0125)]
        assert abs(ratios[-1] - 1.0) < 1e-3, f"N={N}: reduction not certified vs full block: {ratios}"
        assert abs(ratios[-1] - 1.0) < abs(ratios[0] - 1.0) + 1e-9, f"N={N}: ratio not converging: {ratios}"
    # N=4 is the documented exception: the reduction is ~1.5% off (cert. ratio sticks ~0.985).
    g0_4 = gapR(4, 0.28868, ring=True)
    r4 = leb_full_block(4, 0.28868, ring=True, gamma=0.0125) / (0.0125 * g0_4)
    assert 0.97 < r4 < 0.995, f"N=4 reduction ratio expected ~0.985 (K2,2 accidental degeneracy), got {r4}"
    print("[1/4] reduction certified vs full block at N=5,6,7,8 (ratio->1); N=4 stuck ~0.985 (K2,2).  OK")


def check_n4_no_handover():
    """[5a] ring N=4 has NO handover: the FULL half-filling block (not the reduction) is tangent to
    the floor at the XY point Delta=0 (peak ~0.99998*2g, gamma-converged) and a survivor for every
    other Delta -- it never reaches the floor. (The reduction spuriously 'crosses' at Delta~0.289,
    the K2,2 artifact -- so N=4 is read off the full block.)"""
    gamma = 0.0125
    rates = {D: leb_full_block(4, D, ring=True, gamma=gamma) / (2 * gamma)
             for D in (-0.5, -0.289, 0.0, 0.289, 0.5, 1.0, 1.5)}
    assert max(rates.values()) < 1.0, f"ring N=4 should never reach the floor; got {rates}"
    assert rates[0.0] == max(rates.values()) and rates[0.0] > 0.9999, f"N=4 peak not the Delta=0 tangency: {rates}"
    print(f"[5a] ring N=4 no handover: full-block tangent to the floor at Delta=0 (peak {rates[0.0]:.5f}*2g).  OK")


def check_sequence_and_shape():
    """[5b/5c] the recorded ring Delta*(N) is reproduced live (FAST_REVAL), is NON-MONOTONE (both
    parities hump then descend, odd peak at N=9, even peak at N=10), and the ring crosses ABOVE the
    chain (ring < chain for N<=7, ring > chain for N>=8)."""
    for N in FAST_REVAL:
        got = delta_star_ring(N)
        assert abs(got - RING_DSTAR_SEQUENCE[N]) < 1e-3, f"recorded ring Delta*({N})={RING_DSTAR_SEQUENCE[N]} != live {got}"
    odd = {n: RING_DSTAR_SEQUENCE[n] for n in (5, 7, 9, 11, 13)}
    even = {n: RING_DSTAR_SEQUENCE[n] for n in (6, 8, 10, 12, 14)}
    assert max(odd, key=odd.get) == 9, f"odd peak expected at N=9, got {max(odd, key=odd.get)}"
    assert max(even, key=even.get) == 10, f"even peak expected at N=10, got {max(even, key=even.get)}"
    assert odd[13] < odd[9] and even[14] < even[10], "parities should descend after the peak"
    for N in range(5, 8):
        assert RING_DSTAR_SEQUENCE[N] < CHAIN_DSTAR[N], f"N={N}: ring should be below chain"
    for N in range(8, 15):
        assert RING_DSTAR_SEQUENCE[N] > CHAIN_DSTAR[N], f"N={N}: ring should be above chain"
    print("[5b/c] ring Delta*(N) non-monotone (odd peak N=9, even peak N=10), crosses above chain at N=8.  OK")


def check_limit_open():
    """[6b/7] the N->inf limit is OPEN at N<=14: a power-law fit Delta* = L + a*N^(-alpha) to the
    humped sequence DEGENERATES (|a| huge, alpha huge), so it carries no limit. Do not report it."""
    from scipy.optimize import curve_fit
    model = lambda N, L, a, al: L + a * np.power(N, -al)
    for par in (0, 1):
        Ns = np.array([n for n in RING_DSTAR_SEQUENCE if n % 2 == par], float)
        ys = np.array([RING_DSTAR_SEQUENCE[int(n)] for n in Ns])
        popt, _ = curve_fit(model, Ns, ys, p0=[1.3, -1.0, 1.5], maxfev=20000)
        assert abs(popt[1]) > 1e6 or popt[2] > 10, \
            f"parity {par}: power-law fit did NOT degenerate (a={popt[1]:.2e}, alpha={popt[2]:.1f}); re-examine"
    print("[6b/7] power-law fit to the hump degenerates (alpha~33) -> N->inf limit OPEN at N<=14.  OK")


def verdict():
    print("\n=== VERDICT (ring Delta*(N), the periodic twin) ===")
    print("  Ring Delta*(N) is NON-MONOTONE: both parities hump to ~1.31-1.33 near N=9-10, then")
    print("  descend through N=14 (odd 1.331@9 -> 1.295@13; even 1.308@10 -> 1.278@14). The ring")
    print("  crosses ABOVE the chain near N=7-8 (the chain keeps descending to Delta=1; the ring humps).")
    print("  N=4 has NO handover (full-block tangent to the floor at the XY point Delta=0). The N->inf")
    print("  limit is OPEN; EXTENDED to N=15 (Delta*=1.27413, simulations/ring_delta_star_extend.py): the")
    print("  odd descent CONTINUES with non-shrinking steps (-.015,-.021,-.021) -- no plateau deceleration,")
    print("  disfavoring 'settles above 1', but Delta*=1.274 is still far from 1 (resolving 1-vs-settles")
    print("  needs N>>16, infeasible). New evidence shifts the weight; it does not close the question.")
    print("  FRAME (reflections/ON_THE_ONE_DIAGONAL.md): the floor 2g is the one diagonal's first rung")
    print("  -- the light's question, topology-free. Delta*(N) is the Hamiltonian's argument about that")
    print("  fixed floor -- topology-dependent. The diagonal is one; the climb is many. The chain")
    print("  descends to the critical point; the ring's lattice argues a hump about the same floor.")


if __name__ == "__main__":
    check_generator()
    check_momenta()
    check_floor_universal()
    check_sector_hosts_lebensader()
    check_reduction_certified()
    check_n4_no_handover()
    check_sequence_and_shape()
    check_limit_open()
    verdict()
