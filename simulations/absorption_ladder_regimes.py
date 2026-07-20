"""Regime map of the absorption ladder (PROOF_ABSORPTION_THEOREM.md).

Reproduces every quantitative statement the proof makes about WHERE its
corollaries hold, as opposed to the identity Re(lambda) = -2*gamma*<n_XY>
itself (that one is verified per mode in ABSORPTION_THEOREM_DISCOVERY.md and
gated in C# by F8PartnerLightComplementarityTests).

The distinction this script exists to pin: some levels of the spectrum are
exact at every coupling, others are only a J/gamma -> infinity limit, and
several published band edges erode at moderate coupling. Each check below
asserts the number that appears in the proof.

Convention: Pauli J=1, H = J * sum_<kl> (X_k X_l + Y_k Y_l + Z_k Z_l) on an
open chain, D_Z(rho) = Z rho Z - rho. This is NOT the C# spin J/4 convention;
rebuilding there gives 4x the imaginary parts.

Run: python simulations/absorption_ladder_regimes.py
"""

import itertools

import numpy as np
import scipy.linalg as sla

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)

CANONICAL_GAMMA = 0.05
CANONICAL_J = 0.075  # Q = J/gamma = 1.5, the repo's canonical regime


def site_op(site, pauli, n):
    """Pauli operator acting on one site of an n-qubit register."""
    out = np.array([[1]], dtype=complex)
    for k in range(n):
        out = np.kron(out, pauli if k == site else I2)
    return out


def heisenberg(n, coupling):
    h = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for k in range(n - 1):
        for pauli in (X, Y, Z):
            h += coupling * site_op(k, pauli, n) @ site_op(k + 1, pauli, n)
    return h


def xy_chain(n, coupling):
    """XY chain: the Heisenberg bond without its ZZ term."""
    h = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for k in range(n - 1):
        for pauli in (X, Y):
            h += coupling * site_op(k, pauli, n) @ site_op(k + 1, pauli, n)
    return h


def liouvillian(n, coupling, gamma, hamiltonian=heisenberg):
    """Vectorized L = -i[H, .] + sum_k gamma * D_{Z_k}."""
    dim = 2 ** n
    h = hamiltonian(n, coupling)
    lio = -1j * (np.kron(h, np.eye(dim)) - np.kron(np.eye(dim), h.T))
    for k in range(n):
        zk = site_op(k, Z, n)
        lio += gamma * (np.kron(zk, zk.T) - np.eye(dim * dim))
    return lio


def rates(n, coupling, gamma, hamiltonian=heisenberg):
    """Decay rates -Re(lambda), in units of gamma, sorted."""
    return np.sort(
        -np.linalg.eigvals(liouvillian(n, coupling, gamma, hamiltonian)).real / gamma)


def threshold(n, hamiltonian=heisenberg, gamma=CANONICAL_GAMMA):
    """Bisect the smallest Q = J/gamma at which the gap reaches 2*gamma."""
    lo, hi = 0.05, 3.0
    r_lo = rates(n, lo * gamma, gamma, hamiltonian)
    r_hi = rates(n, hi * gamma, gamma, hamiltonian)
    assert not r_lo[r_lo > 1e-9][0] > 2.0 - 1e-9, "bracket: gap already 2*gamma at lo"
    assert r_hi[r_hi > 1e-9][0] > 2.0 - 1e-9, "bracket: gap still below 2*gamma at hi"
    for _ in range(40):
        mid = (lo + hi) / 2
        r = rates(n, mid * gamma, gamma, hamiltonian)
        if r[r > 1e-9][0] > 2.0 - 1e-9:
            hi = mid
        else:
            lo = mid
    return hi


def levels(values, tol=1e-6):
    """Collapse a sorted array into (value, multiplicity) pairs."""
    out = []
    for v in values:
        if not out or abs(v - out[-1][0]) > tol:
            out.append([v, 1])
        else:
            out[-1][1] += 1
    return [(v, m) for v, m in out]


def check_n3_ladder():
    """The four pure-weight rungs are exact at every J; multiplicities close to 64."""
    print("N=3 rate ladder (units of gamma)")
    for coupling, gamma in ((CANONICAL_J, CANONICAL_GAMMA), (1.0, 0.05), (1.0, 0.3), (1000.0, 0.3)):
        lv = levels(rates(3, coupling, gamma))
        total = sum(m for _, m in lv)
        assert total == 64
        for rung in (0.0, 2.0, 4.0, 6.0):
            hit = [(v, m) for v, m in lv if abs(v - rung) < 1e-9]
            assert len(hit) == 1, f"rung {rung} missing at J/gamma={coupling / gamma}"
            expected_mult = 4 if rung in (0.0, 6.0) else 14
            assert hit[0][1] == expected_mult, f"rung {rung}: {hit[0][1]} != {expected_mult}"
        print(f"  J/gamma={coupling / gamma:>8.1f}: "
              + "  ".join(f"{v:.4f}x{m}" for v, m in lv))
    print("  rungs 0(x4), 2(x14), 4(x14), 6(x4) exact at every J -> 36 modes;")
    print("  the two mixed bands carry 14 each -> 36 + 28 = 64\n")


def check_fractional_rates_are_a_limit():
    """8/3 and 10/3 are a J/gamma -> infinity limit; each is a triple at finite J."""
    print("The two fractional rates are a limit, not exact rationals")
    gamma = 0.3
    for ratio, split in ((1.5, True), (20, True), (1000, False)):
        r = rates(3, ratio * gamma, gamma)
        band = [v for v, m in levels(r) if 2.2 < v < 3.0]
        mult = sum(m for v, m in levels(r) if 2.2 < v < 3.0)
        assert mult == 14, mult
        if split:
            assert len(band) == 3, "the band is a triple at finite coupling"
        else:
            # at J/gamma = 1000 the triple has merged to within 1e-6: that IS the limit
            assert len(band) == 1, band
        print(f"  J/gamma={ratio:>7.1f}: " + "  ".join(f"{v:.4f}" for v in band))
    print(f"  limit 8/3 = {8 / 3:.4f}")

    # canonical regime: up to 8% away from 8/3
    r = rates(3, CANONICAL_J, CANONICAL_GAMMA)
    band = [v for v, m in levels(r) if 2.2 < v < 3.0]
    worst = max(abs(v - 8 / 3) / (8 / 3) for v in band)
    assert abs(worst - 0.0772) < 5e-4, worst
    print(f"  canonical Q=1.5 band: " + "  ".join(f"{v:.4f}" for v in band)
          + f"  -> up to {worst:.1%} away from 8/3")

    # the splitting closes as 0.46 * (gamma/J)^2
    print("  splitting constant dev * (J/gamma)^2:")
    for ratio in (10, 100, 1000):
        r = rates(3, ratio * gamma, gamma)
        dev = max(abs(v - 8 / 3) for v, _ in levels(r) if 2.2 < v < 3.0)
        const = dev * ratio ** 2
        assert 0.45 < const < 0.47, const
        print(f"    J/gamma={ratio:>5}: {const:.4f}")
    print()


def check_gap_regimes():
    """Gap = 2*gamma above threshold, 2*J^2/gamma below it (Zeno suppression)."""
    print("Spectral gap has two regimes")
    gamma = 0.3
    for coupling in (0.001, 0.01):
        r = rates(3, coupling, gamma) * gamma
        nonzero = r[r > 1e-12]
        zeno = nonzero[0] / (2 * coupling ** 2 / gamma)
        assert abs(zeno - 1.0) < 0.01, zeno
        print(f"  N=3 J={coupling}: gap/(2*gamma)={nonzero[0] / (2 * gamma):.8f}   "
              f"gap/(2J^2/gamma)={zeno:.6f}")

    for n in (3, 4):
        r = rates(n, CANONICAL_J, CANONICAL_GAMMA)
        nonzero = r[r > 1e-9]
        assert abs(nonzero[0] - 2.0) < 1e-9, (n, nonzero[0])
        print(f"  N={n} canonical: smallest nonzero rate = {nonzero[0]:.6f} gamma")
    print()


def check_erosion_is_population_only():
    """Everything that falls below 2*gamma is a POPULATION mode.

    Sorting by |delta popcount| (the excitation-number difference between a
    coherence's bra and ket labels): the sub-2*gamma modes all sit at 0, while
    every coherence sector keeps its ladder rung at every coupling. This bounds
    the reach of the Q*_gap correction: coherent windows are untouched by it.
    """
    print("The gap erosion is confined to the population sector")

    # The two ingredients that make the coherence floor a THEOREM rather than a
    # sample (PROOF_ABSORPTION_THEOREM.md §4.3): L is exactly block-diagonal in
    # |delta popcount| for a number-conserving H, and inside a |dp|=k block every
    # basis coherence has Hamming >= k, so the convex combination <n_XY> >= k and
    # the rate 2*gamma*<n_XY> >= 2*k*gamma.
    for n in (3, 4):
        dim = 2 ** n
        pc = np.array([bin(i).count("1") for i in range(dim)])
        dp = np.abs(pc[:, None] - pc[None, :]).reshape(-1)
        lio = liouvillian(n, CANONICAL_J, CANONICAL_GAMMA)
        mask = dp[:, None] != dp[None, :]
        assert np.max(np.abs(lio[mask])) == 0.0, f"L mixes |dp| sectors at N={n}"
        ham = np.array([[bin(i ^ j).count("1") for j in range(dim)]
                        for i in range(dim)])
        assert np.all(ham >= np.abs(pc[:, None] - pc[None, :])), "Hamming < |dp|"
    print("  L is exactly |dp|-block-diagonal and Hamming >= |dp|: the floor is derived")

    for n in (3, 4, 5):
        dim = 2 ** n
        popcount = np.array([bin(i).count("1") for i in range(dim)])
        delta = np.abs(popcount[:, None] - popcount[None, :]).reshape(-1)
        for q in (0.2, 0.5, 1.5):
            lio = liouvillian(n, q * CANONICAL_GAMMA, CANONICAL_GAMMA)
            eigvals, eigvecs = np.linalg.eig(lio)
            r = -eigvals.real / CANONICAL_GAMMA
            weight = np.abs(eigvecs) ** 2
            weight = weight / weight.sum(axis=0)
            mins = {}
            for k in (0, 1, 2):
                share = weight[delta == k, :].sum(axis=0)
                inside = r[(share > 0.99) & (r > 1e-9)]
                mins[k] = inside.min() if len(inside) else float("nan")
            # the coherence sectors sit exactly on their rungs, at every Q
            assert abs(mins[1] - 2.0) < 1e-9, (n, q, mins[1])
            assert abs(mins[2] - 4.0) < 1e-9, (n, q, mins[2])
            if q < 1.5 or n == 5:
                assert mins[0] < 2.0 - 1e-9, (n, q, mins[0])
            print(f"  N={n} Q={q}: population {mins[0]:.6f}   "
                  f"|dp|=1 {mins[1]:.6f}   |dp|=2 {mins[2]:.6f}")
    print("  so 'the gap is 2*gamma' is regime-bounded, but")
    print("  'a coherence decays no slower than 2*gamma' is not")
    print()


def check_band_erosion_is_palindromic():
    """At N=5, canonical, the F3 band spills at BOTH ends, in palindromic pairs."""
    print("F3's band edges are not universal")
    r = rates(5, CANONICAL_J, CANONICAL_GAMMA)
    nonzero = r[r > 1e-9]
    lo = nonzero[0]
    hi = max(v for v in r if v < 10.0 - 1e-6)
    assert abs(lo - 1.2754) < 1e-3, lo
    assert abs(hi - 8.7246) < 1e-3, hi
    # The live content: the band spills at BOTH ends. (lo + hi = 2N is then
    # automatic from the F1 palindrome, echoed here as a consistency check,
    # not asserted as a finding.)
    assert lo < 2.0 - 1e-3 and hi > 8.0 + 1e-3, "band must spill at both ends"
    assert abs(lo + hi - 10.0) < 1e-9, "F1 palindrome consistency echo"
    print(f"  N=5 canonical: extreme band rates {lo:.4f} and {hi:.4f} gamma")
    print(f"  outside [2, 8] at both ends; sum = {lo + hi:.6f} = 2N exactly")

    # The F3 caveat's numbers pin at Q = 1.0; the edges sharpen as Q grows.
    expected = {(3, 0.5): 0.539, (3, 1.0): 2.000, (3, 1.5): 2.000, (3, 2.0): 2.000,
                (4, 0.5): 0.272, (4, 1.0): 0.978, (4, 1.5): 2.000, (4, 2.0): 2.000,
                (5, 0.5): 0.177, (5, 1.0): 0.617, (5, 1.5): 1.275, (5, 2.0): 2.000}
    print("  min nonzero rate vs coupling ratio Q = J/gamma:")
    for n in (3, 4, 5):
        row = []
        for q in (0.5, 1.0, 1.5, 2.0):
            r = rates(n, q * CANONICAL_GAMMA, CANONICAL_GAMMA)
            lowest = r[r > 1e-9][0]
            assert abs(lowest - expected[(n, q)]) < 1e-3, (n, q, lowest)
            row.append(f"Q={q}: {lowest:.3f}")
        print(f"    N={n}: " + "   ".join(row))
    print()


def check_immortality_is_about_eigenmodes():
    """The {I,Z} sector has 2^N strings but only N+1 frozen modes."""
    print("Immortality is a property of eigenmodes, not of Pauli content")
    n = 3
    kernel = int(np.sum(np.abs(np.linalg.eigvals(
        liouvillian(n, CANONICAL_J, CANONICAL_GAMMA)).real) < 1e-10))
    assert kernel == n + 1, kernel
    print(f"  N={n}: kernel dimension {kernel} = N+1, "
          f"while the (I,Z) sector holds 2^N = {2 ** n} strings")

    z1 = site_op(0, Z, n).reshape(-1)
    for coupling in (0.0, CANONICAL_J, 1.0):
        evolved = sla.expm(liouvillian(n, coupling, CANONICAL_GAMMA) * 50.0) @ z1
        overlap = (z1.conj() @ evolved).real / (z1.conj() @ z1).real
        print(f"  <Z1(50)|Z1(0)> at J={coupling}: {overlap:.6f}")
        if coupling == 0.0:
            assert abs(overlap - 1.0) < 1e-9, "Z1 is exactly frozen at J=0"
        else:
            # the asymptote is exactly 1/N; t=50 has not fully converged
            late = sla.expm(liouvillian(n, coupling, CANONICAL_GAMMA) * 2000.0) @ z1
            share = (z1.conj() @ late).real / (z1.conj() @ z1).real
            assert abs(share - 1 / n) < 1e-9, (coupling, share)
    print(f"  asymptote at J={coupling}: {share:.9f} = 1/N exactly\n")


def check_per_coherence_rate_needs_an_eigenmode():
    """2*gamma*n_diff is a diagonal entry; a decay rate only if |A><B| is an eigenmode."""
    print("The per-coherence rate is a diagonal entry, not always a rate")
    n, dim = 3, 8
    for coupling, expected in ((0.0, 0.0), (CANONICAL_J, 0.30), (1.0, 4.0)):
        rho = np.zeros((dim, dim), dtype=complex)
        rho[0b001, 0b010] = 1.0
        vec = rho.reshape(-1)
        residual = np.linalg.norm(
            liouvillian(n, coupling, CANONICAL_GAMMA) @ vec - (-2 * CANONICAL_GAMMA * 2) * vec)
        assert abs(residual - expected) < 0.01, (coupling, residual)
        print(f"  |001><010| eigen-residual at J={coupling}: {residual:.4f}")
    print()


def check_depolarizing_stays_diagonal():
    """Depolarizing keeps Herm(L) diagonal and obeys Re(lambda) = -4*gamma*<n_nonI>."""
    print("The boundary is Pauli-string jump operators, not dephasing")
    n, gamma = 3, CANONICAL_GAMMA
    dim = 2 ** n

    basis = []
    weights = []
    for combo in itertools.product([I2, X, Y, Z], repeat=n):
        mat = np.array([[1]], dtype=complex)
        for pauli in combo:
            mat = np.kron(mat, pauli)
        basis.append(mat.reshape(-1) / np.sqrt(dim))
        weights.append(sum(1 for p in combo if p is not I2))
    basis = np.array(basis).T
    weights = np.array(weights)

    dissipator = np.zeros((dim * dim, dim * dim), dtype=complex)
    for k in range(n):
        for pauli in (X, Y, Z):
            pk = site_op(k, pauli, n)
            dissipator += gamma * (np.kron(pk, pk.T) - np.eye(dim * dim))

    in_pauli = basis.conj().T @ dissipator @ basis
    off = np.max(np.abs(in_pauli - np.diag(np.diag(in_pauli))))
    assert off < 1e-12, off
    print(f"  depolarizing: max off-diagonal in the Pauli basis = {off:.3e}")

    h = heisenberg(n, CANONICAL_J)
    full = -1j * (np.kron(h, np.eye(dim)) - np.kron(np.eye(dim), h.T)) + dissipator
    eigvals, eigvecs = np.linalg.eig(full)
    coeffs = np.abs(np.linalg.solve(basis, eigvecs)) ** 2
    coeffs /= coeffs.sum(axis=0)
    predicted = -4 * gamma * (weights @ coeffs)
    err = np.max(np.abs(eigvals.real - predicted))
    assert err < 1e-12, err
    print(f"  Re(lambda) = -4*gamma*<n_nonI>: max error = {err:.3e}\n")


def check_zeno_law_is_n_dependent():
    """Below threshold: Delta = 2*(1-cos(pi/N)) * 2*J^2/gamma. The bare form is N=3 only."""
    print("The Zeno gap law carries an N-dependent factor")
    gamma, coupling = 0.3, 0.001
    for n in (2, 3, 4, 5):
        r = -np.linalg.eigvals(liouvillian(n, coupling, gamma)).real
        lowest = np.sort(r[r > 1e-14])[0]
        measured = lowest / (2 * coupling ** 2 / gamma)
        band = 2 * (1 - np.cos(np.pi / n))
        assert abs(measured - band) < 1e-4, (n, measured, band)
        print(f"  N={n}: measured/(2J^2/gamma)={measured:.6f}   2(1-cos(pi/N))={band:.6f}")

    # The threshold must be BISECTED. A coarse grid reports the first grid point
    # above Q*_gap, which is a property of the grid, not of the physics.
    print("  threshold Q*_gap above which the gap equals 2*gamma (bisected):")
    for n, expected_q in ((2, 0.500000), (3, 0.800243), (4, 1.342243), (5, 1.819350)):
        found = threshold(n)
        assert abs(found - expected_q) < 1e-4, (n, found, expected_q)
        print(f"    Heisenberg N={n}: Q*_gap = {found:.6f}")
    print("  not linear in N: successive gaps 0.300, 0.542, 0.477")

    # Near the threshold the asymptote fails in BOTH directions: at 98% of
    # Q*_gap the measured/predicted ratio (gap / Zeno asymptote) undershoots
    # at small N and overshoots from N=4 on (D06's near-threshold table).
    for n, expected_ratio in ((2, 1.67), (3, 1.43), (4, 0.91), (5, 0.77)):
        q = 0.98 * threshold(n)
        r = rates(n, q * CANONICAL_GAMMA, CANONICAL_GAMMA)
        gap_g = r[r > 1e-9][0]
        asym = 2 * (1 - np.cos(np.pi / n)) * 2 * q ** 2
        ratio = gap_g / asym
        assert abs(ratio - expected_ratio) < 0.01, (n, ratio)
        print(f"    98% of Q*_gap, N={n}: gap/asymptote = {ratio:.2f}")

    # Bisection returns A crossing; it returns THE threshold only if the gap is
    # monotone below it and stays at 2*gamma above. Both checked, so the word
    # "threshold" and the universally-quantified "above Q*_gap" downstream
    # statements are earned rather than assumed.
    for n in (2, 3, 4, 5):
        th = threshold(n)
        below = []
        for q in np.linspace(0.05, th * 0.99, 25):
            r = rates(n, q * CANONICAL_GAMMA, CANONICAL_GAMMA)
            below.append(r[r > 1e-9][0])
        assert all(below[i] <= below[i + 1] + 1e-9 for i in range(len(below) - 1)), \
            f"gap not monotone below threshold at N={n}"
        for q in (th * 1.01, 1.5 * th, 3.0, 6.0, 12.0, 30.0):
            r = rates(n, q * CANONICAL_GAMMA, CANONICAL_GAMMA)
            assert abs(r[r > 1e-9][0] - 2.0) < 1e-9, (n, q, r[r > 1e-9][0])
    print("  monotone below the crossing, exactly 2*gamma up to Q=30 above it")

    # The threshold belongs to the Hamiltonian, not to N alone.
    for n, expected_q in ((3, 0.707107), (4, 0.939271), (5, 1.186087)):
        found = threshold(n, xy_chain)
        assert abs(found - expected_q) < 1e-5, (n, found, expected_q)
        print(f"    XY         N={n}: Q*_gap = {found:.6f}")
    assert abs(threshold(3, xy_chain) - 1 / np.sqrt(2)) < 1e-5
    print("  XY N=3 lands on 1/sqrt(2) = 0.707107")

    # Converted to carrier-clock units (x2, GLOSSARY factor-2 convention), the XY
    # gap threshold meets the coherence horizon Q*(N) = 1, sqrt(2), 1.87874, 2.37367
    # (ANALYTICAL_FORMULAS F2b) exactly at N = 2, 3 only, the clean-2x2 tangency
    # where <n_XY> touches 1 at the EP. From N = 4 on the two are distinct objects.
    assert abs(2 * threshold(2, xy_chain) - 1.0) < 1e-6
    assert abs(2 * threshold(3, xy_chain) - np.sqrt(2)) < 1e-6
    assert abs(2 * threshold(4, xy_chain) - 1.87874) > 1e-4, "should be distinct"
    assert abs(2 * threshold(5, xy_chain) - 2.37367) > 1e-3, "should be distinct"
    print("  x2 (carrier units): meets the coherence horizon exactly at N=2,3 only")

    # The asymptote's error, measured against the TRUE gap ((asym - gap)/gap,
    # N = 5; the proof states all drift percentages in this convention):
    # +8% at Q = 0.5, +35% at Q = 1.5; the 10% line is crossed near Q = 0.57.
    print("  Zeno-asymptote error vs the true gap (N=5):")
    for q, expected_err in ((0.5, 0.081), (1.5, 0.348)):
        r = rates(5, q * CANONICAL_GAMMA, CANONICAL_GAMMA)
        true_gap = r[r > 1e-9][0]
        asym = 2 * (1 - np.cos(np.pi / 5)) * 2 * q ** 2
        err = (asym - true_gap) / true_gap
        assert abs(err - expected_err) < 0.005, (q, err)
        print(f"    Q={q}: (asym - gap)/gap = {err:+.1%}")

    # ... so on the XY chain the canonical Q = 1.5 shows no erosion at N=5.
    r = rates(5, 1.5 * CANONICAL_GAMMA, CANONICAL_GAMMA, xy_chain)
    nonzero = r[r > 1e-9]
    top = max(v for v in r if v < 10.0 - 1e-6)
    assert abs(nonzero[0] - 2.0) < 1e-9 and abs(top - 8.0) < 1e-9, (nonzero[0], top)
    print(f"  XY N=5 at Q=1.5: band exactly [{nonzero[0]:.4f}, {top:.4f}] gamma, no erosion")
    print("  canonical Q=1.5: above Q*_gap at Heisenberg N=3,4, below at N=5\n")


def check_frozen_strings_are_two():
    """Only I^N and Z^N are frozen strings; the (N+1)-dim kernel is superpositions."""
    print("Frozen {I,Z} STRINGS vs kernel DIMENSION")
    for n in (2, 3, 4):
        lio = liouvillian(n, CANONICAL_J, CANONICAL_GAMMA)
        frozen = []
        for combo in itertools.product([("I", I2), ("Z", Z)], repeat=n):
            mat = np.array([[1]], dtype=complex)
            label = ""
            for name, pauli in combo:
                mat = np.kron(mat, pauli)
                label += name
            if np.linalg.norm(lio @ mat.reshape(-1)) < 1e-10:
                frozen.append(label)
        kernel = int(np.sum(np.abs(np.linalg.eigvals(lio).real) < 1e-10))
        assert kernel == n + 1 and frozen == ["I" * n, "Z" * n], (n, kernel, frozen)
        print(f"  N={n}: kernel dim {kernel} = N+1, but only {len(frozen)} frozen strings {frozen}")
    print()


def check_parity_forbids_delta_w_one():
    """[L, (-1)^n_XY] = 0 at every J, so a Delta_w = 1 mixture is impossible.

    Aggregate weight shares CANNOT establish this: all four weights carry
    nonzero share, so a share vector is equally consistent with Delta_w = 1.
    The claim is per-eigenvector and needs the commutator plus the parity split.
    """
    print("Parity forbids Delta_w = 1 (the per-eigenvector statement)")
    n, dim = 3, 8
    signs = []
    for combo in itertools.product([I2, X, Y, Z], repeat=n):
        mat = np.array([[1]], dtype=complex)
        for pauli in combo:
            mat = np.kron(mat, pauli)
        signs.append((mat.reshape(-1), (-1) ** sum(1 for p in combo if p is X or p is Y)))
    basis = np.array([v for v, _ in signs]).T / np.sqrt(dim)
    parity = basis @ np.diag([s for _, s in signs]) @ np.linalg.inv(basis)

    for coupling in (0.0, CANONICAL_J, 1.0, 100.0):
        lio = liouvillian(n, coupling, CANONICAL_GAMMA)
        comm = np.max(np.abs(lio @ parity - parity @ lio))
        assert comm < 1e-12, (coupling, comm)
        print(f"  J={coupling:<6}: max |[L, (-1)^n_XY]| = {comm:.2e}")

    # every eigenvector of the fractional cluster is purely even or purely odd
    gamma, coupling = 0.3, 300.0
    eigvals, eigvecs = np.linalg.eig(liouvillian(n, coupling, gamma))
    r = -eigvals.real / gamma
    sel = np.where((r > 2.2) & (r < 3.0))[0]
    even = 0
    for k in sel:
        vec = eigvecs[:, k]
        share = np.abs((parity @ vec + vec) / 2) ** 2
        frac = share.sum() / np.linalg.norm(vec) ** 2
        assert frac < 1e-10 or frac > 1 - 1e-10, frac
        even += frac > 0.5
    assert (even, len(sel) - even) == (10, 4), (even, len(sel))
    print(f"  cluster splits {even} even + {len(sel) - even} odd, none mixed\n")


def check_cluster_shares_are_a_limit():
    """The (5,5,10,1)/21 shares and <n_XY> = 4/3 are a J -> infinity limit."""
    print("The fractional cluster's weight shares (basis-free)")
    n, gamma, coupling = 3, 0.3, 300.0
    dim = 2 ** n
    basis, weights = [], []
    for combo in itertools.product([I2, X, Y, Z], repeat=n):
        mat = np.array([[1]], dtype=complex)
        for pauli in combo:
            mat = np.kron(mat, pauli)
        basis.append(mat.reshape(-1) / np.sqrt(dim))
        weights.append(sum(1 for p in combo if p is X or p is Y))
    basis, weights = np.array(basis).T, np.array(weights)

    eigvals, eigvecs = np.linalg.eig(liouvillian(n, coupling, gamma))
    r = -eigvals.real / gamma
    sel = np.where((r > 2.2) & (r < 3.0))[0]
    assert len(sel) == 14, len(sel)
    ortho, _ = np.linalg.qr(eigvecs[:, sel])
    coeffs = np.linalg.solve(basis, ortho)
    shares = (np.abs(coeffs) ** 2).sum(axis=1) / len(sel)

    per_w = [shares[weights == k].sum() for k in range(4)]
    for k, expected in zip(range(4), (5 / 21, 5 / 21, 10 / 21, 1 / 21)):
        assert abs(per_w[k] - expected) < 1e-5, (k, per_w[k], expected)
    n_xy = sum(k * per_w[k] for k in range(4))
    assert abs(n_xy - 4 / 3) < 1e-6, n_xy
    print(f"  Q={coupling / gamma:<8.0f} shares x21 = "
          + " ".join(f"{21 * s:.3f}" for s in per_w) + f"   <n_XY> = {n_xy:.6f}")

    # ... but at the canonical coupling they are NOT those values.
    eigvals, eigvecs = np.linalg.eig(liouvillian(n, 1.5 * gamma, gamma))
    r = -eigvals.real / gamma
    sel = np.where((r > 2.2) & (r < 3.0))[0]
    ortho, _ = np.linalg.qr(eigvecs[:, sel])
    coeffs = np.linalg.solve(basis, ortho)
    sh = (np.abs(coeffs) ** 2).sum(axis=1) / len(sel)
    can = [sh[weights == k].sum() for k in range(4)]
    can_nxy = sum(k * can[k] for k in range(4))
    assert abs(can_nxy - 1.274489) < 1e-4, can_nxy
    assert abs(can[0] * 21 - 5.309) < 1e-2, can[0] * 21
    print("  Q=1.5     shares x21 = " + " ".join(f"{21 * s:.3f}" for s in can)
          + f"   <n_XY> = {can_nxy:.6f}")
    print("  so 4/3 is a limit, forced by <n_XY> = rate/(2*gamma)\n")


def check_mirror_cluster():
    """The 10*gamma/3 cluster mirrors the 8*gamma/3 one: (1,10,5,5)/21 and 5/3."""
    print("The mirror cluster")
    n, gamma, coupling = 3, 0.3, 300.0
    dim = 2 ** n
    basis, weights = [], []
    for combo in itertools.product([I2, X, Y, Z], repeat=n):
        mat = np.array([[1]], dtype=complex)
        for pauli in combo:
            mat = np.kron(mat, pauli)
        basis.append(mat.reshape(-1) / np.sqrt(dim))
        weights.append(sum(1 for p in combo if p is X or p is Y))
    basis, weights = np.array(basis).T, np.array(weights)

    eigvals, eigvecs = np.linalg.eig(liouvillian(n, coupling, gamma))
    r = -eigvals.real / gamma
    sel = np.where((r > 3.0) & (r < 3.8))[0]
    assert len(sel) == 14, len(sel)
    ortho, _ = np.linalg.qr(eigvecs[:, sel])
    coeffs = np.linalg.solve(basis, ortho)
    shares = (np.abs(coeffs) ** 2).sum(axis=1) / len(sel)
    per_w = [shares[weights == k].sum() for k in range(4)]
    for k, expected in zip(range(4), (1 / 21, 10 / 21, 5 / 21, 5 / 21)):
        assert abs(per_w[k] - expected) < 1e-5, (k, per_w[k], expected)
    n_xy = sum(k * per_w[k] for k in range(4))
    assert abs(n_xy - 5 / 3) < 1e-6, n_xy
    print("  shares x21 = " + " ".join(f"{21 * s:.3f}" for s in per_w)
          + f"   <n_XY> = {n_xy:.6f} = 5/3 (same limit)\n")


def check_band_spread_coefficient():
    """0.46 is the lowest level's approach; the band's own spread closes as 0.53."""
    print("Two different (gamma/J)^2 coefficients")
    gamma = 0.3
    for ratio in (50, 100, 200, 400):
        r = rates(3, ratio * gamma, gamma)
        band = [v for v, _ in levels(r, tol=1e-9) if 2.2 < v < 3.0]
        spread = (max(band) - min(band)) * ratio ** 2
        approach = (8 / 3 - min(band)) * ratio ** 2
        assert abs(spread - 0.527) < 0.005 and abs(approach - 0.461) < 0.005
        print(f"  Q={ratio:>4}: spread*Q^2={spread:.4f}   (8/3 - lowest)*Q^2={approach:.4f}")
    print()


def check_section4_needs_the_number_conserving_family():
    """Kernel dim, drain, and the 2N*gamma endpoint are H-dependent, not theorem consequences."""
    print("Section 4's shape belongs to the number-conserving family")
    n, gamma, dim = 3, CANONICAL_GAMMA, 8

    def lio_of(h):
        out = -1j * (np.kron(h, np.eye(dim)) - np.kron(np.eye(dim), h.T))
        for k in range(n):
            zk = site_op(k, Z, n)
            out += gamma * (np.kron(zk, zk.T) - np.eye(dim * dim))
        return out

    rng = np.random.default_rng(7)
    a = rng.standard_normal((dim, dim))
    ising = sum(0.3 * site_op(k, Z, n) @ site_op(k + 1, Z, n) for k in range(n - 1))
    cases = {
        "Heisenberg": heisenberg(n, CANONICAL_J),
        "random symm": (a + a.T) / 2,
        "Ising ZZ": ising,
        "H = 0": np.zeros((dim, dim), dtype=complex),
    }
    seen = {}
    for name, h in cases.items():
        r = -np.linalg.eigvals(lio_of(h)).real / gamma
        kernel = int(np.sum(np.abs(r) < 1e-8))
        drain = int(np.sum(np.abs(r - 2 * n) < 1e-8))
        seen[name] = (r.max(), kernel, drain)
        print(f"  {name:<12} max rate {r.max():.4f} gamma (2N = {2 * n})  kernel {kernel}  drain {drain}")
    assert seen["Heisenberg"][1:] == (n + 1, n + 1)
    assert seen["Ising ZZ"][1:] == (2 ** n, 2 ** n)
    assert seen["random symm"][0] < 2 * n - 1 and seen["random symm"][1] == 1
    print("  so N+1 and the 2N*gamma endpoint come from F4/F1, not from the theorem\n")


def check_identity_is_bendixson():
    """Re(lambda) = Rayleigh quotient of Herm(L) holds even where the rate formula fails."""
    print("The Rayleigh identity is universal; only the Pauli diagonal is special")
    n, gamma, gamma_t1, dim = 3, CANONICAL_GAMMA, 0.02, 8
    lowering = np.array([[0, 1], [0, 0]], dtype=complex)
    lio = liouvillian(n, CANONICAL_J, gamma)
    for k in range(n):
        sk = site_op(k, lowering, n)
        lio = lio + gamma_t1 * (
            np.kron(sk, sk.conj())
            - 0.5 * np.kron(sk.conj().T @ sk, np.eye(dim))
            - 0.5 * np.kron(np.eye(dim), (sk.conj().T @ sk).T))
    eigvals, eigvecs = np.linalg.eig(lio)
    herm = (lio + lio.conj().T) / 2
    err = max(
        abs(eigvals[i].real
            - (eigvecs[:, i].conj() @ herm @ eigvecs[:, i]).real / np.linalg.norm(eigvecs[:, i]) ** 2)
        for i in range(len(eigvals)))
    assert err < 1e-12, err
    print(f"  under amplitude damping: max |Re lambda - Rayleigh(Herm L)| = {err:.3e}")

    # ... while the -2*gamma*<n_XY> READING misses by a finite amount there.
    basis, weights = [], []
    for combo in itertools.product([I2, X, Y, Z], repeat=n):
        mat = np.array([[1]], dtype=complex)
        for pauli in combo:
            mat = np.kron(mat, pauli)
        basis.append(mat.reshape(-1) / np.sqrt(dim))
        weights.append(sum(1 for p in combo if p is X or p is Y))
    basis, weights = np.array(basis).T, np.array(weights)
    coeffs = np.abs(np.linalg.solve(basis, eigvecs)) ** 2
    coeffs /= coeffs.sum(axis=0)
    miss = np.max(np.abs(eigvals.real + 2 * gamma * (weights @ coeffs)))
    assert abs(miss - 0.06) < 5e-3, miss
    print(f"  but the -2*gamma*<n_XY> reading misses by {miss:.4f} (gamma_T1 = {gamma_t1})\n")


def check_the_fence_is_one_string_not_locality():
    """The fence is ONE Pauli string per jump operator, not product-vs-correlated.

    Both directions measured: a strictly local (X_k+Z_k)/sqrt2 channel BREAKS
    the reading because each jump is a sum on its own site, while a two-site
    Z_k Z_k+1 channel KEEPS Herm(L) exactly diagonal. Correlated is fine; a sum
    is not.
    """
    print("The fence is one Pauli string per jump, not locality")
    n, gamma, dim = 3, CANONICAL_GAMMA, 8
    h = heisenberg(n, CANONICAL_J)
    basis, weights = [], []
    for combo in itertools.product([I2, X, Y, Z], repeat=n):
        mat = np.array([[1]], dtype=complex)
        for pauli in combo:
            mat = np.kron(mat, pauli)
        basis.append(mat.reshape(-1) / np.sqrt(dim))
        weights.append(sum(1 for p in combo if p is X or p is Y))
    basis, weights = np.array(basis).T, np.array(weights)

    def off_diagonal(jumps):
        lio = -1j * (np.kron(h, np.eye(dim)) - np.kron(np.eye(dim), h.T))
        for jump in jumps:
            lio = lio + gamma * (
                np.kron(jump, jump.conj())
                - 0.5 * np.kron(jump.conj().T @ jump, np.eye(dim))
                - 0.5 * np.kron(np.eye(dim), (jump.conj().T @ jump).T))
        herm = (lio + lio.conj().T) / 2
        in_pauli = basis.conj().T @ herm @ basis
        return np.max(np.abs(in_pauli - np.diag(np.diag(in_pauli))))

    local_sum = [(site_op(k, X, n) + site_op(k, Z, n)) / np.sqrt(2) for k in range(n)]
    two_site = [site_op(k, Z, n) @ site_op(k + 1, Z, n) for k in range(n - 1)]
    local_z = [site_op(k, Z, n) for k in range(n)]

    off_local_sum, off_two_site, off_local_z = (
        off_diagonal(local_sum), off_diagonal(two_site), off_diagonal(local_z))
    assert abs(off_local_sum - 0.05) < 1e-9, off_local_sum
    assert off_two_site < 1e-12, off_two_site
    assert off_local_z < 1e-12, off_local_z
    print(f"  local (X+Z)/sqrt2, a PRODUCT channel: off-diagonal {off_local_sum:.4f} -> breaks")
    print(f"  Z_k Z_k+1, a CORRELATED channel:      off-diagonal {off_two_site:.1e} -> holds")
    print(f"  local Z_k (baseline):                 off-diagonal {off_local_z:.1e} -> holds")
    print()


def check_collective_dephasing_breaks_the_reading():
    """L = sum_k Z_k is pure dephasing, population-preserving, and still breaks it."""
    print("Collective dephasing: pure, Z-only, and outside the fence")
    n, gamma, dim = 3, CANONICAL_GAMMA, 8
    rng = np.random.default_rng(5)
    a = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    h = (a + a.conj().T) / 2
    collective = sum(site_op(k, Z, n) for k in range(n))

    lio = -1j * (np.kron(h, np.eye(dim)) - np.kron(np.eye(dim), h.T))
    lio = lio + gamma * (np.kron(collective, collective.T)
                         - 0.5 * np.kron(collective @ collective, np.eye(dim))
                         - 0.5 * np.kron(np.eye(dim), (collective @ collective).T))

    basis, weights = [], []
    for combo in itertools.product([I2, X, Y, Z], repeat=n):
        mat = np.array([[1]], dtype=complex)
        for pauli in combo:
            mat = np.kron(mat, pauli)
        basis.append(mat.reshape(-1) / np.sqrt(dim))
        weights.append(sum(1 for p in combo if p is X or p is Y))
    basis, weights = np.array(basis).T, np.array(weights)

    herm = (lio + lio.conj().T) / 2
    in_pauli = basis.conj().T @ herm @ basis
    off = np.max(np.abs(in_pauli - np.diag(np.diag(in_pauli))))
    assert abs(off - 0.2) < 1e-9, off
    print(f"  max off-diagonal of Herm(L) in the Pauli basis = {off:.4f} (not diagonal)")

    eigvals, eigvecs = np.linalg.eig(lio)
    coeffs = np.abs(np.linalg.solve(basis, eigvecs)) ** 2
    coeffs /= coeffs.sum(axis=0)
    miss = np.max(np.abs(eigvals.real + 2 * gamma * (weights @ coeffs)))
    assert miss > 1e-2, miss
    print(f"  max |Re lambda + 2*gamma*<n_XY>| = {miss:.4f}: the reading fails\n")


def check_topology_dependence_starts_at_n3():
    """Chain and triangle already differ at N=3, in levels AND rung multiplicities."""
    print("Topology-dependence starts at N=3")

    def lio_bonds(bonds):
        dim = 2 ** 3
        h = np.zeros((dim, dim), dtype=complex)
        for i, j in bonds:
            for pauli in (X, Y, Z):
                h += CANONICAL_J * site_op(i, pauli, 3) @ site_op(j, pauli, 3)
        out = -1j * (np.kron(h, np.eye(dim)) - np.kron(np.eye(dim), h.T))
        for k in range(3):
            zk = site_op(k, Z, 3)
            out += CANONICAL_GAMMA * (np.kron(zk, zk.T) - np.eye(dim * dim))
        return out

    for name, bonds, expected in (
            ("chain", [(0, 1), (1, 2)], (4, 14, 14, 4)),
            ("triangle", [(0, 1), (1, 2), (0, 2)], (4, 16, 16, 4))):
        lv = levels(np.sort(-np.linalg.eigvals(lio_bonds(bonds)).real / CANONICAL_GAMMA))
        rungs = tuple(m for v, m in lv if min(abs(v - t) for t in (0, 2, 4, 6)) < 1e-9)
        assert rungs == expected, (name, rungs)
        print(f"  {name:<9} rung multiplicities {rungs}   levels: "
              + "  ".join(f"{v:.4f}" for v, _ in lv))
    print("  the 2.4607 level exists on the chain and not on the triangle\n")


if __name__ == "__main__":
    check_n3_ladder()
    check_fractional_rates_are_a_limit()
    check_band_spread_coefficient()
    check_parity_forbids_delta_w_one()
    check_cluster_shares_are_a_limit()
    check_mirror_cluster()
    check_gap_regimes()
    check_zeno_law_is_n_dependent()
    check_band_erosion_is_palindromic()
    check_erosion_is_population_only()
    check_immortality_is_about_eigenmodes()
    check_frozen_strings_are_two()
    check_per_coherence_rate_needs_an_eigenmode()
    check_depolarizing_stays_diagonal()
    check_identity_is_bendixson()
    check_the_fence_is_one_string_not_locality()
    check_collective_dephasing_breaks_the_reading()
    check_section4_needs_the_number_conserving_family()
    check_topology_dependence_starts_at_n3()
    print("All regime checks passed.")
