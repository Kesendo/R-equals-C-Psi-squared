"""Which half of the Coulson-Rushbrooke palindrome does the CENTRE carry?

The universal palindrome condition (hypotheses/UNIVERSAL_PALINDROME_CONDITION.md,
Open Question 7, answered 2026-06-01) splits a generator X = A + B against an
involution Q into two independent sub-conditions:

    (A-condition)  Q A Q^-1 = -A                 the coupling part is mirror-odd
    (B-condition)  B + Q B Q^-1 = -2c * I        the bath part pairs to a constant

Together they force Q X Q^-1 = -X - 2c*I, hence a spectrum palindromic about -c.

For Hueckel theory X = B + beta*A_graph with Q the bipartite sublattice sign flip.
In an ALL-CARBON molecule B = alpha*I, so -c = alpha and the B-condition is free:
alpha*I commutes with every invertible Q, and the condition collapses to
2*alpha*I = -2c*I with no reference to the graph. That is the whole content of
"the centre is tr/dim", and it is why the centre identification in
docs/carbon/BENZENE_HUCKEL_FRAMEWORK_LENS.md is graded below the pair sums: for
the molecules in that document it carries no structural information at all.

The B-condition is not empty in general, and this script gates both breaks, which
are PROOF_F103_F87_Z2_CUBED_REFINEMENT.md section 7.2's two ways to lose
bipartiteness:

  * an odd cycle breaks the A-condition (C3, C5, C7 below): the graph is not
    two-colourable, the palindrome goes, and the centre survives untouched.
  * a heteroatom breaks the B-condition (pyridine below): the graph is still
    bipartite and the A-condition still exact, but alpha_N = alpha + h_N*beta
    lifts the diagonal and the palindrome goes anyway.

What the B-condition asks is that the diagonal be ODD UNDER Q about the centre.
"Constant diagonal" is only what that reduces to when Q is diagonal, as the
sublattice flip is, and the push-pull chain below is the case that separates the
two: a diagonal running alpha+0.5*beta down to alpha-0.5*beta, a large B-residual
against the flip, and an exact palindrome under Q = P*K.

Framework side: PROOF_ABSORPTION_THEOREM.md section 4.4 gets the spectral mean
from Tr(L_H) = 0 AND Tr(L_D) = -gamma*N*d^2, and notes that the mean and the
PAIRING coincide only when the spectrum is already palindromic. Sign convention:
section 4.4 works in DECAY RATES, where the centre is +Sigma(gamma); in
eigenvalue coordinates it is -Sigma(gamma).

Referenced from docs/carbon/BENZENE_HUCKEL_FRAMEWORK_LENS.md.
"""

import collections
import itertools

import numpy as np
from scipy.optimize import linear_sum_assignment

ALPHA = -11.4  # on-site Coulomb integral, eV
BETA = -2.4  # resonance integral, eV
H_NITROGEN = 0.5  # textbook Hueckel heteroatom parameter, alpha_N = alpha + h*beta


def ring(n):
    a = np.zeros((n, n))
    for i in range(n):
        a[i, (i + 1) % n] = a[(i + 1) % n, i] = 1.0
    return a


def chain(n):
    a = np.zeros((n, n))
    for i in range(n - 1):
        a[i, i + 1] = a[i + 1, i] = 1.0
    return a


def two_colour(adj):
    """BFS sublattice sign flip Q = diag(+-1); returns (Q, is_bipartite).

    On a bipartite graph this is THE sublattice flip. On an odd cycle no valid
    colouring exists and BFS returns one with a single frustrated edge, which
    happens to be the minimum; min_a_residual() below quantifies over all 2^n
    sign patterns so the A-row does not rest on that accident.
    """
    n = len(adj)
    colour = {0: 1}
    queue = collections.deque([0])
    bipartite = True
    while queue:
        u = queue.popleft()
        for v in range(n):
            if adj[u, v]:
                if v not in colour:
                    colour[v] = -colour[u]
                    queue.append(v)
                elif colour[v] == colour[u]:
                    bipartite = False
    q = np.diag([float(colour.get(i, 1)) for i in range(n)])
    return q, bipartite


def carbon_bath(n, alpha=ALPHA):
    return alpha * np.eye(n)


def hetero_bath(n, h, alpha=ALPHA, beta=BETA, site=0):
    """One heteroatom: alpha_X = alpha + h*beta on a single site."""
    bath = alpha * np.eye(n)
    bath[site, site] = alpha + h * beta
    return bath


def push_pull_bath(n, alpha=ALPHA, beta=BETA):
    """A diagonal that is ODD under the chain reflection about alpha.

    A linear gradient of on-site energies along the chain, the shape a push-pull
    substitution pattern produces.
    It is the counterexample to "the B-condition means a constant diagonal":
    that reading is true only for a DIAGONAL Q like the sublattice flip. Here
    the composite Q = P*K, reflection times sublattice flip, satisfies both
    sub-conditions and the palindrome is exact, with a diagonal that is not
    constant at all. See needs_a_different_q().
    """
    lift = np.linspace(0.5, -0.5, n) * beta
    return np.diag(alpha + lift)


def reflection(n):
    return np.eye(n)[::-1]


def residuals(adj, bath, beta=BETA, involution=None):
    n = len(adj)
    colouring, bipartite = two_colour(adj)
    q = colouring if involution is None else involution
    coupling = beta * adj
    generator = bath + coupling
    # The centre is the spectral mean, i.e. -c = tr(X)/n. For an all-carbon
    # bath this is alpha; for a heteroatom it is not. Where the diagonal is
    # CONSTANT there is an exact route to it that does not divide, and it is
    # taken: tr/n is not bit-exact at every n, see
    # the_centre_is_arithmetic_not_bit_exact().
    diagonal = np.diag(generator)
    centre = diagonal[0] if np.all(diagonal == diagonal[0]) else np.trace(generator) / n
    c = -centre
    q_inv = np.linalg.inv(q)
    a_res = np.linalg.norm(q @ coupling @ q_inv + coupling)
    b_res = np.linalg.norm(bath + q @ bath @ q_inv + 2 * c * np.eye(n))
    spectrum = np.sort(np.linalg.eigvalsh(generator))
    pair_dev = max(abs(spectrum[i] + spectrum[n - 1 - i] - 2 * centre) for i in range(n))
    return bipartite, a_res, b_res, spectrum, pair_dev, centre


def min_a_residual(adj, beta=BETA):
    """min over ALL +-1 sign patterns of ||Q A Q + A||, and the count reached.

    An odd cycle forces an ODD number of frustrated edges, so the minimum is one
    and the residual is 2*sqrt(2)*|beta|, independent of ring size. Any single
    colouring can do worse; this is the statement that does not depend on which
    Q a particular algorithm hands back.
    """
    n = len(adj)
    coupling = beta * adj
    best = None
    for signs in itertools.product((-1.0, 1.0), repeat=n):
        q = np.diag(signs)
        best = min(best, np.linalg.norm(q @ coupling @ q + coupling)) if best is not None \
            else np.linalg.norm(q @ coupling @ q + coupling)
    return best


def ring_spectrum_closed_form(n, alpha=ALPHA, beta=BETA):
    """alpha + 2*beta*cos(2*pi*k/n), the Hueckel ring, eigensolver-free."""
    return np.sort(np.array([alpha + 2 * beta * np.cos(2 * np.pi * k / n) for k in range(n)]))


SYSTEMS = [
    ("benzene C6 ring", ring(6), carbon_bath(6)),
    ("cyclobutadiene C4 ring", ring(4), carbon_bath(4)),
    ("cyclodecapentaene C10", ring(10), carbon_bath(10)),
    ("butadiene C4 chain", chain(4), carbon_bath(4)),
    ("hexatriene C6 chain", chain(6), carbon_bath(6)),
    ("cyclopropenyl C3 ring", ring(3), carbon_bath(3)),
    ("C5 ring", ring(5), carbon_bath(5)),
    ("C7 ring", ring(7), carbon_bath(7)),
    ("pyridine C5N h=0.5", ring(6), hetero_bath(6, H_NITROGEN)),
    ("pyridine-like h=1.0", ring(6), hetero_bath(6, 1.0)),
    ("push-pull 4-chain", chain(4), push_pull_bath(4)),
]

# Which molecules are expected to be palindromic, and why. Written down so that
# the palindrome gate applies to EVERY row rather than only to the ones a
# convenience condition let through: an earlier version gated only the uniform
# baths, which is exactly where a counterexample would have hidden.
PALINDROMIC = {
    "benzene C6 ring", "cyclobutadiene C4 ring", "cyclodecapentaene C10",
    "butadiene C4 chain", "hexatriene C6 chain", "push-pull 4-chain",
}


def main():
    print(f"alpha = {ALPHA} eV, beta = {BETA} eV\n")
    print(f"{'system':24s} {'bipartite':>10s} {'|QAQ+A|':>12s} {'|B+QBQ+2c|':>12s} "
          f"{'pair dev':>12s} {'dev/eps||X||':>13s} {'centre':>10s}")

    for name, adj, bath in SYSTEMS:
        bipartite, a_res, b_res, spectrum, pair_dev, centre = residuals(adj, bath)
        noise = np.finfo(float).eps * abs(spectrum).max()
        print(f"{name:24s} {str(bipartite):>10s} {a_res:12.4e} {b_res:12.4e} "
              f"{pair_dev:12.4e} {pair_dev / noise:13.2e} {centre:10.4f}")

        # (1) The zero diagonal of the adjacency is what makes the centre the
        #     bath's alone. It is a property of the constructors above, not a
        #     measurement, and it is printed rather than gated for that reason;
        #     what a Huckel adjacency cannot have is a self-loop.

        uniform_bath = np.count_nonzero(bath - np.diag(np.diag(bath))) == 0 and \
            np.all(np.diag(bath) == np.diag(bath)[0])

        if uniform_bath:
            # (2) An all-carbon bath is alpha*I, so against the sublattice flip
            #     the B-condition is an IDENTITY, not a measurement. It is
            #     asserted exactly because it cannot break here, which is the
            #     finding rather than evidence for it; b_condition_is_q_free()
            #     says what it rests on, and needs_a_different_q() says what the
            #     B-residual does NOT mean.
            assert b_res == 0.0, name

        if bipartite:
            # (3) Mirror-oddness of the coupling part: an exact rearrangement of
            #     signs, so an exact comparison, not a tolerance.
            assert a_res == 0.0, name
        else:
            # (4) The minimum over ALL sign patterns, not the one BFS returned:
            #     one frustrated edge contributes 2*beta at two symmetric
            #     positions, so ||.||_F = 2*sqrt(2)*|beta| and nothing smaller
            #     is reachable on an odd cycle.
            assert abs(min_a_residual(adj, BETA) - 2 * np.sqrt(2) * abs(BETA)) <= 64 * noise, name
            # The printed A-cell comes from the BFS colouring, the quoted bound
            # from the minimum over all of them. They coincide on a cycle, and
            # that coincidence is asserted rather than assumed.
            assert abs(a_res - min_a_residual(adj, BETA)) <= 64 * noise, name

        # (5) EVERY row is palindrome-gated, against the expectation recorded in
        #     PALINDROMIC. The palindrome has no exact route (an eigensolver),
        #     so it is read against eps*||X||, whose constant is checked across
        #     the beta sweep in scaling_law().
        if name in PALINDROMIC:
            assert pair_dev <= 64 * noise, name
        else:
            assert pair_dev > 64 * noise, name

    odd_ring_mechanism()
    needs_a_different_q()
    framework_side_is_topology_blind()
    b_condition_is_q_free()
    the_centre_is_arithmetic_not_bit_exact()
    scaling_law()
    print("\nAll checks passed.")


def odd_ring_mechanism():
    """The odd-ring break is gated on its closed form, not on "> noise".

    Asserting only that the deviation exceeds machine noise would pass for any
    perturbation whatever. The Hueckel ring has an eigensolver-free spectrum,
    alpha + 2*beta*cos(2*pi*k/n), so the break can be predicted and compared.
    """
    print("\nodd-ring palindrome break against the closed-form prediction:")
    for n in (3, 5, 7, 9, 11):
        _, _, _, spectrum, pair_dev, centre = residuals(ring(n), carbon_bath(n))
        predicted = ring_spectrum_closed_form(n)
        pred_dev = max(abs(predicted[i] + predicted[n - 1 - i] - 2 * ALPHA) for i in range(n))
        noise = np.finfo(float).eps * abs(spectrum).max()
        print(f"  C{n:<3d} measured dev = {pair_dev:.9f}   closed form = {pred_dev:.9f}   "
              f"dev/(eps*||X||) = {pair_dev / noise:.2e}")
        assert abs(pair_dev - pred_dev) <= 64 * noise, n


def needs_a_different_q():
    """The B-condition is Q-RELATIVE, and a constant diagonal is only its
    diagonal-Q reading.

    Against the sublattice flip K, which is diagonal, K B K = B and the
    condition reduces to "the diagonal is constant". That is NOT what the
    sub-condition says. It says the diagonal must be ODD under Q about the
    centre, and a Q that is not diagonal can arrange that for a diagonal that
    varies from site to site.

    The push-pull 4-chain is the case: the on-site energies run in a gradient
    from alpha + 0.5*beta to alpha - 0.5*beta along the chain, so K leaves a
    B-residual of 2*||D|| = 3.578 (D the lift about the centre) and the naive
    reading would call the palindrome broken. It is exact, because Q = P*K
    (reflection times sublattice flip) satisfies BOTH sub-conditions. Pyridine
    is not rescued by any Q at all, for the exact reason given below.

    Naming, so the next reader is not misled: P*K squares to -I here, so it is
    not itself an involution; what is involutive is the conjugation MAP it
    induces, which is all the sub-conditions use.
    """
    print("\nthe B-condition is Q-relative, not 'constant diagonal':")
    n = 4
    adj, bath = chain(n), push_pull_bath(n)
    composite = reflection(n) @ two_colour(adj)[0]
    for label, involution in (("K (sublattice flip)", None), ("P*K (reflection x flip)", composite)):
        _, a_res, b_res, spectrum, pair_dev, _ = residuals(adj, bath, involution=involution)
        noise = np.finfo(float).eps * abs(spectrum).max()
        print(f"  push-pull 4-chain, Q = {label:24s} A-res = {a_res:.3e}  "
              f"B-res = {b_res:.3e}  pair dev = {pair_dev:.3e}")
        if involution is None:
            # Exactly 2*||D||, D the lift about the centre: an exact route, so
            # an exact comparison rather than a threshold.
            lift = np.diag(bath) - np.diag(bath).mean()
            assert abs(b_res - 2 * np.linalg.norm(lift)) <= 64 * noise
        else:
            assert a_res == 0.0 and b_res == 0.0
        # The spectrum does not know which Q was tried: it is palindromic.
        assert pair_dev <= 64 * noise

    # Pyridine is broken for EVERY Q, and that does not need a search. X is
    # diagonalizable, so a Q with Q X Q^-1 = -X - 2c*I exists if and only if the
    # two sides have the same spectrum with multiplicity, i.e.
    # Spec(X) = -2c - Spec(X). Summing that identity over the multiset pins
    # -c = tr(X)/n, which is the centre the code already uses, so there is no
    # freedom left in c either. The pair deviation IS the obstruction, and
    # quoting a best-over-some-family instead would be the weaker claim reached
    # by a route that cannot support the stronger one.
    pyridine_adj, pyridine_bath = ring(6), hetero_bath(6, H_NITROGEN)
    _, _, _, spectrum, pair_dev, _ = residuals(pyridine_adj, pyridine_bath)
    noise = np.finfo(float).eps * abs(spectrum).max()
    print(f"  pyridine, pair dev = {pair_dev:.4e} = {pair_dev / noise:.1e} x the noise model, "
          f"so NO invertible Q rescues it")
    assert pair_dev > 64 * noise

    # And the push-pull chain's rescue is not a lucky Q either: its spectrum IS
    # palindromic, so by the same equivalence some Q had to exist. P*K is an
    # exhibit, not the reason.
    _, _, _, pushpull_spectrum, pushpull_dev, _ = residuals(chain(n), push_pull_bath(n))
    assert pushpull_dev <= 64 * np.finfo(float).eps * abs(pushpull_spectrum).max()


def _pauli_op(single, site, n):
    identity = np.eye(2)
    out = np.array([[1.0]], complex)
    for k in range(n):
        out = np.kron(out, single if k == site else identity)
    return out


def _liouvillian(n, bonds, gammas, coupling=1.0, delta=1.0):
    """Heisenberg XXZ on the given bonds plus per-site Z-dephasing."""
    pauli_x = np.array([[0, 1], [1, 0]], complex)
    pauli_y = np.array([[0, -1j], [1j, 0]])
    pauli_z = np.diag([1.0, -1.0]).astype(complex)
    dim = 2 ** n
    ham = np.zeros((dim, dim), complex)
    for a, b in bonds:
        xx = _pauli_op(pauli_x, a, n) @ _pauli_op(pauli_x, b, n)
        yy = _pauli_op(pauli_y, a, n) @ _pauli_op(pauli_y, b, n)
        zz = _pauli_op(pauli_z, a, n) @ _pauli_op(pauli_z, b, n)
        ham += coupling * (xx + yy) / 4 + delta * zz / 4
    identity = np.eye(dim)
    liouvillian = -1j * (np.kron(ham, identity) - np.kron(identity, ham.T))
    for site, gamma in enumerate(gammas):
        z_site = _pauli_op(pauli_z, site, n)
        liouvillian += gamma * (np.kron(z_site, z_site.T) - np.kron(identity, identity))
    return liouvillian


def _pairing_count(liouvillian, sigma, tolerance):
    """How many eigenvalues find a partner at -2*sigma - lambda: (paired, total).

    A COUNT, in the form MIRROR_SYMMETRY_PROOF uses (64/64, 28/64, 0/64), on the
    complex spectrum rather than the real parts.

    Two earlier versions of this function returned a DISTANCE and both were
    misread. Independent sorts of the two multisets give only an upper bound;
    a minimum-cost assignment gives a max that saturates once the break exceeds
    2*Sigma(gamma), so above that it reports gamma rather than the perturbation.
    A count needs no such reading. Note that it is min-SUM matching, so the
    count is exact where every mode pairs and is an under-count in the partial
    regime, which is the only regime this script does not enter.
    """
    values = np.linalg.eigvals(liouvillian)
    partners = -2 * sigma - values
    cost = np.abs(values[:, None] - partners[None, :])
    rows, columns = linear_sum_assignment(cost)
    return int((cost[rows, columns] <= tolerance).sum()), len(values)


def _worst_pairing_distance(liouvillian, sigma):
    """The matched distance itself, printed so the tolerance is read not trusted."""
    values = np.linalg.eigvals(liouvillian)
    partners = -2 * sigma - values
    cost = np.abs(values[:, None] - partners[None, :])
    rows, columns = linear_sum_assignment(cost)
    return cost[rows, columns].max()


def framework_side_is_topology_blind():
    """The odd cycle has no F1 counterpart, which is the only F1 fact this needs.

    Coulson-Rushbrooke dies on an odd ring. F1 does not notice: the Liouvillian
    palindrome pairs every eigenvalue on the C3 and C5 rings, and it survives an
    arbitrary per-site gamma profile, both of which its registry entry already
    records ("any graph; any N; non-uniform gamma per qubit").

    What this deliberately does NOT test is the other direction, the framework
    counterpart of pyridine's lifted diagonal. That is an on-site field term in
    H, and F138's clause 2 is the law it must satisfy. F138 records clause 2 as
    SPOT-CHECKED, not swept, so it is genuinely open work rather than settled
    elsewhere; what is exhaustively swept there (3^N, zero exceptions) is clause
    1, about dephasing axes. Four successive attempts to state the field side
    here from a handful of C3 rows produced four wrong statements, the last of
    which was an all-or-nothing pairing verdict that is an accident of odd
    rings: at h_z = 0.3 the C4 ring pairs 182 of 256 and the 4-chain 157 of 256,
    not 0. It belongs with F138 and is left there.
    """
    print("\nF1 on the same rings, and under a gamma profile (complex spectrum):")
    rng = np.random.default_rng(20260803)
    triangle = [(0, 1), (1, 2), (2, 0)]
    cases = [
        ("C3 ring, uniform", 3, triangle, np.full(3, 0.05)),
        ("C5 ring, uniform", 5, [(i, (i + 1) % 5) for i in range(5)], np.full(5, 0.05)),
        ("C3 ring, gamma profile", 3, triangle, rng.uniform(0.01, 0.3, 3)),
        ("5-chain, gamma profile", 5, [(i, i + 1) for i in range(4)], rng.uniform(0.01, 0.3, 5)),
    ]
    for name, n, bonds, gammas in cases:
        liouvillian = _liouvillian(n, bonds, gammas)
        # No exact route (a non-normal generator through an eigensolver), so the
        # model is eps*||L|| and the constant is the measured one: the worst
        # matched distance over N = 3..5 x {chain, ring, star, complete} x gamma
        # from 1e-3 to 50 x Delta in {0, 1, 5}, plus random gamma profiles, is
        # 636 * eps*||L||. The gate takes 2048, ~3x that. An earlier version
        # wrote 64 * 4**n, which is ~100x looser again at n = 5 and had no
        # stated model at all.
        scale = np.finfo(float).eps * abs(liouvillian).max()
        paired, total = _pairing_count(liouvillian, gammas.sum(), 2048 * scale)
        worst = _worst_pairing_distance(liouvillian, gammas.sum())
        print(f"  {name:24s} paired {paired:4d} / {total:4d}   "
              f"worst / (eps*||L||) = {worst / scale:6.1f}")
        assert paired == total, name


def b_condition_is_q_free(alpha=ALPHA):
    """Why the all-carbon B-residual is a tautology, stated as one.

    A reader could object that on a non-bipartite ring the sublattice sign flip
    does not exist, so quoting a B-residual there compares against a Q that is
    not available. The objection dissolves, and not by measurement: B = alpha*I
    is a multiple of the identity, so Q B Q^-1 = B for EVERY invertible Q and
    the condition reduces to 2*alpha*I = -2c*I with no reference to the graph,
    the colouring, or Q. The loop below only illustrates it on the class where
    the float arithmetic is exact; the algebra, not the loop, is the argument.
    """
    rng = np.random.default_rng(20260803)
    n = 7
    bath = alpha * np.eye(n)
    c = -alpha
    for _ in range(200):
        perm = rng.permutation(n)
        signs = rng.choice([-1.0, 1.0], size=n)
        q = np.zeros((n, n))
        for i, j in enumerate(perm):
            q[i, j] = signs[i]
        assert np.linalg.norm(bath + q @ bath @ np.linalg.inv(q) + 2 * c * np.eye(n)) == 0.0

    print(f"\nB-condition over 200 signed permutations at N={n}: residual identically 0.0 "
          f"(an identity, illustrated, not measured)")


def the_centre_is_arithmetic_not_bit_exact():
    """alpha = tr(H)/N is exact as arithmetic and NOT bit-exact as floats.

    Read, not gated: the deviation is a function of N through the division and
    of nothing in the chemistry, and it is exactly 0.0 for most N. Quoting
    "alpha = tr(H)/N exactly" without this is an exactness claim on a route
    that carries a rounding channel.
    """
    print("\ntr(H)/N against alpha, all-carbon rings and chains:")
    off = []
    for n in range(3, 13):
        value = np.trace(carbon_bath(n)) / n
        if value != ALPHA:
            off.append(n)
    print(f"  bit-exact for N in 3..12 except N = {off} (deviation "
          f"{np.trace(carbon_bath(7)) / 7 - ALPHA:+.3e}, a float division, not a spectrum)")
    for probe in (-11.2, -11.0):
        exceptions = [n for n in range(3, 13) if np.trace(probe * np.eye(n)) / n != probe]
        print(f"  at alpha = {probe}: {exceptions}")
    print("  the list is a function of alpha alone, nothing chemical. Read, not gated.")


def _mantissa_sweep():
    """The same fit on other grids, to show the slope is grid-dependent."""
    print("  log-log slope on other grids of the same clean system:")
    for mantissa in (2.4, 3.0, 7.0, 1.0):
        betas, ratios = [], []
        for k in range(-3, 10):
            beta = mantissa * 10.0 ** k
            _, _, _, spectrum, dev, _ = residuals(ring(6), carbon_bath(6), beta=beta)
            norm = abs(spectrum).max()
            if norm > 10 * abs(BETA):
                betas.append(beta)
                ratios.append(dev / (np.finfo(float).eps * norm))
        fit = np.polyfit(np.log10(betas), np.log10(ratios), 1)[0]
        print(f"    mantissa {mantissa:4.1f} x 10^k, {len(betas)} points: slope = {fit:+.4f}")


def scaling_law():
    """The bipartite pair deviation is machine noise on ||X||, not a window.

    An eigensolver has no exact route, so the tolerance in main() must be a law
    rather than a number. The law is that dev/(eps*||X||) is CONSTANT, and the
    gate is on the SPREAD of that ratio, which a deviation growing faster than
    linearly in beta would blow up where a fixed ceiling would absorb it for
    several decades first.

    The lever is SIX decades, not the twelve the beta values span: below
    beta ~ 1e3 the norm is alpha-dominated (max|spec| = 11.41 at beta = 2.4e-3),
    so ||X|| stops tracking beta and those rows test nothing about the scaling.
    They are printed because saying which rows are inert is part of the law.
    """
    print("\nbipartite pair deviation against the eps*||X|| error model:")
    scaling_ratios = []
    scaling_betas = []
    for beta in (2.4e-3, 2.4, 2.4e3, 2.4e4, 2.4e5, 2.4e6, 2.4e7, 2.4e8, 2.4e9):
        _, _, _, spectrum, pair_dev, _ = residuals(ring(6), carbon_bath(6), beta=beta)
        norm = abs(spectrum).max()
        ratio = pair_dev / (np.finfo(float).eps * norm)
        tracks_beta = norm > 10 * abs(BETA)
        if tracks_beta:
            scaling_ratios.append(ratio)
            scaling_betas.append(beta)
        print(f"  beta = {beta:10.1e}   ||X|| = {norm:12.4e}   dev = {pair_dev:10.3e}   "
              f"dev/(eps*||X||) = {ratio:7.2f}   {'scaling' if tracks_beta else 'alpha-dominated'}")

    # THE GATE IS THE RATIO ITSELF, per point. A symmetric eigensolver's backward
    # error bounds |lambda~ - lambda| by a modest polynomial in n times eps*||X||,
    # so the ratio is O(1) by construction; the constant is not derived from that
    # bound, it is measured. Worst over a mantissa sweep and over other systems
    # (C6 ring 3.71 on the gated grid, 5.56 over 60 random mantissas, C6 chain
    # 5.61, C40 chain 6.15), so the ceiling 8 is about 1.3x the worst seen and
    # roughly 2x the gated case. It is a measured ceiling, honestly a thin one.
    #
    # The slope and spread below are DIAGNOSTICS, printed and not gated. An
    # earlier version asserted spread <= 4 on three points, which more points
    # failed at 5.45; its replacement asserted a log-log slope <= 0.05, and the
    # mantissa sweep printed underneath shows why that is not a law either: the
    # slope is not even stable in SIGN across grids of the same clean system.
    exponents = np.log10(np.array(scaling_betas))
    slope = np.polyfit(exponents, np.log10(np.array(scaling_ratios)), 1)[0]
    _mantissa_sweep()
    spread = max(scaling_ratios) / min(scaling_ratios)
    print(f"  over the six decades where ||X|| tracks beta: max ratio = "
          f"{max(scaling_ratios):.2f}, spread = {spread:.2f}, log-log slope = {slope:+.4f} "
          f"(both diagnostics, grid-dependent; the gate is the ratio)")
    assert max(scaling_ratios) <= 8.0, max(scaling_ratios)


if __name__ == "__main__":
    main()
