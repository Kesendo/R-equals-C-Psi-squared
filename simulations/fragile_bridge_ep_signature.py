"""The fragile bridge's threshold: where two levels of opposite Krein signature meet.

Two chains of `n_per_chain` qubits, Heisenberg within each, joined by a Heisenberg
bridge; chain A dephases at +gamma, chain B at -gamma, so Sigma-gamma is exactly 0
and the palindrome is centred on zero. Everything below runs at two qubits per chain.
Pi forces inversion lambda <-> -lambda; Hermiticity preservation separately supplies
conjugate pairing. Together they give the generic off-axis quartet
{lambda, lambda*, -lambda, -lambda*}.

Symbols. P_site is the site reversal (F71's mirror, site k <-> N-1-k), G the map
rho -> P_site rho P_site. In the exact search of section 3, x = J_bridge / J with J = 1. R is the one-sided flip rho -> rho X^N, the edge mirror of
F118's mirror group (Pi = R o D there, D the transpose). F is the two-sided flip
rho -> X^N rho X^N, which is Pi^2 (F118's centre).

L conserves the popcounts (p, q) of ket and bra, so it splits into 25 blocks. What holds
an eigenvalue on the imaginary axis inside a block is GK, G composed with complex
conjugation: it maps every block onto itself and sends lambda to -lambda*, so a simple
eigenvalue of a block on the axis is its own partner and cannot leave alone. Leaving
takes two levels meeting, and by Krein's rule two levels of opposite signature under G.
Pi composed with Hermitian conjugation is R composed with complex conjugation
(F118 and F119); it also sends lambda to -lambda* and maps block (p, q) onto (p, N - q).
The product of the two antilinear maps is the linear S = G o R, which commutes with L at
Sigma-gamma = 0 and maps (p, q) onto (p, N - q) with the same spectrum. F maps (p, q)
onto (N - p, N - q) with the same spectrum, Hermitian conjugation onto (q, p) with the
conjugate one. The middle block (2,2) is its own image under all of them; F splits it
into an even and an odd half.

Section 1 is the threshold certificate at J_bridge = 1.0 and 1.9, where the first sector
is a (2,2) half. Inside it the squared gap f = (lambda_a - lambda_b)^2 of the colliding
pair is real for real gamma: -(difference of the two imaginary parts)^2 below the
collision, (2 Re lambda)^2 above it. So f changes sign there and its zero is real;
brentq on that sign change locates gamma*. At gamma* the producer reads:
  - the pair gap (of the order of the square root of machine precision at a located EP2);
  - the smallest singular values of L - lambda* I, inside the sector and in the full L.
    One at rounding level and the next far above it is geometric multiplicity 1 at
    algebraic multiplicity 2: the pair is a 2x2 Jordan block. L is block-diagonal, so
    the singular values of L - lambda* I are those of its blocks together; the output
    names the block each of the three smallest comes from;
  - the pair at gamma*(1 -/+ delta) for delta = 1e-4 and 1e-6: two distinct eigenvalues
    on the imaginary axis below, a mirror pair +/-Re at one Im above, and the same
    half-split / sqrt|delta| on both sides (one Puiseux branch, lambda - lambda* =
    +/- c sqrt(gamma/gamma* - 1)), plus that coefficient above gamma* from delta = 1e-2
    down to 1e-6, where the next Puiseux order shows at the largest delta;
  - the Krein signatures v^+ G v of the two approaching eigenvectors below gamma*.
    G commutes with the Hamiltonian part of L and anticommutes with the gain-loss
    dissipator, and L is complex symmetric here, so i*L is G-pseudo-Hermitian at every
    real gamma. Opposite signatures that shrink to zero at gamma* are what a Krein
    collision is.
Where the bisected gamma_crit sits below the collision, section 1 says why: the bisection
returns a transition point of its own sequence of max Re lambda readings against 1e-12,
and near gamma* the approaching pair is ill-conditioned, so rounding noise in its real
parts can cross that level while the pair is still split on the imaginary axis.

Section 2 is the sector structure: the symmetry checks, then at each sampled coupling the
first block to go unstable, the number of blocks its orbit under F and S holds (the copies
of the collision in the full L), its collision gamma* and frequency, and what the full L
carries at lambda*. It locates the window in which the (0,1) orbit goes first, follows
the two (2,2) halves near the sampled maximum and locates where they cross.

Section 3 is the zeros of the threshold. At gamma = 0 each block of L_H splits into
G-even and G-odd levels (and by F and S where they act inside it). Where an even and an
odd level of one such sector have the same frequency (omega, an eigenvalue of i*L_H),
the gain-loss part, odd under G, couples them at first order:
lambda = -i*omega +/- s*gamma + O(gamma^2), so max Re lambda = s*gamma + O(gamma^2) > 0
for small gamma > 0 and gamma_crit = 0, with a semisimple pair at gamma = 0 and no
exceptional point. The couplings are found exactly (the H spectrum per popcount sector
in closed form, every pair of opposite-G levels of one sector, the radicals eliminated,
the resultant factored over Q) and again numerically (level crossings on a grid);
s is read by first-order degenerate perturbation theory and measured on the full L,
whose max Re lambda is also checked positive over a range of gammas. The section then
shows the threshold around each zero and where the finer sweep's low points sit.

Section 4 keeps the finite-offset diagnostics above the bisected gamma_crit:

  1. Re lambda ~ sqrt(delta), delta = gamma/gamma_crit - 1, the coefficient nearly
     constant across the sampled three decades (four samples).
  2. The finite-offset gap to the across-axis partner -lambda* equals
     2*abs(Re lambda), a derived symmetry identity rather than independent evidence.
  3. The finite-offset Petermann reading grows approximately by a decade per
     sampled decade of delta, the K ~ 1/delta of an EP2.

Measured from the bisected gamma_crit, the smallest-delta rows carry the bisection offset;
section 1 measures the square-root law from gamma* itself.

Section 5 reads the max-real-part selected mode's K across a wide gamma range,
omitting degenerate selections for which a single-vector K is basis-dependent.

Section 6 sweeps the coupling, J_bridge = 0.005 to 20 in steps of 0.005, reading the
first block to go unstable and its threshold at every point: the runs of first blocks,
their copies, the largest threshold on the grid, whether each run rises or falls, and at
every grid coupling the Jordan rank of the collision in the symmetry sector that goes first.

Section 7 prints gamma_crit * J_bridge at six couplings from 20 to 1000, as samples.

Anchors: hypotheses/FRAGILE_BRIDGE.md Sections 2 and 3, experiments/PT_SYMMETRY_ANALYSIS.md,
docs/ANALYTICAL_FORMULAS.md F19 and F40, F118 (PROOF_PI_FACTORS_AS_R_TIMES_D.md).
"""
import ast
import itertools
import os
import sys
from itertools import product

import numpy as np
import scipy.sparse as sp
import sympy
from scipy.linalg import eig, eigvals, svdvals
from scipy.optimize import brentq

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "fragile_bridge_ep_signature.txt")

_LINES = []


def log(msg=""):
    print(msg, flush=True)
    _LINES.append(msg)


I2 = sp.eye(2, format="csr", dtype=complex)
SX = sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
SY = sp.csr_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex))
SZ = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))


def op_at(op, qubit, n_qubits):
    result = sp.eye(1, format="csr", dtype=complex)
    for k in range(n_qubits):
        result = sp.kron(result, op if k == qubit else I2, format="csr")
    return result


def build_L(n_per_chain, gamma, J=1.0, J_bridge=0.1):
    """The gain-loss Liouvillian: +gamma on chain A, -gamma on chain B."""
    n_total = 2 * n_per_chain
    d = 2 ** n_total
    d2 = d * d
    H = sp.csr_matrix((d, d), dtype=complex)
    for i in range(n_per_chain - 1):
        for P in (SX, SY, SZ):
            H = H + J * (op_at(P, i, n_total) @ op_at(P, i + 1, n_total))
    for i in range(n_per_chain, 2 * n_per_chain - 1):
        for P in (SX, SY, SZ):
            H = H + J * (op_at(P, i, n_total) @ op_at(P, i + 1, n_total))
    for P in (SX, SY, SZ):
        H = H + J_bridge * (op_at(P, n_per_chain - 1, n_total)
                            @ op_at(P, n_per_chain, n_total))
    Id = sp.eye(d, format="csr", dtype=complex)
    L = -1j * (sp.kron(H, Id, format="csr") - sp.kron(Id, H.T, format="csr"))
    Id2 = sp.eye(d2, format="csr", dtype=complex)
    for k in range(n_per_chain):
        Zk = op_at(SZ, k, n_total)
        L = L + gamma * (sp.kron(Zk, Zk.conj(), format="csr") - Id2)
    for k in range(n_per_chain, 2 * n_per_chain):
        Zk = op_at(SZ, k, n_total)
        L = L + (-gamma) * (sp.kron(Zk, Zk.conj(), format="csr") - Id2)
    return L


def max_re(n_per_chain, gamma, J_bridge):
    return float(np.max(eigvals(build_L(n_per_chain, gamma, J_bridge=J_bridge).toarray()).real))


def gamma_crit(n_per_chain, J_bridge, lo=0.0, hi=2.0, iterations=60):
    """Bisect on max Re lambda leaving zero. The axis below is clean to ~1e-14."""
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if max_re(n_per_chain, mid, J_bridge) > 1e-12:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def leading_mode(n_per_chain, gamma, J_bridge):
    """Return (lambda, gap to the nearest other eigenvalue, Petermann K)."""
    A = build_L(n_per_chain, gamma, J_bridge=J_bridge).toarray()
    w, vr = eig(A)
    wl, vl = eig(A.conj().T)
    i = int(np.argmax(w.real))
    gap = float(np.sort(np.abs(w - w[i]))[1])
    j = int(np.argmin(np.abs(np.conj(wl) - w[i])))
    r = vr[:, i] / np.linalg.norm(vr[:, i])
    l = vl[:, j] / np.linalg.norm(vl[:, j])
    overlap = abs(np.vdot(l, r)) ** 2
    return w[i], gap, (float("inf") if overlap == 0.0 else 1.0 / overlap)


def generator_parts(n_per_chain, J_bridge):
    """L(gamma) = LH + gamma * D. H is real symmetric, so the Hamiltonian part LH = L(0) is
    purely imaginary and D = Re L(1) is the real diagonal gain-loss dissipator: the split
    is exact."""
    LH = build_L(n_per_chain, 0.0, J_bridge=J_bridge).toarray()
    D = build_L(n_per_chain, 1.0, J_bridge=J_bridge).toarray().real
    return LH, D


def site_reversal(n_total):
    """P_site, the site reversal of the 2n-qubit chain on the computational basis (F71's
    mirror). It maps chain A onto chain B and the bridge onto itself."""
    d = 2 ** n_total
    R = np.zeros((d, d))
    for b in range(d):
        R[int(format(b, "0{0}b".format(n_total))[::-1], 2), b] = 1.0
    return R


def chain_reflection_metric(n_per_chain):
    """G = P_site (x) P_site on the row-major vec basis: rho -> P_site rho P_site."""
    P = site_reversal(2 * n_per_chain)
    return np.kron(P, P)


def symmetry_superoperators(n_total):
    """Superoperators on the row-major vec basis (a*d + b holds |a><b|): G (rho -> P_site rho
    P_site), F (rho -> X^N rho X^N), R (rho -> rho X^N, F118's edge mirror), S = G R, and SWAP,
    the index map of Hermitian conjugation (rho^dagger = SWAP applied to the complex conjugate)."""
    d = 2 ** n_total
    P = site_reversal(n_total)
    XN = np.zeros((d, d))
    for a in range(d):
        XN[d - 1 - a, a] = 1.0
    Id = np.eye(d)
    G = np.kron(P, P)
    F = np.kron(XN, XN)
    T = np.kron(Id, XN)
    SWAP = np.zeros((d * d, d * d))
    for a in range(d):
        for b in range(d):
            SWAP[b * d + a, a * d + b] = 1.0
    return G, F, T, G @ T, SWAP


def palindrome_superoperator(n_total):
    """Pi in the repo convention, per site I <-> X with phase 1, Y -> iZ and Z -> iY, as a
    superoperator on the row-major vec basis (Pauli strings are orthogonal with norm^2 = d)."""
    d = 2 ** n_total
    paulis = {"I": np.eye(2, dtype=complex), "X": SX.toarray(), "Y": SY.toarray(), "Z": SZ.toarray()}
    image = {"I": ("X", 1.0), "X": ("I", 1.0), "Y": ("Z", 1j), "Z": ("Y", 1j)}
    Pi = np.zeros((d * d, d * d), dtype=complex)
    for labels in itertools.product("IXYZ", repeat=n_total):
        source = np.eye(1, dtype=complex)
        target = np.eye(1, dtype=complex)
        phase = 1.0 + 0.0j
        for label in labels:
            source = np.kron(source, paulis[label])
            target = np.kron(target, paulis[image[label][0]])
            phase *= image[label][1]
        Pi += np.outer(phase * target.reshape(-1), source.reshape(-1).conj()) / d
    return Pi


def popcount_sector_frame(n_per_chain):
    """Index sets of the (p, q) popcount blocks of the row-major vec basis (a*d + b holds
    |a><b|, p = popcount a, q = popcount b), and orthonormal bases of the spin-flip even
    and odd halves of the self-conjugate middle block."""
    n_total = 2 * n_per_chain
    d = 2 ** n_total
    ones = [bin(a).count("1") for a in range(d)]
    blocks = {(p, q): np.array([a * d + b for a in range(d) for b in range(d)
                                if ones[a] == p and ones[b] == q])
              for p, q in product(range(n_total + 1), repeat=2)}
    middle = blocks[(n_per_chain, n_per_chain)]
    position = {index: i for i, index in enumerate(middle)}
    even, odd, seen = [], [], set()
    for index in middle:
        if index in seen:
            continue
        a, b = divmod(index, d)
        partner = (d - 1 - a) * d + (d - 1 - b)
        seen.update((index, partner))
        e = np.zeros(len(middle))
        o = np.zeros(len(middle))
        e[position[index]], e[position[partner]] = 1.0, 1.0
        o[position[index]], o[position[partner]] = 1.0, -1.0
        even.append(e / np.sqrt(2.0))
        odd.append(o / np.sqrt(2.0))
    return blocks, np.array(even).T, np.array(odd).T


def sectors(LH, D, frame, n_per_chain):
    """name -> (LH restricted, D restricted) for every block, the middle one split by parity."""
    blocks, v_even, v_odd = frame
    middle = (n_per_chain, n_per_chain)
    out = {}
    for (p, q), index in blocks.items():
        if (p, q) != middle:
            out["({0},{1})".format(p, q)] = (LH[np.ix_(index, index)], D[np.ix_(index, index)])
    index = blocks[middle]
    LHm, Dm = LH[np.ix_(index, index)], D[np.ix_(index, index)]
    for parity, V in (("flip-even", v_even), ("flip-odd", v_odd)):
        out["({0},{0}) {1}".format(n_per_chain, parity)] = (V.T @ LHm @ V, V.T @ Dm @ V)
    return out


def block_orbit(p, q, n_total):
    """The blocks reached from (p, q) by F, S and Hermitian conjugation."""
    same = {(p, q), (n_total - p, n_total - q), (p, n_total - q), (n_total - p, q)}
    return sorted(same | {(b, a) for a, b in same})


def copies_of(name, n_total):
    """How many blocks carry the same spectrum as this one: its orbit under F and S."""
    if "flip" in name:
        return 1
    p, q = ast.literal_eval(name)
    return len({(p, q), (n_total - p, n_total - q), (p, n_total - q), (n_total - p, q)})


def canonical(name, n_total):
    """The smallest member of the block's orbit, so that equal thresholds print one name."""
    if "flip" in name:
        return name
    p, q = ast.literal_eval(name)
    a, b = block_orbit(p, q, n_total)[0]
    return "({0},{1})".format(a, b)


def orbit_label(name, n_total):
    if "flip" in name:
        return name
    p, q = ast.literal_eval(name)
    return "/".join("({0},{1})".format(a, b) for a, b in block_orbit(p, q, n_total))


def sector_max_re(LHs, Ds, gamma):
    return float(np.max(eigvals(LHs + gamma * Ds).real))


def sector_threshold(LHs, Ds, tol=1e-9, gamma_max=2.0, points=300):
    """First gamma at which the sector leaves the axis: a geometric scan to the first
    unstable grid point, then bisection inside that bracket."""
    previous = 0.0
    for gamma in np.geomspace(1e-3, gamma_max, points):
        if sector_max_re(LHs, Ds, gamma) > tol:
            lo, hi = previous, gamma
            for _ in range(55):
                mid = 0.5 * (lo + hi)
                if sector_max_re(LHs, Ds, mid) > tol:
                    hi = mid
                else:
                    lo = mid
            return 0.5 * (lo + hi)
        previous = gamma
    return np.inf


def squared_gap(LHs, Ds, gamma, anchor):
    """f = (lambda_a - lambda_b)^2 of the two sector eigenvalues nearest the anchor, with the
    pair, the midpoint and the sector spectrum."""
    w = eigvals(LHs + gamma * Ds)
    order = np.argsort(np.abs(w - anchor))
    a, b = w[order[0]], w[order[1]]
    return (a - b) ** 2, a, b, w


def locate_in_sector(LHs, Ds, gamma_threshold):
    """The collision inside one sector: brentq on the real squared gap of the pair that
    leaves the axis, bracketed where it is negative (below) and positive (above)."""
    w = eigvals(LHs + gamma_threshold * (1.0 + 1e-4) * Ds)
    anchor = complex(0.0, w[int(np.argmax(w.real))].imag)

    def f(gamma):
        return squared_gap(LHs, Ds, gamma, anchor)[0].real

    lo, hi = gamma_threshold * (1.0 - 1e-6), gamma_threshold * (1.0 + 1e-4)
    k = 0
    while f(lo) >= 0.0 and k < 6:
        k += 1
        lo = gamma_threshold * max(0.0, 1.0 - 1e-6 * 10 ** k)
    gamma_star = brentq(f, lo, hi, xtol=1e-17, rtol=8.9e-16, maxiter=200)
    _, a, b, spectrum = squared_gap(LHs, Ds, gamma_star, anchor)
    return gamma_star, 0.5 * (a + b), abs(a - b), anchor, spectrum


def all_thresholds(LH, D, frame, n_per_chain):
    parts = sectors(LH, D, frame, n_per_chain)
    return parts, {name: sector_threshold(*mats) for name, mats in parts.items()}


def first_of(thresholds):
    first_value = min(thresholds.values())
    first_names = [name for name, value in thresholds.items()
                   if abs(value - first_value) <= 1e-9 * max(first_value, 1e-300)]
    return min(first_names), first_value


def sector_survey(LH, D, frame, n_per_chain):
    """Every sector's threshold, the first sector(s), and the collision located in it."""
    parts, thresholds = all_thresholds(LH, D, frame, n_per_chain)
    first_name, first_value = first_of(thresholds)
    gamma_star, lam_star, gap, anchor, spectrum = locate_in_sector(*parts[first_name], first_value)
    return parts, thresholds, first_name, gamma_star, lam_star, gap, anchor, spectrum


def full_l_at_collision(LH, D, gamma_star, lam_star):
    """What the full L carries at lambda*: eigenvalues within 1e-5, the singular values of
    L - lambda* I (sorted), and the eigenvalue count with Re > 1e-9 at gamma*(1 + 1e-6)."""
    L_star = LH + gamma_star * D
    spectrum = eigvals(L_star)
    within = int(np.sum(np.abs(spectrum - lam_star) < 1e-5))
    sv = np.sort(svdvals(L_star - lam_star * np.eye(L_star.shape[0])))
    above = eigvals(LH + gamma_star * (1.0 + 1e-6) * D)
    return within, sv, int(np.sum(above.real > 1e-9))


def format_eig(z):
    return "({0:+.1e} {1:+.9f}i)".format(z.real, z.imag)


def fmt_poly(expr):
    return str(expr).replace("**", "^").replace("*", "")


# ---------------------------------------------------------------------------------------
# The H spectrum per popcount sector in closed form, and the exact degeneracy search
# ---------------------------------------------------------------------------------------

JB = sympy.Symbol("x", positive=True)
XS = sympy.Symbol("E")


def sector_states(n_total, n_up):
    return [a for a in range(2 ** n_total) if bin(a).count("1") == n_up]


def bits_of(a, n_total):
    return [(a >> (n_total - 1 - s)) & 1 for s in range(n_total)]


def sector_hamiltonian(n_total, n_up, bonds, numeric=None):
    """H on one popcount sector, symbolic in J_bridge (or numeric if a value is given)."""
    states = sector_states(n_total, n_up)
    pos = {a: i for i, a in enumerate(states)}
    M = sympy.zeros(len(states)) if numeric is None else np.zeros((len(states), len(states)))
    for a in states:
        b = bits_of(a, n_total)
        for i, j, coupling in bonds:
            c = coupling if numeric is None else (numeric if coupling == JB else float(coupling))
            M[pos[a], pos[a]] += c * (1 if b[i] == b[j] else -1)
            if b[i] != b[j]:
                b2 = list(b)
                b2[i], b2[j] = b2[j], b2[i]
                a2 = sum(v << (n_total - 1 - s) for s, v in enumerate(b2))
                M[pos[a2], pos[a]] += 2 * c
    return states, M


def symbolic_levels(n_per_chain, j_generic=0.8123):
    """Per popcount sector: [(E(J) closed form, R parity, X parity or 0)], with the closed forms
    from the factored characteristic polynomial and the parities read from numeric eigenpairs
    at a generic coupling."""
    n_total = 2 * n_per_chain
    bonds = [(i, i + 1, sympy.Integer(1)) for i in range(n_per_chain - 1)]
    bonds += [(n_per_chain - 1, n_per_chain, JB)]
    bonds += [(i, i + 1, sympy.Integer(1)) for i in range(n_per_chain, n_total - 1)]
    out = {}
    for n_up in range(n_total + 1):
        states, M = sector_hamiltonian(n_total, n_up, bonds)
        closed = []
        for factor, mult in sympy.factor_list(M.charpoly(XS).as_expr())[1]:
            if sympy.degree(factor, XS) > 2:
                raise NotImplementedError("a sector level needs more than a quadratic")
            for root, m in sympy.roots(sympy.Poly(factor, XS)).items():
                closed.extend([root] * (m * mult))
        _, Mn = sector_hamiltonian(n_total, n_up, bonds, numeric=j_generic)
        w, V = np.linalg.eigh(Mn)
        if len(w) > 1 and np.min(np.diff(w)) < 1e-6:
            raise RuntimeError("the generic coupling is degenerate in sector {0}".format(n_up))
        pos = {a: i for i, a in enumerate(states)}
        R = np.zeros_like(Mn)
        X = np.zeros_like(Mn)
        for a in states:
            R[pos[int(format(a, "0{0}b".format(n_total))[::-1], 2)], pos[a]] = 1.0
            if 2 * n_up == n_total:
                X[pos[(2 ** n_total - 1) ^ a], pos[a]] = 1.0
        levels, unused = [], list(closed)
        for k in range(len(w)):
            v = V[:, k]
            values = [abs(float(e.subs(JB, j_generic)) - w[k]) for e in unused]
            match = unused.pop(int(np.argmin(values)))
            if min(values) > 1e-9:
                raise RuntimeError("no closed form matches a sector eigenvalue")
            levels.append((match, int(round(v @ R @ v)),
                           int(round(v @ X @ v)) if 2 * n_up == n_total else 0))
        out[n_up] = levels
    return out


def block_levels(levels, p, q, n_total):
    """(omega = E_a - E_b, G parity, the parities of F, S or FS where they act inside the block)."""
    half = n_total // 2
    out = []
    for (ea, ra, xa), (eb, rb, xb) in itertools.product(levels[p], levels[q]):
        label = []
        if (p, q) == (half, half):
            label.append(("F", xa * xb))
        if q == half:
            label.append(("S", ra * rb * xb))
        if p == half and q != half:
            label.append(("FS", ra * rb * xa))
        out.append((sympy.expand(ea - eb), ra * rb, tuple(label)))
    return out


def eliminate_radicals(expr):
    """A polynomial in J that vanishes wherever expr does, the square roots eliminated by
    resultants against their defining squares."""
    radicals = sorted({atom for atom in expr.atoms(sympy.Pow) if atom.exp == sympy.Rational(1, 2)},
                      key=str)
    symbols = sympy.symbols("r0:{0}".format(len(radicals))) if radicals else ()
    P = sympy.expand(expr.subs({rad: sym for rad, sym in zip(radicals, symbols)}))
    for rad, sym in zip(radicals, symbols):
        if P.has(sym):
            P = sympy.resultant(P, sym ** 2 - rad.base, sym)
    return sympy.expand(P)


def exact_degeneracy_couplings(levels, n_total, representatives):
    """Every x > 0 where a G-even and a G-odd level of one symmetry sector of a block coincide.
    Returns {J0 float: [(block, sector label, minimal polynomial, exact root)]}, the smallest
    level separation left by a spurious (squaring) root, and whether an opposite-G pair is
    identically equal (which would make the threshold zero at every coupling)."""
    found = {}
    spurious = np.inf
    identical = False
    for blk in representatives:
        lv = block_levels(levels, blk[0], blk[1], n_total)
        for (w1, g1, l1), (w2, g2, l2) in itertools.combinations(lv, 2):
            if g1 == g2 or l1 != l2:
                continue
            if sympy.simplify(w1 - w2) == 0:
                identical = True
                continue
            P = eliminate_radicals(w1 - w2)
            if P == 0:
                identical = True
                continue
            for factor, _ in sympy.factor_list(P)[1]:
                if not factor.has(JB):
                    continue
                for root in sympy.Poly(factor, JB).real_roots():
                    value = float(root.evalf(40))
                    if value <= 0.0:
                        continue
                    residual = abs(complex((w1 - w2).subs(JB, root.evalf(40)).evalf(40)))
                    if residual < 1e-25:
                        found.setdefault(round(value, 12), []).append((blk, l1, factor, root))
                    else:
                        spurious = min(spurious, residual)
    return found, spurious, identical


# numeric sector bases for the first-order slopes and the grid scan

def block_sector_bases(blk, frame, ops, n_total):
    """Orthonormal bases (block-local) of the joint eigenspaces of the involutions acting
    inside the block, each split by G: {sector label: {+1: basis, -1: basis}}."""
    blocks = frame[0]
    index = blocks[blk]
    G, F, S, FS = ops
    half = n_total // 2
    acting = []
    if blk == (half, half):
        acting.append(F)
    if blk[1] == half:
        acting.append(S)
    if blk[0] == half and blk[1] != half:
        acting.append(FS)
    mats = [M[np.ix_(index, index)] for M in acting] + [G[np.ix_(index, index)]]
    rng = np.random.default_rng(11)
    combo = sum(rng.uniform(1.0, 2.0) * 2.0 ** k * M for k, M in enumerate(mats))
    _, V = np.linalg.eigh(combo)
    out = {}
    for k in range(V.shape[1]):
        v = V[:, k]
        label = tuple(int(round(v @ M @ v)) for M in mats)
        out.setdefault(label[:-1], {}).setdefault(label[-1], []).append(k)
    return {key: {g: V[:, cols] for g, cols in parts_.items()} for key, parts_ in out.items()}


def sector_frequencies(LHb, basis):
    return np.sort(np.linalg.eigvalsh(1j * (basis.conj().T @ LHb @ basis)))


def first_order_slope(LH, D, index, basis_pair, omega0, tol=1e-7):
    """s: the largest eigenvalue of D on the degenerate G-even + G-odd subspace at omega0."""
    LHb, Db = LH[np.ix_(index, index)], D[np.ix_(index, index)]
    vectors = []
    for g, B in basis_pair.items():
        w, U = np.linalg.eigh(1j * (B.conj().T @ LHb @ B))
        for k in range(len(w)):
            if abs(w[k] - omega0) < tol:
                vectors.append(B @ U[:, k])
    V = np.array(vectors).T
    M = V.conj().T @ Db @ V
    return float(np.max(np.linalg.eigvalsh(0.5 * (M + M.conj().T)))), len(vectors)


NPC = 2
N_TOTAL = 2 * NPC
BRIDGES = (1.0, 1.9)
DELTAS = (1e-2, 1e-3, 1e-4, 1e-5)
SIMPLE_MODE_GAP_TOL = 1e-10
CERTIFICATE_DELTAS = (1e-4, 1e-6)
SQRT_LAW_DELTAS = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
SECTOR_BRIDGES = (0.1, 0.5, 1.0, 1.2, 1.4, 1.46, 1.5, 1.9, 1.95, 1.96, 2.0, 5.0, 10.0)
CROSSING_BRIDGES = (1.8, 1.85, 1.9, 1.94, 1.95, 1.953, 1.96, 2.0)
REPRESENTATIVES = ((0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
FINER_GRID_POINTS = (0.7, 0.8, 1.1, 1.3, 3.2, 3.6, 4.0, 4.4)
SLOPE_GAMMAS = (1e-2, 1e-4, 1e-6)
V_OFFSETS = (-0.01, -0.005, 0.005, 0.01)

log("=" * 78)
log("The fragile bridge threshold: where two levels of opposite Krein signature meet")
log("=" * 78)
log()
log("Two chains of 2 qubits, Heisenberg, gain-loss dephasing +/- gamma, Sigma-gamma = 0.")
log()

FRAME = popcount_sector_frame(NPC)
G_CRIT = {J_BRIDGE: gamma_crit(NPC, J_BRIDGE) for J_BRIDGE in BRIDGES}
G_METRIC = chain_reflection_metric(NPC)

log("-" * 78)
log("1. Threshold certificate: the colliding pair, located and read at the collision")
log("-" * 78)
log()
log("The collision is located inside the first sector to go unstable by brentq on the pair's")
log("squared gap f = (lambda_a - lambda_b)^2, which is real for real gamma: negative while the")
log("two sit on the imaginary axis, positive once they are a mirror pair. The rank of")
log("L - lambda* I is read from its singular values. The Krein metric is G: rho -> P_site rho P_site,")
log("P_site the site reversal (F71's mirror) that swaps the gain chain and the loss chain.")
log()
GAMMA_STAR = {}
CERTIFIED = []
BELOW = {}
for J_BRIDGE in BRIDGES:
    GC = G_CRIT[J_BRIDGE]
    LH, D = generator_parts(NPC, J_BRIDGE)
    PARTS, _, NAME, g_star, lam_star, pair_gap, anchor, sector_spectrum = sector_survey(LH, D, FRAME, NPC)
    GAMMA_STAR[J_BRIDGE] = g_star
    LHs, Ds = PARTS[NAME]
    f_below = squared_gap(LHs, Ds, g_star * (1.0 - 1e-6), anchor)[0].real
    f_above = squared_gap(LHs, Ds, g_star * (1.0 + 1e-4), anchor)[0].real
    sector_third = np.sort(np.abs(sector_spectrum - lam_star))[2]
    sector_sv = np.sort(svdvals(LHs + g_star * Ds - lam_star * np.eye(LHs.shape[0])))
    within, sv, _ = full_l_at_collision(LH, D, g_star, lam_star)
    ATTRIBUTED = sorted((float(value), name) for name, (lhs, ds) in PARTS.items()
                        for value in svdvals(lhs + g_star * ds - lam_star * np.eye(lhs.shape[0])))
    NEAREST = {name: float(np.min(np.abs(eigvals(PARTS[name][0] + g_star * PARTS[name][1]) - lam_star)))
               for _, name in ATTRIBUTED[1:3]}
    log("J_bridge = {0}   first sector to go unstable: {1} (dimension {2})".format(
        J_BRIDGE, NAME, LHs.shape[0]))
    log("  collision gamma* = {0:.12f}; f(gamma*(1 - 1e-6)) = {1:+.1e} and f(gamma*(1 + 1e-4)) = {2:+.1e},"
        " so the zero of f is real".format(g_star, f_below, f_above))
    log("  lambda* = {0}   pair gap {1:.1e}   next eigenvalue in the sector at {2:.3f}".format(
        format_eig(lam_star), pair_gap, sector_third))
    log("  in the sector: smallest singular values of L - lambda* I: {0:.1e}, then {1:.4f}".format(
        sector_sv[0], sector_sv[1]))
    log("  in the full L: eigenvalues within 1e-5 of lambda*: {0}; smallest singular values"
        " {1:.1e}  {2:.1e}  {3:.1e}".format(within, sv[0], sv[1], sv[2]))
    log("  by block (L is block-diagonal, so its singular values are the blocks' together, here to"
        " {0:.0e}): {1}".format(max(abs(a - b) for (a, _), b in zip(ATTRIBUTED[:3], sv[:3])),
                                 ", ".join("{0:.1e} in {1}".format(v, n) for v, n in ATTRIBUTED[:3])))
    log("  the eigenvalues nearest lambda* in the blocks of the second and third: " + ", ".join(
        "{0} at {1:.2e}".format(n, NEAREST[n]) for _, n in ATTRIBUTED[1:3]))
    log("  the colliding pair on both sides, delta = gamma/gamma* - 1:")
    for delta in CERTIFICATE_DELTAS:
        for sign in (-1.0, 1.0):
            _, pa, pb, _ = squared_gap(LHs, Ds, g_star * (1.0 + sign * delta), lam_star)
            pa, pb = sorted((pa, pb), key=lambda z: (z.real, z.imag))
            log("     delta = {0:+.0e}:  {1}  {2}   half-split/sqrt|delta| = {3:.5f}".format(
                sign * delta, format_eig(pa), format_eig(pb), abs(pa - pb) / 2.0 / np.sqrt(delta)))
    coefficients = []
    for delta in SQRT_LAW_DELTAS:
        _, pa, pb, _ = squared_gap(LHs, Ds, g_star * (1.0 + delta), lam_star)
        coefficients.append("{0:.5f} ({1:.0e})".format(abs(pa - pb) / 2.0 / np.sqrt(delta), delta))
    log("  half-split/sqrt(delta) above gamma*: " + ", ".join(coefficients))
    signatures = []
    for delta in CERTIFICATE_DELTAS:
        values, vectors = eig(LH + g_star * (1.0 - delta) * D)
        order = np.argsort(np.abs(values - lam_star))[:2]
        readings = []
        for index in order:
            v = vectors[:, index] / np.linalg.norm(vectors[:, index])
            readings.append(np.vdot(v, G_METRIC @ v).real)
        signatures.append("{0:+.4f} / {1:+.4f} at delta = -{2:.0e}".format(
            max(readings), min(readings), delta))
    log("  Krein signatures v^+ G v of the pair below gamma*: " + ", ".join(signatures))
    L_star = LH + g_star * D
    log("  metric check: |G LH G - LH| = {0:.1e}, |G D G + D| = {1:.1e}, |L - L^T| = {2:.1e}".format(
        np.max(np.abs(G_METRIC @ LH @ G_METRIC - LH)), np.max(np.abs(G_METRIC @ D @ G_METRIC + D)),
        np.max(np.abs(L_star - L_star.T))))
    at_crit = eigvals(build_L(NPC, GC, J_bridge=J_BRIDGE).toarray())
    pa, pb = at_crit[np.argsort(np.abs(at_crit - lam_star))[:2]]
    BELOW[J_BRIDGE] = ((GC - g_star) / g_star < 0.0 and abs(pa - pb) > 1e3 * pair_gap
                       and float(np.max(at_crit.real)) < 1e-12)
    log("  bisected gamma_crit = {0:.9f}, (gamma_crit - gamma*)/gamma* = {1:.1e}: there the pair"
        " {2} {3} is still split on the axis by {4:.2e}, and max Re lambda = {5:.1e} is rounding"
        " noise".format(GC, (GC - g_star) / g_star, format_eig(pa), format_eig(pb), abs(pa - pb),
                        float(np.max(at_crit.real))))
    log()
    if (within == 2 and sector_sv[0] < 1e-3 * pair_gap and sector_sv[1] > 1e3 * pair_gap
            and sv[0] < 1e-3 * pair_gap and sv[1] > 1e3 * pair_gap and f_below < 0.0 < f_above):
        CERTIFIED.append(J_BRIDGE)
if CERTIFIED == list(BRIDGES):
    log("At J_bridge = 1.0 and 1.9: geometric multiplicity 1 at algebraic multiplicity 2, the colliding pair")
    log("forms a 2x2 Jordan block, a second-order exceptional point (EP2) on the real gamma axis at the threshold.")
    log("Below gamma* two distinct imaginary eigenvalues with opposite Krein signatures approach each other")
    log("(a Krein collision); above it they leave the axis as a mirror pair, with the same")
    log("half-split c*sqrt|delta| on both sides.")
    log("A Hopf crossing moves a simple eigenvalue across Re = 0 with a finite slope. At these two couplings")
    log("GK keeps every simple eigenvalue of a block on the axis, and the departure has infinite slope:")
    log("the threshold is not a Hopf crossing.")
else:
    log("The certificate gates failed at J_bridge = {0}.".format(
        ", ".join(str(j) for j in BRIDGES if j not in CERTIFIED)))
if all(BELOW.values()):
    log("At both couplings the bisected gamma_crit sits slightly below the collision: the approaching pair")
    log("is ill-conditioned, and the rounding noise in its real parts trips the 1e-12 test while the pair is")
    log("still split on the axis.")
else:
    log("The bisected gamma_crit does not sit below a still-split pair at J_bridge = {0}.".format(
        ", ".join(str(j) for j, ok in BELOW.items() if not ok)))
log()

log("-" * 78)
log("2. Sector structure: the first block to go unstable, and its copies in the full L")
log("-" * 78)
log()
G_OP, F_OP, R_OP, S_OP, SWAP_OP = symmetry_superoperators(N_TOTAL)
PI_OP = palindrome_superoperator(N_TOTAL)
PI_DAG = PI_OP @ SWAP_OP
LH, D = generator_parts(NPC, 1.46)
L_TEST = LH + 0.37 * D
log("Symmetry checks at J_bridge = 1.46, gamma = 0.37 (antilinear maps written A o conj, so that")
log("A conj(L) A^-1 = -L means lambda -> -lambda*):")
log("  S = G o R:                  |S L - L S| = {0:.1e}, |S^2 - 1| = {1:.1e}, |S F - F S| = {2:.1e}".format(
    np.max(np.abs(S_OP @ L_TEST - L_TEST @ S_OP)), np.max(np.abs(S_OP @ S_OP - np.eye(S_OP.shape[0]))),
    np.max(np.abs(S_OP @ F_OP - F_OP @ S_OP))))
log("  F = X^N (.) X^N:            |F L F - L| = {0:.1e}, |Pi^2 - F| = {1:.1e} (F118's centre)".format(
    np.max(np.abs(F_OP @ L_TEST @ F_OP - L_TEST)), np.max(np.abs(PI_OP @ PI_OP - F_OP))))
log("  R = (rho -> rho X^N), F118: |R L(gamma) R - L(-gamma)| = {0:.1e}".format(
    np.max(np.abs(R_OP @ L_TEST @ R_OP - (LH - 0.37 * D)))))
log("  GK:                         |G conj(L) G + L| = {0:.1e}".format(
    np.max(np.abs(G_OP @ L_TEST.conj() @ G_OP + L_TEST))))
log("  Pi o dagger:                |A conj(L) A^-1 + L| = {0:.1e}, A = Pi SWAP; |A - R| = {1:.1e},".format(
    np.max(np.abs(PI_DAG @ L_TEST.conj() @ PI_DAG.conj().T + L_TEST)), np.max(np.abs(PI_DAG - R_OP))))
log("                              so Pi o dagger = R o conj")
log("  their product:              |G conj(A) - S| = {0:.1e}, so GK o (Pi o dagger) = S".format(
    np.max(np.abs(G_OP @ PI_DAG.conj() - S_OP))))
BLOCKS = FRAME[0]


def maps_block(M, blk):
    targets = set()
    d = 2 ** N_TOTAL
    for i in BLOCKS[blk]:
        for j in np.nonzero(np.abs(M[:, i]) > 0.5)[0]:
            a, b = divmod(j, d)
            targets.add((bin(a).count("1"), bin(b).count("1")))
    return targets


S_MAP = all(maps_block(S_OP, (p, q)) == {(p, N_TOTAL - q)} for p, q in BLOCKS)
GK_MAP = all(maps_block(G_OP, blk) == {blk} for blk in BLOCKS)
PD_MAP = all(maps_block(PI_DAG, (p, q)) == {(p, N_TOTAL - q)} for p, q in BLOCKS)
log("  block maps over all 25 blocks: S sends (p,q) to (p,4-q): {0}; GK sends every block to itself:"
    " {1}; Pi o dagger sends (p,q) to (p,4-q): {2}".format(
        *("holds" if ok else "fails" for ok in (S_MAP, GK_MAP, PD_MAP))))
log()
log("Per coupling: the first block to go unstable, its orbit under F, S and Hermitian conjugation,")
log("the number of blocks its orbit under F and S holds (copies), the collision gamma* inside it,")
log("and the full L at lambda*: eigenvalues within 1e-5, singular values of L - lambda* I below 1e-10")
log("of the largest, eigenvalues with Re > 1e-9 at gamma*(1 + 1e-6).")
log()
log("  {0:>8}  {1:<16} {2:<50} {3:>6} {4:>12} {5:>9} {6:>12} {7:>14} {8:>12}".format(
    "J_bridge", "first block", "orbit", "copies", "gamma*", "|Im l*|", "within 1e-5", "sv at rounding",
    "Re > 0 above"))
ROWS_HOLD = []
for J_BRIDGE in SECTOR_BRIDGES:
    LH, D = generator_parts(NPC, J_BRIDGE)
    _, _, NAME, g_star, lam_star, _, _, _ = sector_survey(LH, D, FRAME, NPC)
    within, sv, positive = full_l_at_collision(LH, D, g_star, lam_star)
    at_rounding = int(np.sum(sv < 1e-10 * sv[-1]))
    copies = copies_of(NAME, N_TOTAL)
    ROWS_HOLD.append((J_BRIDGE, within == 2 * copies and at_rounding == copies and positive == 2 * copies))
    log("  {0:>8.2f}  {1:<16} {2:<50} {3:>6d} {4:>12.9f} {5:>9.3f} {6:>12d} {7:>14d} {8:>12d}".format(
        J_BRIDGE, canonical(NAME, N_TOTAL), orbit_label(NAME, N_TOTAL), copies, g_star, abs(lam_star.imag),
        within, at_rounding, positive))
log()
FAILED = [j for j, ok in ROWS_HOLD if not ok]
if not FAILED:
    log("In every row the full L holds one Jordan pair per copy: eigenvalues within 1e-5 = 2 x copies,")
    log("singular values at the rounding level = copies, eigenvalues with Re > 0 just above = 2 x copies.")
else:
    log("The count 2 x copies / copies / 2 x copies fails at J_bridge = {0}.".format(
        ", ".join("{0}".format(j) for j in FAILED)))
log()


def block_threshold_at(j_bridge, names):
    lh, d = generator_parts(NPC, j_bridge)
    parts = sectors(lh, d, FRAME, NPC)
    return {name: sector_threshold(*parts[name]) for name in names}


W_LO = brentq(lambda j: np.subtract(*block_threshold_at(j, ["(0,1)", "(1,2)"]).values()),
              1.43, 1.44, xtol=1e-10)
W_HI = brentq(lambda j: np.subtract(*block_threshold_at(j, ["(0,1)", "(2,2) flip-odd"]).values()),
              1.49, 1.50, xtol=1e-10)
INDEX_01 = BLOCKS[(0, 1)]
D_01 = np.real(np.diag(generator_parts(NPC, 1.0)[1][np.ix_(INDEX_01, INDEX_01)]))
SITES_01 = [N_TOTAL - 1 - int(np.log2(i % (2 ** N_TOTAL))) for i in INDEX_01]
log("The (0,1) block holds the coherences between the empty state and one excitation. Its gain-loss part")
log("is diagonal, " + ", ".join("{0:+.0f} gamma on site {1}".format(v, k) for v, k in zip(D_01, SITES_01))
    + " (sites 0 and 1 are chain A),")
log("so it is a four-site hopping chain with loss 2 gamma on chain A and gain 2 gamma on chain B.")
log()
log("The (0,1) orbit goes first for J_bridge between {0:.6f} (where it meets the (1,2) orbit, at gamma ="
    " {1:.6f})".format(W_LO, block_threshold_at(W_LO, ["(0,1)"])["(0,1)"]))
log("and {0:.6f} (where it meets the flip-odd half of (2,2), at gamma = {1:.6f}).".format(
    W_HI, block_threshold_at(W_HI, ["(0,1)"])["(0,1)"]))
log()
log("The two (2,2) halves near the sampled maximum, with the first sector over all blocks:")
log()
log("  {0:>8}  {1:>10} {2:>10}   {3:<16} {4:>10}".format(
    "J_bridge", "flip-even", "flip-odd", "first block", "gamma_c"))
CROSSING_ROWS = []
for J_BRIDGE in CROSSING_BRIDGES:
    LH, D = generator_parts(NPC, J_BRIDGE)
    _, th = all_thresholds(LH, D, FRAME, NPC)
    first, _ = first_of(th)
    CROSSING_ROWS.append((canonical(first, N_TOTAL), th[first]))
    log("  {0:>8.3f}  {1:>10.6f} {2:>10.6f}   {3:<16} {4:>10.6f}".format(
        J_BRIDGE, th["(2,2) flip-even"], th["(2,2) flip-odd"], canonical(first, N_TOTAL), th[first]))


def parity_gap(j_bridge):
    th = block_threshold_at(j_bridge, ["(2,2) flip-even", "(2,2) flip-odd"])
    return th["(2,2) flip-even"] - th["(2,2) flip-odd"]


J_CROSS = brentq(parity_gap, 1.95, 1.953, xtol=1e-9)
LH, D = generator_parts(NPC, J_CROSS)
_, TH = all_thresholds(LH, D, FRAME, NPC)
LOWEST_OTHER = min((name for name in TH if not name.startswith("(2,2)")), key=TH.get)
log()
log("  the halves cross at J_bridge = {0:.6f}, gamma = {1:.6f}; the lowest other sector there is"
    " {2} at {3:.6f}".format(J_CROSS, TH["(2,2) flip-even"], orbit_label(LOWEST_OTHER, N_TOTAL),
                             TH[LOWEST_OTHER]))
log("  at J_bridge = 2.0 the threshold, {0:.6f}, is {1:.2%} below that crossing value".format(
    CROSSING_ROWS[-1][1], (TH["(2,2) flip-even"] - CROSSING_ROWS[-1][1]) / TH["(2,2) flip-even"]))
LEADS = [(lead, value) for lead, value in CROSSING_ROWS]
SWITCH = [lead for lead, _ in LEADS].index("(2,2) flip-even") if "(2,2) flip-even" in [l for l, _ in LEADS] else None
KINK = (SWITCH is not None and all(l == "(2,2) flip-odd" for l, _ in LEADS[:SWITCH])
        and all(l == "(2,2) flip-even" for l, _ in LEADS[SWITCH:])
        and all(b > a for (_, a), (_, b) in zip(LEADS[:SWITCH], LEADS[1:SWITCH]))
        and all(b < a for (_, a), (_, b) in zip(LEADS[SWITCH:], LEADS[SWITCH + 1:]))
        and max(v for _, v in LEADS) < TH["(2,2) flip-even"])
if KINK:
    log("Across these samples the threshold rises with the flip-odd half up to the crossing and falls with the")
    log("flip-even half after it: the sampled maximum is a kink where the first (2,2) half switches spin-flip parity.")
else:
    log("The samples do not show a single switch from the flip-odd to the flip-even half at the maximum.")
log()

log("-" * 78)
log("3. Zeros of the threshold: couplings where an even and an odd level already coincide")
log("-" * 78)
log()
log("At gamma = 0 each block of L_H splits into G-even and G-odd levels, and further by F and S where")
log("they act inside the block. The gain-loss part anticommutes with G, so it couples an even level to")
log("an odd one. Where the two have the same frequency omega (an eigenvalue of i*L_H), it splits them at")
log("first order, lambda = -i*omega +/- s*gamma + O(gamma^2), so max Re lambda = s*gamma + O(gamma^2) > 0")
log("for small gamma > 0: gamma_crit = 0, the pair is semisimple at gamma = 0, and there is no")
log("exceptional point.")
log()
LEVELS = symbolic_levels(NPC)
for n_up in (1, 2):
    log("  H on the popcount-{0} sector: ".format(n_up) + ", ".join(
        "{0} (P_site{1:+d}{2})".format(fmt_poly(e), r, "" if x == 0 else ", X{0:+d}".format(x))
        for e, r, x in LEVELS[n_up]))
log("  (x = J_bridge / J with J = 1; P_site the site reversal, X the spin flip; popcount 0 and 4: x + 2;")
log("  popcount 3 repeats popcount 1)")
EXACT, SPURIOUS, IDENTICAL = exact_degeneracy_couplings(LEVELS, N_TOTAL, REPRESENTATIVES)
log()
log("Exact search over every pair of opposite-G levels in one symmetry sector of the blocks {0}".format(
    ", ".join("({0},{1})".format(*b) for b in REPRESENTATIVES)))
log("(every other block is an image of these under F, S or Hermitian conjugation, except the four 1x1")
log("corner blocks (0,0), (0,4), (4,0), (4,4), which hold a single level and so no pair), radicals eliminated")
log("by resultants, real roots x > 0 of every factor; a root counts where the two levels agree to 1e-25")
log("at 40 digits. The spurious roots squaring adds leave the levels at least {0:.3f} apart; an".format(SPURIOUS))
log("identically equal opposite-G pair: {0}.".format(IDENTICAL))
log()

# the numeric route: level crossings on a grid
OPS = (G_OP, F_OP, S_OP, F_OP @ S_OP)
BASES = {blk: block_sector_bases(blk, FRAME, OPS, N_TOTAL) for blk in REPRESENTATIVES}
LH0, _ = generator_parts(NPC, 0.0)
LH1 = generator_parts(NPC, 1.0)[0] - LH0
GRID = np.round(np.arange(0.002, 12.0 + 1e-9, 0.0025), 10)


def grid_levels(j_bridge):
    lh = LH0 + j_bridge * LH1
    out = {}
    for blk in REPRESENTATIVES:
        index = BLOCKS[blk]
        lhb = lh[np.ix_(index, index)]
        for key, pair in BASES[blk].items():
            if 1 in pair and -1 in pair:
                out[(blk, key)] = (sector_frequencies(lhb, pair[1]), sector_frequencies(lhb, pair[-1]))
    return out


NUMERIC = []
previous = grid_levels(GRID[0])
for j0, j1 in zip(GRID[:-1], GRID[1:]):
    current = grid_levels(j1)
    for key, (wp, wm) in current.items():
        wp0, wm0 = previous[key]
        change = np.sign(wp0[:, None] - wm0[None, :]) * np.sign(wp[:, None] - wm[None, :]) < 0
        for i, k in zip(*np.nonzero(change)):
            def gap(j, key=key, i=i, k=k):
                lv = grid_levels(j)[key]
                return lv[0][i] - lv[1][k]
            root = brentq(gap, j0, j1, xtol=1e-14, rtol=8.9e-16)
            NUMERIC.append((root, key, grid_levels(root)[key][0][i]))
    previous = current
NUMERIC.sort()
NUMERIC_J = []
for j, _, _ in NUMERIC:
    if not NUMERIC_J or j - NUMERIC_J[-1] > 1e-9:
        NUMERIC_J.append(j)
EXACT_J = sorted(EXACT)
log("Numeric route: level crossings of opposite G parity in one sector, on a grid of step 0.0025 over")
log("[0.002, 12], each refined by brentq: J_bridge = " + ", ".join("{0:.12f}".format(j) for j in NUMERIC_J))
log("The exact search finds {0} couplings in all of x > 0; the grid finds {1} in [0.002, 12], at the same"
    " places to {2:.0e}.".format(len(EXACT_J), len(NUMERIC_J), max(
        min(abs(j - e) for e in EXACT_J) for j in NUMERIC_J)))
log()
log("The table gives, per coupling, each sector holding such a pair: its frequency omega, the dimension")
log("of the degenerate subspace, the first-order slope s, and the V slope |d(omega_even - omega_odd)/dJ|/(2s)")
log("at which that sector's threshold leaves zero; then max Re lambda / gamma of the full L at")
log("gamma = 1e-2, 1e-4, 1e-6.")
log()
log("  {0:>15}  {1:<40} {2:<18} {3:>9} {4:>3} {5:>9} {6:>8}   {7}".format(
    "J_bridge", "exact", "block / sector", "omega", "dim", "s", "V slope", "max Re lambda / gamma"))
ZERO_ROWS = []
for j_key in EXACT_J:
    entries = EXACT[j_key]
    factor, root = entries[0][2], entries[0][3]
    j0 = float(root.evalf(30))
    if sympy.degree(factor, JB) <= 2:
        exact = "x = {0}".format(sympy.nsimplify(root) if sympy.degree(factor, JB) == 1 else root)
    else:
        exact = "real root of {0}".format(fmt_poly(factor))
    LH, D = generator_parts(NPC, j0)
    measured = [float(np.max(eigvals(LH + g * D).real)) / g for g in SLOPE_GAMMAS]
    first = True
    s_max = 0.0
    v_min = np.inf
    for blk, label in sorted({(e[0], e[1]) for e in entries}):
        key = tuple(v for _, v in label)
        LHb = LH[np.ix_(BLOCKS[blk], BLOCKS[blk])]
        omegas = [w for w in sector_frequencies(LHb, BASES[blk][key][1])
                  if np.min(np.abs(sector_frequencies(LHb, BASES[blk][key][-1]) - w)) < 1e-7]
        for omega in omegas:
            s, count = first_order_slope(LH, D, BLOCKS[blk], BASES[blk][key], omega)
            s_max = max(s_max, s)
            spread = []
            for h in (-1e-6, 1e-6):
                lhb = generator_parts(NPC, j0 + h)[0][np.ix_(BLOCKS[blk], BLOCKS[blk])]
                even = sector_frequencies(lhb, BASES[blk][key][1])
                odd = sector_frequencies(lhb, BASES[blk][key][-1])
                spread.append(even[np.argmin(np.abs(even - omega))] - odd[np.argmin(np.abs(odd - omega))])
            v_slope = abs(spread[1] - spread[0]) / 2e-6 / (2.0 * s)
            if count == 2:
                v_min = min(v_min, v_slope)
            where = "({0},{1}) ".format(*blk) + " ".join("{0}{1:+d}".format(n, v) for n, v in label)
            log("  {0:>15} {1:<40} {2:<18} {3:>+9.5f} {4:>3d} {5:>9.6f} {6:>8.4f}   {7}".format(
                "{0:.12f}".format(j0) if first else "", exact if first else "", where, omega, count, s,
                v_slope, "  ".join("{0:.6f}".format(m) for m in measured) if first else ""))
            first = False
    ZERO_ROWS.append((j0, exact, s_max, measured, v_min))
log()
LH, D = generator_parts(NPC, 1.0)
CONTROL = [float(np.max(eigvals(LH + g * D).real)) for g in SLOPE_GAMMAS]
log("Control at a generic coupling, J_bridge = 1.0: max Re lambda / gamma = " + "  ".join(
    "{0:.1e}".format(m / g) for m, g in zip(CONTROL, SLOPE_GAMMAS)) + (
    " (max Re lambda below 1e-12 at all three: rounding, no linear onset)" if max(CONTROL) < 1e-12
    else " (max Re lambda is not at the rounding level)"))
log()
LINEAR = all(abs(m[1] / s - 1.0) < 1e-6 and abs(m[2] / s - 1.0) < 1e-6 for _, _, s, m, _ in ZERO_ROWS)
POSITIVITY_GAMMAS = np.geomspace(1e-6, 2.0, 61)
LOWEST = (np.inf, None, None)
for j0, _, _, _, _ in ZERO_ROWS:
    LH, D = generator_parts(NPC, j0)
    for g in POSITIVITY_GAMMAS:
        ratio = float(np.max(eigvals(LH + g * D).real)) / g
        if ratio < LOWEST[0]:
            LOWEST = (ratio, j0, g)
if LINEAR:
    log("At each of these couplings max Re lambda / gamma equals the largest first-order s to 1e-6 from")
    log("gamma = 1e-4 down: gamma_crit = 0 there, and any small gain makes the bridge unstable.")
else:
    log("The linear onset fails at one of these couplings; see the table.")
log("Over 61 gammas from 1e-6 to 2 at each of them, max Re lambda {0} (smallest max Re lambda / gamma:".format(
    "stays positive" if LOWEST[0] > 0 else "does not stay positive"))
log("{0:.4f}, at J_bridge = {1:.6f}, gamma = {2:.3g}).".format(*LOWEST))
log()
V_ROWS = []
for j0, exact, _, _, v_min in ZERO_ROWS:
    values = []
    for h in V_OFFSETS:
        LH, D = generator_parts(NPC, j0 + h)
        _, th = all_thresholds(LH, D, FRAME, NPC)
        values.append(min(th.values()))
    V_ROWS.append((j0, values, v_min))
DOUBLING = all(1.9 < v[0] / v[1] < 2.1 and 1.9 < v[3] / v[2] < 2.1 for _, v, _ in V_ROWS)
log("Around each zero, gamma_crit (all blocks) at J_bridge = J0 + offset; {0}".format(
    "it doubles from 0.005 to 0.01 on both sides, a linear rise:" if DOUBLING
    else "it does not double from 0.005 to 0.01 on every side:"))
log()
log("  {0:>15}  {1:>10} {2:>10} {3:>10} {4:>10}   {5:>11} {6:>11} {7:>13}".format(
    "J_bridge", "-0.01", "-0.005", "+0.005", "+0.01", "slope left", "slope right", "V slope (s)"))
for j0, values, v_min in V_ROWS:
    log("  {0:>15.12f}  {1:>10.6f} {2:>10.6f} {3:>10.6f} {4:>10.6f}   {5:>11.4f} {6:>11.4f} {7:>13.4f}".format(
        j0, values[0], values[1], values[2], values[3], values[1] / 0.005, values[2] / 0.005, v_min))
MATCH = all(abs(0.5 * (v[1] + v[2]) / 0.005 / v_min - 1.0) < 0.02 for _, v, v_min in V_ROWS)
log("The last column is the smallest V slope of the table above; it {0} the measured mean of the two".format(
    "matches" if MATCH else "does not match"))
log("sides to 2 %{0}".format(
    ": near each zero the threshold is the frequency mismatch of the even-odd pair divided by twice"
    if MATCH else "."))
if MATCH:
    log("its coupling, |omega_even - omega_odd|/(2s).")
log()
log("Where the finer sweep's low points and the reversal of its strong-bridge tail sit:")
log()
log("  {0:>8}  {1:>10}   {2:<16} {3:>15} {4:>10}".format(
    "J_bridge", "gamma_c", "first block", "nearest zero", "offset"))
for J_BRIDGE in FINER_GRID_POINTS:
    LH, D = generator_parts(NPC, J_BRIDGE)
    _, th = all_thresholds(LH, D, FRAME, NPC)
    first, value = first_of(th)
    j0 = min((row[0] for row in ZERO_ROWS), key=lambda z: abs(z - J_BRIDGE))
    log("  {0:>8.2f}  {1:>10.6f}   {2:<16} {3:>15.12f} {4:>+10.4f}".format(
        J_BRIDGE, value, canonical(first, N_TOTAL), j0, J_BRIDGE - j0))
TAKE = brentq(lambda j: np.subtract(*block_threshold_at(j, ["(2,2) flip-even", "(1,2)"]).values()),
              4.4, 4.6, xtol=1e-10)
log()
TAKE_GAMMA = block_threshold_at(TAKE, ["(1,2)"])["(1,2)"]
LAST_ZERO = ZERO_ROWS[-1]
log("Beyond the zero near 3.3289 the flip-even half of (2,2) rises until it meets the (1,2) orbit")
log("at J_bridge = {0:.6f}, gamma = {1:.6f}. It leaves the zero at the V slope {2:.4f} and rises on".format(
    TAKE, TAKE_GAMMA, LAST_ZERO[4]))
log("average at {0:.4f} per unit of J_bridge up to that meeting.".format(TAKE_GAMMA / (TAKE - LAST_ZERO[0])))
log()

log("-" * 78)
log("4. Finite-offset diagnostics above the bisected gamma_crit")
log("-" * 78)
log()
for J_BRIDGE in BRIDGES:
    GC = G_CRIT[J_BRIDGE]
    log("J_bridge = {0}   gamma_crit = {1:.9f}".format(J_BRIDGE, GC))
    log()
    log("  {0:>8} {1:>15} {2:>15} {3:>15} {4:>13} {5:>13}".format(
    "delta", "max Re lambda", "Re/sqrt(delta)", "gap to -lambda*", "nearest gap", "Petermann K"))
    coefficients = []
    ks = []
    last_gap = None
    for DELTA in DELTAS:
        lam, gap, K = leading_mode(NPC, GC * (1.0 + DELTA), J_BRIDGE)
        if gap <= SIMPLE_MODE_GAP_TOL:
            raise RuntimeError(
                "near-threshold Petermann row is not a simple mode: "
                "J_bridge={0}, delta={1}, nearest gap={2}".format(J_BRIDGE, DELTA, gap))
        coefficients.append(lam.real / np.sqrt(DELTA))
        ks.append(K)
        last_gap = 2.0 * abs(lam.real)
        log("  {0:>8.0e} {1:>15.8f} {2:>15.4f} {3:>15.3e} {4:>13.3e} {5:>13.4g}".format(
            DELTA, lam.real, lam.real / np.sqrt(DELTA), last_gap, gap, K))
    log("  near-threshold simplicity gate: nearest gap > 1e-10 for every reported K")
    spread = max(coefficients) / min(coefficients)
    decades = np.sqrt(DELTAS[0] / DELTAS[-1])
    offset = (GAMMA_STAR[J_BRIDGE] - GC) / GAMMA_STAR[J_BRIDGE]
    log()
    log("  1. Re lambda / sqrt(delta) spans {0:.4f} over three decades of delta (four samples)."
        " A transversal crossing (Re ~ delta) would make this ratio grow by"
        " {1:.0f} across the same range.".format(spread, decades))
    log("     measured from the bisected gamma_crit, which sits {0:.1e} (relative) below gamma*;"
        " that offset is {1:.2%} of delta at {2:.0e} and {3:.2%} at {4:.0e}.".format(
            offset, offset / DELTAS[-2], DELTAS[-2], offset / DELTAS[-1], DELTAS[-1]))
    log("  2. the gap to the across-axis partner -lambda* is 2*abs(Re lambda) = {0:.2e};"
        " it is derived from max Re and supplies no independent EP evidence."
        .format(last_gap))
    ratios = [ks[i + 1] / ks[i] for i in range(len(ks) - 1)]
    log("  3. Petermann K multiplies by {0} per sampled decade of delta, the K ~ 1/delta"
        " of an EP2; the Jordan rank itself is read in section 1.".format(
            ", ".join("{0:.2f}".format(r) for r in ratios)))
    for below in (-1e-3, -1e-2, -1e-1):
        log("     below threshold, delta = {0:>+.0e}: max Re lambda = {1:.2e}".format(
            below, max_re(NPC, GC * (1.0 + below), J_BRIDGE)))
    log()

log("-" * 78)
log("5. What a coarse scan sees instead")
log("-" * 78)
log()
log("A finite grid reads sampling-dependent Petermann values for simple selected modes.")
log("A single-vector K is omitted wherever the selected eigenvalue is degenerate:")
log()
GC = G_CRIT[1.0]
log("  {0:>17} {1:>43}".format("gamma/gamma_crit", "Petermann reading"))
for RATIO in (1.001, 1.01, 1.05, 1.2, 1.46, 1.5, 2.0, 3.0):
    _, gap, K = leading_mode(NPC, GC * RATIO, 1.0)
    if gap <= SIMPLE_MODE_GAP_TOL:
        reading = "single-vector K not reported for a degenerate eigenspace"
    else:
        reading = "K={0:.4g} (nearest gap {1:.2e})".format(K, gap)
    log("  {0:>17.3f} {1:>43}".format(RATIO, reading))
log()
log("Degenerate rows do not support a basis-invariant single-eigenvector conditioning value.")
log("This coarse scan selects max Re(lambda) independently at each gamma;")
log("the colliding pair itself is followed only in section 1.")


log()
log("-" * 78)
log("6. The first block over a sweep of couplings")
log("-" * 78)
log()


def symmetry_parts(blk):
    """Orthonormal bases (block-local, real) of the joint eigenspaces of the symmetries of L(gamma)
    acting inside the block: F on (2,2), S on (p,2), F S on (2,q). G is not one of them."""
    index = BLOCKS[blk]
    half = N_TOTAL // 2
    acting = []
    if blk == (half, half):
        acting.append(F_OP)
    if blk[1] == half:
        acting.append(S_OP)
    if blk[0] == half and blk[1] != half:
        acting.append(F_OP @ S_OP)
    if not acting:
        return {(): np.eye(len(index))}
    mats = [M[np.ix_(index, index)] for M in acting]
    rng = np.random.default_rng(5)
    combo = sum(rng.uniform(1.0, 2.0) * 2.0 ** k * M for k, M in enumerate(mats))
    _, V = np.linalg.eigh(combo)
    out = {}
    for k in range(V.shape[1]):
        v = V[:, k]
        out.setdefault(tuple(int(round(v @ M @ v)) for M in mats), []).append(k)
    return {key: V[:, cols] for key, cols in out.items()}


D_FULL = generator_parts(NPC, 0.0)[1]
SWEEP_PIECES = {}
for blk in REPRESENTATIVES:
    index = BLOCKS[blk]
    for key, B in symmetry_parts(blk).items():
        name = "({0},{1})".format(*blk)
        if blk == (NPC, NPC):
            name += " flip-even" if key[0] == 1 else " flip-odd"
        SWEEP_PIECES.setdefault(name, []).append(
            (B.T @ LH0[np.ix_(index, index)] @ B, B.T @ LH1[np.ix_(index, index)] @ B,
             B.T @ D_FULL[np.ix_(index, index)] @ B))
SWEEP_GAMMAS = np.geomspace(1e-4, 2.0, 60)


def piece_threshold(A, Dp, tol=1e-9):
    previous = 0.0
    for g in SWEEP_GAMMAS:
        if np.max(eigvals(A + g * Dp).real) > tol:
            lo, hi = previous, g
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if np.max(eigvals(A + mid * Dp).real) > tol:
                    hi = mid
                else:
                    lo = mid
            return 0.5 * (lo + hi)
        previous = g
    return np.inf


SWEEP_J = np.round(np.arange(0.005, 20.0 + 1e-9, 0.005), 6)
SWEEP = []
SWEEP_FIRST_PIECE = []
for j in SWEEP_J:
    best = None
    for name, pieces in SWEEP_PIECES.items():
        for a0, a1, dp in pieces:
            value = piece_threshold(a0 + j * a1, dp)
            if best is None or value < best[0]:
                best = (value, name, a0 + j * a1, dp)
    SWEEP.append((float(j), best[1], best[0]))
    SWEEP_FIRST_PIECE.append(best)
RUNS = []
for j, first, value in SWEEP:
    if RUNS and RUNS[-1][0] == first:
        RUNS[-1][2] = j
        RUNS[-1][3].append((j, value))
    else:
        RUNS.append([first, j, j, [(j, value)]])


def run_shape(points):
    values = [v for _, v in points]
    turns = []
    for k in range(1, len(values) - 1):
        if values[k] < values[k - 1] and values[k] < values[k + 1]:
            turns.append("a minimum {0:.6f} at {1:.3f}".format(values[k], points[k][0]))
        if values[k] > values[k - 1] and values[k] > values[k + 1]:
            turns.append("a maximum {0:.6f} at {1:.3f}".format(values[k], points[k][0]))
    ends = "{0:.6f} to {1:.6f}".format(values[0], values[-1])
    if not turns:
        return "{0} from {1}".format("rises" if values[-1] > values[0] else "falls", ends)
    return "{0}, through {1}".format(ends, " and ".join(turns))


log("Grid: J_bridge = 0.005 to 20 in steps of 0.005 ({0} couplings). Per block the threshold is the first".format(
    len(SWEEP_J)))
log("gamma of a 60-point geometric grid from 1e-4 to 2 at which max Re lambda exceeds 1e-9, refined by 40")
log("bisection steps. The blocks are the representatives (0,1), (0,2), (1,1), (1,2) and the two (2,2)")
log("halves, each split further by the symmetries of L(gamma) acting inside it.")
log()
log("  {0:<16} {1:>7} {2:>7} {3:>6}   {4}".format("first block", "from", "to", "copies", "threshold along the run"))
for name, start_j, end_j, points in RUNS:
    log("  {0:<16} {1:>7.3f} {2:>7.3f} {3:>6d}   {4}".format(
        name, start_j, end_j, copies_of(name, N_TOTAL), run_shape(points)))
log()
FIRSTS = {name for name, _, _, _ in RUNS}
FOUR = [(name, a, b) for name, a, b, _ in RUNS if copies_of(name, N_TOTAL) == 4]
log("On this grid (1,1) is {0} the first block, and the first blocks with four copies are: {1}.".format(
    "never" if "(1,1)" not in FIRSTS else "sometimes",
    ", ".join("{0} from {1:.3f} to {2:.3f}".format(n, a, b) for n, a, b in FOUR) or "none"))
VALUES = np.array([v for _, _, v in SWEEP])
IMAX = int(np.argmax(VALUES))
PEAKS = [(SWEEP_J[k], VALUES[k]) for k in range(1, len(VALUES) - 1)
         if VALUES[k] > VALUES[k - 1] and VALUES[k] > VALUES[k + 1]]
log("The largest threshold on the grid is {0:.6f}, at J_bridge = {1:.3f}; the crossing of the (2,2) halves".format(
    VALUES[IMAX], SWEEP_J[IMAX]))
log("(section 2, gamma = {0:.6f}) {1} every grid value.".format(
    TH["(2,2) flip-even"], "lies above" if TH["(2,2) flip-even"] > VALUES.max() else "does not lie above"))
log("Local maxima of the grid threshold: " + ", ".join(
    "{0:.6f} at {1:.3f}".format(v, j) for j, v in PEAKS) + ".")
log()


def read_collision(A, Dp, threshold):
    """Locate the pair's collision by the sign change of its real squared gap and read it:
    (gamma*, pair gap, sigma_1, sigma_2, distance of the third eigenvalue, relative mismatch of
    half-split / sqrt|delta| below and above at delta = 1e-6); None where f never changes sign."""
    w = eigvals(A + threshold * (1.0 + 1e-4) * Dp)
    anchor = complex(0.0, w[int(np.argmax(w.real))].imag)

    def f(gamma):
        return squared_gap(A, Dp, gamma, anchor)[0].real

    lo, hi = threshold * (1.0 - 1e-6), threshold * (1.0 + 1e-4)
    k = 0
    while f(lo) >= 0.0 and k < 6:
        k += 1
        lo = threshold * max(0.0, 1.0 - 1e-6 * 10 ** k)
    if f(lo) >= 0.0 or f(hi) <= 0.0:
        return None
    g = brentq(f, lo, hi, xtol=1e-17, rtol=8.9e-16, maxiter=200)
    _, a, b, spectrum = squared_gap(A, Dp, g, anchor)
    lam = 0.5 * (a + b)
    sv = np.sort(svdvals(A + g * Dp - lam * np.eye(A.shape[0])))
    split = []
    for sign in (-1.0, 1.0):
        _, pa, pb, _ = squared_gap(A, Dp, g * (1.0 + sign * 1e-6), lam)
        split.append(abs(pa - pb) / 2.0 / 1e-3)
    return (g, abs(a - b), sv[0], sv[1], np.sort(np.abs(spectrum - lam))[2], abs(split[0] / split[1] - 1.0))


READS = [(j, piece[1], read_collision(piece[2], piece[3], piece[0])) for (j, _, _), piece in zip(SWEEP, SWEEP_FIRST_PIECE)]
EP2 = [(j, r) for j, _, r in READS if r is not None and r[2] <= 1e-3 * max(r[1], 1e-300)
       and r[3] >= 1e3 * r[1] and r[5] < 1e-3]
OTHER = [(j, name, r) for j, name, r in READS if (j, r) not in EP2]
log("At every grid coupling the colliding pair of the sector that goes first (a symmetry sector of the")
log("first block) is located by the sign change of its squared gap and read there. It is an EP2 (inside")
log("that sector one singular value of L - lambda* I at the rounding level, the next above 1e3 times the")
log("pair gap, one Puiseux branch) at {0} of the {1} couplings:".format(
    len(EP2), len(READS)))
log("  sigma_1 at most {0:.1e}; sigma_2 at least {1:.1e}; the third eigenvalue at least {2:.2e} away; the".format(
    max(r[2] for _, r in EP2), min(r[3] for _, r in EP2), min(r[4] for _, r in EP2)))
log("  half-split / sqrt|delta| below and above equal to {0:.1e} (relative) at delta = 1e-6.".format(
    max(r[5] for _, r in EP2)))
for j, name, r in OTHER:
    if r is None:
        log("  not read at J_bridge = {0:.3f} ({1}): the squared gap does not change sign.".format(j, name))
    else:
        log("  not an EP2 at J_bridge = {0:.3f} ({1}): gamma* = {2:.1e}, pair gap {3:.1e}, sigma_1 = {4:.1e},"
            " sigma_2 = {5:.1e}{6}".format(j, name, r[0], r[1], r[2], r[3],
                                           ", two null directions at zero gain (a zero of the threshold)"
                                           if r[0] < 1e-12 and r[3] < 1e-12 else ""))

log()
log("-" * 78)
log("7. The product gamma_crit * J_bridge at large couplings")
log("-" * 78)
log()
log("The same per-block threshold as section 6 (first unstable point of a 60-point geometric gamma grid")
log("from 1e-4 to 2, refined by 40 bisection steps), at the sweep's last coupling and five beyond it:")
log()
log("  {0:>9} {1:>16} {2:>22}".format("J_bridge", "first block", "gamma_crit * J_bridge"))
for j in (20.0, 30.0, 50.0, 100.0, 300.0, 1000.0):
    best = None
    for name, pieces in SWEEP_PIECES.items():
        for a0, a1, dp in pieces:
            value = piece_threshold(a0 + j * a1, dp)
            if best is None or value < best[0]:
                best = (value, name)
    log("  {0:>9.0f} {1:>16} {2:>22.4f}".format(j, best[1], best[0] * j))
log("These are samples; they establish no limit.")

with open(OUT_PATH, "w", encoding="utf-8") as handle:
    handle.write("\n".join(_LINES) + "\n")
print("\nResults: " + OUT_PATH)
