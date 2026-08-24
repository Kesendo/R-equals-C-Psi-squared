"""The other candidate for a window-stable functional, and what it turned into.

`experiments/MEDIATOR_NOISE_GATE_LEVEL_THREE.md` asked for a window-stable
functional and named two candidates: MI integrated over time, and steady-state
MI.  `experiments/THE_BLIND_SITE.md` section 7 ruled out the first, because the
window mean dilutes as the window grows.  This script tries the second.

THE ANSWER IS NEGATIVE, AND EXACTLY SO.  Whenever dephasing at a SINGLE seat
leaves no blind subspace, the single-excitation Liouvillian has a one-dimensional
kernel and its steady state is I/N to machine precision, at every seat and every
rate.  A support of SEVERAL seats needs connectivity as well: part 9 prints the
cut-chain counterexample (N = 4, [1,0,1], gamma on {0,2}: nothing jointly blind,
kernel 2), and the sweep's own arms, which carry gamma on ten or eleven seats,
have their one-dimensional kernels certified profile by profile by the exact
GF(p) rank instead.  So the steady-state mutual information is the same number for every
profile the sweep compares, and the span it would report is zero.  It is
window-stable because it has stopped looking.  Its value is a closed form in N
alone:

    I(A:B)|_inf  =  log2 N  -  ((N+1)/N) * log2((N+1)/2)

for halves A and B of size (N-1)/2 with the centre in neither.

The repo already held the premise.  PROOF_ASYMPTOTIC_SECTOR_PROJECTION states
that rho(inf) is a function of the initial sector populations alone, a complete
invariant, with no gamma in the conclusion.  An untracked design spec under
docs/superpowers/ makes the same point about the fixed point being profile
independent; it is not quoted here, because a tracked artifact should not rest
on words a clone cannot read.  What was not written down anywhere is
the consequence for a functional: a steady-state observable cannot rank seats
that all carry gamma > 0.  That is this script's first part.

WHAT THE DIMENSION SAYS INSTEAD.  The steady STATE is blind; the stationary
MANIFOLD is not.  In the single-excitation sector, with dephasing at one seat
alone,

    dim ker L_SE(seat j)  =  1  +  (gcd(2j+1, N) - 1) / 2

which is the committed divisor law as an integer, and which holds under the
wider criterion below for any bond profile WITH NO ZERO BOND.  No window, no
functional, no propagation, and no eigensolver and nothing to tolerate on the
kernel side: it is an exact GF(p) rank on integer inputs, taken at two
independent primes in part 8.  Verified N = 3..13 at every seat.  One number in
part 8 does come off a float eigensolver and says so where it stands; it
certifies nothing, and what it exhibits is a BASIS effect rather than a mechanism
failing, which part 8(a3) prints both ways.  Neither gamma's magnitude nor a uniform
J enters the exact statement at all; what does enter is gamma's SUPPORT and the
couplings' profile.

THE XY CHAIN, which was section 11's first open item.  That page predicted the
node condition becomes m(j+1) = 0 (mod N+1) and that "a different divisor answer
is expected".  It is

    XY:          dim = gcd(j+1, N+1) - 1
    Heisenberg:  dim = (gcd(2j+1, N) - 1) / 2

with no halving on the XY side.  Three routes agree, and only the first two are
about the dynamics: the exact kernel dimension, the node count off the
eigenvectors, and an integer enumeration of the node condition, which is
arithmetic rather than physics and certifies the COUNT rather than the kernel.  The two laws disagree loudly: at N = 6 and
N = 12 Heisenberg has blind seats and XY has none; at N = 11 Heisenberg has one
and XY has seven.  At the reflection-fixed seat of an odd chain they agree, both
giving (N-1)/2, which is why section 7's statement that the exact blindness
there does not need Heisenberg over XY is right for a reason that now follows
from both formulas rather than from an argument.

WHAT THIS DOES NOT CLOSE.  F4 (`docs/ANALYTICAL_FORMULAS.md`) reports that a
single dephased seat gives full-space kernel 4 at the end of an N = 3 chain and
6 at its middle, and says in as many words: "Which seat carries the gamma
decides, and why is an open question rather than a formula here."  The last part
here reproduces both numbers and takes one step into that question without
closing it.  On the PATH the full-space kernel is entirely block-diagonal in
popcount, cross-sector weight 0.000000000 against a diagonal weight of exactly
12 at N = 5.  On the RING it is not: at N = 4 with gamma on seat 0 the kernel
has dimension 13 and carries cross-sector weight exactly 4, so "no cross-sector
term" is a fact about the path and not about single-seat dephasing.  An earlier
version of this script asserted the general form and could not have caught its
own error, because it attributed each kernel vector to its DOMINANT block, which
can never report a mixed one.  On the path, each block's contribution is the
dimension of the commutant of H restricted to that sector, intersected with the
operators the dissipator does not kill.  In the single-excitation sector H is nondegenerate and
that dimension is 1 + blind, which is the law above.  For popcount >= 2 the count
is LARGER than 1 + blind: at N = 5 with gamma on the centre the six sectors give
1, 3, 2, 2, 3, 1 = 12 where 1 + blind would give 1, 3, 1, 1, 3, 1 = 10.  The
reason is NOT degeneracy of the sector Hamiltonian: on the uniform chain popcount
2 and 3 are simple, and part 8(a2) checks it exactly, deg gcd(chi, chi') = 0 at
N = 5, 6 and 7 in both.
So F4's question stays open; what is new is its shape and its first sector.

Usage: python simulations/seat_cut_blindness.py
           [steady | kernel | scope | xy | full | criterion | graphs |
            sector | deleted | all]
"""
import itertools
import math
import sys
from fractions import Fraction
from collections import Counter

import numpy as np

sys.path.insert(0, "simulations")
from bridge_sector import bit, mutual_information

TOL = 1e-8


# ----------------------------------------------------------------------
# the single-excitation sector
# ----------------------------------------------------------------------

def se_basis(n):
    states = sorted(1 << (n - 1 - a) for a in range(n))
    return states, {s: i for i, s in enumerate(states)}


def se_hamiltonian(n, bonds, zz=True):
    """H restricted to the single-excitation sector, ZZ optional."""
    states, index = se_basis(n)
    h = np.zeros((n, n), dtype=complex)
    for col, s in enumerate(states):
        for (a, b, j) in bonds:
            ba, bb = bit(s, a, n), bit(s, b, n)
            if zz:
                h[col, col] += j * (1 - 2 * ba) * (1 - 2 * bb)
            if ba != bb:
                flipped = s ^ (1 << (n - 1 - a)) ^ (1 << (n - 1 - b))
                h[index[flipped], col] += 2.0 * j
    return h


def se_liouvillian(n, bonds, gammas, zz=True):
    """The full superoperator on the sector, built column by column."""
    states, index = se_basis(n)
    h = se_hamiltonian(n, bonds, zz)
    site = [0] * n
    for a in range(n):
        site[index[1 << (n - 1 - a)]] = a
    mask = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                mask[i, j] = -2.0 * (gammas[site[i]] + gammas[site[j]])
    lio = np.zeros((n * n, n * n), dtype=complex)
    for a in range(n):
        for b in range(n):
            e = np.zeros((n, n), dtype=complex)
            e[a, b] = 1.0
            lio[:, a * n + b] = (-1j * (h @ e - e @ h) + mask * e).reshape(-1)
    return lio


# ----------------------------------------------------------------------
# the exact route: no eigensolver, no tolerance, no floating point
# ----------------------------------------------------------------------

MODP = (1 << 61) - 1


def _rank_modp(rows, ncols, p=MODP):
    """Rank of an integer matrix over GF(p), by elimination."""
    rows = [r[:] for r in rows]
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], p - 2, p)
        rows[r] = [(x * inv) % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] % p:
                f = rows[i][c]
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows):
            break
    return r


def _as_int_coupling(j):
    """J must be an exact integer here; a float J has no GF(p) image."""
    if isinstance(j, int):
        return j
    if isinstance(j, float) and j.is_integer():
        return int(j)
    raise ValueError(
        f"the exact route needs integer couplings, got J = {j!r}.  A "
        "non-integer J is not a limitation of the physics: the kernel is a "
        "commutant and a uniform rescaling of J cannot move it, so run the "
        "integer chain instead.  For a genuinely rational PROFILE, clear the "
        "denominators first (the commutant is scale-invariant, so [1, 2, 3] "
        "and [0.1, 0.2, 0.3] give the same answer)."
    )


def se_hamiltonian_int(n, bonds, zz=True):
    """The same H as `se_hamiltonian`, in INTEGER arithmetic, for integer J.

    INDEXING, and the two conventions differ.  This matrix is SITE-indexed:
    row/column s is site s.  `se_hamiltonian` is built on `se_basis`, whose
    ascending state order makes its row/column k the site n-1-k, so the two
    are related by the site reversal j -> n-1-j and NOT equal entry by entry.
    They agree on every number this script prints only because both divisor
    laws are symmetric under that reversal and the chains here are uniform or
    reflection-symmetric.  An asymmetric-J extension must fix a convention
    first; `node_count` below reverses explicitly for this reason.
    """
    bonds = [(a, b, _as_int_coupling(j)) for (a, b, j) in bonds]
    h = [[0] * n for _ in range(n)]
    for (a, b, j) in bonds:
        h[a][b] += 2 * j
        h[b][a] += 2 * j
    if zz:
        for s in range(n):
            h[s][s] += sum(-j if (s == a or s == b) else j for (a, b, j) in bonds)
    return h


def exact_kernel_dim(n, bonds, seats, zz=True):
    """dim ker L_SE, exactly, with no eigensolver and nothing to tolerate.

    `seats` is gamma's SUPPORT, not a rate list, and the unstated precondition
    is that gamma > 0 on every seat named and gamma = 0 on every seat not named.
    Passing a seat that carries no dephasing gives the wrong answer, not an
    error: at N = 9 with `seats = [4]` this returns 5, which is the kernel of
    the dephased chain, while the honest gamma = 0 Liouvillian has kernel 9.

    WHY THIS IS THE SAME KERNEL.  L(rho) = -i[H, rho] + M o rho, with M real,
    non-positive, and zero exactly on the entries whose dephased-seat bits
    agree.  For Hermitian rho, Tr(rho . (-i)[H, rho]) = 0, so
    0 = Tr(rho . L(rho)) = sum_ij M_ij |rho_ij|^2 <= 0, forcing rho_ij = 0
    wherever M_ij < 0, and then [H, rho] = 0 as well.  L preserves
    Hermiticity, so its kernel is spanned by Hermitian elements and

        ker L = { rho : [H, rho] = 0, and rho_ij = 0 wherever a dephased
                  seat's bit differs between i and j }.

    Neither gamma's MAGNITUDE nor J's appears: gamma enters only through which
    entries M kills, and J only scales H, which does not move a commutant.  So
    the two "independence" claims are not measurements here, they are visible
    in the statement.

    The float route this replaces is not merely less elegant, it is wrong in
    reachable corners: an SVD rank of the same operator returns 21 instead of 6
    at N = 11, J = 1e-5 while reporting a singular-value gap of 5.95e3, and
    101 at gamma = 1e9.  Both are reproduced by `run_kernel` below, so the
    numbers live in the run file and not only here.  MirrorWorld's `Divisor`
    records the same trap: "a floating-point rank silently miscounts it once
    the coupling is small and the chain long, where the other eigenvalues crowd
    the root at spacing J^(2d)" (compute/MirrorWorld/Divisor.cs).
    """
    h = se_hamiltonian_int(n, bonds, zz)
    allowed = [(i, j) for i in range(n) for j in range(n)
               if all((i == s) == (j == s) for s in seats)]
    index = {c: k for k, c in enumerate(allowed)}
    rows = []
    for a in range(n):
        for b in range(n):
            row = [0] * len(allowed)
            for (i, j) in allowed:
                v = 0
                if j == b:
                    v += h[a][i]
                if i == a:
                    v -= h[j][b]
                if v:
                    row[index[(i, j)]] = v % MODP
            if any(row):
                rows.append(row)
    return len(allowed) - _rank_modp(rows, len(allowed))


def svd_kernel_dim(matrix, tol=TOL):
    """The float route this script does NOT use, kept to show what it does.

    It is here for one purpose: `run_kernel` calls it on two settings where it
    is demonstrably wrong, so the two numbers that justify the exact route are
    produced by the run rather than asserted in prose.  Nothing else calls it.
    """
    sv = np.sort(np.linalg.svd(matrix, compute_uv=False))
    k = int((sv < tol * max(1.0, sv[-1])).sum())
    gap = (sv[k] / sv[k - 1]) if 0 < k < len(sv) else float("inf")
    return k, gap


def chain(n, j=1.0):
    return [(i, i + 1, j) for i in range(n - 1)]


def chain_int(n, j=1):
    """The same chain with integer couplings, for the exact route."""
    return [(i, i + 1, j) for i in range(n - 1)]


def one_seat(n, seat, gamma):
    g = [0.0] * n
    g[seat] = gamma
    return g


def heisenberg_blind(n, j):
    return (math.gcd(2 * j + 1, n) - 1) // 2


def xy_blind(n, j):
    return math.gcd(j + 1, n + 1) - 1


# ----------------------------------------------------------------------
# part 1: the steady state says nothing
# ----------------------------------------------------------------------

def unique_steady(n, bonds, gammas):
    """The stationary state, taken as the last right-singular vector.

    That step is only legitimate when the kernel is ONE-dimensional, and this
    page faults two `simulations/` helpers for taking it without checking.  So
    it is checked here, and exactly: `exact_kernel_dim` certifies the dimension
    over GF(p) on integer inputs before the float solver is allowed to pick a
    single direction.  All 23 profiles the sweep compares pass.
    """
    support = [k for k, g in enumerate(gammas) if g > 0.0]
    dim = exact_kernel_dim(n, bonds, support)
    if dim != 1:
        raise ValueError(
            f"kernel is {dim}-dimensional for support {support}; a single "
            "singular vector is not the steady state there"
        )
    lio = se_liouvillian(n, bonds, gammas)
    _, _, vh = np.linalg.svd(lio)
    rho = vh[-1].conj().reshape(n, n)
    rho = 0.5 * (rho + rho.conj().T)
    return rho / np.trace(rho).real


def run_steady(n=11):
    states, _ = se_basis(n)
    m = (n - 1) // 2
    a_sites, b_sites = list(range(m)), list(range(m + 1, n))
    print(f"STEADY STATE: N = {n}, uniform chain J = 1, halves A = {a_sites},")
    print(f"B = {b_sites}, the centre {m} in neither.")
    print()
    print("The sweep of THE_BLIND_SITE section 7 moves one seat from gamma = 0")
    print("to 0.5 with a baseline of 0.05 on the rest.  Both arms, every seat:")
    print()
    uniform = np.eye(n, dtype=complex) / n
    profiles = [("uniform 0.05", [0.05] * n)]
    for j in range(n):
        profiles.append((f"seat {j} at 0", [0.0 if k == j else 0.05
                                            for k in range(n)]))
        profiles.append((f"seat {j} at 0.5", [0.5 if k == j else 0.05
                                              for k in range(n)]))
    mis, worst = [], 0.0
    for _, g in profiles:
        rho = unique_steady(n, chain(n), g)
        worst = max(worst, float(np.abs(rho - uniform).max()))
        mis.append(mutual_information(rho, n, states, a_sites, b_sites))
    print(f"  profiles compared                     {len(profiles)}")
    print(f"  worst |rho_inf - I/N| over all of them {worst:.3e}")
    print(f"  I(A:B) at rho_inf, spread              {max(mis) - min(mis):.3e}")
    print(f"  the one value                          {mis[0]:.15f}")
    closed = math.log2(n) - ((n + 1) / n) * math.log2((n + 1) / 2)
    print(f"  log2 N - ((N+1)/N) log2((N+1)/2)       {closed:.15f}")
    print()
    print("the span the sweep would report from this functional:")
    for j in range(n):
        lo = unique_steady(n, chain(n), [0.0 if k == j else 0.05
                                         for k in range(n)])
        hi = unique_steady(n, chain(n), [0.5 if k == j else 0.05
                                         for k in range(n)])
        a = mutual_information(lo, n, states, a_sites, b_sites)
        b = mutual_information(hi, n, states, a_sites, b_sites)
        print(f"  seat {j:>2}: {100.0 * (a - b) / a:+.3e} %")
    print()
    print("Zero at every seat.  The functional is window-stable because it has")
    print("stopped looking: the destination is fixed by the conserved sector")
    print("populations, and no rate appears in it.  The repo's premise for this")
    print("is PROOF_ASYMPTOTIC_SECTOR_PROJECTION, Consequence 2.")
    print()
    print("the closed form, checked across N:")
    print(f"{'N':>4} {'measured on I/N':>20} {'closed form':>20} {'diff':>11}")
    for m2 in (5, 7, 9, 11, 13, 15):
        st, _ = se_basis(m2)
        k = (m2 - 1) // 2
        got = mutual_information(np.eye(m2, dtype=complex) / m2, m2, st,
                                 list(range(k)), list(range(k + 1, m2)))
        cf = math.log2(m2) - ((m2 + 1) / m2) * math.log2((m2 + 1) / 2)
        print(f"{m2:>4} {got:20.15f} {cf:20.15f} {got - cf:11.2e}")


# ----------------------------------------------------------------------
# part 2: the dimension says everything
# ----------------------------------------------------------------------

def run_kernel():
    print("KERNEL DIMENSION, computed exactly: dephasing at ONE seat only.")
    print("Claim: dim ker L_SE = 1 + (gcd(2j+1, N) - 1)/2, the committed")
    print("divisor law as an integer.")
    print()
    print("Integer arithmetic and a GF(p) rank for the table below; see the note")
    print("on `exact_kernel_dim` for why this is the same kernel, and for the two")
    print("float answers it replaces.  Nothing in the table is tolerated.  The")
    print("SVD demonstration and the Q table further down are the exceptions and")
    print("are marked where they stand.")
    print()
    print(f"{'N':>3} {'seat':>5} {'blind(j)':>9} {'dim ker':>8} {'ker - 1':>8} "
          f"{'match':>7}")
    bad = 0
    for n in range(3, 14):
        for j in range(n):
            blind = heisenberg_blind(n, j)
            k = exact_kernel_dim(n, chain_int(n), [j])
            ok = (k - 1) == blind
            bad += not ok
            if blind > 0 or j == 0:
                print(f"{n:>3} {j:>5} {blind:>9} {k:>8} {k - 1:>8} "
                      f"{'yes' if ok else 'NO':>7}")
    print()
    print(f"mismatches over N = 3..13, every seat: {bad}")
    print()
    print("WHY THE EXACT ROUTE, measured rather than asserted.  The same")
    print("operator through an SVD rank, at N = 11 with gamma on the centre,")
    print("where the exact answer is 6 at every setting below:")
    print()
    print(f"{'setting':<34} {'SVD rank':>9} {'sv gap':>10} {'exact':>6}")
    for label, bonds, gammas in (
            ("J = 1,    gamma = 0.5", chain(11, 1.0), one_seat(11, 5, 0.5)),
            ("J = 1e-5, gamma = 0.5", chain(11, 1e-5), one_seat(11, 5, 0.5)),
            ("J = 1,    gamma = 1e9", chain(11, 1.0), one_seat(11, 5, 1e9))):
        k, gap = svd_kernel_dim(se_liouvillian(11, bonds, gammas))
        print(f"{label:<34} {k:>9} {gap:>10.2e} "
              f"{exact_kernel_dim(11, chain_int(11), [5]):>6}")
    print()
    print("The middle row is the one that matters: a gap of 5.95e3 and the")
    print("wrong answer, so 'the gap is the evidence, not the tolerance' is")
    print("false as a general defence.  An earlier version of this script made")
    print("exactly that argument.  MirrorWorld's Divisor already recorded the")
    print("trap: 'a floating-point rank silently miscounts it once the coupling")
    print("is small and the chain long, where the other eigenvalues crowd the")
    print("root at spacing J^(2d)' (compute/MirrorWorld/Divisor.cs).")
    print()
    print("WHAT THE KERNEL ACTUALLY CONTAINS, because 'blind' is a DIMENSION and")
    print("not a description of a state.  Plenty of stationary states carry")
    print("amplitude AT the seat.  The kernel is spanned by Q, the sum of the")
    print("projectors onto the modes that do NOT vanish at the seat, plus one")
    print("projector per node-mode; blind(j) counts the second group only.")
    print("N = 7, seat 3, gamma = 0.5, residual = max|L rho| entrywise:")
    print()
    n, seat = 7, 3
    bonds = [(i, i + 1, 1.0) for i in range(n - 1)]
    lio = se_liouvillian(n, bonds,
                         {a: (0.5 if a == seat else 0.0) for a in range(n)})
    h = se_hamiltonian(n, bonds)
    states, index = se_basis(n)
    site = [0] * n
    for a in range(n):
        site[index[1 << (n - 1 - a)]] = a
    row = site.index(seat)
    evals, vecs = np.linalg.eigh(h)
    node = [k for k in range(n) if abs(vecs[row, k]) < 1e-9]
    other = [k for k in range(n) if k not in node]
    proj = lambda k: np.outer(vecs[:, k], vecs[:, k])
    big = sum(proj(k) for k in other)

    def show(label, rho):
        r = rho / np.trace(rho)
        res = np.abs(lio @ r.reshape(-1)).max()
        print(f"   {label:<38} seat amplitude {r[row, row].real:>8.5f}   "
              f"residual {res:.1e}")

    show("Q, the non-node sum (unnormalised: 1.0)", big)
    show("Q + one node projector", big + proj(node[0]))
    show("Q + all node projectors  ( = I/N )", big + sum(proj(k) for k in node))
    show("the node projectors alone", sum(proj(k) for k in node))
    print()
    print(f"   Q sits on the seat with weight {big[row, row].real:.4f} and is")
    print("   stationary; so is every combination above.  Only the last row has")
    print("   no amplitude at the seat, and its dimension is what blind(j)")
    print("   counts.  An earlier version of the page said a state survives")
    print("   ONLY by having no amplitude at the seat, which the second and")
    print("   third rows refute.")


def run_scope():
    print("SCOPE of the kernel dimension.")
    print()
    print("(a) gamma does not enter at all, and that is visible rather than")
    print("measured: the exact kernel depends on which entries the dissipator")
    print("kills, i.e. on the SUPPORT of gamma, and not on any rate.  The same")
    print("holds for a uniform J, which only scales H and cannot move a")
    print("commutant.  The float route needed a five-decade sweep here and")
    print("still had corners where it was wrong.")
    print()
    print("(b) the couplings, where there IS something to measure.  Every seat")
    print("of N = 11 under three bond profiles, exactly:")
    print()
    cases = [
        ("uniform J = 1", [1] * 10),
        ("reflection-symmetric, non-uniform", [1, 2, 3, 1, 2, 2, 1, 3, 2, 1]),
        ("asymmetric", [1, 2, 3, 4, 5, 1, 1, 1, 1, 1]),
    ]
    print(f"{'bond profile':<36} " + " ".join(f"{j:>3}" for j in range(11)))
    for label, js in cases:
        bonds = [(i, i + 1, js[i]) for i in range(10)]
        row = [exact_kernel_dim(11, bonds, [j]) - 1 for j in range(11)]
        print(f"{label:<36} " + " ".join(f"{v:>3}" for v in row))
    print()
    print("Read the rows, not a single cell.  Under the uniform chain only the")
    print("centre is blind at N = 11.  A reflection-symmetric but non-uniform")
    print("chain keeps the centre and nothing else, so what the centre needs is")
    print("the reflection and not uniformity.  Under an asymmetric chain the")
    print("centre goes to 0.  Whether the OTHER seats' blindness (which at")
    print("N = 11 is empty anyway) needs the reflection is not decided by this")
    print("table and is not claimed.")


# ----------------------------------------------------------------------
# part 3: the XY chain, section 11's first open item
# ----------------------------------------------------------------------

def node_count(n, seat, zz):
    """Count H-eigenmodes with an exact node at the seat, from eigenvectors.

    `se_hamiltonian` is built on `se_basis`, whose row k is the site n-1-k, so
    the row to read for site `seat` is n-1-seat.  Both laws are symmetric under
    that reversal, so this correction moves no number here; it is made anyway
    because the next non-symmetric chain would be silently wrong.
    """
    h = se_hamiltonian(n, chain(n), zz)
    _, v = np.linalg.eigh(h)
    return int((np.abs(v[n - 1 - seat, :]) < 1e-9).sum())


def run_xy():
    print("THE XY CHAIN: section 11's first open item.")
    print()
    print("That page predicted the node condition becomes m(j+1) = 0 (mod N+1)")
    print("and that a different divisor answer was to be expected.  It is")
    print("gcd(j+1, N+1) - 1, with no halving.  Route 1, the kernel dimension:")
    print()
    print(f"{'N':>3} {'seat':>5} {'Heis ker-1':>11} {'(gcd(2j+1,N)-1)/2':>18} "
          f"{'XY ker-1':>9} {'gcd(j+1,N+1)-1':>15}")
    bad = 0
    for n in (6, 7, 9, 11, 12, 13):
        for j in range(n):
            kh = exact_kernel_dim(n, chain_int(n), [j], zz=True)
            kx = exact_kernel_dim(n, chain_int(n), [j], zz=False)
            ph, px = heisenberg_blind(n, j), xy_blind(n, j)
            bad += (kh - 1 != ph) + (kx - 1 != px)
            if ph or px:
                print(f"{n:>3} {j:>5} {kh - 1:>11} {ph:>18} {kx - 1:>9} "
                      f"{px:>15}")
    print(f"\n  mismatches: {bad}")
    print()
    print("Route 2, counting nodes off the eigenvectors, N = 3..20:")
    bad_h = bad_x = 0
    for n in range(3, 21):
        for j in range(n):
            bad_h += node_count(n, j, True) != heisenberg_blind(n, j)
            bad_x += node_count(n, j, False) != xy_blind(n, j)
    print(f"  Heisenberg mismatches: {bad_h}    XY mismatches: {bad_x}")
    print()
    print("Route 3, EXACT: no floating point anywhere.  The XY modes are")
    print("sin(pi m (j+1)/(N+1)) for m = 1..N, so a node at seat j means")
    print("(N+1) divides m(j+1).  With d = gcd(j+1, N+1) those m are exactly")
    print("the multiples of (N+1)/d, of which d-1 lie in 1..N.  Counted by")
    print("integer arithmetic, N = 2..200, every seat:")
    bad_h = bad_x = 0
    for n in range(2, 201):
        for j in range(n):
            ex_x = sum(1 for m in range(1, n + 1) if (m * (j + 1)) % (n + 1) == 0)
            ex_h = sum(1 for m in range(0, n)
                       if (m * (2 * j + 1)) % (2 * n) == n % (2 * n))
            bad_x += ex_x != xy_blind(n, j)
            bad_h += ex_h != heisenberg_blind(n, j)
    print(f"  Heisenberg mismatches: {bad_h}    XY mismatches: {bad_x}")
    print()
    print("Where the two laws disagree:")
    print(f"{'N':>4} {'Heisenberg blind seats':<36} {'XY blind seats':<36}")
    for n in (6, 7, 9, 11, 12, 13, 15, 16):
        h = [f"{j}:{heisenberg_blind(n, j)}" for j in range(n)
             if heisenberg_blind(n, j)]
        x = [f"{j}:{xy_blind(n, j)}" for j in range(n) if xy_blind(n, j)]
        print(f"{n:>4} {' '.join(h) or 'none':<36} {' '.join(x) or 'none':<36}")
    print()
    print("At the reflection-fixed seat of an odd chain both give (N-1)/2:")
    print("gcd(N, N) = N on one side and gcd((N+1)/2, N+1) = (N+1)/2 on the")
    print("other.  That is section 7's claim, now following from both formulas")
    print("instead of from an argument about which Hamiltonian is used.")
    print()
    print("PARITY-FORCED BLINDNESS, which the XY book carries and the ZZ book")
    print("does not.  On XY the H_SE diagonal is zero, so a principal block of")
    print("odd size m is singular: expanding along its last row gives")
    print("det T_m = -b^2 det T_{m-2}, b the block's own last bond, with")
    print("det T_1 = 0.  At an ODD seat of an")
    print("ODD chain both leftover blocks have odd size, so both share the")
    print("root 0 and blind(j) >= 1 for EVERY zero-free profile, however")
    print("irregular.  Swept with the exact criterion over pseudo-random")
    print("integer profiles carrying BOTH signs:")
    print()
    import random
    rng = random.Random(20260824)
    alphabet = [1, 2, 3, 5, 7, -1, -2, -3, -5, -7]
    for n in (5, 7, 9):
        odd_hits = odd_total = 0
        heis_blind_pairs = xy_even_hits = even_total = 0
        for _ in range(20):
            js = [rng.choice(alphabet) for _ in range(n - 1)]
            for j in range(n):
                bx = blind_by_gcd(n, _path(n, js), j, zz=False)
                bh = blind_by_gcd(n, _path(n, js), j, zz=True)
                heis_blind_pairs += bh >= 1
                if j % 2 == 1:
                    odd_total += 1
                    odd_hits += bx >= 1
                else:
                    even_total += 1
                    xy_even_hits += bx >= 1
        print(f"  N = {n}: XY blind at odd seats {odd_hits}/{odd_total}, "
              f"at even seats {xy_even_hits}/{even_total}; "
              f"Heisenberg blind anywhere {heis_blind_pairs}/{20 * n}")
    print()
    print("So on XY a third kind of blindness exists beside the mirror-forced")
    print("and the met: PARITY-forced, generic rather than exceptional.  The")
    print("sentence that a genuinely irregular chain will usually have no blind")
    print("seat at all is a ZZ-book sentence; on the swept profiles above the")
    print("Heisenberg count says it, and the XY count refuses it at every odd")
    print("seat of every odd chain.")
    print()
    print("WHAT THE BLIND SUBSPACE CAN HIDE.  Every blind MODE of a zero-free")
    print("profile carries nonzero amplitude at both ends (Jacobi: a vanishing")
    print("end drags the whole vector to zero).  A blind STATE need not: the")
    print("corner entry P[0, N-1] of the exact projector onto the blind")
    print("subspace, I minus the projector onto the seat's Krylov space, all")
    print("arithmetic in Fraction with no eigensolver, is")
    print()
    for label, n, seat, zz in (
            ("XY uniform N = 8, seat 2 (blind 2)", 8, 2, False),
            ("XY uniform N = 11, seat 3 (blind 3), N ODD", 11, 3, False),
            ("XY uniform N = 14, seat 4 (blind 4)", 14, 4, False),
            ("XY uniform N = 7, seat 3 (blind 3)", 7, 3, False),
            ("Heisenberg uniform N = 7, seat 3 (blind 3)", 7, 3, True)):
        corner = _blind_projector_corner(n, chain_int(n), seat, zz)
        print(f"  {label:<44} P[0, N-1] = {corner}")
    print()
    print("The zeros are exact.  P/blind is then a stationary state inside the")
    print("blind subspace with NO end-to-end entry at all, so the modes' end")
    print("coherences can cancel in a blind state.  The criterion is m-PARITY,")
    print("not the parity of N: the end product psi_m(0) psi_m(N-1) is")
    print("(-1)^(m+1) sin^2(pi m/(N+1)) up to normalisation, and with")
    print("d = gcd(j+1, N+1) the node modes are the multiples of M = (N+1)/d,")
    print("whose parity alternates with the multiplier exactly when M is odd.")
    print("So the corner is 0 iff the node set carries BOTH m-parities, that")
    print("is iff M is odd and d >= 3, at even and odd N alike; when M is")
    print("even every node mode is even and the products share one sign, which")
    print("is the N = 7 XY control (M = 2).  Swept exactly, every uniform-XY")
    print("(N, seat) with at least two node modes over N = 3..21:")
    print()
    zero_hits = par_pred = both = tested = 0
    for n in range(3, 22):
        for j in range(n):
            d = math.gcd(j + 1, n + 1)
            if d < 3:
                continue
            tested += 1
            cz = _blind_projector_corner(n, chain_int(n), j, False) == 0
            pp = ((n + 1) // d) % 2 == 1
            zero_hits += cz
            par_pred += pp
            both += (cz == pp)
    print(f"  cases {tested}, corner exactly 0 on {zero_hits}, M odd on "
          f"{par_pred}, agreement {both}/{tested}")
    print()
    print("On the uniform HEISENBERG chain the node modes of any one seat all")
    print("share one m-parity, checked by the integer node condition below, so")
    print("every end product shares one sign and no cancellation is available:")
    print()
    mixed = 0
    for n in range(3, 22):
        for j in range(n):
            ms = [m for m in range(n)
                  if (m * (2 * j + 1)) % (2 * n) == n % (2 * n)]
            if len({m % 2 for m in ms}) > 1:
                mixed += 1
    print(f"  Heisenberg N = 3..21, every seat: node sets with mixed m-parity:"
          f" {mixed}")


# ----------------------------------------------------------------------
# part 4: one step into F4's open question, which stays open
# ----------------------------------------------------------------------

def full_hamiltonian(n, bonds, j=1.0, zz=True):
    d = 1 << n
    h = np.zeros((d, d), dtype=complex)
    for s in range(d):
        for (a, c) in bonds:
            ba, bc = (s >> (n - 1 - a)) & 1, (s >> (n - 1 - c)) & 1
            if zz:
                h[s, s] += j * (1 - 2 * ba) * (1 - 2 * bc)
            if ba != bc:
                h[s ^ (1 << (n - 1 - a)) ^ (1 << (n - 1 - c)), s] += 2.0 * j
    return h


def path_bonds(n):
    return [(i, i + 1) for i in range(n - 1)]


def ring_bonds(n):
    return [(i, (i + 1) % n) for i in range(n)]


def full_kernel_basis(n, bonds, seat, gamma=0.5, zz=True):
    """An orthonormal basis of ker L on the full 4^N space."""
    d = 1 << n
    h = full_hamiltonian(n, bonds, zz=zz)
    mask = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            if ((i >> (n - 1 - seat)) & 1) != ((j >> (n - 1 - seat)) & 1):
                mask[i, j] = -2.0 * gamma
    lio = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d):
        for b in range(d):
            e = np.zeros((d, d), dtype=complex)
            e[a, b] = 1.0
            lio[:, a * d + b] = (-1j * (h @ e - e @ h) + mask * e).reshape(-1)
    _, sv, vh = np.linalg.svd(lio)
    k = int((sv < TOL * max(1.0, sv[0])).sum())
    gap = (sv[-k - 1] / sv[-k]) if 0 < k < len(sv) else float("inf")
    return k, [vh[-(r + 1)].conj().reshape(d, d) for r in range(k)], d, gap


def block_weights(vectors, d):
    """Weight of the kernel in each (bra popcount, ket popcount) block.

    Summed over an ORTHONORMAL basis of the kernel, so this is the trace of the
    kernel projector restricted to the block and does not depend on which basis
    the solver happened to return.  The earlier version of this part attributed
    each basis vector to its dominant block, which cannot report a cross-sector
    direction at all: a gate that could not fail.
    """
    pc = [bin(i).count("1") for i in range(d)]
    w = Counter()
    for m in vectors:
        for i in range(d):
            for j in range(d):
                a = abs(m[i, j]) ** 2
                if a > 0.0:
                    w[(pc[i], pc[j])] += a
    return w


def run_full():
    print("F4's OPEN QUESTION, one step in and not closed.")
    print()
    print("docs/ANALYTICAL_FORMULAS.md, F4: 'on the N = 3 open chain at J = 1,")
    print("gamma = 0.5 on the end seat alone already gives kernel 4, while the")
    print("same gamma on the middle seat gives 6.  Which seat carries the gamma")
    print("decides, and why is an open question rather than a formula here.'")
    print()
    print("This part alone uses the float route, on spaces small enough that")
    print("its answers were reproduced independently; the exact route above is")
    print("built for the single-excitation sector and is not extended here.")
    print()
    print(f"{'N':>3} {'seat':>5} {'full dim ker':>13} {'N+1':>5} {'excess':>8} "
          f"{'sv gap':>10}")
    worst = float("inf")
    for n in (3, 4, 5, 6):
        for j in range(n):
            k, _, _, gap = full_kernel_basis(n, path_bonds(n), j)
            worst = min(worst, gap)
            print(f"{n:>3} {j:>5} {k:>13} {n + 1:>5} {k - n - 1:>8} "
                  f"{gap:>10.2e}")
    print()
    print(f"The count is a bare TOL = {TOL:g} here, which is the pattern this")
    print("page faults elsewhere, so the gap it rests on is printed rather than")
    print(f"assumed: the narrowest deciding RATIO over the table is {worst:.2e}.")
    print("A ratio cannot be counted in decades against an absolute TOL, so it")
    print("is printed and not converted.  It is still a threshold, and this is")
    print("not the only one in the script: run_graphs takes the same TOL, and")
    print("three eigh sites decide a node or blind partition at a bare 1e-9.")
    print("That is five; the sixth is the component count, which is a float")
    print("eigensolver by necessity, being a claim about an eigenBASIS, and")
    print("which certifies nothing.  The page enumerates the same six and says")
    print("of FOUR of them that they print no gap at all.")
    print()
    print("F4's own two numbers reproduce: 4 at the end, 6 at the middle.")
    print()
    print("WHERE THE KERNEL LIVES, and this is where an earlier version of this")
    print("page was WRONG.  It claimed there is no cross-sector stationary")
    print("coherence, full stop.  That holds on the PATH and fails on the RING;")
    print("but the ring is not the reason, and part 8(b) measures why: drop the ZZ")
    print("term and the PATH carries it too.  The qualifier is the ZZ term, not")
    print("the topology on its own.")
    print()
    for label, n, bonds, seat in (("N = 5 path, centre", 5, path_bonds(5), 2),
                                  ("N = 4 path, seat 0", 4, path_bonds(4), 0),
                                  ("N = 4 RING, seat 0", 4, ring_bonds(4), 0)):
        k, vecs, d, _ = full_kernel_basis(n, bonds, seat)
        w = block_weights(vecs, d)
        diag = sum(v for (a, b), v in w.items() if a == b)
        cross = sum(v for (a, b), v in w.items() if a != b)
        cells = ", ".join(f"{w[(a, a)]:.6g}" for a in range(n + 1)
                          if (a, a) in w)
        print(f"  {label:<20} dim ker = {k:>3}   diagonal weight {diag:12.9f}   "
              f"cross-sector weight {cross:12.9f}")
        print(f"  {'':<20} per sector: {cells}")
    print()
    print("On the path the cross-sector weight is zero to machine precision and")
    print("the per-sector weights are exact integers.  On the ring it is not:")
    print("the ring's kernel carries genuine cross-sector stationary coherence,")
    print("so 'no cross-sector term' is a fact about the path and NOT about")
    print("single-seat dephasing in general.  The F4 bullet this page annotates")
    print("lists ring among its topologies, so the distinction has to be kept.")
    print()
    print("What survives, on the path: the seat dependence is a sum of")
    print("per-sector contributions, each of them the commutant of that")
    print("sector's Hamiltonian intersected with the operators the dissipator")
    print("does not kill.  In the single-excitation sector H is nondegenerate")
    print("and that is 1 + blind, the law above; the N = 5 centre's 3 is 1 + 2.")
    print("For popcount >= 2 the count is LARGER: the six sectors give")
    print("1, 3, 2, 2, 3, 1 = 12 where 1 + blind would give 10.  The reason is")
    print("NOT that the sector Hamiltonian is degenerate: on the uniform chain")
    print("popcount 2 and 3 are SIMPLE, and part 8(a2) below checks it exactly,")
    print("deg gcd(chi, chi') = 0 at N = 5, 6 and 7 in both.  So the obvious")
    print("sum rule is false and F4's question stays open.")


# ----------------------------------------------------------------------
# part 5: what the divisor law is a special case OF
# ----------------------------------------------------------------------

def _charpoly(mat):
    """Characteristic polynomial coefficients, exactly, low degree first.

    Expanded by the recursion det(xI - A) over the leading principal minors is
    not available for a general matrix, so this uses Faddeev-LeVerrier with
    Fractions: exact for the integer H_SE this script builds, no eigensolver.
    """
    n = len(mat)
    if n == 0:
        return [Fraction(1)]
    a = [[Fraction(x) for x in row] for row in mat]
    coeffs = [Fraction(1)]
    m = [[Fraction(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        for i in range(n):
            m[i][i] += coeffs[0]
        m = [[sum(a[i][t] * m[t][j] for t in range(n)) for j in range(n)]
             for i in range(n)]
        c = -sum(m[i][i] for i in range(n)) / k
        coeffs.insert(0, c)
    return coeffs


def _poly_gcd_degree(f, g):
    """Degree of gcd(f, g) over Q, by the Euclidean algorithm on Fractions."""
    def trim(q):
        while len(q) > 1 and q[-1] == 0:
            q.pop()
        return q

    f, g = trim(list(f)), trim(list(g))
    while not (len(g) == 1 and g[0] == 0):
        if len(f) < len(g):
            f, g = g, f
            continue
        shift = len(f) - len(g)
        factor = f[-1] / g[-1]
        for i, c in enumerate(g):
            f[shift + i] -= factor * c
        trim(f)
        if len(f) < len(g) or (len(f) == 1 and f[0] == 0):
            f, g = g, f
    return len(f) - 1


def blind_by_gcd(n, bonds, seat, zz=True):
    """blind(seat) = deg gcd(charpoly(H_left), charpoly(H_right)).

    H_SE on an open chain is a Jacobi matrix, hence always simple, and an
    eigenvector vanishes at a site exactly when its eigenvalue is shared by
    the two principal submatrices that site leaves behind.  The blind
    count at a seat is therefore the degree of that shared factor.  The
    divisor law (gcd(2j+1, N) - 1)/2 is this quantity EVALUATED ON THE
    UNIFORM CHAIN; it is not the general answer.
    """
    h = se_hamiltonian_int(n, bonds, zz)
    left = [row[:seat] for row in h[:seat]]
    right = [row[seat + 1:] for row in h[seat + 1:]]
    if not left or not right:
        return 0
    return _poly_gcd_degree(_charpoly(left), _charpoly(right))


def _blind_projector_corner(n, bonds, seat, zz):
    """Corner entry P[0, n-1] of the exact projector onto the blind subspace.

    The blind subspace is the orthogonal complement of the Krylov space the
    seat generates (the same object `blind_truth` counts), so
    P = I - K (K^T K)^{-1} K^T over an independent column subset of the
    Krylov matrix K = [e_seat, H e_seat, ...].  All arithmetic in Fraction;
    no eigensolver, nothing to tolerate.
    """
    h = se_hamiltonian_int(n, bonds, zz)
    cols = []
    vec = [Fraction(1) if s == seat else Fraction(0) for s in range(n)]
    for _ in range(n):
        cols.append(vec[:])
        vec = [sum(Fraction(h[a][b]) * vec[b] for b in range(n))
               for a in range(n)]
    basis, echelon = [], []
    for c in cols:
        v = c[:]
        for r in echelon:
            piv = next(i for i, x in enumerate(r) if x != 0)
            if v[piv] != 0:
                f = v[piv] / r[piv]
                v = [a - f * b for a, b in zip(v, r)]
        if any(x != 0 for x in v):
            echelon.append(v)
            basis.append(c)
    k = len(basis)
    gram = [[sum(basis[i][t] * basis[j][t] for t in range(n))
             for j in range(k)] for i in range(k)]
    rhs = [basis[i][n - 1] for i in range(k)]
    m = [gram[i][:] + [rhs[i]] for i in range(k)]
    for col in range(k):
        piv = next(r for r in range(col, k) if m[r][col] != 0)
        m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]
        m[col] = [x / pv for x in m[col]]
        for r in range(k):
            if r != col and m[r][col] != 0:
                f = m[r][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[col])]
    x = [m[r][k] for r in range(k)]
    proj = sum(basis[i][0] * x[i] for i in range(k))
    return (Fraction(1) if n == 1 else Fraction(0)) - proj


CRITERION_CASES = [
    ("uniform N = 7", 7, [1] * 6),
    ("N = 7 palindrome [1,2,1,1,2,1]", 7, [1, 2, 1, 1, 2, 1]),
    ("N = 9 palindrome [1,2,3,4,4,3,2,1]", 9, [1, 2, 3, 4, 4, 3, 2, 1]),
    ("N = 5 asymmetric [1,4,2,2]", 5, [1, 4, 2, 2]),
    ("N = 4 asymmetric [1,3,2]", 4, [1, 3, 2]),
]


def run_criterion():
    print("WHAT THE DIVISOR LAW IS A SPECIAL CASE OF.")
    print()
    print("An earlier version of this script and of its page said the law needs")
    print("REFLECTION-SYMMETRIC couplings.  That is false in both directions,")
    print("and the three rows of the `scope` part could not see it, because at")
    print("N = 11 only one of the eleven cells is ever non-zero.  Below, every")
    print("seat of five bond profiles, against the law:")
    print()
    for label, n, js in CRITERION_CASES:
        bonds = [(i, i + 1, js[i]) for i in range(n - 1)]
        got = [exact_kernel_dim(n, bonds, [j]) - 1 for j in range(n)]
        law = [heisenberg_blind(n, j) for j in range(n)]
        crit = [blind_by_gcd(n, bonds, j) for j in range(n)]
        print(f"  {label}")
        print(f"    measured  {got}")
        print(f"    law       {law}")
        print(f"    criterion {crit}")
    print()
    print("The N = 7 palindrome is reflection-symmetric and has MORE blindness")
    print("than the law, at seats 1 and 5.  The N = 9 palindrome has LESS, at")
    print("seats 1 and 7.  And the N = 5 asymmetric chain reaches the law's own")
    print("value 2 at the centre with no reflection anywhere, so the reflection")
    print("is not necessary even there.  The law is a UNIFORM-chain law.")
    print()
    print("THE CRITERION, which reproduces all of it.  H_SE on an open chain is")
    print("a Jacobi matrix and therefore simple; an eigenvector vanishes at a")
    print("site exactly when its eigenvalue is shared by the two principal")
    print("submatrices that site leaves behind.  So")
    print()
    print("    blind(j) = deg gcd( charpoly(H_left), charpoly(H_right) )")
    print()
    print("computed here in exact Fractions, no eigensolver.  Swept against the")
    print("exact GF(p) kernel over uniform, ramp, palindromic and pseudo-random")
    print("integer profiles:")
    print()
    seed = 12345
    # BOTH BOOKS.  An earlier version swept zz=True only, and all four surfaces
    # that carry this count stated it with no book named, while part 4 hands the
    # criterion to XY; that is the same defect this script's own fence paragraph
    # was rewritten to remove.  The sweep prints one line per book now, and the
    # page, the arc and this file all say "on each book".
    stats = {True: [0, 0], False: [0, 0]}
    for n in range(3, 9):
        profiles = [[1] * (n - 1), list(range(1, n)),
                    [((i * 7) % 5) + 1 for i in range(n - 1)]]
        # A PALINDROME ON THE BONDS, which is n-1 long, not n.  An earlier
        # version built it n+1 long and truncated, so not one of the six was a
        # palindrome and four surfaces said otherwise; the assertion below is
        # what keeps that from coming back silently.
        m = n - 1
        half = [((i * 3) % 4) + 1 for i in range(m // 2)]
        pal = half + ([3] if m % 2 else []) + half[::-1]
        assert len(pal) == m and pal == pal[::-1], (n, pal)
        profiles.append(pal)
        for _ in range(6):
            seed = (1103515245 * seed + 12345) % (1 << 31)
            profiles.append([((seed >> (3 * i)) % 5) + 1 for i in range(n - 1)])
        for js in profiles:
            bonds = [(i, i + 1, js[i]) for i in range(n - 1)]
            for j in range(n):
                for zz in (True, False):
                    stats[zz][0] += 1
                    stats[zz][1] += (exact_kernel_dim(n, bonds, [j], zz=zz) - 1
                                     != blind_by_gcd(n, bonds, j, zz=zz))
    for zz, name in ((True, "Heisenberg"), (False, "XY        ")):
        t, b = stats[zz]
        print(f"  {name}   (profile, seat) pairs tested: {t} over N = 3..8"
              f"    mismatches: {b}")
    print()
    print("Two corollaries, read off the criterion rather than measured.  On")
    print("the UNIFORM chain the shared factor is the cosine coincidence the")
    print("divisor law counts, which is why that law is exactly right there.")
    print("And a chain reflection-symmetric ABOUT SEAT j makes the two blocks")
    print("equal, so the gcd is the whole characteristic polynomial and")
    print("blind(j) = j.  That is sufficient at the ONE seat it speaks about")
    print("and says nothing about any other, so global reflection symmetry buys")
    print("the centre of an odd chain and nothing else.")
    print()
    print("THE ZERO-BOND FENCE, swept exhaustively rather than exhibited.  Every")
    print("profile in {0,1,2}^(N-1) for N = 3..6, every seat, criterion against")
    print("the exact kernel, split by whether the profile has a zero bond:")
    print()
    print(f"{'book':>11} {'N':>3} {'zero-bond pairs':>16} {'of those wrong':>15} "
          f"{'zero-free pairs':>16} {'of those wrong':>15}")
    grand = {}
    # The simplicity reading below used to be ASSERTED.  Count it: the four
    # cells of (spectrum simple?) x (criterion right?) over the zero-bond half.
    simple_cells = {}
    for zz in (True, False):
        tot = [0, 0, 0, 0]
        cells = {(True, True): 0, (True, False): 0,
                 (False, True): 0, (False, False): 0}
        for n in (3, 4, 5, 6):
            zb = zbw = zf = zfw = 0
            for combo in itertools.product((0, 1, 2), repeat=n - 1):
                bonds = [(i, i + 1, combo[i]) for i in range(n - 1)]
                has_zero = 0 in combo
                simple = _spectrum_is_simple(n, bonds, zz) if has_zero else None
                for seat in range(n):
                    got = exact_kernel_dim(n, bonds, [seat], zz) - 1
                    want = blind_by_gcd(n, bonds, seat, zz)
                    if has_zero:
                        zb += 1
                        zbw += (got != want)
                        cells[(simple, got == want)] += 1
                    else:
                        zf += 1
                        zfw += (got != want)
            tot = [tot[0] + zb, tot[1] + zbw, tot[2] + zf, tot[3] + zfw]
            print(f"{'Heisenberg' if zz else 'XY':>11} {n:>3} {zb:>16} "
                  f"{zbw:>15} {zf:>16} {zfw:>15}")
        print(f"{'':>11} {'all':>3} {tot[0]:>16} {tot[1]:>15} {tot[2]:>16} "
              f"{tot[3]:>15}")
        simple_cells[zz] = cells
        grand['Heisenberg' if zz else 'XY'] = tot
        print()
    h, y = grand['Heisenberg'], grand['XY']
    print("WHAT THE FENCE GUARDS IS BOOK-SPECIFIC: the condition, no zero")
    print("bond, is the same on both books; what a zero bond DOES is not,")
    print("which is why both books are swept.")
    print(f"On HEISENBERG the failure is TOTAL, not partial: all {h[0]} zero-bond")
    print(f"pairs are wrong and all {h[2]} zero-free pairs are right.  On XY the")
    print(f"same zero-bond set has {y[0] - y[1]} pairs the criterion still gets")
    print(f"RIGHT, while the {y[2]} zero-free pairs are right on both books.")
    print()
    print("The reason is that the fence is not really about the zero bond; it is")
    print("about the DEGENERACY a zero bond forces.  On Heisenberg it always")
    print("forces one: each component contributes the one-magnon descendant of")
    print("its own ferromagnetic vacuum at the same eigenvalue, so two")
    print("components always repeat a level.  On XY there is no such term and a")
    print("cut chain can keep a simple spectrum.  Simplicity is NECESSARY there")
    print("and not sufficient, and that is COUNTED rather than asserted; the")
    print("four cells of the XY zero-bond half, and the Heisenberg ones beside")
    print("them because a column of zeros is the whole Heisenberg story:")
    print()
    print(f"{'book':>11} {'simple & right':>15} {'simple & wrong':>15} "
          f"{'degen & right':>14} {'degen & wrong':>14}")
    for zz, name in ((True, "Heisenberg"), (False, "XY")):
        c = simple_cells[zz]
        print(f"{name:>11} {c[(True, True)]:>15} {c[(True, False)]:>15} "
              f"{c[(False, True)]:>14} {c[(False, False)]:>14}")
    print()
    print("Read the XY row: the 'degenerate and right' cell is EMPTY, so every")
    print("pair the criterion gets right is simple, which is necessity; the")
    print("'simple and wrong' cell is not, which is why it is not sufficient.")
    print("On Heisenberg the two simple cells are BOTH empty, which is the")
    print("no-simple-spectrum claim above, counted.")
    print()
    print("The smallest case, seat by seat, so that 'the end seat' is not used")
    print("as a name for it:")
    print()
    _cut = [(0, 1, 0), (1, 2, 1)]
    for zz, name in ((True, "Heisenberg"), (False, "XY        ")):
        print(f"   N = 3, bonds [0, 1], {name}:  measured "
              f"{[exact_kernel_dim(3, _cut, [s], zz=zz) - 1 for s in range(3)]}"
              f"   criterion "
              f"{[blind_by_gcd(3, _cut, s, zz=zz) for s in range(3)]}")
    print()
    print("   On Heisenberg the two ends do not even agree with each other, so")
    print("   the cut is what the criterion cannot see, not the seat's position.")
    print("   The XY row is the smallest of the 60: at SEAT 1 the cut leaves two")
    print("   1x1 blocks with the same characteristic polynomial, deg gcd = 1,")
    print("   and that is the true count, so a zero bond leaves the criterion")
    print("   standing there.  It is the same arithmetic on both books; what")
    print("   differs is the diagonal the ZZ term adds.")
    print()
    print("WHAT THE DEPHASING ACTUALLY KILLS, because an earlier version of")
    print("this work said it 'kills every entry of rho that crosses the seat'.")
    print("It does not.  Z-dephasing at seat j kills rho_ab exactly when the")
    print("SEAT'S OWN bit differs between a and b, which in this sector means")
    print("exactly one of a, b IS the seat: it clears the seat's row and column")
    print("and touches nothing else.  The decay mask at N = 7, seat 3, in units")
    print("of gamma (0 = untouched, -2 = killed), exactly and by construction:")
    print()
    n, seat = 7, 3
    for a in range(n):
        row = [(-2 if (a == seat) != (b == seat) else 0) for b in range(n)]
        print("   " + " ".join(f"{v:>3}" for v in row))
    print()
    print("   entry (0, 6), which CROSSES the seat:   untouched")
    print("   entry (2, 3), which TOUCHES the seat:   killed")
    print()
    print("So the seat does not decohere the two sides from each other, and a")
    print("surviving state is free to be coherent straight across it.  It is:")
    evals, vecs = np.linalg.eigh(np.array(se_hamiltonian_int(n, chain_int(n)),
                                          dtype=float))
    blind = [k for k in range(n) if abs(vecs[seat, k]) < 1e-9]
    print(f"   the {len(blind)} blind modes at N = 7, seat 3, each as rho = "
          f"|psi><psi|, sorted by energy:")
    for k in blind:
        psi = vecs[:, k]
        print(f"     seat amplitude {abs(psi[seat]):.2e}   "
              f"rho[0, 6] = {psi[0] * psi[6]:+.4f}")
    print("   Every one of them is non-zero across the full length of the chain,")
    print("   and that is a THEOREM rather than these three numbers: on a Jacobi")
    print("   matrix no eigenvector can vanish at an END, since the three-term")
    print("   recursion would then drag the whole vector to zero.  So every blind")
    print("   mode of every zero-free profile carries a non-zero end-to-end")
    print("   coherence.  H_SE is simple here, so each value is fixed up to a sign")
    print("   the product cancels.  NOTE THE SCOPE: this is about the BLIND modes,")
    print("   not about every surviving state.  I/N survives too and its")
    print("   end-to-end entry is exactly zero.  One wave with a node at the")
    print("   watched place, not two worlds that cannot see each other.")
    print()
    print("AND WHAT THE TWO HALVES ARE NOT.  It is tempting to read the")
    print("criterion as 'the two pieces can each ring at that pitch ON THEIR")
    print("OWN', and for the XY chain that is exactly what it says.  For")
    print("Heisenberg it is FALSE: H_left and H_right are PRINCIPAL SUBMATRICES,")
    print("and the ZZ term puts a boundary term on the cut site and a different")
    print("shift on each side, so neither block is the free-standing subchain's")
    print("own Hamiltonian.  Below, 'isolated' rebuilds each half as a chain in")
    print("its own right and takes the same gcd:")
    print()
    print(f"{'chain':<18} {'criterion (principal)':<30} "
          f"{'isolated subchains':<30} {'same':>5}")
    for n, zz, label in ((5, True, "Heisenberg N = 5"),
                         (9, True, "Heisenberg N = 9"),
                         (11, True, "Heisenberg N = 11"),
                         (5, False, "XY N = 5"),
                         (11, False, "XY N = 11")):
        bonds = chain_int(n)
        principal = [blind_by_gcd(n, bonds, j, zz=zz) for j in range(n)]
        isolated = [_blind_isolated_halves(n, bonds, j, zz=zz)
                    for j in range(n)]
        print(f"{label:<18} {str(principal):<30} {str(isolated):<30} "
              f"{('yes' if principal == isolated else 'NO'):>5}")
    print()
    print("The two agree at every seat of the XY chain and disagree at every")
    print("Heisenberg row here; at N = 11 the 'on their own' reading invents two")
    print("blind seats where the chain has none.  So that phrasing belongs to")
    print("the XY case, and the general statement is about the two principal")
    print("submatrices the seat leaves behind.")


# ----------------------------------------------------------------------
# part 6: which graphs carry a cross-sector stationary coherence
# ----------------------------------------------------------------------


def shared_neighbourhood_pairs(n, bonds):
    """Non-adjacent vertex pairs with identical neighbourhoods.

    Graph theory calls these FALSE TWINS.  The phrase "twin pair" is not used
    for them anywhere in this repository, because it is already spent on the
    sigma_T pair of eigenvalue strands in the Galois and monodromy work.
    """
    nbr = {v: set() for v in range(n)}
    for (a, b) in bonds:
        nbr[a].add(b)
        nbr[b].add(a)
    return [(u, v) for u in range(n) for v in range(u + 1, n)
            if v not in nbr[u] and nbr[u] == nbr[v]]


GRAPHS = [
    ("P4   path 0-1-2-3", 4, [(0, 1), (1, 2), (2, 3)]),
    ("C3   ring 0-1-2", 3, [(0, 1), (1, 2), (2, 0)]),
    ("C4   ring 0-1-2-3", 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    ("C5   ring 0-1-2-3-4", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
    ("C6   ring 0-1-2-3-4-5", 6,
     [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
    ("S4   star hub 0, leaves 1-3", 4, [(0, 1), (0, 2), (0, 3)]),
    ("S5   star hub 0, leaves 1-4", 5, [(0, 1), (0, 2), (0, 3), (0, 4)]),
    ("C4+p ring 0-1-2-3, pendant 3-4", 5,
     [(0, 1), (1, 2), (2, 3), (3, 0), (3, 4)]),
    ("K4-e  0-1-2-3 minus edge 2-3", 4,
     [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]),
    ("NO FALSE TWIN  {01,02,03,12,13,24}", 5,
     [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 4)]),
]


def run_graphs():
    print("WHICH GRAPHS CARRY A CROSS-SECTOR STATIONARY COHERENCE.")
    print()
    print("The full-space kernel of part 4 is block-diagonal in popcount on the")
    print("open chain and NOT on the N = 4 ring.  This part asks every seat of")
    print("every graph below, so the reading is off the run and not off a")
    print("sentence.  The vertex numbering is printed with each graph, because")
    print("the seat labels below mean nothing without it.")
    print()
    print("The weight is the trace of the kernel projector restricted to the")
    print("off-diagonal popcount blocks; it is basis-independent.  The float")
    print("route of part 4 is used, on spaces of at most 2^5.")
    print()
    for label, n, bonds in GRAPHS:
        pairs = shared_neighbourhood_pairs(n, bonds)
        row = []
        for seat in range(n):
            k, vecs, d, _ = full_kernel_basis(n, bonds, seat)
            w = block_weights(vecs, d)
            cross = sum(v for (a, b), v in w.items() if a != b)
            row.append((seat, k, cross))
        print(f"{label}")
        print(f"    non-adjacent shared-neighbourhood pairs: "
              f"{pairs if pairs else 'none'}")
        for seat, k, cross in row:
            print(f"    seat {seat}: dim ker = {k:>3}   "
                  f"cross-sector weight {cross:12.9f}")
        print()
    print("THIS TABLE IS A NULL RESULT ABOUT THE PAIRS, and it took three")
    print("tries to say so.  The reading under test was that a non-adjacent")
    print("shared-neighbourhood pair (graph theory's FALSE TWINS) is what")
    print("carries the cross-sector weight.  The table refuses it on both")
    print("sides, and the first two summaries written for it were wrong.")
    print()
    print("NOT SUFFICIENT: the star has three such pairs at N = 4 and six at")
    print("N = 5, and carries no cross-sector weight at any seat.")
    print()
    print("NOT NECESSARY: the last row has no NON-ADJACENT shared-neighbourhood")
    print("pair at all (vertices 0 and 1 are ADJACENT with the same closed")
    print("neighbourhood, which is a different thing and is not what is counted)")
    print("yet it carries 6, 12, 6 at seats 2, 3, 4.  That is the row that ends")
    print("the reading, and it was found by a reviewer sweeping every connected")
    print("graph on five vertices, not by this script's own case list.")
    print()
    print("AND THE SEAT RULE POINTS BOTH WAYS, which is the sharpest of the")
    print("three.  On C4+pendant, one pair, the weight is ZERO at exactly that")
    print("pair's seats.  On K4 minus an edge, also one pair, the weight is")
    print("non-zero at exactly that pair's seats and zero everywhere else.  And")
    print("on C4 every seat lies in a pair and every seat carries 4.  No rule")
    print("in terms of the pairs survives all three.")
    print()
    print("What survives is only the open-chain statement of part 4: on the")
    print("HEISENBERG path the kernel is block-diagonal in popcount, and")
    print("elsewhere it need not be.  Not 'off the path', which part 8(b)")
    print("refutes by dropping the ZZ term and keeping the path.  Which systems")
    print("carry cross-sector stationary coherence is OPEN, it is not a question")
    print("about topology alone, and this table is the record of one reading")
    print("that failed.")


# ----------------------------------------------------------------------
# part 7: the exact per-popcount-sector kernel, and what breaks off it
# ----------------------------------------------------------------------

def sector_states(n, w):
    return [s for s in range(1 << n) if bin(s).count("1") == w]


def sector_hamiltonian_int(n, bonds, w, zz=True):
    """H restricted to the popcount-w sector, in INTEGER arithmetic."""
    states = sector_states(n, w)
    index = {s: i for i, s in enumerate(states)}
    h = [[0] * len(states) for _ in states]
    for col, s in enumerate(states):
        for (a, b, j) in bonds:
            ba = (s >> (n - 1 - a)) & 1
            bb = (s >> (n - 1 - b)) & 1
            if zz:
                h[col][col] += j * (1 - 2 * ba) * (1 - 2 * bb)
            if ba != bb:
                flipped = s ^ (1 << (n - 1 - a)) ^ (1 << (n - 1 - b))
                h[index[flipped]][col] += 2 * j
    return states, h


def exact_sector_kernel_dim(n, bonds, w, seats, zz=True):
    """dim of the kernel's popcount-(w,w) block, exactly, by a GF(p) rank.

    Same argument as `exact_kernel_dim`, one sector up: L preserves each (w, w)
    block, the dissipator kills every entry whose dephased-seat bits differ, and
    on what is left the condition is [H_w, rho] = 0.  No eigensolver, no
    tolerance.  A float rank of the same block miscounts for the reason part 2
    demonstrates.
    """
    bonds = [(a, b, _as_int_coupling(j)) for (a, b, j) in bonds]
    states, h = sector_hamiltonian_int(n, bonds, w, zz)
    m = len(states)

    def bits(s):
        return tuple((s >> (n - 1 - q)) & 1 for q in seats)

    allowed = [(i, j) for i in range(m) for j in range(m)
               if bits(states[i]) == bits(states[j])]
    index = {c: k for k, c in enumerate(allowed)}
    rows = []
    for a in range(m):
        for b in range(m):
            row = [0] * len(allowed)
            for (i, j) in allowed:
                v = 0
                if j == b:
                    v += h[a][i]
                if i == a:
                    v -= h[j][b]
                if v:
                    row[index[(i, j)]] = v % MODP
            if any(row):
                rows.append(row)
    return len(allowed) - _rank_modp(rows, len(allowed))


def _blind_isolated_halves(n, bonds, seat, zz=True):
    """blind(seat) computed from the two halves as FREE-STANDING chains.

    Not the criterion.  This is the reading the criterion is often mistaken
    for, kept so the difference can be measured rather than argued: each side
    of the seat is rebuilt as a chain in its own right, carrying only its own
    bonds, and the same gcd degree is taken.  For the XY chain it coincides
    with the criterion because there is no diagonal to differ in; with the ZZ
    term it does not, and part 5 prints both.
    """
    bonds = [(a, b, _as_int_coupling(j)) for (a, b, j) in bonds]
    couplings = {}
    for (a, b, j) in bonds:
        couplings[(min(a, b), max(a, b))] = j
    left, right = seat, n - 1 - seat
    if left == 0 or right == 0:
        return 0

    def half(size, offset):
        sub = [(i, i + 1, couplings[(offset + i, offset + i + 1)])
               for i in range(size - 1)]
        return se_hamiltonian_int(size, sub, zz=zz) if size > 1 else [[0]]

    hl = half(left, 0)
    hr = half(right, seat + 1)
    return _poly_gcd_degree(_charpoly(hl), _charpoly(hr))


MODP32 = (1 << 31) - 1


def _rank_modp_np(matrix, p=MODP32):
    """Rank over GF(p) for a word-sized prime, vectorised.

    `_rank_modp` above is exact but eliminates in Python big integers, which is
    why the cross-block route below was out of reach when this script only had
    it.  Here p is small enough that every product of two reduced entries fits
    in an int64 (p^2 < 2^62), so the elimination is numpy row arithmetic and
    stays exact.  It is a SECOND prime as well as a faster one: the two routes
    are compared head to head in part (a) below, which is worth more than the
    speed, because a GF(p) rank can only be too small and one prime is a
    certificate with a probability attached rather than a proof.
    """
    a = np.mod(np.asarray(matrix, dtype=np.int64), p)
    nrows, ncols = a.shape
    r = 0
    for c in range(ncols):
        nz = np.nonzero(a[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            a[[r, i]] = a[[i, r]]
        a[r] = np.mod(a[r] * pow(int(a[r, c]), p - 2, p), p)
        col = a[:, c].copy()
        col[r] = 0
        hit = np.nonzero(col)[0]
        if hit.size:
            a[hit] = np.mod(a[hit] - np.outer(col[hit], a[r]), p)
        r += 1
        if r == nrows:
            break
    return r


def exact_block_kernel_dim(n, bonds, w_ket, w_bra, seats, zz=True):
    """dim of the kernel's popcount-(w_ket, w_bra) block, exactly.

    The same two conditions as `exact_sector_kernel_dim`, with the two sectors
    allowed to differ, which is what makes it able to measure CROSS-sector
    stationary coherence rather than assume its absence.  On the block the
    commutator term is the Sylvester map rho -> H_ket rho - rho H_bra, built
    here as kron(H_ket, I) - kron(I, H_bra^T); the dissipator keeps exactly the
    entries whose dephased-seat bits agree.  For w_ket == w_bra it reduces to
    `exact_sector_kernel_dim` and part (a) checks that it does.
    """
    bonds = [(a, b, _as_int_coupling(j)) for (a, b, j) in bonds]
    states_k, h_k = sector_hamiltonian_int(n, bonds, w_ket, zz)
    states_b, h_b = sector_hamiltonian_int(n, bonds, w_bra, zz)
    nk, nb = len(states_k), len(states_b)
    hk = np.array(h_k, dtype=np.int64)
    hb = np.array(h_b, dtype=np.int64)
    op = (np.kron(hk, np.eye(nb, dtype=np.int64))
          - np.kron(np.eye(nk, dtype=np.int64), hb.T))

    def bits(s):
        return tuple((s >> (n - 1 - q)) & 1 for q in seats)

    allowed = [a * nb + b for a in range(nk) for b in range(nb)
               if bits(states_k[a]) == bits(states_b[b])]
    if not allowed:
        return 0
    return len(allowed) - _rank_modp_np(op[:, allowed])


def _sum_m_from_invariants(d, dim_a, dim_ker, blocks):
    """sum m_i, solved from the four Wedderburn invariants by search.

    sum n_i m_i = d, sum n_i^2 = dim_a, sum m_i^2 = dim_ker, over exactly
    `blocks` blocks.  Small integer search; raises if the data do not pin a
    unique answer, so the printed number is never a guess.
    """
    from itertools import product
    hits = set()
    rng = range(1, d + 1)
    for ns in product(rng, repeat=blocks):
        if sum(n * n for n in ns) != dim_a:
            continue
        for ms in product(rng, repeat=blocks):
            if (sum(n * m for n, m in zip(ns, ms)) == d
                    and sum(m * m for m in ms) == dim_ker):
                hits.add(sum(ms))
    if len(hits) != 1:
        raise AssertionError(f"invariants do not pin sum m_i: {sorted(hits)}")
    return hits.pop()


def _algebra_dim(mats, dim, primes=(MODP, MODP32)):
    """Dimension of the unital algebra generated by `mats`, exactly.

    Words are grown until the flattened span stops growing; the span's rank is
    taken over GF(p) at two primes, which must agree.  Integer entries only,
    so no eigensolver and nothing to tolerate.
    """
    ident = [[1 if a == b else 0 for b in range(dim)] for a in range(dim)]

    def mul(x, y):
        return [[sum(x[a][t] * y[t][b] for t in range(dim)) for b in range(dim)]
                for a in range(dim)]

    def flat(x):
        return [x[a][b] for a in range(dim) for b in range(dim)]

    span, frontier = [ident], [ident]
    for _ in range(2 * dim):
        grown = []
        for w in frontier:
            for g in mats:
                p = mul(w, g)
                rows = [flat(x) for x in span] + [flat(p)]
                if all(_rank_modp(rows, dim * dim, q) > len(span) for q in primes):
                    span.append(p)
                    grown.append(p)
        if not grown:
            break
        frontier = grown
    return len(span)


def _centre_dim(mats, dim, primes=(MODP, MODP32)):
    """Number of Wedderburn blocks of the algebra generated by `mats`.

    The centre of the commutant is A intersect A', and for a *-algebra its
    dimension is the block count.  Computed as the nullity of the commutator
    map restricted to the algebra's own span; exact GF(p) rank at two primes.
    """
    ident = [[1 if a == b else 0 for b in range(dim)] for a in range(dim)]

    def mul(x, y):
        return [[sum(x[a][t] * y[t][b] for t in range(dim)) for b in range(dim)]
                for a in range(dim)]

    def flat(x):
        return [x[a][b] for a in range(dim) for b in range(dim)]

    span, frontier = [ident], [ident]
    for _ in range(2 * dim):
        grown = []
        for w in frontier:
            for g in mats:
                p = mul(w, g)
                rows = [flat(x) for x in span] + [flat(p)]
                if all(_rank_modp(rows, dim * dim, q) > len(span) for q in primes):
                    span.append(p)
                    grown.append(p)
        if not grown:
            break
        frontier = grown
    rows = []
    for basis_el in span:
        row = []
        for g in mats:
            c = [[sum(g[a][t] * basis_el[t][b] - basis_el[a][t] * g[t][b]
                      for t in range(dim)) for b in range(dim)]
                 for a in range(dim)]
            row.extend(flat(c))
        rows.append(row)
    dims = {len(span) - _rank_modp(rows, len(rows[0]), q) for q in primes}
    if len(dims) != 1:
        raise AssertionError(f"primes disagree on the centre: {sorted(dims)}")
    return dims.pop()


def _sector_components(n, bonds, w, seat, rebases=0, seed=20260824,
                       adapt=False):
    """Connected components of the sector's eigenmodes under the seat's n_seat.

    THE ONE FLOAT IN THIS PART, and deliberately so: the object is a claim ABOUT
    an eigenbasis, so it cannot be phrased without one.  It is never used to
    certify a kernel dimension; those are the exact ranks beside it.  Returns
    the SET of counts found, the first entry from the solver's own eigenbasis
    and the rest from `rebases` random orthogonal re-bases of each degenerate
    eigenspace.

    WHAT THE RE-BASES CAN AND CANNOT SHOW, and an earlier reading of this
    function got it backwards.  On a DEGENERATE sector the count depends on the
    basis, and the basis the mechanism's own argument entitles you to is the one
    that diagonalises n_seat INSIDE each degenerate eigenspace (`adapt=True`),
    because the argument is that rho is diagonal there.  That basis is
    measure-zero in the eigenspace, so a random re-base misses it with
    probability one and returns the generic, maximally connected value every
    time.  A one-element set from random re-bases is therefore NOT evidence
    that the count is basis-free; on a degenerate sector it is evidence of
    nothing at all.  Compare the two: `adapt=True` against the exact rank.
    """
    bonds = [(a, b, _as_int_coupling(j)) for (a, b, j) in bonds]
    states, h = sector_hamiltonian_int(n, bonds, w)
    evals, vecs = np.linalg.eigh(np.array(h, dtype=float))
    occupied = np.diag([1.0 if (s >> (n - 1 - seat)) & 1 else 0.0
                        for s in states])
    scale = max(1.0, float(np.max(np.abs(evals))))
    rng = np.random.default_rng(seed)
    found = set()
    for attempt in range(rebases + 1):
        v = vecs.copy()
        if adapt and not attempt:
            i = 0
            while i < len(evals):
                k = i
                while k + 1 < len(evals) and evals[k + 1] - evals[i] < 1e-9 * scale:
                    k += 1
                if k > i:
                    blk = v[:, i:k + 1].T @ occupied @ v[:, i:k + 1]
                    _, u = np.linalg.eigh(blk)
                    v[:, i:k + 1] = v[:, i:k + 1] @ u
                i = k + 1
        if attempt:
            i = 0
            while i < len(evals):
                k = i
                while k + 1 < len(evals) and evals[k + 1] - evals[i] < 1e-9 * scale:
                    k += 1
                if k > i:
                    q, _ = np.linalg.qr(rng.standard_normal((k - i + 1,
                                                             k - i + 1)))
                    v[:, i:k + 1] = v[:, i:k + 1] @ q
                i = k + 1
        coupling = v.T @ occupied @ v
        d = len(evals)
        seen = [-1] * d
        count = 0
        for i in range(d):
            if seen[i] >= 0:
                continue
            stack = [i]
            seen[i] = count
            while stack:
                a = stack.pop()
                for b in range(d):
                    if seen[b] < 0 and abs(coupling[a][b]) > 1e-8:
                        seen[b] = count
                        stack.append(b)
            count += 1
        found.add(count)
    return found


def full_kernel_by_block(n, bonds, seats, zz=True):
    """The WHOLE 4^n kernel, exactly, as (diagonal sum, cross sum).

    Every (p, q) block is ranked, so the split is measured rather than assumed.
    No eigensolver and no tolerance; the 4^n SVD of part 4 cannot reach n = 7
    and this does.
    """
    diag = cross = 0
    for p in range(n + 1):
        for q in range(n + 1):
            k = exact_block_kernel_dim(n, bonds, p, q, seats, zz)
            if p == q:
                diag += k
            else:
                cross += k
    return diag, cross


def run_sector():
    print("THE FULL KERNEL BY SECTOR, EXACTLY, and two things it settles.")
    print()
    print("Part 4 took the full-space kernel by an SVD and stopped at N = 6.")
    print("Every (p, q) popcount block is ranked here instead, the diagonal ones")
    print("and the CROSS ones, so the whole 4^N kernel is measured and its split")
    print("is not assumed.  Every number in (a), (b) and (c) is an exact rank")
    print("with no eigensolver and nothing to tolerate; (a3) is the exception")
    print("and says so where it stands.")
    print()
    print("(a) The middle seat (N-1)//2 of the open chain, uniform J.  At ODD N")
    print("that seat is reflection-fixed and blind; at EVEN N it is neither, and")
    print("both are printed because the even rows are what decide WHY the sum")
    print("rule below survives where it does.")
    print()
    print(f"{'N':>3} {'seat':>4} {'b':>2} {'per sector, p = 2^61-1':<34} "
          f"{'per sector, p = 2^31-1':<34} {'cross':>6} {'FULL':>5} "
          f"{'(N+1)+b*b+b':>12} {'3(N-1)':>8}")
    for n in (3, 4, 5, 6, 7):
        seat = (n - 1) // 2
        per61 = [exact_sector_kernel_dim(n, chain_int(n), w, [seat])
                 for w in range(n + 1)]
        per31 = [exact_block_kernel_dim(n, chain_int(n), w, w, [seat])
                 for w in range(n + 1)]
        diag, cross = full_kernel_by_block(n, chain_int(n), [seat])
        b = heisenberg_blind(n, seat)
        print(f"{n:>3} {seat:>4} {b:>2} {str(per61):<34} {str(per31):<34} "
              f"{cross:>6} {diag + cross:>5} {n + 1 + b * b + b:>12} "
              f"{(str(3 * (n - 1)) if n % 2 else '-'):>8}")
    print()
    print("The two middle lists are the same numbers by two independent routes:")
    print("the pure-Python elimination at p = 2^61 - 1 and the numpy one at")
    print("p = 2^31 - 1.  THE DIAGONAL BLOCKS ONLY: the cross column is ranked")
    print("at p = 2^31 - 1 alone, and saying otherwise would overstate it.  Two")
    print("primes matter on the diagonal because a GF(p) rank can only be too")
    print("SMALL, so a kernel dimension read off it can only be too LARGE, and")
    print("one prime is a certificate with a probability attached.  They do not")
    print("matter on the cross blocks, for a better reason: a nullity can only")
    print("be reported too LARGE, so a measured nullity of ZERO forces the")
    print("rational nullity to be zero.  One prime settles an empty block.")
    print()
    print("THE CROSS COLUMN IS THE POINT.  The claim below is about the FULL")
    print("kernel, and a sum over diagonal blocks is not a full kernel unless")
    print("the cross blocks are empty.  They are, measured, on every row above.")
    print("Part 4 could confirm that only through N = 6, by an SVD; this reaches")
    print("N = 7, which is the first N where the two formulas differ at all.")
    print()
    print("The 'b^2+b' column is docs/CAUGHT_ERRORS.md's committed parenthetical")
    print("'kernel excess is b^2+b, not gcd', evaluated.  It agrees at N = 3 and")
    print("N = 5 and is WRONG at N = 7, where it says 20 and the exact answer is")
    print("18.  The even rows do not discriminate: b = 0 there, so the formula")
    print("says N+1 and every sector carries 1.  3(N-1) fits the three ODD rows,")
    print("and is printed only for them; three points are not a law and this is")
    print("not offered as one.  What is settled is that the committed formula is")
    print("false, at the smallest N that can tell.")
    print()
    print("AND THE EVEN ROWS ANSWER THE SUM RULE'S 'WHY'.  The natural guess is")
    print("(N+1) + sum_w blind_w, which holds at N = 3, 4 and 6 and fails at")
    print("N = 5 and 7.  An earlier version of the page said the surviving N")
    print("have no in-between sector.  They do: N = 4 has w = 2 and N = 6 has")
    print("w = 2, 3, 4, and the table shows each of them carrying 1.  What the")
    print("even rows lack is not the sector, it is the BLINDNESS: b = 0 at an")
    print("even chain's middle seat, and where nothing is blind every sector")
    print("carries 1 and any sum rule is trivially right.  The rule fails at")
    print("exactly the N where a blind middle seat and an in-between sector")
    print("exist together, which is N = 5 and N = 7 here.")
    print()
    print("(a2) AND THE MECHANISM IS NOT DEGENERACY.  An earlier version of")
    print("this script said the popcount >= 2 sector Hamiltonians are degenerate")
    print("and that this is where the excess comes from.  They are simple, and")
    print("that is checked EXACTLY rather than by a smallest gap: a matrix has a")
    print("simple spectrum iff its characteristic polynomial is squarefree, iff")
    print("deg gcd(chi, chi') = 0, computed here in Fractions.")
    print()
    print(f"{'N':>3} {'sector w':>9} {'dim':>5} {'deg gcd(chi, chi_prime)':>24} "
          f"{'simple':>7}")
    for n in (5, 6, 7):
        for w in (2, 3):
            states, h = sector_hamiltonian_int(n, chain_int(n), w)
            c = _charpoly(h)
            dc = [c[k] * k for k in range(1, len(c))]
            g = _poly_gcd_degree(c, dc)
            print(f"{n:>3} {w:>9} {len(states):>5} {g:>24} "
                  f"{'yes' if g == 0 else 'NO':>7}")
    print()
    print("Heisenberg degeneracy is ACROSS popcount sectors, in the SU(2)")
    print("multiplets, not inside one.  That is a UNIFORM-chain statement, and")
    print("(a3) below carries reflection-symmetric profiles whose popcount-2 and")
    print("popcount-3 sectors DO have a doubled eigenvalue, so it is the PROFILE")
    print("that puts a degeneracy inside a sector.  On a SIMPLE sector the")
    print("count is the number of connected components of the sector's")
    print("eigenmodes under the dephased seat's occupation operator: rho commutes")
    print("with H_w, so it is diagonal in the eigenbasis, and the seat's")
    print("condition forces one common coefficient along every edge.  In the")
    print("single-excitation sector that operator is a rank-one projector, every")
    print("mode not vanishing at the seat fuses into one component and each")
    print("node-mode stands alone, which is part 2's 1 + blind.  A component is")
    print("NOT a blind state: at N = 5, popcount 2, nothing is blind and the")
    print("count is still 2.")
    print()
    print("(a3) THE HYPOTHESIS IS LOAD-BEARING.  Off a simple sector no single")
    print("eigenbasis need diagonalise every rho in the kernel at once, and the")
    print("component count then depends on the basis: it is a LOWER bound in")
    print("every basis, and the uniform chain above never shows the gap while a")
    print("palindromic profile at the same reflection-fixed centre does.  The")
    print("kernel and the simplicity degree below are exact; the component")
    print("counts are float, and what they exhibit is a basis effect.")
    print()
    print(f"{'chain':<28} {'w':>2} {'deg gcd':>8} {'simple':>7} "
          f"{'exact ker':>10} {'components':>11} {'re-bases':>9} "
          f"{'adapted':>8}")
    for label, n, js, w in (("N = 5 uniform", 5, [1] * 4, 2),
                            ("N = 5 palindrome [4,3,3,4]", 5, [4, 3, 3, 4], 2),
                            ("N = 7 uniform", 7, [1] * 6, 3),
                            ("N = 7 palindrome [4,3,4,4,3,4]", 7,
                             [4, 3, 4, 4, 3, 4], 3)):
        seat = (n - 1) // 2
        bonds = [(i, i + 1, js[i]) for i in range(n - 1)]
        states, h = sector_hamiltonian_int(n, bonds, w)
        c = _charpoly(h)
        dc = [c[k] * k for k in range(1, len(c))]
        g = _poly_gcd_degree(c, dc)
        k = exact_sector_kernel_dim(n, bonds, w, [seat])
        seen = _sector_components(n, bonds, w, seat, rebases=50)
        comp = sorted(seen)[0] if len(seen) == 1 else sorted(seen)
        fitted = sorted(_sector_components(n, bonds, w, seat, adapt=True))[0]
        print(f"{label:<28} {w:>2} {g:>8} {'yes' if g == 0 else 'NO':>7} "
              f"{k:>10} {str(comp):>11} "
              f"{('one value' if len(seen) == 1 else 'VARIES'):>9} "
              f"{fitted:>8}")
    print()
    print("Both palindromes carry a doubled sector eigenvalue, and in the")
    print("SOLVER'S basis the component count comes out below the exact kernel")
    print("there.  Neither is a counterexample to the mechanism, and the last")
    print("two columns are why.")
    print()
    print("The 're-bases' column re-runs the count on 50 random re-bases of each")
    print("degenerate eigenspace and returns one value.  That value carries no")
    print("information about the mechanism.  The count is a LOWER bound in every")
    print("basis, so the question is what its MAXIMUM over eigenbases is, and a")
    print("maximum is attained on a measure-zero set of bases: a random draw")
    print("lands in the generic stratum, where the mode graph is maximally")
    print("connected and the count minimal, with probability one.  The test is")
    print("blind in exactly the direction the question runs.")
    print()
    print("The 'adapted' column fits the basis to n_seat inside each degenerate")
    print("eigenspace and reaches the exact kernel on both palindromes.  That")
    print("recipe is NOT the general prescription: it can fall short of the best")
    print("basis, and it is undefined when n_seat is itself degenerate on a")
    print("degenerate eigenspace.  The general statement is algebraic.  The")
    print("kernel is the commutant of A = <H_w, n_seat>; writing A's Wedderburn")
    print("blocks as (n_i, m_i) with sum n_i m_i = dim, the kernel is sum m_i^2")
    print("while the best achievable component count is sum m_i.  They agree")
    print("exactly when every m_i = 1, i.e. when the commutant is abelian.")
    print()
    print("That settles both rows with no solver at all: a multiplicity 2 costs")
    print("4, so a kernel of 1, 2 or 3 forces every m_i = 1 and the best basis")
    print("then realises the kernel as its component count.  Both exhibited")
    print("kernels are 2.  FOUR is the first kernel dimension at which a")
    print("multiplicity 2 fits, so the first at which the mechanism can fail")
    print("at all; the smallest failure exhibited here sits one above it.  The")
    print("star at")
    print("N = 4 with the HUB dephased has blocks (2,1) and (1,2) in the")
    print("single-excitation sector, so its kernel is 5 while no basis can give")
    print("more than 3 components.  The gap is sum m_i^2 - sum m_i, and it is")
    print("open only for the zero-free open chain, where none is exhibited")
    print("here either way.")
    print()
    print("The star witness, measured rather than asserted.  Four exact numbers")
    print("fix the block structure and a reader can redo the arithmetic:")
    print()
    print(f"{'case':<28} {'d':>3} {'dim A':>6} {'dim ker':>8} {'blocks':>7} "
          f"{'best count':>11}")
    for label, nn, ww in (("star N = 4, hub, w = 1", 4, 1),
                          ("star N = 4, hub, w = 2", 4, 2)):
        sb = [(0, s, 1) for s in range(1, nn)]
        states, hh = sector_hamiltonian_int(nn, sb, ww)
        dd = len(states)
        nop = [[1 if (a == b and (states[a] >> (nn - 1 - 0)) & 1) else 0
                for b in range(dd)] for a in range(dd)]
        da = _algebra_dim([hh, nop], dd)
        kk = exact_sector_kernel_dim(nn, sb, ww, [0])
        bb = _centre_dim([hh, nop], dd)
        best = _sum_m_from_invariants(dd, da, kk, bb)
        print(f"{label:<28} {dd:>3} {da:>6} {kk:>8} {bb:>7} {best:>11}")
    print()
    print("sum n_i m_i = d, sum n_i^2 = dim A, sum m_i^2 = dim ker, and the")
    print("block count b pin the (n_i, m_i) uniquely here: (2,1) and (1,2) at")
    print("w = 1, (2,2) and (2,1) at w = 2.  Both give sum m_i = 3 against a")
    print("kernel of 5, so no eigenbasis reaches the kernel and the undercount")
    print("is basis-free.  'best count' is that sum, solved from the four")
    print("numbers to its left rather than read off any basis.")
    print()
    print("These are reflection-symmetric profiles read at the reflection-fixed")
    print("centre, which is the setting of table (a), so the hypothesis has to")
    print("travel with the sentence.")
    print()
    print("(b) WHAT THE 'OPEN CHAIN' QUALIFIER OF PART 4 REALLY IS.  Part 4 says")
    print("the kernel is block-diagonal in popcount on the open chain.  Drop the")
    print("ZZ term and keep the same open chain, N = 5, uniform J:")
    print()
    print(f"{'seat':>5} {'Heis: sector dims':<26} {'Heis total':>11} "
          f"{'XY: sector dims':<26} {'XY total':>9}")
    for seat in range(5):
        ph = [exact_sector_kernel_dim(5, chain_int(5), w, [seat], zz=True)
              for w in range(6)]
        px = [exact_sector_kernel_dim(5, chain_int(5), w, [seat], zz=False)
              for w in range(6)]
        print(f"{seat:>5} {str(ph):<26} {sum(ph):>11} {str(px):<26} "
              f"{sum(px):>9}")
    print()
    print("Now the same XY chain on the WHOLE 4^N space, so the comparison is")
    print("measured here and not asserted.  Two routes: the exact block ranks,")
    print("and part 4's SVD with its projector-trace weight, which is where this")
    print("counterexample was first seen.")
    print()
    print(f"{'seat':>5} {'exact FULL':>11} {'exact diag':>11} "
          f"{'exact cross':>12} {'SVD full':>9} {'SVD cross':>10}")
    for seat in range(5):
        diag, cross = full_kernel_by_block(5, chain_int(5), [seat], zz=False)
        k, vecs, d, _ = full_kernel_basis(5, path_bonds(5), seat, zz=False)
        w = block_weights(vecs, d)
        svd_cross = sum(v for (a, b), v in w.items() if a != b)
        print(f"{seat:>5} {diag + cross:>11} {diag:>11} {cross:>12} "
              f"{k:>9} {svd_cross:>10.6f}")
    print()
    print("The exact route and the SVD agree on every row, which is worth having")
    print("because the exact one needs no threshold and the SVD is what part 4")
    print("uses everywhere else.  The XY chain carries cross-sector stationary")
    print("coherence at the three interior seats and none at the two ends.")
    print()
    print("So the load-bearing qualifier is not the topology, it is the open")
    print("chain WITH the ZZ term, and part 4's dichotomy needs that word.")
    print()
    print("(c) Further Heisenberg chains, every seat, now with the cross blocks")
    print("measured rather than left open.  An earlier version of this part")
    print("printed diagonal sums only and had to say in its own text that a")
    print("diagonal sum is not evidence of block-diagonality; the cross column")
    print("is what turns that caveat into a measurement.")
    print()
    for label, n, js in (("N = 5 uniform", 5, [1] * 4),
                         ("N = 5 [1,4,2,2]", 5, [1, 4, 2, 2]),
                         ("N = 6 [2,1,1,2,1]", 6, [2, 1, 1, 2, 1]),
                         ("N = 7 uniform", 7, [1] * 6)):
        bonds = [(i, i + 1, js[i]) for i in range(n - 1)]
        diag, cross, full = [], [], []
        for seat in range(n):
            d, c = full_kernel_by_block(n, bonds, [seat])
            diag.append(d)
            cross.append(c)
            full.append(d + c)
        print(f"  {label:<20} full kernel by seat: {full}")
        print(f"  {'':<20} of which cross-sector: {cross}")
    print()
    print("Every cross entry is zero, so on the open Heisenberg chain the")
    print("diagonal sum IS the full kernel, at every seat of every row here and")
    print("not only where part 4's SVD could reach.  The N = 5 uniform row also")
    print("reproduces part 4 exactly, 6, 6, 12, 6, 6, so the two routes agree")
    print("where they overlap.  Part 2's law is NOT re-read off this table,")
    print("which prints no per-sector cells at all; for that, see part 2.")
    print("What is NOT claimed here: that this holds for every N, or off the")
    print("chain, or with the ZZ term dropped, where part (b) has already")
    print("measured that it FAILS.")


# ----------------------------------------------------------------------
# part 9: the criterion without the chain, and without the fence
# ----------------------------------------------------------------------

def blind_by_deleted(n, bonds, seat, zz=True):
    """blind(seat) = deg gcd(charpoly(H), charpoly(H with row/column seat cut)).

    The two-halves form `blind_by_gcd` compares the seat's two leftover blocks
    WITH EACH OTHER, which presumes the seat leaves exactly two of them: that
    is a chain statement, and a zero bond breaks it.  This form compares the
    FULL spectrum against the spectrum of the principal submatrix, which is
    defined on any graph and says nothing about pieces.  On an open chain with
    no zero bond the two agree, since there the cut matrix is the direct sum of
    the halves.

    PRINCIPAL SUBMATRIX, not induced subgraph: with the ZZ term on, cutting the
    row and column keeps the diagonal shift the seat's own bonds put on its
    neighbours, and rebuilding the rest as a free-standing graph does not.
    Part 5 measures those two apart at every Heisenberg row.
    """
    h = se_hamiltonian_int(n, bonds, zz)
    keep = [s for s in range(n) if s != seat]
    cut = [[h[a][b] for b in keep] for a in keep]
    return _poly_gcd_degree(_charpoly(h), _charpoly(cut))


def blind_truth(n, bonds, seat, zz=True, primes=(MODP, MODP32)):
    """blind(seat) from the DEFINITION, with no gcd and no eigensolver.

    The blind space is the largest H-invariant subspace annihilated at the
    seat, whose dimension is n minus the rank of the Krylov matrix
    [e_seat, H e_seat, ..., H^n e_seat].  Exact integer entries, ranked over
    GF(p).  A GF(p) rank can only be too SMALL, so this count can only be too
    LARGE; two primes are taken and both must agree.
    """
    h = se_hamiltonian_int(n, bonds, zz)
    vec = [1 if s == seat else 0 for s in range(n)]
    rows = [vec[:]]
    for _ in range(n):
        vec = [sum(h[a][b] * vec[b] for b in range(n)) for a in range(n)]
        rows.append(vec[:])
    dims = {n - _rank_modp(rows, n, p) for p in primes}
    if len(dims) != 1:
        raise AssertionError(f"primes disagree on blind({seat}): {sorted(dims)}")
    return dims.pop()


def _joint_blind(n, bonds, seats, zz=True, primes=(MODP, MODP32)):
    """Largest H-invariant subspace annihilated at EVERY seat in `seats`.

    n minus the rank of the block Krylov matrix generated by all the seats'
    unit vectors at once.  Exact integers, ranked over GF(p) at two primes,
    which must agree; a GF(p) rank can only be too small, so this can only be
    too large.  The single-seat case is `blind_truth`.
    """
    h = se_hamiltonian_int(n, bonds, zz)
    rows = []
    for seat in seats:
        vec = [1 if s == seat else 0 for s in range(n)]
        rows.append(vec[:])
        for _ in range(n):
            vec = [sum(h[a][b] * vec[b] for b in range(n)) for a in range(n)]
            rows.append(vec[:])
    dims = {n - _rank_modp(rows, n, p) for p in primes}
    if len(dims) != 1:
        raise AssertionError(f"primes disagree on joint blind: {sorted(dims)}")
    return dims.pop()


def _spectrum_is_simple(n, bonds, zz=True):
    """deg gcd(chi, chi') == 0, exactly.  Used only to LABEL the rows."""
    c = _charpoly(se_hamiltonian_int(n, bonds, zz))
    deriv = [c[k] * k for k in range(1, len(c))]
    return _poly_gcd_degree(c, deriv) == 0


def _path(n, js):
    return [(s, s + 1, js[s]) for s in range(n - 1)]


def _ring(n, j=1):
    return [(s, (s + 1) % n, j) for s in range(n)]


def _star(n, j=1):
    return [(0, s, j) for s in range(1, n)]


def _complete(n, j=1):
    return [(a, b, j) for a in range(n) for b in range(a + 1, n)]


DELETED_GRAPHS = [
    ("path N = 5 uniform", 5, _path(5, [1] * 4)),
    ("path N = 6 [1,1,0,1,1]  ZERO BOND", 6, _path(6, [1, 1, 0, 1, 1])),
    ("path N = 6 [2,0,3,0,2]  TWO ZERO BONDS", 6, _path(6, [2, 0, 3, 0, 2])),
    ("path N = 4 [1,0,1]  ZERO BOND", 4, _path(4, [1, 0, 1])),
    ("path N = 5 [1,-2,3,-1] NEGATIVE J", 5, _path(5, [1, -2, 3, -1])),
    ("ring N = 4", 4, _ring(4)),
    ("ring N = 5", 5, _ring(5)),
    ("ring N = 6", 6, _ring(6)),
    ("star N = 5", 5, _star(5)),
    ("star N = 6", 6, _star(6)),
    ("complete K4", 4, _complete(4)),
    ("complete K5", 5, _complete(5)),
    ("Y-tree N = 5", 5, [(0, 1, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1)]),
    ("T-tree N = 6", 6, [(0, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1), (4, 5, 1)]),
    ("triangle + pendant N = 4", 4,
     [(0, 1, 1), (1, 2, 1), (2, 0, 1), (2, 3, 1)]),
    ("two triangles bridged N = 6", 6,
     [(0, 1, 1), (1, 2, 1), (2, 0, 1), (3, 4, 1), (4, 5, 1), (5, 3, 1),
      (2, 3, 1)]),
    ("K4 minus an edge N = 4", 4,
     [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1)]),
    ("weighted ring N = 5 [1,2,3,2,1]", 5,
     [(0, 1, 1), (1, 2, 2), (2, 3, 3), (3, 4, 2), (4, 0, 1)]),
    ("two triangles, DISCONNECTED N = 6", 6,
     [(0, 1, 1), (1, 2, 1), (2, 0, 1), (3, 4, 1), (4, 5, 1), (5, 3, 1)]),
    ("isolated seat N = 4", 4, [(0, 1, 1), (1, 2, 1)]),
]


def run_deleted():
    """Where the chain phrasing stops and the criterion does not."""
    print("THE CRITERION WITHOUT THE CHAIN.  Part 6's criterion compares the")
    print("seat's two leftover blocks WITH EACH OTHER, which presumes the seat")
    print("leaves exactly two.  A zero bond breaks that, which is the fence")
    print("part 6 sweeps; a general graph breaks it differently, since cutting")
    print("one vertex usually leaves something connected and there is no second")
    print("spectrum to compare.")
    print()
    print("The form that survives compares the FULL spectrum against the")
    print("principal submatrix's:")
    print()
    print("    blind(j) = deg gcd( chi(H), chi(H with row and column j cut) )")
    print()
    print("On an open chain with no zero bond the cut matrix IS the direct sum")
    print("of the two halves and the two forms coincide.  Below, both are taken")
    print("against the DEFINITION, n - rank of the Krylov matrix on the seat,")
    print("exactly, at two primes.  No eigensolver anywhere.")
    print()

    def _connected(n, bonds):
        # over the bonds that actually couple: a J = 0 entry is not an edge.
        parent = list(range(n))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for (a, b, j) in bonds:
            if j != 0:
                parent[find(a)] = find(b)
        return len({find(v) for v in range(n)}) == 1

    span_fails = {}
    for zz in (True, False):
        label = "Heisenberg (ZZ on)" if zz else "XY (ZZ off)"
        span_fails[zz] = []
        print(f"--- {label}")
        print(f"{'graph':40s} {'simple':>6} {'blind per seat':>20} "
              f"{'deleted':>8} {'halves':>7} {'1+bl':>6}")
        ok_del = span_rows = 0
        paths = wrong_half = degenerate = 0
        for name, n, bonds in DELETED_GRAPHS:
            truth = [blind_truth(n, bonds, s, zz) for s in range(n)]
            dele = [blind_by_deleted(n, bonds, s, zz) for s in range(n)]
            ok_del += (truth == dele)
            ker = [exact_kernel_dim(n, bonds, [s], zz) for s in range(n)]
            ker = [k[0] if isinstance(k, (list, tuple)) else k for k in ker]
            span_ok = all(k == 1 + b for k, b in zip(ker, truth))
            span_rows += span_ok
            is_path = len(bonds) == n - 1 and \
                sorted(a for (a, b, _) in bonds) == list(range(n - 1)) and \
                all(b == a + 1 for (a, b, _) in bonds)
            hs = "n/a"
            if is_path:
                js = [0] * (n - 1)
                for (a, b, j) in bonds:
                    js[a] = j
                half = [blind_by_gcd(n, _path(n, js), s, zz) for s in range(n)]
                paths += 1
                if truth != half:
                    wrong_half += 1
                hs = "match" if truth == half else "WRONG"
            is_simple = _spectrum_is_simple(n, bonds, zz)
            degenerate += (not is_simple)
            simple = "yes" if is_simple else "NO"
            if not span_ok:
                span_fails[zz].append((name, _connected(n, bonds)))
            print(f"{name:40s} {simple:>6} {str(truth):>20} "
                  f"{'match' if truth == dele else 'WRONG':>8} {hs:>7} "
                  f"{'holds' if span_ok else 'BREAKS':>6}")
        print(f"    deleted form: {ok_del}/{len(DELETED_GRAPHS)} graphs match "
              f"the definition at every seat.")
        print(f"    halves form: defined on {paths} paths, WRONG on "
              f"{wrong_half} of them.")
        print(f"    {degenerate} of the {len(DELETED_GRAPHS)} graphs have a "
              f"DEGENERATE spectrum, so the COUNT carries no simplicity "
              f"hypothesis.")
        print(f"    the kernel identity dim ker L_SE = 1 + blind holds on "
              f"{span_rows}/{len(DELETED_GRAPHS)} graphs: that one, unlike "
              f"the count, DOES carry a hypothesis, named below.")
        print()

    print("Read the table this way.  The deleted form matches the definition on")
    print("every row, the rows whose spectrum is NOT simple included: the two")
    print("triangle graphs, the star, the complete graphs, and every zero-bond")
    print("path, whose halves repeat each other.  The halves form is wrong on")
    print("exactly the zero-bond paths, which is part 6's fence.")
    print()
    print("So for the COUNT the fence belongs to the PHRASING and not to the")
    print("phenomenon.  Stated on the principal submatrix the criterion needs")
    print("no fence, no simplicity hypothesis and no chain; the chain is where")
    print("the cut matrix happens to fall into two pieces, which is what let")
    print("the criterion be written as a gcd of two halves in the first place.")
    print()
    nz = len(span_fails[True])
    nx = len(span_fails[False])
    conn_z = sorted(nm for nm, c in span_fails[True] if c)
    conn_x = sorted(nm for nm, c in span_fails[False] if c)
    assert conn_z == conn_x, (conn_z, conn_x)
    print("The SPAN is a different story, and the last column is what the")
    print("chain's 'load-bearing twice' really covers: a zero bond does break")
    print("the span, but the zero bond is NOT the discriminator.  The counts")
    print("are book-specific and the tables ABOVE carry them:")
    print()
    print(f"    span failures    Heisenberg {nz} of {len(DELETED_GRAPHS)}"
          f"    XY {nx} of {len(DELETED_GRAPHS)}")
    print()
    print("and the honest reading of them is CONNECTIVITY, not the bond list.")
    print("Counting failures 'with no zero bond' turns on whether a missing")
    print(f"edge is written as J = 0 or left out: {len(conn_z)} of the failures")
    print("are on CONNECTED graphs, the SAME ones on both books,")
    print()
    for nm in conn_z:
        print(f"        {nm}")
    print()
    print("while the zero-bond path [1,0,1], which is disconnected, HOLDS.  So")
    print("neither the zero bond nor the disconnection decides it.  What does:")
    print("dim ker L_SE(j) = 1 + blind(j) can BREAK only where the spectrum is")
    print("degenerate, and plenty of degenerate rows hold anyway: the rings, on")
    print("BOTH books; K4 minus an edge, degenerate and holding on the ZZ book")
    print("and simple on the XY one; and the isolated seat, degenerate on both,")
    print("holding on ZZ and BREAKING on XY.  So simplicity is")
    print("SUFFICIENT for the span and NOT necessary.  That is the hypothesis")
    print("the summary lines above point at.  On [1,1,0,1,1] at seats 1 and 4")
    print("the count is 4 and the kernel is 7.  The count generalises and the")
    print("kernel identity does not, so a reader carrying the fence-free")
    print("criterion to the kernel would be carrying it past the hypothesis it")
    print("needs.")
    print()
    js110 = [1, 1, 0, 1, 1]
    per = []
    for s in range(6):
        kk = exact_kernel_dim(6, _path(6, js110), [s], True)
        kk = kk[0] if isinstance(kk, (list, tuple)) else kk
        per.append((s, blind_truth(6, _path(6, js110), s, True), kk))
    print("Per seat on [1,1,0,1,1], ZZ book (blind | dim ker):")
    print("    " + "   ".join(f"seat {s}: {b} | {k}" for s, b, k in per))
    print("Seats 1 and 4 carry the excess; the chain's other four seats give")
    print("a blind count of 3, hence a predicted span of 4, against a kernel")
    print("of 4.")
    print()
    print("ONE SEAT AGAINST A SUPPORT.  For a single seat, blind(j) = 0 forces")
    print("a one-dimensional kernel on ANY graph and needs no simplicity: it")
    print("says e_j is cyclic for H, so anything commuting with H and with the")
    print("seat's projector sends e_j to a multiple of itself and is therefore a")
    print("multiple of the identity.  For a SUPPORT of several seats that fails")
    print("without connectivity, and the smallest witness is a cut chain:")
    print()
    for label, nn, js, sup in (("N = 4 [1,0,1], gamma on {0, 2}", 4, [1, 0, 1], [0, 2]),
                               ("N = 4 [1,0,1], gamma on {0, 3}", 4, [1, 0, 1], [0, 3]),
                               ("N = 4 [1,1,1], gamma on {0, 2}", 4, [1, 1, 1], [0, 2])):
        bb = [(i, i + 1, js[i]) for i in range(nn - 1)]
        jb = _joint_blind(nn, bb, sup)
        kk = exact_kernel_dim(nn, bb, sup)
        kk = kk[0] if isinstance(kk, (list, tuple)) else kk
        print(f"    {label:<32} joint blind = {jb}, dim ker = {kk}")
    print()
    print("The connected chain gives 1 + 0, as the single-seat statement would;")
    print("the two disconnected pairs each keep their own scale and give 2 with")
    print("nothing blind to both seats at all.  So the one-seat form is a")
    print("theorem and the several-seat form carries a connectivity clause.")
    print()
    print("What this does NOT do: it does not decide which systems carry")
    print("cross-sector stationary coherence, which is what the open item asks.")
    print("It removes the typecheck failure that stood in the way of asking the")
    print("question on a graph at all.")

PARTS = [("steady", run_steady), ("kernel", run_kernel),
         ("scope", run_scope), ("xy", run_xy), ("full", run_full),
         ("criterion", run_criterion), ("graphs", run_graphs),
         ("sector", run_sector), ("deleted", run_deleted)]


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    names = [n for n, _ in PARTS]
    if mode not in names + ["all"]:
        print(f"unknown part {mode!r}; known: {' | '.join(names)} | all")
        return
    for i, (name, fn) in enumerate(PARTS):
        if mode in (name, "all"):
            fn()
            if mode == "all" and i < len(PARTS) - 1:
                print()
                print("-" * 72)
                print()


if __name__ == "__main__":
    main()
