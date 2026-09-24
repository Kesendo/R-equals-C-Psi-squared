"""
The reflection keeps past and future apart: the Born shadow's cross term in the reflection sectors
================================================================================================
experiments/BORN_RULE_SHADOW.md splits rho(t) of a Z-dephased chain at the rate Σγ into a slow
part (the past) and a fast part (the rest) and reads the purity cross term 2·Tr(rho_slow·rho_fast).
Books: H = J Σ_bonds (XX + YY + Δ·ZZ) with Pauli matrices and uniform bonds, and the dissipator
γ Σ_l (Z_l rho Z_l − rho), so a coherence |i><j| decays at 2γ times the number of sites where i
and j differ (as in simulations/born_rule_shadow.py). At N = 3 on the open chain, a state with the
chain's reflection symmetry (R rho R = rho, R the site reversal 0 <-> 2, F71's spatial mirror)
carries no cross term at any time exactly when (4 − Δ²)·J² > γ² (for J ≠ 0). This script gates
the exact premises of that statement and of its companions, and prints numerical readings beside
them (experiments/THE_REFLECTION_KEEPS_PAST_AND_FUTURE_APART.md).

Exact part (sympy, symbolic in J, γ and a real Δ; every check is an identity and can fail):
  E1  of the sixteen joint-popcount blocks (p, q) only (1,1), (2,2), (1,2), (2,1) carry more than
      one dephasing rate, so only they can be split by the cut;
  E2  the reflection-even part of (1,1) has characteristic polynomial λ(λ + 4γ)·p3(λ),
      p3 = λ³ + 8γλ² + (16γ² + (32 + 4Δ²)J²)λ + 96γJ²; its root −4γ belongs to
      X_Δ = Δ(|100><001| + |001><100|) + (the one-excitation hopping matrix), which commutes with H,
      and L and L† both send X_Δ to −4γ·X_Δ and the identity of the sector to 0; the premises
      that make the remaining triple irreducible: p3(0) = 96γJ², p3(−4γ) = −16γJ²(Δ² + 2), and
      the one level-0 direction d = |100><100| + |001><001| − 2|010><010| of the triple is moved by
      the commutator;
  E3  the reflection-even part of (1,2) has characteristic polynomial (λ + 2γ)(λ + 6γ)·(−p3(−λ − 6γ)),
      and its modes at −2γ and −6γ are eigenvectors of L and of L† at symbolic J and Δ;
  E4  Routh-Hurwitz on p3 at the cut: the cubic in ν = −(λ + 3γ) has coefficients 1, γ,
      (32 + 4Δ²)J² − 5γ², 3γ(γ² + 4Δ²J²), and Hurwitz determinant 8γ((4 − Δ²)J² − γ²);
  E5  the reflection-odd parts: each has dimension 4 and trace −3γ per mode; (2,2) and (2,1) have
      the characteristic polynomials of (1,1) and (1,2), and (1,2) is the palindromic partner of
      (1,1); the dissipator is diagonal with one level-0 and three −4γ directions on (1,1), one
      −6γ and three −2γ directions on (1,2); for (1,1), the resultants in ω for a root on the
      lines Re λ = −4γ and Re λ = 0 are 2²⁴·Δ⁴J¹²γ⁴ and a product of factors no real J ≠ 0 can
      make zero (a nonzero resultant means no root on the line; the partner relation carries
      this to the lines −2γ and −6γ of (1,2)); in μ = λ + 3γ the (1,1) polynomial has μ³
      coefficient 0 and μ¹ coefficient −8γ(Δ²J² + γ²); at Δ = 0 (the XY chain) the odd (1,1)
      part factors as (λ² + 4γλ + 8J²)((λ + 4γ)² + 8J²), whose second factor sits on the edge,
      and the first two roots are slower than Σγ exactly when 8J² > 3γ²;
  E6  the W law: N·|W><W| − P_1 commutes with the Heisenberg model with unit couplings on the
      chain, ring, star and complete graph and with random integer couplings on the complete
      graph, N = 3, 4, 5, and fails to commute with the XXZ chain at Δ = 2, the control;
  E7  on the N = 3 ring (the triangle) every mixed part, of either reflection sector, factors into
      window-edge modes times one triple q3 = λ³ + 8γλ² + (16γ² + 36J²)λ + 96γJ² or its palindromic
      image, for every Δ.

Numerical readings (printed, not thresholded; the exact part carries the claims): the sector-
resolved non-orthogonality ||P_s^H P_f|| on the chain at Δ = 1 and Δ = 0, on the ring and at N = 4
(with, per table, the sectors that sit wholly on the cut, the largest distance of an on-cut mode from
Σγ and the gap to the nearest mode off it, which is what makes the cut band harmless),
the threshold bracketed at 0.98 and 1.02 of the exact values, the rates of every mixed part at
J/γ = 1.5 and 20, the purity split of |+0+> and |++0> at their folds and of |+0+> below the
threshold, and the W law against expm.

Import-inert; writes only through main(output_path).
Output: simulations/results/born_shadow_reflection_gate.txt
"""

import argparse
import itertools
from pathlib import Path
import sys

import numpy as np
import sympy as sp
from scipy.linalg import expm, schur, solve_sylvester

DEFAULT_OUTPUT = Path(__file__).parent / "results" / "born_shadow_reflection_gate.txt"
J_PAGE, GAMMA_PAGE = 1.0, 0.05          # the page's parameters
CUT_BAND = 1e-9                         # the page's cut rule: modes within this of Σγ go to the fast part
SEED = 20260924

PX = [[0, 1], [1, 0]]
PY = [[0, -1j], [1j, 0]]
PZ = [[1, 0], [0, -1]]
CHAIN3 = [(0, 1), (1, 2)]
RING3 = [(0, 1), (1, 2), (0, 2)]


# ---------------------------------------------------------------- shared combinatorics

def popcount(b):
    return bin(b).count("1")


def reverse_bits(b, n):
    return int(format(b, f"0{n}b")[::-1], 2)


def hamming(a, b):
    return popcount(a ^ b)


# ---------------------------------------------------------------- exact (sympy) builders

def sym_hamiltonian(n, bonds, J, Delta):
    """H = J Σ_bonds (XX + YY + Δ ZZ) as a sympy matrix; site 0 is the leftmost tensor factor."""
    I2 = sp.eye(2)
    mats = {"X": sp.Matrix([[0, 1], [1, 0]]), "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]),
            "Z": sp.Matrix([[1, 0], [0, -1]])}          # exact entries: no float enters the symbolic part
    H = sp.zeros(2**n)
    for i, j in bonds:
        for letter, c in (("X", 1), ("Y", 1), ("Z", Delta)):
            ops = [I2] * n
            ops[i] = mats[letter]
            ops[j] = mats[letter]
            term = ops[0]
            for o in ops[1:]:
                term = sp.kronecker_product(term, o)
            H += J * c * term
    return H


def sym_block(H, n, p, q, gamma):
    """L on the joint-popcount block (p, q): cells |i><j| with popcount(i) = p, popcount(j) = q.
    L(|i><j|) = −i(H|i><j| − |i><j|H) − 2γ·hamming(i, j)·|i><j|."""
    kets = [b for b in range(2**n) if popcount(b) == p]
    bras = [b for b in range(2**n) if popcount(b) == q]
    cells = [(i, j) for i in kets for j in bras]
    pos = {c: k for k, c in enumerate(cells)}
    M = sp.zeros(len(cells))
    for (i, j), col in pos.items():
        for i2 in kets:                                   # −i H|i><j|
            if H[i2, i] != 0:
                M[pos[(i2, j)], col] += -sp.I * H[i2, i]
        for j2 in bras:                                   # +i |i><j| H
            if H[j, j2] != 0:
                M[pos[(i, j2)], col] += sp.I * H[j, j2]
        M[col, col] += -2 * gamma * hamming(i, j)
    return M, cells


def reflection_isometry(cells, n, sign):
    """Orthonormal basis (columns) of the cells' span with R(ρ) = sign·ρ, R(|i><j|) = |Ri><Rj|."""
    pos = {c: k for k, c in enumerate(cells)}
    cols, seen = [], set()
    for c in cells:
        if c in seen:
            continue
        rc = (reverse_bits(c[0], n), reverse_bits(c[1], n))
        seen.update({c, rc})
        v = sp.zeros(len(cells), 1)
        if rc == c:
            if sign > 0:
                v[pos[c]] = 1
                cols.append(v)
            continue
        v[pos[c]] = 1 / sp.sqrt(2)
        v[pos[rc]] = sign / sp.sqrt(2)
        cols.append(v)
    return sp.Matrix.hstack(*cols) if cols else sp.zeros(len(cells), 0)


def cell_vector(cells, entries):
    """Column vector over `cells` from a dict {(i, j): value}."""
    return sp.Matrix([entries.get(c, 0) for c in cells])


def edge_resultant(poly, lam, r):
    """Resultant in ω of the real and imaginary parts of poly(r + iω); a nonzero value means that no
    root has Re λ = r (a zero may also come from two roots mirrored through that line)."""
    w = sp.Symbol("omega", real=True)
    e = sp.expand(poly.subs(lam, r + sp.I * w))
    return sp.factor(sp.resultant(sp.re(e), sp.im(e), w))


# ---------------------------------------------------------------- numeric builders (float)

def num_hamiltonian(n, bonds, J=1.0, Delta=1.0):
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P, c in ((PX, 1.0), (PY, 1.0), (PZ, Delta)):
            ops = [np.eye(2, dtype=complex)] * n
            ops[i] = np.array(P, dtype=complex)
            ops[j] = np.array(P, dtype=complex)
            term = ops[0]
            for o in ops[1:]:
                term = np.kron(term, o)
            H += J * c * term
    return H


def num_liouvillian(H, n, gamma):
    """Row-major vec (rho.ravel()), as in simulations/born_rule_shadow.py."""
    d = 2**n
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(n):
        ops = [np.eye(2, dtype=complex)] * n
        ops[k] = np.array(PZ, dtype=complex)
        Zk = ops[0]
        for o in ops[1:]:
            Zk = np.kron(Zk, o)
        L += gamma * (np.kron(Zk, Zk) - np.kron(Id, Id))
    return L


def spectral_projector(M, select):
    """Riesz projector onto the eigenvalues select() picks, along the others (sorted Schur + Sylvester)."""
    T, Q, k = schur(M, output="complex", sort=select)
    n = T.shape[0]
    if k in (0, n):
        return (np.eye(n, dtype=complex) if k == n else np.zeros((n, n), dtype=complex)), k
    X = solve_sylvester(T[:k, :k], -T[k:, k:], -T[:k, k:])
    S = np.eye(n, dtype=complex); S[:k, k:] = X
    Si = np.eye(n, dtype=complex); Si[:k, k:] = -X
    E = np.zeros((n, n), dtype=complex); E[:k, :k] = np.eye(k)
    return Q @ S @ E @ Si @ Q.conj().T, k


def slow_select(sigma):
    return lambda z: (-z.real) < sigma - CUT_BAND


def num_sector(n, bonds, gamma, J, Delta, p, q, sign, L=None):
    d = 2**n
    if L is None:
        L = num_liouvillian(num_hamiltonian(n, bonds, J, Delta), n, gamma)
    kets = [b for b in range(d) if popcount(b) == p]
    bras = [b for b in range(d) if popcount(b) == q]
    cells = [(i, j) for i in kets for j in bras]
    idx = [i * d + j for i, j in cells]
    W = np.array(reflection_isometry(cells, n, sign).evalf(), dtype=complex)
    if W.shape[1] == 0:
        return None
    return W.conj().T @ L[np.ix_(idx, idx)] @ W


def num_sector_table(n, bonds, gamma, J, Delta, with_cut=False):
    """(p, q, sign, dim, slow count, ||P_s^H P_f||) for every sector that the cut splits; with_cut also
    returns (sectors wholly on the cut, largest on-cut distance from Σγ, gap to the nearest off-cut mode)."""
    L = num_liouvillian(num_hamiltonian(n, bonds, J, Delta), n, gamma)
    sigma = n * gamma
    rows, wholly, on_res, gap = [], [], None, np.inf
    for p in range(n + 1):
        for q in range(n + 1):
            for sign in (+1, -1):
                M = num_sector(n, bonds, gamma, J, Delta, p, q, sign, L)
                if M is None:
                    continue
                dist = np.abs(-np.linalg.eigvals(M).real - sigma)
                on = dist <= CUT_BAND
                if on.any():
                    on_res = max(on_res or 0.0, dist[on].max())
                if (~on).any():
                    gap = min(gap, dist[~on].min())
                if on.all():
                    wholly.append((p, q, sign))
                P, k = spectral_projector(M, slow_select(sigma))
                if 0 < k < M.shape[0]:
                    Pf = np.eye(M.shape[0]) - P
                    rows.append((p, q, sign, M.shape[0], k, np.linalg.norm(P.conj().T @ Pf, 2)))
    return (rows, wholly, on_res, gap) if with_cut else rows


def purity_split(L, P, rho0, t):
    d = rho0.shape[0]
    v = expm(L * t) @ rho0.ravel()
    rp = (P @ v).reshape(d, d); rp = (rp + rp.conj().T) / 2
    rho = v.reshape(d, d); rho = (rho + rho.conj().T) / 2
    rf = rho - rp
    return (np.real(np.trace(rp @ rp)), np.real(np.trace(rf @ rf)), 2 * np.real(np.trace(rp @ rf)))


def cpsi(rho):
    """The page's book: purity times the sum of the off-diagonal moduli over d - 1."""
    d = rho.shape[0]
    return np.real(np.trace(rho @ rho)) * (np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))) / (d - 1)


def first_downward_quarter(L, rho0, t_max=60.0, steps=12000):
    """First time CΨ falls through 1/4 from above, refined by bisection (as in born_rule_shadow.py)."""
    d = rho0.shape[0]
    dt = t_max / steps
    U = expm(L * dt)
    v = rho0.ravel().copy()
    c_prev = cpsi(rho0)
    for k in range(1, steps + 1):
        v = U @ v
        c = cpsi(v.reshape(d, d))
        if c_prev > 0.25 >= c:
            lo, hi = (k - 1) * dt, k * dt
            for _ in range(60):
                mid = (lo + hi) / 2
                if cpsi((expm(L * mid) @ rho0.ravel()).reshape(d, d)) > 0.25:
                    lo = mid
                else:
                    hi = mid
            return hi
        c_prev = c
    return None


def product_state(letters):
    vecs = {"0": np.array([1, 0], dtype=complex), "1": np.array([0, 1], dtype=complex),
            "+": np.array([1, 1], dtype=complex) / np.sqrt(2)}
    psi = vecs[letters[0]]
    for c in letters[1:]:
        psi = np.kron(psi, vecs[c])
    return np.outer(psi, psi.conj())


def reflect_operator(rho, n):
    d = 2**n
    perm = [reverse_bits(b, n) for b in range(d)]
    return rho[np.ix_(perm, perm)]


def integer_xxz(m, weighted_bonds, Delta=1):
    """Integer matrix of Σ w·(XX + YY + Δ ZZ) over (i, j, w): XX + YY hops with amplitude 2, ZZ is ±1."""
    d = 2**m
    Hi = np.zeros((d, d), dtype=np.int64)
    for i, j, w in weighted_bonds:
        for b in range(d):
            bi, bj = (b >> (m - 1 - i)) & 1, (b >> (m - 1 - j)) & 1
            if bi == bj:
                Hi[b, b] += w * Delta
            else:
                Hi[b, b] -= w * Delta
                Hi[b ^ (1 << (m - 1 - i)) ^ (1 << (m - 1 - j)), b] += 2 * w
    return Hi


# ---------------------------------------------------------------- the run

def main(output_path=DEFAULT_OUTPUT):
    destination = Path(output_path)
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    out, results = [], []

    def log(msg=""):
        print(msg)
        out.append(msg)

    def check(label, ok):
        results.append((label, bool(ok)))
        log(f"   [{'PASS' if ok else 'FAIL'}] {label}")

    J, g = sp.symbols("J gamma", positive=True)
    D = sp.Symbol("Delta", real=True)
    lam = sp.Symbol("lambda")
    n = 3
    H3 = sym_hamiltonian(n, CHAIN3, J, D)
    sigma = 3 * g
    p3 = lam**3 + 8 * g * lam**2 + (16 * g**2 + (32 + 4 * D**2) * J**2) * lam + 96 * g * J**2

    log("=" * 96)
    log("THE BORN SHADOW'S CROSS TERM IN THE REFLECTION SECTORS (N = 3, uniform γ)")
    log("=" * 96)
    log()
    log("EXACT PART (sympy; J, γ, Δ symbolic, Δ real)")
    log()

    # E1
    mixed = sorted((p, q) for p in range(n + 1) for q in range(n + 1)
                   if len({hamming(i, j) for i in range(8) for j in range(8)
                           if popcount(i) == p and popcount(j) == q}) > 1)
    log(f"E1  blocks carrying more than one dephasing rate: {mixed}")
    check("E1  exactly (1,1), (1,2), (2,1), (2,2) carry more than one rate", mixed == [(1, 1), (1, 2), (2, 1), (2, 2)])

    # E2
    M11, cells11 = sym_block(H3, n, 1, 1, g)
    W11e = reflection_isometry(cells11, n, +1)
    M11e = sp.simplify(W11e.T * M11 * W11e)
    cp11e = sp.expand((lam * sp.eye(M11e.shape[0]) - M11e).det())
    log(f"E2  reflection-even (1,1): dimension {M11e.shape[0]}, characteristic polynomial = {sp.factor(cp11e)}")
    check("E2  det(λ − L) on the reflection-even (1,1) part equals λ(λ + 4γ)·p3(λ)",
          sp.expand(cp11e - lam * (lam + 4 * g) * p3) == 0)
    xdict = {(4, 2): 1, (2, 4): 1, (2, 1): 1, (1, 2): 1, (4, 1): D, (1, 4): D}
    xvec = cell_vector(cells11, xdict)
    kets1 = [b for b in range(8) if popcount(b) == 1]
    X = sp.Matrix(3, 3, lambda a, b: xdict.get((kets1[a], kets1[b]), 0))
    H1 = sp.Matrix(3, 3, lambda a, b: H3[kets1[a], kets1[b]])
    check("E2  X_Δ commutes with H on the one-excitation sector, identically in J and Δ",
          sp.simplify(H1 * X - X * H1) == sp.zeros(3))
    check("E2  L(X_Δ) = −4γ·X_Δ and L†(X_Δ) = −4γ·X_Δ on the full (1,1) block",
          sp.simplify(M11 * xvec + 4 * g * xvec) == sp.zeros(len(cells11), 1)
          and sp.simplify(M11.H * xvec + 4 * g * xvec) == sp.zeros(len(cells11), 1))
    ivec = cell_vector(cells11, {(b, b): 1 for b in kets1})
    check("E2  the identity of the sector: L and L† both send it to 0",
          sp.simplify(M11 * ivec) == sp.zeros(len(cells11), 1) and sp.simplify(M11.H * ivec) == sp.zeros(len(cells11), 1))
    check("E2  p3(0) = 96γJ² and p3(−4γ) = −16γJ²(Δ² + 2): the triple has no root at either edge value",
          sp.expand(p3.subs(lam, 0) - 96 * g * J**2) == 0
          and sp.expand(p3.subs(lam, -4 * g) + 16 * g * J**2 * (D**2 + 2)) == 0)
    dvec = cell_vector(cells11, {(4, 4): 1, (1, 1): 1, (2, 2): -2})
    in_triple = (ivec.T * dvec)[0] == 0 and (xvec.T * dvec)[0] == 0
    hopped = sp.simplify(M11 * dvec)                       # dissipator gives 0 on diagonal cells
    check("E2  d = |100><100| + |001><001| − 2|010><010| lies in the triple's space, and the commutator sends it to "
          "an operator with the entry 6iJ at |100><010|",
          in_triple and sp.simplify(hopped[cells11.index((4, 2))] - 6 * sp.I * J) == 0)

    # E3
    M12, cells12 = sym_block(H3, n, 1, 2, g)
    W12e = reflection_isometry(cells12, n, +1)
    M12e = sp.simplify(W12e.T * M12 * W12e)
    cp12e = sp.expand((lam * sp.eye(M12e.shape[0]) - M12e).det())
    log(f"E3  reflection-even (1,2): dimension {M12e.shape[0]}, characteristic polynomial = {sp.factor(cp12e)}")
    check("E3  det(λ − L) on the reflection-even (1,2) part equals (λ + 2γ)(λ + 6γ)·(−p3(−λ − 6γ))",
          sp.expand(cp12e - (lam + 2 * g) * (lam + 6 * g) * sp.expand(-p3.subs(lam, -lam - 6 * g))) == 0)
    v2 = cell_vector(cells12, {(1, 3): D, (4, 6): D, (1, 5): 1, (2, 3): 1, (2, 6): 1, (4, 5): 1})
    v6 = cell_vector(cells12, {(1, 6): 1, (2, 5): 1, (4, 3): 1})
    ok = True
    for vec, root in ((v2, -2 * g), (v6, -6 * g)):
        ok &= sp.simplify(M12 * vec - root * vec) == sp.zeros(len(cells12), 1)
        ok &= sp.simplify(M12.H * vec - root * vec) == sp.zeros(len(cells12), 1)
    check("E3  the modes at −2γ and −6γ are eigenvectors of L and of L† at symbolic J and Δ", ok)
    check("E3  the reflection-even (1,1) and (1,2) parts are palindromic partners (λ -> −λ − 6γ)",
          sp.expand(cp11e.subs(lam, -lam - 6 * g) + cp12e) == 0)

    # E4
    nu = sp.Symbol("nu")
    r = sp.Poly(sp.expand(-p3.subs(lam, -sigma - nu)), nu)
    c3, c2, c1, c0 = r.all_coeffs()
    log(f"E4  shifted cubic in ν = −(λ + 3γ): coefficients {c3}, {c2}, {sp.factor(c1)}, {sp.factor(c0)}")
    check("E4  coefficients 1, γ, (32 + 4Δ²)J² − 5γ², 3γ(γ² + 4Δ²J²)",
          sp.simplify(c3 - 1) == 0 and sp.simplify(c2 - g) == 0
          and sp.expand(c1 - ((32 + 4 * D**2) * J**2 - 5 * g**2)) == 0
          and sp.expand(c0 - 3 * g * (g**2 + 4 * D**2 * J**2)) == 0)
    check("E4  Hurwitz determinant c2·c1 − c0 = 8γ((4 − Δ²)J² − γ²)",
          sp.expand(c2 * c1 - c0 - 8 * g * ((4 - D**2) * J**2 - g**2)) == 0)

    # E5
    odd = {}
    for (p, q) in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        Mb, cells = sym_block(H3, n, p, q, g)
        Wo = reflection_isometry(cells, n, -1)
        Mo = sp.simplify(Wo.T * Mb * Wo)
        odd[(p, q)] = (Mo, sp.expand((lam * sp.eye(Mo.shape[0]) - Mo).det()), Wo, cells)
        check(f"E5  block ({p},{q}) reflection-odd: dimension 4, trace −3γ per mode",
              Mo.shape[0] == 4 and sp.simplify(Mo.trace() + sigma * 4) == 0)
    cpo11, cpo12 = odd[(1, 1)][1], odd[(1, 2)][1]
    log(f"E5  reflection-odd (1,1): characteristic polynomial = {cpo11}")
    check("E5  the reflection-odd (1,1) and (1,2) parts are palindromic partners",
          sp.expand(cpo11.subs(lam, -lam - 6 * g) - cpo12) == 0)
    check("E5  the reflection-odd (2,2) and (2,1) parts have the characteristic polynomials of (1,1) and (1,2)",
          sp.expand(odd[(2, 2)][1] - cpo11) == 0 and sp.expand(odd[(2, 1)][1] - cpo12) == 0)
    Wo12, c12 = odd[(1, 2)][2], odd[(1, 2)][3]
    D12 = sp.simplify(Wo12.T * sp.diag(*[-2 * g * hamming(i, j) for i, j in c12]) * Wo12)
    check("E5  the dissipator on the reflection-odd (1,2) part is diagonal with one −6γ direction and three at −2γ",
          D12.is_diagonal() and sorted([D12[k, k] for k in range(4)], key=lambda v: v.subs(g, 1)) == [-6 * g, -2 * g, -2 * g, -2 * g])
    Wo11, c11 = odd[(1, 1)][2], odd[(1, 1)][3]
    levels = sp.diag(*[-2 * g * hamming(i, j) for i, j in c11])
    Dodd = sp.simplify(Wo11.T * levels * Wo11)
    check("E5  the dissipator on the reflection-odd (1,1) part is diagonal with one level-0 direction and three at −4γ",
          Dodd.is_diagonal() and sorted([Dodd[k, k] for k in range(4)], key=lambda v: v.subs(g, 1)) == [-4 * g, -4 * g, -4 * g, 0])
    res_edge = edge_resultant(cpo11, lam, -4 * g)
    res_zero = edge_resultant(cpo11, lam, 0)
    log(f"E5  resultant for a root on Re λ = −4γ: {res_edge}")
    log(f"E5  resultant for a root on Re λ = 0:   {res_zero}")
    check("E5  the Re λ = −4γ resultant equals 2²⁴·Δ⁴J¹²γ⁴ (zero only at Δ = 0)",
          sp.expand(res_edge - 2**24 * D**4 * J**12 * g**4) == 0)
    check("E5  the Re λ = 0 resultant equals 2²⁴·J²γ⁴(J² + 2γ²)(Δ⁴J⁴ + 9Δ²J⁴ + 20Δ²J²γ² + 72J²γ² + 64γ⁴)², "
          "which no real J ≠ 0 and Δ can make zero",
          sp.expand(res_zero - 2**24 * J**2 * g**4 * (J**2 + 2 * g**2)
                    * (D**4 * J**4 + 9 * D**2 * J**4 + 20 * D**2 * J**2 * g**2 + 72 * J**2 * g**2 + 64 * g**4)**2) == 0)
    mu = sp.Symbol("mu")
    qmu = sp.Poly(sp.expand(cpo11.subs(lam, mu - sigma)), mu)
    m1 = sp.factor(qmu.coeff_monomial(mu))
    log(f"E5  in μ = λ + 3γ: μ³ coefficient {qmu.coeff_monomial(mu**3)}, μ¹ coefficient {m1}")
    check("E5  μ³ coefficient 0 and μ¹ coefficient −8γ(Δ²J² + γ²) ≠ 0: the cut splits every odd part",
          qmu.coeff_monomial(mu**3) == 0 and sp.expand(m1 + 8 * g * (D**2 * J**2 + g**2)) == 0)
    xy = sp.expand(cpo11.subs(D, 0))
    check("E5  at Δ = 0 the odd (1,1) part factors as (λ² + 4γλ + 8J²)((λ + 4γ)² + 8J²)",
          sp.expand(xy - (lam**2 + 4 * g * lam + 8 * J**2) * ((lam + 4 * g)**2 + 8 * J**2)) == 0)
    quad = sp.Poly(sp.expand((lam**2 + 4 * g * lam + 8 * J**2).subs(lam, -sigma - nu)), nu)
    check("E5  its first two roots, shifted to the cut, solve ν² + 2γν + (8J² − 3γ²) = 0: slower than Σγ exactly when 8J² > 3γ²",
          quad.all_coeffs() == [1, 2 * g, 8 * J**2 - 3 * g**2])

    # E6
    log()
    graphs = {
        "chain": lambda m: [(i, i + 1) for i in range(m - 1)],
        "ring": lambda m: [(i, (i + 1) % m) for i in range(m)],
        "star": lambda m: [(0, i) for i in range(1, m)],
        "complete": lambda m: list(itertools.combinations(range(m), 2)),
    }
    rng_int = np.random.default_rng(SEED)
    for m in (3, 4, 5):
        d = 2**m
        ones = [b for b in range(d) if popcount(b) == 1]
        Xw = np.zeros((d, d), dtype=np.int64)
        for a in ones:
            for b in ones:
                if a != b:
                    Xw[a, b] = 1                               # N|W><W| − P_1
        for name, bonds in graphs.items():
            Hi = integer_xxz(m, [(i, j, 1) for i, j in bonds(m)])
            check(f"E6  N = {m}, {name}, unit couplings: [H, N|W><W| − P_1] = 0 exactly", not (Hi @ Xw - Xw @ Hi).any())
        weighted = [(i, j, int(w)) for (i, j), w in zip(graphs["complete"](m), rng_int.integers(1, 10, size=m * (m - 1) // 2))]
        Hw = integer_xxz(m, weighted)
        check(f"E6  N = {m}, complete graph, random integer couplings {[w for _, _, w in weighted]}: "
              "[H, N|W><W| − P_1] = 0 exactly", not (Hw @ Xw - Xw @ Hw).any())
    Hc = integer_xxz(3, [(0, 1, 1), (1, 2, 1)], Delta=2)
    X3 = np.zeros((8, 8), dtype=np.int64)
    for a in (1, 2, 4):
        for b in (1, 2, 4):
            if a != b:
                X3[a, b] = 1
    check("E6  control: on the XXZ chain at Δ = 2 the commutator is nonzero (the check can fail)",
          (Hc @ X3 - X3 @ Hc).any())

    # E7
    log()
    Hr = sym_hamiltonian(n, RING3, J, D)
    q3 = lam**3 + 8 * g * lam**2 + (16 * g**2 + 36 * J**2) * lam + 96 * g * J**2
    q3i = sp.expand(-q3.subs(lam, -lam - 6 * g))
    expected = {((1, 1), +1): lam * (lam + 4 * g) * q3, ((1, 1), -1): (lam + 4 * g) * q3,
                ((2, 2), +1): lam * (lam + 4 * g) * q3, ((2, 2), -1): (lam + 4 * g) * q3,
                ((1, 2), +1): (lam + 2 * g) * (lam + 6 * g) * q3i, ((1, 2), -1): (lam + 2 * g) * q3i,
                ((2, 1), +1): (lam + 2 * g) * (lam + 6 * g) * q3i, ((2, 1), -1): (lam + 2 * g) * q3i}
    for ((p, q), s), target in expected.items():
        Mb, cells = sym_block(Hr, n, p, q, g)
        Ws = reflection_isometry(cells, n, s)
        Ms = sp.simplify(Ws.T * Mb * Ws)
        cp = sp.expand((lam * sp.eye(Ms.shape[0]) - Ms).det())
        check(f"E7  ring, block ({p},{q}) reflection-{'even' if s > 0 else 'odd'}: window-edge factors times "
              f"{'q3' if (p, q) in ((1, 1), (2, 2)) else 'its palindromic image'}, for every Δ",
              sp.expand(cp - sp.expand(target)) == 0)

    # ------------------------------------------------------------ numerical readings
    log()
    log("NUMERICAL READINGS (J = 1, γ = 0.05 unless stated; ||P_s^H P_f|| = 0 means the two sides are orthogonal)")
    for label, m, bonds, dl in (("chain N=3, Δ = 1", 3, CHAIN3, 1.0), ("chain N=3, Δ = 0", 3, CHAIN3, 0.0),
                                ("ring N=3, Δ = 1", 3, RING3, 1.0), ("chain N=4, Δ = 1", 4, graphs["chain"](4), 1.0)):
        rows, wholly, on_res, gap = num_sector_table(m, bonds, GAMMA_PAGE, J_PAGE, dl, with_cut=True)
        log(f"   {label}: sectors with modes on both sides of the cut")
        for p, q, s, dim, k, nono in rows:
            log(f"      ({p},{q}) reflection-{'even' if s > 0 else 'odd '}: dim {dim:2d}, slow {k:2d}, ||P_s^H P_f|| = {nono:.2e}")
        wl = ", ".join(f"({p},{q}){'+' if s > 0 else '-'}" for p, q, s in wholly if (p, q) not in ((0, 0), (m, m)))
        on_text = "no mode on the cut" if on_res is None else f"largest on-cut distance from Σγ {on_res:.1e}"
        log(f"      sectors wholly on the cut: {wl or 'none'}; {on_text}; gap to the nearest off-cut mode {gap:.2e}")
    log()
    log("   rates −Re λ/γ of every mixed part, chain N=3, Δ = 1 (compare PROOF_ABSORPTION_THEOREM §4.6)")
    for q_ratio in (1.5, 20.0):
        log(f"      J/γ = {q_ratio}:")
        for (p, q) in [(1, 1), (1, 2)]:
            for sgn in (+1, -1):
                M = num_sector(3, CHAIN3, GAMMA_PAGE, q_ratio * GAMMA_PAGE, 1.0, p, q, sgn)
                rates = np.sort(-np.linalg.eigvals(M).real / GAMMA_PAGE)
                log(f"         ({p},{q}) reflection-{'even' if sgn > 0 else 'odd '}: " + "  ".join(f"{x:.4f}" for x in rates))
    log()
    log("   the thresholds, bracketed: worst ||P_s^H P_f|| at 0.98 and 1.02 of the exact value")
    for Dv in (0.0, 1.0, 1.5, 1.9):
        jc = GAMMA_PAGE / np.sqrt(4 - Dv**2)
        worst = [max([x[5] for x in num_sector_table(3, CHAIN3, GAMMA_PAGE, f * jc, Dv) if x[2] > 0] or [0.0])
                 for f in (0.98, 1.02)]
        log(f"      reflection-even, Δ = {Dv:3.1f}: J_c/γ = {jc / GAMMA_PAGE:.4f}   0.98·J_c: {worst[0]:.2e}   1.02·J_c: {worst[1]:.2e}")
    jo = GAMMA_PAGE * np.sqrt(3 / 8)
    worst = [max([x[5] for x in num_sector_table(3, CHAIN3, GAMMA_PAGE, f * jo, 0.0) if x[2] < 0] or [0.0])
             for f in (0.98, 1.02)]
    log(f"      reflection-odd, Δ = 0:   J_o/γ = {jo / GAMMA_PAGE:.4f}   0.98·J_o: {worst[0]:.2e}   1.02·J_o: {worst[1]:.2e}")
    log()
    L3 = num_liouvillian(num_hamiltonian(3, CHAIN3), 3, GAMMA_PAGE)
    P3, _ = spectral_projector(L3, slow_select(3 * GAMMA_PAGE))
    for letters in ("+0+", "++0"):
        rho0 = product_state(letters)
        oddpart = (rho0 - reflect_operator(rho0, 3)) / 2
        fold = first_downward_quarter(L3, rho0)
        log(f"   |{letters}>, Δ = 1: purity split slow / fast / cross, and the reflection-odd part's cross term "
            f"(fold, CΨ = 1/4 in the page's book, at t = {fold:.3f})")
        for t in (0.0, 1.0, 3.0, fold):
            ps, pf, pc = purity_split(L3, P3, rho0, t)
            oc = purity_split(L3, P3, oddpart, t)[2]
            log(f"      t = {t:5.3f}: {ps:.4f} / {pf:.4f} / {pc:+.2e} ({100 * pc / (ps + pf + pc):+.4f}% of the purity)"
                f"   reflection-odd part: {oc:+.2e}")
    jlow = 0.5 * GAMMA_PAGE / np.sqrt(3)
    Llow = num_liouvillian(num_hamiltonian(3, CHAIN3, J=jlow), 3, GAMMA_PAGE)
    Plow, _ = spectral_projector(Llow, slow_select(3 * GAMMA_PAGE))
    rho0 = product_state("+0+")
    log(f"   |+0+>, Δ = 1, below the threshold (J = 0.5·γ/√3 = {jlow:.5f}): cross term at t = 0, 3, 10: "
        + "  ".join(f"{purity_split(Llow, Plow, rho0, t)[2]:+.2e}" for t in (0.0, 3.0, 10.0)))
    log()
    rng = np.random.default_rng(SEED)
    log("   the W law against expm: largest |<W|rho(t)|W> − (P_1/N + e^(−4γt)(<W|rho(0)|W> − P_1/N))| over t in [0, 20]")
    for m, name in ((3, "chain"), (4, "chain"), (4, "star")):
        d = 2**m
        Lm = num_liouvillian(num_hamiltonian(m, graphs[name](m)), m, GAMMA_PAGE)
        a = rng.normal(size=d) + 1j * rng.normal(size=d)
        a /= np.linalg.norm(a)
        rho0 = np.outer(a, a.conj())
        w = np.zeros(d, dtype=complex)
        for b in range(d):
            if popcount(b) == 1:
                w[b] = 1 / np.sqrt(m)
        p1 = sum(np.real(rho0[b, b]) for b in range(d) if popcount(b) == 1)
        fw0 = np.real(w.conj() @ rho0 @ w)
        dev = max(abs(np.real(w.conj() @ (expm(Lm * t) @ rho0.ravel()).reshape(d, d) @ w)
                      - (p1 / m + np.exp(-4 * GAMMA_PAGE * t) * (fw0 - p1 / m))) for t in np.linspace(0, 20, 41))
        log(f"      N = {m} {name}, a random pure state: {dev:.1e}")

    log()
    passed = sum(ok for _, ok in results)
    log(f"{passed} of {len(results)} checks pass.")
    for label, ok in results:
        if not ok:
            log(f"   FAILED: {label}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    sys.exit(main(parser.parse_args().output))
