"""Where do the corner block's exceptional couplings actually sit?

`docs/proofs/PROOF_R90_FROZEN_DIVISOR` section 9: on a NON-uniform locus
profile the frozen divisor gains an extra algebraic dimension at
isolated NONZERO real couplings, "one small coupling away from
anywhere", and section 12 leaves the COUNT open, reporting that at N = 6
two generic profiles on the same locus give four real couplings and
eight. An earlier corner-beat pass (2026-08-18) reported "no exceptional
couplings" after seeing multiplicity jumps only at Q <= 0.25 and reading
them as "the J -> 0 root section 9 excludes". Section 9 excludes J = 0
EXACTLY ("Exceptional coupling below always means a nonzero root"), not
small J, so that reading had to be checked. It was wrong twice over.

This script produces every number Amendment 2 section 2.4b cites, for
the two FLOWN lighting profiles C = [2,0,0,2,2,0] and
C' = [2,2,2,0,0,0] (both on the R90 locus, both gamma-bar != 0, so the
taxed stratum section 9 is about) and for the uniform control:

  1. the profiles are on the locus, and the generic algebraic
     multiplicity of the frozen root is exactly floor(N/2) = 3 and
     semisimple;
  2. the polynomial whose nonzero real roots ARE the exceptional
     couplings, built exactly, and its real roots by Sturm;
  3. the Jordan structure AT each exceptional coupling, exactly;
  4. the conditioning at the FLOWN point Q = 10;
  5. the weight the frozen modes carry on the population cells.

Steps 1-3 use no floating point at all: characteristic polynomials mod
p over verified primes, lifted by CRT, then exact real-root isolation
over the integers. A rank test at a tuned degeneracy is the least
trustworthy measurement there is, which is why the defect is reached
through a polynomial and not through a numerical rank. Steps 4 and 5 are
floating point by nature (they are conditioning statements) and are
labelled as such.

Conventions, pinned:
  - one-magnon sector, H per bond = XX + YY + ZZ, hopping element 2 and
    a ZZ diagonal: `bond_H_1m` here is verbatim the committed gate's
    (`simulations/corner_beat_gate.py`). The `--xy` variant drops the ZZ
    diagonal and is a model-sensitivity check, not the flown model.
  - Z-dephasing in the LINDBLAD book: a one-magnon coherence (a,b) with
    a != b decays at 2(gamma_a + gamma_b), and the diagonal cells a = b
    pay NOTHING. That even defect is what the frozen divisor lives on.
  - gamma-bar is scaled to 1, so gamma IS the integer profile and the
    coupling variable IS Q = J/gamma-bar; the generator is homogeneous
    of degree 1 in (J, gamma) jointly, so only the ratio matters and
    Q = 10 here is the flight's Q = 10. The frozen root is then
    lambda_0 = -4, which is F140's -4*gamma-bar.
  - substitution u = i*J makes every matrix entry an INTEGER; a real
    nonzero coupling is a root with u^2 < 0, i.e. Q = sqrt(-u^2).

Run: python simulations/_corner_beat_exceptional_verify.py
"""

import sympy as sp

N = 6
D6 = 1 << N
FLOWN_Q = 10
LAM0 = -4                      # -4 * gamma-bar at gamma-bar = 1
NPTS = 41                      # interpolation nodes in u; degree is 30

# Verified primes just below 2^61. Guessed moduli are the fastest way to
# make every CRT lift silently wrong: the first run of this script used
# four composite ones and returned degree 40 where the a-priori bound is
# 33, which is how it was caught.
PRIMES = [int(p) for p in (sp.prevprime(1 << 61),
                           sp.prevprime(sp.prevprime(1 << 61)),
                           sp.prevprime(sp.prevprime(sp.prevprime(1 << 61))),
                           sp.prevprime(sp.prevprime(sp.prevprime(
                               sp.prevprime(1 << 61)))))]
assert all(sp.isprime(p) for p in PRIMES)

PROFILES = [
    ("C  (maximizing, flown)",     [2, 0, 0, 2, 2, 0]),
    ("C' (non-maximizing, flown)", [2, 2, 2, 0, 0, 0]),
    ("U  (uniform, control)",      [1, 1, 1, 1, 1, 1]),
]


def bond_H_1m(bonds, zz=True):
    """The committed gate's one-magnon Hamiltonian block, verbatim."""
    H = [[0] * D6 for _ in range(D6)]
    for a in range(D6):
        for l in bonds:
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a][a] += za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a][a ^ (1 << l) ^ (1 << (l + 1))] += 2
    cfg = [1 << l for l in range(N)]
    M = [[H[a][b] for b in cfg] for a in cfg]
    if not zz:
        for a in range(N):
            M[a][a] = 0
    return M


IDX = [(a, b) for a in range(N) for b in range(N)]
DIM = len(IDX)
DIAG_CELLS = [i for i, (a, b) in enumerate(IDX) if a == b]


def corner_block(h, gamma, u):
    """M on the (1,1) coherence block, DIM x DIM, in the variable u = i*J.

    -i*J*[H, rho] is -u*(H rho - rho H), an integer matrix whenever u is;
    the dissipator is the diagonal -2(gamma_a + gamma_b) off the
    population cells and zero on them."""
    M = [[0] * DIM for _ in range(DIM)]
    for r, (a, b) in enumerate(IDX):
        for c, (ap, bp) in enumerate(IDX):
            v = 0
            if b == bp:
                v -= u * h[a][ap]
            if a == ap:
                v += u * h[bp][b]
            M[r][c] = v
        if a != b:
            M[r][r] -= 2 * (gamma[a] + gamma[b])
    return M


# --------------------------------------------------------------- GF(p)

def det_mod(mat, p):
    """Determinant mod p by Gaussian elimination. `mat` is consumed."""
    n = len(mat)
    det = 1
    for col in range(n):
        piv = next((r for r in range(col, n) if mat[r][col] % p), None)
        if piv is None:
            return 0
        if piv != col:
            mat[col], mat[piv] = mat[piv], mat[col]
            det = -det
        pv = mat[col][col] % p
        det = det * pv % p
        inv = pow(pv, p - 2, p)
        row = mat[col]
        for r in range(col + 1, n):
            f = mat[r][col] * inv % p
            if f:
                rr = mat[r]
                for c in range(col, n):
                    rr[c] = (rr[c] - f * row[c]) % p
    return det % p


def rank_mod(mat, p):
    n, m = len(mat), len(mat[0])
    mat = [row[:] for row in mat]
    rank = row = 0
    for col in range(m):
        piv = next((r for r in range(row, n) if mat[r][col] % p), None)
        if piv is None:
            continue
        mat[row], mat[piv] = mat[piv], mat[row]
        inv = pow(mat[row][col], p - 2, p)
        pr = mat[row]
        for r in range(row + 1, n):
            f = mat[r][col] * inv % p
            if f:
                rr = mat[r]
                for c in range(col, m):
                    rr[c] = (rr[c] - f * pr[c]) % p
        rank += 1
        row += 1
        if row == n:
            break
    return rank


def matmul_mod(A, B, p):
    Bt = list(zip(*B))
    return [[sum(x * y for x, y in zip(row, col)) % p for col in Bt]
            for row in A]


def lagrange_mod(xs, ys, p):
    """Coefficients (ascending) of the interpolating polynomial mod p."""
    n = len(xs)
    coeffs = [0] * n
    for i in range(n):
        num, deg, den = [1] + [0] * n, 0, 1
        for j in range(n):
            if j == i:
                continue
            new = [0] * (n + 1)
            for k in range(deg + 1):
                new[k + 1] = (new[k + 1] + num[k]) % p
                new[k] = (new[k] - num[k] * xs[j]) % p
            num, deg = new, deg + 1
            den = den * (xs[i] - xs[j]) % p
        f = ys[i] * pow(den, p - 2, p) % p
        for k in range(n):
            coeffs[k] = (coeffs[k] + f * num[k]) % p
    return coeffs


def eps_coeffs_mod(A, p, kmax):
    """Coefficients of eps^0..eps^(kmax-1) in det(A - eps*I) mod p."""
    n = len(A)
    xs = list(range(n + 1))
    ys = [det_mod([[(A[r][c] - (e if r == c else 0)) % p
                    for c in range(n)] for r in range(n)], p) for e in xs]
    return lagrange_mod(xs, ys, p)[:kmax]


def crt(residues, primes):
    M = 1
    for p in primes:
        M *= p
    x = 0
    for r, p in zip(residues, primes):
        Mi = M // p
        x = (x + r * Mi * pow(Mi % p, p - 2, p)) % M
    return (x - M if x > M // 2 else x), M


def shifted(h, gamma, u, p):
    A = corner_block(h, gamma, u)
    return [[(A[r][c] - (LAM0 if r == c else 0)) % p for c in range(DIM)]
            for r in range(DIM)]


def kernel_dims(h, gamma, u, p, kmax=4):
    """Kernel dimensions of (M - lambda_0)^k, exactly, mod p."""
    A = shifted(h, gamma, u, p)
    dims, P = [], [row[:] for row in A]
    for _ in range(kmax):
        dims.append(DIM - rank_mod(P, p))
        P = matmul_mod(P, A, p)
    return dims


def exceptional_polynomial(h, gamma):
    """The eps^3 coefficient of det(M(u) - (lambda_0 + eps) I), exactly.

    Returns (integer coefficients ascending in u, low-order report, the
    CRT modulus). The eps^0..eps^2 coefficients must vanish identically
    when the generic multiplicity is exactly 3; they are returned so the
    caller checks rather than assumes."""
    us = list(range(NPTS))
    per_prime, low = [], None
    for p in PRIMES:
        v3, v012 = [], []
        for u in us:
            co = eps_coeffs_mod(shifted(h, gamma, u, p), p, kmax=4)
            v012.append(co[:3])
            v3.append(co[3])
        per_prime.append(lagrange_mod(us, v3, p))
        report = [all(v[k] == 0 for v in v012) for k in range(3)]
        low = report if low is None else [a and b for a, b in zip(low, report)]
    coeffs, mod = [], None
    for k in range(len(per_prime[0])):
        c, mod = crt([pp[k] for pp in per_prime], PRIMES)
        coeffs.append(c)
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs, low, mod


def exceptional_couplings(coeffs):
    """Real nonzero couplings Q > 0, exactly, from the u-polynomial."""
    w = sp.Symbol("w")
    assert not [c for k, c in enumerate(coeffs) if k % 2 and c], \
        "the cofactor must be even in u (a polynomial in J^2)"
    gw = sp.Poly([coeffs[k] for k in range(len(coeffs) - 1, -1, -2)], w)
    zeros = 0
    while gw.degree() > 0 and gw.eval(0) == 0:
        gw = sp.Poly(sp.div(gw.as_expr(), w, w)[0], w)
        zeros += 1
    if gw.degree() == 0:
        return [], gw, zeros
    neg = [r for r in gw.real_roots() if r < 0]
    return neg, gw, zeros


def jordan_at(h, gamma, wval):
    """Kernel dims at the exceptional coupling w = u^2, exactly.

    Picks a prime in which w is a square, so u has an integer
    representative and the same integer machinery applies; no float and
    no tuned rank test."""
    wq = sp.Rational(sp.nsimplify(wval)) if wval.is_rational else None
    q = 1 << 61
    for _ in range(400):
        q = int(sp.prevprime(q))
        if wq is not None:
            target = (wq.p * pow(wq.q, q - 2, q)) % q
        else:
            # w is an algebraic number: use its minimal polynomial's
            # roots mod q instead
            mp = sp.Poly(sp.minimal_polynomial(wval, sp.Symbol("x")),
                         sp.Symbol("x"), modulus=q)
            lin = [(-c.as_expr().subs(sp.Symbol("x"), 0)) % q
                   for c, _ in mp.factor_list()[1] if c.degree() == 1]
            target = next((t for t in lin
                           if sp.sqrt_mod(int(t) % q, q) is not None), None)
            if target is None:
                continue
            target = int(target) % q
        r = sp.sqrt_mod(target, q)
        if r is None:
            continue
        u = int(r)
        return q, kernel_dims(h, gamma, u, q), kernel_dims(h, gamma,
                                                           (u + 1) % q, q)
    raise RuntimeError("no usable prime found")


def jordan_blocks(dims):
    """Block sizes from the kernel dimensions of the powers.

    increments m_k = dim_k - dim_{k-1} count the blocks of size >= k, so
    the blocks of size exactly k are m_k - m_{k+1}. Writing the size out
    of (algebraic - geometric + 1) happens to be right when there is a
    single nontrivial block and is wrong the moment there are two."""
    m = [dims[0]] + [dims[k] - dims[k - 1] for k in range(1, len(dims))]
    m.append(0)
    return {k: m[k - 1] - m[k] for k in range(1, len(m)) if m[k - 1] - m[k]}


def conditioning_at_flown_point(h, gamma):
    """Float, and labelled so: conditioning is not an exact question."""
    import numpy as np
    M = np.zeros((DIM, DIM), dtype=complex)
    for r, (a, b) in enumerate(IDX):
        for c, (ap, bp) in enumerate(IDX):
            v = 0j
            if b == bp:
                v += -1j * FLOWN_Q * h[a][ap]
            if a == ap:
                v += +1j * FLOWN_Q * h[bp][b]
            M[r, c] = v
        if a != b:
            M[r, r] += -2 * (gamma[a] + gamma[b])
    lam, V = np.linalg.eig(M)
    frozen = [i for i in range(DIM) if abs(lam[i] - LAM0) < 1e-8]
    Vi = np.linalg.inv(V)
    kappa = [float(np.linalg.norm(V[:, i]) * np.linalg.norm(Vi[i, :]))
             for i in frozen]
    pop = [float(np.linalg.norm(V[DIAG_CELLS, i] / np.linalg.norm(V[:, i])))
           for i in frozen]
    gap = sorted(abs(lam - LAM0))[len(frozen)]
    return len(frozen), kappa, pop, float(gap)


def main():
    for zz, model in ((True, "Heisenberg (the flown H)"), (False, "XY (no ZZ)")):
        h = bond_H_1m(list(range(N - 1)), zz=zz)
        print(f"\n{'=' * 68}\nMODEL: {model}\n{'=' * 68}")
        for name, gamma in PROFILES:
            gbar = sum(gamma) // N
            assert sum(gamma) % N == 0 and gbar == 1
            rev = list(reversed(gamma))
            on_locus = all(2 * gbar - rev[i] == gamma[i] for i in range(N))
            print(f"\n--- {name}   gamma = {gamma}, gamma-bar = {gbar}")
            print(f"    on the R90 locus (x -> 2*avg - reverse): {on_locus}")

            for probe in (3, 7, 13):
                print(f"    kernel dims of (M - lam0)^k at u = {probe}: "
                      f"{kernel_dims(h, gamma, probe, PRIMES[0])}"
                      f"   (a control point, not a real coupling: u = i*J)")

            coeffs, low, mod = exceptional_polynomial(h, gamma)
            print(f"    eps^0, eps^1, eps^2 vanish identically, over every "
                  f"prime: {low}")
            print(f"    cofactor degree in u: {len(coeffs) - 1} "
                  f"(the proof's N(N-1) = {N * (N - 1)})")
            peak = max(abs(c) for c in coeffs)
            print(f"    largest |coefficient| {peak} against the CRT modulus "
                  f"{mod} (ratio {float(2 * peak / mod):.3e}, so no wrap)")

            neg, gw, zeros = exceptional_couplings(coeffs)
            print(f"    J = 0 root multiplicity in w = u^2: {zeros}; "
                  f"reduced g(w) degree {gw.degree()}")
            print(f"    g(w) = {sp.factor(gw.as_expr())}")
            if not neg:
                print("    NO real nonzero exceptional coupling")
                continue
            for r in neg:
                Q = sp.sqrt(-r)
                print(f"    exceptional coupling: Q* = {sp.nsimplify(Q)} "
                      f"= {float(Q):.9f}   (distance to the flown Q = "
                      f"{FLOWN_Q}: {abs(float(Q) - FLOWN_Q):.6f})")
                q, at, near = jordan_at(h, gamma, r)
                blocks = jordan_blocks(at)
                shape = ", ".join(f"{n} of size {k}"
                                  for k, n in sorted(blocks.items()))
                print(f"      kernel dims AT it (exact, mod {q}): {at}"
                      f"  -> geometric {at[0]}, algebraic {at[-1]}, "
                      f"Jordan blocks: {shape}")
                print(f"      kernel dims at a neighbour: {near}")

            if zz:
                nf, kappa, pop, gap = conditioning_at_flown_point(h, gamma)
                print(f"    AT THE FLOWN POINT Q = {FLOWN_Q} (float, since "
                      f"conditioning is not an exact question):")
                print(f"      frozen modes {nf}; eigenvalue condition "
                      f"numbers {[round(k, 3) for k in kappa]}")
                print(f"      separation to the nearest other eigenvalue "
                      f"{gap:.6f}")
                print(f"      weight of the frozen modes on the POPULATION "
                      f"cells: {['%.2e' % x for x in pop]}")


if __name__ == "__main__":
    main()
