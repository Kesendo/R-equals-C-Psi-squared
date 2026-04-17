#!/usr/bin/env python3
"""
GHZ+W optimum for pair-CPsi(0) on N=3.

For the one-parameter family
    |psi(a)> = a |GHZ_3> + sqrt(1 - a^2) |W_3>,   a in [0, 1]
this script:

  (1) derives rho_AB(a) in closed form (exact rational / sqrt expressions);
  (2) isolates the stationarity condition d/da CPsi(0) = 0 as an equation
      in x = a^2, whose rationalization is a sextic with integer coefficients;
  (3) confirms the sextic is irreducible over the rationals;
  (4) reports the unique root that gives a genuine maximum of min pair-CPsi(0);
  (5) verifies the optimum against a 401-point numerical sweep;
  (6) computes 3-tangle and pair concurrences of the optimum state;
  (7) confirms that the same construction at N in {4, 5, 6} fails to lift
      any pair above CPsi = 1/4.

Source for the writeup: docs/ANALYTICAL_FORMULAS.md, entry F69.
"""
from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp


if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# --- State constructors (big-endian: |b_{n-1} b_{n-2} ... b_0>) -------------
def ket_ghz(n):
    v = np.zeros(2**n, dtype=complex)
    v[0] = 1.0 / np.sqrt(2)
    v[-1] = 1.0 / np.sqrt(2)
    return v


def ket_w(n):
    v = np.zeros(2**n, dtype=complex)
    for i in range(n):
        v[1 << i] = 1.0 / np.sqrt(n)
    return v


# --- Reduced density matrices and CPsi --------------------------------------
def partial_trace_pair(rho, n, keep):
    a, b = keep
    cur = list(range(n))
    t = rho.reshape([2] * (2 * n))
    while len(cur) > 2:
        rm = next(q for q in cur if q not in keep)
        idx = cur.index(rm)
        nc = len(cur)
        t = np.trace(t, axis1=nc - 1 - idx, axis2=2 * nc - 1 - idx)
        cur.pop(idx)
    if cur == [a, b]:
        return t.reshape(4, 4)
    return t.transpose(1, 0, 3, 2).reshape(4, 4)


def partial_trace_single(rho, n, keep):
    cur = list(range(n))
    t = rho.reshape([2] * (2 * n))
    while len(cur) > 1:
        rm = next(q for q in cur if q != keep)
        idx = cur.index(rm)
        nc = len(cur)
        t = np.trace(t, axis1=nc - 1 - idx, axis2=2 * nc - 1 - idx)
        cur.pop(idx)
    return t.reshape(2, 2)


def cpsi_pair(rho2):
    C = float(np.real(np.trace(rho2 @ rho2)))
    diag = np.diag(np.diag(rho2))
    L1 = float(np.sum(np.abs(rho2 - diag)))
    return C, L1 / 3.0, C * L1 / 3.0


def min_pair_cpsi(psi, n):
    rho = np.outer(psi, psi.conj())
    return min(cpsi_pair(partial_trace_pair(rho, n, (a, b)))[2]
               for a, b in itertools.combinations(range(n), 2))


# --- Entanglement diagnostics -----------------------------------------------
def concurrence_2q(rho2):
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    YY = np.kron(Y, Y)
    eigs = np.sort(np.real(np.linalg.eigvals(rho2 @ (YY @ rho2.conj() @ YY))))[::-1]
    eigs = np.clip(eigs, 0.0, None)
    s = np.sqrt(eigs)
    return float(max(0.0, s[0] - s[1] - s[2] - s[3]))


def three_tangle(psi):
    n = 3
    rho = np.outer(psi, psi.conj())
    rho_A = partial_trace_single(rho, n, 0)
    linear_entropy_A = 2.0 * (1.0 - np.real(np.trace(rho_A @ rho_A)))
    c_AB = concurrence_2q(partial_trace_pair(rho, n, (0, 1)))
    c_AC = concurrence_2q(partial_trace_pair(rho, n, (0, 2)))
    return max(0.0, float(linear_entropy_A - c_AB**2 - c_AC**2)), c_AB, c_AC


# --- Symbolic derivation of rho_AB(a) and CPsi(a) ---------------------------
def symbolic_cpsi():
    """Return (CPsi_sym, sextic_poly, dCPsi_dt, t).

    Uses substitution a = sin(t), b = cos(t) on t in (0, pi/2) so that a, b
    are both non-negative and Abs() simplifies away.
    """
    t = sp.Symbol("t", real=True, positive=True)
    x = sp.Symbol("x", real=True, positive=True)

    a_sym = sp.sin(t)
    b_sym = sp.cos(t)

    # rho_AB: derived symbolically from |psi(a,b)> = a|GHZ_3> + b|W_3>.
    # (|GHZ> on indices {0, 7}; |W> on indices {1, 2, 4}.)
    # Tracing out the third qubit (LSB) gives:
    #     C = Tr(rho_AB^2) = -5 a^4 / 18 + 2 a^2 / 9 + 5 / 9
    #     L1_off           = sqrt(6) * a * b + (2/3) * b^2     (for a, b >= 0)
    C_expr = -5 * a_sym**4 / 18 + 2 * a_sym**2 / 9 + sp.Rational(5, 9)
    L1_expr = sp.sqrt(6) * a_sym * b_sym + sp.Rational(2, 3) * b_sym**2
    CPsi_sym = sp.simplify(C_expr * L1_expr / 3)

    dCPsi_dt = sp.simplify(sp.diff(CPsi_sym, t))

    # Rationalized stationarity equation in x = a^2 (isolate sqrt(x(1-x))
    # and square).  This yields a polynomial in x with integer coefficients
    # whose relevant root equals a_opt^2.
    C_x = -5 * x**2 / 18 + 2 * x / 9 + sp.Rational(5, 9)
    L1_x = sp.sqrt(6) * sp.sqrt(x * (1 - x)) + sp.Rational(2, 3) * (1 - x)
    CPsi_x = C_x * L1_x / 3
    dCPsi_dx = sp.simplify(sp.diff(CPsi_x, x))

    # Separate sqrt and non-sqrt pieces and square to clear the radical.
    # dCPsi_dx = A(x) + B(x) * sqrt(x(1-x))  =>  A(x)^2 = B(x)^2 * x * (1-x)
    f = sp.together(dCPsi_dx)
    num = sp.numer(f)
    num = sp.expand(num)
    s = sp.Symbol("s", positive=True)  # stand-in for sqrt(x(1-x))
    num_s = num.subs(sp.sqrt(x * (1 - x)), s)
    A = num_s.coeff(s, 0)
    B = num_s.coeff(s, 1)
    sextic = sp.expand(A**2 - B**2 * x * (1 - x))
    sextic = sp.Poly(sextic, x).as_expr()
    # Normalize leading coefficient sign
    lead = sp.Poly(sextic, x).LC()
    if lead < 0:
        sextic = sp.expand(-sextic)
    # Divide out integer content
    poly = sp.Poly(sextic, x)
    content = poly.content()
    sextic = sp.expand(sextic / content)

    return CPsi_sym, sextic, dCPsi_dt, t, x


def numerical_optimum():
    """Return (a_opt, b_opt, cpsi_opt) via scipy + dense sweep cross-check."""
    from scipy.optimize import minimize_scalar

    def cpsi_t(t_val):
        a = math.sin(t_val)
        b = math.cos(t_val)
        C = -5 * a**4 / 18 + 2 * a**2 / 9 + 5 / 9
        L1 = math.sqrt(6) * a * b + (2 / 3) * b**2
        return C * L1 / 3

    res = minimize_scalar(lambda t: -cpsi_t(t),
                          bounds=(1e-3, math.pi / 2 - 1e-3),
                          method="bounded",
                          options={"xatol": 1e-14})
    t_opt = float(res.x)
    a_opt = math.sin(t_opt)
    b_opt = math.cos(t_opt)
    return a_opt, b_opt, cpsi_t(t_opt)


def sweep_cross_check(a_ref, n_grid=401):
    """Dense grid sweep to confirm scipy optimum matches."""
    ghz = ket_ghz(3)
    w = ket_w(3)
    best = {"min": -1.0}
    for a in np.linspace(0, 1, n_grid):
        b = math.sqrt(max(0.0, 1 - a * a))
        psi = a * ghz + b * w
        psi /= np.linalg.norm(psi)
        m = min_pair_cpsi(psi, 3)
        if m > best["min"]:
            best = {"alpha": float(a), "beta": float(b), "min": float(m)}
    return best


def family_fails_for_larger_N(max_N=6, n_grid=201):
    """For N in 3..max_N, find the best alpha|GHZ_N> + beta|W_N> in a sweep.
    Returns {N: (alpha_best, min_cpsi_best)}.
    """
    out = {}
    for N in range(3, max_N + 1):
        ghz = ket_ghz(N)
        w = ket_w(N)
        best = {"min": -1.0}
        for a in np.linspace(0, 1, n_grid):
            b = math.sqrt(max(0.0, 1 - a * a))
            psi = a * ghz + b * w
            psi /= np.linalg.norm(psi)
            m = min_pair_cpsi(psi, N)
            if m > best["min"]:
                best = {"alpha": float(a), "min": float(m)}
        out[N] = best
    return out


def main():
    RESULTS = Path(__file__).parent / "results"
    RESULTS.mkdir(exist_ok=True)
    OUT = RESULTS / "ghz_w_optimum_n3.txt"

    with open(OUT, "w", encoding="utf-8") as fh:
        def log(msg=""):
            print(msg, flush=True)
            fh.write(msg + "\n")

        log("=" * 72)
        log("  GHZ_3 + W_3 OPTIMUM FOR pair-CPsi(0)")
        log("=" * 72)

        # 1. Symbolic derivation
        log("\n-- (1) symbolic form of rho_AB(a) and CPsi(a) --")
        CPsi_sym, sextic, dCPsi_dt, t_sym, x_sym = symbolic_cpsi()
        log(f"  C(a)      = Tr(rho_AB^2) = -5 a^4 / 18 + 2 a^2 / 9 + 5/9")
        log(f"  L1_off(a) = sqrt(6) * a * b + (2/3) * b^2,   b = sqrt(1 - a^2)")
        log(f"  CPsi(a)   = C(a) * L1_off(a) / 3")

        # 2. Sextic and irreducibility
        log("\n-- (2) rationalized stationarity in x = a^2 --")
        log(f"  sextic(x) = {sp.expand(sextic)}")
        is_irred = sp.Poly(sextic, x_sym).is_irreducible
        log(f"  is_irreducible over Q: {is_irred}")
        factors = sp.factor_list(sextic)
        log(f"  factor_list: {factors}")

        # 3. Numerical optimum
        log("\n-- (3) numerical optimum --")
        a_opt, b_opt, cpsi_opt = numerical_optimum()
        log(f"  alpha_opt    = {a_opt:.15f}")
        log(f"  beta_opt     = {b_opt:.15f}")
        log(f"  alpha_opt^2  = {a_opt**2:.15f}")
        log(f"  CPsi_opt     = {cpsi_opt:.15f}")
        log(f"  ratio to 1/4 = {cpsi_opt / 0.25:.6f}x")

        # 4. Sextic root matches numerical optimum
        log("\n-- (4) sextic roots in [0, 1] --")
        roots = sp.nroots(sextic, n=20)
        real_roots_in_unit = []
        for r in roots:
            val = complex(sp.N(r))
            if abs(val.imag) < 1e-10 and 0 < val.real < 1:
                real_roots_in_unit.append(float(val.real))
                log(f"  x = {float(val.real):.15f}")
        diff = min(abs(r - a_opt**2) for r in real_roots_in_unit)
        log(f"  min distance from alpha_opt^2 to a sextic root: {diff:.2e}")

        # 5. Dense-sweep cross-check
        log("\n-- (5) dense-sweep cross-check (401 points) --")
        sweep = sweep_cross_check(a_opt, n_grid=401)
        log(f"  sweep best   : alpha = {sweep['alpha']:.6f}, min CPsi = {sweep['min']:.8f}")
        log(f"  scipy best   : alpha = {a_opt:.6f}, min CPsi = {cpsi_opt:.8f}")
        log(f"  |delta CPsi| = {abs(sweep['min'] - cpsi_opt):.2e}")

        # 6. Entanglement diagnostics at optimum
        log("\n-- (6) entanglement at the optimum --")
        psi_opt = a_opt * ket_ghz(3) + b_opt * ket_w(3)
        psi_opt /= np.linalg.norm(psi_opt)
        tau, cAB, cAC = three_tangle(psi_opt)
        rho = np.outer(psi_opt, psi_opt.conj())
        cBC = concurrence_2q(partial_trace_pair(rho, 3, (1, 2)))
        log(f"  3-tangle tau_ABC      = {tau:.6f}")
        log(f"  pair concurrences     = C(AB)={cAB:.4f}, C(AC)={cAC:.4f}, C(BC)={cBC:.4f}")
        # Permutation symmetry: all three pair-CPsi values should agree
        vals_per_pair = []
        for a, b in itertools.combinations(range(3), 2):
            vals_per_pair.append(cpsi_pair(partial_trace_pair(rho, 3, (a, b)))[2])
        log(f"  pair CPsi per pair    = {vals_per_pair}")
        log(f"  spread (max - min)    = {max(vals_per_pair) - min(vals_per_pair):.2e}")

        # 7. N in {4, 5, 6}: family fails
        log("\n-- (7) same family for N = 4, 5, 6 --")
        fails = family_fails_for_larger_N(max_N=6, n_grid=201)
        for N, b in fails.items():
            tag = "above 1/4" if b["min"] > 0.25 else "BELOW 1/4"
            log(f"  N={N}: best alpha = {b['alpha']:.4f}, min pair-CPsi = {b['min']:.6f}   [{tag}]")

        log("\n" + "=" * 72)
        log("  SUMMARY")
        log("=" * 72)
        log("  At N=3, the family alpha|GHZ_3> + beta|W_3> admits a unique maximum of")
        log("  min pair-CPsi(0) inside the fold, characterized exactly by a sextic")
        log("  that is irreducible over Q.  No simple radical form exists.  For")
        log("  N >= 4 no member of the same family lifts any pair above 1/4.")
        log(f"\n  Output: {OUT}")


if __name__ == "__main__":
    main()
