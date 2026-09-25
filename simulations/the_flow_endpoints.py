#!/usr/bin/env python3
"""The flow's two singularities, the toy 2x2 EP kept apart, and where the Kingston handover sits.

Four separately scoped objects:

  - the flow's own EP (parameter space): the single-excitation (1,1) block, H = (J/2) sum(XX+YY)
    (hop J), jump sqrt(gamma) Z on every site, Q = J/gamma. At N = 2, 3 the block's characteristic
    polynomial is factored exactly (sympy); every eigenvalue coincidence at Q > 0 is either a zero of
    one factor's discriminant or a common root of two factors (a zero of their resultant). The only
    defective one is the coherence horizon Q*(N) (Q*(2) = 1, Q*(3) = sqrt 2), certified by exact
    ranks: rank(L - lambda* I) = n^2 - 1 and rank((L - lambda* I)^2) = n^2 - 2, algebraic 2,
    geometric 1. At N = 3 the only other coincidence (Q = 4/3) is semisimple. At N = 2..5 the
    block's slowest non-kernel mode is read at 0.99 Q* and 1.01 Q* (Q* from the F2b corollary):
    real below, a complex pair above.
  - the full flow's target (state space): L is singular and the 1/N single-excitation state sits in
    its lambda = 0 kernel. The global zero eigenvalue is N+1-fold semisimple (geometric = algebraic),
    one stationary state per particle-number sector; it is not a simple eigenvalue globally.
  - the toy 2x2 EP (ExceptionalPointClock): DEFECTIVE at its own Q_EP, a separate object from the
    block; no branch continuation connects it to the flow.
  - the Kingston placement (N = 3, J = 1.5 rad/us, dt = 0.5 us, runner label Q_label = J/Gamma,
    gamma = Gamma/2, so Q_Lindblad = 2 Q_label): the continuous model's revival at the flown labels,
    the probe-time dependence of the "n0 > 0.45 at every larger Q" statistic at the k-th coherent
    return, the Q_label where the flown Trotter map's slowest mode turns oscillatory, and the Trotter
    map's own revival (exact twirl average) with the probe-time drift of the same statistic.

stdout only; no file is written.
"""
import numpy as np
import sympy as sp

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)

QSTAR_F2B = {2: 1.0, 3: float(np.sqrt(2.0)), 4: 1.87874, 5: 2.37367}   # ANALYTICAL_FORMULAS F2b corollary


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def bond_op(N, b, P, Qop):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == b else (Qop if i == b + 1 else I2))
    return o


def liouvillian_dimensionless(N, Q):
    """L' at canonical Q=J/gamma with H=(Q/2)sum(XX+YY), in column-stack vec."""
    d = 2 ** N
    Id = np.eye(d)
    H1 = sum(bond_op(N, b, X, X) + bond_op(N, b, Y, Y) for b in range(N - 1))
    L = -1j * (Q / 2.0) * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += np.kron(Zl, Zl) - np.kron(Id, Id)
    return L


def se_block(N, hop, gamma, module=np):
    """The (1,1) block on rho_ij (row-major i*N+j): -i hop [A, rho] - 4 gamma (rho - diag rho)."""
    if module is sp:
        L = sp.zeros(N * N, N * N)
    else:
        L = np.zeros((N * N, N * N), complex)
    for i in range(N):
        for j in range(N):
            r = i * N + j
            for k in (i - 1, i + 1):
                if 0 <= k < N:
                    L[r, k * N + j] += -module.I * hop if module is sp else -1j * hop
            for k in (j - 1, j + 1):
                if 0 <= k < N:
                    L[r, i * N + k] += module.I * hop if module is sp else 1j * hop
            if i != j:
                L[r, r] += -4 * gamma
    return L


def flow_ep_endpoint():
    """The flow's own EP: exact factorization, discriminants, resultants and ranks at N = 2, 3."""
    print("the flow's own EP, the (1,1)-block coherence horizon Q*(N) (gamma = 1, exact):")
    lam = sp.symbols("lambda")
    q = sp.symbols("Q", positive=True)
    for N, qstar in ((2, sp.Integer(1)), (3, sp.sqrt(2))):
        L = se_block(N, q, 1, sp)
        cp = sp.factor((lam * sp.eye(N * N) - L).det())
        factors = [f for f, _ in sp.factor_list(cp)[1]]
        print(f"  N={N}: charpoly = {cp}")
        coincidences = set()
        for f in factors:
            if sp.degree(f, lam) >= 2:
                disc = sp.factor(sp.discriminant(f, lam))
                roots = sp.solve(disc, q)
                print(f"    disc({f}) = {disc}; zeros at Q>0: {roots}")
                coincidences |= set(roots)
        for a in range(len(factors)):
            for b in range(a + 1, len(factors)):
                res = sp.factor(sp.resultant(factors[a], factors[b], lam))
                roots = sp.solve(res, q) if res.has(q) else []
                if roots:
                    print(f"    common root of {factors[a]} and {factors[b]} at Q = {roots}")
                coincidences |= set(roots)
        for qc in sorted(coincidences, key=float):
            Lq = se_block(N, qc, 1, sp)
            eig = sp.roots(sp.Poly(cp.subs(q, qc), lam))
            for ev, mult in eig.items():
                if mult < 2:
                    continue
                M = Lq - ev * sp.eye(N * N)
                geo = N * N - M.rank()
                r1, r2 = M.rank(), (M * M).rank()
                kind = "DEFECTIVE (EP2)" if geo < mult else "semisimple"
                print(f"    Q={qc}, lambda={ev}: alg {mult}, geo {geo}, rank(L-l) = {r1}, "
                      f"rank((L-l)^2) = {r2} of {N*N}: {kind}")
                if qc == qstar:
                    assert ev == -2 and mult == 2 and r1 == N * N - 1 and r2 == N * N - 2
                else:
                    assert geo == mult
        assert qstar in coincidences
    print("  slowest non-kernel mode of the block across Q* (N=2..5):")
    for N, qs in QSTAR_F2B.items():
        row = []
        for fac in (0.99, 1.01):
            w = np.linalg.eigvals(se_block(N, fac * qs, 1.0))
            w = w[np.argsort(-w.real)][1:]          # drop the kernel (the 1/N state)
            row.append(w[0])
        print(f"    N={N}: at 0.99 Q*: {row[0].real:+.4f}{row[0].imag:+.4f}i   "
              f"at 1.01 Q*: {row[1].real:+.4f}{row[1].imag:+.4f}i")


def target_endpoint():
    """The target: lambda=0 stationary state in an N+1-fold semisimple global kernel."""
    print("\nfull-flow stationary kernel:")
    for N in [3, 4, 5, 6]:
        d = 2 ** N
        L = liouvillian_dimensionless(N, 2.0)               # Q = 2
        w = np.linalg.eigvals(L)
        n_zero = int(np.sum(np.abs(w) < 1e-9))              # algebraic multiplicity of lambda=0
        ker = d * d - int(np.linalg.matrix_rank(L, tol=1e-9))  # geometric multiplicity
        rho = np.zeros((d, d), complex)
        for b in range(N):
            psi = np.array([1], complex)
            for i in range(N):
                psi = np.kron(psi, np.array([0, 1], complex) if i == b else np.array([1, 0], complex))
            rho += np.outer(psi, psi.conj()) / N
        res = float(np.linalg.norm(L @ rho.flatten(order="F")))
        kind = "SEMISIMPLE (non-defective)" if ker == n_zero else "defective"
        print(f"  N={N}: ker L = {ker} = N+1 (one fixed point per number sector);  "
              f"L.vec(1/N) = {res:.0e} (fixed point);  lambda=0 {kind} (geom {ker} = alg {n_zero})")


def toy_ep():
    """The toy 2x2 L_eff(k=1) at its own Q_EP is DEFECTIVE (rank 1 = a Jordan block)."""
    g0, k, g_eff = 1.0, 1, 4.0 / 3.0
    J = (2.0 / g_eff) * g0                                  # Q_EP = 2/g_eff
    L = np.array([[-2 * g0 * (2 * k - 1), 1j * J * g_eff],
                  [1j * J * g_eff,        -2 * g0 * (2 * k + 1)]], dtype=complex)
    lam = -4.0 * g0 * k
    M = L - lam * np.eye(2)
    rank = int(np.linalg.matrix_rank(M, tol=1e-9))
    print("\ntoy 2x2 EP (ExceptionalPointClock, a separate object):")
    print(f"  rank(L_eff - lambda_EP*I) = {rank}  -> {'DEFECTIVE (Jordan block)' if rank == 1 else 'not rank 1'}")


def kingston_placement():
    """Where the Kingston Part B handover sits relative to the walk's EP (N = 3)."""
    J, dt = 1.5, 0.5
    print("\nKingston Part B placement (N=3, J=1.5 rad/us, dt=0.5 us, Q_Lindblad = 2 Q_label):")
    print(f"  the walk's EP: Q_Lindblad* = sqrt 2, i.e. Q_label* = {np.sqrt(2.0) / 2:.4f}")

    def n0(q_lindblad, t):
        L = se_block(3, J, J / q_lindblad)
        w, V = np.linalg.eig(L)
        r0 = np.zeros(9, complex)
        r0[0] = 1.0
        return float((V @ (np.exp(w * t) * np.linalg.solve(V, r0)))[0].real)

    q_labels = [0.5, 1.0, 1.5, 2.5, 5.0, 20.0]
    rev = [max(n0(2 * ql, 2.0), n0(2 * ql, 3.0)) for ql in q_labels]
    print("  continuous-model revival max(n0(2 us), n0(3 us)) at the flown Q_label:")
    print("    " + ", ".join(f"{ql:g}: {r:.3f}" for ql, r in zip(q_labels, rev)))
    t_r = 2 * np.pi / (np.sqrt(2.0) * J)
    print(f"  first coherent return t_r = 2 pi/(sqrt2 J) = {t_r:.3f} us")
    grid = np.round(np.arange(1.0, 40.0, 0.01), 2)
    for k in (1, 2, 3):
        ok = np.array([n0(q, k * t_r) > 0.45 for q in grid])
        onset = grid[np.where(~ok)[0][-1] + 1]
        print(f"    'n0(t_{k}) > 0.45 at every larger Q' sets in at Q_Lindblad = {onset:.2f} (grid step 0.01)")

    def trotter_step(q_label):
        A = [np.zeros((3, 3)) for _ in range(2)]
        for b in range(2):
            A[b][b, b + 1] = A[b][b + 1, b] = 1.0
        U = np.eye(3, dtype=complex)
        for b in range(2):                                  # bond (0,1) first, then (1,2)
            w, V = np.linalg.eigh(A[b])
            U = (V @ np.diag(np.exp(-1j * J * dt * w)) @ V.T) @ U
        gamma_cap = J / q_label
        D = np.diag([1.0 if i == j else np.exp(-2 * gamma_cap * dt) for i in range(3) for j in range(3)])
        return D @ np.kron(U, U.conj())

    def slow_imag(q_label):
        w = np.linalg.eigvals(trotter_step(q_label))
        w = w[np.argsort(-np.abs(w))][1:]                   # drop the unit multiplier (1/N state)
        return abs(w[0].imag), w[0]

    lo, hi = 0.35, 0.37
    print(f"  Trotter map slowest multiplier at Q_label 0.35: {slow_imag(lo)[1]:.4f}, "
          f"at 0.37: {slow_imag(hi)[1]:.4f}")
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        # Close to the onset rounding can leave a real pair with |Im| up to ~sqrt(eps) ~ 1e-8; above it
        # Im grows as c*sqrt(Q - Q_c) with c ~ 0.5 (Im = 0.014 at Q_label 0.36), so a cut anywhere in
        # [1e-9, 1e-7] moves the located onset by (cut/c)^2 < 1e-13 in Q_label, below the printed 6 digits.
        if slow_imag(mid)[0] > 1e-9:
            hi = mid
        else:
            lo = mid
    print(f"  the flown Trotter map's slowest mode turns oscillatory at Q_label = {hi:.6f}")

    def trotter_n0(q_label, steps):
        r0 = np.zeros(9, complex)
        r0[0] = 1.0
        return float((np.linalg.matrix_power(trotter_step(q_label), steps) @ r0)[0].real)

    trev = [max(trotter_n0(ql, 4), trotter_n0(ql, 6)) for ql in q_labels]
    print("  flown Trotter map (exact twirl average) revival max(n0(step 4), n0(step 6)) at the flown Q_label:")
    print("    " + ", ".join(f"{ql:g}: {r:.3f}" for ql, r in zip(q_labels, trev)))
    grid_l = np.round(np.arange(0.3, 40.0, 0.01), 2)
    for steps in (6, 12, 18):
        ok = np.array([trotter_n0(q, steps) > 0.45 for q in grid_l])
        onset = grid_l[np.where(~ok)[0][-1] + 1]
        print(f"    Trotter map: 'n0(step {steps}, t = {steps * dt:g} us) > 0.45 at every larger Q' sets in at "
              f"Q_label = {onset:.2f}, Q_Lindblad = {2 * onset:.2f} (grid step 0.01)")


def main():
    print("=" * 74)
    print("THE FLOW'S TWO SINGULARITIES, THE TOY EP, AND THE KINGSTON PLACEMENT")
    print("=" * 74)
    flow_ep_endpoint()
    target_endpoint()
    toy_ep()
    kingston_placement()
    print("\nThe flow's EP is defective (alg 2, geo 1); its zero eigenvalue is N+1-fold semisimple.")
    print("The toy EP is a separate object: no branch continuation connects it to the flow.")


if __name__ == "__main__":
    main()
