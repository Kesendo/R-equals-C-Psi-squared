"""Large-N Coherence Horizon slope via SPARSE shift-invert (clock_hand_ladder thread a).

The dense SE-EP solver tops out near N=30; the asymptotic slope (dQ* still rising, ~0.602 at N=30) needs
larger N to decide exactly-2/pi vs a nearby constant. The SE Liouvillian is N^2 x N^2 but very sparse
(commutator of a tridiagonal h + a diagonal dephasing), so shift-invert near the gap (~ -2g) finds the
slowest cluster cheaply, reaching N ~ 150.

GATE FIRST: the sparse Q* must reproduce the dense ladder (N=8,12,16,20: 3.96162, 6.20222, 8.51261,
10.86307) before any large-N value is trusted. A firing gate means the sparse slow-mode capture is wrong,
not the physics."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

TWO_PI = 2.0 / np.pi

DENSE = {8: 3.961618, 12: 6.202221, 16: 8.512605, 20: 10.863071}   # from the validated dense sweep


def L_se_sparse(N, J, g):
    """Sparse N^2 x N^2 single-excitation Liouvillian, identical to coherence_horizon_se_block.L_se."""
    off = J * np.ones(N - 1)
    h = sp.diags([off, off], [1, -1], format="csc", dtype=complex)
    I = sp.identity(N, format="csc", dtype=complex)
    L = -1j * (sp.kron(h, I, format="csc") - sp.kron(I, h.T, format="csc"))
    deph = np.array([(-4.0 * g if i != j else 0.0) for i in range(N) for j in range(N)])
    return (L + sp.diags(deph, format="csc")).tocsc()


def slow_im(N, J, g, k=8):
    """max|Im| over the slowest non-zero cluster, via shift-invert. The resummed-ladder dispersion
    lambda^2 + 8g lambda + 4J^2 q^2 = 0 puts the slow mode at Re = -4g, so target sigma = -4g."""
    L = L_se_sparse(N, J, g)
    sigma = -2.0 * g   # between the kernel (0) and the slow mode (-4g): near it, but not ON it (singular)
    kk = min(k, N * N - 2)
    ncv = min(N * N - 1, max(2 * kk + 2, 40))
    vals = None
    for _ in range(4):
        try:
            vals = spla.eigs(L, k=kk, sigma=sigma, which="LM", ncv=ncv,
                             return_eigenvectors=False, maxiter=20000, tol=1e-9)
            break
        except spla.ArpackNoConvergence as e:
            if len(e.eigenvalues) >= 2:   # use the partial set that did converge
                vals = e.eigenvalues
                break
            ncv = min(N * N - 1, ncv * 2)
    if vals is None or len(vals) == 0:
        return 0.0
    nz = vals[vals.real < -1e-9]
    if len(nz) == 0:
        return 0.0
    gap = nz.real.max()
    band = nz[np.abs(nz.real - gap) < 1e-6 * max(1.0, abs(gap))]
    return float(np.abs(band.imag).max())


def qstar_sparse(N, J=1.0, iters=40):
    """Bisection on g within a Q in [0.4N, 0.85N] bracket (g* ~ J/(0.6N) sits inside)."""
    lo, hi = J / (0.85 * N), J / (0.40 * N)     # g range: small g (high Q) oscillates, large g frozen
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if slow_im(N, J, m) > 1e-7:
            lo = m                               # still oscillating -> EP at larger g
        else:
            hi = m                               # frozen -> EP at smaller g
    return J / (0.5 * (lo + hi))


def main():
    print("-- GATE: sparse vs dense --")
    ok = True
    for N, q_dense in DENSE.items():
        q = qstar_sparse(N)
        good = abs(q - q_dense) < 2e-3
        ok = ok and good
        print(f"  N={N:>2}: sparse Q*={q:.6f}  dense={q_dense:.6f}  {'ok' if good else 'MISMATCH'}")
    assert ok, "sparse solver does not reproduce the dense ladder -- fix before trusting large N"
    print("  gate passed.\n")

    print("-- large-N sweep (Q*, Q*/N, discrete slope dQ* over the step) --")
    Ns = [30, 40, 50, 60, 80, 100, 120]
    qs = {}
    prev_N, prev_q = None, None
    for N in Ns:
        q = qstar_sparse(N)
        qs[N] = q
        dq = (q - prev_q) / (N - prev_N) if prev_q is not None else float("nan")
        print(f"  N={N:>3}: Q*={q:.5f}  Q*/N={q / N:.5f}  dQ*/dN={dq:.5f}", flush=True)
        prev_N, prev_q = N, q

    # secant slope between the two largest N, and a 1/N extrapolation of Q*/N
    big = Ns[-2:]
    sec = (qs[big[1]] - qs[big[0]]) / (big[1] - big[0])
    print(f"\n  secant slope N={big[0]}->{big[1]}: {sec:.5f}")
    print(f"  2/pi = {TWO_PI:.5f}   (|gap| = {abs(sec - TWO_PI):.4f})")


if __name__ == "__main__":
    main()
