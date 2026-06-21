"""F86a retraction review, part 2: where is the genuine EP, and what are the
real-axis numbers the plan calls 'diabolic'?

Confirms the CORRECT physics (vs the plan's 'diabolic / eig-artifact' story):
  (A) On the real Q axis at the Petermann peak: print the SLOWEST-mode Petermann
      (plan cites '=1.0') AND the eigenvector-matrix condition number (plan cites
      '~150') AND the gap -- to show these are the fingerprints of a NEAR-EP
      (strong, genuine non-normality concentrated at Re=-4 gamma0), not a
      diabolic (normal, ||P||~1) crossing.
  (B) Complexify Q (= complex J at fixed gamma0) and FIND the genuine defective
      EP off the real axis: the two eigenvalues nearest Re=-4 gamma0 coalesce at
      complex Q, with a true Jordan block (departure-from-normality large). This
      is the 'the genuine EP is in complex gamma' piece -- made concrete.
"""
from __future__ import annotations
import sys
import numpy as np

sys.path.insert(0, 'simulations')
from _review_f86a_diabolic_vs_defective import (  # noqa: E402
    coherence_block, riesz_projector, max_petermann)

np.set_printoptions(precision=4, suppress=True, linewidth=140)


def per_mode_petermann(L):
    w, VR = np.linalg.eig(L)
    VR = VR / np.linalg.norm(VR, axis=0, keepdims=True)
    Rinv = np.linalg.inv(VR)
    K = (np.linalg.norm(VR, axis=0) ** 2) * (np.linalg.norm(Rinv, axis=1) ** 2)
    cond = np.linalg.cond(VR)
    return w, K, cond


def real_axis_fingerprint(N, q_peak, gamma0=1.0):
    L, _ = coherence_block(N, J=q_peak * gamma0, gamma0=gamma0, n=1)
    w, K, cond = per_mode_petermann(L)
    order = np.argsort(-w.real)
    slow = order[0]
    nmax = int(np.argmax(K))
    print(f"  N={N}, Q={q_peak:.3f}:")
    print(f"     slowest mode  lam={w[slow]:.4f}  Petermann K={K[slow]:.3f}   (plan: 'slowest-mode Petermann=1.0')")
    print(f"     MAX-K  mode   lam={w[nmax]:.4f}  Petermann K={K[nmax]:.1f}   <- the non-normal near-EP mode at Re=-4")
    print(f"     eigenvector-matrix condition number = {cond:.1f}   (plan: 'cond bounded ~150'; cond>>1 => strongly NON-normal)")
    print(f"     => a diabolic crossing would have K~1 for ALL modes and cond~1; this is a NEAR-EP.")


def _min_pair_gap(N, qc, gamma0):
    L, _ = coherence_block(N, J=qc * gamma0, gamma0=gamma0, n=1)
    w = np.linalg.eigvals(L)
    D = np.abs(w[:, None] - w[None, :]) + np.eye(len(w)) * 1e9
    return float(D.min())


def find_offaxis_ep(N, q_center, gamma0=1.0):
    """Locate the nearest complex-Q point where two eigenvalues coalesce
    (Nelder-Mead on the min pairwise gap, several seeds)."""
    from scipy.optimize import minimize
    best = (1e9, None)
    seeds = [q_center + dr + 1j * di
             for dr in (-0.4, 0.0, 0.4) for di in (-1.0, -0.4, 0.4, 1.0)]
    for s in seeds:
        res = minimize(lambda x: _min_pair_gap(N, x[0] + 1j * x[1], gamma0),
                       x0=[s.real, s.imag], method='Nelder-Mead',
                       options={'xatol': 1e-7, 'fatol': 1e-12, 'maxiter': 2000})
        if res.fun < best[0]:
            best = (res.fun, res.x[0] + 1j * res.x[1])
    return best[1], best[0]


def classify_singularity(N, q_ep, gamma0=1.0):
    """Enclose the ACTUAL coalescing pair with a tight contour; report whether
    the 2x2 compression is a DEFECTIVE Jordan block (dep>0) or DIABOLIC (dep~0)."""
    L, _ = coherence_block(N, J=q_ep * gamma0, gamma0=gamma0, n=1)
    w = np.linalg.eigvals(L)
    D = np.abs(w[:, None] - w[None, :]) + np.eye(len(w)) * 1e9
    i, j = np.unravel_index(np.argmin(D), D.shape)
    lam0 = 0.5 * (w[i] + w[j])
    coalescence = abs(w[i] - w[j])
    third = np.sort(np.abs(w - lam0))[2]            # distance to nearest OTHER eigenvalue
    r = max(3 * coalescence, 0.2 * third)
    r = min(r, 0.6 * third)
    P = riesz_projector(L, lam0, r=r)
    m = np.trace(P).real
    U, S, _ = np.linalg.svd(P)
    V = U[:, :2]
    A = V.conj().T @ L @ V
    eigA = np.linalg.eigvals(A)
    dep = np.sqrt(max(0.0, np.linalg.norm(A, 'fro') ** 2 - np.sum(np.abs(eigA) ** 2)))
    pn = np.linalg.norm(P, 2)
    kind = "DEFECTIVE EP (Jordan block)" if dep > 1e-2 else "diabolic (normal crossing)"
    print(f"  N={N}: nearest coalescence at  Q = {q_ep.real:+.4f} {q_ep.imag:+.4f}i   "
          f"(|Im Q| = {abs(q_ep.imag):.3f}{'  OFF real axis' if abs(q_ep.imag) > 1e-3 else ''})")
    print(f"        coalescing pair gap = {coalescence:.2e}   enclosed alg.mult = {m:.2f}   ||P|| = {pn:.2f}")
    print(f"        departure-from-normality of the 2x2 = {dep:.4f}   =>  {kind}")


if __name__ == '__main__':
    print("=" * 78)
    print("(A) REAL-AXIS FINGERPRINT at the Petermann peak  -- sourcing the plan's 'diabolic' numbers")
    print("=" * 78)
    for N, q in [(5, 2.150), (6, 0.925), (7, 1.850)]:
        real_axis_fingerprint(N, q)
        print()

    print("=" * 78)
    print("(B) THE GENUINE EP IS OFF THE REAL AXIS (complex Q = complex J / gamma0)")
    print("=" * 78)
    for N, qc in [(5, 2.15), (6, 0.925), (7, 1.85)]:
        q_ep, gap = find_offaxis_ep(N, qc)
        classify_singularity(N, q_ep)
        print()
