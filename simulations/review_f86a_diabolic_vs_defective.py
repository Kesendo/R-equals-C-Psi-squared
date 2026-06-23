"""F86a retraction review: independent diabolic-vs-defective test of the
(n,n+1)-coherence-block real-Q degeneracy.

GATE-FIRST. The load-bearing claim under review (plan 2026-06-21-f86a-ep-retraction):
the full (1,2)-coherence block's degeneracy on the REAL Q axis is DIABOLIC
(symmetry-protected; eigenvalues touch, eigenvectors stay independent;
geometric multiplicity = algebraic multiplicity), NOT a defective exceptional
point (Jordan block, eigenvectors coalesce, Petermann -> inf). The "6x Petermann"
evidence is alleged to be an `eig`-on-degenerate-eigenspace artifact.

This script settles it WITHOUT trusting either the proof or the plan, using
diagnostics that never read `eig` eigenvectors:
  (D1) geometric multiplicity  = # singular values of (L_block - lam0 I) ~ 0
  (D2) algebraic multiplicity  = trace of the Riesz spectral projector (contour integral)
  (D3) departure-from-normality of the eigenspace compression  (0 => diabolic)
  (D4) the Petermann factor computed BOTH ways: raw `eig` (artifact-prone) vs
       the stable projector route -- to reproduce/refute the alleged artifact.
diabolic  <=>  D1 == D2  and  D3 ~ 0  and  D4-stable ~ 1
defective <=>  D1 <  D2  and  D3 large and  D4 -> inf

GATE 0 validates the test on the KNOWN-defective toy 2x2 (must read DEFECTIVE),
so a 'diabolic' verdict on the full block cannot be the test being blind.
"""
from __future__ import annotations

import sys
import numpy as np

sys.path.insert(0, 'simulations')
from framework.pauli import _build_bilinear  # noqa: E402

np.set_printoptions(precision=4, suppress=True, linewidth=140)


def popcount(x):
    return bin(x).count('1')


def hamming(a, b):
    return popcount(a ^ b)


# ----------------------------------------------------------------------
# Build the (n, n+1) coherence block of L = -i[H,.] + gamma0 * sum_l (Z_l . Z_l - .)
# directly (H is 2^N x 2^N; the 4^N Liouvillian is never formed).
# ----------------------------------------------------------------------
def coherence_block(N, J, gamma0, n=1, kind='XY'):
    bonds = [(i, i + 1) for i in range(N - 1)]
    terms = {'XY': [('X', 'X', 0.5), ('Y', 'Y', 0.5)],
             'heisenberg': [('X', 'X', 1.0), ('Y', 'Y', 1.0), ('Z', 'Z', 1.0)]}[kind]
    H = J * _build_bilinear(N, bonds, terms)             # 2^N x 2^N
    A = [a for a in range(2 ** N) if popcount(a) == n]
    B = [b for b in range(2 ** N) if popcount(b) == n + 1]
    basis = [(a, b) for a in A for b in B]
    idx = {ab: i for i, ab in enumerate(basis)}
    D = len(basis)
    L = np.zeros((D, D), dtype=complex)
    for i, (a, b) in enumerate(basis):
        # -i[H, rho]:  -i H_{a'a} |a'><b|  +i H_{b b'} |a><b'|
        for ap in A:
            if abs(H[ap, a]) > 1e-15:
                L[idx[(ap, b)], i] += -1j * H[ap, a]
        for bp in B:
            if abs(H[b, bp]) > 1e-15:
                L[idx[(a, bp)], i] += 1j * H[b, bp]
        # dephasing (diagonal): -2 gamma0 * HammingDistance(a,b)
        L[i, i] += -2.0 * gamma0 * hamming(a, b)
    return L, basis


# ----------------------------------------------------------------------
# Artifact-free degeneracy diagnostics
# ----------------------------------------------------------------------
def riesz_projector(L, lam0, r, nq=400):
    """P = (1/2pi i) oint (zI - L)^{-1} dz, circle radius r about lam0."""
    D = L.shape[0]
    Id = np.eye(D, dtype=complex)
    P = np.zeros((D, D), dtype=complex)
    for k in range(nq):
        th = 2 * np.pi * (k + 0.5) / nq
        z = lam0 + r * np.exp(1j * th)
        dz = 1j * r * np.exp(1j * th) * (2 * np.pi / nq)
        P += np.linalg.solve(z * Id - L, Id) * dz
    return P / (2j * np.pi)


def geometric_multiplicity(L, lam0, tol_rel=1e-6):
    s = np.linalg.svd(L - lam0 * np.eye(L.shape[0]), compute_uv=False)
    return int(np.sum(s < tol_rel * np.max(s))), s


def characterize(L, lam0, r, label=""):
    """Return dict of artifact-free diagnostics around eigenvalue cluster lam0."""
    P = riesz_projector(L, lam0, r)
    m_alg = np.trace(P).real
    m = int(round(m_alg))
    # orthonormal basis of range(P)
    U, S, _ = np.linalg.svd(P)
    V = U[:, :max(m, 1)]
    A = V.conj().T @ L @ V                       # compression onto the eigenspace
    eigA = np.linalg.eigvals(A)
    dep = np.sqrt(max(0.0, np.linalg.norm(A, 'fro') ** 2 - np.sum(np.abs(eigA) ** 2)))
    g, svals = geometric_multiplicity(L, lam0)
    return {
        'label': label, 'lam0': lam0, 'r': r,
        'm_alg': m_alg, 'm_alg_round': m, 'm_geo': g,
        'departure_from_normality': dep,
        'A_norm': np.linalg.norm(A, 'fro'),
        'eig_spread': float(np.max(np.abs(eigA - eigA.mean()))) if m > 0 else 0.0,
        'smallest_svals': svals[-min(4, len(svals)):],
    }


def max_petermann(L):
    """MAX Petermann K over all modes, the way F86PetermannProbe does it:
    K_n = ||r_n||^2 * ||l_n||^2 with r_n = eig columns, l_n = rows of inv(R).
    Returns (K_max, eigenvalue of the argmax mode, all eigenvalues, all K)."""
    w, VR = np.linalg.eig(L)
    VR = VR / np.linalg.norm(VR, axis=0, keepdims=True)
    Rinv = np.linalg.inv(VR)
    K = (np.linalg.norm(VR, axis=0) ** 2) * (np.linalg.norm(Rinv, axis=1) ** 2)
    nmax = int(np.argmax(K))
    return float(K[nmax]), w[nmax], w, K


def projector_norm(L, lam0, r, nq=400):
    """Spectral norm of the Riesz projector onto eigenvalues within radius r of
    lam0, and the enclosed count (algebraic multiplicity). Artifact-free: it is
    a contour integral of the resolvent, never an eig eigenvector."""
    P = riesz_projector(L, lam0, r, nq)
    return np.linalg.norm(P, 2), float(np.trace(P).real)


def nearest_gap(w, lam):
    d = np.sort(np.abs(w - lam))
    return float(d[1])              # d[0] = 0 (itself)


# ----------------------------------------------------------------------
# GATE 0: the known-defective toy 2x2 MUST read DEFECTIVE.
# ----------------------------------------------------------------------
def toy_Leff(J, g_eff, gamma0, k=1):
    return np.array([[-2 * gamma0 * (2 * k - 1), 1j * J * g_eff],
                     [1j * J * g_eff, -2 * gamma0 * (2 * k + 1)]], dtype=complex)


def gate0():
    print("=" * 78)
    print("GATE 0  toy 2x2 L_eff at its EP (J*g_eff = 2 gamma0)  -- MUST be DEFECTIVE")
    print("=" * 78)
    gamma0, g_eff = 1.0, 1.0
    J = 2.0 * gamma0 / g_eff                      # exactly the EP
    L = toy_Leff(J, g_eff, gamma0, k=1)
    lam0 = -4 * gamma0 * 1                        # degenerate eigenvalue
    c = characterize(L, lam0, r=0.5, label="toy EP")
    K_eig, lam_star, _, _ = max_petermann(L)
    pn, _ = projector_norm(L, lam0, r=0.5)
    print(f"  eigenvalues          = {np.linalg.eigvals(L)}")
    print(f"  m_alg (proj trace)   = {c['m_alg']:.4f}  -> {c['m_alg_round']}")
    print(f"  m_geo (svd nullity)  = {c['m_geo']}")
    print(f"  departure-from-norm  = {c['departure_from_normality']:.4e}   (A_norm={c['A_norm']:.3f})")
    print(f"  Petermann via eig    = {K_eig:.4e}")
    print(f"  Riesz ||P|| (r=0.5)  = {pn:.4e}   (artifact-free; LARGE => genuine EP)")
    verdict = "DEFECTIVE" if (c['m_geo'] < c['m_alg_round'] and
                              c['departure_from_normality'] > 1e-3) else "diabolic"
    print(f"  >>> test reads: {verdict}   (expected DEFECTIVE)")
    assert verdict == "DEFECTIVE", "GATE 0 FAILED: test cannot detect a real EP!"
    print("  GATE 0 PASSED: the diagnostic correctly flags a genuine EP.\n")


# ----------------------------------------------------------------------
# MAIN: reproduce the F86PetermannProbe peak, then dissect it artifact-free.
# Probe's reported peaks (LocalGlobalEpLink frozen table): (N, Q, K)
# ----------------------------------------------------------------------
PROBE_PEAKS = {5: (1.288, 1333.6), 6: (0.938, 337.9), 7: (1.842, 2384.7), 8: (2.046, 795.4)}


def sweep_maxK(N, gamma0, qs, n=1):
    Ks = []
    for q in qs:
        L, _ = coherence_block(N, J=q * gamma0, gamma0=gamma0, n=n)
        Ks.append(max_petermann(L)[0])
    return np.array(Ks)


def analyze_N(N, gamma0=1.0, n=1):
    print("=" * 78)
    print(f"(n={n},{n+1}) COHERENCE BLOCK, N={N}, XY chain  -- reproduce probe peak, then dissect")
    print("=" * 78)
    # (1) reproduce the probe: sweep Q in [0.5, 4], find the max-Petermann peak.
    qs = np.linspace(0.5, 4.0, 141)
    Ks = sweep_maxK(N, gamma0, qs, n=n)
    ip = int(np.argmax(Ks))
    q_peak, K_peak = qs[ip], Ks[ip]
    ref = PROBE_PEAKS.get(N)
    refstr = f"   [probe: Q={ref[0]}, K={ref[1]}]" if ref else ""
    print(f"  K-peak on real Q axis:  Q_peak = {q_peak:.3f},  max eig-Petermann K = {K_peak:.1f}{refstr}")

    # (2) at the peak, identify the argmax-K mode and its eigenvalue lam*
    L, _ = coherence_block(N, J=q_peak * gamma0, gamma0=gamma0, n=n)
    K_eig, lam_star, w, Kall = max_petermann(L)
    gap = nearest_gap(w, lam_star)
    print(f"  argmax-K eigenvalue lam* = {lam_star:.4f}   nearest-neighbour gap = {gap:.3e}")
    print(f"    (tiny gap => lam* sits in a near-DEGENERATE cluster => eig-K is the artifact-prone regime)")

    # (3) ARTIFACT-FREE: spectral projector norm onto the cluster around lam*, several radii.
    print(f"  Riesz ||P|| onto cluster around lam* (artifact-free; O(1)=diabolic, >>1=genuine EP):")
    for r in (0.05, 0.1, 0.2, 0.4):
        pn, m = projector_norm(L, lam_star, r)
        print(f"      r={r:.2f}:  ||P|| = {pn:7.3f}   enclosed alg.mult = {m:5.2f}")

    # (4) departure-from-normality of the cluster compression (basis-independent)
    lam0 = w[np.abs(w - lam_star) < 0.15].mean()
    c = characterize(L, lam0, r=0.15, label=f"N={N}")
    relA = c['departure_from_normality'] / max(1.0, c['A_norm'])
    print(f"  cluster: m_alg={c['m_alg_round']}  m_geo(svd)={c['m_geo']}  "
          f"departure-from-normality={c['departure_from_normality']:.3e} (rel {relA:.2e})")

    # (5) grid + perturbation stability of the eig-K (artifact <=> wildly unstable)
    Kg = []
    for dq in (-1e-3, 0.0, 1e-3):
        Lq, _ = coherence_block(N, J=(q_peak + dq) * gamma0, gamma0=gamma0, n=n)
        Kg.append(max_petermann(Lq)[0])
    # tiny random perturbation at fixed Q:
    rng = np.random.default_rng(0)
    Pert = (rng.standard_normal(L.shape) + 1j * rng.standard_normal(L.shape)) * 1e-9
    K_pert = max_petermann(L + Pert)[0]
    print(f"  eig-K stability:  Q_peak-1e-3/Q_peak/+1e-3 = {Kg[0]:.1f}/{Kg[1]:.1f}/{Kg[2]:.1f};  "
          f"+1e-9 random noise -> K={K_pert:.1f}")

    pn_big, _ = projector_norm(L, lam_star, 0.2)
    artifact = (gap < 5e-2 and pn_big < 10 and relA < 1e-2)
    verdict = ("DIABOLIC near-degeneracy; eig-K is an ARTIFACT (||P||~O(1), normal eigenspace)"
               if artifact else "GENUINE non-normality (||P|| large or normality departs) -- NOT an artifact")
    print(f"  >>> verdict: {verdict}")
    print()
    return dict(N=N, q_peak=q_peak, K_peak=K_peak, gap=gap, pn02=pn_big,
                relA=relA, m_alg=c['m_alg_round'], m_geo=c['m_geo'], artifact=artifact)


if __name__ == '__main__':
    gate0()
    rows = [analyze_N(N) for N in (4, 5, 6, 7)]
    print("=" * 78)
    print("SUMMARY  (does the real-axis Petermann '6x' survive an artifact-free test?)")
    print("=" * 78)
    print(f"  {'N':>2} {'Q_peak':>7} {'K_eig':>9} {'nn-gap':>9} {'||P||(0.2)':>10} {'dep_rel':>9}  verdict")
    for r in rows:
        print(f"  {r['N']:>2} {r['q_peak']:>7.3f} {r['K_peak']:>9.1f} {r['gap']:>9.2e} "
              f"{r['pn02']:>10.3f} {r['relA']:>9.2e}  {'ARTIFACT' if r['artifact'] else 'GENUINE'}")
