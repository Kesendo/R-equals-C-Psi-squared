"""F130 in the time domain: the collision decoupling as a q^4-vs-q^2 power law.

The law (docs/proofs/PROOF_F130_COLLISION_DECOUPLING.md): two distinct mode triples with
equal levels S(tau) = S(sigma) have vanishing second-order cross block B(tau, sigma) = 0.
Every prior verification is STATIC (Gram matrices). This gate pins the law's DYNAMICAL
face on the (1,2) coherence block of the full Lindbladian (XX chain, uniform local
Z-dephasing), which is exactly invariant (number-conserving H on both sides, dephasing
diagonal on |a><b|): in gamma units the block generator is

    M(q) = A + i q K,   A = -2 n_diff (rungs -2 / -6),   K = -(H2 (x) I - I (x) H1),

q = J/gamma. The basis-free observable is the exact effective coupling on the span of the
two Slater-lift multiplets [W_tau, W_sigma] (Schur complement at the cluster point
z = -6 - 2iqS; no eigenvector ambiguity inside the exactly degenerate 6-space, no
rung-decay separation). Lab = the documented off-resonance pair at n = 12 (N = 11):
tau = (1,2,10), sigma = (3,5,6), S = cos(pi/12), disjoint, eps = -1 (cell 2, the
code-trust cell; this gate is its first dynamics-side check, independent of the static
Gram chain). Conventions = resonant_n_twinning.blocks + y_zero_and_level_law lifts.

Gates (all two-sided; exit 0 iff all pass):
  G1  static pin: equal levels (1e-12); K66 W = -2S W on both multiplets (1e-12);
      cross-Gram |B(tau,sigma)| < 1e-12; generic control cross-Gram > 1e-3.
  G2  the power law: the tau-sigma Schur off-block scales as q^4 (log-log slope over
      q in [0.02, 0.32], |slope - 4| < 0.05) while the generic control scales as q^2
      (|slope - 2| < 0.05); coefficients stable (max/min < 1.05 across the window).
      The q^2 AND q^3 orders of the protected coupling are absent (the pure q^4 fit
      is the statement; a q^3 admixture would bend the slope).
  G3  the persistent degeneracy is protected to ALL orders: the X-kernel combination
      (totally antisymmetric lift) of EACH triple is an exact eigenvector of the FULL
      hop K (residual < 1e-12), while the orthogonal (3a, a) branch combinations are
      NOT (residual > 0.1) -- so the check discriminates.
  G4  the broken comb: a defect bond (delta = 0.15) destroys the eigenmultiplet
      property (residual > 0.05) -- the protection is the comb's, not generic.
  G5  self-block consistency: the tau self-block deviation from the cluster scalar
      scales as q^2 (|slope - 2| < 0.05): the coincidence is shifted apart at second
      order, so the q^4 mixing is never resonant (no avoided crossing).

Runtime a few seconds (a dozen 599-dim complex solves). Scope, honestly: N = 11 only (the
documented pair); the power law is measured, its exponent-4 mechanism (which weighted
sums beyond B vanish) is not derived here. The scout that found the law and the
decay-separation trap of naive late-time transfer ratios stayed local (WIP policy).
"""

import sys
import numpy as np
from itertools import combinations

sys.path.insert(0, "simulations")
from resonant_n_twinning import blocks, _hop  # noqa: E402
from y_zero_and_level_law import modes, slater_norm, lift  # noqa: E402

N = 11
n = N + 1
TAU = (1, 2, 10)
SIGMA = (3, 5, 6)
CTRL = (3, 4, 9)
QS = (0.02, 0.04, 0.08, 0.16, 0.32)

FAIL = []


def gate(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'}  ({detail})")
    if not ok:
        FAIL.append(name)


def level(tau):
    return sum(np.cos(k * np.pi / n) for k in tau)


def build():
    U = modes(N)
    basis = [(a, b) for a in combinations(range(N), 2) for b in combinations(range(N), 1)]
    ndiff = np.array([len(set(a) ^ set(b)) for a, b in basis])
    m6 = ndiff == 3
    H2, H1 = _hop(N, 2), _hop(N, 1)
    K = -(np.kron(H2, np.eye(N)) - np.kron(np.eye(H2.shape[0]), H1))
    A = -2.0 * ndiff
    return U, m6, A, K


def multiplet(U, tau, m6):
    nrm = slater_norm(U, tau, N)
    cols = []
    for s in range(3):
        w = np.zeros(m6.size)
        w[m6] = lift(U, tau, s, N, nrm)
        cols.append(w)
    return np.array(cols).T


def schur_off(A, K, Wt, Wp, S, q):
    Q6, _ = np.linalg.qr(np.hstack([Wt, Wp]))
    Qfull = np.linalg.qr(np.hstack([Q6, np.eye(Q6.shape[0])]))[0]
    P = Qfull[:, 6:]
    M = np.diag(A.astype(complex)) + 1j * q * K
    z = -6.0 + 1j * q * (-2.0 * S)
    Meff = (Q6.T @ M @ Q6) + (Q6.T @ M @ P) @ np.linalg.solve(
        z * np.eye(P.shape[1]) - P.T @ M @ P, P.T @ M @ Q6)
    off = np.linalg.norm(Meff[:3, 3:])
    self_dev = np.linalg.norm(Meff[:3, :3] - Meff[0, 0] * np.eye(3))
    return off, self_dev


def slope(qs, ys):
    return np.polyfit(np.log(qs), np.log(ys), 1)[0]


def main():
    U, m6, A, K = build()
    K66, K26, _ = blocks(N)
    S = level(TAU)

    print(f"n = {n}, N = {N}: tau = {TAU}, sigma = {SIGMA} (S = cos 15 deg), ctrl = {CTRL}")

    # ---- G1: the static pin
    Wt, Ws, Wc = (multiplet(U, t, m6) for t in (TAU, SIGMA, CTRL))
    dS = abs(level(TAU) - level(SIGMA))
    r_t = np.linalg.norm(K66 @ Wt[m6] - (-2.0 * S) * Wt[m6])
    r_s = np.linalg.norm(K66 @ Ws[m6] - (-2.0 * S) * Ws[m6])
    xg = np.linalg.norm((K26 @ Wt[m6]).T @ (K26 @ Ws[m6]))
    xg_c = np.linalg.norm((K26 @ Wt[m6]).T @ (K26 @ Wc[m6]))
    gate("G1", dS < 1e-12 and r_t < 1e-12 and r_s < 1e-12 and xg < 1e-12 and xg_c > 1e-3,
         f"dS = {dS:.1e}, eigen residuals {r_t:.1e}/{r_s:.1e}, "
         f"|B| = {xg:.1e}, |B_ctrl| = {xg_c:.3e}")

    # ---- G2 + G5: the power laws
    offs, offc, selfs = [], [], []
    for q in QS:
        o_s, sd = schur_off(A, K, Wt, Ws, S, q)
        o_c, _ = schur_off(A, K, Wt, Wc, S, q)
        offs.append(o_s)
        offc.append(o_c)
        selfs.append(sd)
    offs, offc, selfs = map(np.array, (offs, offc, selfs))
    qs = np.array(QS)
    sl_s, sl_c, sl_d = slope(qs, offs), slope(qs, offc), slope(qs, selfs)
    coef_s = offs / qs ** 4
    coef_c = offc / qs ** 2
    print("    q       off(protected)  off(generic)   self-dev")
    for i, q in enumerate(QS):
        print(f"    {q:5.2f}   {offs[i]:.3e}      {offc[i]:.3e}     {selfs[i]:.3e}")
    gate("G2", abs(sl_s - 4.0) < 0.05 and abs(sl_c - 2.0) < 0.05
         and coef_s.max() / coef_s.min() < 1.05 and coef_c.max() / coef_c.min() < 1.05,
         f"protected slope = {sl_s:.4f} (q^4), generic slope = {sl_c:.4f} (q^2), "
         f"coefficients {coef_s.mean():.3e} q^4 / {coef_c.mean():.3e} q^2")
    gate("G5", abs(sl_d - 2.0) < 0.05,
         f"self-block deviation slope = {sl_d:.4f} (q^2: split before the q^4 mixing)")

    # ---- G3: the all-orders protection of the antisymmetric branch (two-sided)
    ok3, det3 = True, []
    for name, tau, W in (("tau", TAU, Wt), ("sigma", SIGMA, Ws)):
        X = (K26 @ W[m6]).T @ (K26 @ W[m6])
        wX, vX = np.linalg.eigh(X)
        w_anti = W @ vX[:, 0]
        w_top = W @ vX[:, 2]
        r_anti = np.linalg.norm(K @ w_anti - (-2.0 * level(tau)) * w_anti) / np.linalg.norm(w_anti)
        r_top = np.linalg.norm(K @ w_top - (-2.0 * level(tau)) * w_top) / np.linalg.norm(w_top)
        ok3 &= r_anti < 1e-12 and r_top > 0.1
        det3.append(f"{name}: anti {r_anti:.1e}, 3a-branch {r_top:.2f}")
    gate("G3", ok3, "; ".join(det3))

    # ---- G4: the defect bond breaks the comb's protection
    delta = 0.15

    def hop_defect(k):
        st = list(combinations(range(N), k))
        idx = {s: i for i, s in enumerate(st)}
        H = np.zeros((len(st), len(st)))
        for s in st:
            occ = set(s)
            for i in range(N - 1):
                if (i in occ) ^ (i + 1 in occ):
                    H[idx[tuple(sorted(occ ^ {i, i + 1}))], idx[s]] += 1.0 + (delta if i == 0 else 0.0)
        return H

    Kd = -(np.kron(hop_defect(2), np.eye(N)) - np.kron(np.eye(len(list(combinations(range(N), 2)))), hop_defect(1)))
    K66d = Kd[np.ix_(m6.nonzero()[0], m6.nonzero()[0])]
    rd = np.linalg.norm(K66d @ Wt[m6] - (-2.0 * S) * Wt[m6])
    gate("G4", rd > 0.05, f"defect delta = {delta}: eigenmultiplet residual {rd:.3e}")

    print(f"\n{'ALL GATES PASS' if not FAIL else 'FAILED: ' + ', '.join(FAIL)}")
    sys.exit(0 if not FAIL else 1)


if __name__ == "__main__":
    main()
