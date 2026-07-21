"""Gamma-channel rank probe: is the full-rank gamma-Jacobian caused by the palindrome?

Discriminating test for the mechanism claim in experiments/GAMMA_AS_SIGNAL.md
("the palindromic pairing creates a full-rank response matrix"). Strategy: keep
the knobs (per-site dephasing rates) and the feature set fixed, and destroy the
spectral palindrome with a spectator amplitude-damping channel (T1 breaks the
mirror; see simulations/direct_sum_scope_probe.py). If the Jacobian stays full
rank with comparable conditioning, full rank is a generic property of local
rate perturbations, not a consequence of the palindrome.

Also included: the transverse-field control (hx does NOT break the spectral
palindrome; it breaks the parity wall, not the mirror), so the field rows are
expected to be near-identical to baseline, which is a consistency check rather
than a discriminator.

Setup matches gamma_signal_analysis.py: N=5 Heisenberg chain (Pauli convention,
J=1), dissipator sum_k gamma_k (Z_k rho Z_k - rho), initial |+>^5. Features:
5 single-qubit purities + 4 nearest-neighbour CPsi + end-to-end MI, at
t = 0.5, 1.0, ..., 3.0. Jacobian d(features)/d(gamma_k) by central differences
around uniform gamma = 0.05.

Expected output (verified 2026-07-21):
    baseline Z-dephasing      rank 5/5, cond ~18, palindrome residual ~1e-14
    + transverse field hx=0.3 rank 5/5, cond ~18, palindrome residual ~1e-14
    + amplitude damping 0.05  rank 5/5, cond ~17, palindrome residual ~1e-1

Verdict: full rank survives palindrome breaking. The palindrome is not the
mechanism behind the channel's full rank.
"""

import numpy as np
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SM = (X + 1j * Y) / 2  # sigma_-

N = 5
D = 2 ** N
GAMMA0 = 0.05
EPS = 1e-4
DT = 0.5
N_TIMES = 6  # t = 0.5 .. 3.0


def op(single, site):
    ops = [I2] * N
    ops[site] = single
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def build_H(hx=0.0):
    H = np.zeros((D, D), dtype=complex)
    for i in range(N - 1):
        for P in (X, Y, Z):
            H += op(P, i) @ op(P, i + 1)
    for i in range(N):
        H += hx * op(X, i)
    return H


def build_L(H, jumps):
    Id = np.eye(D, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for r, A in jumps:
        AA = A.conj().T @ A
        L += r * (np.kron(A, A.conj())
                  - 0.5 * np.kron(AA, Id) - 0.5 * np.kron(Id, AA.T))
    return L


def reduced(rho, keep):
    dims = [2] * N
    r = rho.reshape(dims + dims)
    cur = list(range(N))
    while len(cur) > len(keep):
        t = next(i for i, s in enumerate(cur) if s not in keep)
        r = np.trace(r, axis1=t, axis2=t + len(cur))
        cur.pop(t)
    m = 2 ** len(keep)
    return r.reshape(m, m)


def vn_entropy(r):
    ev = np.linalg.eigvalsh(r)
    ev = ev[ev > 1e-12]
    return -np.sum(ev * np.log2(ev))


def features(rho_vec):
    rho = rho_vec.reshape(D, D)
    f = []
    for i in range(N):
        ri = reduced(rho, [i])
        f.append(np.real(np.trace(ri @ ri)))
    for i in range(N - 1):
        rij = reduced(rho, [i, i + 1])
        pur = np.real(np.trace(rij @ rij))
        l1 = np.sum(np.abs(rij)) - np.sum(np.abs(np.diag(rij)))
        f.append(pur * l1 / 3.0)  # CPsi with d=4
    r0, r4 = reduced(rho, [0]), reduced(rho, [N - 1])
    r04 = reduced(rho, [0, N - 1])
    f.append(vn_entropy(r0) + vn_entropy(r4) - vn_entropy(r04))
    return np.array(f)


def feature_trajectory(gammas, hx, ad):
    H = build_H(hx)
    jumps = [(g, op(Z, k)) for k, g in enumerate(gammas)]
    jumps += [(ad, op(SM, k)) for k in range(N)]
    L = build_L(H, jumps)
    E = expm(L * DT)
    plus = np.ones(2, dtype=complex) / np.sqrt(2)
    psi = plus
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    rho = np.outer(psi, psi.conj()).reshape(-1)
    out = []
    for _ in range(N_TIMES):
        rho = E @ rho
        out.append(features(rho))
    return np.concatenate(out)


def palindrome_residual(hx, ad):
    H = build_H(hx)
    jumps = [(GAMMA0, op(Z, k)) for k in range(N)]
    jumps += [(ad, op(SM, k)) for k in range(N)]
    L = build_L(H, jumps)
    re = np.sort(np.linalg.eigvals(L).real)
    # palindromic pairing: sorted real parts satisfy re[i] + re[-1-i] = -2 Sigma gamma
    return np.max(np.abs(re + re[::-1] + 2 * N * GAMMA0))


def jacobian_svd(hx, ad):
    base = np.full(N, GAMMA0)
    cols = []
    for k in range(N):
        gp = base.copy()
        gp[k] += EPS
        gm = base.copy()
        gm[k] -= EPS
        cols.append((feature_trajectory(gp, hx, ad)
                     - feature_trajectory(gm, hx, ad)) / (2 * EPS))
    return np.linalg.svd(np.column_stack(cols), compute_uv=False)


def main():
    print("Gamma-channel rank probe (N=5 Heisenberg chain, knobs = per-site dephasing)")
    print(f"{'system':<28s} {'palin.resid':>12s} {'sv_max':>8s} {'sv_min':>8s} "
          f"{'rank':>5s} {'cond':>7s}")
    for name, hx, ad in [("baseline Z-dephasing", 0.0, 0.0),
                         ("+ transverse field hx=0.3", 0.3, 0.0),
                         ("+ amplitude damping 0.05", 0.0, 0.05)]:
        resid = palindrome_residual(hx, ad)
        s = jacobian_svd(hx, ad)
        rank = int(np.sum(s > 1e-6 * s[0]))
        print(f"{name:<28s} {resid:12.3e} {s[0]:8.3f} {s[-1]:8.3f} "
              f"{rank:>4d}/5 {s[0] / s[-1]:7.1f}")
    print()
    print("Reading: amplitude damping destroys the palindrome (residual ~1e-1)")
    print("yet the gamma-Jacobian stays full rank with comparable conditioning.")
    print("Full rank is generic for independent local rate perturbations; the")
    print("palindrome is not its mechanism. The transverse-field row is a")
    print("control: hx does not break the spectral palindrome (wall, not mirror).")


if __name__ == "__main__":
    main()
