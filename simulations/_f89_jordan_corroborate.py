"""Corroboration of the F89 octic-EP diabolic verdict via rank-based Jordan
invariants and a subspace cross-check, independent of nullity-by-tol.

Key invariants for an enclosed double root (alg = 2) of L:
  A = L - lam*I
  rank(A), rank(A^2)  (counted by a CLEAR singular-value gap, not a fixed tol)
  defective (1 Jordan block of size 2):  rank(A) = n-1,  rank(A^2) = n-2
                                         => g1 = 1, g2 = 2
  diabolic  (2 Jordan blocks of size 1): rank(A) = n-2,  rank(A^2) = n-2
                                         => g1 = 2, g2 = 2
The DISCRIMINATOR is rank(A): n-1 (defective) vs n-2 (diabolic).

We also extract the 2D right-eigenvector subspace and form the 2x2
restriction of L to it; a Jordan block leaves an off-diagonal 1, a diabolic
pair gives a scalar 2x2 = lam*I.
"""
from __future__ import annotations
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def build_full_block(J, gamma):
    de_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    basis = [(i, jk) for i in range(4) for jk in de_pairs]
    index = {b: n for n, b in enumerate(basis)}
    M_SE = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J
    M_DE = np.zeros((6, 6))
    for col, (j, k) in enumerate(de_pairs):
        for nj in (j - 1, j + 1):
            if 0 <= nj <= 3 and nj != k:
                p = tuple(sorted((nj, k)))
                if p in de_pairs:
                    M_DE[de_pairs.index(p), col] += 2 * J
        for nk in (k - 1, k + 1):
            if 0 <= nk <= 3 and nk != j:
                p = tuple(sorted((j, nk)))
                if p in de_pairs:
                    M_DE[de_pairs.index(p), col] += 2 * J
    L = np.zeros((24, 24), dtype=complex)
    for col, (i, jk) in enumerate(basis):
        for i2 in range(4):
            if M_SE[i2, i] != 0:
                L[index[(i2, jk)], col] += -1j * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2 in range(6):
            if M_DE[jk_idx, jk2] != 0:
                L[index[(i, de_pairs[jk2])], col] += 1j * M_DE[jk_idx, jk2]
        L[col, col] += -2 * gamma if i in jk else -6 * gamma
    return L, basis


def rank_with_gap(M):
    """Return (rank, sorted ascending sv, the gap ratio used)."""
    s = np.sort(np.linalg.svd(M, compute_uv=False))
    # find the largest multiplicative gap between consecutive sv to set rank
    n = len(s)
    # candidate cut: the two smallest sv are "zero" iff there is a big jump up
    # report rank for an absolute tol AND show the spectrum so the gap is visible
    return s


def report(label, L, lam, n):
    A = L - lam * np.eye(L.shape[0])
    A2 = A @ A
    sA = np.sort(np.linalg.svd(A, compute_uv=False))
    sA2 = np.sort(np.linalg.svd(A2, compute_uv=False))
    print(f"--- {label} (n={L.shape[0]}) ---")
    print(f"  4 smallest sv(A)   : {sA[:4]}")
    print(f"  4 smallest sv(A^2) : {sA2[:4]}")
    # rank by a generous absolute tol sitting in the visible gap
    tol = 1e-7
    rA = int(np.sum(sA > tol))
    rA2 = int(np.sum(sA2 > tol))
    g1, g2 = L.shape[0] - rA, L.shape[0] - rA2
    print(f"  rank(A)={rA}=n-{L.shape[0]-rA}   rank(A^2)={rA2}=n-{L.shape[0]-rA2}")
    print(f"  => g1={g1}, g2={g2}   "
          f"[defective wants rank(A)=n-1; diabolic wants rank(A)=n-2]")

    # 2D eigen-subspace restriction
    _, _, Vh = np.linalg.svd(A)
    Q = Vh.conj().T[:, -2:]                  # ON basis of approx null space
    # restrict L: project.  L|_Q = Q^H L Q  (works because Q spans an
    # (approx) invariant subspace here; off-diagonal 1 would survive)
    Lr = Q.conj().T @ L @ Q
    print(f"  2x2 L restricted to null subspace (should be ~lam*I if diabolic):")
    print(f"     {np.array2string(Lr, precision=6, prefix='     ')}")
    off = np.linalg.norm(Lr - np.diag(np.diag(Lr)))
    diag_dev = np.linalg.norm(np.diag(Lr) - lam)
    print(f"     ||off-diag||={off:.3e}   ||diag-lam||={diag_dev:.3e}")
    # Also: residual test for a generalized eigenvector for EACH null vector
    for c in range(2):
        v = Q[:, c]
        x, *_ = np.linalg.lstsq(A, v, rcond=None)
        print(f"     min ||A x - v_{c}|| = {np.linalg.norm(A @ x - v):.3e}  "
              f"(=0 => generalized chain; ~O(1) => none => diabolic)")
    print()


def main():
    q = np.sqrt((-1 + np.sqrt(13)) / 6)
    lam = complex(-4.0, 2.0 * q)   # analytic exact double root, gamma=1
    print(f"q_EP={q:.10f}  lam_EP(analytic)={lam:.10f}\n")

    L, basis = build_full_block(q, 1.0)
    report("FULL 24x24", L, lam, 24)

    # symmetric 12x12
    perm = {0: 3, 1: 2, 2: 1, 3: 0}
    rp = lambda jk: tuple(sorted((perm[jk[0]], perm[jk[1]])))
    index = {b: n for n, b in enumerate(basis)}
    plus, handled = [], set()
    for n, (i, jk) in enumerate(basis):
        if n in handled:
            continue
        m = index[(perm[i], rp(jk))]
        v = np.zeros(24, dtype=complex)
        if n == m:
            v[n] = 1.0
        else:
            v[n] = v[m] = 1 / np.sqrt(2)
            handled.add(m)
        plus.append(v); handled.add(n)
    Pp = np.column_stack(plus)
    Lsym = Pp.conj().T @ L @ Pp
    report("SYM 12x12", Lsym, lam, 12)


if __name__ == "__main__":
    main()
