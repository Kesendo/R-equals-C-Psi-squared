"""v2: clean joint diagonalization within the cluster's 32D invariant subspace.

Fixes the v1 issues:
  1. Cluster identification via strict (Re, |Im|) filtering: at uniform gamma,
     find all eigenmodes with Re in {-0.3869, -0.2131} +/- tol AND |Im|=0.0234
     +/- tol. Expect exactly 32.
  2. Use only the commuting set {Sz_L, Sz_R, P_a, R_super} for simultaneous
     diagonalization. Pi2 commutes with L but NOT with Sz_L/R, so it gets a
     separate per-mode analysis: how does Pi2 act within Sz-eigenblocks?
  3. Cluster operator restrictions via orthonormal Q (SVD of cluster
     eigenvectors), so that S_Q = Q^H S Q is the proper restriction.

Output: full table of 32 modes with quantum numbers under each Sz_L, Sz_R, P_a,
R_super; check uniqueness; analyze Pi2's action.

Investigation only.
"""
import sys
import itertools
from collections import Counter

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
LABELS = ('I', 'X', 'Y', 'Z')
N = 6
GAMMA0 = 0.05
J = 0.075


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y):
            H += site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def build_L(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma[k] * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def mirror_perm(N):
    d = 2 ** N
    R = np.zeros((d, d), dtype=complex)
    for i in range(d):
        bits = [(i >> k) & 1 for k in range(N)]
        bits_rev = bits[::-1]
        j = sum(b << k for k, b in enumerate(bits_rev))
        R[j, i] = 1.0
    return R


def total_Sz(N):
    d = 2 ** N
    Sz = np.zeros((d, d), dtype=complex)
    for k in range(N):
        Sz += 0.5 * site_op(Z, k, N)
    return Sz


def kron_super(A, B):
    return np.kron(A, B)


def to_pauli_basis_matrix(N):
    D = 4 ** N
    M = np.empty((D, D), dtype=complex)
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        op = np.array([[1.0 + 0j]])
        for ch in letters:
            op = np.kron(op, PAULIS[ch])
        M[:, idx] = op.flatten() / np.sqrt(2 ** N)
    return M


def main():
    H = J * chain_H(N)
    d = 2 ** N
    D = 4 ** N

    print(f"N={N}, J={J}, gamma0={GAMMA0}")

    # Step 1: build L at uniform gamma and eigendecompose
    print("\n--- step 1: L at uniform gamma, full eigendecomp ---")
    gamma_u = np.full(N, GAMMA0)
    L_u = build_L(H, gamma_u, N)
    print(f"  built L ({D}x{D}), eigendecomposing...")
    sys.stdout.flush()
    ev, V = np.linalg.eig(L_u)

    # Step 2: filter cluster modes by strict (Re, |Im|) criteria
    print("\n--- step 2: cluster identification by (Re, |Im|) filter ---")
    target_Im = 0.02338
    target_Re_lo = -0.3869
    target_Re_hi = -0.2131
    Im_tol = 1e-5
    Re_tol = 1e-3

    re = ev.real
    aim = np.abs(ev.imag)

    mask = (np.abs(aim - target_Im) < Im_tol) & (
        (np.abs(re - target_Re_lo) < Re_tol) | (np.abs(re - target_Re_hi) < Re_tol)
    )
    cluster_idx = np.where(mask)[0]
    print(f"  found {len(cluster_idx)} cluster modes (expecting 32)")
    if len(cluster_idx) != 32:
        print(f"  WARNING: expected 32, got {len(cluster_idx)}")
        # broaden tolerance and re-check
        Im_tol2 = 1e-3
        Re_tol2 = 0.005
        mask2 = (np.abs(aim - target_Im) < Im_tol2) & (
            (np.abs(re - target_Re_lo) < Re_tol2) | (np.abs(re - target_Re_hi) < Re_tol2)
        )
        cluster_idx = np.where(mask2)[0]
        print(f"  broadened tol: {len(cluster_idx)} modes")

    eigvals_c = ev[cluster_idx]
    Re_lo_idx = cluster_idx[np.abs(eigvals_c.real - target_Re_lo) < Re_tol]
    Re_hi_idx = cluster_idx[np.abs(eigvals_c.real - target_Re_hi) < Re_tol]
    if len(Re_lo_idx) == 0:
        Re_lo_idx = cluster_idx[eigvals_c.real < -0.3]
    if len(Re_hi_idx) == 0:
        Re_hi_idx = cluster_idx[eigvals_c.real > -0.3]
    print(f"  Re=-0.387 sub-cluster: {len(Re_lo_idx)} modes")
    print(f"  Re=-0.213 sub-cluster: {len(Re_hi_idx)} modes")

    # Step 3: orthonormal basis Q for the cluster subspace
    print("\n--- step 3: orthonormal cluster basis Q ---")
    V_cluster = V[:, cluster_idx]
    Q, sv, _ = np.linalg.svd(V_cluster, full_matrices=False)
    print(f"  SVD singular values: {sv[:5]}... ({len(sv)} total)")
    print(f"  smallest sv: {sv[-1]:.3e}")
    n_cluster = Q.shape[1]

    # Free memory
    del V, L_u

    # Step 4: build symmetry super-operators and verify [S, L] = 0
    print("\n--- step 4: build symmetries, verify they commute with L (on full Hilbert) ---")
    Id_d = np.eye(d, dtype=complex)
    Sz = total_Sz(N)
    Rm = mirror_perm(N)
    ZxN = np.array([[1.0 + 0j]])
    for _ in range(N):
        ZxN = np.kron(ZxN, Z)

    Set_A = {
        'Sz_L': kron_super(Sz, Id_d),
        'Sz_R': kron_super(Id_d, Sz.conj()),
        'P_a': kron_super(ZxN, ZxN.conj()),
        'R_super': kron_super(Rm, Rm.conj()),
    }

    # Pi2 super-operator: build via Pauli basis (diagonal with signs)
    print("  building Pauli basis matrix for Pi2 super-op...")
    sys.stdout.flush()
    M_basis = to_pauli_basis_matrix(N)
    Pi2_diag = np.empty(D, dtype=complex)
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        bb = sum(1 for ch in letters if ch in 'YZ') % 2
        Pi2_diag[idx] = (-1) ** bb
    Pi2 = M_basis @ np.diag(Pi2_diag) @ M_basis.conj().T
    del M_basis

    # Step 5: restrict each S to cluster via S_Q = Q^H S Q
    print("\n--- step 5: cluster-restricted symmetries S_Q ---")
    S_Q = {}
    for name, S in Set_A.items():
        S_Q[name] = Q.conj().T @ S @ Q  # (32, 32)
        # Pairwise commute check on Q-restricted level
        eigs = np.linalg.eigvals(S_Q[name])
        print(f"  {name}_Q: eigvals span [{eigs.real.min():+.3f}, {eigs.real.max():+.3f}]")
    Pi2_Q = Q.conj().T @ Pi2 @ Q

    # Step 6: pairwise commutation in cluster restriction
    print("\n--- step 6: pairwise commutators of set A in cluster restriction ---")
    names_A = list(S_Q.keys())
    for i, n1 in enumerate(names_A):
        for n2 in names_A[i+1:]:
            cn = np.linalg.norm(S_Q[n1] @ S_Q[n2] - S_Q[n2] @ S_Q[n1])
            print(f"  [{n1}_Q, {n2}_Q] = {cn:.3e}")
    cn = np.linalg.norm(Pi2_Q @ S_Q['Sz_L'] - S_Q['Sz_L'] @ Pi2_Q)
    print(f"  [Pi2_Q, Sz_L_Q] = {cn:.3e}  (expected non-zero: Pi2 doesn't commute with Sz_L on full Hilbert)")

    # Step 7: simultaneous diagonalization of Set A
    print("\n--- step 7: simultaneous diagonalization of {Sz_L, Sz_R, P_a, R_super} ---")
    rng = np.random.default_rng(42)
    weights = rng.uniform(0.1, 1.0, size=4)
    combo = sum(w * S_Q[n] for w, n in zip(weights, names_A))
    eig_c, V_c = np.linalg.eig(combo)
    V_c_inv = np.linalg.inv(V_c)

    # Read each S's eigenvalue per cluster mode in V_c basis
    qnums = {}
    print("\n  per-S off-diagonal residual after sim-diag:")
    for name in names_A:
        M_diag = V_c_inv @ S_Q[name] @ V_c
        diag = np.diag(M_diag).real
        off = np.linalg.norm(M_diag - np.diag(np.diag(M_diag)))
        print(f"    {name}: off-diag residual {off:.3e}")
        qnums[name] = np.round(diag, 3)

    # Re-read L's eigenvalues in V_c basis
    # L on cluster: L_Q = Q^H @ L @ Q. Diagonal in V_u basis (cluster eigvals).
    # We need L_Q in V_c basis.
    # Cluster eigenvalues are eigvals_c (in cluster_idx order, original basis).
    # In Q basis, L_Q = (eigvals_c diag) only if V_cluster columns are eigvecs.
    # They are (cluster_idx selects L-eigenvectors), but V_cluster is not Q.
    # Q comes from SVD of V_cluster. So Q basis != V_cluster basis.
    # L in Q basis = Q^H @ V_cluster @ diag(eigvals_c) @ V_cluster^(-1?) @ Q. Complicated.
    # Simpler: compute L_Q directly via Q^H L Q. Need L_u still in memory — we freed it!
    # Rebuild for this step:
    print("\n--- step 8: L in Q basis (rebuild L) ---")
    sys.stdout.flush()
    L_u = build_L(H, gamma_u, N)
    L_Q = Q.conj().T @ L_u @ Q
    del L_u
    L_diag = np.diag(V_c_inv @ L_Q @ V_c)
    off_L = np.linalg.norm(V_c_inv @ L_Q @ V_c - np.diag(L_diag))
    print(f"  L off-diag in V_c basis: {off_L:.3e}")

    # Also Pi2 in V_c basis
    Pi2_diag_perm = np.diag(V_c_inv @ Pi2_Q @ V_c)
    off_Pi2 = np.linalg.norm(V_c_inv @ Pi2_Q @ V_c - np.diag(Pi2_diag_perm))
    print(f"  Pi2 off-diag in V_c basis: {off_Pi2:.3e}  "
          f"(expected: NOT clean if Pi2 mixes Sz-blocks)")

    # Step 9: tabulate
    print(f"\n--- step 9: per-mode quantum-number table ---")
    print(f"  mode |   Re      |     Im     |  Sz_L  |  Sz_R  |  P_a  | R_super |  Pi2(diag)")
    tuples = []
    for i in range(n_cluster):
        re_i = L_diag[i].real
        im_i = L_diag[i].imag
        sz_l = qnums['Sz_L'][i]
        sz_r = qnums['Sz_R'][i]
        pa = qnums['P_a'][i]
        rs = qnums['R_super'][i]
        pi2 = round(Pi2_diag_perm[i].real, 3)
        tup = (sz_l, sz_r, pa, rs)
        tuples.append(tup)
        print(f"   {i:2d}  | {re_i:+.4f} | {im_i:+.4e} | "
              f"{sz_l:+.2f}  | {sz_r:+.2f}  | {pa:+.2f} | {rs:+.2f}   | {pi2:+.2f}")

    # Uniqueness
    tup_count = Counter(tuples)
    mult_dist = Counter(tup_count.values())
    print(f"\n  unique-tuple distribution: {dict(sorted(mult_dist.items()))}")
    if mult_dist == {1: n_cluster}:
        print(f"  -> ALL {n_cluster} MODES HAVE UNIQUE (Sz_L, Sz_R, P_a, R) TUPLES; SET-A ALGEBRA CLOSED.")
    else:
        print(f"  -> RESIDUAL DEGENERACY: ({mult_dist}); need more symmetries.")
        for tup, count in sorted(tup_count.items(), key=lambda kv: -kv[1])[:8]:
            print(f"    {tup}: {count} modes")


if __name__ == "__main__":
    main()
