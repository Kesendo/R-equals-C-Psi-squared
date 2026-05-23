"""Open question #2: which symmetry protects the 16-fold (N=6) / 12-fold (N=5)
degeneracy of each sub-cluster?

From _atmosphere_cluster_modes.py we know:
- N=6, eps=-0.83: 32 cluster modes split into two 16-mode sub-clusters by Re:
  Re=-0.3906 (XY-wt 5) and Re=-0.2094 (XY-wt 1), F1-mirror partners summing to -2sigma.
- N=5, eps=-0.3357: 24 cluster modes, two 12-mode sub-clusters at Re=-0.2697 / -0.2303.

Each sub-cluster has 16 (or 12) algebraically-degenerate modes with full geometric
rank (not an EP). The question: what symmetry algebra commutes with L and acts
within each sub-cluster, making the 16D (12D) subspace its own invariant?

This script:
  1. Build L at the dip (N=5 then N=6).
  2. Find all candidate super-operator symmetries: bit_a parity P_a, bit_b parity
     (Pi^2) P_b, spatial mirror R_super, Sz_total commutator Sz_comm, separate
     left/right Sz Sz_L / Sz_R, chiral K_super.
  3. Test [S, L] = 0 numerically: take Frobenius norm of the commutator and L.
  4. For each S that commutes with L, restrict S to the cluster subspace (32D for
     N=6, 24D for N=5) and compute eigenvalues; find the simultaneous-eigenbasis
     of all commuting Ss and L within the cluster. The dimensions of joint eigen-
     subspaces tell us the algebra responsible for the protection.

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


def sym_shape(N):
    i = np.arange(N, dtype=float)
    u = (i - (N - 1) / 2.0) ** 2
    u = u - u.mean()
    return u / np.max(np.abs(u))


def mirror_perm(N):
    """Site-mirror permutation matrix on Hilbert space (N-qubit, dim 2^N).
    Maps basis vector |b_0 b_1 ... b_{N-1}> -> |b_{N-1} ... b_1 b_0>."""
    d = 2 ** N
    R = np.zeros((d, d), dtype=complex)
    for i in range(d):
        bits = [(i >> k) & 1 for k in range(N)]
        bits_rev = bits[::-1]
        j = sum(b << k for k, b in enumerate(bits_rev))
        R[j, i] = 1.0
    return R


def chiral_op(N):
    """Chiral operator K = Z on even sites, I on odd (KHK = -H for nearest-
    neighbour XX+YY chain). Returns the Hilbert-space matrix."""
    K = np.array([[1.0 + 0j]])
    for k in range(N):
        K = np.kron(K, Z if k % 2 == 0 else I2)
    return K


def total_Sz(N):
    """Total Sz operator on 2^N Hilbert space: Sz = sum_k (1/2) Z_k."""
    d = 2 ** N
    Sz = np.zeros((d, d), dtype=complex)
    for k in range(N):
        Sz += 0.5 * site_op(Z, k, N)
    return Sz


def pauli_basis(N):
    """Returns (P shape (4^N, 4^N) with vec(sigma_alpha)/sqrt(2^N) columns
    [orthonormal in vec-space], labels)."""
    D = 4 ** N
    P = np.empty((D, D), dtype=complex)
    labels = []
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        op = np.array([[1.0 + 0j]])
        for ch in letters:
            op = np.kron(op, PAULIS[ch])
        P[:, idx] = op.flatten() / np.sqrt(2 ** N)
        labels.append(''.join(letters))
    return P, labels


def xy_weight(lbl):
    return sum(1 for c in lbl if c in 'XY')


def yz_weight(lbl):
    return sum(1 for c in lbl if c in 'YZ')


def commutator_norm(S, L):
    """||SL - LS||_F / ||L||_F. ~0 means symmetry."""
    return np.linalg.norm(S @ L - L @ S) / np.linalg.norm(L)


def build_super(A_L, B_R, d):
    """Build superoperator (A_L) o (B_R) acting as vec(rho) -> vec(A_L rho B_R^T)
    via Kron: A_L (x) B_R^conj. For real symmetries A=B=hermitian: this is
    A (x) A* in vec(rho) -> vec(A rho A^dagger) convention.
    """
    return np.kron(A_L, B_R.conj())


def main():
    for N, dip_eps, cluster_size in ((5, -0.3357, 24), (6, -0.83, 32)):
        print(f"\n========== N={N}, eps={dip_eps}, cluster_size={cluster_size} ==========\n")
        d = 2 ** N
        D = 4 ** N
        H = J * chain_H(N)
        u = sym_shape(N)
        gamma = GAMMA0 * (1.0 + dip_eps * u)

        print(f"building L ({D}x{D})...")
        L = build_L(H, gamma, N)

        print("computing eigendecomposition...")
        ev, V = np.linalg.eig(L)
        aim = np.abs(ev.imag)
        order = np.argsort(aim)
        cluster_idx = [i for i in order if aim[i] >= 1e-13][:cluster_size]
        Vc = V[:, cluster_idx]
        eigvals_c = ev[cluster_idx]
        print(f"cluster Re values: min={eigvals_c.real.min():+.4f} "
              f"max={eigvals_c.real.max():+.4f}, |Im| range: "
              f"{aim[cluster_idx[0]]:.3e} ... {aim[cluster_idx[-1]]:.3e}")
        # split into sub-clusters by Re
        re_round = np.round(eigvals_c.real, 4)
        unique_re = sorted(set(re_round.tolist()))
        sub_indices = {r: [i for i in range(cluster_size) if re_round[i] == r]
                       for r in unique_re}
        for r in unique_re:
            print(f"  Re={r:+.4f}: {len(sub_indices[r])} modes")

        # ---------- build symmetry-operator candidates ----------
        Id_d = np.eye(d, dtype=complex)
        Sz = total_Sz(N)
        R = mirror_perm(N)
        K = chiral_op(N)
        Z_total = sum(site_op(Z, k, N) for k in range(N))
        ZxN = np.array([[1.0 + 0j]])
        for _ in range(N):
            ZxN = np.kron(ZxN, Z)  # product Z⊗N

        candidates = {}

        # 1. Sz commutator: [Sz, .]
        candidates['Sz_comm'] = build_super(Sz, Id_d, d) - build_super(Id_d, Sz, d)
        # 2. Sz_L only: Sz (x) I (acts as Sz·rho)
        candidates['Sz_L'] = build_super(Sz, Id_d, d)
        # 3. Sz_R only: I (x) Sz* (acts as rho·Sz)
        candidates['Sz_R'] = build_super(Id_d, Sz, d)
        # 4. P_a = Z⊗N⊗(Z⊗N)* — bit_a parity
        candidates['P_a (ZxN⊗ZxN)'] = build_super(ZxN, ZxN, d)
        # 5. Spatial mirror super-op (R rho R^†): R⊗R*
        candidates['R_super'] = build_super(R, R, d)
        # 6. Chiral K super-op (K rho K^†): K⊗K*
        candidates['K_super'] = build_super(K, K, d)
        # 7. Chiral K-conjugation (rho -> K rho K^†): same as K_super (K hermitian)
        # 8. Z_total commutator
        candidates['Ztot_comm'] = (build_super(Z_total, Id_d, d) -
                                    build_super(Id_d, Z_total, d))
        # 9. ZxN_L only (rho -> ZxN·rho): not necessarily a symmetry
        candidates['ZxN_L'] = build_super(ZxN, Id_d, d)
        # 10. ZxN_R only
        candidates['ZxN_R'] = build_super(Id_d, ZxN, d)
        # 11. K_L only
        candidates['K_L'] = build_super(K, Id_d, d)
        # 12. K_R only
        candidates['K_R'] = build_super(Id_d, K, d)

        print(f"\n--- commutator [S, L] / |L|_F (small = symmetry) ---")
        commuting = {}
        for name, S in candidates.items():
            cn = commutator_norm(S, L)
            mark = '<-- SYM' if cn < 1e-9 else ('~zero?' if cn < 1e-6 else '')
            print(f"  {name:24s}  ||SL-LS||/|L| = {cn:.3e}  {mark}")
            if cn < 1e-6:
                commuting[name] = S

        print(f"\n--- symmetry action on cluster subspace ---")
        # Project each commuting S into the cluster eigenspace
        # S restricted to cluster: Vc^† S Vc (need Vc dual basis since V is not unitary)
        # Use pseudo-inverse: a^cluster_ij = (V^{-1})_i,c_idx · L V_c -> for symmetry
        # restriction, compute Vc^+ S Vc where Vc^+ is the left-eigenvector dual.
        # Simpler: since V is generally not unitary for non-normal L, use V_inv.
        V_inv = np.linalg.inv(V)
        # rows of V_inv corresponding to cluster_idx
        Vc_left = V_inv[cluster_idx, :]   # (cluster_size, D)
        # S restricted to cluster: M_ij = (V^{-1})_i S V_c_j  (i, j in cluster index)
        for name, S in commuting.items():
            M = Vc_left @ S @ Vc  # cluster_size x cluster_size
            # Hermitian-symmetrize for diagonalization stability
            eig_s = np.linalg.eigvals(M)
            print(f"\n  {name}: eigenvalues on cluster ({cluster_size} total):")
            # Round and count
            eig_round = np.round(eig_s, 4)
            counts = Counter([(float(np.real(x)), float(np.imag(x))) for x in eig_round])
            sorted_eigs = sorted(counts.items(),
                                 key=lambda kv: (kv[0][0], kv[0][1]))
            for (re, im), c in sorted_eigs:
                if abs(im) < 1e-6:
                    print(f"    eig = {re:+.4f}  multiplicity {c}")
                else:
                    print(f"    eig = {re:+.4f}{im:+.4f}j  multiplicity {c}")

        # ---------- simultaneous-block dimensions ----------
        # Tag each cluster mode by its eigenvalues under all commuting Ss, find
        # joint-eigenspace dimensions.
        print(f"\n--- joint-eigenspace block structure ---")
        if commuting:
            # Diagonalize each S on the cluster and assign labels.
            # For non-normal L the V may not orthogonalize S, so the eigenvalues
            # of S restricted to cluster may not be exact — round them.
            mode_tags = [[] for _ in range(cluster_size)]
            for name, S in commuting.items():
                M = Vc_left @ S @ Vc
                # Get S eigenvalues evaluated at each cluster eigenmode:
                # if mode i is a true eigenvector of S restricted to cluster,
                # then M is diagonal in that basis. In general, take diagonal
                # of M as approximate eigenvalue per cluster mode.
                diag_vals = np.diag(M)
                for i in range(cluster_size):
                    mode_tags[i].append((name, float(round(np.real(diag_vals[i]), 3)),
                                          float(round(np.imag(diag_vals[i]), 3))))
            # Cluster by tag-tuple
            tag_tuples = [tuple(t) for t in mode_tags]
            block_count = Counter(tag_tuples)
            print(f"  joint-tag block sizes: {dict(Counter(block_count.values()))}")
            print(f"  first few blocks:")
            for tag, count in list(block_count.items())[:8]:
                tag_str = ', '.join(f"{n}={r:+.2f}" for n, r, im in tag)
                print(f"    [{tag_str}]  count={count}")

    print(f"\n--- done ---")


if __name__ == "__main__":
    main()
