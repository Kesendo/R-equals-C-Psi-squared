"""OQ#2 follow-up: joint diagonalization within each sub-cluster.

The first sweep (_atmosphere_symmetry_protection.py) found:
- N=6: all 32 cluster modes have P_a = -1; the Re-split (16+16 at -0.3906/-0.2094)
  is NOT bit_a-parity-based.
- N=5: 12+12 P_a-split exists; the relation to the Re-split (12+12 at -0.2697/-0.2303)
  is unknown.
- R_super splits N=5 cluster 16/8 (asymmetric); N=6 cluster 16/16.

This script splits the cluster by L-eigenvalue first (Re-sub-clusters), then within
each Re-sub-cluster:
  1. Compute the cluster-restricted action of each commuting symmetry.
  2. Diagonalize all of them simultaneously (find joint-eigenspace).
  3. Report dimensions of joint-eigenspaces.

If the joint-eigenspace dimensions are 1 (or close to 1) for every (sub-cluster,
symmetry-tuple) combination, then the listed symmetries fully classify the modes
and explain the degeneracy. If joint blocks are larger than 1, there's a residual
symmetry not in our list.

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
    """vec(rho) -> vec(A rho B^T) is A (x) B."""
    return np.kron(A, B)


def main():
    for N, dip_eps, cluster_size in ((5, -0.3357, 24), (6, -0.83, 32)):
        print(f"\n========== N={N}, eps={dip_eps}, cluster_size={cluster_size} ==========\n")
        d = 2 ** N
        D = 4 ** N
        H = J * chain_H(N)
        u = sym_shape(N)
        gamma = GAMMA0 * (1.0 + dip_eps * u)
        print(f"building L ({D}x{D}) and diagonalizing...")
        L = build_L(H, gamma, N)
        ev, V = np.linalg.eig(L)
        aim = np.abs(ev.imag)
        order = np.argsort(aim)
        cluster_idx = np.array([i for i in order if aim[i] >= 1e-13][:cluster_size])
        eigvals_c = ev[cluster_idx]
        # split into sub-clusters by Re
        re_round = np.round(eigvals_c.real, 4)
        sub_re_values = sorted(set(re_round.tolist()))
        print(f"Re sub-clusters: {sub_re_values}")
        sub_groups = {r: cluster_idx[re_round == r] for r in sub_re_values}
        for r, idxs in sub_groups.items():
            print(f"  Re={r:+.4f}: {len(idxs)} modes")

        # candidate symmetries (those commuting with L)
        Id_d = np.eye(d, dtype=complex)
        Sz = total_Sz(N)
        R = mirror_perm(N)
        ZxN = np.array([[1.0 + 0j]])
        for _ in range(N):
            ZxN = np.kron(ZxN, Z)
        sym_ops = {
            'Sz_L': kron_super(Sz, Id_d),
            'Sz_R': kron_super(Id_d, Sz.conj()),
            'P_a': kron_super(ZxN, ZxN.conj()),
            'R_super': kron_super(R, R.conj()),
        }
        # Verify these commute with L
        Lnorm = np.linalg.norm(L)
        print("verifying [S, L] = 0:")
        for name, S in sym_ops.items():
            cn = np.linalg.norm(S @ L - L @ S) / Lnorm
            print(f"  {name}: {cn:.3e}")

        # For each Re sub-cluster, do joint diagonalization of all symmetries
        for r, idxs in sub_groups.items():
            print(f"\n-- Re sub-cluster {r:+.4f} ({len(idxs)} modes) --")
            # Extract V columns and V_inv rows for this sub-cluster
            V_sub = V[:, idxs]  # (D, n)
            V_inv = np.linalg.inv(V)
            Vleft_sub = V_inv[idxs, :]  # (n, D)
            n = len(idxs)

            # For each symmetry, restrict to sub-cluster: M_S = Vleft_sub @ S @ V_sub
            M_ops = {}
            for name, S in sym_ops.items():
                M_ops[name] = Vleft_sub @ S @ V_sub  # (n, n)

            # Simultaneous diagonalization: diagonalize M_R_super first (small commuting
            # algebra; pick one that splits the most). Then within each block, diagonalize
            # M_Sz_L, etc.
            # Random linear combination to break ties:
            random_combo = (0.7 * M_ops['Sz_L'] + 0.5 * M_ops['Sz_R']
                            + 0.3 * M_ops['P_a'] + 0.2 * M_ops['R_super'])
            eig_combo, V_combo = np.linalg.eig(random_combo)
            # Now in V_combo basis, each M should be (approximately) diagonal
            # Tag each mode by its eigenvalues under each symmetry:
            tags_per_mode = []
            for i in range(n):
                tag = []
                for name in ('P_a', 'R_super', 'Sz_L', 'Sz_R'):
                    M_diag = V_combo @ np.diag(np.linalg.solve(V_combo, M_ops[name] @ V_combo[:, i] / np.linalg.norm(V_combo[:, i]))).real if False else (
                        np.linalg.solve(V_combo, M_ops[name] @ V_combo)
                    )
                    # That formula is wrong; do the right thing:
                    pass
                # Better: explicitly transform each M into the V_combo basis
                # and read the diagonal
                # (move this outside the per-mode loop)
            # Cleaner: transform M into V_combo basis once, take diagonal
            for name in ('P_a', 'R_super', 'Sz_L', 'Sz_R'):
                M_in_combo = np.linalg.solve(V_combo, M_ops[name] @ V_combo)
                diag = np.diag(M_in_combo)
                # Round to handle numerical noise
                diag_round = np.round(diag.real, 3) + 1j * np.round(diag.imag, 3)
                if name == 'P_a':
                    print(f"  {name}: {sorted(Counter(diag_round.real).items())}")
                else:
                    rd = sorted(Counter(diag_round.real).items())
                    print(f"  {name}: {rd}")

            # Joint tagging
            joint_tags = []
            tag_per_mode_dict = {name: [] for name in ('P_a', 'R_super', 'Sz_L', 'Sz_R')}
            for name in ('P_a', 'R_super', 'Sz_L', 'Sz_R'):
                M_in_combo = np.linalg.solve(V_combo, M_ops[name] @ V_combo)
                diag = np.diag(M_in_combo)
                diag_round = tuple(np.round(d_val.real, 2) for d_val in diag)
                tag_per_mode_dict[name] = diag_round
            for i in range(n):
                joint_tags.append(tuple(tag_per_mode_dict[name][i]
                                        for name in ('P_a', 'R_super', 'Sz_L', 'Sz_R')))
            block_count = Counter(joint_tags)
            block_sizes = Counter(block_count.values())
            print(f"  joint (P_a, R_super, Sz_L, Sz_R) block sizes: "
                  f"{dict(sorted(block_sizes.items()))}")
            for tag, c in sorted(block_count.items()):
                print(f"    (P_a={tag[0]:+.0f}, R={tag[1]:+.2f}, "
                      f"Sz_L={tag[2]:+.2f}, Sz_R={tag[3]:+.2f}): count={c}")

    print(f"\n--- done ---")


if __name__ == "__main__":
    main()
