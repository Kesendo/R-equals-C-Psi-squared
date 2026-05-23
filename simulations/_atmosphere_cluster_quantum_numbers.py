"""Tests beta + gamma: simultaneous diagonalization within the cluster at uniform
gamma, with full identification of each mode's symmetry quantum numbers.

Goal: assign each of the 32 cluster modes a unique tuple (Sz_L, Sz_R, P_a, R)
or similar. If the assignment is bijective (32 distinct tuples), the listed
symmetries fully classify the cluster. If multiple modes share a tuple, there
is residual unexplained degeneracy.

Method:
  1. Build L at uniform gamma (N=6). Eigendecompose.
  2. Identify 32 cluster modes by max-overlap with the dip-cluster (eps=-0.83).
     This re-uses the protocol from cluster_modes_at_uniform.
  3. Build candidate symmetry super-operators: Sz_L, Sz_R, P_a (Z(x)N tensor),
     R_super (spatial mirror), Pi_squared (bit_b parity super-op).
  4. Verify each commutes with L (numerical [S, L] = 0).
  5. Verify each pair commutes ([S_i, S_j] = 0 within the cluster).
  6. Restrict each S to cluster: M_S = V_left @ S @ V_right.
  7. Take a random linear combination of M_S's, diagonalize -> get a basis
     in which all M_S are simultaneously (almost) diagonal.
  8. Read each cluster mode's eigenvalue under each S.
  9. Count unique tuples vs degenerate tuples.

If we find residual degeneracy, we need an additional symmetry; that becomes
the next question. If we find unique tuples, the algebra is closed.

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
CLUSTER_SIZE = 32


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
    return np.kron(A, B)


def pi_squared_super(N):
    """Pi^2 super-operator: acts on Pauli string sigma_alpha as (-1)^bit_b(alpha)
    where bit_b = #Y + #Z mod 2. This is NOT a simple kron product.

    For each computational-basis transition |i><j|: bit_b changes by ... complex.
    Instead, build via Pauli basis: P_b = sum_alpha (-1)^bit_b(alpha) |sigma_alpha><sigma_alpha|

    On vec(rho) basis: P_b vec(sigma_alpha) = (-1)^bit_b(alpha) vec(sigma_alpha).
    So in vec basis, P_b is diagonal in the Pauli-vector basis.
    """
    D = 4 ** N
    diag_signs = np.empty(D, dtype=complex)
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        bb = sum(1 for ch in letters if ch in 'YZ') % 2
        diag_signs[idx] = (-1) ** bb
    return diag_signs  # length-D vector representing diagonal in Pauli basis


def to_pauli_basis_matrix(N):
    """Returns matrix M of shape (4^N, 4^N) such that M[:, alpha] = vec(sigma_alpha)/sqrt(2^N).
    M is unitary (Pauli basis orthonormal in trace inner product / 2^N)."""
    D = 4 ** N
    M = np.empty((D, D), dtype=complex)
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        op = np.array([[1.0 + 0j]])
        for ch in letters:
            op = np.kron(op, PAULIS[ch])
        M[:, idx] = op.flatten() / np.sqrt(2 ** N)
    return M


def commutator_norm(S, L):
    return np.linalg.norm(S @ L - L @ S) / max(np.linalg.norm(L), 1e-30)


def main():
    H = J * chain_H(N)
    u = sym_shape(N)
    d = 2 ** N
    D = 4 ** N

    print(f"N={N}, J={J}, gamma0={GAMMA0}")

    # Step 1: cluster at eps=-0.83
    print("\n--- step 1: cluster at eps=-0.83 (dip), |Im|-smallest 32 modes ---")
    gamma_dip = GAMMA0 * (1.0 + (-0.83) * u)
    L_dip = build_L(H, gamma_dip, N)
    ev_dip, V_dip = np.linalg.eig(L_dip)
    aim = np.abs(ev_dip.imag)
    nonreal = [i for i in np.argsort(aim) if aim[i] >= 1e-13]
    idx_dip = np.array(nonreal[:CLUSTER_SIZE])
    Q_dip, _, _ = np.linalg.svd(V_dip[:, idx_dip], full_matrices=False)
    print(f"  cluster <Re>={np.mean(ev_dip[idx_dip].real):+.4f}, "
          f"<|Im|>={np.mean(np.abs(ev_dip[idx_dip].imag)):.4e}")
    del L_dip, V_dip

    # Step 2: cluster at uniform gamma via overlap with Q_dip
    print("\n--- step 2: cluster at uniform gamma via Q_dip overlap ---")
    gamma_u = np.full(N, GAMMA0)
    L_u = build_L(H, gamma_u, N)
    ev_u, V_u = np.linalg.eig(L_u)
    overlaps = np.linalg.norm(Q_dip.conj().T @ V_u, axis=0)
    V_norms = np.linalg.norm(V_u, axis=0)
    overlaps_n = overlaps / np.maximum(V_norms, 1e-30)
    idx_u = np.argsort(-overlaps_n)[:CLUSTER_SIZE]
    eigvals_c = ev_u[idx_u]
    print(f"  cluster <Re>={np.mean(eigvals_c.real):+.4f}, "
          f"<|Im|>={np.mean(np.abs(eigvals_c.imag)):.4e}")
    Vc_right = V_u[:, idx_u]

    # Get left-eigenvectors (rows of V_u^-1) for restriction
    print(f"  computing V^-1 (D={D})...")
    sys.stdout.flush()
    V_u_inv = np.linalg.inv(V_u)
    Vc_left = V_u_inv[idx_u, :]
    del V_u, V_u_inv

    # Step 3: build candidate symmetry super-operators
    print("\n--- step 3: building symmetry candidates ---")
    Id_d = np.eye(d, dtype=complex)
    Sz = total_Sz(N)
    Rm = mirror_perm(N)
    ZxN = np.array([[1.0 + 0j]])
    for _ in range(N):
        ZxN = np.kron(ZxN, Z)

    candidates = {
        'Sz_L': kron_super(Sz, Id_d),
        'Sz_R': kron_super(Id_d, Sz.conj()),
        'P_a': kron_super(ZxN, ZxN.conj()),
        'R_super': kron_super(Rm, Rm.conj()),
    }
    # Pi_squared super-op: diagonal in Pauli basis. Convert to vec(rho) basis.
    print(f"  building Pauli basis matrix M_basis ({D}x{D})...")
    sys.stdout.flush()
    M_basis = to_pauli_basis_matrix(N)
    Pi2_diag_pauli = pi_squared_super(N)
    # Pi2 in vec basis: M_basis @ diag(Pi2_diag_pauli) @ M_basis^H
    Pi2 = M_basis @ np.diag(Pi2_diag_pauli) @ M_basis.conj().T
    candidates['Pi2'] = Pi2
    del M_basis

    # Step 4: verify each commutes with L
    print("\n--- step 4: [S, L] commutator norms ---")
    Lnorm = np.linalg.norm(L_u)
    for name, S in candidates.items():
        cn = np.linalg.norm(S @ L_u - L_u @ S) / Lnorm
        mark = '<-- SYM' if cn < 1e-9 else '<-- NOT SYM'
        print(f"  [{name}, L] / |L| = {cn:.3e}  {mark}")
    del L_u

    # Step 5: pairwise commutators among candidates
    print("\n--- step 5: pairwise commutator norms among candidates ---")
    names = list(candidates.keys())
    for i, n1 in enumerate(names):
        for n2 in names[i+1:]:
            S1, S2 = candidates[n1], candidates[n2]
            cn = np.linalg.norm(S1 @ S2 - S2 @ S1) / max(np.linalg.norm(S1), 1)
            print(f"  [{n1}, {n2}] / |{n1}| = {cn:.3e}")

    # Step 6: restrict each S to cluster subspace
    print("\n--- step 6: cluster-restricted symmetries M_S = V_left @ S @ V_right ---")
    M_ops = {}
    for name, S in candidates.items():
        M_ops[name] = Vc_left @ S @ Vc_right  # (32, 32)
        print(f"  {name}: M_S shape {M_ops[name].shape}, "
              f"eigvals span [{np.linalg.eigvals(M_ops[name]).real.min():.2f}, "
              f"{np.linalg.eigvals(M_ops[name]).real.max():.2f}]")

    # Step 7: simultaneous diagonalization via random linear combination
    print("\n--- step 7: simultaneous diagonalization via random combination ---")
    rng = np.random.default_rng(42)
    # Random weights, each in [0.1, 1.0]
    weights = rng.uniform(0.1, 1.0, size=len(M_ops))
    combo = sum(w * M for w, M in zip(weights, M_ops.values()))
    # Make Hermitian-symmetrize for better numerical behavior
    print(f"  diagonalizing random combo (32x32)...")
    eig_combo, V_combo = np.linalg.eig(combo)
    # V_combo columns are now simultaneous-eigenvectors (approximately) of all M_S

    # Step 8: read off each cluster mode's eigenvalue under each S
    print(f"\n--- step 8: per-mode quantum numbers ---")
    V_combo_inv = np.linalg.inv(V_combo)
    quantum_numbers = {}  # name -> array of 32 eigenvalues
    for name, M_S in M_ops.items():
        M_diag = V_combo_inv @ M_S @ V_combo
        diag = np.diag(M_diag).real
        # Round to handle numerical noise
        diag_round = np.round(diag, 2)
        quantum_numbers[name] = diag_round
        # Also check off-diagonal residual to verify simultaneity
        off_diag = np.linalg.norm(M_diag - np.diag(np.diag(M_diag))) / max(np.linalg.norm(M_diag), 1e-30)
        print(f"  {name}: off-diagonal residual after sim-diag = {off_diag:.3e}")

    # Step 9: tabulate per-mode tuples and count uniqueness
    print(f"\n--- step 9: per-mode quantum-number tuples ---")
    names_ordered = list(M_ops.keys())
    tuples = []
    for i in range(CLUSTER_SIZE):
        # Also get L-eigenvalue (Re only, since |Im| degenerate)
        # L-eigenvalue in combo basis: diag(V_combo^-1 L_cluster V_combo)
        tup = tuple(quantum_numbers[name][i] for name in names_ordered)
        tuples.append(tup)

    # Print all 32 modes' tuples
    print(f"  (per-mode) Re | " + " | ".join(f"{n:>8s}" for n in names_ordered))
    # Get L-eigenvalues restricted to cluster basis
    L_eigs_combo = np.diag(V_combo_inv @ (Vc_left @ Vc_right) @ V_combo)  # this is just identity restricted; we need L on cluster
    # Actually L-eigenvalue per cluster mode (in V_combo basis):
    # L_cluster (in V_u-basis) is diagonal with eigvals_c on cluster modes
    # Translate to V_combo basis: M_L = V_combo^-1 @ diag(eigvals_c) @ V_combo (since cluster is already diagonal in original V basis)
    L_cluster_diag = np.diag(eigvals_c)
    M_L = V_combo_inv @ L_cluster_diag @ V_combo
    L_per_mode = np.diag(M_L)
    for i in range(CLUSTER_SIZE):
        tup_str = " | ".join(f"{v:+8.3f}" for v in tuples[i])
        print(f"  mode {i:2d}: Re={L_per_mode[i].real:+.4f} (Im={L_per_mode[i].imag:+.4e}) | {tup_str}")

    # Count unique tuples
    tup_count = Counter(tuples)
    multiplicity_dist = Counter(tup_count.values())
    print(f"\n  unique-tuple distribution: {dict(sorted(multiplicity_dist.items()))}")
    if multiplicity_dist == {1: CLUSTER_SIZE}:
        print(f"  -> ALL 32 MODES HAVE UNIQUE QUANTUM-NUMBER TUPLES; ALGEBRA CLOSED.")
    else:
        print(f"  -> RESIDUAL DEGENERACY: some tuples appear multiple times.")
        print(f"  Most-degenerate tuples:")
        for tup, count in sorted(tup_count.items(), key=lambda kv: -kv[1])[:5]:
            tup_str = ", ".join(f"{n}={v:+.2f}" for n, v in zip(names_ordered, tup))
            print(f"    [{tup_str}]: {count} modes")


if __name__ == "__main__":
    main()
