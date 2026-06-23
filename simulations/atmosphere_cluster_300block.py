"""Direct examination of L's structure in the (Sz_L=0, Sz_R=+1) sector at
uniform gamma, N=6.

In H-eigenstate basis, the operator basis sigma_{ij} = |phi_i><phi_j| is 20*15
= 300 dimensional (Sz=0 H-eigenstates × Sz=+1 H-eigenstates).

The H-commutator part of L is diagonal with values -i(lambda_i - lambda_j).
The dissipator D is off-diagonal in this basis but preserves Sz.

Aim: find L's eigenvalues in this block, identify the 4 cluster modes at
|Im|=0.0234, and explicitly see what symmetry pairs the doublets.

Investigation only.
"""
import sys
import math
import itertools
from collections import Counter

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
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


def total_Sz(N):
    d = 2 ** N
    Sz = np.zeros((d, d), dtype=complex)
    for k in range(N):
        Sz += 0.5 * site_op(Z, k, N)
    return Sz


def mirror_perm(N):
    d = 2 ** N
    R = np.zeros((d, d), dtype=complex)
    for i in range(d):
        bits = [(i >> k) & 1 for k in range(N)]
        bits_rev = bits[::-1]
        j = sum(b << k for k, b in enumerate(bits_rev))
        R[j, i] = 1.0
    return R


def main():
    d = 2 ** N
    H = J * chain_H(N)
    Sz = total_Sz(N)
    R_op = mirror_perm(N)

    # H and Sz_total commute (XX+YY conserves Sz). Diagonalize H within each Sz sector.
    print(f"N={N}, J={J}, gamma0={GAMMA0}")
    print(f"\n--- step 1: diagonalize H in each Sz sector ---")
    # First diagonalize Sz to get sector projectors
    sz_eig, V_sz = np.linalg.eigh(Sz)
    # Each Sz eigenvalue (half-integer) corresponds to a sector
    sz_round = np.round(2 * sz_eig).astype(int) / 2  # round to half-integers
    sectors = {}
    for sz in sorted(set(sz_round.tolist())):
        idx = np.where(sz_round == sz)[0]
        sectors[sz] = V_sz[:, idx]  # basis of this Sz sector
        print(f"  Sz={sz:+.1f}: dim {len(idx)}")

    # Restrict H to each sector and diagonalize
    H_eigs = {}  # Sz -> (energies, eigvecs in full Hilbert space)
    for sz, V_sec in sectors.items():
        H_sec = V_sec.conj().T @ H @ V_sec
        # Make hermitian (should already be)
        H_sec = 0.5 * (H_sec + H_sec.conj().T)
        eigs_sec, U_sec = np.linalg.eigh(H_sec)
        full_vecs = V_sec @ U_sec  # H-eigenstates in full Hilbert basis
        H_eigs[sz] = (eigs_sec, full_vecs)
        print(f"  Sz={sz:+.1f}: H eigenvalues span [{eigs_sec.min():+.4f}, "
              f"{eigs_sec.max():+.4f}]")

    # Step 2: build operator basis sigma_{ij} = |phi_i_Sz0><phi_j_Sz+1| for
    # bra in Sz=0 (20 states), ket in Sz=+1 (15 states). 20*15 = 300 operators.
    print(f"\n--- step 2: build sigma_{{ij}} basis for (Sz_L=0, Sz_R=+1) ---")
    eigs_b, phi_b = H_eigs[0.0]  # bra Sz=0: 20 states
    eigs_k, phi_k = H_eigs[1.0]  # ket Sz=+1: 15 states
    nb = len(eigs_b)
    nk = len(eigs_k)
    n_ops = nb * nk
    print(f"  bra dim {nb}, ket dim {nk}, operator basis dim {n_ops}")
    assert n_ops == 300

    # Build sigma_{ij} as 64x64 matrices, then vectorize
    # Actually let's just work in (i, j) index space directly for the block matrix.

    # L block structure:
    #   H-commutator part: L_H[ij,kl] = -i(lambda_i - lambda_j) delta_ik delta_jl
    #   Dissipator part: D[ij, kl] = gamma sum_m <phi_i|Z_m|phi_k> <phi_l|Z_m|phi_j> - delta_ik delta_jl
    # where Z_m acts on bra-side (gives <phi_i|Z_m|phi_k> on bra index) and ket-side
    # (gives <phi_l|Z_m|phi_j> on ket index).

    # Compute Z_m matrix elements in H-eigenstate basis
    print(f"  computing Z_m matrix elements in H-eigenstate basis...")
    Z_bra_bra = []  # list of 20x20 matrices, one per site
    Z_ket_ket = []  # list of 15x15 matrices, one per site
    for m in range(N):
        Z_m = site_op(Z, m, N)
        Z_bb = phi_b.conj().T @ Z_m @ phi_b
        Z_kk = phi_k.conj().T @ Z_m @ phi_k
        Z_bra_bra.append(Z_bb)
        Z_ket_ket.append(Z_kk)

    # Build L block
    print(f"  building L block ({n_ops}x{n_ops})...")
    L_block = np.zeros((n_ops, n_ops), dtype=complex)
    # H-commutator part (diagonal)
    for i in range(nb):
        for j in range(nk):
            idx = i * nk + j
            L_block[idx, idx] += -1j * (eigs_b[i] - eigs_k[j])
    # Dissipator part: for each site m, add gamma * (Z_m sigma Z_m - sigma)
    # Z_m sigma_{ij} Z_m = sum_{kl} Z_bb[k,i] * Z_kk[j,l] * sigma_{kl}^T_or_similar?
    # Careful: sigma_{ij} = |phi_i^b><phi_j^k|. Z_m sigma Z_m = Z_m |phi_i^b><phi_j^k| Z_m
    # = (Z_m |phi_i^b>)(<phi_j^k| Z_m) = (sum_k Z_bb[k,i] |phi_k^b>)(sum_l Z_kk[j,l]* <phi_l^k|)
    # = sum_{k,l} Z_bb[k,i] Z_kk[j,l]^* sigma_{kl}
    # But Z is real (diagonal {+/-1}), so Z_kk[j,l]^* = Z_kk[j,l]. Z is hermitian, so
    # Z_kk[j,l] = Z_kk[l,j]^*. For Z hermitian + real, Z_kk[j,l] = Z_kk[l,j].
    # In matrix form: D_m[sigma]_(kl) = (Z_bb @ sigma @ Z_kk^T)_(kl) where sigma is the
    # 20x15 matrix of coefficients. Equivalently in vectorized form: D_m sigma_vec =
    # (Z_kk^T (x) Z_bb) sigma_vec (Kron with bra on outer index).
    for m in range(N):
        # In the basis ordering idx = i * nk + j (bra-major), the action is
        # L[i' j', i j] += gamma * (Z_bb[i', i] * Z_kk[j, j'] - delta_ii' delta_jj')
        # Wait need to be careful with index conventions.
        # sigma_{ij} = sum_{i', j'} D_m[(i'j'), (ij)] * sigma_{i' j'}
        # where Z_m sigma_{ij} Z_m = sum_{i', j'} Z_bb[i', i] Z_kk[j, j'] sigma_{i' j'}
        # (using <phi_l^k| Z_m = sum_l Z_kk[j, l] <phi_l^k| -- since Z_kk[j,l] = <phi_j| Z |phi_l>)
        # So D_m[(i'j'),(ij)] = Z_bb[i',i] * Z_kk[j, j']
        # Index: row = i'*nk + j', col = i*nk + j
        # As outer product: D_m = Z_bb (x_row) ZkT_along_col? Let me just use kron.
        # D_m = Z_bb (kron) Z_kk^T (with Z_kk^T meaning transpose along ket axis)
        # In Kron convention: A(x)B has [A(x)B][i*p+k, j*q+l] = A[i,j]*B[k,l]
        # We want: D_m[i'*nk + j', i*nk + j] = Z_bb[i', i] * Z_kk[j, j']
        # = Z_bb[i', i] * (Z_kk^T)[j', j]
        # So D_m = Z_bb (x) Z_kk^T (using standard kron).
        D_m = np.kron(Z_bra_bra[m], Z_ket_ket[m].T)
        L_block += GAMMA0 * D_m
        # Subtract gamma * delta_ii' delta_jj' (the -sigma part of D[sigma])
        L_block -= GAMMA0 * np.eye(n_ops, dtype=complex)

    # Step 3: diagonalize L_block
    print(f"\n--- step 3: diagonalize L_block ---")
    sys.stdout.flush()
    ev, V_block = np.linalg.eig(L_block)
    print(f"  L_block eigenvalue range: Re in [{ev.real.min():+.4f}, {ev.real.max():+.4f}], "
          f"|Im| in [{abs(ev.imag).min():.4e}, {abs(ev.imag).max():+.4e}]")

    # Step 4: find eigenvalues at |Im| ~ 0.0234
    target_Im = 0.02338
    Im_tol = 1e-4
    cluster_mask = (np.abs(np.abs(ev.imag) - target_Im) < Im_tol)
    cluster_idx = np.where(cluster_mask)[0]
    print(f"\n--- step 4: eigenvalues at |Im|={target_Im} +/- {Im_tol} ---")
    print(f"  found {len(cluster_idx)} cluster modes (expecting 4)")
    if len(cluster_idx) > 0:
        for i in cluster_idx:
            print(f"    Re={ev[i].real:+.4f}, Im={ev[i].imag:+.5e}")

    if len(cluster_idx) != 4:
        # broaden tolerance
        cluster_mask = (np.abs(np.abs(ev.imag) - target_Im) < 1e-3)
        cluster_idx = np.where(cluster_mask)[0]
        print(f"  broadened to tol 1e-3: {len(cluster_idx)} modes")
        for i in cluster_idx:
            print(f"    Re={ev[i].real:+.4f}, Im={ev[i].imag:+.5e}")

    # Step 5: look at cluster eigenvectors in the (i, j) index space
    print(f"\n--- step 5: structure of cluster eigenvectors ---")
    if len(cluster_idx) >= 2:
        # Sort by Im to identify R=+1 vs R=-1 partner pairs
        cluster_ev = ev[cluster_idx]
        cluster_V = V_block[:, cluster_idx]

        # Compute R action on cluster eigenvectors (R acts on Hilbert space, but
        # in operator basis sigma_{ij}, R sigma R^-1 changes (i, j) by R's action
        # on |phi_i^b>, |phi_j^k>. We need R_bb and R_kk matrices.
        R_bb = phi_b.conj().T @ R_op @ phi_b
        R_kk = phi_k.conj().T @ R_op @ phi_k

        # R super-op on sigma_{ij}: R sigma_{ij} R^-1 -> sum_{i'j'} R_bb[i',i] R_kk[j',j] sigma_{i' j'}
        # (since R is real and orthogonal, R^-1 = R^T = R^* if real)
        # = (R_bb (x) R_kk) sigma_vec  in convention idx = i*nk + j
        # but wait, R should act symmetrically: R rho R^dagger. R is real orthogonal, so R^dagger = R^T.
        R_super_block = np.kron(R_bb, R_kk.conj())

        for idx, ce in enumerate(cluster_ev):
            v = cluster_V[:, idx]
            # Apply R super-op: should give +/- v (R-eigenstate)
            Rv = R_super_block @ v
            # Compute overlap <v | Rv>
            ov = np.vdot(v, Rv) / np.linalg.norm(v) ** 2
            print(f"  cluster {idx}: Re={ce.real:+.4f} Im={ce.imag:+.5e}, "
                  f"<v|R|v>={ov.real:+.3f}{ov.imag:+.3f}j")

        # Identify R=+1 and R=-1 cluster modes; within each R-sector check if 2 modes share L-eigval
        print(f"\n  R-classified cluster modes:")
        R_signs = []
        for idx in range(len(cluster_idx)):
            v = cluster_V[:, idx]
            Rv = R_super_block @ v
            ov = float(np.vdot(v, Rv).real / np.linalg.norm(v) ** 2)
            R_signs.append(round(ov, 2))
        for r_val in [1.0, -1.0]:
            same_R = [i for i in range(len(cluster_idx)) if R_signs[i] == r_val]
            print(f"    R={r_val:+.1f}: {len(same_R)} modes, "
                  f"Im values {[float(ev[cluster_idx[i]].imag) for i in same_R]}")

        # Step 6: for the doublet (same R, same L-eigenvalue), examine eigenvector
        # structure. What pairs them?
        print(f"\n--- step 6: doublet eigenvector structure ---")
        # Group cluster modes by (R-sign, Im-sign)
        grouped = {}
        for idx in range(len(cluster_idx)):
            v = cluster_V[:, idx]
            Rv = R_super_block @ v
            ov = float(np.vdot(v, Rv).real / np.linalg.norm(v) ** 2)
            r_sign = +1 if ov > 0 else -1
            im_sign = +1 if ev[cluster_idx[idx]].imag > 0 else -1
            key = (r_sign, im_sign)
            grouped.setdefault(key, []).append(idx)

        for key, idx_list in grouped.items():
            r_sign, im_sign = key
            print(f"  (R={r_sign}, Im_sign={im_sign}): {len(idx_list)} modes")
            if len(idx_list) >= 2:
                # Take the two modes in the doublet
                v1 = cluster_V[:, idx_list[0]]
                v2 = cluster_V[:, idx_list[1]]
                # Reshape from 300 to (nb=20, nk=15)
                m1 = v1.reshape(nb, nk)
                m2 = v2.reshape(nb, nk)
                # Show top |coefficient| entries
                norms_m1 = np.abs(m1).flatten() ** 2
                top_idx = np.argsort(-norms_m1)[:6]
                print(f"    doublet member 1: top (bra_i, ket_j) by coefficient amplitude^2:")
                for ti in top_idx:
                    i_b, j_k = divmod(ti, nk)
                    val = m1[i_b, j_k]
                    print(f"      (bra {i_b} E={eigs_b[i_b]:+.4f}, "
                          f"ket {j_k} E={eigs_k[j_k]:+.4f}): "
                          f"|c|^2 = {abs(val)**2:.4f}")
                norms_m2 = np.abs(m2).flatten() ** 2
                top_idx2 = np.argsort(-norms_m2)[:6]
                print(f"    doublet member 2: top (bra_i, ket_j) by coefficient amplitude^2:")
                for ti in top_idx2:
                    i_b, j_k = divmod(ti, nk)
                    val = m2[i_b, j_k]
                    print(f"      (bra {i_b} E={eigs_b[i_b]:+.4f}, "
                          f"ket {j_k} E={eigs_k[j_k]:+.4f}): "
                          f"|c|^2 = {abs(val)**2:.4f}")

                # Overlap?
                overlap = np.vdot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                print(f"    <member 1 | member 2> = {overlap:+.4f}")


if __name__ == "__main__":
    main()
