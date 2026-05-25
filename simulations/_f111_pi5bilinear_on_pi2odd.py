"""
F111 exploration: does Pi_5bilinear (Z-deph variant) COMMUTE with Pi^2-Z-odd
2-bilinear commutator superops? F108 says anti-commute=0 for Pi^2-Z-EVEN.
Maybe it COMMUTES for Pi^2-Z-ODD.

For 2-qubit XY bilinear: compute Q [H, .] Q^-1 vs +[H, .] and -[H, .].

This would mean: for Pi^2-D-odd Hamiltonians, conjugation by Pi_5bilinear
gives Q L_H Q^-1 = +L_H. Then combined with Y^N (which gives -L_H),
the product (Pi_5bilinear) (Y^N-superop) achieves a constraint.

Also check: for the canonical Pi (NOT 5bilinear), what's the action on Pi^2-D-odd bilinears?
"""

import os
from itertools import product as iprod

import numpy as np


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": sx, "Y": sy, "Z": sz}
LABELS = ["I", "X", "Y", "Z"]


def build_pauli_op_full(letters):
    mat = PAULI[letters[0]]
    for k in range(1, len(letters)):
        mat = np.kron(mat, PAULI[letters[k]])
    return mat


def build_Q_Nsite(per_site, N):
    d = 2 ** N
    d2 = d * d
    basis_N = list(iprod(LABELS, repeat=N))
    label_to_idx = {bl: i for i, bl in enumerate(basis_N)}
    mats_N = []
    for bl in basis_N:
        mat = PAULI[bl[0]]
        for k in range(1, N):
            mat = np.kron(mat, PAULI[bl[k]])
        mats_N.append(mat)
    vecs_N = [m.flatten() for m in mats_N]
    Q = np.zeros((d2, d2), dtype=complex)
    for idx, bl in enumerate(basis_N):
        phase = 1.0
        tgt_labels = []
        for sl in bl:
            ph, t = per_site[sl]
            phase *= ph
            tgt_labels.append(t)
        tgt_idx = label_to_idx[tuple(tgt_labels)]
        Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
    return Q


def pi5bilinear_per_site(dephase):
    if dephase == "Z":
        return {"I": (1, "X"), "X": (-1, "I"), "Y": (1j, "Z"), "Z": (-1j, "Y")}
    elif dephase == "X":
        return {"I": (1, "Z"), "Z": (-1, "I"), "X": (-1j, "Y"), "Y": (1j, "X")}
    elif dephase == "Y":
        return {"I": (1, "X"), "X": (-1, "I"), "Y": (-1j, "Z"), "Z": (1j, "Y")}
    raise ValueError(dephase)


def canonical_pi_per_site(dephase):
    if dephase == "Z":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}
    elif dephase == "X":
        return {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}
    elif dephase == "Y":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}
    raise ValueError(dephase)


def commutator_superop(H):
    d = H.shape[0]
    return np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T)


def main():
    print("=" * 80)
    print("F111: Pi_5bilinear action on Pi^2-D-odd bilinears")
    print("=" * 80)

    bilinears = {
        "Z_even": ["XX", "YY", "YZ", "ZY", "ZZ"],
        "Z_odd": ["XY", "XZ", "YX", "ZX"],
    }

    Q_5bi = build_Q_Nsite(pi5bilinear_per_site("Z"), 2)
    Q_can = build_Q_Nsite(canonical_pi_per_site("Z"), 2)
    Q_5bi_inv = np.linalg.inv(Q_5bi)
    Q_can_inv = np.linalg.inv(Q_can)

    print(f"\n{'Bilinear':<10} {'||Q5bi*C-(-C)*Q5bi||':>22} {'||Q5bi*C-(+C)*Q5bi||':>22} {'||Qcan*C-(-C)*Qcan||':>22}")
    for parity, bils in bilinears.items():
        print(f"\n{parity} bilinears:")
        for bil in bils:
            H = build_pauli_op_full(bil)
            C = commutator_superop(H)
            # Q anti-commutes with C: ||{Q, C}|| = ||QC + CQ||
            anti_5bi = np.linalg.norm(Q_5bi @ C + C @ Q_5bi)
            comm_5bi = np.linalg.norm(Q_5bi @ C - C @ Q_5bi)
            anti_can = np.linalg.norm(Q_can @ C + C @ Q_can)
            print(f"  {bil:<8} anti_5bi={anti_5bi:>10.3e}, "
                  f"comm_5bi={comm_5bi:>10.3e}, "
                  f"anti_can={anti_can:>10.3e}")


    print()
    print("=" * 80)
    print("Y^N-conjugation on Pi^2-D-odd: H -> -H means commutator C -> -C")
    print("=" * 80)
    print("Y^N as 4x4 = Y o Y; Y^N H Y^N gives sign flip per non-I letter X,Z")
    print()
    Y_2 = np.kron(sy, sy)
    Y2_inv = np.linalg.inv(Y_2)
    print(f"{'Bilinear':<10} {'Y^2 H Y^-2 ?= H':>20}")
    for parity, bils in bilinears.items():
        print(f"\n{parity} bilinears:")
        for bil in bils:
            H = build_pauli_op_full(bil)
            conj = Y_2 @ H @ Y2_inv
            ratio = "+H" if np.allclose(conj, H) else "-H" if np.allclose(conj, -H) else "other"
            n_x = bil.count("X")
            n_z = bil.count("Z")
            n_y = bil.count("Y")
            sign_expect = (-1)**(n_x + n_z)
            sign_str = "+1" if sign_expect > 0 else "-1"
            print(f"  {bil:<8} #X={n_x}, #Y={n_y}, #Z={n_z}, sign expected={sign_str}, observed={ratio}")


if __name__ == "__main__":
    main()
