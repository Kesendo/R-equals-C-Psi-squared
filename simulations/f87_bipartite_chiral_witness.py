#!/usr/bin/env python3
"""Computational witness for the F87 diagonal-cell bipartite mechanism.

For a diagonal-cell pair, soft <=> H's hopping graph (in the dephasing letter's eigenbasis)
is bipartite, i.e. a chiral K (diagonal there, KHK=-H) exists. Proof of bipartite => soft,
verified link by link:

  link 1 (F80 one-sidedness):  M := Pi L Pi^-1 + L + 2sigma  =  -2i (H (x) I)
  link 2 (algebra, from 1):    Pi L Pi^-1  =  -i{H,.} - D - 2sigma
  link 3 (chiral mirror):      W=K(x)I with KHK=-H sends (-i{H,.}-D) -> (i[H,.]-D) = -L
  => Spec(L) = Spec(Pi L Pi^-1) = Spec(-i{H,.}-D-2sigma) = Spec(-L-2sigma) : SOFT.

For a hard pair, link 3's K does not exist (odd hopping cycle), and Spec(L) != Spec(-L-2sigma)
(optimal lambda<->-lambda-2sigma matching leaves a finite residual).
"""
from __future__ import annotations
import sys
from collections import deque
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.lindblad import lindbladian_pauli_dephasing, palindrome_residual
from framework.pauli import _build_kbody_chain, _vec_to_pauli_basis_transform


def two_coloring(H, tol=1e-9):
    """Return diag(+-1) K with KHK=-H if H's hopping graph is bipartite (zero diagonal),
    else None."""
    d = H.shape[0]
    if np.max(np.abs(np.diag(H))) > tol:
        return None
    adj = [np.where(np.abs(H[a]) > tol)[0] for a in range(d)]
    color = {}
    for s in range(d):
        if s in color:
            continue
        color[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v == u:
                    continue
                if v not in color:
                    color[v] = color[u] ^ 1
                    q.append(v)
                elif color[v] == color[u]:
                    return None
    return np.diag([1.0 if color[a] == 0 else -1.0 for a in range(d)]).astype(complex)


def optimal_pairing_residual(L, Sigma):
    ev = np.linalg.eigvals(L)
    tgt = -ev - 2 * Sigma
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def marginal_diffs(L, Sigma):
    """(Re-marginal symmetric about -sigma, Im-marginal symmetric about 0) diffs.

    Im always pairs (Lindbladian spectrum is conjugation-closed); for hard pairs the Re
    marginal (decay rates) no longer pairs about -sigma. This is the §7.3 observation."""
    ev = np.linalg.eigvals(L)
    re_diff = float(np.max(np.abs(np.sort(ev.real) - np.sort(-2 * Sigma - ev.real))))
    im_diff = float(np.max(np.abs(np.sort(ev.imag) - np.sort(-ev.imag))))
    return re_diff, im_diff


def witness(label, pair, N=4):
    chain = fw.ChainSystem(N=N)
    g0 = chain.gamma_0
    Sigma = N * g0
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    T = _vec_to_pauli_basis_transform(N)

    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    L = lindbladian_pauli_dephasing(H, [g0] * N, dephase_letter='Z')

    # link 1: M = -2i (H (x) I)
    M_pauli = palindrome_residual(L, Sigma, N, dephase_letter='Z')
    M_vec = (T @ M_pauli @ T.conj().T) / d
    link1 = np.allclose(M_vec, -2j * np.kron(H, Id), atol=1e-9)

    # dissipator D and the two H-superoperators
    comm = np.kron(H, Id) - np.kron(Id, H.T)
    anti = np.kron(H, Id) + np.kron(Id, H.T)
    D = L + 1j * comm                       # since L = -i*comm + D
    G_minus = -1j * anti - D                # -i{H,.} - D
    minus_L = 1j * comm - D                 # i[H,.] - D = -L

    K = two_coloring(H)
    print(f"{label}:")
    print(f"   link1  M = -2i(H(x)I)                 : {link1}")
    if K is None:
        print(f"   link2/3  chiral K (KHK=-H)            : DOES NOT EXIST (non-bipartite)")
    else:
        link2 = np.allclose(K @ H @ K + H, 0, atol=1e-9)   # KHK = -H
        W = np.kron(K, Id)                                  # W^-1 = W (K^2=I)
        link3 = np.allclose(W @ G_minus @ W, minus_L, atol=1e-9)
        print(f"   link2  KHK = -H                       : {link2}")
        print(f"   link3  (K(x)I)(-i{{H,.}}-D)(K(x)I) = -L : {link3}")
    re_diff, im_diff = marginal_diffs(L, Sigma)
    print(f"   marginals   Re about -sigma diff = {re_diff:.3e}   Im about 0 diff = {im_diff:.3e}")
    print(f"   conclusion  optimal |lambda<->-lambda-2s| residual = {optimal_pairing_residual(L, Sigma):.3e}")
    print()


def main():
    witness("SOFT  XXZ+ZXX (bipartite)", [('X', 'X', 'Z'), ('Z', 'X', 'X')])
    witness("HARD  XXZ+XZX (odd cycle)", [('X', 'X', 'Z'), ('X', 'Z', 'X')])
    witness("HARD  ZZZ+XXZ (pure-D template lifts diagonal)", [('Z', 'Z', 'Z'), ('X', 'X', 'Z')])


if __name__ == "__main__":
    main()
