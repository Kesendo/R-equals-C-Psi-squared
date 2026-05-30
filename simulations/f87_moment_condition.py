#!/usr/bin/env python3
"""Resolve the first-order soft condition (moment form).

Claim: soft <=> Tr((L+sigma)^(2k+1)) = 0 for all k, and to first order in gamma this is
   S_k := sum_ab (E_a-E_b)^(2k) <z_a,z_b> = 0  (k>=1),  z_a^l = <E_a|Z_l|E_a>.
A hand expansion suggested bipartite forces S_2 = 6*sum_l Tr(H^2 Z_l)^2 != 0, which would
contradict bipartite=>soft. So there is an error. Compute, for a SOFT and a HARD pair:
  - Tr(H^2 Z_l) per site
  - S_k (k=1,2,3) directly
  - the ACTUAL first-order odd moments Tr((L+sigma)^(2k+1))/gamma at small gamma
and see what is actually zero.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain, site_op


def analyze(pair, N=4, gamma=1e-4):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    Z = [site_op(N, l, 'Z') for l in range(N)]
    z = np.array([[np.real(V[:, a].conj() @ (Zl @ V[:, a])) for Zl in Z] for a in range(len(E))])
    zz = z @ z.T                                    # <z_a,z_b>
    omega = E[:, None] - E[None, :]
    # Tr(H^2 Z_l)
    H2 = H @ H
    trH2Z = [np.real(np.trace(H2 @ Zl)) for Zl in Z]
    # S_k = sum_ab omega^(2k) <z_a,z_b>
    Sk = [float(np.sum((omega ** (2 * k)) * zz)) for k in (1, 2, 3)]
    # actual first-order odd moments
    d = 2 ** N
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter='Z')
    sigma = N * gamma
    M = L + sigma * np.eye(d * d)
    moms = {}
    P = M.copy()
    for p in range(1, 6):
        moms[p] = np.trace(P) if p == 1 else moms.get(p)
        P = P  # placeholder
    # compute Tr(M^p) for p=1,3,5
    tr = {}
    Mp = np.eye(d * d)
    for p in range(1, 6):
        Mp = Mp @ M
        tr[p] = complex(np.trace(Mp))
    return trH2Z, Sk, tr, gamma


def main():
    for label, pair in [("SOFT  XXZ+ZXX", [('X', 'X', 'Z'), ('Z', 'X', 'X')]),
                        ("HARD  XXZ+XZX", [('X', 'X', 'Z'), ('X', 'Z', 'X')])]:
        trH2Z, Sk, tr, g = analyze(pair)
        print(f"{label}:")
        print(f"   Tr(H^2 Z_l) = {[round(x,4) for x in trH2Z]}")
        print(f"   S_1={Sk[0]:.4f}  S_2={Sk[1]:.4f}  S_3={Sk[2]:.4f}")
        print(f"   Tr((L+σ)^3)/γ = {tr[3].real/g:.4f}{tr[3].imag/g:+.4f}i   "
              f"(predict -3·S_1 = {-3*Sk[0]:.4f})")
        print(f"   Tr((L+σ)^5)/γ = {tr[5].real/g:.4f}{tr[5].imag/g:+.4f}i   "
              f"(predict 5·S_2 = {5*Sk[1]:.4f})")
        print()


if __name__ == "__main__":
    main()
