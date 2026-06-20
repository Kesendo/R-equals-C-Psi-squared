"""N=9 d_real(2) via block decomposition: the 4th odd point (N=3,5,7,9).

Verifies the sparse linear forms (chain 4N-6 -> 30, ring 7N-13 -> 50) and gives star/complete a 4th
odd point. Block functions are a copy of _dreal_block_decomposed.py (already gate-validated against
the full CSVs through N=7). N=9's biggest blocks are C(9,4)=126 -> 15876^2 (feasible, ~minutes each).

Run:  python simulations/_dreal_block_n9.py
"""
from collections import defaultdict
import numpy as np

GAMMA = 0.05
GQ = 2 * GAMMA
TOL = 1e-7
TOPOS = ("chain", "ring", "star", "complete")


def bonds_of(topo, N):
    if topo == "chain":
        return [(i, i + 1) for i in range(N - 1)]
    if topo == "ring":
        return [(i, (i + 1) % N) for i in range(N)]
    if topo == "star":
        return [(0, i) for i in range(1, N)]
    if topo == "complete":
        return [(i, j) for i in range(N) for j in range(i + 1, N)]
    raise ValueError(topo)


def sector_H(N, p, bonds, J=1.0):
    states = [x for x in range(1 << N) if bin(x).count("1") == p]
    idx = {s: i for i, s in enumerate(states)}
    m = len(states)
    H = np.zeros((m, m))
    for ai, a in enumerate(states):
        d = 0.0
        for (i, j) in bonds:
            bi, bj = (a >> i) & 1, (a >> j) & 1
            d += J if bi == bj else -J
            if bi != bj:
                a2 = a ^ (1 << i) ^ (1 << j)
                H[idx[a2], ai] += 2.0 * J
        H[ai, ai] += d
    return H, states


def block_dreal2(N, topo):
    bonds = bonds_of(topo, N)
    Hs = {p: sector_H(N, p, bonds) for p in range(N + 1)}
    count = 0
    for p in range(N + 1):
        Hp, A = Hs[p]
        na = len(A)
        for q in range(N + 1):
            Hq, B = Hs[q]
            nb = len(B)
            if na == 0 or nb == 0:
                continue
            adH = np.kron(Hp, np.eye(nb)) - np.kron(np.eye(na), Hq)
            ham = np.array([bin(A[a] ^ B[b]).count("1") for a in range(na) for b in range(nb)], float)
            L = -1j * adH - 2 * GAMMA * np.diag(ham)
            ev = np.linalg.eigvals(L)
            real = ev[np.abs(ev.imag) < TOL].real
            count += int(np.sum(np.abs(real - (-2 * GQ)) < TOL))
    return count


print("N=9 d_real(2) (the 4th odd point):", flush=True)
res = {}
for topo in TOPOS:
    v = block_dreal2(9, topo)
    res[topo] = v
    print(f"  {topo:9} = {v}", flush=True)

print("\nVERIFY sparse linear forms at N=9:")
print(f"  chain: {res['chain']} vs 4N-6 = {4*9-6}  -> {'OK' if res['chain']==4*9-6 else 'FAIL'}")
print(f"  ring:  {res['ring']} vs 7N-13 = {7*9-13}  -> {'OK' if res['ring']==7*9-13 else 'FAIL'}")
print(f"\nDense odd 4 points (N=3,5,7,9): star=[6,30,136,{res['star']}]  complete=[8,54,216,{res['complete']}]")
