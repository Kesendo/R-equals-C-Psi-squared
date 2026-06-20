"""Block-decomposed d_real(k): reach N=8,9,(10) for the d_real(2) closed-form hunt without the
full 4^N Liouvillian. The Heisenberg Liouvillian under Z-dephasing is EXACTLY block-diagonal in the
(popcount_a, popcount_b) coherence sectors (H conserves magnetization; dephasing is diagonal). So
d_real(k) = sum over (p,q) of the purely-real, on-grid eigenvalues of each small block L_{pq}.

L_{pq} = -i * ad_H|_{pq} - 2*gamma * diag(hamming(a,b)),  ad_H = kron(H_p, I) - kron(I, H_q),
H_p = Heisenberg H (XX+YY swap + ZZ diagonal) restricted to the popcount-p Hilbert sector.

GATE (gate-first): the block sum must reproduce the full-CSV d_real(2) for every (topology, N<=7)
that exists on disk. A mismatch means the block model is wrong -> diagnose, do not loosen. Only after
the gate passes are the N=8,9 numbers trustworthy.

Run:  python simulations/_dreal_block_decomposed.py
"""
from pathlib import Path
from collections import defaultdict
from math import comb
import numpy as np

RESULTS = Path(__file__).parent / "results"
GAMMA = 0.05          # gamma_phys; grid quantum gamma_q = 2*gamma = 0.1
GQ = 2 * GAMMA
TOL = 1e-7


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
    """Heisenberg H = J*sum_bonds (X_iX_j + Y_iY_j + Z_iZ_j) on the popcount-p sector."""
    states = [x for x in range(1 << N) if bin(x).count("1") == p]
    idx = {s: i for i, s in enumerate(states)}
    m = len(states)
    H = np.zeros((m, m))
    for ai, a in enumerate(states):
        d = 0.0
        for (i, j) in bonds:
            bi, bj = (a >> i) & 1, (a >> j) & 1
            d += J if bi == bj else -J          # Z_iZ_j
            if bi != bj:                          # X_iX_j + Y_iY_j swap, amplitude 2J
                a2 = a ^ (1 << i) ^ (1 << j)
                H[idx[a2], ai] += 2.0 * J
        H[ai, ai] += d
    return H, states


def block_dreal(N, topo):
    """d_real(k) for k=0..N via the (p,q)-block decomposition."""
    bonds = bonds_of(topo, N)
    Hs = {p: sector_H(N, p, bonds) for p in range(N + 1)}
    counts = defaultdict(int)
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
            for k in range(N + 1):
                counts[k] += int(np.sum(np.abs(real - (-k * GQ)) < TOL))
    return [counts[k] for k in range(N + 1)]


def csv_dreal2(topo, N):
    name = f"rmt_eigenvalues_N{N}.csv" if topo == "chain" else f"rmt_eigenvalues_{topo}_N{N}.csv"
    p = RESULTS / name
    if not p.exists():
        return None
    data = np.loadtxt(p, delimiter="\t", skiprows=1)
    ev = data[:, 0] + 1j * data[:, 1]
    real = ev[np.abs(ev.imag) < TOL].real
    return int(np.sum(np.abs(real - (-2 * GQ)) < TOL))


TOPOS = ("chain", "ring", "star", "complete")

# ---- GATE: block d_real(2) must equal full-CSV d_real(2) for all available N<=7 ----
print("=" * 84)
print("GATE: block-decomposed d_real(2) vs full-CSV d_real(2) (must match where the CSV exists)")
print("=" * 84)
gate_ok = True
for topo in TOPOS:
    for N in range(3, 8):
        csv = csv_dreal2(topo, N)
        if csv is None:
            continue
        blk = block_dreal(N, topo)[2]
        ok = blk == csv
        gate_ok &= ok
        print(f"  {topo:9} N={N}: block={blk:4d}  csv={csv:4d}  {'OK' if ok else 'MISMATCH'}")
        assert ok, f"GATE FIRED: {topo} N={N} block {blk} != csv {csv} (block model wrong)"
print(f"\nGATE {'PASS' if gate_ok else 'FAIL'}: the block decomposition reproduces the full Liouvillian.\n")

# ---- reach beyond the CSVs: N=8 (even) and N=9 (odd) ----
print("=" * 84)
print("EXTENSION: d_real(2) at N=8 (even, the missing even point) and N=9 (odd)")
print("=" * 84)
import sys
for N in (7, 8):
    line = f"  N={N}: "
    for topo in TOPOS:
        # N=7 cross-checks against the running full-dense export; N=8 is new (the even point)
        d2 = block_dreal(N, topo)[2]
        line += f"{topo}={d2}  "
        print(line, flush=True) if False else None
    print(line, flush=True)
print("\n(Feeds _dreal2_closed_form_hunt.py: even N now has 4,6,8(,10) and odd has 3,5,7,9.)")
