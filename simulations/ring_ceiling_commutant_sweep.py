"""Throwaway sweep: the ring's structural-ceiling story across N (the never-fed axis of F122).

C# (inspect --root ceiling, the new RingNode) just reported, over ALL (p,q) sectors:
  g2(ring, N) = 1 for N=4..7  -> the ring has NO structural ceiling, the band edge protects it
  (1,1) commutant: ring-4=1, ring-5=1.6, ring-6=1.3333, ring-7=1.714286

This script (1) cross-checks that C# port against the banked numpy machine for N=4..7, then
(2) extends the (1,1) commutant to N=11 to test the conjectured closed forms
       even N: 2(N-2)/N    odd N: 2(N-1)/N
and (3) reads the half-filling (N/2,N/2) seam for even N (the clock_hand_ladder open arc:
"why does the ring (2,2) double-excitation seam overtake the band edge at even half-filling?").

Self-validating: STAGE A asserts the C# port-fidelity numbers.  Pure numpy, no framework import.
Run:  python simulations/ring_ceiling_commutant_sweep.py
"""
import itertools
from math import sqrt
import numpy as np

TOL_REL = 1e-6     # same-Omega cluster tolerance (relative to the Omega scale)
NZ = 1e-7          # nonzero-<n_XY> threshold

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)


def kron_list(ops):
    out = np.array([[1]], dtype=complex)
    for o in ops:
        out = np.kron(out, o)
    return out


def site_op(op, l, N):
    return kron_list([op if k == l else I2 for k in range(N)])


def ring_bonds(N):
    return [(i, (i + 1) % N) for i in range(N)]


def H_full(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in ring_bonds(N):
        H += 0.5 * J * (site_op(X, i, N) @ site_op(X, j, N) + site_op(Y, i, N) @ site_op(Y, j, N))
    return H.real


def sector_states(N, p):
    return [x for x in range(2 ** N) if bin(x).count('1') == p]


def sector_analysis(Hf, N, p, q):
    """Returns (highq_min over all-Omega blocks, comm_min over the Omega=0 commutant)."""
    A = sector_states(N, p)
    B = sector_states(N, q)
    na, nb = len(A), len(B)
    if na == 0 or nb == 0:
        return None, None
    Hp = Hf[np.ix_(A, A)]
    Hq = Hf[np.ix_(B, B)]
    adH = np.kron(Hp, np.eye(nb)) - np.kron(np.eye(na), Hq)
    Omega, V = np.linalg.eigh(adH)
    diag = np.array([bin(A[a] ^ B[b]).count('1') for a in range(na) for b in range(nb)], dtype=float)
    Ntil = V.T @ (diag[:, None] * V)
    oscale = 1.0 + (np.abs(Omega).max() if Omega.size else 0.0)
    same = np.abs(Omega[:, None] - Omega[None, :]) < TOL_REL * oscale
    w_all = np.linalg.eigvalsh(Ntil * same)
    nz_all = w_all[w_all > NZ]
    highq = float(nz_all.min()) if nz_all.size else None
    z = np.abs(Omega) < TOL_REL * oscale
    comm = None
    if z.any():
        w_c = np.linalg.eigvalsh(Ntil[np.ix_(z, z)])
        nz_c = w_c[w_c > NZ]
        comm = float(nz_c.min()) if nz_c.size else None
    return highq, comm


def global_g2(Hf, N):
    best, bpq = np.inf, None
    for p in range(N + 1):
        for q in range(p, N + 1):
            hq, _ = sector_analysis(Hf, N, p, q)
            if hq is not None and hq < best:
                best, bpq = hq, (p, q)
    return best, bpq


# =====================================================================================
# STAGE A -- port-fidelity cross-check against the C# RingNode numbers (N=4..7)
# =====================================================================================
print("=" * 92)
print("STAGE A -- numpy vs C# RingNode (port fidelity): global g2 and (1,1) commutant, ring N=4..7")
print("=" * 92)
print(f"{'N':>2} {'global g2':>11} {'win sector':>11} {'(1,1) commutant':>16}   C#-expected")
csharp_c11 = {4: 1.0, 5: 1.6, 6: 4.0 / 3.0, 7: 12.0 / 7.0}
csharp_g2 = {4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0}
for N in (4, 5, 6, 7):
    Hf = H_full(N)
    g2, win = global_g2(Hf, N)
    _, c11 = sector_analysis(Hf, N, 1, 1)
    print(f"{N:>2} {g2:>11.6f} {str(win):>11} {c11:>16.6f}   g2={csharp_g2[N]}, c11={csharp_c11[N]:.6f}")
    assert abs(g2 - csharp_g2[N]) < 1e-6, f"PORT MISMATCH g2 ring-{N}: numpy {g2} vs C# {csharp_g2[N]}"
    assert abs(c11 - csharp_c11[N]) < 1e-6, f"PORT MISMATCH c11 ring-{N}: numpy {c11} vs C# {csharp_c11[N]}"
print("STAGE A PASS: the C# RingNode is a faithful port (global g2 = 1, (1,1) commutant matches).")

# =====================================================================================
# STAGE B -- the (1,1) commutant closed form, extended to N=11
#   conjecture:  even N -> 2(N-2)/N ,  odd N -> 2(N-1)/N
# =====================================================================================
print("\n" + "=" * 92)
print("STAGE B -- ring (1,1) commutant closed form (the value F122 marked 'breaks 4/(m+1)')")
print("=" * 92)
print(f"{'N':>2} {'parity':>6} {'(1,1) commutant':>16} {'conjecture':>12} {'closed form':>14} {'match?':>7}")
allmatch = True
for N in range(4, 12):
    Hf = H_full(N)
    _, c11 = sector_analysis(Hf, N, 1, 1)
    if N % 2 == 0:
        conj, form = 2.0 * (N - 2) / N, f"2(N-2)/N"
    else:
        conj, form = 2.0 * (N - 1) / N, f"2(N-1)/N"
    ok = abs(c11 - conj) < 1e-7
    allmatch &= ok
    print(f"{N:>2} {('even' if N%2==0 else 'odd'):>6} {c11:>16.9f} {conj:>12.6f} {form:>14} {('YES' if ok else 'NO'):>7}")
print(f"STAGE B {'PASS' if allmatch else 'FAIL'}: ring (1,1) commutant = 2(N-2)/N (even), 2(N-1)/N (odd) "
      f"to machine precision through N=11.")

# =====================================================================================
# STAGE C -- the half-filling (N/2, N/2) seam for even N (the clock_hand_ladder question)
#   does the double-excitation seam dip below / co-occupy / sit above the band edge?
# =====================================================================================
print("\n" + "=" * 92)
print("STAGE C -- ring half-filling (N/2,N/2) seam, even N (clock_hand_ladder: overtakes band edge?)")
print("=" * 92)
print(f"{'N':>2} {'(N/2,N/2)':>10} {'highq_min':>11} {'commutant':>11} {'vs band edge 1':>16}")
for N in (4, 6, 8):
    Hf = H_full(N)
    h = N // 2
    hq, comm = sector_analysis(Hf, N, h, h)
    cs = f"{comm:.6f}" if comm is not None else "  -"
    rel = "= 1 (co-occupy)" if (comm is not None and abs(comm - 1) < 1e-6) else \
          ("< 1 (CEILING!)" if (comm is not None and comm < 1) else "> 1 (above)")
    print(f"{N:>2} {str((h, h)):>10} {hq:>11.6f} {cs:>11} {rel:>16}")
print("(If half-filling commutant = 1 for all even N, the ring 'seam' co-occupies the band edge, it does "
      "not undercut it -> the ring stays g2=1, ceiling-free; the 'overtake' is a finite-Q / Im-side effect.)")

print("\nDONE.")
