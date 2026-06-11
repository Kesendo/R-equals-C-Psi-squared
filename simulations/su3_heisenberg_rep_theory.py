#!/usr/bin/env python3
"""
_su3_heisenberg_rep_theory.py - WIP scout: derive the SU(3)-Heisenberg
interacting count (60) and the <Q> distribution {0:6, 1:36, 1.5:12, 2:27}
from representation theory.

Setup (N=2 qutrits, full-Cartan dephasing, H = sum_a lambda_a x lambda_a):
the full L = L_H + L_D has real parts -2g*<Q> exactly (the symmetric case of
F121 sec.4); <Q> in {0,1,1.5,2}. WHY those multiplicities?

The 2-qutrit Hilbert space C3 x C3 = 6 (+) 3bar (sym (+) antisym); H is constant
on each (a Casimir). The Liouville space End(C3 x C3) splits into four L_H-blocks
by energy difference:
  A = 6 x 6bar   L_H=0     dim 36   = 1 (+) 8 (+) 27
  B = 3bar x 3   L_H=0     dim  9   = 1 (+) 8
  C = 6 x 3bar   L_H=-iD   dim 18   = 8 (+) 10
  D = 3bar x 6bar L_H=+iD  dim 18   = 8 (+) 10bar
The global SU(3) Casimir on Liouville space identifies the irreps; the dephasing
Q breaks SU(3) but commutes with the COMMON symmetry of H and L_D (the global
level-permutations S3 and the site swap). This scout measures how <Q> distributes
over the SU(3) blocks/irreps and whether it is a good quantum number on the
common-symmetry sectors.

Blocks (this scout pins, exact):
  [1] SU(3) irrep multiplicities of the 81-dim Liouville space {1:2, 8:4, 27:1, 10:1, 10bar:1}
  [2] L_H spectrum: 0 (x45), +/- i*Delta (x18 each)
  [3] <Q> per L_H-block and per SU(3) irrep -> the {0:6,1:36,1.5:12,2:27} decomposition
  [4] the 60 pairing read off the blocks
"""

import numpy as np
from itertools import product as iprod

# ---- Gell-Mann + qutrit infrastructure (matches the committed verifiers) ----
gm = [
    np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
    np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3),
]
I3 = np.eye(3, dtype=complex)
I9 = np.eye(9, dtype=complex)
GAMMA = 0.05


def kron(a, b):
    return np.kron(a, b)


def main():
    np.set_printoptions(suppress=True, precision=4)
    d, N = 3, 2
    g = GAMMA

    # ---- Hamiltonian and Liouvillian pieces ----
    H = sum(kron(m, m) for m in gm)                 # 9x9 SU(3) Heisenberg, one bond
    LH = -1j * (kron(H, I9) - kron(I9, H.T))        # 81x81 commutator superoperator
    # full-Cartan dephasing dissipator with jumps lambda_3, lambda_8 on each site
    LD = np.zeros((81, 81), dtype=complex)
    for site in (0, 1):
        for m in (gm[2], gm[7]):
            Mk = kron(m, I3) if site == 0 else kron(I3, m)
            MdM = Mk.conj().T @ Mk
            LD += g * (kron(Mk, Mk.conj()) - 0.5 * kron(MdM, I9) - 0.5 * kron(I9, MdM.T))
    L = LH + LD

    # Q operator (Hamming distance) on the 81 coherence indices, diagonal
    states = list(iprod(range(d), repeat=N))
    Qdiag = np.array([sum(1 for a, b in zip(i, j) if a != b)
                      for i in states for j in states], dtype=float)

    # ---- [1] SU(3) Casimir on Liouville space -> irrep multiplicities ----
    G = [kron(m, I3) + kron(I3, m) for m in gm]      # global SU(3) generators, 9x9
    Cas = np.zeros((81, 81), dtype=complex)
    for Ga in G:
        ad = kron(Ga, I9) - kron(I9, Ga.T)           # adjoint generator on Liouville space
        Cas += ad @ ad
    cas_eigs = np.linalg.eigvalsh(Cas)               # Hermitian
    print("=" * 72)
    print("SU(3)-Heisenberg rep theory of the 81-dim Liouville space (N=2 qutrits)")
    print("=" * 72)
    print("\n[1] SU(3) Casimir spectrum (distinct value : multiplicity):")
    vals, counts = np.unique(np.round(cas_eigs.real, 3), return_counts=True)
    # identify irreps by dimension: 1,8,27,10/10bar (10 & 10bar share a Casimir)
    irrep_by_cas = {}
    for v, cnt in zip(vals, counts):
        print(f"    C2 = {v:8.3f}  total dim {cnt}")
        irrep_by_cas[v] = cnt
    print("    expected: dims {2 (=1x2), 32 (=8x4), 20 (=10+10bar), 27}")

    # ---- [2] L_H spectrum (energy differences) ----
    lh_eigs = np.linalg.eigvals(LH)
    print("\n[2] L_H spectrum (imaginary energy differences):")
    lh_round = np.round(lh_eigs.imag / g, 2)
    lvals, lcounts = np.unique(lh_round, return_counts=True)
    for v, cnt in zip(lvals, lcounts):
        print(f"    Im(L_H)/g = {v:8.2f}  mult {cnt}")
    # H eigenvalues on sym(6) vs antisym(3bar)
    hev = np.linalg.eigvalsh(H)
    hv, hc = np.unique(np.round(hev, 3), return_counts=True)
    print(f"    H eigenvalues (Hilbert): {dict(zip(hv, hc))}  (6 = sym, 3 = antisym)")

    # ---- [3] sym/antisym block structure and <Q> per block ----
    # symmetric / antisymmetric projectors on C3 x C3 (site swap)
    SWAP = np.zeros((9, 9), dtype=complex)
    for a, (a0, a1) in enumerate(states):
        b = states.index((a1, a0))
        SWAP[b, a] = 1.0
    Psym = (I9 + SWAP) / 2     # rank 6
    Panti = (I9 - SWAP) / 2    # rank 3
    print(f"\n[3] sym/antisym dims: 6={int(round(np.trace(Psym).real))}, "
          f"3bar={int(round(np.trace(Panti).real))}")

    # Dephasing BREAKS SU(3) (single-site jumps do not commute with the global
    # SU(3) nor the site swap), so <Q> is NOT an SU(3) class function: compressing
    # L_D into an irrep gives a spread. The clean split is by H's ENERGY structure
    # (which DOES commute with L_H): intra-sector (L_H = 0, the 6x6bar (+) 3barx3
    # coherences) vs inter-sector (L_H = +/- i*Delta, the 6<->3bar coherences).
    # Use the FULL-L eigenmodes directly, classified by |Im(lambda)|.
    P_inter = kron(Psym, Panti.T) + kron(Panti, Psym.T)   # projector onto C (+) D
    ev, Vr = np.linalg.eig(L)
    W = np.linalg.inv(Vr)
    qfull = np.einsum('ki,i,ik->k', W, Qdiag, Vr).real
    inter_wt = np.real(np.einsum('ki,ij,jk->k', W, P_inter, Vr))   # biorthogonal inter weight
    print("\n    full-L modes split by H-energy sector (Im(lambda) large = inter 6<->3bar):")
    intra = np.abs(ev.imag) < 1e-6
    inter = ~intra

    def qhist(mask):
        qv, qc = np.unique(np.round(qfull[mask], 3), return_counts=True)
        return {float(a): int(b) for a, b in zip(qv, qc)}

    qi, qo = qhist(intra), qhist(inter)
    print(f"      intra-sector (L_H=0,   {int(intra.sum()):2d} modes): <Q> = {qi}")
    print(f"      inter-sector (L_H=+-iD, {int(inter.sum()):2d} modes): <Q> = {qo}")
    # the half-integer <Q>=1.5 lives ONLY in the inter-sector
    half = np.abs(qfull - 1.5) < 1e-6
    print(f"\n    the 12 <Q>=1.5 modes: all inter-sector? "
          f"{bool(np.all(inter[half]))}; inter-weight ~ 1? "
          f"min={inter_wt[half].min():.3f}")
    print(f"    intra-sector <Q> values are integers only: "
          f"{sorted(set(np.round(qfull[intra]).astype(int)))}")
    assert np.all(inter[half]), "a <Q>=1.5 mode is NOT inter-sector"
    qhist_total = qhist(intra | inter)
    assert {round(k, 1): v for k, v in qhist_total.items()} == {0.0: 6, 1.0: 36, 1.5: 12, 2.0: 27}
    print(f"\n    reconstructed full-L <Q>: {qhist_total}  (target {{0:6,1:36,1.5:12,2:27}})")

    # ---- [4] the 60: palindrome pairs (<Q>, Im) -> (3 - <Q>, -Im) about center -3g ----
    # the reflection lambda -> -2*Sg - lambda at center -3g sends <Q> -> 3 - <Q>
    # AND Im -> -Im, so a mode (q, omega) pairs with (3 - q, -omega).
    from collections import Counter
    cells = Counter((round(q, 1), int(np.sign(round(w / g)))) for q, w in zip(qfull, ev.imag))
    print("\n[4] census by (<Q>, sign Im) cell:")
    for (q, s), n in sorted(cells.items()):
        tag = "intra" if s == 0 else ("inter +iD" if s > 0 else "inter -iD")
        print(f"    <Q>={q:<4} {tag:<10}: {n}")
    # reflect (q, s) -> (3 - q, -s); the <Q>=1.5 cells map (+iD)<->(-iD) into each
    # other (never onto themselves, since the sign flips), so the plain loop pairs
    # everything once with no self-mirror special case.
    paired, used = 0, set()
    for (q, s), n in cells.items():
        if (q, s) in used:
            continue
        pq, ps = round(3.0 - q, 1), -s
        m = min(n, cells.get((pq, ps), 0))
        paired += 2 * m
        used.add((q, s)); used.add((pq, ps))
    print(f"\n    palindrome pairing from the (<Q>, Im) census = {paired}/81")
    print("    decomposition: intra Q1<->Q2 = 2*min(18,21) = 36;  inter Q1<->Q2 across +-iD")
    print("    = 2*(min(9,3)+min(9,3)) = 12;  inter Q=1.5 self across +-iD = 2*min(6,6) = 12;")
    print("    total 36 + 12 + 12 = 60.")
    assert paired == 60, f"rep-theory pairing census gives {paired}, expected 60"

    print("\n" + "=" * 72)
    print("READING: SU(3) gives the operator-space SKELETON ({1:2,8:4,27:1,10:1,10bar:1})")
    print("and, through H's two energy levels (6 at +4/3, 3bar at -8/3), the L_H blocks")
    print("(intra L_H=0 x45, inter 6<->3bar L_H=+-4J x36). But the dephasing BREAKS SU(3),")
    print("so <Q> is a common-symmetry count, not an SU(3) class function. The clean")
    print("rep-theory reading: the -3g rung (<Q>=1.5) is EXACTLY the inter-sector 6<->3bar")
    print("coherences - the half-integer Hamming average is the signature of crossing the")
    print("symmetric/antisymmetric divide; intra-sector modes carry integer <Q> only.")
    print("=" * 72)


if __name__ == "__main__":
    main()
