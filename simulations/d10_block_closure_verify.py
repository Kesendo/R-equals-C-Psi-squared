#!/usr/bin/env python3
"""D10's (0,1) block, verified past the range a full 4^N eigendecomposition can reach.

D10 (docs/proofs/derivations/D10_W1_DISPERSION.md) Step 3 derives, for the Heisenberg chain
under uniform Z-dephasing,

    L|(0,1) = -2i*J*Laplacian - 2*gamma*Id

on the N-dimensional block spanned by the |0><j| between the ferromagnet and the single
excitations, and Step 4 reads off omega_k = 4J(1 - cos(pi*k/N)).

D10's Numerical Verification section compares against a FULL Liouvillian eigendecomposition, and
simulations/analytical_spectrum_verify.py runs that comparison to N=6. This script does not need
it. The block is claimed INVARIANT, so applying L to the N basis elements one at a time settles
everything at O(N * 4^N) memory instead of O(16^N):

  (1) does the image of each |0><j| stay inside the block (leak = 0),
  (2) is the block operator entry-wise the closed form, and
  (3) do its frequencies match omega_k.

If (1) holds, the block's eigenvalues ARE eigenvalues of L, so (3) is a statement about the full
spectrum, not about a truncation. That distinction is the point: D10 Step 6 warns that the Pauli
w=1 SECTOR is not invariant, and a spectrum read off a non-invariant subspace is the compression's,
not the system's.

What this script does NOT do: it does not take D10's closed form on trust for the frequencies (it
diagonalizes the measured block), and it does not assume the graph is a chain (the Laplacian is
built from the bond list, and ring, star and complete are exercised below alongside the path, so
a wrong topology cannot pass). It does NOT pin the other half of Step 6: nothing here would fail
if the claim about the w=1 PAULI sector (non-invariant, 160 against 5 at N=5, its twenty
frequencies an artifact of the compression) were wrong. That half is measured in the arc
`site_resolved_vacuum_block`, not here.

Two things it deliberately checks on the side, both scope fences D10 itself states:
  - the ZZ term is what makes the operator a LAPLACIAN: dropping it returns the adjacency form,
    which is the XY sibling F2b's answer, not this one's;
  - uniform gamma is a scope, not a feature: with a gamma PROFILE the block is still exactly
    closed and the operator is -2i*J*Laplacian - 2*diag(gamma), but the modes no longer share a
    decay rate, which is where D10 Step 4's "all decaying at 2gamma" stops.
"""
from __future__ import annotations
import sys
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

J = 1.0
GAMMA = 0.05
N_MAX = 10                 # the full-4^N route stops near 6; this one is limited by 4^N, not 16^N
TOL = 1e-12


def op_at(op, site, n):
    m = np.array([[1.0 + 0j]])
    for k in range(n):
        m = np.kron(m, op if k == site else I2)
    return m


def hamiltonian(n, bonds, j=J, letters=(X, Y, Z)):
    d = 2 ** n
    h = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        for p in letters:
            h += j * op_at(p, a, n) @ op_at(p, b, n)
    return h


def graph_laplacian(n, bonds):
    a = np.zeros((n, n))
    for (u, v) in bonds:
        a[u, v] = a[v, u] = 1.0
    return np.diag(a.sum(axis=1)) - a


def chain(n):
    return [(i, i + 1) for i in range(n - 1)]


def ring(n):
    return [(i, (i + 1) % n) for i in range(n)]


def star(n):
    return [(0, i) for i in range(1, n)]


def complete(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def dissipator_mask(n, gammas):
    """The full Z-dephasing dissipator on the whole d x d coherence space, as an elementwise
    factor: sum_l gamma_l (Z_l rho Z_l - rho)[a,b] = sum_l gamma_l (z_l(a) z_l(b) - 1) rho[a,b].
    Built over ALL of d x d and not on the block, so that a dissipator which did NOT preserve
    the block would show up in the leak below rather than being assumed away."""
    d = 2 ** n
    z = np.array([[1 - 2 * ((s >> (n - 1 - l)) & 1) for l in range(n)] for s in range(d)],
                 dtype=float)
    return (z * gammas) @ z.T - float(np.sum(gammas))


def measure_block(n, bonds, gammas, j=J, letters=(X, Y, Z)):
    """Apply L to each |0><j| and read the (0,1) block off the images.

    Site 0 sits on the MOST significant bit (op_at's kron order), so the single excitation at
    site k is the computational index 1 << (n-1-k).  Returns (block, leak)."""
    d = 2 ** n
    h = hamiltonian(n, bonds, j, letters)
    zmask = dissipator_mask(n, gammas)
    cols = [1 << (n - 1 - k) for k in range(n)]
    block = np.zeros((n, n), dtype=complex)
    leak = 0.0
    inside = np.zeros((d, d), dtype=bool)
    inside[0, cols] = True
    for k, c in enumerate(cols):
        e = np.zeros((d, d), dtype=complex)
        e[0, c] = 1.0
        img = -1j * (h @ e - e @ h) + zmask * e
        block[:, k] = img[0, cols]
        leak = max(leak, float(np.abs(img[~inside]).max()))
    return block, leak


def measure_named_block(n, bonds, gammas, rows, cols, j=J):
    """Same technique on an arbitrary (rows x cols) coherence block: apply L to each |a><b|
    and read the image back on the block. Returns (block, leak) in the given basis order."""
    d = 2 ** n
    h = hamiltonian(n, bonds, j)
    zmask = dissipator_mask(n, gammas)
    idx = [(a, b) for a in rows for b in cols]
    inside = np.zeros((d, d), dtype=bool)
    for (a, b) in idx:
        inside[a, b] = True
    block = np.zeros((len(idx), len(idx)), dtype=complex)
    leak = 0.0
    for c, (a, b) in enumerate(idx):
        e = np.zeros((d, d), dtype=complex)
        e[a, b] = 1.0
        img = -1j * (h @ e - e @ h) + zmask * e
        block[:, c] = [img[p, q] for (p, q) in idx]
        leak = max(leak, float(np.abs(img[~inside]).max()))
    return block, leak


def main():
    ok = True

    print("D10 Step 3 + Step 4 on the (0,1) block, uniform gamma, Heisenberg chain")
    print(f"{'N':>3} {'leak out of block':>18} {'|block - closed form|':>22} {'|omega - D10|':>14}")
    for n in range(2, N_MAX + 1):
        bonds = chain(n)
        block, leak = measure_block(n, bonds, np.full(n, GAMMA))
        pred = -2j * J * graph_laplacian(n, bonds) - 2 * GAMMA * np.eye(n)
        d_entry = float(np.abs(block - pred).max())
        omega = np.sort(np.abs(np.linalg.eigvals(block).imag))
        omega_d10 = np.sort([4 * J * (1 - np.cos(np.pi * k / n)) for k in range(n)])
        d_omega = float(np.abs(omega - omega_d10).max())
        good = leak <= TOL and d_entry <= TOL and d_omega <= 1e-10
        ok &= good
        print(f"{n:>3} {leak:>18.2e} {d_entry:>22.2e} {d_omega:>14.2e}"
              f"{'' if good else '   FAIL'}")

    print("\nStep 3's operator identity is graph-general, which D10's Steps 4-5 are not: the bond")
    print("count cancels between the ferromagnet and the one-magnon diagonal, so only E_ferro*Id")
    print("- H_1 = 2J*Laplacian survives off the chain (the path spectrum and the Neumann waves")
    print("do not). Checked on the topologies the repo builds:")
    for name, gen in (("ring", ring), ("star", star), ("complete", complete)):
        for n in (4, 5, 6):
            bonds = gen(n)
            block, leak = measure_block(n, bonds, np.full(n, GAMMA))
            pred = -2j * J * graph_laplacian(n, bonds) - 2 * GAMMA * np.eye(n)
            d_entry = float(np.abs(block - pred).max())
            good = leak <= TOL and d_entry <= TOL
            ok &= good
            print(f"  {'PASS' if good else 'FAIL'}  {name:<8} N={n}: leak {leak:.1e}, "
                  f"entry-wise {d_entry:.1e}")

    print("\nfour blocks carry this generator, not one: the (0,1) block, the ferromagnet's mirror")
    print("(N,N-1) in the hole basis, and their two transposes, which carry its conjugate. That is")
    print("why the decay shell's frequencies come four-fold (cavity-modes experiment, Result 4).")
    for n in (4, 5):
        bonds = chain(n)
        gammas = np.full(n, GAMMA)
        M = -2j * J * graph_laplacian(n, bonds) - 2 * GAMMA * np.eye(n)
        top, ones = 2 ** n - 1, [1 << (n - 1 - k) for k in range(n)]
        holes = [top ^ s for s in ones]          # popcount N-1, same site order
        cases = [("(0,1)", [0], ones, M), ("(1,0)", ones, [0], M.conj()),
                 ("(N,N-1)", [top], holes, M), ("(N-1,N)", holes, [top], M.conj())]
        for name, rows, cols, want in cases:
            block, leak = measure_named_block(n, bonds, gammas, rows, cols)
            dev = float(np.abs(block - want).max())
            good = leak <= TOL and dev <= TOL
            ok &= good
            print(f"  {'PASS' if good else 'FAIL'}  N={n} {name:<8} leak {leak:.1e}, deviation "
                  f"from {'M' if want is M else 'conj(M)'} {dev:.1e}")

    print("\nthe ZZ term is what makes it a Laplacian (drop it and the adjacency form returns,")
    print("which is F2b's answer, not this one's)")
    for n in (4, 6):
        bonds = chain(n)
        block, leak = measure_block(n, bonds, np.full(n, GAMMA), letters=(X, Y))
        lap = graph_laplacian(n, bonds)
        adj = np.diag(np.diag(lap)) - lap
        d_adj = float(np.abs(block - (2j * J * adj - 2 * GAMMA * np.eye(n))).max())
        d_lap = float(np.abs(block - (-2j * J * lap - 2 * GAMMA * np.eye(n))).max())
        good = leak <= TOL and d_adj <= TOL and d_lap > 1e-6
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  N={n} XY only: matches adjacency ({d_adj:.1e}), "
              f"differs from Laplacian ({d_lap:.2f})")

    print("\nuniform gamma is a scope, not a feature: with a PROFILE the block is still closed")
    print("and the operator is -2iJ*Laplacian - 2*diag(gamma), but the decay rates come apart")
    rng = np.random.default_rng(20260802)
    for n in (4, 5, 6, 7):
        bonds = chain(n)
        gammas = rng.uniform(0.01, 2.5, size=n)
        block, leak = measure_block(n, bonds, gammas)
        pred = -2j * J * graph_laplacian(n, bonds) - 2 * np.diag(gammas)
        d_entry = float(np.abs(block - pred).max())
        spread = float(np.ptp(np.linalg.eigvals(block).real))
        flat = float(np.ptp(np.linalg.eigvals(
            -2j * J * graph_laplacian(n, bonds) - 2 * GAMMA * np.eye(n)).real))
        good = leak <= TOL and d_entry <= TOL and spread > 1e-3 and flat <= 1e-9
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  N={n} profile: leak {leak:.1e}, entry-wise "
              f"{d_entry:.1e}, Re spread {spread:.4f} against {flat:.1e} at uniform gamma")

    print("\nall checks pass" if ok else "\nFAILURES above")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
