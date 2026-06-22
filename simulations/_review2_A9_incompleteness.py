"""
Independent from-below verifier for review finding A9 (INCOMPLETENESS_PROOF.md, Section 2,
Candidate 1: "[Pi^2, L] = 0  =>  noise cannot originate internally").

The dispute is NOT whether [Pi^2,L]=0 (everyone agrees it does). It is whether that
TRUE algebraic identity ENTAILS the physical conclusion "noise must be external."

Refute-first plan: show [Pi^2,L]=0 is
  (1) ORIGIN-agnostic: it holds for the STANDARD EXTERNAL-noise L we always use, so it
      is satisfied by external origin and cannot rule out internal origin; and
  (2) AXIS-agnostic: it holds equally for X-, Y-, and Z-dephasing, so it does not even
      single out the noise axis -> underdetermination, not impossibility.
Both => it is a SYMMETRY/FORM constraint on L, carrying no origin information. Non-sequitur.

Uses the framework's CANONICAL Pi (symmetry.build_pi_full), not a hand-rolled one.
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from framework.symmetry import build_pi_full
from framework.pauli import _vec_to_pauli_basis_transform, _k_to_indices, PAULI_LABELS

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
PMAT = {'X': X, 'Y': Y, 'Z': Z}

def op_at(Pmat, l, N):
    out = np.array([[1]], complex)
    for i in range(N):
        out = np.kron(out, Pmat if i == l else I2)
    return out

def heisenberg_H(N, J=1.0):
    H = np.zeros((2**N, 2**N), complex)
    for l in range(N - 1):
        for Pm in (X, Y, Z):
            H += J * op_at(Pm, l, N) @ op_at(Pm, l + 1, N)
    return H

def vec_adH(H):
    d = H.shape[0]; I = np.eye(d)
    return -1j * (np.kron(I, H) - np.kron(H.T, I))   # vec(-i[H,rho]), column-stacking

def vec_dephasing(letter, gamma, N):
    d = 2**N; I = np.eye(d)
    D = np.zeros((d*d, d*d), complex)
    Pm = PMAT[letter]
    for l in range(N):
        Pl = op_at(Pm, l, N)
        D += gamma * (np.kron(Pl.T, Pl) - np.kron(I, I))   # gamma (P rho P - rho)
    return D

N = 3
gamma = 0.05
H = heisenberg_H(N)
M_basis = _vec_to_pauli_basis_transform(N)
Pi = build_pi_full(N, 'Z')         # canonical palindrome operator (Z convention)
Pi2 = Pi @ Pi
def to_pauli(Lvec):
    return (M_basis.conj().T @ Lvec @ M_basis) / (2**N)

print("=" * 70)
print("FACT (premise is TRUE): [Pi^2, L] = 0 for Heisenberg H + dephasing")
print("=" * 70)
print("Pi^2 diagonal in Pauli basis (a sign per string)?",
      bool(np.allclose(Pi2, np.diag(np.diag(Pi2)))))
adH = vec_adH(H)
results = {}
for letter in ['Z', 'X', 'Y']:
    D = vec_dephasing(letter, gamma, N)
    Lp = to_pauli(adH + D)
    Dp = to_pauli(D)
    comm = Pi2 @ Lp - Lp @ Pi2
    diss_diag = bool(np.allclose(Dp, np.diag(np.diag(Dp))))
    # immune count: damped Pauli strings have nonzero diagonal; immune have ~0
    diss_diag_vals = np.real(np.diag(Dp))
    immune = int(np.sum(np.abs(diss_diag_vals) < 1e-9))
    results[letter] = (np.linalg.norm(comm), diss_diag, immune)
    print(f"  {letter}-dephasing: ||[Pi^2, L]|| = {np.linalg.norm(comm):.2e}   "
          f"dissipator Pauli-diagonal: {diss_diag}   immune strings: {immune}/{4**N}")

print("\n" + "=" * 70)
print("THE NON-SEQUITUR, demonstrated from below")
print("=" * 70)
print("(1) ORIGIN-agnostic: [Pi^2,L]=0 holds for the STANDARD EXTERNAL Z-dephasing L")
print("    (the L we always use, noise admittedly external). A property satisfied by")
print("    external-origin noise cannot, by itself, forbid internal origin.")
print(f"(2) AXIS-agnostic: it holds for Z AND X AND Y dephasing "
      f"(all ||[Pi^2,L]||={results['Z'][0]:.1e}/{results['X'][0]:.1e}/{results['Y'][0]:.1e}).")
print("    Parity does not even fix the noise AXIS, let alone its origin.")
print("MECHANISM: every dephasing dissipator is DIAGONAL in the Pauli basis, and Pi^2 is")
print("    a sign on each Pauli string -> two diagonal operators commute trivially.")
print("    So [Pi^2,L]=0 is a STATEMENT ABOUT L's FORM (parity block-structure), not")
print("    about where the noise comes from. The origin conclusion is unentailed.")

print("\n" + "=" * 70)
print("SLIPPAGE check: per-letter C=1/2 vs GLOBAL immune fraction (1/2)^N")
print("=" * 70)
print("The 'C=1/2' reading is the PER-LETTER 2:2 Klein split (immune {I,Z} : damped {X,Y}).")
print("The GLOBAL Liouville-space immune fraction is (1/2)^N (strings all-letters-in-{I,Z}):")
for n in range(1, 7):
    frac = 0.5 ** n
    marker = "  <- N=2 is 0.25, NOT 0.5" if n == 2 else ""
    print(f"  N={n}: immune fraction = {frac:.4f}{marker}")
# cross-check against the measured immune count above for Z at N=3:
print(f"  (measured Z immune at N={N}: {results['Z'][2]}/{4**N} = {results['Z'][2]/4**N:.4f} "
      f"= (1/2)^{N} = {0.5**N:.4f})")
print("\nDONE.")
