"""The Hamiltonian dynamics on the rungs of the one diagonal: L_H is an EVEN-STEP ladder.

The dynamics complement to the one-diagonal SYMMETRY (mirror_inventory_d4.py / the just-landed
ThreeDephasingDiagonalsOrbitClaim). The one diagonal is the disagreement count k = popcount(i^j)
(L_D = gamma*(Q - N*I), PROOF_ABSORPTION_THEOREM sect.4.7); its integer levels are the "rungs".

THE QUESTION: how does L_H = -i[H,.] connect rung k <-> k'? A single hop is X_a X_b + Y_a Y_b -- it
flips TWO bits (an excitation moves a->b), so on a coherence |i><j| it changes k = popcount(i^j) by
0 or +-2 (each flipped bit toggles its XOR-membership: +-1 each); ZZ is diagonal (Delta k = 0). So:

  HYPOTHESIS (Stage-0 gate): L_H connects rung k <-> {k, k+-2} -- EVEN steps, NEVER k+-1. Hence the
  PARITY of the disagreement count is CONSERVED by the whole Lindbladian (L_D diagonal, L_H even).

Consequences verified:
  Stage 1: the band edge |vac><magnon| sits at k=1 (ODD); the {0,2}-coherence handover/survivor is
           EVEN -- DIFFERENT parity sectors, so they cannot hybridize (a level crossing, not an EP-
           coalescence): the structural reason the handover is a crossing.
  Stage 2: the mirror R (ket-flip, R|i><j| = |i><j_bar|) maps rung k -> N-k -- the palindrome is "the
           disagreement count read from the other end" (ON_THE_ONE_DIAGONAL), linking the even-step
           L_H dynamics to the mirror group (the symmetry weld) and the F87 hard/soft palindrome.

Convention: row-stacking (C-order) vec, |i><j| -> e_i (x) e_j, kron(A,B): rho -> A rho B^T (matches
framework.lindblad + the one-diagonal verifier). Self-validating; a Stage-0 failure (a k+-1 block, or
parity broken) is THE FINDING. H = J*sum(XX+YY) + Delta*sum(ZZ), open chain (U(1)-conserving).
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import site_op  # noqa: E402

TOL = 1e-12


def popcount(x):
    return bin(x).count("1")


def build_H(N, Delta, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for b in range(N - 1):
        H += (J * (site_op(N, b, "X") @ site_op(N, b + 1, "X")
                   + site_op(N, b, "Y") @ site_op(N, b + 1, "Y"))
              + Delta * (site_op(N, b, "Z") @ site_op(N, b + 1, "Z")))
    return H


def rung_labels(N):
    """k = popcount(i^j) for coherence basis index m = i*d + j (row-stacking)."""
    d = 2 ** N
    return np.array([popcount(i ^ j) for i in range(d) for j in range(d)])


def L_H_super(H, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))   # vec(Hrho - rhoH); kron(A,B): rho -> A rho B^T


def L_D_super(N, gamma):
    """L_D|i><j| = -2*gamma*popcount(i^j)*|i><j| (uniform Z-dephasing; diagonal, Delta k = 0)."""
    return np.diag(-2.0 * gamma * rung_labels(N).astype(complex))


# ============================ Stage 0: the gate ============================
def stage0_even_step_ladder(Ns=(3, 4), Delta=0.7, J=1.0, gamma=0.05):
    print("Stage 0 (gate): L_H connects rungs by EVEN steps only; k-parity CONSERVED:")
    for N in Ns:
        H = build_H(N, Delta, J)
        LH = L_H_super(H, N)
        r = rung_labels(N)
        nz = np.argwhere(np.abs(LH) > TOL)
        dks = sorted({abs(int(r[m]) - int(r[n])) for m, n in nz})
        assert set(dks) <= {0, 2}, (
            f"N={N}: L_H connects rungs by Delta k in {dks} -- expected only {{0,2}}. THE FINDING.")
        assert 1 not in dks and 3 not in dks, f"N={N}: L_H has an ODD rung step {dks} -- parity broken!"
        # parity conservation: the FULL L = L_H + L_D is block-diagonal in k mod 2
        L = LH + L_D_super(N, gamma)
        par = r % 2
        cross = np.abs(L[np.ix_(par == 0, par == 1)]).max()
        assert cross < TOL, f"N={N}: even<->odd k blocks mix (max {cross:.1e}) -- parity NOT conserved."
        print(f"   N={N}: L_H rung steps = {dks} (even only); even<->odd block max = {cross:.1e}; "
              f"k-parity CONSERVED. OK")
    print("   [0] L_H is an even-step rung ladder (k<->k,k+-2); disagreement-count parity is conserved. OK")


# ============================ Stage 1: the sector split ============================
def stage1_sector_split(N=3):
    print(f"\nStage 1: band edge (k=1, odd) vs the {{0,2}} survivor (even) -- different parity sectors (N={N}):")
    d = 2 ** N
    r = rung_labels(N)
    # band edge |vac><magnon|: i=0 (vacuum), j = a single excitation (bit 0) -> m = 0*d + 1
    k_be = int(r[0 * d + 1])
    # {0,2} sector: population |0><0| (k=0); a 2-disagreement coherence |0><3| (j = bits 0,1) (k=2)
    k_pop = int(r[0])
    k_two = int(r[0 * d + 3])
    print(f"   band edge |vac><1magnon| : k = {k_be} ({'odd' if k_be % 2 else 'even'})")
    print(f"   population |vac><vac|     : k = {k_pop} ({'odd' if k_pop % 2 else 'even'})")
    print(f"   2-disagreement coherence : k = {k_two} ({'odd' if k_two % 2 else 'even'})")
    assert k_be == 1 and k_pop == 0 and k_two == 2
    assert k_be % 2 == 1 and k_pop % 2 == 0 and k_two % 2 == 0, "expected band edge odd, {0,2}-survivor even"
    print("   [1] the band edge lives in the ODD sector, the {0,2} survivor in the EVEN sector -- they")
    print("       cannot hybridize (parity-protected), so the handover is a level CROSSING, not a coalescence. OK")


# ============================ Stage 2: the mirror reads the rungs backwards ============================
def stage2_mirror_reads_backwards(N=3):
    print(f"\nStage 2: the mirror R (ket-flip) maps rung k -> N-k (the palindrome, read from the other end) (N={N}):")
    d = 2 ** N
    r = rung_labels(N)
    F = np.array([[1.0]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    for _ in range(N):
        F = np.kron(F, X)                         # F = X^(x)N
    R = np.kron(np.eye(d, dtype=complex), F)       # R: |i><j| -> |i><j|F = |i><j_bar|
    # R permutes coherence index m -> m'; verify r[m'] = N - r[m] for every m
    worst = 0
    for m in range(d * d):
        mp = int(np.argmax(np.abs(R[:, m])))       # R sends basis m to basis m' (a permutation)
        worst = max(worst, abs(int(r[mp]) - (N - int(r[m]))))
    print(f"   max | rung(R m) - (N - rung(m)) | over all {d * d} coherences = {worst}")
    assert worst == 0, "R does not map k -> N-k"
    print("   [2] R maps every rung k -> N-k: the palindrome IS the disagreement count read from the other")
    print("       end. The even-step L_H dynamics and the mirror group meet on the one diagonal. OK")


if __name__ == "__main__":
    stage0_even_step_ladder()
    stage1_sector_split()
    stage2_mirror_reads_backwards()
    print("\nALL STAGES PASS: L_H is an even-step rung ladder (k<->k,k+-2); disagreement-count parity is")
    print("conserved; band edge (odd) and the {0,2} survivor (even) are parity-split (-> level crossing);")
    print("the mirror reads the rungs backwards (k -> N-k). The dynamics complement to the symmetry weld.")
