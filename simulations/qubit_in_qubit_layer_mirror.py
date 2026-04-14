"""
Qubit-in-Qubit: Layer Mirror Structure in the Minimal Nest
==========================================================

Minimal two-qubit Lindblad system where dephasing acts only on the outer
qubit (B), and the inner qubit (S) decoheres solely through S-B coupling.
Investigates whether the full Liouvillian's eigenmodes exhibit a
three-class partial-trace structure aligned with the palindromic spectral
symmetry from MIRROR_SYMMETRY_PROOF.

Setup:
    H_SB = J * 0.5 * (X (x) X + Y (x) Y)        (XX+YY coupling)
    L_B  = sqrt(gamma_B) * I (x) Z               (dephasing on B only)
    No direct dephasing on S.

Three separate analyses:
    (1) Non-Markovian signature in reduced S dynamics (coherence rebound)
    (2) Eigenvalue classes and partial-trace weights per mode
    (3) SWAP-operator expectation values per eigenmode

Scope of this script:
    Single parameter point (J=1.0, gamma_B=0.1). Reproduces the numerical
    observations recorded in hypotheses/NESTED_MIRROR_STRUCTURE.md.
    Sanity checks (SWAP-artifact robustness, N=3 scaling, alternative
    couplings) are pending and tracked as open items in the hypothesis doc.

Authors: Tom and Claude (chat)
Date: 2026-04-14
"""

import numpy as np
from scipy.linalg import expm

np.set_printoptions(precision=4, suppress=True)

# Pauli basis
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def liouvillian(H, jumps):
    """Vectorized Liouvillian using column-stacking vec convention."""
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def evolve(rho0, L, t):
    """Evolve rho0 under Liouvillian L to time t."""
    v = rho0.reshape(-1, order='F')
    vt = expm(L * t) @ v
    d = rho0.shape[0]
    return vt.reshape((d, d), order='F')


def partial_trace_B(rho_AB):
    """Trace out B (second qubit) from 2-qubit density matrix."""
    r = rho_AB.reshape(2, 2, 2, 2)
    return np.einsum('ibjb->ij', r)


def partial_trace_S(rho_AB):
    """Trace out S (first qubit) from 2-qubit density matrix."""
    r = rho_AB.reshape(2, 2, 2, 2)
    return np.einsum('aiaj->ij', r)



# Parameters
J = 1.0
gamma_B = 0.1

# Hamiltonian and dissipator
H = J * 0.5 * (kron(X, X) + kron(Y, Y))
L_jump = np.sqrt(gamma_B) * kron(I2, Z)
L_super = liouvillian(H, [L_jump])

# Initial state: S in |+>, B maximally mixed (no prior information in bath)
rho_plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
rho_mix = 0.5 * np.eye(2, dtype=complex)
rho0 = np.kron(rho_plus, rho_mix)


def analysis_1_non_markovian_rebound():
    """Observation 1: rho_S coherence is non-monotonic.

    A single-qubit Lindblad with any fixed gamma_eff produces strictly
    monotonic exponential decay of the off-diagonal element. The two-layer
    reduced dynamics shows rebound, which no single-layer model can
    reproduce.
    """
    print("=" * 72)
    print("Analysis 1: Non-Markovian rebound in rho_S")
    print("=" * 72)
    ts = np.linspace(0, 20, 201)
    coh_2layer = []
    for t in ts:
        rho_t = evolve(rho0, L_super, t)
        rho_S = partial_trace_B(rho_t)
        coh_2layer.append(abs(rho_S[0, 1]))
    # Fit single-layer gamma_eff to match at t=1
    i1 = 10  # ts[10] = 1.0
    gamma_eff = -0.5 * np.log(coh_2layer[i1] / 0.5)
    print(f"  Fitted single-layer gamma_eff at t=1: {gamma_eff:.4f}")
    print()
    print(f"  {'t':>6} {'|S01| two-layer':>18} {'|S01| one-layer':>18} {'diff':>10}")
    for i in [0, 5, 10, 20, 50, 100, 150, 200]:
        c2 = coh_2layer[i]
        c1 = 0.5 * np.exp(-2 * gamma_eff * ts[i])
        print(f"  {ts[i]:6.2f} {c2:18.4f} {c1:18.4f} {c2-c1:+10.4f}")
    print()
    print("  Key observation: at t=10, two-layer |S01|=0.1684, one-layer=0.0018.")
    print("  The rebound is the signature of information flowing from B back to S.")
    print("  No Markovian single-layer model reproduces this trajectory.")
    print()



def analysis_2_eigenmode_classes():
    """Observation 2+3: three eigenvalue classes with distinct partial-trace weights.

    Diagonalize the Liouvillian. Eigenvalues fall into three classes:
      - Re(lambda) = 0     (steady-state manifold)
      - Re(lambda) = -gamma  (coupled/mirror modes)
      - Re(lambda) = -2*gamma  (correlation modes)

    Partial-trace weights per class reveal the layer structure:
      - 0-modes:       |M_S| = |M_B| = 1.0   (shared conserved quantities)
      - -gamma-modes:  |M_S| = |M_B| = 0.707 (equal split, mirror axis)
      - -2*gamma-modes: |M_S| = |M_B| = 0.0  (traceless on both single layers,
                                              pure correlation content)
    """
    print("=" * 72)
    print("Analysis 2+3: Eigenvalue classes and partial-trace weights")
    print("=" * 72)
    print(f"  Expected palindrome pairing: Re(lam_a) + Re(lam_b) = -2*Sum(gamma) = {-2*gamma_B}")
    print()
    eigvals, eigvecs = np.linalg.eig(L_super)
    order = sorted(range(16), key=lambda i: (eigvals[i].real, eigvals[i].imag))
    print(f"  {'Re(lam)':>9} {'Im(lam)':>10}  {'|M_S|':>7} {'|M_B|':>7}  class")
    print("  " + "-" * 66)
    for idx in order:
        lam = eigvals[idx]
        vec = eigvecs[:, idx]
        M = vec.reshape((4, 4), order='F')
        nTotal = np.linalg.norm(M)
        M_S = partial_trace_B(M)
        M_B = partial_trace_S(M)
        nS = np.linalg.norm(M_S) / nTotal
        nB = np.linalg.norm(M_B) / nTotal
        if abs(lam.real) < 1e-6:
            cls = "conserved (0)"
        elif abs(lam.real + 2 * gamma_B) < 1e-6:
            cls = "correlation (-2g)"
        else:
            cls = "mirror (-g)"
        print(f"  {lam.real:9.4f} {lam.imag:10.4f}  {nS:7.4f} {nB:7.4f}  {cls}")
    print()


def analysis_3_swap_structure():
    """Observation 4: eigenmodes are simultaneous SWAP eigenmodes (numerical).

    Caveat: [L, SWAP] != 0 (gamma is only on B, which breaks strict S-B
    exchange symmetry). The ±1 SWAP expectation values below are therefore
    a property of the numpy-selected basis within each degenerate eigenspace
    of L, not a commuting-operator relationship. Verifying whether this is
    a basis artifact or a genuine structural property requires breaking the
    degeneracy (e.g., asymmetric gamma on both qubits) and checking whether
    the ±1 pattern survives. Open item in the hypothesis doc.
    """
    print("=" * 72)
    print("Analysis 3: SWAP expectation values per eigenmode")
    print("=" * 72)
    # SWAP operator on 2 qubits
    SWAP = np.zeros((4, 4), dtype=complex)
    for i in range(2):
        for j in range(2):
            SWAP[2 * j + i, 2 * i + j] = 1
    SWAP_super = np.kron(SWAP.T, SWAP)
    comm = L_super @ SWAP_super - SWAP_super @ L_super
    print(f"  ||[L, SWAP]|| = {np.linalg.norm(comm):.6f}  (nonzero: strict symmetry broken)")
    print()
    eigvals, eigvecs = np.linalg.eig(L_super)
    order = sorted(range(16), key=lambda i: (eigvals[i].real, eigvals[i].imag))
    print(f"  {'Re(lam)':>9} {'Im(lam)':>10}  {'<SWAP>':>10}")
    print("  " + "-" * 40)
    for idx in order:
        lam = eigvals[idx]
        v = eigvecs[:, idx]
        v = v / np.linalg.norm(v)
        swap_exp = v.conj() @ SWAP_super @ v
        print(f"  {lam.real:9.4f} {lam.imag:10.4f}  {swap_exp.real:+7.4f}")
    print()


if __name__ == "__main__":
    analysis_1_non_markovian_rebound()
    analysis_2_eigenmode_classes()
    analysis_3_swap_structure()
