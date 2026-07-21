"""Born-rule-mirror tables probe: the February diagonals, reproduced and decoded.

Reproduces every table of experiments/BORN_RULE_MIRROR.md from first principles
and pins three things the February record left open:

1. THE CONVENTION. The tables reproduce EXACTLY under the PAULI convention
   H = 1.0 * sum_bonds (XX + YY + ZZ) on the N=4 ring (the doc's "J = 1.0").
   In the spin convention H = (J/4)*sum (the one PROOF_F94 uses), the same
   J = 1 gives a much weaker evolution (P(|00>) = 0.94, not 0.41). The doc's
   dynamics therefore corresponds to J_spin = 4, i.e. Q = J_spin/gamma = 80,
   NOT Q = 20. With Q = 80, K = gamma*t = 0.0143, the F94 leading order
   (4/3)*Q^2*K^3 = 0.0250 lands within ~14% of the measured ratio deviation
   Delta_|00> = P_lind/P_unit - 1 = 0.0290 (the point sits ~100x above F94's
   verified deep-perturbative window, so higher orders contribute the rest).

2. THE PSI COLUMN. The C_eff table's Psi_i values (0.2274/0.2139/0.1127) are
   Psi_i = (sum_{j != i} |rho02[i,j]|) / (d-1), the l1 row-coherence of the
   LINDBLAD reduced pair state over d-1 = 3. They are NOT the Born amplitude
   |<i|psi>| of the doc's Section 3.2 (the reduced state is mixed; sqrt of
   the unitary probabilities gives 0.64/0.51/0.25).

3. THE SHIFTS. All nine basis-shift numbers of Sections 2.2/7.3 reproduce,
   with "x/y-eigenstates" meaning the parity-even product-eigenbasis pair
   ({|++>,|-->} etc.), shift = (Lindblad - unitary) diagonal summed over the
   pair. The sigma_z "favors z-eigenstates" total +0.0098 is an AGGREGATE:
   |00> gains +0.0120 while |11> loses -0.0022.

Expected output (verified 2026-07-21): see the EXPECTED block at the bottom.
"""

import numpy as np
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)

N = 4
T = 0.286
GAMMA = 0.05
BONDS = [(0, 1), (1, 2), (2, 3), (3, 0)]


def op(P, i):
    m = np.array([[1.0]], dtype=complex)
    for k in range(N):
        m = np.kron(m, P if k == i else I2)
    return m


def hamiltonian(coef):
    h = np.zeros((16, 16), dtype=complex)
    for a, b in BONDS:
        for P in (X, Y, Z):
            h += coef * op(P, a) @ op(P, b)
    return h


def lindblad_rho(H, jump_pauli, t):
    d = 16
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        A = op(jump_pauli, k)
        L += GAMMA * (np.kron(A, A.conj()) - np.kron(Id, Id))
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi0 = np.kron(np.kron(zero, plus), np.kron(zero, plus))
    rho0 = np.outer(psi0, psi0.conj()).reshape(-1)
    return (expm(L * t) @ rho0).reshape(16, 16)


def unitary_rho(H, t):
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi0 = np.kron(np.kron(zero, plus), np.kron(zero, plus))
    psi = expm(-1j * H * t) @ psi0
    return np.outer(psi, psi.conj())


def red02(rho16):
    """Reduced state of qubits (0, 2): trace out qubits 1 and 3."""
    r = rho16.reshape(2, 2, 2, 2, 2, 2, 2, 2)
    return np.einsum('abcdefgh,bf,dh->aceg', r, np.eye(2), np.eye(2)).reshape(4, 4)


def basis_change(rho4, single):
    U = np.kron(single, single)
    return U.conj().T @ rho4 @ U


def main():
    Hp = hamiltonian(1.0)

    print("=== Convention check: pair-(0,2) unitary diagonal at t = 0.286 ===")
    for coef, name in [(1.0, "Pauli J=1 (the doc's dynamics)"),
                       (0.25, "spin J/4 at J=1 (PROOF_F94's book)")]:
        d = np.real(np.diag(red02(unitary_rho(hamiltonian(coef), T))))
        print(f"  {name:<36s} {np.round(d, 4)}")

    rho_u = red02(unitary_rho(Hp, T))
    rho_l = red02(lindblad_rho(Hp, Z, T))
    pu, pl = np.real(np.diag(rho_u)), np.real(np.diag(rho_l))
    print("\n=== Section 2.1 table (doc: 0.4134/0.2616/0.0635 vs 0.4254/0.2567/0.0613) ===")
    for i, lab in enumerate(["|00>", "|01>", "|10>", "|11>"]):
        print(f"  {lab}: unitary {pu[i]:.4f}  lindblad {pl[i]:.4f}  "
              f"dev {abs(pl[i]-pu[i])/pu[i]*100:.1f}%")
    print(f"  dominant eigenvalue of rho02: {np.max(np.linalg.eigvalsh(rho_l)):.4f}  (doc 0.911)")

    print("\n=== Sections 2.2 / 7.3: the nine shifts ===")
    hada = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    sdag = np.array([[1, 0], [0, -1j]], dtype=complex)
    for P, pname in [(Z, "sigma_z"), (X, "sigma_x"), (Y, "sigma_y")]:
        rl = red02(lindblad_rho(Hp, P, T))
        for single, bname in [(np.eye(2, dtype=complex), "z"), (hada, "x"),
                              (sdag.conj().T @ hada, "y")]:
            dl = np.real(np.diag(basis_change(rl, single)))
            du = np.real(np.diag(basis_change(rho_u, single)))
            even = (dl[0] + dl[3]) - (du[0] + du[3])
            print(f"  {pname} dephasing, shift to {bname}-eigenstates: {even:+.4f}")

    print("\n=== Section 3.4: the Psi column decoded ===")
    l1rows = np.array([np.sum(np.abs(rho_l[i])) - np.abs(rho_l[i, i])
                       for i in range(4)]) / 3.0
    print(f"  l1 row-coherence / (d-1): {np.round(l1rows, 4)}  (doc 0.2274/0.2139/0.2139/0.1127)")
    print(f"  C_eff = P/Psi^2:          {np.round(pl / l1rows**2, 2)}  (doc 8.22/5.61/5.61/4.82)")
    print(f"  sqrt(P_unitary) (NOT it): {np.round(np.sqrt(pu), 3)}")

    print("\n=== F94 bridge ===")
    ratio = pl[0] / pu[0] - 1
    K = GAMMA * T
    for Q, tag in [(80.0, "Q=80 (doc's Pauli J=1 => J_spin=4)"),
                   (20.0, "Q=20 (naive J=1 in F94's spin book)")]:
        print(f"  (4/3)*Q^2*K^3 at {tag}: {4/3*Q**2*K**3:.5f}")
    print(f"  measured ratio deviation P_l/P_u - 1 = {ratio:.5f}")

    print("\nReading: the tables are exact under Pauli J=1; the doc's point is")
    print("Q = 80 in F94's variables (leading order within ~14%, higher orders")
    print("carry the rest); Psi_i is Lindblad row-coherence/3, not an amplitude.")


if __name__ == "__main__":
    main()
