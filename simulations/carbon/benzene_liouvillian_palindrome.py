"""Selected C6 spin-ring comparison for the F1 palindrome.

The Hamiltonian below is an XX+YY ring and the two dissipators are deliberately
chosen model jumps: local ``Z_l`` dephasing and bond-operator ``B_b`` dephasing.
If a site occupation ``n_l = (I-Z_l)/2`` and its density jump have first been
selected, then ``D[n_l] = D[Z_l]/4``.  That algebra does not identify a carbon
degree of freedom, a beta-to-J conversion, a bath, gamma, T2, or Q.

F1 covers the selected XX+YY plus all-site-Z model, so this script checks its
palindrome there.  The bond-jump comparison lies outside that F1 premise.  Neither
calculation classifies a molecular vibrational environment.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
d = 2 ** N

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
eye_d = np.eye(d)


def site(op, l):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == l else I2)
    return m


bonds = [(l, (l + 1) % N) for l in range(N)]    # ring, incl. the wrap-around bond


def bond_op(a, b):                              # XX + YY on one C-C bond
    return site(X, a) @ site(X, b) + site(Y, a) @ site(Y, b)


# Hueckel pi-ring Hamiltonian: XX+YY ring (J = 1), the JW image of free-fermion hopping
H = sum(bond_op(a, b) for a, b in bonds)


def commutator(A):                              # -i[A, .] as a superoperator
    return -1j * (np.kron(A, eye_d) - np.kron(eye_d, A.T))


def dissipator(jump):                           # D[L] = L.L+ - 1/2 {L+ L, .}
    ld_l = jump.conj().T @ jump
    return (np.kron(jump, jump.conj())
            - 0.5 * np.kron(ld_l, eye_d)
            - 0.5 * np.kron(eye_d, ld_l.T))


gamma = 1.0
L_H = commutator(H)
D_site = sum(dissipator(np.sqrt(gamma) * site(Z, l)) for l in range(N))
D_bond = sum(dissipator(np.sqrt(gamma) * bond_op(a, b)) for a, b in bonds)


def palindrome_residual(L, f1_centre=None):
    """Residual for Spec(L) closed under lambda -> 2*centre - lambda. When f1_centre
    is given (the F1 prediction -Sigma_gamma) the strict F1 involution
    lambda -> -lambda - 2*Sigma_gamma is tested directly. When it is None the
    spectrum mean is used: the most generous "is it palindromic about anything"
    test, since any palindromic spectrum has mean = centre. The residual is the
    largest distance from a reflected eigenvalue to the nearest actual one."""
    ev = np.linalg.eigvals(L)
    centre = ev.mean() if f1_centre is None else complex(f1_centre)
    reflected = 2.0 * centre - ev
    resid = max(np.min(np.abs(ev - r)) for r in reflected)
    return centre.real, resid


print(f"=== Q2: selected C{N} XX+YY-ring Liouvillian palindrome (d={d}, L is {d * d}x{d * d}) ===\n")

h_ev = np.linalg.eigvalsh(H)
h_resid = float(np.max(np.abs(np.sort(h_ev) + np.sort(h_ev)[::-1])))
print(f"closed-system sanity: many-body H spectrum palindromic about 0, "
      f"residual {h_resid:.2e}  (Coulson-Rushbrooke / truly-class check)\n")

sigma_gamma = N * gamma          # F1 centre is -Sigma_gamma for the Z-dephasing case

for name, L, f1_centre in [
    ("(A) selected local-Z jump       D[Z_l]   F1 premise",
     L_H + D_site, -sigma_gamma),
    ("(B) selected bond jump          D[B_b]   outside F1 premise",
     L_H + D_bond, None),
]:
    centre, resid = palindrome_residual(L, f1_centre)
    verdict = "PALINDROME HOLDS" if resid < 1e-6 else "palindrome BROKEN"
    if f1_centre is None:
        where = f"reflected about the spectrum mean {centre:.4f} (no F1 centre)"
    else:
        where = f"strict F1 involution about -Sigma_gamma = {centre:.4f}"
    print(f"{name}")
    print(f"    {where}    residual = {resid:.2e}    -> {verdict}\n")
