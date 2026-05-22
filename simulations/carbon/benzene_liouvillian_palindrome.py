"""Question 2 (docs/carbon, open question): does benzene's open-system Liouvillian
show the F1 palindrome under vibrational dephasing?

BENZENE_HUCKEL_FRAMEWORK_LENS.md settled the closed system (the Hueckel MO spectrum
is palindromic, Coulson-Rushbrooke = F1 on H). Question 2 is the OPEN system: build
the benzene-pi + vibrational-dephasing Liouvillian and check the F1 palindrome
Spec(L) closed under lambda -> -lambda - 2*Sigma_gamma.

Benzene C6 pi-system in the framework qubit picture: the Hueckel pi-hopping is the
XX+YY ring Hamiltonian (the Jordan-Wigner image of free-fermion hopping). The
carbon-substrate subtlety is the vibrational coupling, and Q2's answer bifurcates
on it:

  (A) Holstein / on-site coupling: the phonon couples to the local pi-density n_l.
      The dephasing dissipator is D[n_l]; since D[a*Z + b*I] = a^2 * D[Z], we get
      D[n_l] = D[(I - Z_l)/2] = (1/4) * D[Z_l]. That IS the framework's Z-dephasing,
      so F1 (proven for XY-on-any-graph + Z-dephasing) forces the palindrome.

  (B) Peierls / SSH / bond coupling: the phonon couples to the C-C bond and
      modulates the hopping. The jump operator is the bond operator B_b itself.
      This is NOT Z-dephasing; F1 does not cover it. Whether the palindrome
      survives is the open part of Question 2 and is what this script tests.

Tom 2026-05-22: try to answer Q2; a null result is also fine.
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


print(f"=== Q2: benzene-ring Liouvillian palindrome, N={N} (d={d}, L is {d * d}x{d * d}) ===\n")

h_ev = np.linalg.eigvalsh(H)
h_resid = float(np.max(np.abs(np.sort(h_ev) + np.sort(h_ev)[::-1])))
print(f"closed-system sanity: many-body H spectrum palindromic about 0, "
      f"residual {h_resid:.2e}  (Coulson-Rushbrooke / truly-class check)\n")

sigma_gamma = N * gamma          # F1 centre is -Sigma_gamma for the Z-dephasing case

for name, L, f1_centre in [
    ("(A) Holstein on-site dephasing  D[Z_l]   = framework Z-dephasing, F1 applies",
     L_H + D_site, -sigma_gamma),
    ("(B) Peierls bond dephasing      D[B_b]   not Z-dephasing, F1 does not cover",
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
