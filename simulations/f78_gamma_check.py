"""
Gate-first check of the F78 single-body M eigenvalue formula.

PROOF_SVD_CLUSTER_STRUCTURE.md Theorem F78 (lines 41/43) states, for a single-body
H = sum_l c_l P_l with P in {Y,Z} under uniform Z-dephasing gamma_l:

    M_l eigenvalues = +/- 2 c_l gamma_l i        (full M: sum_l eps_l 2 c_l gamma_l i)

The Abstract (added this session) claims instead, citing the Master Lemma that M is
gamma-independent:

    M_l eigenvalues = +/- 2 c_l i                 (full M: sum_l eps_l 2 c_l i)

These differ by a factor gamma_l. The two coincide only at gamma = 1, which would hide
the discrepancy. This verifier builds M = Pi L Pi^-1 + L + 2*sigma*I from first
principles in the Pauli basis at gamma != 1 and reads the eigenvalues directly.

STAGE 0 GATE: at gamma != 1, the numeric eigenvalues must match exactly one of the two
candidate formulas. The one they match is the correct one; the other's gate fires.
"""

import itertools
import numpy as np

I2 = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
P1 = [I2, X, Y, Z]  # index 0,1,2,3 = I,X,Y,Z

# Per-site Pi action as an operator image (with phase):
#   I -> X (1),  X -> I (1),  Y -> iZ (i),  Z -> iY (i)
PI_IMG = [X, I2, 1j * Z, 1j * Y]


def kron_list(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def pauli_string(idx_tuple):
    return kron_list([P1[i] for i in idx_tuple])


def build_basis(N):
    return list(itertools.product(range(4), repeat=N))


def project(op, basis, d):
    """Coefficients of op in the (orthogonal, HS) Pauli basis."""
    coeffs = np.empty(len(basis), dtype=complex)
    for j, idx in enumerate(basis):
        Pb = pauli_string(idx)
        coeffs[j] = np.trace(Pb.conj().T @ op) / d  # Tr(Pa Pb) = d * delta
    return coeffs


def superop_in_pauli_basis(action, basis, d):
    """Matrix of a superoperator (given as a function op->op) in the Pauli basis."""
    n = len(basis)
    Mmat = np.empty((n, n), dtype=complex)
    for col, idx in enumerate(basis):
        Pa = pauli_string(idx)
        Mmat[:, col] = project(action(Pa), basis, d)
    return Mmat


def single_site_op(letter_op, site, N):
    mats = [I2] * N
    mats[site] = letter_op
    return kron_list(mats)


def run(N, c, gamma):
    """c: list of length N (coefficient of Y on each site). gamma: list of length N."""
    d = 2 ** N
    basis = build_basis(N)
    H = sum(c[l] * single_site_op(Y, l, N) for l in range(N))
    Zops = [single_site_op(Z, l, N) for l in range(N)]
    sigma = sum(gamma)

    def L(op):
        out = -1j * (H @ op - op @ H)
        for l in range(N):
            out = out + gamma[l] * (Zops[l] @ op @ Zops[l] - op)
        return out

    def Pi(op):
        # Pi(P_a) = tensor of per-site images (with phase)
        # build as superoperator on the operator directly
        return _pi_apply(op, N)

    Lmat = superop_in_pauli_basis(L, basis, d)
    Pimat = superop_in_pauli_basis(Pi, basis, d)
    Mmat = Pimat @ Lmat @ np.linalg.inv(Pimat) + Lmat + 2 * sigma * np.eye(len(basis))

    eig = np.linalg.eigvals(Mmat)
    nz = np.sort_complex(eig[np.abs(eig) > 1e-9])

    # candidate spectra
    gi = [c[l] for l in range(N)]
    gd = [c[l] * gamma[l] for l in range(N)]
    pred_gamma_indep = sorted(
        (2j * sum(s * gi[l] for l, s in enumerate(signs))
         for signs in itertools.product([1, -1], repeat=N)),
        key=lambda z: (z.imag, z.real))
    pred_gamma_dep = sorted(
        (2j * sum(s * gd[l] for l, s in enumerate(signs))
         for signs in itertools.product([1, -1], repeat=N)),
        key=lambda z: (z.imag, z.real))

    return eig, nz, pred_gamma_indep, pred_gamma_dep


def _pi_apply(op, N):
    # Express op in Pauli basis, map each string by PI_IMG, recombine.
    d = 2 ** N
    basis = build_basis(N)
    out = np.zeros((d, d), dtype=complex)
    for idx in basis:
        Pb = pauli_string(idx)
        coeff = np.trace(Pb.conj().T @ op) / d
        if abs(coeff) < 1e-15:
            continue
        img = kron_list([PI_IMG[i] for i in idx])
        out = out + coeff * img
    return out


def report(N, c, gamma):
    eig, nz, pi, pd = run(N, c, gamma)
    print(f"--- N={N}, c={c}, gamma={gamma} (gamma != 1 so the two formulas differ) ---")
    re_max = np.max(np.abs(eig.real))
    print(f"max |Re(eig)|                = {re_max:.2e}   (Master Lemma: M purely imaginary => ~0)")
    nz_im = np.sort(np.round(nz.imag, 6))
    print(f"numeric nonzero eig (Im)     = {nz_im}")
    print(f"pred  +/-2 c_l i      (Im)   = {np.sort(np.round([z.imag for z in pi], 6))}")
    print(f"pred  +/-2 c_l gamma_l i(Im) = {np.sort(np.round([z.imag for z in pd], 6))}")
    # gate
    set_num = set(np.round(nz.imag, 5))
    set_indep = set(np.round([z.imag for z in pi if abs(z) > 1e-9], 5))
    set_dep = set(np.round([z.imag for z in pd if abs(z) > 1e-9], 5))
    match_indep = set_num == set_indep
    match_dep = set_num == set_dep
    print(f"matches gamma-INDEPENDENT (+/-2c_l i)      : {match_indep}")
    print(f"matches gamma-DEPENDENT   (+/-2c_l gamma i): {match_dep}")
    assert re_max < 1e-9, "GATE: M is not purely imaginary -> Master Lemma broken"
    assert match_indep, "GATE: numeric eig does NOT match +/-2c_l i"
    assert not match_dep, "GATE: numeric eig matches the gamma-dependent formula (proof body)"
    print("=> Abstract (gamma-independent +/-2c_l i) CONFIRMED; proof-body gamma_l is spurious.\n")


if __name__ == "__main__":
    report(1, [0.37], [0.71])
    report(2, [0.37, 0.19], [0.71, 0.53])
    report(3, [0.37, 0.19, 0.41], [0.71, 0.53, 0.29])
    print("All gates passed: F78 single-body M eigenvalues are +/-2 c_l i (gamma-independent).")
