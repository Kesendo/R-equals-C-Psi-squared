#!/usr/bin/env python3
"""Compare three distinct diagnostics at canonical Q=3 and Q=2000.

The slow real-rate tolerance cluster can be multidimensional, so arbitrary eigenvector overlaps
and per-vector commutator residuals are not invariant. This producer compares orthogonal projectors
onto the selected right-invariant subspaces and a Frobenius aggregate of the Hamiltonian commutator
over an orthonormal subspace basis. Spectral-edge rate drift, subspace overlap, and
Hamiltonian-commutator residual are reported side by side; no exact equivalence or causal shortcut
criterion is inferred from this finite profile list.
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def H_xy_unit(N):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        for P in (X, Y):
            t = np.array([[1]], complex)
            for i in range(N):
                t = np.kron(t, P if i in (b, b + 1) else I2)
            H += t
    return H


def L_dimless(N, Q, profile, H1):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (Q / 2.0) * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def slow_subspace(N, Q, profile, H1, cluster_tolerance=1e-6):
    L = L_dimless(N, Q, profile, H1)
    w, V = np.linalg.eig(L)
    nonkernel = np.abs(w) > 1e-7
    edge = float(np.max(np.where(nonkernel, w.real, -np.inf)))
    cluster = np.where(nonkernel & (np.abs(w.real - edge) <= cluster_tolerance))[0]
    selected = V[:, cluster]
    U, singular, _ = np.linalg.svd(selected, full_matrices=False)
    rank_tol = max(selected.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.sum(singular > rank_tol))
    return -edge, U[:, :rank]


def subspace_overlap(left, right):
    """Mean squared canonical correlation; basis invariant, in [0,1]."""
    return float(np.linalg.norm(left.conj().T @ right, "fro") ** 2 / min(left.shape[1], right.shape[1]))


def commutator_residual(basis, H1, N):
    """Basis-invariant aggregate ||C|_V||_F / ||H_left|_V||_F for C(rho)=[H,rho]."""
    d = 2 ** N
    numerator = denominator = 0.0
    for j in range(basis.shape[1]):
        rho = basis[:, j].reshape(d, d, order="F")
        comm = H1 @ rho - rho @ H1
        numerator += np.linalg.norm(comm) ** 2
        denominator += np.linalg.norm(H1 @ rho) ** 2
    return float(np.sqrt(numerator / (denominator + 1e-30)))


def main():
    N = 5
    H1 = H_xy_unit(N)
    profiles = {
        "uniform":           [1, 1, 1, 1, 1],
        "V peaked (edges.25)": [0.25, 0.75, 3.0, 0.75, 0.25],
        "V peaked strong":   [0.1, 0.5, 3.8, 0.5, 0.1],
        "edge flat-bulk":    [0.25, 1.5, 1.5, 1.5, 0.25],
        "center-quiet":      list(np.array([1.333, 1.333, 0.5, 1.333, 0.5]) * 5 / 4.999),
        "inverse-V":         list(np.array([2.4, 0.6, 0.2, 0.6, 2.4]) * 5 / 6.2),
    }
    print(f"N={N}. Distinct rate-drift, slow-subspace-overlap, and commutator diagnostics.\n")
    print(f"  {'profile':>20}  {'edge Q=3':>10}  {'edge Q=2000':>11}  {'subspace ov2':>12}  {'||C|V||/||H|V||':>16}")
    for name, p in profiles.items():
        p = list(np.array(p, float) * N / np.sum(p))
        r_lo, basis_lo = slow_subspace(N, 3.0, p, H1)
        r_hi, basis_hi = slow_subspace(N, 2000.0, p, H1)
        ov = subspace_overlap(basis_lo, basis_hi)
        res = commutator_residual(basis_hi, H1, N)
        print(f"  {name:>20}  {r_lo:>10.5f}  {r_hi:>11.5f}  {ov:>12.6f}  {res:>16.5f}")
    print("\n  reading: these are separate finite-scan observables; their correlation is not an iff theorem.")
    print("  This finite table does not isolate the cause of spectral-edge rate drift.")


if __name__ == "__main__":
    main()
