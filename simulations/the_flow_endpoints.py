#!/usr/bin/env python3
"""Two separately scoped objects: a toy 2x2 EP and the full flow's semisimple stationary kernel.

This probe verifies both objects without asserting that one is an endpoint of the other:

  - the EP (parameter-space, at Q_EP): DEFECTIVE. The 2-level L_eff - lambda_EP*I has rank 1
    (a Jordan block; the two eigenvectors coalesce; the Petermann factor diverges). Where the
    memory is born.
  - the full flow's target (state-space): L is singular and the 1/N single-excitation state sits
    in its lambda=0 kernel. The global zero eigenvalue is N+1-fold semisimple (geometric = algebraic),
    one stationary state per particle-number sector; it is not a simple eigenvalue globally.
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


def bond_op(N, b, P, Qop):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == b else (Qop if i == b + 1 else I2))
    return o


def liouvillian_dimensionless(N, Q):
    """L' at canonical Q=J/gamma with H=(Q/2)sum(XX+YY), in column-stack vec."""
    d = 2 ** N
    Id = np.eye(d)
    H1 = sum(bond_op(N, b, X, X) + bond_op(N, b, Y, Y) for b in range(N - 1))
    L = -1j * (Q / 2.0) * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += np.kron(Zl, Zl) - np.kron(Id, Id)
    return L


def ep_endpoint():
    """The EP: the 2-level L_eff(k=1) at Q_EP is DEFECTIVE (rank 1 = a Jordan block)."""
    g0, k, g_eff = 1.0, 1, 4.0 / 3.0
    J = (2.0 / g_eff) * g0                                  # Q_EP = 2/g_eff
    L = np.array([[-2 * g0 * (2 * k - 1), 1j * J * g_eff],
                  [1j * J * g_eff,        -2 * g0 * (2 * k + 1)]], dtype=complex)
    lam = -4.0 * g0 * k
    M = L - lam * np.eye(2)
    rank = int(np.linalg.matrix_rank(M, tol=1e-9))
    print("EP (birth)  , parameter-space singularity:")
    print(f"  rank(L_eff - lambda_EP*I) = {rank}  -> {'DEFECTIVE (Jordan block, eigenvectors coalesce)' if rank == 1 else 'not rank 1'}")


def target_endpoint():
    """The target: lambda=0 stationary state in an N+1-fold semisimple global kernel."""
    print("\nfull-flow stationary kernel:")
    for N in [3, 4, 5, 6]:
        d = 2 ** N
        L = liouvillian_dimensionless(N, 2.0)               # Q = 2
        w = np.linalg.eigvals(L)
        n_zero = int(np.sum(np.abs(w) < 1e-9))              # algebraic multiplicity of lambda=0
        ker = d * d - int(np.linalg.matrix_rank(L, tol=1e-9))  # geometric multiplicity
        rho = np.zeros((d, d), complex)
        for b in range(N):
            psi = np.array([1], complex)
            for i in range(N):
                psi = np.kron(psi, np.array([0, 1], complex) if i == b else np.array([1, 0], complex))
            rho += np.outer(psi, psi.conj()) / N
        res = float(np.linalg.norm(L @ rho.flatten(order="F")))
        kind = "SEMISIMPLE (non-defective)" if ker == n_zero else "defective"
        print(f"  N={N}: ker L = {ker} = N+1 (one fixed point per number sector);  "
              f"L.vec(1/N) = {res:.0e} (fixed point);  lambda=0 {kind} (geom {ker} = alg {n_zero})")


def main():
    print("=" * 74)
    print("TWO SEPARATELY SCOPED SPECTRAL OBJECTS")
    print("=" * 74)
    ep_endpoint()
    target_endpoint()
    print("\nThe toy EP is defective. The full flow's zero eigenvalue is N+1-fold semisimple.")
    print("This probe establishes no branch continuation or endpoint relation between them.")


if __name__ == "__main__":
    main()
