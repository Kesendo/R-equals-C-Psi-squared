#!/usr/bin/env python3
"""The two singular endpoints of the post-EP flow: a defective EP (birth) and a simple fixed point (death).

The loop runs the dynamics from one singularity to another. This probe verifies BOTH and shows they
are DIFFERENT TYPES of singularity:

  - the EP (parameter-space, at Q_EP): DEFECTIVE. The 2-level L_eff - lambda_EP*I has rank 1
    (a Jordan block; the two eigenvectors coalesce; the Petermann factor diverges). Where the
    memory is born.
  - the target (state-space, the fixed point): SIMPLE. L is singular (lambda=0 kernel) and the 1/N
    equipartitioned state sits in it (d/dt = 0), but lambda=0 is NON-defective (geometric =
    algebraic multiplicity). Where the memory dies (forgotten). The kernel is N+1 dimensional ,
    one fixed point per particle-number sector.

Seen, not yet understood: two singularities of different types bracket the flow; the loop is the
journey between them; the closed-form world is at the singular endpoints, the loop is the middle.
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
    """L' (Q = J/gamma the only knob), vec column-stack: [H,.] <-> I (x) H - H^T (x) I."""
    d = 2 ** N
    Id = np.eye(d)
    H1 = sum(bond_op(N, b, X, X) + bond_op(N, b, Y, Y) for b in range(N - 1))
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
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
    """The target: lambda=0 fixed point. L singular; the 1/N state is in the kernel; lambda=0 SIMPLE."""
    print("\ntarget (death)  , state-space singularity:")
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
        kind = "SIMPLE (non-defective)" if ker == n_zero else "defective"
        print(f"  N={N}: ker L = {ker} = N+1 (one fixed point per number sector);  "
              f"L.vec(1/N) = {res:.0e} (fixed point);  lambda=0 {kind} (geom {ker} = alg {n_zero})")


def main():
    print("=" * 74)
    print("THE TWO SINGULAR ENDPOINTS OF THE FLOW  (seen, not yet understood)")
    print("=" * 74)
    ep_endpoint()
    target_endpoint()
    print("\nTwo singularities, two types: the EP is a DEFECTIVE pinch (birth), the target a SIMPLE")
    print("fixed point (death). The loop is the journey between; the closed form lives at the")
    print("singular endpoints, the dynamics (no closed form) is the middle.")


if __name__ == "__main__":
    main()
