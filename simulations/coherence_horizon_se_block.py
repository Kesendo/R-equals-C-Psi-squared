"""The coherence horizon Q*(N) via the single-excitation (Haken-Strobl) Liouvillian.
Q*(N) is the EP of the N-site dephased tight-binding chain's rho^(1) dynamics, an N^2-dim
object that reproduces the full 4^N horizon. This script derives its closed form / exact
condition. Self-validating: run it, every assert must pass."""
import numpy as np


def h_single(N, J):
    """N x N single-particle Hamiltonian: tridiagonal, off-diagonal J (the (J/2)(XX+YY) hop)."""
    h = np.zeros((N, N), complex)
    for i in range(N - 1):
        h[i, i + 1] = J
        h[i + 1, i] = J
    return h


def L_se(N, J, g):
    """The N^2-dim single-excitation Liouvillian in the real-space basis. Coherent part
    -i[h, rho]; Z-dephasing makes coherence rho[i][j] (i!=j) decay at 4g, populations untouched
    (a single-excitation coherence has popcount(i^j)=2, Re = -2g*2 = -4g, the Absorption Theorem)."""
    h = h_single(N, J)
    I = np.eye(N)
    L = -1j * (np.kron(h, I) - np.kron(I, h.T))
    deph = np.array([(-4.0 * g if i != j else 0.0) for i in range(N) for j in range(N)])
    return L + np.diag(deph)


def qstar_se(N, J=1.0, lo=0.02, hi=6.0):
    """Q* = J/g* where g* is the largest g (smallest Q) at which the slowest non-zero SE mode
    stops oscillating (the EP), by bisection on g."""
    for _ in range(70):
        m = 0.5 * (lo + hi)
        ev = np.linalg.eigvals(L_se(N, J, m))
        nz = ev[ev.real < -1e-7]
        if len(nz) == 0:
            hi = m
            continue
        gap = nz.real.max()
        band = nz[np.abs(nz.real - gap) < 1e-7]
        if np.abs(band.imag).max() > 1e-7:
            lo = m
        else:
            hi = m
    return J / (0.5 * (lo + hi))


# Anchors: the live-witness ladder.
LADDER = {2: 1.0, 3: np.sqrt(2.0), 4: 1.87850, 5: 2.37220}


def _assert_ladder():
    for N, q in LADDER.items():
        got = qstar_se(N)
        assert abs(got - q) < 2e-3, f"SE Q*({N})={got:.5f} != ladder {q:.5f}"
    print("[Task1] SE builder reproduces the witness ladder:",
          {N: round(qstar_se(N), 5) for N in LADDER})


def _L_full(N, J, g):
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1]).astype(complex)

    def site(op, l):
        m = np.array([[1.0 + 0j]])
        for k in range(N):
            m = np.kron(m, op if k == l else I2)
        return m

    d = 2 ** N
    H = np.zeros((d, d), complex)
    for b in range(N - 1):
        H += (J / 2) * (site(X, b) @ site(X, b + 1) + site(Y, b) @ site(Y, b + 1))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site(Z, l)
        L += g * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def _assert_full_reduction():
    # SE spectrum must be a sub-spectrum of the full L (the single-excitation block is invariant).
    for N in (2, 3):
        g = 0.5
        full = np.linalg.eigvals(_L_full(N, 1.0, g))
        sub = np.linalg.eigvals(L_se(N, 1.0, g))
        for s in sub:
            assert np.min(np.abs(full - s)) < 1e-9, f"N={N}: SE eig {s} not in full spectrum"
    print("[Task1] SE spectrum is an exact sub-spectrum of the full Liouvillian (N=2,3).")


if __name__ == "__main__":
    _assert_ladder()
    _assert_full_reduction()
