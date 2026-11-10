"""Evidence, not a gate (docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md, the general triangle). Search the
triangle's rate space for a paired spectrum that no ratio carries: minimise
f = sum_i min_j |lam_i + lam_j + c|^2 / c^2 (the spectrum against its mirror about -c/2), divided by the spread of the
ratios g^2 = sum_{l<k} (log r_l - log r_k)^2, r = gu/gd. Near one ratio f ~ g^2, so f/g^2 stays finite there; a
pairing off one ratio is a zero. Control: the equal triangle (equal |t|, one Delta), where two equal rate points (the mirror)
are zeros with g > 0. Log rates are held in [-2, 2] (unbounded, the search runs off to rates 0 or infinity, where
everything is small). Numpy alone, its own Liouvillian.

Run: python simulations/triangle_pairing_search.py equal|generic STARTS   (the proof's numbers: STARTS = 24)"""
import sys
import numpy as np
from scipy.optimize import minimize
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1])
SM = np.array([[0, 1], [0, 0]], complex); SP = SM.T.copy()

def op(o, l, n=3):
    m = np.eye(1)
    for q in range(n): m = np.kron(m, o if q == l else I2)
    return m

def ham(bonds, n=3):
    H = np.zeros((2 ** n, 2 ** n), complex)
    for a, b, J, D, Dz in bonds:
        H += J * (op(X, a) @ op(X, b) + op(Y, a) @ op(Y, b)) + D * (op(X, a) @ op(Y, b) - op(Y, a) @ op(X, b)) + Dz * op(Z, a) @ op(Z, b)
    return H

def lind(H, gd, gu, n=3):
    d = 2 ** n; Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(n):
        for g, c in ((gd[l], op(SM, l)), (gu[l], op(SP, l))):
            cd = c.conj().T; cc = cd @ c
            L = L + g * (np.kron(c, c.conj()) - 0.5 * (np.kron(cc, Id) + np.kron(Id, cc.T)))
    return L

def objective(q, H):
    p = 2 * np.tanh(q)                     # log rates in a box [-2, 2]
    gd, gu = np.exp(p[:3]), np.exp(p[3:])
    ev = np.linalg.eigvals(lind(H, gd, gu)); c = gd.sum() + gu.sum()
    dist = np.abs(ev[:, None] + ev[None, :] + c).min(axis=1)
    f = (dist ** 2).sum() / c ** 2
    lr = p[3:] - p[:3]
    g2 = sum((lr[i] - lr[j]) ** 2 for i in range(3) for j in range(i + 1, 3))
    return f / g2, f, g2

def search(H, starts, seed):
    rng = np.random.default_rng(seed); best = []
    for _ in range(starts):
        p0 = rng.normal(0, 1, 6)
        r = minimize(lambda p: objective(p, H)[0], p0, method="Nelder-Mead", options={"maxiter": 3000, "xatol": 1e-10, "fatol": 1e-16})
        q, f, g2 = objective(r.x, H)
        best.append((q, f, g2, r.x))
    best.sort(key=lambda t: t[0])
    return best

if __name__ == "__main__":
    which = sys.argv[1]; starts = int(sys.argv[2])
    if which == "equal":
        H = ham([(0, 1, 1.0, 0.4, 0.3), (1, 2, 1.0, 0.4, 0.3), (2, 0, 1.0, 0.4, 0.3)])
    else:
        H = ham([(0, 1, 1.0, 0.0, 0.3), (1, 2, 1.2, 0.5, -0.5), (2, 0, 0.2, 0.8, 0.2)])
    rng = np.random.default_rng(2)
    rand = [objective(np.arctanh(rng.uniform(-2, 2, 6) / 2), H)[0] for _ in range(200)]
    print(f"median of f/g2 at 200 random points of the box: {np.median(rand):.2e}")
    best = search(H, starts, 1)
    edge = sum(1 for _, _, _, x in best if np.abs(2 * np.tanh(x)).max() > 1.999)
    print(f"minima with a log rate at the box edge (|log rate| > 1.999): {edge} of {len(best)}")
    for q, f, g2, x in best:
        p = 2 * np.tanh(x)
        gd, gu = np.exp(p[:3]), np.exp(p[3:])
        print(f"f/g2 {q:.3e}  f {f:.2e}  g2 {g2:.2e}  max|log rate| {np.abs(p).max():.4f} (box edge 2)  "
              f"gd {np.round(gd, 4)}  gu {np.round(gu, 4)}")
