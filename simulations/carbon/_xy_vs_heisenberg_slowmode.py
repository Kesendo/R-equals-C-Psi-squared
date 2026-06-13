"""Decisive reconciliation: is the carbon/benzene XY slow mode the same object as the
CHAIN_GAP_SECTOR_DIAGNOSTIC half-filling magnon-admixture? NO - that diagnostic is HEISENBERG (XXX,
with ZZ); the carbon work is XY (free-fermion, the Hueckel hopping model, no ZZ).

This script builds BOTH Hamiltonians through the SAME canonical framework lindbladian_z_dephasing
and lists the eigenvalues nearest zero, so the ONLY difference is the ZZ coupling. It answers two
questions at once:

  1. Does the Heisenberg ring N=6 reproduce CHAIN_GAP's slow mode at Re=-0.230 (sector (3,3))?
     (confirms the discrepancy is the ZZ term, not a bug in my builder.)
  2. Where is the XY ring N=6 slowest non-kernel mode? Is it the band-edge coherence / frozen
     double-excitation at the -2gamma scale (my committed benzene finding), or is there a slower
     near-stationary diffusion mode I missed?

Run it; it prints, it does not assert (this is a diagnostic, not a claim verifier)."""
import sys
import numpy as np

sys.path.insert(0, "simulations")
from framework.lindblad import lindbladian_z_dephasing

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)


def site(op, l, N):
    m = np.array([[1.0 + 0j]])
    for k in range(N):
        m = np.kron(m, op if k == l else I2)
    return m


def H_xy(N, J, ring):
    """(J/2) sum (XX + YY) - the carbon/benzene convention (hopping amplitude J)."""
    bonds = [(b, b + 1) for b in range(N - 1)] + ([(N - 1, 0)] if ring else [])
    H = np.zeros((2 ** N, 2 ** N), complex)
    for a, b in bonds:
        H += (J / 2) * (site(X, a, N) @ site(X, b, N) + site(Y, a, N) @ site(Y, b, N))
    return H


def H_heis(N, J, ring):
    """(J/4) sum (XX + YY + ZZ) - the CHAIN_GAP / framework heisenberg_graph_h convention."""
    bonds = [(b, b + 1) for b in range(N - 1)] + ([(N - 1, 0)] if ring else [])
    H = np.zeros((2 ** N, 2 ** N), complex)
    for a, b in bonds:
        H += (J / 4) * (site(X, a, N) @ site(X, b, N)
                        + site(Y, a, N) @ site(Y, b, N)
                        + site(Z, a, N) @ site(Z, b, N))
    return H


def nearest_zero(L, N, k=12, tol=1e-7):
    """The k non-kernel eigenvalues nearest the imaginary axis, with their dominant filling sector."""
    w, V = np.linalg.eig(L)
    nz = np.where(w.real < -tol)[0]
    order = nz[np.argsort(-w[nz].real)]  # closest to 0 first
    d = 2 ** N
    pc = np.array([bin(i).count("1") for i in range(d)])
    idxp = [np.where(pc == a)[0] for a in range(N + 1)]
    out = []
    for j in order[:k]:
        wt = np.abs(V[:, j].reshape(d, d)) ** 2
        blocks = np.zeros((N + 1, N + 1))
        for a in range(N + 1):
            for b in range(N + 1):
                blocks[a, b] = wt[np.ix_(idxp[a], idxp[b])].sum()
        blocks /= blocks.sum()
        a, b = np.unravel_index(np.argmax(blocks), blocks.shape)
        out.append((w[j].real, abs(w[j].imag), (int(a), int(b)), blocks[a, b]))
    return out


def report(N, J, g, ring):
    Q = J / g
    topo = "ring" if ring else "chain"
    print(f"\n{'='*78}\n{topo} N={N}, J={J}, gamma={g}  (Q=J/gamma={Q:.4g})\n{'='*78}")
    for name, H in (("XY  (J/2)(XX+YY)   [carbon/benzene]", H_xy(N, J, ring)),
                    ("HEIS (J/4)(XX+YY+ZZ) [CHAIN_GAP]   ", H_heis(N, J, ring))):
        L = lindbladian_z_dephasing(H, [g] * N)
        rows = nearest_zero(L, N)
        re0, im0, sec0, w0 = rows[0]
        print(f"\n  {name}")
        print(f"    SLOWEST non-kernel: Re={re0:+.5f}  |Im|={im0:.4f}  "
              f"dominant sector={sec0} ({w0*100:.0f}%)   [Q-AbsThm <n_XY>={-re0/(2*g):.4f}]")
        print("    nearest-zero spectrum:")
        for re, im, sec, wmax in rows[:8]:
            print(f"       Re={re:+.5f}  |Im|={im:.4f}  sector={sec} ({wmax*100:.0f}%)")


if __name__ == "__main__":
    # The CHAIN_GAP anchor: ring N=6, J=1, gamma=0.5 (their Q=2). Heisenberg row says (3,3) Re=-0.230.
    report(6, 1.0, 0.5, ring=True)
    # The benzene V-Effect-seam anchor: ring N=6 below the beat (Q=1.6 -> g=0.625).
    report(6, 1.0, 1.0 / 1.6, ring=True)
