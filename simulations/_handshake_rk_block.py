"""handshake R_k on the (1,1) single-excitation block (2026-06-19). The big unlock.

The decoder's BondingMode carrier |psi_1> is pure single excitation; <X_i>=<Y_i>=0, so the per-site purity
P_a = 1/2(1+<Z_a>^2) is a function of the POPULATION n_a(t) ALONE. The dynamics is the N^2-dim (1,1)
Liouvillian (single-particle Haken-Strobl): L = -i[H_1, .] + dephasing(-4g on single-exc off-diagonals).
This reproduces the full-4^N painter for this carrier EXACTLY (cross-checked vs _handshake_carrier_compare
N=4,5) but costs N^2, so it runs to N=20+. That kills the small-sample (2-distinct-bond) caveat on the
strength channel and gives the location channel at high N.

CHANNELS (both falsifiable):
  STRENGTH = mean_a f_a(b) (uniform level). Hypothesis: proportional to <psi_1|V_b|psi_1>=2 psi_1(b)psi_1(b+1)
             = d eps_1/dJ_b (Hellmann-Feynman energy sensitivity). corr over ALL bonds (N>=6 -> 3+ distinct).
  LOCATION = f_a(b) - mean (deviation). Hypothesis under test: the eigenvector-perturbation footprint
             sum_{k>=2} c_k(b) R_k(a), R_k(a) ~ psi_1(a)psi_k(a)/(eps_1-eps_k). Reported, not assumed.
"""
import importlib.util
import sys

import numpy as np
from scipy.linalg import eig

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


st = _load("stone", "simulations/_stone_survivor_alpha_closure.py")


def psi(k, a, N):
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (a + 1) / (N + 1))


def eps(k, N, J):
    return 2.0 * J * np.cos(np.pi * k / (N + 1))


def hopping(N, J, bnds):
    H = np.zeros((N, N))
    for a, b in bnds:
        H[a, b] = H[b, a] = J
    return H


def L11(H1, g):
    """(1,1)-sector Liouvillian for a given N x N single-excitation hopping H1. Basis position a = site a
    (bit a). Off-diagonal single-exc coherences a!=b decay at -4g (Hamming weight 2); populations don't."""
    N = H1.shape[0]
    s1 = [1 << a for a in range(N)]
    L = -1j * (np.kron(H1, np.eye(N)) - np.kron(np.eye(N), H1.T))
    L += np.diag([-2.0 * g * bin(s1[a] ^ s1[b]).count("1") for a in range(N) for b in range(N)])
    return L


def populations(w, R, c, tg, N):
    """n_a(t) = rho(t)[a,a] for the (1,1) block trajectory."""
    n = np.zeros((N, len(tg)))
    for ti, t in enumerate(tg):
        rho = (R @ (np.exp(w * t) * c)).reshape(N, N)
        n[:, ti] = np.real(np.diag(rho))
    return n


def paint_block(N, J, g, dJ, bnds):
    """Per-site f_a=(alpha_a-1)/dJ for every bond defect, via (1,1)-block population dynamics."""
    psi1 = np.array([psi(1, a, N) for a in range(N)])
    rho0 = np.outer(psi1, psi1).astype(complex)
    v0 = rho0.flatten()
    H1 = hopping(N, J, bnds)
    LA = L11(H1, g)
    wA, RA = eig(LA)
    cA = np.linalg.inv(RA) @ v0
    ds = st.dominant_slow(wA, cA)
    re = ds.real if ds is not None else -4.0 * g
    tg = np.linspace(0.0, 5.0 / max(abs(re), 1e-3), 60)
    nA = populations(wA, RA, cA, tg, N)
    PA = 0.5 * (1 + (1 - 2 * nA) ** 2)

    nb = len(bnds)
    f = np.full((nb, N), np.nan)
    rel = np.zeros((nb, N), bool)
    for bi, b in enumerate(bnds):
        def f_at(frac):
            Hb = H1.copy()
            Hb[b[0], b[1]] += dJ * frac
            Hb[b[1], b[0]] += dJ * frac
            wB, RB = eig(L11(Hb, g))
            cB = np.linalg.inv(RB) @ v0
            nB = populations(wB, RB, cB, tg, N)
            PB = 0.5 * (1 + (1 - 2 * nB) ** 2)
            return np.array([st.fit_alpha(PA[a], PB[a], tg) for a in range(N)]), PB

        a1, _ = f_at(1.0)
        a2, _ = f_at(0.5)
        f1, f2 = (a1 - 1) / dJ, (a2 - 1) / (dJ / 2)
        f[bi] = f1
        rel[bi] = (np.abs(f1) <= 10.0) & (np.abs(f1 - f2) <= 0.5)
    return f, rel, float(re)


def main():
    Q, J, dJ = 1.5, 1.0, 0.02
    g = 1.0 / Q
    print("=== handshake R_k, (1,1)-block painter: strength + location channels vs N (Q=1.5 chain) ===\n",
          flush=True)
    print(f"{'N':>3} {'#bond':>5} {'rate':>7} {'strength~<v|V|v> corr':>21} {'loc~eigvec corr':>16} "
          f"{'reliable':>9}", flush=True)
    for N in (4, 5, 6, 7, 8, 9):
        bnds = bonds(N, "chain")
        nb = len(bnds)
        f, rel, re = paint_block(N, J, g, dJ, bnds)
        # strength channel
        strength = np.array([f[bi][rel[bi]].mean() if rel[bi].any() else np.nan for bi in range(nb)])
        eshift = np.array([2 * psi(1, b[0], N) * psi(1, b[1], N) for b in bnds])
        sb = ~np.isnan(strength)
        cs = float(np.corrcoef(strength[sb], eshift[sb])[0, 1]) if sb.sum() > 2 else float("nan")
        n_distinct = len(set(np.round(eshift, 6)))
        # location channel: eigenvector-perturbation footprint, per-bond centered
        loc = np.zeros((nb, N))
        for bi, b in enumerate(bnds):
            for a in range(N):
                loc[bi, a] = sum((psi(k, b[0], N) * psi(1, b[1], N) + psi(k, b[1], N) * psi(1, b[0], N))
                                 * psi(1, a, N) * psi(k, a, N) / (eps(1, N, J) - eps(k, N, J))
                                 for k in range(2, N))
        m = rel.copy()

        def center(M):
            out = np.zeros_like(M, float)
            for bi in range(nb):
                if m[bi].any():
                    out[bi] = M[bi] - M[bi][m[bi]].mean()
            return out

        fd, ld = center(f), center(loc)
        cl = float(np.corrcoef(fd[m], ld[m])[0, 1]) if m.sum() > 2 else float("nan")
        print(f"{N:>3} {nb:>5} {re:>7.3f} {cs:>+18.3f}({n_distinct}d) {cl:>+16.3f} {int(rel.sum()):>6}/{nb*N}",
              flush=True)
    print("\nstrength corr uses ALL bonds; (Nd) = # distinct |eshift| values (need >=3 for a non-trivial "
          "corr). loc corr = per-bond-centered f vs the eigenvector footprint (the OPEN channel).", flush=True)


if __name__ == "__main__":
    main()
