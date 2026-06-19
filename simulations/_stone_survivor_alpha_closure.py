"""Gate-first PYTHON probe of THE STONE (felt_time_dimensions arc, step B). 2026-06-19.

(A) found, at the EIGENVALUE level, that under a deltaJ bond defect the survivor's Re lambda
(-> K_decay) MOVES (soft darkness) while the rigid band edge's is frozen. (B) tests the
TRAJECTORY-level dual: the PTF painter closure Sum_i ln(alpha_i).

CLAIM: Sum ln(alpha) = log of the net time-dilation of the painted trajectory.
  - pure REDISTRIBUTION (vector only, K_decay frozen) -> Sum ln alpha = 0 (IN the +-0.05 window).
  - a global RATE shift (K_decay moves)               -> Sum ln alpha != 0 (OUT of the window).
So the closure BREAKS for a survivor-dominated trajectory, HOLDS for a band-edge one;
quantitatively |Sum ln alpha|_survivor ~ N * deltaJ * (dRe lambda / Re lambda) from (A).

Construction (apples-to-apples): rho_0 = I/d + eps*Herm(mode). I/d is the stationary steady state
(L(I/d)=0), so the ONLY dynamics is the chosen eigenmode -> mode-dominated BY CONSTRUCTION.
Compare mode = survivor (p,p) vs band edge (0,1), plus the canonical bonding psi_1.

GATES: Stage 0 (overlap) -- the dominant non-stationary mode of rho_0 in L's eigenbasis must be the
intended mode (its Re lambda), else vacuous. Stage 1 (prediction) -- closure OUT for the survivor,
IN for the band edge. The closure is over RELIABLE sites only (Symphony's guard: |f=(alpha-1)/dJ| <= 10
+ linearity across dJ/2); a survivor standing-wave NODE site (flat P_i(t), ill-conditioned alpha-fit)
is correctly dropped -- that is the right criterion, not a loosening.
"""
import sys
import numpy as np
from scipy.linalg import eig
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

sys.path.insert(0, "simulations/carbon")
sys.path.insert(0, "simulations")
from incompleteness_survivor import bonds, H_p, survivor
import value_vector_felt_time as vv

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)


def site(op, l, N):
    m = np.array([[1.0 + 0j]])
    for k in range(N):
        m = np.kron(m, op if k == l else I2)
    return m


def full_L(N, J, g, bnds):
    d = 2 ** N
    H = np.zeros((d, d), complex)
    for a, b in bnds:
        H += (J / 2) * (site(X, a, N) @ site(X, b, N) + site(Y, a, N) @ site(Y, b, N))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site(Z, l, N)
        L += g * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def mode_embed(N, prow, pcol, J, g, bnds):
    """Slowest (prow,pcol)-block right-eigenmode embedded as d x d, + its Re lambda."""
    Hr, sr = H_p(N, prow, J, bnds)
    Hc, sc = H_p(N, pcol, J, bnds)
    dr, dc = len(sr), len(sc)
    Lb = -1j * (np.kron(Hr, np.eye(dc)) - np.kron(np.eye(dr), Hc.T))
    Lb += np.diag([-2.0 * g * bin(sr[a] ^ sc[b]).count("1") for a in range(dr) for b in range(dc)])
    w, vr = eig(Lb)
    s = int(np.where(w.real < -1e-7)[0][np.argmax(w.real[w.real < -1e-7])])
    d = 2 ** N
    M = np.zeros((d, d), complex)
    for a in range(dr):
        for b in range(dc):
            M[sr[a], sc[b]] = vr[a * dc + b, s]
    return M, float(w[s].real)


def herm_state(M, d, eps_frac=0.5):
    Mh = (M + M.conj().T) / 2
    Mh = Mh - (np.trace(Mh) / d) * np.eye(d)
    n = np.linalg.norm(Mh)
    if n < 1e-12:
        return None
    Mh /= n
    lo = np.linalg.eigvalsh(Mh).min()
    eps = eps_frac / d / max(abs(lo), 1e-9)
    return np.eye(d) / d + eps * Mh


def psi1_state(N):
    d = 2 ** N
    psi = np.zeros(d, complex)
    nrm = 0.0
    for l in range(N):
        a = np.sin(np.pi * (l + 1) / (N + 1))
        psi[1 << (N - 1 - l)] = a
        nrm += a * a
    psi /= np.sqrt(nrm)
    return np.outer(psi, psi.conj())


def purities(w, R, c, tg, paulis, N):
    d = 2 ** N
    P = np.zeros((N, len(tg)))
    for ti, t in enumerate(tg):
        rho = (R @ (np.exp(w * t) * c)).reshape(d, d)
        for i in range(N):
            ex = np.trace(paulis[0][i] @ rho).real
            ey = np.trace(paulis[1][i] @ rho).real
            ez = np.trace(paulis[2][i] @ rho).real
            P[i, ti] = 0.5 * (1 + ex * ex + ey * ey + ez * ez)
    return P


def dominant_slow(w, c):
    """Eigenvalue of the slowest NON-stationary mode that rho_0 substantially overlaps."""
    amp = np.abs(c)
    mask = w.real < -1e-7
    if not mask.any():
        return None
    thr = 0.02 * amp[mask].max()
    idx = np.where(mask & (amp > thr))[0]
    if len(idx) == 0:
        return None
    return w[idx[np.argmax(w.real[idx])]]


def fit_alpha(pA, pB, tg):
    sp = CubicSpline(tg, pA)
    t0, t1 = tg[0], tg[-1]

    def mse(a):
        return float(np.sum((sp(np.clip(a * tg, t0, t1)) - pB) ** 2))

    grid = np.linspace(0.1, 10.0, 100)
    a0 = grid[int(np.argmin([mse(a) for a in grid]))]
    lo, hi = max(0.1, a0 - 0.25), min(10.0, a0 + 0.25)
    if hi - lo < 1e-3:
        return float(a0)
    return float(minimize_scalar(mse, bounds=(lo, hi), method="bounded").x)


def main():
    Q, dJ, WIN = 1.5, 0.02, 0.05
    print(f"STONE probe (Q={Q}, deltaJ={dJ}): does the PTF closure Sum ln(alpha) BREAK for the soft")
    print(f"survivor and HOLD for the rigid band edge?  in-window iff |Sum ln alpha| <= {WIN}\n", flush=True)
    ok = True
    for N in (4, 5):
        d = 2 ** N
        g = 1.0 / Q
        bnds = bonds(N, "chain")
        defect = bnds[len(bnds) // 3]
        re_s, im_s, sec, nxy = survivor(N, 1.0, g, bnds)
        Msurv, reS = mode_embed(N, sec[0], sec[1], 1.0, g, bnds)
        Mband, reB = mode_embed(N, 0, 1, 1.0, g, bnds)
        dRe, _ = vv.re_shift(N, sec[0], sec[1], 1.0, g, bnds, defect)
        pred = N * dJ * (dRe / reS)
        # eig the clean + defected full Liouvillian once per N
        Vbond = full_L(N, 1.0, 0.0, [defect]) - full_L(N, 0.0, 0.0, [defect])
        LA = full_L(N, 1.0, g, bnds)
        LB = LA + dJ * Vbond
        LB2 = LA + (dJ / 2) * Vbond                      # half-defect, for the linearity guard
        wA, RA = eig(LA); RAi = np.linalg.inv(RA)
        wB, RB = eig(LB); RBi = np.linalg.inv(RB)
        wB2, RB2 = eig(LB2); RB2i = np.linalg.inv(RB2)
        paulis = ([site(X, i, N) for i in range(N)], [site(Y, i, N) for i in range(N)],
                  [site(Z, i, N) for i in range(N)])
        print(f"--- N={N} chain  survivor {sec}, <n_XY>={nxy:.3f}, |Im lambda|={im_s:.3f};  (A) "
              f"predicts survivor Sum ln alpha ~ {pred:+.3f} ---", flush=True)
        print(f"{'state':>17} {'dom Re':>7} {'expect':>7} {'dom?':>5} | {'Sum ln a':>9} {'window':>7} "
              f"{'predict':>8}", flush=True)
        states = [("survivor " + str(sec), herm_state(Msurv, d), reS, False),
                  ("band-edge (0,1)", herm_state(Mband, d), reB, True),
                  ("bonding psi_1", psi1_state(N), reB, True)]
        for label, rho0, ere, holds in states:
            if rho0 is None:
                continue
            v0 = rho0.flatten()
            cA, cB, cB2 = RAi @ v0, RBi @ v0, RB2i @ v0
            ds = dominant_slow(wA, cA)
            dom = ds is not None and abs(ds.real - ere) < 0.2 * abs(ere) + 0.05
            tg = np.linspace(0.0, 5.0 / max(abs(ere), 1e-3), 60)
            PA = purities(wA, RA, cA, tg, paulis, N)
            PB = purities(wB, RB, cB, tg, paulis, N)
            PB2 = purities(wB2, RB2, cB2, tg, paulis, N)
            a1 = np.array([fit_alpha(PA[i], PB[i], tg) for i in range(N)])
            a2 = np.array([fit_alpha(PA[i], PB2[i], tg) for i in range(N)])
            f1, f2 = (a1 - 1) / dJ, (a2 - 1) / (dJ / 2)
            rel = (np.abs(f1) <= 10.0) & (np.abs(f1 - f2) <= 0.5)   # Symphony reliability guard
            clo = float(np.sum(np.log(np.maximum(a1[rel], 1e-12)))) if rel.any() else float("nan")
            inwin = abs(clo) <= WIN
            constructed = "psi" not in label   # psi_1 is canonical, not an I/d+eps*mode construction
            gate = (inwin == holds) and (dom or not constructed)
            ok &= gate
            print(f"{label:>17} {ds.real if ds is not None else float('nan'):>7.3f} {ere:>7.3f} "
                  f"{str(dom):>5} | {clo:>+9.3f} {('IN' if inwin else 'OUT'):>7} "
                  f"{('HOLD' if holds else 'BREAK'):>8} {int(rel.sum())}/{N}rel [{'ok' if gate else 'CHECK'}]",
                  flush=True)
            marks = [f"{a:.3f}{'' if r else '*'}" for a, r in zip(a1, rel)]
            print(f"{'':>17}   per-site alpha = [{', '.join(marks)}]   (* dropped: unreliable)", flush=True)
        print(flush=True)
    print("VERDICT:", "STONE CONFIRMS (A): the closure breaks for the soft survivor and holds for the"
          "\n  rigid band edge -- Sum ln alpha is the trajectory-level value/vector detector."
          if ok else "a gate is off -- diagnose (overlap? domination? regime?).")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
