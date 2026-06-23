"""Gate-first PYTHON probe of THE STONE (felt_time_dimensions arc, step B). 2026-06-19.

(A) found, at the EIGENVALUE level, that under a deltaJ bond defect the survivor's Re lambda
(-> K_decay) MOVES (soft darkness) while the rigid band edge's is frozen. (B) tests the
TRAJECTORY-level dual: the PTF painter closure Sum_i ln(alpha_i).

CLAIM (precise, two-lens reviewed 2026-06-19): for the mode-isolating probe rho_0 = I/d + eps*Herm(mode),
the closure Sum_i ln(alpha_i) reads the chosen mode's first-order RATE shift -- OUT of the +-0.05 window
for the soft survivor (Re moves) and IN for the rigid band edge (Re frozen). So (B) is a CONSTRUCTIVE
confirmation of (A) FOR THIS PROBE, not a universal trajectory law.

SCOPE (review-pinned, do not overstate): (i) the break is PROBE-STATE-DEPENDENT -- a survivor-dominated
but POLARIZED state HOLDS (the closure collapses ~0.16->0.007 under added single-site polarization while
the survivor stays the dominant slow mode); the closure measures how cleanly the probe isolates the mode
in the SINGLE-SITE PURITY, which the near-stationary I/d+eps probe maximizes. (ii) Sum ln alpha != 0 does
NOT by itself imply a rate shift -- an asymmetric per-site redistribution also breaks it; the break
certifies a K_decay shift only when the reliable per-site f=(alpha-1)/dJ are SIGN-COHERENT (asserted below
as `coh`). (iii) the magnitude is SCALING + SIGN, not a quantitative law (measured/predicted drifts ~2.5x
toward the handover as biorth->0.28 makes the first-order dRe overshoot). (iv) psi_1 is a fast multi-mode
BASELINE (dom Re ~ -1.7), not a band-edge contrast; the genuine band-edge HOLD is the I/d+eps band-edge row.

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
        f_pred = dRe / reS   # the per-site first-order f if the rate shift were uniform
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
        print(f"--- N={N} chain  survivor {sec}, <n_XY>={nxy:.3f}, |Im lambda|={im_s:.3f};  (A): "
              f"f_pred=dRe/reS={f_pred:+.3f} per site (Sum ~ {pred:+.3f}) ---", flush=True)
        print(f"{'state':>21} {'dom Re':>7} {'dom?':>5} | {'Sum ln a':>9} {'win':>4} {'coh':>5} "
              f"{'kind':>13}  gate", flush=True)
        states = [("survivor " + str(sec), herm_state(Msurv, d), reS, False),
                  ("band-edge (0,1)", herm_state(Mband, d), reB, True),
                  ("psi_1 (fast baseline)", psi1_state(N), reB, True)]
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
            # REVIEW fix (math lens): sign-coherence of the reliable per-site f. coh~1 = all one sign
            # (a genuine global RATE shift); coh~0 = sign-mixed (a mere redistribution, which ALSO breaks
            # Sum ln alpha). The break certifies a K_decay shift only when coh is high.
            rf = f1[rel]
            denom = float(np.sum(np.abs(rf)))
            coh = float(abs(rf.sum()) / denom) if rel.any() and denom > 0 else float("nan")
            inwin = abs(clo) <= WIN
            kind = "frozen" if inwin else ("rate-shift" if coh > 0.8 else "redistrib.")
            constructed = "psi" not in label
            # survivor must BREAK *and* be sign-coherent (a real rate shift); the others must HOLD
            gate_clo = inwin if holds else (not inwin and coh > 0.8)
            gate = gate_clo and (dom or not constructed)
            ok &= gate
            print(f"{label:>21} {ds.real if ds is not None else float('nan'):>7.3f} {str(dom):>5} | "
                  f"{clo:>+9.3f} {('IN' if inwin else 'OUT'):>4} {coh:>5.2f} {kind:>13}  "
                  f"[{'ok' if gate else 'CHECK'}]", flush=True)
            marks = [f"{a:.3f}{'' if r else '*'}" for a, r in zip(a1, rel)]
            print(f"{'':>21}   per-site alpha = [{', '.join(marks)}]   (* dropped: unreliable)", flush=True)
        print(flush=True)
    print("VERDICT:", "STONE CONFIRMS (A) FOR THIS PROBE: for the mode-isolating I/d+eps*Herm(mode) state,"
          "\n  the closure reads the mode's first-order rate shift -- OUT + sign-coherent (rate-shift) for the"
          "\n  soft survivor, IN for the rigid band edge. NOT a universal trajectory law (see SCOPE docstring)."
          if ok else "a gate is off -- diagnose (overlap? domination? sign-coherence? regime?).")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
