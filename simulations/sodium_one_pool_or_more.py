"""One sodium pool or more: the mean-rate identity against published agarose data.

Reads experiments/SODIUM_ONE_POOL_OR_MORE.md. For a single sodium pool under secular Redfield
quadrupolar relaxation, at any spectral density, static splitting, offset or dynamic shift,

    1/T2s = (1/T1s + 1/T1f)/2                       (the mean-rate identity of the sodium page)

Data, read by hand from the open-access PDFs (9.4 T, 154 mM NaCl in agarose, w/w):
  NBM: Reichert, Schepkin, Kleimaier, Zöllner, Schad, NMR Biomed 37 (2024) e5106, Table 2:
       T1s, T1f and T1m for four phase cycles (fit amplitudes fixed to 0.8 / 0.2) and the T1-TQ/SQ ratio;
       fitted with its Eq. 13 (single and triple quantum jointly, shared T1s, T1f).
  MRM: Reichert, Schepkin, Kleimaier, Zöllner, Schad, Magn Reson Med (2024), doi 10.1002/mrm.30107,
       Table 1: T2s, T2f, As, Af from fixed-evolution-time TQTPPI with a refocusing pulse at 26 evolution
       times in 0.4 to 200 ms, "determined in the same way" as its Eq. 11, read here as the single-quantum
       signal alone fitted as a biexponential plus an offset.
Part 1 checks the identity and splits its gap. Part 2 checks the signal model. Parts 3 to 5 ask which
pictures of more than one pool reproduce the published numbers when their signals go through the authors'
fits. Bound pools have J0 >= J1 >= J2 >= 0, as one spectral-density curve of any distribution of correlation
times gives. Times in ms, rates in 1/ms inside the code. Runtime: about seven minutes; --global adds a
differential-evolution check of Part 3.
"""
import sys
import numpy as np
from scipy.optimize import least_squares, differential_evolution

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")

EPS = np.finfo(float).eps
NBM_T1 = {   # (T1s, sT1s, T1f, sT1f) per phase cycle, NBM Table 2, 9.4 T
    "PC0": {2: (49.44, .96, 29.92, .31), 4: (43.79, 1.18, 25.67, .37), 6: (37.50, 1.70, 23.94, .50)},
    "PC1": {2: (49.55, .90, 30.19, .29), 4: (43.62, 1.22, 26.32, .38), 6: (37.51, 1.66, 23.75, .49)},
    "PC2": {2: (49.22, .73, 30.87, .23), 4: (43.92, .98, 25.65, .31), 6: (37.92, .90, 22.90, .27)},
    "PC3": {2: (49.21, .81, 30.87, .26), 4: (43.78, .98, 25.83, .31), 6: (37.78, .86, 23.09, .26)},
}
NBM_T1M_PC2 = {2: (44.10, .48), 4: (38.65, .75), 6: (33.91, .58)}   # NBM Table 2, T1m (their Eq. 17), PC2
NBM_TQ1_PC2 = {2: 0.0532, 4: 0.0867, 6: 0.1587}                       # NBM Table 2, T1-TQ/SQ, PC2
MRM_T2 = {2: dict(T2s=(45.03, 1.46), T2f=(9.04, .13), As=(0.43, .01), Af=(0.61, .01)),
          4: dict(T2s=(39.17, .73), T2f=(5.12, .06), As=(0.42, .01), Af=(0.62, .01)),
          6: dict(T2s=(33.05, .47), T2f=(3.54, .04), As=(0.43, .01), Af=(0.63, .01))}
NBM_FIG6A_T2S = {2: 44.1, 4: 38.1, 6: 32.2}             # NBM Fig. 6A, read off the figure (no error bars)
T_FREE = 54.19   # MRM Table 1, the 0 % sample: T2 = 54.19 +- 0.95 ms, taken as T1 = T2 of free sodium
S_TFREE = 0.95
KEYS = ["T1s", "T1f", "T2s", "T2f", "slow"]


def data_numbers(c):
    """The five independent published numbers per sample: T1s, T1f (PC2), T2s, T2f and the slow
    amplitude fraction As/(As + Af) with its error propagated (As + Af is not fixed by the fit)."""
    T1s, sT1s, T1f, sT1f = NBM_T1["PC2"][c]
    m = MRM_T2[c]
    As, sAs = m["As"]
    Af, sAf = m["Af"]
    slow = As / (As + Af)
    s_slow = ((Af * sAs)**2 + (As * sAf)**2)**0.5 / (As + Af)**2
    return dict(T1s=(T1s, sT1s), T1f=(T1f, sT1f), T2s=m["T2s"], T2f=m["T2f"], slow=(slow, s_slow))


# ---------------- Part 1: the identity and where its gap sits ----------------
def gap(T1s, sT1s, T1f, sT1f, T2s, sT2s, correlated=False):
    R = 1e3 * (1 / T1s + 1 / T1f) / 2
    a, b = 1e3 * sT1s / T1s**2 / 2, 1e3 * sT1f / T1f**2 / 2
    dR = (a + b) if correlated else (a**2 + b**2)**0.5
    R2, dR2 = 1e3 / T2s, 1e3 * sT2s / T2s**2
    return R, R2, R - R2, (dR**2 + dR2**2)**0.5


def part1():
    print("Part 1. The mean-rate identity, 1/T2s = (1/T1s + 1/T1f)/2, against the measured T2s (MRM)")
    for c in (2, 4, 6):
        T1s, sT1s, T1f, sT1f = NBM_T1["PC2"][c]
        T2s, sT2s = MRM_T2[c]["T2s"]
        R, R2, g, s = gap(T1s, sT1s, T1f, sT1f, T2s, sT2s)
        _, _, _, sc = gap(T1s, sT1s, T1f, sT1f, T2s, sT2s, correlated=True)
        pcs = [gap(*NBM_T1[pc][c], T2s, sT2s)[2] for pc in NBM_T1]
        g6 = 1e3 * (1 / T1s + 1 / T1f) / 2 - 1e3 / NBM_FIG6A_T2S[c]
        print(f"  {c}% agarose: T2s for one pool {1e3 / R:6.2f} ms, measured {T2s:5.2f} +- {sT2s:.2f} ms; rates {R:5.2f} vs "
              f"{R2:5.2f} 1/s, gap {g:4.2f} +- {s:4.2f} 1/s = {g / s:3.1f} sigma ({g / sc:3.1f} sigma with the T1s, T1f errors "
              f"fully correlated)")
        print(f"      over the four phase cycles {min(pcs):4.2f} to {max(pcs):4.2f} 1/s; with NBM's own T2s read off "
              f"Fig. 6A ({NBM_FIG6A_T2S[c]} ms) {g6:4.2f} 1/s")
        T1m_eq17 = 1 / (0.8 / T1s + 0.2 / T1f)
        part_a = 1e3 / T1m_eq17 - 1e3 / T2s
        part_b = 0.3 * (1e3 / T1f - 1e3 / T1s)
        T1m, sT1m = NBM_T1M_PC2[c]
        z = (T2s - T1m) / (sT2s**2 + sT1m**2)**0.5
        print(f"      split of the gap: (1/T1m - 1/T2s) = {part_a:+5.2f} 1/s with T1m from Eq. 17 ({T1m_eq17:.2f} ms), "
              f"0.3 (R1f - R1s) = {part_b:5.2f} 1/s;")
        print(f"      T2s - T1m with their tabulated T1m {T1m} +- {sT1m}: {T2s - T1m:+5.2f} +- {(sT2s**2 + sT1m**2)**0.5:.2f} ms "
              f"= {z:+4.1f} sigma (one pool with J1 >= J2 has T1m/2.5 <= T2s <= T1m)")


# ---------------- the signal model: a pool label beside one spin 3/2 ----------------
MS = [1.5, 0.5, -0.5, -1.5]
SPIN = 1.5
Iz = np.diag(MS).astype(complex)
Ip = np.zeros((4, 4), complex)
for i in range(1, 4):
    Ip[i - 1, i] = np.sqrt(SPIN * (SPIN + 1) - MS[i] * (MS[i] + 1))
Im = Ip.conj().T
I4 = np.eye(4)
TQUAD = {0: (3 * Iz @ Iz - SPIN * (SPIN + 1) * I4) / np.sqrt(6), 1: -(Iz @ Ip + Ip @ Iz) / 2,
         -1: (Iz @ Im + Im @ Iz) / 2, 2: Ip @ Ip / 2, -2: Im @ Im / 2}


def superop(jumps, d):
    """Lindblad dissipator, column stacking: vec(A rho B) = (B^T kron A) vec(rho)."""
    Id = np.eye(d)
    L = np.zeros((d * d, d * d), complex)
    for c in jumps:
        cdc = c.conj().T @ c
        L += np.kron(c.conj(), c) - 0.5 * np.kron(Id, cdc) - 0.5 * np.kron(cdc.T, Id)
    return L


def idx(X, a, b):
    return 8 * (4 * X + b) + (4 * X + a)          # the column-stacked index of |X,a><X,b|


ORDER0 = [idx(X, a, a) for X in range(2) for a in range(4)]
ORDER1 = [idx(X, a, a + 1) for X in range(2) for a in range(3)]
T10 = np.array(MS)                                 # rank 1, order 0 (I_z)
T30 = np.array([1, -3, 3, -1.0])                   # rank 3, order 0
T11 = np.array([np.sqrt(3), 2, np.sqrt(3)])        # rank 1, order 1 (I_+)
T31 = np.array([1, -np.sqrt(3), 1.0])              # rank 3, order 1


def voxel(f0, J_a, J_b, k):
    """Two pools: pool 0 (fraction f0, densities J_a = (J0, J1, J2)) and pool 1 (1 - f0, J_b), exchange with
    detailed balance f0 k01 = (1 - f0) k10, k01 + k10 = k. Returns the generator on the 8-level space."""
    P = [np.diag([1, 0]).astype(complex), np.diag([0, 1]).astype(complex)]
    jumps = [np.sqrt(J[abs(m)]) * np.kron(P[X], TQUAD[m])
             for X, J in ((0, J_a), (1, J_b)) for m in (-2, -1, 0, 1, 2) if J[abs(m)] > 0]
    if k > 0:
        up = np.array([[0, 0], [1, 0]], complex)     # |1><0|, out of pool 0
        jumps += [np.sqrt((1 - f0) * k) * np.kron(up, I4), np.sqrt(f0 * k) * np.kron(up.T, I4)]
    return superop(jumps, 8)


def propagate(L, ids, x1, f0, times):
    G = L[np.ix_(ids, ids)].real
    w, V = np.linalg.eig(G)
    x0 = np.concatenate([f0 * x1, (1 - f0) * x1])
    X = (V @ (np.linalg.solve(V, x0)[:, None] * np.exp(np.outer(w, times)))).real
    n = len(ids) // 2
    return [X[s * n:(s + 1) * n] for s in range(2)], np.linalg.cond(V)


def envelopes(f0, J_a, J_b, k, t1, t2):
    """T1 pathway after a deviation along T10: sq1 = 1 - 2 <T10>/<T10>(0), tq1 = <T30>/4; T2 pathway after
    excitation: sq2 = <T11>/<T11>(0), tq2 = <T31>/(2 sqrt 3). One pool gives tq = e^(-R_fast t) - e^(-R_slow t)."""
    L = voxel(f0, J_a, J_b, k)
    (p0, p1), _ = propagate(L, ORDER0, T10, f0, t1)
    sq1 = 1 - 2 * (T10 @ (p0 + p1)) / (T10 @ T10)
    tq1 = (T30 @ (p0 + p1)) / 4.0
    (q0, q1), _ = propagate(L, ORDER1, T11, f0, t2)
    sq2 = (T11 @ (q0 + q1)) / (T11 @ T11)
    tq2 = (T31 @ (q0 + q1)) / (2 * np.sqrt(3))
    return sq1, tq1, sq2, tq2


def hubbard(J, t1, t2, f=1.0):
    """The same envelopes for one pool without exchange, written from Hubbard's rates (the second route)."""
    J0, J1, J2 = J
    R1s, R1f, R2f, R2s = 6 * J2, 6 * J1, 3 * (J0 + J1), 3 * (J1 + J2)
    return (f * (0.8 * np.exp(-R1s * t1) + 0.2 * np.exp(-R1f * t1)), f * (np.exp(-R1f * t1) - np.exp(-R1s * t1)),
            f * (0.6 * np.exp(-R2f * t2) + 0.4 * np.exp(-R2s * t2)), f * (np.exp(-R2f * t2) - np.exp(-R2s * t2)))


def free_J(T=T_FREE):
    return (1 / (6 * T),) * 3


# ---------------- the authors' fits ----------------
T1_GRID = np.concatenate([0.1 * np.arange(1, 651), 65 + 0.8 * np.arange(1, 301)])   # NBM: dense to 65 ms (as their
                                                                                    # text says), then 8x steps
T2_GRIDS = {"linear": np.linspace(0.4, 200, 26), "geometric": np.geomspace(0.4, 200, 26)}   # MRM: 26 times, spacing unstated
TOL = dict(ftol=1e-12, xtol=1e-12, gtol=1e-12)
B_FLOOR = 0.01   # below 1 % of one pool's longitudinal TQ the imposed TQ weight fades with the model's own signal


def fit_eq13(sq1, tq1, guess):
    """NBM Eq. 13 with the physical plus sign between the two longitudinal terms (the printed equation shows a
    minus there): SQ = 1 - 2 A (0.8 e^(-t/T1s) + 0.2 e^(-t/T1f)), TQ = B (e^(-t/T1f) - e^(-t/T1s)). Both orderings
    of the two times are tried and the lower residual kept (the fit has a second minimum with the times swapped)."""
    def res(p):
        A, T1s, T1f, B = p
        return np.concatenate([1 - 2 * A * (0.8 * np.exp(-T1_GRID / T1s) + 0.2 * np.exp(-T1_GRID / T1f)) - sq1,
                               B * (np.exp(-T1_GRID / T1f) - np.exp(-T1_GRID / T1s)) - tq1])
    fits = [least_squares(res, [1.0, a, b, 0.05], bounds=([0, 1, 1, -10], [2, 500, 500, 10]), **TOL)
            for a, b in ((guess[0], guess[1]), (guess[1], guess[0]))]
    return min(fits, key=lambda r: r.cost).x


def fit_mrm(sq2, grid, guess):
    """MRM Eq. 11: Y = A1 e^(-t/T2s) + A2 e^(-t/T2f) + c on the 26 evolution times; both orderings tried."""
    def res(p):
        A1, A2, Ts, Tf, c = p
        return A1 * np.exp(-grid / Ts) + A2 * np.exp(-grid / Tf) + c - sq2
    fits = [least_squares(res, [0.4, 0.6, a, b, 0.0], bounds=([0, 0, 0.5, 0.5, -1], [2, 2, 500, 500, 1]), **TOL)
            for a, b in ((guess[0], guess[1]), (guess[1], guess[0]))]
    A1, A2, Ts, Tf, c = min(fits, key=lambda r: r.cost).x
    if Ts < Tf:
        A1, A2, Ts, Tf = A2, A1, Tf, Ts
    return Ts, Tf, A1 / (A1 + A2)


def evaluate(f0, J_a, J_b, k, c, grid="linear", tqw=1.0):
    """The five numbers the authors' fits return for this voxel, and B1, its longitudinal TQ against one pool with
    the same fitted times. The TQ part of Eq. 13 carries the reported T1-TQ/SQ weight (times tqw), since the measured
    signal had it; below B_FLOOR of one pool's TQ that weight fades with the model's own signal."""
    t2 = T2_GRIDS[grid]
    sq1, tq1, sq2, _ = envelopes(f0, J_a, J_b, k, T1_GRID, t2)
    d = data_numbers(c)
    guess = (d["T1s"][0], d["T1f"][0])
    _, Ts0, Tf0, B = fit_eq13(sq1, tq1, guess)
    scale = tqw * NBM_TQ1_PC2[c] / max(abs(B), B_FLOOR)
    _, T1s, T1f, _ = fit_eq13(sq1, tq1 * scale, guess)
    T2s, T2f, slow = fit_mrm(sq2, t2, (d["T2s"][0], d["T2f"][0]))
    u = np.exp(-T1_GRID / T1f) - np.exp(-T1_GRID / T1s)
    B1 = (u @ tq1) / (u @ u)
    return dict(T1s=T1s, T1f=T1f, T2s=T2s, T2f=T2f, slow=slow), B1


def pulls(ap, c):
    d = data_numbers(c)
    return np.array([(ap[key] - d[key][0]) / d[key][1] for key in KEYS])


# ---------------- Part 2: the checks ----------------
def part2():
    print("\nPart 2. Checks of the signal model")
    rng = np.random.default_rng(20261114)
    import sympy as sp
    sys.path.insert(0, __file__.replace("\\", "/").rsplit("/", 1)[0])
    from sodium_two_mirrors import generator as gated
    J0, J1, J2, wQ, w = sp.Rational(7, 10), sp.Rational(9, 20), sp.Rational(1, 5), sp.Rational(9, 10), sp.Rational(7, 20)
    Jf = (float(J0), float(J1), float(J2))
    # (a) the jump operators and the dissipator with a Hamiltonian, against the gated generator (row stacking there)
    _, _, Lg = gated(sp.Rational(3, 2), J0, J1, J2, wQ, w)
    Lg = np.array(Lg.evalf(30).tolist(), dtype=complex)
    H = float(wQ) * TQUAD[0] + float(w) * Iz
    Lm = superop([np.sqrt(Jf[abs(m)]) * TQUAD[m] for m in (-2, -1, 0, 1, 2)], 4) - 1j * (np.kron(I4, H) - np.kron(H.T, I4))
    rho = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    a = (Lg @ rho.reshape(-1)).reshape(4, 4)
    b = (Lm @ rho.reshape(-1, order="F")).reshape(4, 4, order="F")
    scale = np.abs(Lg).max() * np.abs(rho).max() * 16
    print(f"  (a) one pool with a Hamiltonian against simulations/sodium_two_mirrors.generator, action on a random rho: "
          f"{np.abs(a - b).max():.1e} = {np.abs(a - b).max() / (EPS * scale):.2f} x eps x (max|L| max|rho| d^2)")
    assert np.abs(a - b).max() <= EPS * scale, "the jump operators or the stacking do not match the gated generator"
    # (a2) voxel() itself: its pool-0 block (fraction 1) against the gated generator without a Hamiltonian
    _, _, Lg0 = gated(sp.Rational(3, 2), J0, J1, J2, 0, 0)
    Lg0 = np.array(Lg0.evalf(30).tolist(), dtype=complex)
    Lv = voxel(1.0, Jf, free_J(), 0.0)
    pool0 = [idx(0, aa, bb) for bb in range(4) for aa in range(4)]
    a = (Lg0 @ rho.reshape(-1)).reshape(4, 4)
    b = (Lv[np.ix_(pool0, pool0)] @ rho.reshape(-1, order="F")).reshape(4, 4, order="F")
    print(f"  (a2) voxel(), pool 0, against the gated generator: {np.abs(a - b).max():.1e} = "
          f"{np.abs(a - b).max() / (EPS * scale):.2f} x eps x (max|L| max|rho| d^2)")
    assert np.abs(a - b).max() <= EPS * scale, "voxel() does not match the gated generator"
    # (b) two routes without exchange; error model: n modes summed, each to eps x cond(V)
    Ja, Jb, f0 = (0.031, 0.0055, 0.0034), free_J(), 0.83
    t2 = T2_GRIDS["linear"]
    g = envelopes(f0, Ja, Jb, 0.0, T1_GRID, t2)
    ha, hb = hubbard(Ja, T1_GRID, t2, f0), hubbard(Jb, T1_GRID, t2, 1 - f0)
    cf = (1 - 2 * (ha[0] + hb[0]), ha[1] + hb[1], ha[2] + hb[2], ha[3] + hb[3])
    L = voxel(f0, Ja, Jb, 0.0)
    kappa = max(propagate(L, ORDER0, T10, f0, t2[:1])[1], propagate(L, ORDER1, T11, f0, t2[:1])[1])
    diff = max(np.abs(x - y).max() for x, y in zip(g, cf))
    print(f"  (b) generator against Hubbard's envelopes, no exchange: {diff:.1e} = {diff / (EPS * kappa):.2f} x eps x "
          f"cond(eigenvectors) ({kappa:.2f}); bound 8 (modes summed)")
    assert diff <= 8 * EPS * kappa, "the two routes to the envelopes disagree"
    # (c) exchange: the central-transition cells against the closed two-site matrix; error model: a 2x2 eigensolve
    k = 0.137
    L = voxel(f0, Ja, Jb, k)
    ct = [idx(0, 1, 2), idx(1, 1, 2)]
    others = [i for i in range(64) if i not in ct]
    leak = np.abs(L[np.ix_(others, ct)]).max() + np.abs(L[np.ix_(ct, others)]).max()
    c1a, c1b = 3 * (Ja[1] + Ja[2]), 3 * (Jb[1] + Jb[2])
    M = np.array([[-c1a - (1 - f0) * k, f0 * k], [(1 - f0) * k, -c1b - f0 * k]])
    ev_g = np.sort(np.linalg.eigvals(L[np.ix_(ct, ct)]).real)
    ev_c = np.sort(np.linalg.eigvals(M).real)
    print(f"  (c) exchange: central-transition cells decoupled from the rest exactly ({leak}); eigenvalues {ev_g} against "
          f"the two-site matrix {ev_c}, difference {np.abs(ev_g - ev_c).max():.1e}; bound 4 x eps x max|M|")
    assert leak == 0.0 and np.abs(ev_g - ev_c).max() <= 4 * EPS * np.abs(M).max(), "exchange disagrees"
    # (d) a single pool comes back from the authors' fits (J1 >= J2); bound two orders above the fits' xtol 1e-12
    J = (0.025, 1 / 180, 1 / 240)
    one, _ = evaluate(1.0, J, free_J(), 0.0, 4)
    target = dict(T1s=1 / (6 * J[2]), T1f=1 / (6 * J[1]), T2s=1 / (3 * (J[1] + J[2])), T2f=1 / (3 * (J[0] + J[1])), slow=0.4)
    worst = max(abs(one[key] / v - 1) for key, v in target.items())
    print(f"  (d) one pool through Eq. 13 and Eq. 11 returns its own times and the 0.4 amplitude, worst relative "
          f"deviation {worst:.1e}; bound 1e-10")
    assert worst < 1e-10, "a single pool must come back from the fits"
    # (e) the swapped minimum of Eq. 13: started from the swapped guess, the fit still returns the pool's times
    sq1, tq1, _, _ = hubbard(J, T1_GRID, t2)
    _, T1s, T1f, _ = fit_eq13(1 - 2 * sq1, tq1 * NBM_TQ1_PC2[4], (target["T1f"], target["T1s"]))
    dev = max(abs(T1s / target["T1s"] - 1), abs(T1f / target["T1f"] - 1))
    print(f"  (e) Eq. 13 started from the swapped times returns {T1s:.6f} / {T1f:.6f} ms for {target['T1s']:.6f} / "
          f"{target['T1f']:.6f}; bound 1e-10")
    assert dev < 1e-10, "the fit kept the swapped minimum"


# ---------------- Parts 3 and 4: free + bound ----------------
def bound_J(J2, a, b):
    """Monotone densities: J1 = J2 + a, J0 = J1 + b with J2, a, b >= 0."""
    return (J2 + a + b, J2 + a, J2)


def start_bound(c):
    d = data_numbers(c)
    J2, J1 = 1 / (6 * d["T1s"][0]), 1 / (6 * d["T1f"][0])
    return [J2, max(J1 - J2, 0.0), max(1 / (3 * d["T2f"][0]) - 2 * J1, 0.0)]


def free_bound_chi(q, c, k, grid="linear", tqw=1.0, T=T_FREE, prior=False):
    if prior:
        fF, J2, a, b, T = q
    else:
        fF, J2, a, b = q
    ap, _ = evaluate(1 - fF, bound_J(J2, a, b), free_J(T), k, c, grid, tqw)
    p = pulls(ap, c)
    return np.append(p, (T - T_FREE) / S_TFREE) if prior else p


LO, HI = [0.0, 1e-5, 0.0, 0.0], [0.9, 1.0, 1.0, 5.0]


def fit_free_bound(c, k=0.0, grid="linear", tqw=1.0, starts=None, T=T_FREE, prior=False, fixed_free=None):
    if starts is None:
        starts = [[f] + start_bound(c) for f in (0.0, 0.02, 0.05, 0.1, 0.2, 0.35)]
    best = None
    for s in starts:
        if fixed_free is not None:
            fun = lambda q: free_bound_chi([fixed_free] + list(q), c, k, grid, tqw, T)
            r = least_squares(fun, s[1:4], bounds=(LO[1:], HI[1:]))
            x = [fixed_free] + list(r.x)
        elif prior:
            r = least_squares(free_bound_chi, list(s[:4]) + [T_FREE], bounds=(LO + [40.0], HI + [80.0]),
                              args=(c, k, grid, tqw, T_FREE, True))
            x = list(r.x)
        else:
            r = least_squares(free_bound_chi, s[:4], bounds=(LO, HI), args=(c, k, grid, tqw, T))
            x = list(r.x)
        if best is None or r.cost < best[1]:
            best = (x, r.cost)
    return best[0], 2 * best[1]


def part3(do_global=False):
    print(f"\nPart 3. Free + bound, no exchange (free sodium T1 = T2 = {T_FREE} ms; 5 numbers, 4 parameters, 1 degree "
          f"of freedom; multi-start local search)")
    best = {}
    for c in (2, 4, 6):
        _, chi_one = fit_free_bound(c, fixed_free=0.0)
        q, chi2 = fit_free_bound(c)
        best[c] = q
        J = bound_J(*q[1:])
        print(f"  {c}% agarose: one pool alone chi2 = {chi_one:5.1f}; free + bound: free fraction {q[0]:.3f}, bound pool "
              f"T1s {1 / (6 * J[2]):.2f}, T1f {1 / (6 * J[1]):.2f}, T2f {1 / (3 * (J[0] + J[1])):.2f}, T2s (identity) "
              f"{1 / (3 * (J[1] + J[2])):.2f} ms; chi2 = {chi2:5.1f}")
        print("      pulls: " + "  ".join(f"{key} {v:+.2f}" for key, v in zip(KEYS, free_bound_chi(q, c, 0.0))))
        near = [q] + [[f] + start_bound(c) for f in (0.02, 0.1)]
        for grid, tqw in (("geometric", 1.0), ("linear", 0.5), ("linear", 2.0)):
            _, chi2v = fit_free_bound(c, 0.0, grid, tqw, starts=near)
            print(f"      {grid} evolution times, TQ weight x{tqw}: chi2 = {chi2v:5.1f}")
        for T in (58.0, 65.0):
            _, chi2t = fit_free_bound(c, starts=near, T=T)
            print(f"      free sodium at T1 = T2 = {T} ms: chi2 = {chi2t:5.1f}")
        qp, chi2p = fit_free_bound(c, starts=near, prior=True)
        print(f"      free time fitted with the prior {T_FREE} +- {S_TFREE} ms (6 numbers, 5 parameters): chi2 = {chi2p:5.1f} "
              f"at {qp[4]:.2f} ms")
        if do_global:
            r = differential_evolution(lambda q: float(np.sum(free_bound_chi(q, c, 0.0)**2)),
                                       list(zip(LO, [0.6, 0.02, 0.01, 0.2])), seed=1, maxiter=60, tol=1e-8, polish=True)
            print(f"      differential evolution: chi2 = {r.fun:5.1f} at free fraction {r.x[0]:.3f}")
    return best


def part4(best):
    print("\nPart 4. Free + bound with exchange: chi2 profile over the exchange rate k (the other four parameters refitted "
          "at each k; a search, not a proof of the minimum; fast exchange makes one averaged pool, so the profile "
          "approaches the one-pool chi2 of Part 3)")
    for c in (2, 4, 6):
        row = []
        for kk in (1, 3, 10, 30, 100, 300, 1000):
            starts = [best[c]] + [[f] + start_bound(c) for f in (0.2,)]
            q, chi2 = fit_free_bound(c, kk / 1e3, starts=starts)
            row.append(f"{kk}: {chi2:5.1f} ({q[0]:.2f})")
        print(f"  {c}% agarose, k in 1/s: chi2 (free fraction)  " + "   ".join(row))


# ---------------- Part 5: two bound pools ----------------
def two_bound(q, c):
    p, J2a, aa, ba, J2b, ab, bb = q
    return evaluate(p, bound_J(J2a, aa, ba), bound_J(J2b, ab, bb), 0.0, c)


def part5():
    print("\nPart 5. Two bound pools, no exchange (7 parameters for 5 numbers, not determined), searched with the "
          "longitudinal TQ fraction B1 held at or above a floor (B1 = 1 is one pool's; a search, not the solution set)")
    lo = [0.05, 1e-5, 0, 0, 1e-5, 0, 0]
    hi = [0.95, 1, 1, 5, 1, 1, 5]
    chi_at = {}
    for c in (2, 4, 6):
        d = data_numbers(c)
        Ts, Tf, T2f = d["T1s"][0], d["T1f"][0], d["T2f"][0]
        guided = []
        for p0, TA2, TA1, TB2, TB1 in ((0.88, 0.98 * Ts, 0.98 * Ts, 0.85 * Tf, 0.85 * Tf),
                                       (0.90, 0.99 * Ts, 0.90 * Ts, 0.87 * Tf, 0.63 * Tf),
                                       (0.75, 0.97 * Ts, 0.90 * Ts, 1.0 * Tf, 0.6 * Tf)):
            guided.append([p0, 1 / (6 * TA2), max(1 / (6 * TA1) - 1 / (6 * TA2), 0.0), 0.8 / (3 * T2f),
                           1 / (6 * TB2), max(1 / (6 * TB1) - 1 / (6 * TB2), 0.0), 0.8 / (3 * T2f)])
        row = []
        for floor in (0.0, 0.25, 0.5, 0.75):
            def res(q):
                ap, B1 = two_bound(q, c)
                return np.append(pulls(ap, c), 50 * max(0.0, floor - B1))
            best = None
            for s in guided:
                r = least_squares(res, s, bounds=(lo, hi), x_scale="jac", max_nfev=300)
                if best is None or r.cost < best.cost:
                    best = r
            ap, B1 = two_bound(best.x, c)
            chi2 = float(np.sum(pulls(ap, c)**2))
            Ja, Jb = bound_J(*best.x[1:4]), bound_J(*best.x[4:7])
            tag = "  (no longitudinal TQ: excluded by the measured TQ peak)" if B1 < B_FLOOR else ""
            row.append(f"      B1 >= {floor:.2f}: chi2 {chi2:6.2f} at B1 {B1:.3f}; pool A {best.x[0]:.2f} (T1 parts "
                       f"{1 / (6 * Ja[2]):.1f} / {1 / (6 * Ja[1]):.1f} ms), pool B {1 - best.x[0]:.2f} "
                       f"({1 / (6 * Jb[2]):.1f} / {1 / (6 * Jb[1]):.1f} ms){tag}")
            chi_at[(c, floor)] = chi2
            if floor == 0.25:
                ap2, _ = two_bound(best.x, c)
                tfine = 0.2 * np.arange(1, 1001)
                _, _, _, tq2 = envelopes(best.x[0], Ja, Jb, 0.0, T1_GRID, tfine)
                u2 = np.exp(-tfine / ap2["T2f"]) - np.exp(-tfine / ap2["T2s"])
                row.append(f"      (this solution's transverse TQ against one pool with its fitted T2 times: {(u2 @ tq2) / (u2 @ u2):.3f})")
        print(f"  {c}% agarose:")
        print("\n".join(row))
    return chi_at


# ---------------- Part 6: the longitudinal TQ across the three samples ----------------
NBM_TQ1 = {"PC1": {2: (3.09, .50), 4: (6.00, .93), 6: (11.47, 2.25)},     # NBM Table 2, T1-TQ/SQ in %
           "PC2": {2: (5.32, .48), 4: (8.67, .89), 6: (15.87, 1.50)}}
NBM_TQ2_FIG6C = {2: 23.1, 4: 27.0, 6: 33.8}                                # NBM Fig. 6C, T2-TQ/SQ in %, read off


def part6(chi_at):
    print("\nPart 6. The longitudinal TQ across the samples. For one pool the ratio of the longitudinal to the transverse "
          "TQ amplitude (each over its single-quantum amplitude) is the same in every sample: the times enter only the "
          "fitted shapes, and the pulses after the evolution and the acquisition are shared. Measured (NBM Table 2 and "
          "Fig. 6C, read off):")
    for pc in ("PC1", "PC2"):
        r = {c: NBM_TQ1[pc][c][0] / NBM_TQ2_FIG6C[c] for c in (2, 4, 6)}
        dr = {c: NBM_TQ1[pc][c][1] / NBM_TQ2_FIG6C[c] for c in (2, 4, 6)}
        q = r[6] / r[2]
        dq = q * ((dr[6] / r[6])**2 + (dr[2] / r[2])**2)**0.5
        print(f"  {pc}: ratio {r[2]:.3f}, {r[4]:.3f}, {r[6]:.3f}; relative to 6 %: {r[2] / r[6]:.2f} : {r[4] / r[6]:.2f} : 1; "
              f"6 % over 2 %: {q:.2f} +- {dq:.2f} (one pool: 1)")
    edges = {}
    for c in (2, 4, 6):
        pts = [(f, chi_at[(c, f)]) for f in (0.0, 0.25, 0.5, 0.75)]
        edge = 0.0
        for (f1, x1), (f2, x2) in zip(pts, pts[1:]):
            if x1 <= 2.0 < x2:
                edge = f1 + (f2 - f1) * (2.0 - x1) / (x2 - x1)
        edges[c] = edge
    print(f"  two bound pools (Part 5): longitudinal TQ fraction where chi2 reaches 2, interpolated between the floors: "
          f"{edges[2]:.2f}, {edges[4]:.2f}, {edges[6]:.2f}; relative to 6 %: {edges[2] / edges[6]:.2f} : "
          f"{edges[4] / edges[6]:.2f} : 1 (soft: a local search with a soft floor)")


# ---------------- Part 7: the size a pathway along z would need ----------------
def part7():
    print("\nPart 7. A signal stored along z by imperfect pulses, added to one pool's transverse signal (one pool from the "
          "fitted T1s, T1f with the identity, T2f as measured) and fitted with Eq. 11 on the 26 evolution times:")
    grid = T2_GRIDS["linear"]
    for c in (2, 4, 6):
        T1s, _, T1f, _ = NBM_T1["PC2"][c]
        T2s_one = 2 / (1 / T1s + 1 / T1f)
        T2s_meas, sT2s = MRM_T2[c]["T2s"]
        T2f = MRM_T2[c]["T2f"][0]
        base = 0.6 * np.exp(-grid / T2f) + 0.4 * np.exp(-grid / T2s_one)
        out = []
        for label, Tx in (("at T1s", T1s), ("at T1m", 1 / (0.8 / T1s + 0.2 / T1f))):
            lo, hi = 0.0, 20.0
            for _ in range(60):
                a = 0.5 * (lo + hi)
                Ts, _, slow = fit_mrm(base + a * np.exp(-grid / Tx), grid, (T2s_meas, T2f))
                lo, hi = (a, hi) if Ts < T2s_meas else (lo, a)
            Ts, _, slow = fit_mrm(base + a * np.exp(-grid / Tx), grid, (T2s_meas, T2f))
            out.append(f"{label}: admixture {a:.2f} of the transverse signal, slow fraction then {slow:.2f}")
        T_half = 2 / (1 / T2s_one + 1 / T1s)
        print(f"  {c}% agarose (one pool {T2s_one:.2f} ms, measured {T2s_meas} +- {sT2s} ms): " + "; ".join(out))
        print(f"      stored along z for half the evolution, the whole slow signal would decay at {T_half:.2f} ms "
              f"= {(T_half - T2s_meas) / sT2s:+.1f} sigma from the measured T2s (measured slow fraction "
              f"{data_numbers(c)['slow'][0]:.2f})")


if __name__ == "__main__":
    do_global = "--global" in sys.argv
    part1()
    part2()
    best = part3(do_global)
    part4(best)
    chi_at = part5()
    part6(chi_at)
    part7()
