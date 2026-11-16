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
fits. Part 6 compares the longitudinal TQ across the samples, Part 7 the size a pathway along z would need.
Part 8 runs the two sequences pulse by pulse for one pool, against the hand count with Wigner coefficients,
and reads the measured TQ amplitudes against one pool's. Part 9 reads the first paper's sampling of T1 off its
Fig. 2A and off the same figure in the first author's thesis, drawn there as a vector path (the points in
data/reichert2024_sodium_agarose, read by simulations/read_reichert_thesis_figures.py), and checks the thesis's
Fig. 4.11, the paper's Fig. 5 drawn as vector paths, against Table 2. Part 10 sets the thesis's second evaluation
(its Eq. 3.16, Table 4.7 and Fig. 4.14, read by the same script) and its acquisition at fixed evolution times beside
the count and reading (i), and fits the measured FID itself: the paper's fit identified on it, its plateau, the
weights, a control of known split, and the triple-quantum amplitude held. Bound pools have J0 >= J1 >= J2 >= 0, as
one spectral-density curve of any distribution of correlation times gives. Times in ms, rates in 1/ms inside the
code. Runtime: about eight minutes; --global adds a differential-evolution check of Part 3.
"""
import os
import sys
from collections import Counter
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


THESIS_T47 = {   # (T1s, sT1s, T1f, sT1f, T1m, sT1m, TQ/SQ in %, its error) per sequence, the first author's thesis
    # (Heidelberg 2024), Table 4.7, 9.4 T, read by hand
    "standard": {2: (46.05, .26, 41.79, .10, 45.13, .21, 18.52, .32), 4: (41.38, .42, 33.62, .17, 39.55, .32, 19.13, .24),
                 6: (36.33, .24, 28.50, .11, 34.44, .18, 23.36, .20)},   # the paper's data refitted with the TQ amplitude
    # held at the thesis's Eq. 3.16
    "fixed": {2: (46.68, .89, 42.75, 3.22, 45.84, .07, 19.90, 21.61), 4: (40.28, .67, 35.60, 2.34, 39.25, .06, 31.52, 21.17),
              6: (36.04, .38, 29.73, 1.27, 34.57, .04, 29.29, 8.18)},   # a separate acquisition at 37 fixed evolution
    # times, the SQ and TQ curves separated and fitted with the TQ amplitude free
}


def thesis_t47(seq, c, T2s, sT2s):
    """The identity on one row of the thesis's Table 4.7. Its T1m is Eq. 17 of its two times (checked to the table's
    last digit), so its error constrains the correlation rho of the two rate errors, which the split and the gap carry;
    a value beyond -1 or 1, the rounding of the tabulated errors, is printed and then held at the limit."""
    a, sa, f, sf, m, sm, tq, stq = THESIS_T47[seq][c]
    assert abs(1 / (0.8 / a + 0.2 / f) - m) <= 0.005, "T1m is Eq. 17 of the table's times"
    Rs, Rf, Rm = 1e3 / a, 1e3 / f, 1e3 / m
    ss, sf_, sm_ = 1e3 * sa / a**2, 1e3 * sf / f**2, 1e3 * sm / m**2
    rho_raw = (sm_**2 - 0.64 * ss**2 - 0.04 * sf_**2) / (0.32 * ss * sf_)
    rho = float(np.clip(rho_raw, -1, 1))
    R2, sR2 = 1e3 / T2s, 1e3 * sT2s / T2s**2
    split, s_split = Rf - Rs, (ss**2 + sf_**2 - 2 * rho * ss * sf_)**0.5
    gap_, s_gap = (Rs + Rf) / 2 - R2, ((ss**2 + sf_**2 + 2 * rho * ss * sf_) / 4 + sR2**2)**0.5
    pa, spa, pf, spf = NBM_T1["PC2"][c]          # Table 2 (PC2), its errors taken as independent
    split2, s_split2 = 1e3 / pf - 1e3 / pa, 1e3 * ((spa / pa**2)**2 + (spf / pf**2)**2)**0.5
    return dict(T1s=a, T1f=f, rho=rho, rho_raw=rho_raw, split=split, s_split=s_split, T2s_pool=2e3 / (Rs + Rf), gap=gap_, s_gap=s_gap,
                tq=tq, stq=stq, below=(split2 - split) / (s_split**2 + s_split2**2)**0.5)


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
        for seq, label in (("standard", "the paper's data refitted with the TQ amplitude held at the thesis's Eq. 3.16"),
                           ("fixed", "the thesis's acquisition at fixed evolution times, the TQ amplitude free")):
            r = thesis_t47(seq, c, T2s, sT2s)
            print(f"      thesis Table 4.7, {label}: T1s {r['T1s']:.2f}, T1f {r['T1f']:.2f} ms, rho {r['rho_raw']:+.3f} from its "
                  f"T1m error (used: {r['rho']:+.2f}); split {r['split']:4.2f} +- {r['s_split']:4.2f} 1/s; T2s for one pool {r['T2s_pool']:6.2f} ms, gap "
                  f"{r['gap']:4.2f} +- {r['s_gap']:4.2f} 1/s = {r['gap'] / r['s_gap']:3.1f} sigma; the split {r['below']:3.1f} "
                  f"sigma below Table 2's (PC2); TQ/SQ {r['tq']:.2f} +- {r['stq']:.2f} %")


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
T1_GRID = np.concatenate([0.2 * np.arange(1, 321), 64 + 1.8 * np.arange(1, 321)])   # the sampling of NBM Fig. 2A,
# drawn as a vector path in the thesis's Fig. 3.6a (Part 9): 320 points 0.2 ms apart to 64 ms, then 320 points 1.8 ms
# apart to 640 ms; the text's step and counts (100 us, 8 steps of it, 300 coarse points) give neither
T2_GRIDS = {"linear": np.linspace(0.4, 200, 26), "geometric": np.geomspace(0.4, 200, 26)}   # MRM: 26 times, spacing unstated
TOL = dict(ftol=1e-12, xtol=1e-12, gtol=1e-12)
B_FLOOR = 0.01   # below 1 % of one pool's longitudinal TQ the imposed TQ weight fades with the model's own signal


def fit_eq13(sq1, tq1, guess, t=None):
    """NBM Eq. 13 with the physical plus sign between the two longitudinal terms (the printed equation shows a
    minus there): SQ = 1 - 2 A (0.8 e^(-t/T1s) + 0.2 e^(-t/T1f)), TQ = B (e^(-t/T1f) - e^(-t/T1s)). Both orderings
    of the two times are tried and the lower residual kept (the fit has a second minimum with the times swapped)."""
    # (the printed Eq. 14 shows a minus between its two transverse terms as well; this script does not fit Eq. 14)
    t = T1_GRID if t is None else t

    def res(p):
        A, T1s, T1f, B = p
        return np.concatenate([1 - 2 * A * (0.8 * np.exp(-t / T1s) + 0.2 * np.exp(-t / T1f)) - sq1,
                               B * (np.exp(-t / T1f) - np.exp(-t / T1s)) - tq1])
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


def evaluate(f0, J_a, J_b, k, c, grid="linear", tqw=1.0, t1=None):
    """The five numbers the authors' fits return for this voxel, and B1, its longitudinal TQ against one pool with
    the same fitted times. The TQ part of Eq. 13 carries the reported T1-TQ/SQ weight (times tqw), since the measured
    signal had it; below B_FLOOR of one pool's TQ that weight fades with the model's own signal."""
    t2 = T2_GRIDS[grid]
    t1 = T1_GRID if t1 is None else t1
    sq1, tq1, sq2, _ = envelopes(f0, J_a, J_b, k, t1, t2)
    d = data_numbers(c)
    guess = (d["T1s"][0], d["T1f"][0])
    _, Ts0, Tf0, B = fit_eq13(sq1, tq1, guess, t1)
    scale = tqw * NBM_TQ1_PC2[c] / max(abs(B), B_FLOOR)
    _, T1s, T1f, _ = fit_eq13(sq1, tq1 * scale, guess, t1)
    T2s, T2f, slow = fit_mrm(sq2, t2, (d["T2s"][0], d["T2f"][0]))
    u = np.exp(-t1 / T1f) - np.exp(-t1 / T1s)
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
def two_bound(q, c, t1=None):
    p, J2a, aa, ba, J2b, ab, bb = q
    return evaluate(p, bound_J(J2a, aa, ba), bound_J(J2b, ab, bb), 0.0, c, t1=t1)


TWO_LO = [0.05, 1e-5, 0, 0, 1e-5, 0, 0]
TWO_HI = [0.95, 1, 1, 5, 1, 1, 5]


def guided_starts(c):
    """Three starts for two bound pools near the published times (used by Parts 5 and 8)."""
    d = data_numbers(c)
    Ts, Tf, T2f = d["T1s"][0], d["T1f"][0], d["T2f"][0]
    guided = []
    for p0, TA2, TA1, TB2, TB1 in ((0.88, 0.98 * Ts, 0.98 * Ts, 0.85 * Tf, 0.85 * Tf),
                                   (0.90, 0.99 * Ts, 0.90 * Ts, 0.87 * Tf, 0.63 * Tf),
                                   (0.75, 0.97 * Ts, 0.90 * Ts, 1.0 * Tf, 0.6 * Tf)):
        guided.append([p0, 1 / (6 * TA2), max(1 / (6 * TA1) - 1 / (6 * TA2), 0.0), 0.8 / (3 * T2f),
                       1 / (6 * TB2), max(1 / (6 * TB1) - 1 / (6 * TB2), 0.0), 0.8 / (3 * T2f)])
    return guided


def part5():
    print("\nPart 5. Two bound pools, no exchange (7 parameters for 5 numbers, not determined), searched with the "
          "longitudinal TQ fraction B1 held at or above a floor (B1 = 1 is one pool's; a search, not the solution set)")
    lo, hi = TWO_LO, TWO_HI
    chi_at = {}
    for c in (2, 4, 6):
        guided = guided_starts(c)
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
            row.append(f"      B1 >= {floor:.2f}: chi2 {chi2:7.3f} at B1 {B1:.3f}; pool A {best.x[0]:.2f} (T1 parts "
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
    return edges


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


# ---------------- Part 8: the two sequences, pulse by pulse ----------------
IX = (Ip + Im) / 2
IY = (Ip - Im) / 2j
TPPI = np.deg2rad(90 + 45 * np.arange(1, 9))      # NBM Fig. 1, PC0: phase of the first 90 pulse, 90 + 45 n, n = 1..8
NBM_T1M_PC1 = {2: (44.18, .59), 4: (38.92, .79), 6: (34.02, 1.10)}   # NBM Table 2, T1m (their Eq. 17), PC1
MRM_T2STAR = {2: (16.89, 5.35), 4: (21.26, 4.65), 6: (17.48, 3.07)}  # MRM Table 1: T2s*, T2f* of the single-pulse FID
# The hand count, per |d| (d: what the acquisition makes of a rank-3 coherence of order -1 against a rank-1 one).
# Longitudinal: the inverted populations carry rank-3 order 2 f30 = (4/5)(e^-R1s t - e^-R1f t); the first 90 pulse turns
# it into order +-3 with |d3_{3,0}|, the second back to order -1 with |d3_{-1,3}|, against |d1_{1,0}| |d1_{-1,1}| for the
# single quantum: (4/5)(5 sqrt6/16) = sqrt6/4. Transverse: f31 = (sqrt6/5)(...) turned with |d3_{3,1}| |d3_{-1,3}| against
# |d1_{1,1}| |d1_{-1,1}|: (sqrt6/5)(15/16) = 3 sqrt6/16. The longitudinal rank-3 order also reaches order +-1, with
# |d3_{1,0}| |d3_{-1,1}|, at the single-quantum frequency: the slow weight becomes 0.8 + (2/5)(sqrt6/16)|d| = 0.8 + sqrt6/40 |d|
# (its sign, toward the slow weight, is the simulation's). The relaxation factors 4/5 = 2 x 2/5 and sqrt6/5 are the square
# roots of the weights' products, 2/5 = sqrt(4/5 x 1/5) and sqrt6/5 = sqrt(3/5 x 2/5).
LONG, TRANS, SHIFT = np.sqrt(6) / 4, 3 * np.sqrt(6) / 16, np.sqrt(6) / 40
_PULSES = {}


def hard(theta, phi):
    key = (theta, phi % (2 * np.pi))     # the floats themselves: a rounded key could hand back a neighbour's pulse
    if key not in _PULSES:
        from scipy.linalg import expm
        _PULSES[key] = expm(-1j * theta * (np.cos(phi) * IX + np.sin(phi) * IY))
    return _PULSES[key]


def turn(U, rho):
    return U @ rho @ U.conj().T


class Spin:
    """One pool, d rho/dt = L (rho - rho_eq), rho_eq = Iz (high temperature, on resonance). Without a Hamiltonian the
    generator of these jumps is real symmetric, exactly (asserted), so it is diagonalized once in an orthonormal basis;
    recover=False sets rho_eq = 0 during the evolution, which leaves the transverse pathway alone. Detection at the line:
    the integral of Tr(I+ rho(t)) over the acquisition (the Fourier transform at the resonance), in closed form on the
    eigenmodes; Iz carries no signal, so the acquisition needs no rho_eq."""

    def __init__(self, J, recover=True):
        L = superop([np.sqrt(J[abs(m)]) * TQUAD[m] for m in (-2, -1, 0, 1, 2) if J[abs(m)] > 0], 4)
        assert np.array_equal(L, L.conj().T) and not np.any(L.imag), "the generator is not real symmetric"
        self.w, self.V = np.linalg.eigh(L.real)
        self.Vi = self.V.T
        self.eq = Iz if recover else 0 * Iz
        self.obs = Ip.T.reshape(-1, order="F") @ self.V

    def evolve(self, rho, t):
        x = self.Vi @ (rho - self.eq).reshape(-1, order="F")
        return (self.V @ (np.exp(self.w * t) * x)).reshape(4, 4, order="F") + self.eq

    def detect(self, rho):
        a = self.obs * (self.Vi @ rho.reshape(-1, order="F"))
        live = np.abs(self.w) > 1e-12
        assert np.all(np.abs(a[~live]) < 1e-12), "signal on the trace mode"
        return np.sum(-a[live] / self.w[live])


def ir_scan(sp_, tau, phi1, phi2, s=1.0, invert=True):
    """NBM Fig. 1: 180_phi1 - tau - 90_phi2 - 90_0 - acquisition at receiver phase 0; s scales every flip angle."""
    rho = Iz.astype(complex)
    if invert:
        rho = sp_.evolve(turn(hard(np.pi * s, phi1), rho), tau)
    rho = turn(hard(np.pi / 2 * s, phi2), rho)
    return sp_.detect(turn(hard(np.pi / 2 * s, 0.0), rho))


def tq_scan(sp_, tau, phi, s=1.0, alternate=False, first=np.pi / 2, refocus=np.pi / 2):
    """TQTPPI with refocusing: 90_{phi+first} - tau/2 - 180_{phi+refocus} - tau/2 - 90_phi - 90_0 - acquisition. The
    papers do not give these phases; one pool's TQ/SQ does not depend on them while the single quantum survives
    (checked in Part 8 at 2 refocus - first = 90, 30 and -60 degrees). alternate=True adds the scan with the first pulse turned
    by 180 and the receiver inverted, which cancels everything that does not start from the first pulse's coherence: the
    magnetization recovering during the evolution, which the second 90 pulse turns into the same orders with the same
    TPPI labels, and what an imperfect first pulse leaves along z."""
    def one(f1):
        rho = sp_.evolve(turn(hard(np.pi / 2 * s, phi + f1), Iz.astype(complex)), tau / 2)
        rho = sp_.evolve(turn(hard(np.pi * s, phi + refocus), rho), tau / 2)
        rho = turn(hard(np.pi / 2 * s, phi), rho)
        return sp_.detect(turn(hard(np.pi / 2 * s, 0.0), rho))
    return (one(first) - one(first + np.pi)) / 2 if alternate else one(first)


def orders(vals):
    """The values over the eight phase steps, split into the parts turning as e^(-i p phi), p = -3..3."""
    return {p: np.mean(vals * np.exp(1j * p * TPPI)) for p in range(-3, 4)}


def ir_vals(sp_, tau, cycle="PC0", s=1.0, invert=True):
    k = np.deg2rad(135)
    if cycle == "PC0":
        vals = [ir_scan(sp_, tau, 0.0, f, s, invert) for f in TPPI]
    elif cycle == "PC1":                          # the inversion pulse at 0 and 180
        vals = [(ir_scan(sp_, tau, 0.0, f, s, invert) + ir_scan(sp_, tau, np.pi, f, s, invert)) / 2 for f in TPPI]
    else:                                         # PC2: the first 90 pulse at +-135 about its TPPI phase
        vals = [(ir_scan(sp_, tau, 0.0, f + k, s, invert) + ir_scan(sp_, tau, 0.0, f - k, s, invert)) / 2 for f in TPPI]
    return np.array(vals)


def ir_orders(sp_, tau, cycle="PC0", s=1.0, invert=True):
    return orders(ir_vals(sp_, tau, cycle, s, invert))


def tq_vals(sp_, tau, s=1.0, alternate=False, first=np.pi / 2, refocus=np.pi / 2):
    return np.array([tq_scan(sp_, tau, f, s, alternate, first, refocus) for f in TPPI])


def amp(o, p):
    return abs(o[p]) + abs(o[-p])


def project(o, ref, p):
    """Signed amplitude at order p: the pair (o[p], o[-p]) projected onto a reference pair."""
    return ((o[p] * np.conj(ref[p]) + o[-p] * np.conj(ref[-p])) / (abs(ref[p]) ** 2 + abs(ref[-p]) ** 2)).real


def free(olist, u, p):
    """The amplitude a fit with a free phase reads at order p: each of o[p], o[-p] over the evolution times projected onto
    the real shape u with a complex amplitude, magnitudes added."""
    return sum(abs(sum(o[q] * x for o, x in zip(olist, u))) for q in (p, -p)) / (u @ u)


def d_line(T2s, T2f):
    """d read at the line: the integral of f13 over the integral of f11, with the times the acquisition decays with."""
    return (np.sqrt(6) / 5) * (T2f - T2s) / (0.6 * T2f + 0.4 * T2s)


def pool_J(T1s, T1f, T2f):
    return 1 / (3 * T2f) - 1 / (6 * T1f), 1 / (6 * T1f), 1 / (6 * T1s)     # J0, J1, J2 for R1s, R1f, R2f


SHIFT_T = np.sqrt(6) / 80     # the transverse twin: the slow weight of the TQTPPI single quantum, 0.4 + sqrt6/80 |d|
MRM_TQTPPI_MAX = {2: 12.46, 4: 19.84, 6: 23.93}   # MRM Table 2: maximum TQ/SQ in %, TQTPPI with the refocusing pulse
def fig5_first():
    """NBM Fig. 5 (2 %, PC2) with no coarse points added, T1s and T1f in ms and the T1-TQ/SQ ratio in %, as the thesis's
    Fig. 4.11 draws it (data/reichert2024_sodium_agarose, read in Part 9)."""
    names, fit = read_csv("thesis_fig4_11_fit_vs_points_added_2pct_pc2.csv")
    col = {k: fit[:, i] for i, k in enumerate(names)}
    assert col["points_added"][0] == 0
    return col["T1s_ms"][0], col["T1f_ms"][0], col["tq_sq_percent"][0]


NBM_FIG5_FIRST = (47.5, 35.8, 8.6)   # NBM Fig. 5 (2 %, PC2), no coarse points added: T1s, T1f in ms, T1-TQ/SQ in %, read off


def d_variants(c):
    """|d| read at the line with the times the acquisition may decay with: the homogeneous T2 (no broadening, the largest
    |d|), the second paper's T2* put into the formula, and one Lorentzian broadening shared by both components, matched
    to the slow or to the fast T2* (the two components of its Table 1 do not share one at 2 and 6 %)."""
    T2s, T2f = MRM_T2[c]["T2s"][0], MRM_T2[c]["T2f"][0]
    Ss, Sf = MRM_T2STAR[c]
    rs, rf = 1 / Ss - 1 / T2s, 1 / Sf - 1 / T2f
    return {"T2": abs(d_line(T2s, T2f)), "T2* formula": abs(d_line(Ss, Sf)),
            "shared, slow": abs(d_line(Ss, 1 / (1 / T2f + rs))), "shared, fast": abs(d_line(1 / (1 / T2s + rf), Sf))}


def peak(Ts, Tf):
    """The peak of e^(-t/Tf) - e^(-t/Ts) in size, Ts > Tf."""
    t = np.log(Ts / Tf) * Ts * Tf / (Ts - Tf)
    return np.exp(-t / Ts) - np.exp(-t / Tf)


def part8(edges):
    print("\nPart 8. The two sequences pulse by pulse for one pool (hard pulses, relaxation between them, the eight TPPI "
          "steps, detection at the line), against the hand count with Wigner coefficients")
    import sympy as sp
    from sympy.physics.quantum.spin import Rotation
    from scipy.optimize import brentq
    from scipy.stats import chi2 as chi2_dist

    def wd(j, m, mp):
        return sp.Abs(sp.simplify(Rotation.d(j, m, mp, sp.pi / 2).doit()))
    hand = {"long": sp.Rational(4, 5) * wd(3, 3, 0) * wd(3, -1, 3) / (wd(1, 1, 0) * wd(1, -1, 1)),
            "trans": sp.sqrt(6) / 5 * wd(3, 3, 1) * wd(3, -1, 3) / (wd(1, 1, 1) * wd(1, -1, 1)),
            "shift": sp.Rational(2, 5) * wd(3, 1, 0) * wd(3, -1, 1) / (wd(1, 1, 0) * wd(1, -1, 1)),
            "twin": sp.sqrt(6) / 5 * wd(3, 1, 1) * wd(3, -1, 1) / (wd(1, 1, 1) * wd(1, -1, 1))}
    exact = {"long": sp.sqrt(6) / 4, "trans": 3 * sp.sqrt(6) / 16, "shift": sp.sqrt(6) / 40, "twin": sp.sqrt(6) / 80}
    ok = all(sp.simplify(hand[k] - exact[k]) == 0 for k in hand)
    print(f"  (a) sympy's Wigner elements at 90 degrees give exactly: T1-TQ/SQ = sqrt6/4 |d|, T2-TQ/SQ = 3 sqrt6/16 |d|, "
          f"ratio 4/3; the rank-3 order at the SQ frequency moves the slow weight to 0.8 + sqrt6/40 |d| (inversion recovery) "
          f"and 0.4 + sqrt6/80 |d| (the transverse sequence): {ok}")
    assert ok
    # (b) the two routes. The simulation reads an order out of eight scans whose values are larger than it, and evolves
    # through rounded eigenvalues, so its error is of order eps x (the dynamic range of the scans, largest value over the
    # order's size, + the largest rate times the evolution time); the gate is that law
    route, law, sim_w = 0.0, 0.0, {}
    thetas = ((np.pi / 2, np.pi / 2), (0.0, np.pi / 12), (0.0, -np.pi / 6))    # 2 refocus - first = 90, 30, -60 degrees
    for c in (2, 4, 6):
        T1s, _, T1f, _ = NBM_T1["PC2"][c]
        T2f = MRM_T2[c]["T2f"][0]
        J = pool_J(T1s, T1f, T2f)
        R1s, R1f, R2f, R2s = 6 * J[2], 6 * J[1], 3 * (J[0] + J[1]), 3 * (J[1] + J[2])
        assert max(abs(R1s * T1s - 1), abs(R1f * T1f - 1), abs(R2f * T2f - 1)) <= 8 * EPS, "pool_J"
        d = d_line(1 / R2s, 1 / R2f)
        one = Spin(J)
        wmax = np.max(np.abs(one.w))
        taus = np.array([5.0, 20.0, 60.0, 150.0])
        eq = ir_orders(one, 0.0, invert=False)
        val, sq, dyn = [], [], 0.0
        for t in taus:
            v = ir_vals(one, t)
            o = orders(v)
            val.append(amp(o, 3) / (np.exp(-R1s * t) - np.exp(-R1f * t)) / amp(eq, 1))
            dev = abs(val[-1] / (LONG * abs(d)) - 1)
            route, law = max(route, dev), max(law, dev / (EPS * (np.max(np.abs(v)) / amp(o, 3) + wmax * t)))
            sq.append(project(o, eq, 1))
            dyn = max(dyn, np.max(np.abs(v)) / amp(eq, 1))
        M = np.column_stack([np.exp(-R1s * taus), np.exp(-R1f * taus)])
        ws, wf = np.linalg.lstsq(M, (1 - np.array(sq)) / 2, rcond=None)[0]
        sim_w[c] = (abs(d), ws)
        dev = max(abs(ws - 0.8 - SHIFT * abs(d)), abs(wf - 0.2 + SHIFT * abs(d)))
        # the weights come out of a 4 x 2 least squares on exponentials, which amplifies by its condition number
        route, law = max(route, dev), max(law, dev / (EPS * (dyn + wmax * taus[-1]) * np.linalg.cond(M)))
        tr = Spin(J, recover=False)
        s2 = []
        t2 = np.array([0.5, 5.0, 15.0, 40.0])
        for f1, rf in thetas:
            o0 = orders(tq_vals(tr, 0.0, first=f1, refocus=rf))
            sqt, dyn_t = [], 0.0
            for t in t2:
                v = tq_vals(tr, t, first=f1, refocus=rf)
                o = orders(v)
                s2.append(amp(o, 3) / abs(np.exp(-R2f * t) - np.exp(-R2s * t)) / amp(o0, 1))
                dev = abs(s2[-1] / (TRANS * abs(d)) - 1)
                route, law = max(route, dev), max(law, dev / (EPS * (np.max(np.abs(v)) / amp(o, 3) + wmax * t)))
                sqt.append(project(o, o0, 1))
                dyn_t = max(dyn_t, np.max(np.abs(v)) / amp(o0, 1))
            Mt = np.column_stack([np.exp(-R2s * t2), np.exp(-R2f * t2)])
            wts, wtf = np.linalg.lstsq(Mt, np.array(sqt), rcond=None)[0]
            dev = max(abs(wts - 0.4 - SHIFT_T * abs(d)), abs(wtf - 0.6 + SHIFT_T * abs(d)))
            route, law = max(route, dev), max(law, dev / (EPS * (dyn_t + wmax * t2[-1]) * np.linalg.cond(Mt)))
        print(f"  {c}% pool (published T1s, T1f; T2f measured; T2s by the relation; d at the line {d:.4f}): T1-TQ/SQ "
              f"{val[0]:.5f} (hand {LONG * abs(d):.5f}); T2-TQ/SQ {s2[0]:.5f} (hand {TRANS * abs(d):.5f}); ratio "
              f"{val[0] / s2[0]:.5f}; slow weight at the SQ frequency {ws:.5f} (hand {0.8 + SHIFT * abs(d):.5f}), in the "
              f"transverse sequence {wts:.5f} (hand {0.4 + SHIFT_T * abs(d):.5f})")
    print(f"  (b) the two routes over three pools, four evolution times each, the weights of both sequences, and three choices "
          f"of the transverse phases (2 refocus - first = 90, 30, -60 degrees): largest relative deviation {route:.1e}, at "
          f"most {law:.0f} x eps x (the dynamic range of the scans it is read from + the largest rate times the evolution time)")
    assert law < 100, law
    # the phase cycles, order by order: the scan is a trigonometric polynomial of degree 3 in a pulse's phase, so adding
    # the two scans at +-alpha multiplies order p by cos(alpha p) for any flip angle; PC1 changes nothing for an ideal pi
    T1s, _, T1f, _ = NBM_T1["PC2"][2]
    J2 = pool_J(T1s, T1f, MRM_T2[2]["T2f"][0])
    one2 = Spin(J2)
    rel, scale = 0.0, None
    for s in (0.90, 1.00, 1.10):
        o0, o2 = ir_orders(one2, 20.0, "PC0", s), ir_orders(one2, 20.0, "PC2", s)
        top = max(abs(o0[p]) for p in range(-3, 4))
        rel = max(rel, max(abs(o2[p] - np.cos(np.deg2rad(135) * p) * o0[p]) for p in range(-3, 4)) / top)
        if s == 0.90:
            scale = (abs(o2[1]) / abs(o0[1]), abs(o2[3]) / abs(o0[3]), abs(o2[2]) / top, abs(o0[2]) / top)
    o0, o1 = ir_orders(one2, 20.0, "PC0"), ir_orders(one2, 20.0, "PC1")
    rel1 = max(abs(o1[p] - o0[p]) for p in range(-3, 4)) / max(abs(o0[p]) for p in range(-3, 4))
    print(f"      PC2 against cos(135 p) x PC0, order by order, at every flip angle 0.9, 1.0 and 1.1 times: largest deviation "
          f"{rel:.1e} of the largest order; at 0.9 the single and the triple quantum keep {scale[0]:.6f} and {scale[1]:.6f}, "
          f"the double quantum {scale[2]:.1e} of the largest order (PC0: {scale[3]:.4f}); PC1 against PC0 with ideal pulses: "
          f"{rel1:.1e}")
    assert rel < 100 * EPS and rel1 < 100 * EPS, (rel, rel1)

    # (c) the full dynamics (recovery included) on each pool, the TQ amplitudes read as a fit with a free phase reads them,
    # the first transverse pulse alternated with the receiver, which cancels what does not start from that pulse
    print("  (c) full dynamics (recovery included), the first transverse pulse alternated with the receiver, every flip angle "
          "scaled 0.90, 0.95, 1.05, 1.10 times, three choices of the transverse phases (2 refocus - first = 90, 30, -60):")
    allr = {5: [], 10: []}
    grid = T2_GRIDS["linear"]
    for c in (2, 4, 6):
        T1s, _, T1f, _ = NBM_T1["PC2"][c]
        T2f = MRM_T2[c]["T2f"][0]
        Jc = pool_J(T1s, T1f, T2f)
        onec = Spin(Jc)
        R1s, R1f, R2f, R2s = 6 * Jc[2], 6 * Jc[1], 3 * (Jc[0] + Jc[1]), 3 * (Jc[1] + Jc[2])
        g1, g2 = T1_GRID[::10], np.arange(2.0, 122.0, 2.0)
        u1, u2 = np.exp(-R1s * g1) - np.exp(-R1f * g1), np.exp(-R2s * g2) - np.exp(-R2f * g2)

        def long_tq(s, cyc):
            return free([ir_orders(onec, t, cyc, s) for t in g1], u1, 3) / amp(ir_orders(onec, 0.0, cyc, s, invert=False), 1)

        def trans_tq(s, f1, rf):
            sq0 = amp(orders(tq_vals(onec, 0.0, s, True, f1, rf)), 1)
            return free([orders(tq_vals(onec, t, s, True, f1, rf)) for t in g2], u2, 3) / sq0
        lt = {s: long_tq(s, "PC0") for s in (0.90, 0.95, 1.00, 1.05, 1.10)}
        pc1 = [long_tq(s, "PC1") / lt[s] for s in (0.90, 0.95, 1.05, 1.10)]
        # the same ideal values as (b), read through the free phase: far above the rounding (b) gates, a check of the reading
        assert abs(lt[1.00] / (LONG * abs(d_line(1 / R2s, 1 / R2f))) - 1) < 1e-10, lt[1.00]
        r = {5: [], 10: []}
        dT2, slows = {5: [], 10: []}, {5: [], 10: []}
        for f1, rf in thetas:
            assert abs(lt[1.00] / trans_tq(1.0, f1, rf) - 4 / 3) < 1e-9
            fits = {}
            for s in (0.90, 0.95, 1.00, 1.05, 1.10):
                ref = orders(tq_vals(onec, 0.0, s, True, f1, rf))
                sq = np.array([project(orders(tq_vals(onec, t, s, True, f1, rf)), ref, 1) for t in grid])
                fits[s] = fit_mrm(sq, grid, (1 / R2s, T2f))
            for k, ss in ((5, (0.95, 1.05)), (10, (0.90, 1.10))):
                r[k] += [lt[s] / trans_tq(s, f1, rf) for s in ss]
                dT2[k] += [fits[s][0] - fits[1.00][0] for s in ss]
                slows[k] += [fits[s][2] for s in ss]
        allr[5] += r[5]
        allr[10] += r[10]
        gap_ms = MRM_T2[c]["T2s"][0] - 1 / R2s
        s_meas, s_err = data_numbers(c)["slow"]
        a60, r60 = orders(tq_vals(onec, 60.0, alternate=True)), orders(tq_vals(onec, 60.0))
        a0 = orders(tq_vals(onec, 0.0, alternate=True))
        lk3 = (abs(r60[3] - a60[3]) + abs(r60[-3] - a60[-3])) / amp(a60, 3)
        lk1 = abs(r60[1] - a60[1]) + abs(r60[-1] - a60[-1])
        print(f"    {c}%: ratio {min(r[5]):.2f} to {max(r[5]):.2f} at 5 %, {min(r[10]):.2f} to {max(r[10]):.2f} at 10 %; PC1/PC0 "
              + ", ".join(f"{x:.3f}" for x in pc1) + f"; the transverse SQ fitted with Eq. 11 on the 26 times: T2s longer by "
              f"{min(dT2[5]):+.2f} to {max(dT2[5]):+.2f} ms at 5 %, {min(dT2[10]):+.2f} to {max(dT2[10]):+.2f} ms at 10 % "
              f"({max(dT2[10]) / gap_ms:.2f} of the gap at most), slow fraction {min(slows[5]):.3f} to {max(slows[5]):.3f} at 5 % "
              f"({(min(slows[5]) - s_meas) / s_err:+.1f} to {(max(slows[5]) - s_meas) / s_err:+.1f} sigma) and "
              f"{min(slows[10]):.3f} to {max(slows[10]):.3f} at 10 %; without the alternation, ideal pulses, the recovering "
              f"magnetization at 60 ms is {lk3:.2f} of the transverse TQ and {lk1 / amp(a0, 1):.2f} of the SQ at the start, "
              f"{lk1 / amp(a60, 1):.1f} times the transverse SQ there")
    print(f"      over the pools and phases: ratio {min(allr[5]):.2f} to {max(allr[5]):.2f} at 5 %, {min(allr[10]):.2f} to "
          f"{max(allr[10]):.2f} at 10 %")
    rlo, rhi = min(allr[10]), max(allr[10])

    # (d) d, what the acquisition makes of a rank-3 coherence against a rank-1 one, read at the line
    dv = {c: d_variants(c) for c in (2, 4, 6)}
    print("  (d) |d| read at the line, by the times the acquisition decays with: " + "; ".join(
        f"{c}%: " + ", ".join(f"{k} {v:.3f}" for k, v in dv[c].items()) for c in (2, 4, 6)))
    dmin = {c: min(dv[c].values()) for c in (2, 4, 6)}
    dmax = {c: max(dv[c].values()) for c in (2, 4, 6)}

    # (e) TQ1 through the ratio, two bound pools held at it, and TQ1 from the size alone
    print("  (e) TQ1, the longitudinal TQ against one pool's with the fitted times, = (3/4) x measured ratio (errors from "
          "the T1-TQ/SQ alone; Fig. 6C carries none); two bound pools allow up to "
          f"{edges[2]:.2f}, {edges[4]:.2f}, {edges[6]:.2f} (Part 6, a search from below); held at TQ1 itself and fitted, as in "
          "Part 5, to the five numbers of PC2:")
    for c in (2, 4, 6):
        out = []
        for pc in ("PC1", "PC2"):
            target = 0.75 * NBM_TQ1[pc][c][0] / NBM_TQ2_FIG6C[c]

            def res(q):
                ap, B1 = two_bound(q, c)
                return np.append(pulls(ap, c), 50 * (B1 - target))
            best = min((least_squares(res, s0, bounds=(TWO_LO, TWO_HI), x_scale="jac", max_nfev=300)
                        for s0 in guided_starts(c)), key=lambda r: r.cost)
            ap, B1 = two_bound(best.x, c)
            size_lo = NBM_TQ1[pc][c][0] / (100 * LONG * dmax[c])
            size_hi = NBM_TQ1[pc][c][0] / (100 * LONG * dmin[c])
            out.append(f"{pc} TQ1 {target:.4f}: chi2 {np.sum(pulls(ap, c) ** 2):.2f} at TQ1 {B1:.3f}; from the size alone "
                       f"{size_lo:.3f} to {size_hi:.3f} over the range of d")
        print(f"    {c}%: " + "; ".join(out))
    top = max(0.75 * NBM_TQ1[pc][c][0] / NBM_TQ2_FIG6C[c] for pc in ("PC1", "PC2") for c in (2, 4, 6))
    print(f"      the largest TQ1 is {top:.4f}: TQ1 = 1 would need one pool's ratio {1 / top:.2f} times smaller than 4/3, "
          f"against {rlo:.2f} to {rhi:.2f} over the flip angles 10 % off, the phases and the pools of (c)")

    # (f) reading (i), each phase cycle with its own times and T1m
    split_i = {}
    for pc, t1m in (("PC1", NBM_T1M_PC1), ("PC2", NBM_T1M_PC2)):
        chi, chi_x, rows = 0.0, 0.0, []
        for c in (2, 4, 6):
            T1s, _, T1f, _ = NBM_T1[pc][c]
            T1m, sT1m = t1m[c]
            T2s, sT2s = MRM_T2[c]["T2s"]
            TQ1 = 0.75 * NBM_TQ1[pc][c][0] / NBM_TQ2_FIG6C[c]
            sTQ1 = 0.75 * NBM_TQ1[pc][c][1] / NBM_TQ2_FIG6C[c]
            u = np.exp(-T1_GRID / T1s) - np.exp(-T1_GRID / T1f)

            def rates(delta):
                return 1 / T1m - 0.2 * delta, 1 / T1m + 0.8 * delta      # 0.8 R1s + 0.2 R1f = 1/T1m, R1f - R1s = delta

            def proj(delta):
                Rs, Rf = rates(delta)
                return (u @ (np.exp(-T1_GRID * Rs) - np.exp(-T1_GRID * Rf))) / (u @ u)
            dl = brentq(lambda x: proj(x) - TQ1, 1e-9, 0.05)
            split_i[(pc, c)] = dl
            dlo = brentq(lambda x: proj(x) - max(TQ1 - sTQ1, 1e-6), 1e-12, 0.05)
            dhi = brentq(lambda x: proj(x) - (TQ1 + sTQ1), 1e-9, 0.05)
            Rs, Rf = rates(dl)
            assert abs((0.8 * Rs + 0.2 * Rf) * T1m - 1) < 8 * EPS and abs(((Rs + Rf) / 2 - 1 / T1m) / (0.3 * dl) - 1) < 1e-9
            out = []
            for k in (0.3, 0.3 + SHIFT * dmax[c]):
                pred = 1 / (1 / T1m + k * dl)
                sp_pred = np.hypot(sT1m * (pred / T1m) ** 2,
                                   abs(1 / (1 / T1m + k * dhi) - 1 / (1 / T1m + k * dlo)) / 2)
                out.append((pred, sp_pred, (T2s - pred) / np.hypot(sT2s, sp_pred)))
            chi += out[0][2] ** 2
            chi_x += out[1][2] ** 2
            rows.append(f"    {c}%: TQ1 {TQ1:.4f} +- {sTQ1:.4f}; one pool with a true split of {1e3 * dl:.2f} 1/s (published "
                        f"{1e3 * (1 / T1f - 1 / T1s):.2f}) gives it through the published times' shape and predicts T2s "
                        f"{out[0][0]:.2f} +- {out[0][1]:.2f} ms, measured {T2s} +- {sT2s}: {out[0][2]:+.2f} sigma (with the "
                        f"pathway at the SQ frequency at the largest |d|: {out[1][2]:+.2f}); against NBM's own T2s (Fig. 6A, no "
                        f"error bars) {out[0][0] - NBM_FIG6A_T2S[c]:+.2f} ms; its TQ peaks at {np.log(Rf / Rs) / (Rf - Rs):.1f} ms, "
                        f"the published shape at {np.log(T1s / T1f) * T1s * T1f / (T1s - T1f):.1f} ms")
        print(f"  (f) reading (i), {pc} (its own T1s, T1f and tabulated T1m):")
        print("\n".join(rows))
        print(f"      chi2 {chi:.2f} for 3, p = {chi2_dist.sf(chi, 3):.3f} ({chi_x:.2f}, p = {chi2_dist.sf(chi_x, 3):.3f}, "
              f"with the pathway)")
    def chi_i(pc, t1m, scale):
        """reading (i)'s chi2 with TQ1 scaled by 4/3 over one pool's actual ratio"""
        tot = 0.0
        for c in (2, 4, 6):
            T1s, _, T1f, _ = NBM_T1[pc][c]
            T1m, sT1m = t1m[c]
            T2s, sT2s = MRM_T2[c]["T2s"]
            TQ1 = scale * 0.75 * NBM_TQ1[pc][c][0] / NBM_TQ2_FIG6C[c]
            sTQ1 = scale * 0.75 * NBM_TQ1[pc][c][1] / NBM_TQ2_FIG6C[c]
            u = np.exp(-T1_GRID / T1s) - np.exp(-T1_GRID / T1f)

            def proj(delta):
                return (u @ (np.exp(-T1_GRID * (1 / T1m - 0.2 * delta)) - np.exp(-T1_GRID * (1 / T1m + 0.8 * delta)))) / (u @ u)
            dl = brentq(lambda x: proj(x) - TQ1, 1e-9, 0.05)
            dlo = brentq(lambda x: proj(x) - max(TQ1 - sTQ1, 1e-6), 1e-12, 0.05)
            dhi = brentq(lambda x: proj(x) - (TQ1 + sTQ1), 1e-9, 0.05)
            pred = 1 / (1 / T1m + 0.3 * dl)
            sp_pred = np.hypot(sT1m * (pred / T1m) ** 2, abs(1 / (1 / T1m + 0.3 * dhi) - 1 / (1 / T1m + 0.3 * dlo)) / 2)
            tot += ((T2s - pred) / np.hypot(sT2s, sp_pred)) ** 2
        return tot
    for pc, t1m in (("PC1", NBM_T1M_PC1), ("PC2", NBM_T1M_PC2)):
        lo_chi, hi_chi = chi_i(pc, t1m, (4 / 3) / rhi), chi_i(pc, t1m, (4 / 3) / rlo)
        print(f"      {pc}: with one pool's ratio anywhere in {rlo:.2f} to {rhi:.2f} instead of 4/3, TQ1 scales by "
              f"{(4 / 3) / rhi:.3f} to {(4 / 3) / rlo:.3f} and reading (i)'s chi2 runs from {lo_chi:.2f} (p = "
              f"{chi2_dist.sf(lo_chi, 3):.3f}) to {hi_chi:.2f} (p = {chi2_dist.sf(hi_chi, 3):.3f})")
    T1s, _, T1f, _ = NBM_T1["PC2"][2]
    a, b, r0 = fig5_first()
    print(f"      NBM Fig. 5 (2 %, PC2; the thesis's Fig. 4.11): with no coarse points the fit reads T1s {a}, T1f {b} ms and "
          f"TQ/SQ {r0} %, with all of "
          f"them {T1s}, {T1f} and {NBM_TQ1['PC2'][2][0]} %; the coefficient times the peak of its own shape, the height of the "
          f"triple quantum against the single quantum, is {r0 * peak(a, b):.3f} and {NBM_TQ1['PC2'][2][0] * peak(T1s, T1f):.3f} %")

    # (g) the pathway Eq. 13 leaves out, for one pool with the published split; (h) its transverse twin; (i) the transverse
    # amplitude; (j) PC2/PC1
    print("  (g) one pool with the published split, its SQ with the rank-3 order at the SQ frequency (weights "
          "0.8 + sqrt6/40 |d|, gated in (b)), fitted with Eq. 13 (weights fixed at 0.8 / 0.2), at the smallest and the "
          "largest |d| of (d):")
    for c in (2, 4, 6):
        T1s, _, T1f, _ = NBM_T1["PC2"][c]
        shape = np.exp(-T1_GRID / T1f) - np.exp(-T1_GRID / T1s)
        dd0, ws0 = sim_w[c]
        assert abs(0.8 + SHIFT * dd0 - ws0) < 1e-12, "the closed form against the simulated weight of (b)"
        cols = []
        for dd in (dmin[c], dmax[c]):
            sq = 1 - 2 * ((0.8 + SHIFT * dd) * np.exp(-T1_GRID / T1s) + (0.2 - SHIFT * dd) * np.exp(-T1_GRID / T1f))
            res = []
            for wlab, w in (("one pool's own TQ", LONG * dd), ("the reported TQ", NBM_TQ1_PC2[c])):
                _, Ts, Tf, _ = fit_eq13(sq, w * shape, (T1s, T1f))
                res.append(f"{wlab}: split {100 * ((1 / Tf - 1 / Ts) / (1 / T1f - 1 / T1s) - 1):+.1f} %, T1m "
                           f"{1 / (0.8 / Ts + 0.2 / Tf) - 1 / (0.8 / T1s + 0.2 / T1f):+.2f} ms")
            cols.append(f"[|d| {dd:.3f}, weights moved by {SHIFT * dd:.3f}; " + "; ".join(res) + "]")
        print(f"    {c}%: " + " ".join(cols))
    rows = []
    for c in (2, 4, 6):
        s_meas, s_err = data_numbers(c)["slow"]
        lo, hi = 0.4 + SHIFT_T * dmin[c], 0.4 + SHIFT_T * dmax[c]
        rows.append(f"{c}%: {lo:.3f} to {hi:.3f} against the measured {s_meas:.3f} +- {s_err:.3f} "
                    f"({(lo - s_meas) / s_err:+.1f} to {(hi - s_meas) / s_err:+.1f} sigma; without the twin "
                    f"{(0.4 - s_meas) / s_err:+.1f})")
    print("  (h) the transverse twin, which Parts 3 to 5 do not carry: one pool's slow fraction of the transverse single quantum, "
          "0.4 + sqrt6/80 |d|, over the range of d: " + "; ".join(rows))
    rows = []
    for c in (2, 4, 6):
        th_lo, th_hi = TRANS * dmin[c], TRANS * dmax[c]
        rows.append(f"{c}%: {100 * th_lo:.1f} to {100 * th_hi:.1f} %, the measured {NBM_TQ2_FIG6C[c]} % is "
                    f"{NBM_TQ2_FIG6C[c] / (100 * th_hi):.2f} to {NBM_TQ2_FIG6C[c] / (100 * th_lo):.2f} of it")
    print("  (i) the transverse TQ/SQ of one pool read at the line, over the range of d: " + "; ".join(rows))
    rows = []
    for c in (2, 4, 6):
        T2s, T2f = MRM_T2[c]["T2s"][0], MRM_T2[c]["T2f"][0]
        amp2 = MRM_TQTPPI_MAX[c] / peak(T2s, T2f)
        rows.append(f"{c}%: {MRM_TQTPPI_MAX[c]} % / {peak(T2s, T2f):.3f} = {amp2:.1f} % against Fig. 6C's {NBM_TQ2_FIG6C[c]} % "
                    f"({amp2 / NBM_TQ2_FIG6C[c]:.2f}); over the range of d {amp2 / (100 * TRANS * dmax[c]):.2f} to "
                    f"{amp2 / (100 * TRANS * dmin[c]):.2f} of one pool's")
    print("      a second reading of the measured transverse amplitude: the second paper's transverse maximum (its Table 2, "
          "TQTPPI with the refocusing pulse) over one pool's shape peak at its T2: " + "; ".join(rows))
    print("  (j) PC2 over PC1, measured: " + ", ".join(f"{NBM_TQ1['PC2'][c][0] / NBM_TQ1['PC1'][c][0]:.2f}" for c in (2, 4, 6))
          + " (in the model PC2 equals PC0 for any flip angle, (b); a flip-angle error raises PC1, (c))")
    rows = []
    for dlt in (-5.0, -10.0):
        k = np.deg2rad(135 + dlt)
        rows.append(f"{dlt:+.0f} degrees: PC2's TQ/SQ x {abs(np.cos(3 * k)) / abs(np.cos(k)):.2f}, double quantum left "
                    f"{abs(np.cos(2 * k)):.2f} of PC0's")
    print("      a phase error in PC2's +-135 degree step multiplies order p by cos(p (135 + error)): " + "; ".join(rows))

    # (k) the first order in time. One pool has as many unknowns there (J0, J1, J2, d) as the data have numbers: the two
    # initial rates, 1/T1m and 0.6/T2f + 0.4/T2s, and the two triple quanta's initial slopes, each its amplitude times its
    # split. So the first order measures d and tests the pool only through inequalities. J in the page's convention
    # (R1f = 2 J1, R1s = 2 J2, R2f = J0 + J1, R2s = J1 + J2), three times pool_J's.
    J0, J1, J2, dd = sp.symbols("J0 J1 J2 d")
    E = sp.symbols("E1:5", positive=True)
    first = [sp.Rational(2, 5) * J1 + sp.Rational(8, 5) * J2, sp.Rational(3, 5) * J0 + J1 + sp.Rational(2, 5) * J2,
             exact["long"] * dd * 2 * (J1 - J2), exact["trans"] * dd * (J0 - J2)]
    comb = sp.expand(sp.Rational(6, 5) * first[2] + sp.Rational(16, 5) * first[3] - sp.sqrt(6) * dd * (first[1] - first[0]))
    sols = sp.solve([f - e for f, e in zip(first, E)], [J0, J1, J2, dd], dict=True)
    ok = comb == 0 and len(sols) == 1 and all(sp.simplify(f.subs(sols[0]) - e) == 0 for f, e in zip(first, E))
    print("  (k) the first order in time, one pool: 1/T1m = 0.4 J1 + 1.6 J2, 0.6/T2f + 0.4/T2s = 0.6 J0 + J1 + 0.4 J2, the "
          "triple quanta's initial slopes (amplitude times split) sqrt6/2 |d| (J1 - J2) and 3 sqrt6/16 |d| (J0 - J2), J in "
          "the page's convention; 1.2 and 3.2 times the slopes add to sqrt6 |d| times the difference of the two initial "
          f"rates, and the four numbers have exactly one solution in J0, J1, J2 and d: {ok}")
    assert ok
    # the second route for the rates: the generator's own slopes at t = 0 for one pool, against the formulas in J (the
    # page's J are three times the generator's); each slope is a weighted sum of G x, rounded to at most 64 eps times
    # the sum of the magnitudes of its terms
    checks = []
    rng = np.random.default_rng(7)
    for Jp in [(0.031, 0.0047, 0.0039)] + [tuple(10 ** rng.uniform(-6, 2, 3)) for _ in range(200)]:
        Lg = voxel(1.0, Jp, Jp, 0.0)
        Jg = [3 * x for x in Jp]
        for ids, x1, (w_sq, n_sq), (w_tq, n_tq), r_sq, s_tq in (
                (ORDER0, T10, (T10, T10 @ T10), (T30, 4.0), 0.4 * Jg[1] + 1.6 * Jg[2], -2 * (Jg[1] - Jg[2])),
                (ORDER1, T11, (T11, T11 @ T11), (T31, 2 * np.sqrt(3)), 0.6 * Jg[0] + Jg[1] + 0.4 * Jg[2], -(Jg[0] - Jg[2]))):
            G = Lg[np.ix_(ids, ids)].real
            x = np.concatenate([x1, 0 * x1])
            n = len(x1)
            gx, ga = G @ x, np.abs(G) @ np.abs(x)
            for w, norm, want, sign in ((w_sq, n_sq, r_sq, -1), (w_tq, n_tq, s_tq, 1)):
                got = sign * (w @ gx[:n] + w @ gx[n:]) / norm
                bound = 64 * EPS * (np.abs(w) @ ga[:n] + np.abs(w) @ ga[n:]) / abs(norm) + 4 * EPS * abs(want)
                checks.append(abs(got - want) <= bound)
    print("      second route: the generator's slopes at t = 0 for one pool give the two initial rates and the two "
          f"triple-quantum rate differences of the formulas, each within its rounding bound, at one density set and 200 "
          f"drawn over eight decades: {sum(checks)} of {len(checks)}")
    assert all(checks)
    for pc, t1m in (("PC1", NBM_T1M_PC1), ("PC2", NBM_T1M_PC2)):
        rows = []
        for c in (2, 4, 6):
            T1s, _, T1f, _ = NBM_T1[pc][c]
            T2s, T2f = MRM_T2[c]["T2s"][0], MRM_T2[c]["T2f"][0]
            num = [1 / t1m[c][0], 0.6 / T2f + 0.4 / T2s, NBM_TQ1[pc][c][0] / 100 * (1 / T1f - 1 / T1s),
                   NBM_TQ2_FIG6C[c] / 100 * (1 / T2f - 1 / T2s)]
            sol = {k: float(v.subs(dict(zip(E, num)))) for k, v in sols[0].items()}
            j0, j1, j2, dv1 = 1e3 * sol[J0], 1e3 * sol[J1], 1e3 * sol[J2], abs(sol[dd])
            split = 2 * (j1 - j2)
            rows.append(f"{c}%: |d| {dv1:.3f} (largest {dmax[c]:.3f}); J0, J1, J2 = {j0:.1f}, {j1:.2f}, {j2:.2f} 1/s "
                        f"(J1/J2 {j1 / j2:.3f}, J0/J1 {j0 / j1:.2f}), "
                        f"monotone {j0 >= j1 >= j2}, 3 J0 + J2 - 4 J1 = {3 * j0 + j2 - 4 * j1:.1f} and J2 - J1/4 = "
                        f"{j2 - j1 / 4:.2f} (both >= 0 for any sum of Lorentzians), the exact region of such sums "
                        f"J0 >= J1 >= J2 >= J0 J1/(4 J0 - 3 J1), margin above its lower edge J2/J0 - (J1/J0)/(4 - 3 J1/J0) = "
                        f"{j2 / j0 - (j1 / j0) / (4 - 3 * j1 / j0):.3f}; split 2 (J1 - J2) {split:.2f} 1/s against "
                        f"reading (i)'s {1e3 * split_i[(pc, c)]:.2f} (x {split / (1e3 * split_i[(pc, c)]):.3f})")
        print(f"    {pc}: " + "; ".join(rows))
    T1s, _, T1f, _ = NBM_T1["PC2"][2]
    a, b, r0 = fig5_first()
    s0, s1 = r0 / 100 * 1e3 * (1 / b - 1 / a), NBM_TQ1["PC2"][2][0] / 100 * 1e3 * (1 / T1f - 1 / T1s)
    h0, h1 = r0 * peak(a, b), NBM_TQ1["PC2"][2][0] * peak(T1s, T1f)
    print(f"      NBM Fig. 5 (2 %, PC2; the thesis's Fig. 4.11): the longitudinal slope, amplitude times split, goes from {s0:.3f} to {s1:.3f} 1/s "
          f"({100 * (s1 / s0 - 1):+.1f} %) as the late points are added, the height from {h0:.3f} to {h1:.3f} % "
          f"({100 * (h1 / h0 - 1):+.1f} %)")
    return split_i


# ---------------- Part 9: the first paper's sampling of T1, read off the thesis's Fig. 3.6a and the paper's Fig. 2A ----------------
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "reichert2024_sodium_agarose")
NBM_TEXT_SAMPLING = (0.1, 40, 8, 300)   # NBM text: dense step 100 µs, 40 phase cycles of 8 steps, then 300 steps of 8 dense steps
NBM_FIG2A_PEAKS = ((72.6, 0.604), (87.1, 0.711), (101.5, 0.784), (116.0, 0.846), (130.4, 0.885))   # NBM Fig. 2A (2 %,
# PC0), a raster image: the first five coarse peaks, read off the trace's colour at the figure's resolution, calibrated on
# its gridlines
NBM_FIG2A_DENSE_END = (63.6, 64.5)   # the last column of the dense band that is filled and the first one that is not, ms
NBM_FIG2A_FIRST, NBM_FIG2A_LAST, NBM_FIG2A_PERIODS = 72.6, 635.8, 39   # the first and the last coarse peak, ms, and the
# periods between them, each peak read on its own
THESIS_FIG411_CAL = {"ms": 1.30e-4, "percent": 2.23e-5}   # the largest residuals of Fig. 4.11's two value axes on their tick
# marks (data/reichert2024_sodium_agarose/reader_output.txt)


def read_csv(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        rows = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    return rows[0].split(","), np.array([[float(v) for v in r.split(",")] for r in rows[1:]])


def part9():
    print("\nPart 9. The first paper's sampling of T1, read off the thesis's Fig. 3.6a and the paper's Fig. 2A (2 %, PC0)")
    _, fid = read_csv("thesis_fig3_6a_irtqtppi_fid_2pct_pc0.csv")
    t, y = fid[:, 0], fid[:, 1]
    big = np.diff(t) > 1.0                       # steps of more than a millisecond
    jump = np.flatnonzero(np.diff(big.astype(int)))
    assert len(jump) == 1 and not big[0] and big[-1], "one change of step, from short to long"
    bands = (np.arange(jump[0] + 2), np.arange(jump[0] + 2, len(t)))   # the first long step follows index jump + 1
    assert len(T1_GRID) == len(t), "T1_GRID has the path's number of points"
    print(f"  the thesis's Fig. 3.6a, the same FID drawn as a vector path (data/reichert2024_sodium_agarose): {len(t)} "
          f"vertices, the step changing once")
    dense, coarse = T1_GRID[:len(bands[0])], T1_GRID[len(bands[0]):]
    for name, k, grid in zip(("dense", "coarse"), bands, (dense, coarse)):
        slope, icept = np.polyfit(k, t[k], 1)
        line_res = t[k] - (slope * k + icept)
        scatter = np.abs(line_res).max()        # the band's vertices around their own evenly stepped line
        n = len(k)
        offset = np.mean(t[k] - grid)
        # each end is off the band's line by at most the scatter, so the mean step from the two ends is known to
        # 2 scatter / (n - 1), and the band's mean place against the grid to the scatter
        assert len(grid) == n and abs((t[k][-1] - t[k][0]) - (grid[-1] - grid[0])) / (n - 1) <= 2 * scatter / (n - 1), \
            f"T1_GRID's {name} step is the path's"
        assert abs(offset) <= scatter, f"T1_GRID's {name} band sits where the path's does"
        mx = [i for i in k[1:-1] if y[i] > y[i - 1] and y[i] >= y[i + 1] and y[i] > 0.1]
        print(f"      {name}: {n} samples from {t[k][0]:.4f} to {t[k][-1]:.4f} ms, mean step "
              f"{(t[k][-1] - t[k][0]) / (n - 1):.6f} ms against T1_GRID's {(grid[-1] - grid[0]) / (n - 1):g}; around one "
              f"evenly stepped line within {scatter:.4f} ms (rms {np.sqrt(np.mean(line_res ** 2)):.4f}), off T1_GRID by "
              f"{offset:+.4f} ms on average; maxima above 0.1 spaced "
              f"{dict(sorted(Counter(int(v) for v in np.diff(mx)).items()))} samples (spacing: count)")
    dev = np.abs(t - T1_GRID)
    print(f"      T1_GRID, the grid every refit on this page uses: {len(dense)} points {dense[1] - dense[0]:g} ms apart to "
          f"{dense[-1]:g} ms and {len(coarse)} points {coarse[1] - coarse[0]:g} ms apart to {coarse[-1]:g} ms; every vertex "
          f"lies within {dev.max():.4f} ms of its point, each band's mean step and mean place within what that band's "
          f"scatter allows")
    co = bands[1]
    maxima = [i for i in co[1:-1] if y[i] > y[i - 1] and y[i] >= y[i + 1]]
    tq, hq = t[maxima], y[maxima]
    tr = np.array([p[0] for p in NBM_FIG2A_PEAKS])
    hr = np.array([p[1] for p in NBM_FIG2A_PEAKS])
    step = (NBM_FIG2A_LAST - NBM_FIG2A_FIRST) / (8 * NBM_FIG2A_PERIODS)   # a TPPI period is eight phase steps
    print(f"  the paper's Fig. 2A, read off its raster as a second route: the dense band ends between {NBM_FIG2A_DENSE_END[0]} "
          f"and {NBM_FIG2A_DENSE_END[1]} ms; its coarse peaks, {NBM_FIG2A_PERIODS} periods of eight phase steps from "
          f"{NBM_FIG2A_FIRST} to {NBM_FIG2A_LAST} ms on the raster's own gridline calibration, give {step:.3f} ms per coarse "
          f"point; against the path's {len(tq)} coarse maxima the first five lie {np.min(tq[:5] - tr):.2f} to "
          f"{np.max(tq[:5] - tr):.2f} ms early and within {np.max(np.abs(hq[:5] - hr)):.4f} of their heights, the last "
          f"{NBM_FIG2A_LAST - tq[-1]:.2f} ms late: the raster's calibration stretched by "
          f"{100 * (step / (coarse[1] - coarse[0]) - 1):+.2f} % against the path")
    h, n_cyc, n_step, n_coarse = NBM_TEXT_SAMPLING
    print(f"  the text's step of {h * 1e3:g} µs, {n_cyc} cycles of {n_step} steps and {n_coarse} coarse points of {n_step} "
          f"steps would end the two bands at {h * n_cyc * n_step:g} and {h * n_cyc * n_step + n_coarse * n_step * h:g} ms")
    T1s, _, T1f, _ = NBM_T1["PC0"][2]

    def recovery(tt):
        return 1 - 1.6 * np.exp(-tt / T1s) - 0.4 * np.exp(-tt / T1f)
    t_text = dense[-1] + (tq[:5] - dense[-1]) * n_step * h / (coarse[1] - coarse[0])
    print("  the path's first five coarse maxima against one pool with the published times (PC0, inversion complete) on "
          f"the drawn axis, and if the coarse step were the text's {n_step * h:g} ms: " + "; ".join(
              f"{tt:.2f} ms: {hh:.3f} against {recovery(tt):.3f} and {recovery(t8):.3f}"
              for tt, hh, t8 in zip(tq[:5], hq[:5], t_text)))
    print(f"      largest deviation {np.max(np.abs(hq[:5] - recovery(tq[:5]))):.3f} on the drawn axis, "
          f"{np.max(np.abs(hq[:5] - recovery(t_text))):.3f} with the text's step")
    names, fit = read_csv("thesis_fig4_11_fit_vs_points_added_2pct_pc2.csv")
    col = {k: fit[:, i] for i, k in enumerate(names)}
    T1s2, _, T1f2, _ = NBM_T1["PC2"][2]
    table2 = {"T1s_ms": T1s2, "T1f_ms": T1f2, "T1m_ms": NBM_T1M_PC2[2][0], "tq_sq_percent": 100 * NBM_TQ1_PC2[2]}
    for key, value in table2.items():           # Table 2 gives two decimals: half a unit of the last, plus the axis
        cal = THESIS_FIG411_CAL["percent" if key == "tq_sq_percent" else "ms"]
        assert abs(col[key][-1] - value) <= 0.005 + cal, f"Fig. 4.11's last {key} is Table 2's"
    print(f"  the thesis's Fig. 4.11, the paper's Fig. 5 drawn as vector paths (2 %, PC2): points added "
          f"{', '.join(f'{v:g}' for v in col['points_added'])}; at the last two, "
          + " and ".join(f"T1s {col['T1s_ms'][i]:.3f}, T1f {col['T1f_ms'][i]:.3f}, T1m {col['T1m_ms'][i]:.3f} ms, ratio "
                         f"{col['tq_sq_percent'][i]:.3f} %" for i in (-2, -1))
          + f"; Table 2: {table2['T1s_ms']:.2f}, {table2['T1f_ms']:.2f}, {table2['T1m_ms']:.2f} ms, "
          f"{table2['tq_sq_percent']:.2f} %, met by the last point within the table's rounding; with no points added "
          f"{col['T1s_ms'][0]:.3f}, {col['T1f_ms'][0]:.3f} ms, {col['tq_sq_percent'][0]:.3f} % (read off the paper's raster: "
          f"{NBM_FIG5_FIRST[0]}, {NBM_FIG5_FIRST[1]} ms, {NBM_FIG5_FIRST[2]} %)")


# ---------------- Part 10: the thesis's second evaluation, its acquisition at fixed times, and the measured FID ----------------
NBM_TQ1_PC0_2 = (3.14, .53)    # NBM Table 2, T1-TQ/SQ in %, PC0, 2 %
NBM_TABLE1_PC2_2 = dict(T1s=99.98, A1s=5.45, T1f=43.06, A1f=94.17, tq=2.93)   # NBM Table 1, the amplitudes unbound, PC2, 2 %
# (times in ms, amplitudes and ratio in %): one time at the fit's upper bound of 100 ms
THESIS_T46_T1M_PC2_4 = (38.73, .62)   # the thesis's Table 4.6, its reprint of NBM Table 2: T1m at 4 %, PC2, where NBM prints
# 38.65 +- 0.75, the digits of its own PC0 cell
FIG414_CAL = {"ms": 1.77e-4, "percent": 3.25e-5}   # the largest residuals of Fig. 4.14's value axes on their tick marks
FID_NAMES = ["s", "A", "Ts", "Tf", "ph1", "b", "ph2", "dc", "D", "TD", "ph3", "w"]
FID_START = [1.0, 0.98, 49.0, 30.0, 0.0, 0.05, 0.0, 0.0, 0.0, 10.0, 0.0, 0.2]
FID_LO = [0.5, 0.3, 1, 1, -10, -1, -10, -1, -1, 0.5, -10, 0.0]
FID_HI = [1.5, 1.5, 500, 500, 20, 1, 20, 1, 1, 500, 20, 1.0]
NO_DQ = dict(D=0.0, TD=10.0, ph3=0.0)
PLATEAU_1 = dict(s=1.0)


def fid_model(p, t, car):
    """The IRTQTPPI FID along the evolution time: NBM Eq. 13 with the physical sign, s times its SQ part (s = 1 is the
    printed equation, whose plateau is 1) and the fast weight w (0.2 printed), the TQ part s b times its shape (b the TQ
    over the plateau), a constant DC and a double quantum at twice the carrier with an exponential envelope (D = 0
    printed). car is the TPPI phase by point."""
    s, A, Ts, Tf, ph1, b, ph2, dc, D, TD, ph3, w = p
    sq = s * (1 - 2 * A * ((1 - w) * np.exp(-t / Ts) + w * np.exp(-t / Tf)))
    tq = s * b * (np.exp(-t / Tf) - np.exp(-t / Ts))
    return np.sin(car + ph1) * sq + np.sin(3 * car + ph2) * tq + dc + D * np.sin(2 * car + ph3) * np.exp(-t / TD)


def fid_fit(y, t, car, fixed, starts, physical=True, Tmax=500.0):
    """Least squares over the parameters not held in `fixed` (name -> value), from each start (name -> value); with
    physical=True only solutions with T1s > T1f are kept (the printed fit's second minimum has the times swapped)."""
    free = [i for i, k in enumerate(FID_NAMES) if k not in fixed]
    lo, hi = np.array(FID_LO, float), np.array(FID_HI, float)
    hi[2] = hi[3] = Tmax
    best = None
    for st in starts:
        full = np.array([fixed.get(k, st.get(k, v)) for k, v in zip(FID_NAMES, FID_START)], float)

        def res(q, full=full):
            p = full.copy()
            p[free] = q
            return fid_model(p, t, car) - y
        x0 = np.clip(full[free], lo[free] + 1e-9, hi[free] - 1e-9)
        r = least_squares(res, x0, bounds=(lo[free], hi[free]), x_scale="jac", max_nfev=4000)
        p = full.copy()
        p[free] = r.x
        if (not physical or p[2] > p[3]) and (best is None or r.cost < best[1].cost):
            best = (p, r)
    p, r = best
    s2 = 2 * r.cost / (len(y) - len(free))
    cov = np.full((len(FID_NAMES), len(FID_NAMES)), np.nan)
    cov[np.ix_(free, free)] = np.linalg.pinv(r.jac.T @ r.jac) * s2
    err = np.sqrt(np.diag(cov))
    return dict(p=p, err=err, cost=r.cost, s2=s2, cov=cov)


def split_err(f):
    """The split's error from a fit's covariance of T1s and T1f."""
    p, C = f["p"], f["cov"][np.ix_([2, 3], [2, 3])]
    g = np.array([1e3 / p[2] ** 2, -1e3 / p[3] ** 2])
    return float(np.sqrt(g @ C @ g))


def phase_starts(pairs, n1=8, ph=None):
    """Starts over pairs of times, from given phases (ph1, ph2) or from a grid of them."""
    grid = [ph] if ph is not None else [(a, b) for a in np.linspace(0, 2 * np.pi, n1, endpoint=False)
                                       for b in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2)]
    return [dict(Ts=a, Tf=b, ph1=x, ph2=z) for a, b in pairs for x, z in grid]


def fid_split(p):
    return 1e3 * (1 / p[3] - 1 / p[2])


def fid_t1m(p):
    return 1 / (0.8 / p[2] + 0.2 / p[3])


def reading_i_split(pc, t1m, c, T1m=None):
    """Part 8 (f)'s reading (i): one pool at the table's T1m whose TQ, read through the published shape, gives TQ1 (the
    measured ratio over Fig. 6C's, times 3/4); its split in 1/s with the ends from the T1-TQ/SQ error."""
    from scipy.optimize import brentq
    T1s, _, T1f, _ = NBM_T1[pc][c]
    T1m = t1m[c][0] if T1m is None else T1m
    TQ1 = 0.75 * NBM_TQ1[pc][c][0] / NBM_TQ2_FIG6C[c]
    sTQ1 = 0.75 * NBM_TQ1[pc][c][1] / NBM_TQ2_FIG6C[c]
    u = np.exp(-T1_GRID / T1s) - np.exp(-T1_GRID / T1f)

    def proj(delta):
        return (u @ (np.exp(-T1_GRID * (1 / T1m - 0.2 * delta)) - np.exp(-T1_GRID * (1 / T1m + 0.8 * delta)))) / (u @ u)
    dl = brentq(lambda x: proj(x) - TQ1, 1e-9, 0.05)
    dlo = brentq(lambda x: proj(x) - max(TQ1 - sTQ1, 1e-6), 1e-12, 0.05)
    dhi = brentq(lambda x: proj(x) - (TQ1 + sTQ1), 1e-9, 0.05)
    return 1e3 * dl, 1e3 * dlo, 1e3 * dhi


def part10(split_i=None):
    import sympy as sp
    from sympy.physics.quantum.spin import Rotation
    print("\nPart 10. The thesis's second evaluation of the same data (its Eq. 3.16, Table 4.7 and Fig. 4.14), its acquisition "
          "at fixed evolution times, and the measured FID of its Fig. 3.6a (2 %, PC0)")

    # (a) Eq. 3.16 against the count of Part 8 (a), term by term
    def wd(j, m, mp):
        return sp.Abs(sp.simplify(Rotation.d(j, m, mp, sp.pi / 2).doit()))
    W = wd(3, 3, 0) * wd(3, -1, 3) / (wd(1, 1, 0) * wd(1, -1, 1))
    eq316, count = sp.Rational(2, 5) * W, sp.Rational(4, 5) * W   # Eq. 3.16 without its leading sqrt2, and the count
    assert sp.simplify(count - sp.sqrt(6) / 4) == 0, "the count of Part 8 (a)"
    c316 = float(eq316)
    print(f"  (a) Eq. 3.16 is sqrt2 x 2/5 x the Wigner bracket of the count ({sp.nsimplify(W)}) x d read at the line with T2* "
          f"and the fitted amplitudes, where the count has 4/5 and no sqrt2: per |d|, without its sqrt2 (as the thesis reports "
          f"it) {eq316}, with it {sp.nsimplify(sp.sqrt(2) * eq316)}, against the count's {count} (the count, gated against the "
          f"sequences simulated in Part 8 (b), is the TQ over the single quantum at equilibrium, the plateau of Eq. 13): "
          f"{sp.simplify(eq316 / count)} and {sp.nsimplify(sp.sqrt(2) * eq316 / count)} of it, by the 2/5 transcribed from the "
          f"thesis")
    dstar = {}
    for c in (2, 4, 6):
        rep = THESIS_T47["standard"][c][6]
        dstar[c] = rep / (100 * c316)
    dv = {c: d_variants(c) for c in (2, 4, 6)}
    print("      |d*| behind the refit's reported TQ/SQ (Table 4.7, the sqrt2 removed after the fit): " + ", ".join(
          f"{c}% {dstar[c]:.4f} ({THESIS_T47['standard'][c][6] / (100 * np.sqrt(2) * c316):.4f} if kept)" for c in (2, 4, 6))
          + "; the range of Part 8 (d): " + ", ".join(f"{min(dv[c].values()):.3f} to {max(dv[c].values()):.3f}" for c in (2, 4, 6))
          + "; Eq. 3.16 with the second paper's T2* and the weights 0.4 / 0.6 of the thesis's Eq. 3.17 gives " + ", ".join(
          f"{100 * c316 * dv[c]['T2* formula']:.2f}" for c in (2, 4, 6)) + " % where the table reports " + ", ".join(
          f"{THESIS_T47['standard'][c][6]:.2f}" for c in (2, 4, 6)) + " %")
    print("      Part 8 (e)'s size route at |d*|, TQ1 = measured T1-TQ/SQ / (sqrt6/4 |d*|): " + "; ".join(
          f"{c}% " + ", ".join(f"{pc} {NBM_TQ1[pc][c][0] / (100 * LONG * dstar[c]):.3f} +- {NBM_TQ1[pc][c][1] / (100 * LONG * dstar[c]):.3f} "
                               f"(ratio route {0.75 * NBM_TQ1[pc][c][0] / NBM_TQ2_FIG6C[c]:.3f})" for pc in ("PC1", "PC2"))
          for c in (2, 4, 6)))

    # (b) the splits: Table 2, reading (i), the refit and the acquisition at fixed times
    if split_i is not None:
        for c in (2, 4, 6):
            assert reading_i_split("PC2", NBM_T1M_PC2, c)[0] == 1e3 * split_i[("PC2", c)], "the same reading (i) as Part 8 (f)"
    print("  (b) splits in 1/s (PC2); the refit's and the acquisition's errors carry the correlation their T1m error allows (Part 1):")
    for c in (2, 4, 6):
        T2s, sT2s = MRM_T2[c]["T2s"]
        fx, st = thesis_t47("fixed", c, T2s, sT2s), thesis_t47("standard", c, T2s, sT2s)
        di, dlo, dhi = reading_i_split("PC2", NBM_T1M_PC2, c)
        si = (dhi - dlo) / 2
        print(f"    {c}%: reading (i) {di:.2f} ({dlo:.2f} to {dhi:.2f}); refit {st['split']:.2f} +- {st['s_split']:.2f}; acquisition "
              f"{fx['split']:.2f} +- {fx['s_split']:.2f}; acquisition - (i) {(fx['split'] - di) / fx['s_split']:+.2f} sigma (its error), "
              f"refit - acquisition {(st['split'] - fx['split']) / np.hypot(st['s_split'], fx['s_split']):+.2f}, refit - (i) "
              f"{(st['split'] - di) / np.hypot(st['s_split'], si):+.2f} (both errors)")
    a, b_ = reading_i_split("PC2", NBM_T1M_PC2, 4)[0], reading_i_split("PC2", NBM_T1M_PC2, 4, THESIS_T46_T1M_PC2_4[0])[0]
    print(f"      at 4 % NBM's T1m in PC2 ({NBM_T1M_PC2[4][0]} +- {NBM_T1M_PC2[4][1]}) carries the digits of its PC0 cell; the "
          f"thesis's reprint gives {THESIS_T46_T1M_PC2_4[0]} +- {THESIS_T46_T1M_PC2_4[1]}: reading (i)'s split {a:.3f} and {b_:.3f}")

    # (c) the acquisition's T1-TQ/SQ: its normalization is not stated (Eq. 3.18 writes the SQ as 1 - A_SQ(...), A_SQ about 2)
    print("  (c) the acquisition's T1-TQ/SQ against one pool's amplitude at |d*|, in the two normalizations and by the two counts:")
    for c in (2, 4, 6):
        v, e = THESIS_T47["fixed"][c][6:8]
        cells = []
        for lab, k in (("over the plateau", 1.0), ("over A_SQ = 2", 2.0)):
            for clab, cc in (("count", LONG), ("Eq. 3.16", c316)):
                cells.append(f"{lab}, {clab} {100 * cc * dstar[c]:.2f} %: {(k * v - 100 * cc * dstar[c]) / (k * e):+.2f} sigma")
        print(f"    {c}%: {v} +- {e} %; " + "; ".join(cells))

    # (d) the measured FID: the paper's fit identified on it
    names, fid = read_csv("thesis_fig3_6a_irtqtppi_fid_2pct_pc0.csv")
    y, t, car = fid[:, 1], T1_GRID, np.pi / 4 * np.arange(1, 641)
    wt = 2 * np.pi * (45 / 360) / 0.2 * t        # omega t of Eq. 13 with the dense step of 0.2 ms
    dev = np.max(np.abs(np.exp(1j * car) - np.exp(1j * wt)))
    assert dev <= 4 * EPS * np.max(wt), "on the drawn axis omega t of Eq. 13 is the phase stepped by point"
    print(f"  (d) the FID (data/reichert2024_sodium_agarose): the TPPI phase steps by 45 degrees per point in both bands, and "
          f"omega t of Eq. 13 meets it on the drawn axis to {dev:.1e} (a coarse step of nine dense steps is 405 = 45 degrees); "
          f"the phase is the pulses', by point, so the carrier does not see the times")
    pairs = ((49.44, 29.92), (46.0, 40.0), (44.5, 43.5))
    printed = "Eq. 13 as printed (plateau 1, no DQ)"
    forms = {printed: dict(PLATEAU_1, **NO_DQ, w=0.2),
             "plateau 1 with a DQ term": dict(PLATEAU_1, w=0.2),
             "plateau free, no DQ": dict(NO_DQ, w=0.2),
             "plateau free with a DQ term": dict(w=0.2)}
    fits = {k: fid_fit(y, t, car, fx, phase_starts(pairs)) for k, fx in forms.items()}
    fits["Eq. 13 as printed, either order"] = fid_fit(y, t, car, forms[printed], phase_starts(pairs + ((30.0, 49.0),)), physical=False)
    T1s0, sT1s0, T1f0, sT1f0 = NBM_T1["PC0"][2]
    devs = {k: max(abs(f["p"][2] - T1s0), abs(f["p"][3] - T1f0)) for k, f in fits.items()}
    assert devs[printed] < min(v for k, v in devs.items() if k != printed), "the printed form is the one that meets the cell"
    for k, f in fits.items():
        p, e = f["p"], f["err"]
        print(f"      {k}: T1s {p[2]:.3f} +- {e[2]:.3f}, T1f {p[3]:.3f} +- {e[3]:.3f} ms, T1m {fid_t1m(p):.3f}, split {fid_split(p):+.2f} 1/s, "
              f"TQ {100 * abs(p[5]):.2f} +- {100 * e[5]:.2f} % of the plateau, plateau {p[0]:.4f} of the largest point, rms {np.sqrt(f['s2']):.5f}")
    base = fits[printed]
    p = base["p"]
    ph0 = (p[4], p[6])
    drawn_t = fid[:, 0]
    pd = fid_fit(y, drawn_t, car, forms[printed], phase_starts(pairs, ph=ph0))["p"]
    tq0 = fid_fit(y, t, car, dict(forms[printed], b=0.0), phase_starts(pairs, ph=ph0))["p"]
    tqt = fid_fit(y, t, car, dict(forms[printed], b=np.sign(p[5]) * NBM_TQ1_PC0_2[0] / 100), phase_starts(pairs, ph=ph0))["p"]
    gl = fits["Eq. 13 as printed, either order"]
    print(f"    NBM Table 2, PC0, 2 %: T1s {T1s0} +- {sT1s0}, T1f {T1f0} +- {sT1f0} ms (split {1e3 * (1 / T1f0 - 1 / T1s0):.2f} 1/s), T1-TQ/SQ "
          f"{NBM_TQ1_PC0_2[0]} +- {NBM_TQ1_PC0_2[1]} %: the printed form meets the times to {abs(p[2] - T1s0):.3f} and "
          f"{abs(p[3] - T1f0):.3f} ms ({abs(pd[2] - T1s0):.3f} and {abs(pd[3] - T1f0):.3f} on the drawn vertices' own times), every "
          f"other form tried misses by at least {min(v for k, v in devs.items() if k != printed):.3f}; its TQ, {100 * abs(p[5]):.2f} %, "
          f"is {100 * abs(p[5]) / NBM_TQ1_PC0_2[0]:.2f} times the table's, not explained here. The times barely see the TQ: held at 0 "
          f"they are {tq0[2]:.3f} and {tq0[3]:.3f}; held at the table's {NBM_TQ1_PC0_2[0]} % they move to {tqt[2]:.3f} and "
          f"{tqt[3]:.3f}. The table's cell is a local minimum: the printed form's lower one has the times swapped, by delta chi2 "
          f"{2 * (base['cost'] - gl['cost']) / gl['s2']:.0f}")

    # (e) the refit, from Fig. 4.14
    rows414 = {}
    with open(os.path.join(DATA, "thesis_fig4_14_refit_tq_held_2_4_6pct.csv"), encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("PC"):
                v = line.strip().split(",")
                rows414[(v[0], int(v[1]))] = tuple(float(x) for x in v[2:])
    csv_round = 5e-5      # the CSV's last digit
    for c in (2, 4, 6):
        s_, ss, f_, sf, m_, sm, tq, stq = THESIS_T47["standard"][c]
        T1s, T1m, T1f, r414 = rows414[("PC2", c)]
        assert max(abs(T1s - s_), abs(T1m - m_), abs(T1f - f_)) <= 0.005 + FIG414_CAL["ms"] + csv_round and \
            abs(r414 - tq) <= 0.005 + FIG414_CAL["percent"] + csv_round, "Fig. 4.14's PC2 markers are Table 4.7's standard row"
    T1s4, T1m4, T1f4, r4 = rows414[("PC0", 2)]

    def held_fit(k, fx, tt=t):
        return fid_fit(y, tt, car, dict(fx, b=np.sign(p[5]) * k * r4 / 100), phase_starts(pairs, ph=ph0))

    def solve_k(idx, target, tt=t):
        from scipy.optimize import brentq
        return brentq(lambda k: held_fit(k, forms[printed], tt)["p"][idx] - target, 1.2, 1.7, xtol=1e-7)
    ks, kf = solve_k(2, T1s4), solve_k(3, T1f4)
    ksd, kfd = solve_k(2, T1s4, drawn_t), solve_k(3, T1f4, drawn_t)
    h = 1e-3
    dTs = (held_fit(ks + h, forms[printed])["p"][2] - held_fit(ks - h, forms[printed])["p"][2]) / (2 * h)
    dTf = (held_fit(kf + h, forms[printed])["p"][3] - held_fit(kf - h, forms[printed])["p"][3]) / (2 * h)
    tol_k = (FIG414_CAL["ms"] + csv_round) * (1 / abs(dTs) + 1 / abs(dTf))
    assert abs(ks - kf) <= tol_k, "one held amplitude gives the refit's two times, within the figure's reading"
    alt = {"sqrt2, plateau 1": held_fit(np.sqrt(2), forms[printed]), "1, plateau 1": held_fit(1.0, forms[printed]),
           "sqrt2, plateau free": held_fit(np.sqrt(2), dict(NO_DQ, w=0.2)), "1, plateau free": held_fit(1.0, dict(NO_DQ, w=0.2))}
    dev414 = {k: max(abs(f["p"][2] - T1s4), abs(f["p"][3] - T1f4)) for k, f in alt.items()}
    kfit = held_fit((ks + kf) / 2, forms[printed])
    print(f"  (e) Fig. 4.14 (the refit, read off its markers; PC2 is Table 4.7's standard row within its rounding and the axes), "
          f"PC0 at 2 %: T1s {T1s4:.4f}, T1m {T1m4:.4f}, T1f {T1f4:.4f} ms, TQ/SQ {r4:.4f} %. The printed form with the TQ held at k "
          f"times that, on the sampling of T1_GRID: T1s asks k = {ks:.5f}, T1f asks {kf:.5f}, the same within {tol_k:.1e} (the "
          f"figure's reading carried through), {(ks + kf) / 2 / np.sqrt(2):.4f} sqrt2; on the drawn vertices' own times they ask "
          f"{ksd:.5f} and {kfd:.5f}; at k = sqrt2 the times miss by {dev414['sqrt2, plateau 1']:.3f} ms, at k = 1 or with the plateau "
          f"free by at least {min(v for k, v in dev414.items() if k != 'sqrt2, plateau 1'):.3f}. Eq. 3.16 with its sqrt2 is "
          f"{np.sqrt(2) * c316 / LONG:.3f} of the count at the d its reported value implies; delta chi2 of the held amplitude against "
          f"the free fit of the same form {2 * (kfit['cost'] - base['cost']) / base['s2']:.0f} (on that fit's residual variance)")

    # (f) the plateau, and the free plateau on the drawn axis
    fr = fits["plateau free with a DQ term"]
    pf = fr["p"]
    tail = t > 430
    X = np.column_stack([fn(k * car[tail]) for k in (1, 2, 3) for fn in (np.sin, np.cos)] + [np.ones(tail.sum()), np.cos(4 * car[tail])])
    coef = np.linalg.lstsq(X, y[tail], rcond=None)[0]
    floor = np.sqrt(np.sum((y[tail] - X @ coef) ** 2) / (tail.sum() - X.shape[1]))
    plat = float(np.hypot(coef[0], coef[1]))
    t_first = t[tail][0]

    def recovery(q, tt):
        return 2 * q[1] * q[0] * (0.8 * np.exp(-tt / q[2]) + 0.2 * np.exp(-tt / q[3]))
    print(f"  (f) read without a model on the tail (t > 430 ms, where the recovery left is {recovery(pf, t_first):.1e} by the free "
          f"plateau's fit and {recovery(p, t_first):.1e} by the printed one at the first point, {t_first:.1f} ms): the single "
          f"quantum's amplitude {plat:.4f} of the largest point, the DC {coef[6]:+.4f}, the noise floor {floor:.5f}; Eq. 13 holds that "
          f"amplitude at 1. With the plateau free (a scale on the single quantum) on the drawn axis: plateau {pf[0]:.4f}")
    res = y - fid_model(pf, t, car)
    bands = ((0, 64), (64, 82), (82, 100), (100, 150), (150, 430), (430, 640))
    print("    its residual's rms by band (ms): " + ", ".join(
        f"{lo} to {hi}: {np.sqrt(np.mean(res[(t > lo) & (t <= hi)] ** 2)):.5f}" for lo, hi in bands)
        + f": at the floor except just after the bands join")

    def harmonics(idx):
        """The residual (data - model) over consecutive points, on the single quantum in phase with the fitted SQ term and
        out of it, the constant, and the double and triple quantum in size, with the 1-sigma of the noise floor."""
        B = np.column_stack([np.sin(car[idx] + pf[4]), np.cos(car[idx] + pf[4]), np.ones(len(idx)),
                             np.sin(2 * car[idx]), np.cos(2 * car[idx]), np.sin(3 * car[idx]), np.cos(3 * car[idx])])
        c_, *_ = np.linalg.lstsq(B, res[idx], rcond=None)
        sgn = np.sign(pf[0])
        return sgn * c_[0], max(abs(c_[1]), abs(c_[2]), np.hypot(c_[3], c_[4]), np.hypot(c_[5], c_[6])), floor * np.sqrt(2 / len(idx))
    hc, ho, hs = harmonics(np.arange(320, 328))
    hd, hdo, hds = harmonics(np.arange(304, 320))
    print(f"    the residual (data - model) on the carriers: over the first eight coarse points {hc:+.4f} on the single quantum in "
          f"phase (negative: less recovered than the model), at most {ho:.4f} out of phase, at order 0 and in the double and the "
          f"triple quantum (1 sigma {hs:.4f}); over the last sixteen dense points {hd:+.4f} in phase, at most {hdo:.4f} on the rest "
          f"(1 sigma {hds:.4f})")
    prof = []
    for w in (0.10, 0.15, 0.20, 0.25, 0.30, 0.40):
        g = fid_fit(y, t, car, dict(w=w), phase_starts(pairs, ph=(pf[4], pf[6])))
        prof.append(f"{w:.2f}: split {fid_split(g['p']):.2f}, delta chi2 {2 * (g['cost'] - fr['cost']) / fr['s2']:+.1f}")
    gw = fid_fit(y, t, car, {}, phase_starts(pairs, ph=(pf[4], pf[6])), physical=False)
    gu = fid_fit(y, t, car, dict(PLATEAU_1, **NO_DQ), phase_starts(((99.0, 43.0), (49.0, 30.0)), ph=(pf[4], pf[6])), physical=False, Tmax=100.0)
    print(f"    the fast weight held (delta chi2 on the w = 0.2 fit's residual variance): " + "; ".join(prof)
          + f"; free: {gw['p'][11]:.3f} at {gw['p'][3]:.2f} ms beside {gw['p'][2]:.2f} ms")
    slot = sorted(((1 - gu["p"][11], gu["p"][2]), (gu["p"][11], gu["p"][3])), key=lambda z: z[1])
    tb = NBM_TABLE1_PC2_2
    print(f"    plateau 1, the weights unbound and the times bounded at 100 ms (as the paper's Table 1): {100 * slot[1][0]:.2f} % at "
          f"{slot[1][1]:.2f} ms, {100 * slot[0][0]:.2f} % at {slot[0][1]:.2f} ms, TQ {100 * abs(gu['p'][5]):.2f} %; the paper's Table 1 "
          f"(PC2): {tb['A1s']} % at {tb['T1s']} ms, {tb['A1f']} % at {tb['T1f']} ms, TQ {tb['tq']} %")

    def part_terms(q):
        s, A, Ts, Tf, ph1, b, ph2 = q[:7]
        return (np.sin(car + ph1) * s * (1 - 2 * A * (0.8 * np.exp(-t / Ts) + 0.2 * np.exp(-t / Tf))),
                np.sin(3 * car + ph2) * s * b * (np.exp(-t / Tf) - np.exp(-t / Ts)))

    def jac(k, idx):
        cols = []
        for i in idx:
            dp = np.zeros(12)
            dp[i] = 1e-6 * max(1.0, abs(pf[i]))
            cols.append((part_terms(pf + dp)[k] - part_terms(pf - dp)[k]) / (2 * dp[i]))
        return np.array(cols).T
    gvec = np.array([1e3 / pf[2] ** 2, -1e3 / pf[3] ** 2])
    sd = []
    for k, idx, at in ((0, [0, 1, 2, 3, 4], [2, 3]), (1, [2, 3, 5, 6], [0, 1])):
        J = jac(k, idx)
        C = np.linalg.pinv(J.T @ J) * fr["s2"]
        sd.append(np.sqrt(gvec @ C[np.ix_(at, at)] @ gvec))
    print(f"    the split's spread from the SQ term alone {sd[0]:.2f}, from the TQ term alone {sd[1]:.1f} 1/s: the single quantum "
          f"sets it")

    # (g) the join of the two bands
    fixed_dq = dict(w=0.2)
    phf = (pf[4], pf[6])

    def shifted(dt):
        tt = t.copy()
        tt[320:] += dt
        return tt

    def offset_fit(yy, phs, fixed=fixed_dq):
        from scipy.optimize import minimize_scalar
        r = minimize_scalar(lambda dt: fid_fit(yy, shifted(dt), car, fixed, phase_starts(pairs[:2], ph=phs))["cost"],
                            bounds=(-3.0, 2.0), method="bounded", options=dict(xatol=1e-4))
        return r.x, fid_fit(yy, shifted(r.x), car, fixed, phase_starts(pairs[:2], ph=phs))
    dt0, fo = offset_fit(y, phf)
    po = fo["p"]
    curv = [2 * fid_fit(y, shifted(dt0 + s), car, fixed_dq, phase_starts(pairs[:2], ph=phf))["cost"] / fo["s2"] for s in (-0.05, 0.05)]
    sdt = 0.05 / np.sqrt((curv[0] + curv[1] - 2 * 2 * fo["cost"] / fo["s2"]) / 2)

    def band_step(kind):
        """The coarse band changed against the dense one by one factor: its part still to recover ('unrecovered', the
        inversion term times 1 + e), its whole signal ('gain'), or both its times ('T1')."""
        free = [i for i, k in enumerate(FID_NAMES) if k != "w"]
        lo, hi = np.array(FID_LO, float)[free], np.array(FID_HI, float)[free]

        def model(q):
            full = pf.copy()
            full[free] = q[:-1]
            e = q[-1]
            m = fid_model(full, t, car)
            c2 = full.copy()
            if kind == "unrecovered":
                c2[1] *= 1 + e
            elif kind == "gain":
                c2[0] *= 1 + e
            else:
                c2[2] *= 1 + e
                c2[3] *= 1 + e
            m[320:] = fid_model(c2, t, car)[320:]
            return m
        r = least_squares(lambda q: model(q) - y, np.append(pf[free], 0.0), bounds=(np.append(lo, -0.5), np.append(hi, 0.5)),
                          x_scale="jac", max_nfev=4000)
        q = pf.copy()
        q[free] = r.x[:-1]
        return r.x[-1], q, r.cost
    steps = {k: band_step(k) for k in ("unrecovered", "gain", "T1")}
    srow = "; ".join(f"{k}: {100 * e:+.2f} %, delta chi2 {2 * (fr['cost'] - cst) / fo['s2']:.0f}, split {fid_split(q):.2f}"
                     for k, (e, q, cst) in steps.items())
    prow = [f"{d:+.1f}: split {fid_split(fid_fit(y, shifted(d), car, fixed_dq, phase_starts(pairs[:2], ph=phf))['p']):.2f}, delta chi2 "
            f"{2 * (fid_fit(y, shifted(d), car, fixed_dq, phase_starts(pairs[:2], ph=phf))['cost'] - fo['cost']) / fo['s2']:+.1f}"
            for d in (-1.8, -1.6, -1.4, -1.2)]
    print(f"  (g) the join, plateau free with a DQ term, delta chi2 gained over the unmoved fit on the variance at the floor: the "
          f"coarse band moved in time by {dt0:+.3f} +- {sdt:.3f} ms gains {2 * (fr['cost'] - fo['cost']) / fo['s2']:.0f} (rms "
          f"{np.sqrt(fo['s2']):.5f} against {np.sqrt(fr['s2']):.5f}, the floor {floor:.5f}), split {fid_split(po):.2f} +- "
          f"{split_err(fo):.2f}, T1m {fid_t1m(po):.2f}; the coarse band changed by one factor instead, {srow} (a move by dt scales "
          f"the part still to recover by exp(-dt/T), so the first two are one step); at other moves " + "; ".join(prow))
    from scipy.optimize import minimize

    def timing(free_off, free_cs):
        def make_t(q):
            off = q[0] if free_off else 0.0
            cs = q[-1] if free_cs else 1.8
            return np.concatenate([0.2 * np.arange(1, 321), 64.0 + off + cs * np.arange(1, 321)])
        q0 = [dt0] * free_off + [1.8] * free_cs
        r = minimize(lambda q: fid_fit(y, make_t(q), car, fixed_dq, phase_starts(pairs[:2], ph=phf))["cost"], q0,
                     method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-12, maxiter=600))
        return r.x, 2 * (fr["cost"] - r.fun) / fo["s2"]
    q1, g1 = timing(0, 1)
    q2, g2 = timing(1, 1)
    print(f"    the time scale: the coarse step free, {q1[0]:.4f} ms, gains {g1:.0f} (the dense step free is the same fit, the "
          f"model being unchanged when all times and both T1 scale together); the move with the coarse step free, "
          f"{q2[0]:+.4f} ms and {q2[1]:.4f} ms, gains {g2:.0f}")
    rng = np.random.default_rng(2)
    nulls = []
    for split in (6.0, 9.5):
        Rs = 1 / fid_t1m(pf) - 0.2 * split / 1e3
        truth = pf.copy()
        truth[0], truth[1], truth[2], truth[3] = 1.0, 1.0, 1 / Rs, 1 / (Rs + split / 1e3)
        for _ in range(4):
            yy = fid_model(truth, t, car) + floor * rng.standard_normal(640)
            nulls.append(offset_fit(yy / np.max(np.abs(yy)), phf)[0])
    nulls = np.array(nulls)
    assert np.max(np.abs(nulls)) <= 3 * sdt, "on controls without a step the fit finds none, within three of its own errors"
    print(f"    controls without a step (one pool, split 6.0 and 9.5, the curve's DC, DQ, TQ and noise floor, normalized to the "
          f"largest point, eight draws): the fitted move {nulls.mean():+.3f}, at most {np.max(np.abs(nulls)):.3f} ms in size "
          f"against its error {sdt:.3f}")
    keep = ~((t > 64) & (t <= 150))
    dense = np.arange(640) < 320
    routes = {("the dense band alone, the plateau held at the tail's", True):
              fid_fit(y[dense], t[dense], car[dense], dict(fixed_dq, s=plat), phase_starts(pairs, ph=phf)),
              ("64 to 150 ms left out, the plateau free", True): fid_fit(y[keep], t[keep], car[keep], fixed_dq, phase_starts(pairs, ph=phf)),
              ("the coarse band moved, the plateau free", True): fo,
              ("the dense band alone, the plateau held at the tail's", False):
              fid_fit(y[dense], t[dense], car[dense], dict(fixed_dq, s=plat, **NO_DQ), phase_starts(pairs, ph=phf)),
              ("64 to 150 ms left out, the plateau free", False): fid_fit(y[keep], t[keep], car[keep], dict(fixed_dq, **NO_DQ),
                                                                       phase_starts(pairs, ph=phf))}
    with_dq = [fid_split(f["p"]) for k, f in routes.items() if k[1]]
    print("    routes past the join, the fast weight held at 0.2: " + "; ".join(
        f"{k[0]}, {'with' if k[1] else 'no'} DQ: split {fid_split(f['p']):.2f}, T1m {fid_t1m(f['p']):.2f}, rms {np.sqrt(f['s2']):.5f}"
        for k, f in routes.items()) + f"; with the DQ term, at the floor, {min(with_dq):.2f} to {max(with_dq):.2f} 1/s")
    dfree = fid_fit(y[dense], t[dense], car[dense], fixed_dq, phase_starts(pairs, ph=phf))
    gws = fid_fit(y, shifted(dt0), car, {}, phase_starts(pairs, ph=phf), physical=False)
    wsplit = 1e3 * abs(1 / gws["p"][3] - 1 / gws["p"][2])
    print(f"    the inversion of the moved fit, A = {po[1]:.4f} +- {fo['err'][1]:.4f} (1 is complete), so a coarse band 3.36 % further "
          f"from recovery with the same gain would be inverted beyond its equilibrium")
    print(f"    the dense band alone with its plateau free: split {fid_split(dfree['p']):.2f} at plateau {dfree['p'][0]:.4f} +- "
          f"{dfree['err'][0]:.4f}; on the moved axis with the fast weight free: {gws['p'][11]:.3f} at {gws['p'][3]:.2f} ms beside "
          f"{gws['p'][2]:.2f} ms, split {wsplit:.2f}, and 0.2 costs delta chi2 {2 * (fo['cost'] - gws['cost']) / gws['s2']:.1f}")

    # (h) the late points: the curve against Fig. 4.11, and loops of known truth
    nf, fit411 = read_csv("thesis_fig4_11_fit_vs_points_added_2pct_pc2.csv")
    col = {k: fit411[:, i] for i, k in enumerate(nf)}
    pts = col["points_added"].astype(int)
    rows, curve_series = [], {}
    for i, m in enumerate(pts):
        mask = np.arange(640) < 320 + m
        a1 = fid_fit(y[mask], t[mask], car[mask], forms[printed], phase_starts(pairs, ph=ph0))["p"]
        a2 = fid_fit(y[mask], t[mask], car[mask], forms["plateau free, no DQ"], phase_starts(pairs, ph=ph0))["p"]
        curve_series[m] = fid_split(a1)
        rows.append(f"{m}: {a1[2]:.2f} / {a1[3]:.2f} / {100 * abs(a1[5]):.2f} (split {fid_split(a1):.2f}), free {a2[2]:.2f} / {a2[3]:.2f}, "
                    f"Fig. 4.11 {col['T1s_ms'][i]:.2f} / {col['T1f_ms'][i]:.2f} / {col['tq_sq_percent'][i]:.2f}")
    print("  (h) coarse points added: the FID (PC0) with Eq. 13 as printed, T1s / T1f ms / TQ %, with the plateau free, and the "
          "thesis's Fig. 4.11 (PC2): " + "; ".join(rows))
    wins = (0, 30, 100, 320)
    rng = np.random.default_rng(3)

    def loop(split, tt, noise, draws, src):
        Rs = 1 / fid_t1m(src) - 0.2 * split / 1e3
        truth = src.copy()
        truth[0], truth[1], truth[2], truth[3] = 1.0, 1.0, 1 / Rs, 1 / (Rs + split / 1e3)
        clean = fid_model(truth, tt, car)
        out = {m: [] for m in wins}
        for _ in range(draws):
            yy = clean + noise * rng.standard_normal(640)
            yy = yy / np.max(np.abs(yy))
            for m in wins:
                mask = np.arange(640) < 320 + m
                out[m].append(fid_split(fid_fit(yy[mask], t[mask], car[mask], forms[printed], phase_starts(pairs, ph=ph0))["p"]))
        return {m: np.array(v) for m, v in out.items()}
    split_mid = float(np.median(with_dq))
    arms = {"moved and noisy": loop(split_mid, shifted(dt0), floor, 16, po), "moved, no noise": loop(split_mid, shifted(dt0), 0.0, 1, po),
            "unmoved and noisy": loop(split_mid, t, floor, 16, po), "unmoved, no noise": loop(split_mid, t, 0.0, 1, po),
            "the drawn axis's split, unmoved and noisy": loop(fid_split(pf), t, floor, 16, pf)}
    jump = {k: a[30] - a[0] for k, a in arms.items()}
    cj = curve_series[30] - curve_series[0]
    a = arms["moved and noisy"]
    for m in wins:
        assert abs(a[m].mean() - curve_series[m]) <= 3 * a[m].std(), "the step and the plateau give the curve's series"
    assert abs(jump["moved and noisy"].mean() - cj) <= 3 * jump["moved and noisy"].std(), "with the step, the curve's jump at 30"
    assert jump["unmoved and noisy"].mean() + 3 * jump["unmoved and noisy"].std() < cj, "without the step, no such jump"
    print(f"    loops: one pool of split {split_mid:.2f} 1/s (the median of the routes at the floor) at T1m {fid_t1m(po):.2f} ms with the "
          f"curve's DC, DQ and TQ from its moved fit, generated with and without the move of {dt0:+.3f} ms and the noise floor, "
          f"normalized to the largest point, fitted with Eq. 13 as printed on the drawn axis (+- the spread of sixteen draws): "
          + "; ".join(f"{k}: " + ", ".join(f"{m}: {v[m].mean():.2f}" + (f" +- {v[m].std():.2f}" if len(v[m]) > 1 else "") for m in wins)
                      for k, v in arms.items())
          + f"; the curve: " + ", ".join(f"{m}: {curve_series[m]:.2f}" for m in wins) + f". The jump from no coarse points to 30: "
          + ", ".join(f"{k} {v.mean():.2f}" for k, v in jump.items()) + f", the curve {cj:.2f}. The pool is the curve's own fit, so "
          f"the loops show what the step and the plateau do, not which split is the curve's")

    # (i) controls of known split for the plateau alone
    rng = np.random.default_rng(4)
    out = []
    for split in (2.0, fid_split(fits["plateau free, no DQ"]["p"])):
        Rs = 1 / fid_t1m(pf) - 0.2 * split / 1e3
        truth = pf.copy()
        truth[0], truth[1], truth[2], truth[3], truth[8] = 1.0, 1.0, 1 / Rs, 1 / (Rs + split / 1e3), 0.0
        clean = fid_model(truth, t, car)
        cn = clean / np.max(np.abs(clean))
        g1 = fid_split(fid_fit(cn, t, car, forms[printed], phase_starts(pairs, ph=ph0))["p"])
        g0 = fid_split(fid_fit(cn, t, car, forms["plateau free, no DQ"], phase_starts(pairs, ph=ph0))["p"])
        assert abs(g0 - split) <= 1e-8 * split, "with the plateau free the noiseless control returns its split, to the solver's xtol"
        noisy = {True: [], False: []}
        for _ in range(16):
            yy = clean + floor * rng.standard_normal(640)
            yy = yy / np.max(np.abs(yy))
            for pl in (True, False):
                fx = forms[printed] if pl else forms["plateau free, no DQ"]
                noisy[pl].append(fid_split(fid_fit(yy, t, car, fx, phase_starts(pairs, ph=ph0))["p"]))
        n1, n0 = np.array(noisy[True]), np.array(noisy[False])
        out.append(f"split {split:.2f}: without noise, plateau {1 / np.max(np.abs(clean)):.4f}, printed {g1:.2f}, free {g0:.4f} "
                   f"(off by {abs(g0 - split) / split:.0e} of it); with the noise floor printed {n1.mean():.2f} +- {n1.std() / 4:.2f}, "
                   f"free {n0.mean():.2f} +- {n0.std() / 4:.2f} (+- the error of the mean of sixteen draws)")
    print("  (i) controls for the plateau alone, one pool on the drawn axis with the curve's DC and TQ, normalized to its largest "
          "point, all coarse points: " + "; ".join(out))

    # (j) the held amplitude, the relation and the acquisition, on both axes
    T2s, sT2s = MRM_T2[2]["T2s"]
    fxa = thesis_t47("fixed", 2, T2s, sT2s)
    dv2 = d_variants(2)
    for axis, tt, ref in (("drawn axis", t, fr), ("moved axis", shifted(dt0), fo)):
        cells = []
        for lab, b in (("the refit's", (ks + kf) / 2 * r4 / 100), ("one pool's by the count at the smallest |d| of Part 8 (d)", LONG * min(dv2.values())),
                       ("at its |d*|", LONG * r4 / (100 * c316)), ("at the largest |d|", LONG * max(dv2.values()))):
            hf = fid_fit(y, tt, car, dict(w=0.2, b=np.sign(pf[5]) * b), phase_starts(pairs, ph=phf))
            hq, he = hf["p"], hf["err"]
            _, _, hg, hs = gap(hq[2], he[2], hq[3], he[3], T2s, sT2s)
            cells.append(f"{lab} ({100 * b:.2f} %): split {fid_split(hq):.2f}, delta chi2 {2 * (hf['cost'] - ref['cost']) / ref['s2']:.0f}, "
                         f"gap {hg:.2f} +- {hs:.2f} 1/s = {hg / hs:.1f} sigma")
        q, e = ref["p"], ref["err"]
        R, R2, gp, sg = gap(q[2], e[2], q[3], e[3], T2s, sT2s)
        print(f"  (j) {axis}, plateau free with a DQ term and the TQ free: split {fid_split(q):.2f}, T1m {fid_t1m(q):.2f}, TQ {100 * abs(q[5]):.2f} % "
              f"(one pool's by the count at |d| {abs(q[5]) / LONG:.3f}); the TQ held at " + "; ".join(cells)
              + f"; the relation against the second paper's T2s: gap {gp:.2f} +- {sg:.2f} 1/s = {gp / sg:.1f} sigma; the acquisition "
              f"at fixed times (PC2, a separate session, Eq. 3.18 frees the plateau) {fxa['split']:.2f} +- {fxa['s_split']:.2f} at T1m "
              f"{THESIS_T47['fixed'][2][4]} ms: {(fid_split(q) - fxa['split']) / fxa['s_split']:.1f} sigma on its error")
    for k, f in routes.items():
        q, e = f["p"], f["err"]
        R, R2, gp, sg = gap(q[2], e[2], q[3], e[3], T2s, sT2s)
        print(f"      {k[0]}, {'with' if k[1] else 'no'} DQ: gap {gp:.2f} +- {sg:.2f} 1/s = {gp / sg:.1f} sigma")
    R, R2, gp, sg = gap(T1s0, sT1s0, T1f0, sT1f0, T2s, sT2s)
    print(f"    with NBM Table 2's PC0 cell the gap is {gp:.2f} +- {sg:.2f} 1/s")

if __name__ == "__main__":
    do_global = "--global" in sys.argv
    part1()
    part2()
    best = part3(do_global)
    part4(best)
    chi_at = part5()
    edges = part6(chi_at)
    part7()
    split_i = part8(edges)
    part9()
    part10(split_i)
