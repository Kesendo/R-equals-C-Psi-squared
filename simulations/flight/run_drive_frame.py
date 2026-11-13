"""The drive-frame flight, with the X^3 window arms riding in its Marrakesh and Kingston jobs.

Pre-registrations: experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md (arms A, RO, the forest; verdicts P, V1, VF, VS) and
experiments/X3_WINDOW_HARDWARE_PREDICTION.md (arms B, Q; verdicts PB, WB, VFB, VQ). Circuits are frozen per chip by
drive_frame_freeze.py; the runner never re-transpiles them.

    python run_drive_frame.py --gate CHIP                 # offline: structure, logical identities, the zero, the turn
                                                          # law in the replay, one simulated flight file through the
                                                          # analysis (written under results/, gitignored)
    python run_drive_frame.py --nulls CHIP [N]            # N simulated flights through the analysis; only N = 10 000
                                                          # writes the ruling file drive_frame_<chip>_nulls.json
    python run_drive_frame.py --freeze-numbers CHIP       # the per-chip numbers of the X^3 page's section 4
    python run_drive_frame.py --hardware CHIP             # day-of conditions, submit, persist, analyze
    python run_drive_frame.py --recover STUB.json         # rebuild the counts file from a submitted job id
    python run_drive_frame.py --analyze RESULT.json [--reach]   # the committed rules on measured counts

The replay is the October flight's (reality_window_letters.replay, patched with the drive-crosstalk letter by
simulations/reality_window_crosstalk_fit.py), run on logical qubits 0, 1, 2 with the role labels [4, 5, 6] for a, b, c:
every letter is keyed by role, so a letter of the October path is carried to another path by role. T1 and T2 are the
freeze snapshot's on the chip's path; the gate durations are its target's.

Covariances (section 5 of both pages). For an arm (A or B) and a set of its readings, Sigma is the multinomial
covariance of one circuit's parities at its shots, computed from the measured parities pooled over the arm's standing
variants (logically one circuit). The drift: each standing variant's half difference (h1 - h2) / 2 has covariance Sigma
when nothing drifts; in Sigma's whitened basis their mean outer product S, kappa = max(1, trace S / rank), and
K = S with every eigenvalue below kappa raised to kappa; "kappa Sigma" is Sigma^1/2 K Sigma^1/2, formed block by block
over the basis circuits (DriftCov). With S diagonal this is max(kappa, kappa_i) per reading; as a matrix it also takes a
drift that moves the readings of one circuit together, which the per-reading form, on A's correlated readings, does not.
Credentials: IBM_QUANTUM_TOKEN from the environment or a .env here (gitignored), else the saved account.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import re
import subprocess
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular
from scipy.stats import chi2 as CHI2

import drive_frame_freeze as FZ
import reality_window_letters as rl

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "simulations"))
RESULTS = HERE / "results"
OUTDIR = REPO / "simulations" / "results" / "drive_frame"
PAGE = REPO / "experiments" / "DRIVE_FRAME_HARDWARE_PREDICTION.md"

ROLE = [4, 5, 6]
SHOTS_HALF = 4096
N_NULLS = 10000
BASES = [(b0, b2) for b0 in "XYZ" for b2 in "XYZ"]
BI = {b: i for i, b in enumerate(BASES)}
MI = {m: i for i, m in enumerate(rl.MASKS)}
SG = np.array([[(-1.0) ** (m[0] * (k >> 2 & 1) + m[1] * (k >> 1 & 1) + m[2] * (k & 1)) for k in range(8)]
               for m in rl.MASKS])                       # parity signs, rows in rl.MASKS order, columns k = 4 q0 + 2 q2 + q1
XOR = np.array([[MI.get(tuple(a ^ b for a, b in zip(m1, m2)), -1) for m2 in rl.MASKS] for m1 in rl.MASKS])

# Section 5 of the drive-frame page and of the X^3 window page
P_TURNS, P_BLIND = 0.999, 0.99                           # TURNS above the first, FRAME-BLIND at or below the second
P_OPEN, P_SHUT = 0.999, 0.95
SIM_RULE = 1.05                                          # a simulated point more than 5 % above its chi^2 point rules
P_GATE_MARGIN = 0.05
PB_GATE = (0.7, 0.7, -0.6)
RO_GUARD_P1G0 = 0.01
RO_GUARD_ERR = 0.03
Q_BOUNDS = ([-1.5, -10.0, 0.0, 0.0], [1.5, 10.0, 1.2, 1.0])
Q_DETUNED_Z = 4.0                                        # the phase wander rests on six half differences: a t-like tail
REACH_DRAWS = 2000

# Section 6
MIN_REMAINING_S = 600
MAX_EXEC_S = 600


def say(lines, s=""):
    print(s, flush=True)
    lines.append(s)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def forbidden_b(b, m):
    """X^3 window page section 1: the parity's string carries an odd number of Y and Z letters on q0 and q2, or an
    odd number of Y (q1 read in X carries neither)."""
    letters = [b[0] if m[0] else "I", b[1] if m[1] else "I"]
    return sum(c in "YZ" for c in letters) % 2 == 1 or sum(c == "Y" for c in letters) % 2 == 1


def no_x(b, m):
    return m[2] == 0 and not (m[0] and b[0] == "X") and not (m[1] and b[1] == "X")


A_IDX = [(b, m) for b in BASES for m in rl.MASKS if rl.odd_y(b, m)]                 # 36 odd-Y readings of an A variant
B_IDX = [(b, m) for b in BASES for m in rl.MASKS if forbidden_b(b, m) and no_x(b, m)]   # the 14 of the X^3 page
B_OTHER = [(b, m) for b in BASES for m in rl.MASKS if forbidden_b(b, m) and not no_x(b, m)]   # the other 22
S_POS = A_IDX.index((("Y", "Z"), (1, 1, 0)))
B_KINDS = {"one Z, no Y": [i for i, (b, m) in enumerate(B_IDX) if not any(m[l] and b[l] == "Y" for l in (0, 1))],
           "one Y, no Z": [i for i, (b, m) in enumerate(B_IDX) if not any(m[l] and b[l] == "Z" for l in (0, 1))],
           "one Y and one Z": [i for i, (b, m) in enumerate(B_IDX) if sum(m[:2]) == 2]}
POOLED = {f"{P}{site}": [i for i, (b, m) in enumerate(B_IDX) if sum(m) == 1 and m[l] and b[l] == P]
          for P in "YZ" for site, l in (("0", 0), ("2", 1))}
assert len(A_IDX) == 36 and len(B_IDX) == 14 and len(B_OTHER) == 22
assert sorted(sum(B_KINDS.values(), [])) == list(range(14)) and all(len(v) == 3 for v in POOLED.values())

FIVE = {("det", 4): 0.001, ("det", 5): 0.201, ("det", 6): 0.180, ("zz", (4, 5)): 0.002, ("zz", (5, 6)): -0.111}
XT_PHASE = {(5, 4): -2.24, (5, 6): 0.83}                  # the October crosstalk-alone fit's phases on b's two pairs
XT_FIT = {("xtc", (5, 4)): 0.0079 * math.cos(-2.24), ("xts", (5, 4)): 0.0079 * math.sin(-2.24),
          ("xtc", (5, 6)): 0.0090 * math.cos(0.83), ("xts", (5, 6)): 0.0090 * math.sin(0.83),
          ("xtc", (4, 5)): 0.0106 * math.cos(-2.56), ("xts", (4, 5)): 0.0106 * math.sin(-2.56),
          ("xtc", (6, 5)): 0.0124 * math.cos(2.69), ("xts", (6, 5)): 0.0124 * math.sin(2.69)}


def xt_b(eps):
    """The reach's letter: eps on (b -> a) and (b -> c) at the October phases, the other pairs at zero."""
    L = {}
    for pair, ph in XT_PHASE.items():
        L[("xtc", pair)] = eps * math.cos(ph)
        L[("xts", pair)] = eps * math.sin(ph)
    return L


_REPLAY = {}


def replay_fn(L):
    if L and any(k[0] in ("xtc", "xts", "stx") for k in L):
        if "patched" not in _REPLAY:
            import reality_window_crosstalk_fit as CF
            _REPLAY["patched"] = CF.patched_replay()
        return _REPLAY["patched"]
    return rl.replay


def pinv_rank(S):
    w, V = np.linalg.eigh(S)
    keep = w > 1e-9 * w.max()
    return (V[:, keep] / w[keep]) @ V[:, keep].T, int(keep.sum())


# ── a chip: its frozen file, its freeze record, its replay ────────────────

class Chip:
    def __init__(self, name):
        if name not in FZ.CHIPS:
            raise SystemExit(f"unknown chip {name}")
        self.name = name
        self.qpy, self.meta = FZ.files(name)
        self.freeze = json.loads(self.meta.read_text(encoding="utf-8"))
        self.freeze_sha = FZ.canonical_sha(self.freeze)
        if self.freeze["sha256"] != sha256(self.qpy):
            raise SystemExit(f"{name}: the frozen file does not match its freeze record")
        from qiskit import qpy
        with open(self.qpy, "rb") as f:
            self.circs = {c.name: c for c in qpy.load(f)}
        self.path, self.forest = self.freeze["path"], self.freeze["forest"]
        self.nf = len(self.forest)
        self.rider = name in FZ.RIDERS
        self.dt_us = self.freeze["dt_s"] * 1e6
        snap = self.freeze["snapshot"]["path"]
        self.T1 = [snap[str(q)]["t1_us"] for q in self.path]
        self.T2 = [snap[str(q)]["t2_us"] for q in self.path]
        self.ro_snap = [(snap[str(q)]["p01"], snap[str(q)]["p10"]) for q in (self.path[0], self.path[2], self.path[1])]
        self.A_snap = rl.assignment_product(self.ro_snap)
        self.forest_p01 = [self.freeze["snapshot"]["forest"][str(f)].get("p01", 0.0) for f in self.forest]
        self.taus = list(self.freeze["taus_us"])
        self.ops = {nm: self.extract(c) for nm, c in self.circs.items() if not nm.startswith("RO_")}
        self._ptrue = {}

    def extract(self, tc):
        """(name, logical qubits, params, duration_us) on the path in program order; the forest's delays dropped."""
        lg = {p: i for i, p in enumerate(self.path)}
        d = self.freeze["durations_dt"]
        ops = []
        for inst in tc.data:
            name = inst.operation.name
            qs = [tc.find_bit(q).index for q in inst.qubits]
            if name in ("barrier", "measure"):
                continue
            if not all(q in lg for q in qs):
                if name == "delay" and qs[0] in self.forest:
                    continue
                raise ValueError(f"{tc.name}: {name} on {qs} leaves the path")
            prm = [float(p) for p in inst.operation.params]
            if name == "delay":
                if inst.operation.unit != "dt":
                    raise ValueError(f"{tc.name}: delay in {inst.operation.unit}")
                du = prm[0] * self.dt_us
            elif name == "rz":
                du = 0.0
            elif name in ("sx", "x"):
                du = d[name][str(qs[0])] * self.dt_us
            elif name == "cz":
                du = d["cz"][f"{min(qs)}-{max(qs)}"] * self.dt_us
            else:
                raise ValueError(name)
            ops.append((name, [lg[q] for q in qs], prm, du))
        return ops

    def p_true(self, nm, L=None, ops=None):
        """Outcome probabilities before readout, k = 4 q0 + 2 q2 + q1, under the letters L (role-keyed)."""
        key = (nm, None if not L else tuple(sorted((str(k), v) for k, v in L.items())))
        if ops is None and key in self._ptrue:
            return self._ptrue[key]
        rho = replay_fn(L)(ops if ops is not None else self.ops[nm], ROLE, self.T1, self.T2, L)
        p = np.clip(rl.probs_rho(rho), 0, None)
        p = p / p.sum()
        if ops is None:
            self._ptrue[key] = p
        return p

    def names(self, prefix):
        return sorted(n for n in self.circs if n.startswith(prefix))

    def variants(self, arm):
        return [f"{arm}{r}{k}" for r, k in FZ.TURNS]


# ── counts ────────────────────────────────────────────────────────────────

def compact(counts_halves, nf):
    """{name: [h1, h2]} of IBM count dicts (key bits: forest ... q1 q2 q0, clbit 0 rightmost) ->
    {name: (path counts (2, 8), forest ones (2, nf), shots (2,))}."""
    out = {}
    for nm, halves in counts_halves.items():
        pc = np.zeros((2, 8))
        fo = np.zeros((2, nf))
        n = np.zeros(2)
        for h, cnt in enumerate(halves):
            for key, c in cnt.items():
                key = key.replace(" ", "")
                if len(key) != 3 + nf:
                    raise ValueError(f"{nm}: count key {key!r} has {len(key)} bits, expected {3 + nf}")
                q0, q2, q1 = int(key[-1]), int(key[-2]), int(key[-3])
                pc[h, 4 * q0 + 2 * q2 + q1] += c
                for i in range(nf):
                    fo[h, i] += c * int(key[-4 - i])
                n[h] += c
        out[nm] = (pc, fo, n)
    return out


def readout_matrix(comp):
    A = np.zeros((8, 8))
    nro = np.zeros(8)
    for j in range(8):
        pc, _, n = comp[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"]
        A[:, j] = pc.sum(0) / n.sum()
        nro[j] = n.sum()
    return A, nro


def shot_cov(par0, sel, n):
    """Multinomial covariance of the parities sel = [(basis, mask)] at n shots per circuit from parities par0
    ({basis: 7-vector}); parities of different circuits are independent."""
    bi = np.array([BI[b] for b, _ in sel])
    mi = np.array([MI[m] for _, m in sel])
    P = np.array([par0[b] for b in BASES])               # (9, 7)
    x = XOR[mi[:, None], mi[None, :]]
    full = np.where(x < 0, 1.0, P[bi[:, None], np.maximum(x, 0)])
    S = (full - P[bi, mi][:, None] * P[bi, mi][None, :]) / n
    return S * (bi[:, None] == bi[None, :])


def ro_kappa(comp):
    """The RO arm's drift factor: over its eight circuits, the chi^2 of the difference of the two halves' outcome
    distributions against its multinomial variance, over its degrees of freedom; max(1, that)."""
    chi, dof = 0.0, 0
    for j in range(8):
        pc, _, n = comp[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"]
        p = pc.sum(0) / n.sum()
        k = p > 0
        if k.sum() < 2:
            continue
        chi += float(np.sum((pc[0][k] / n[0] - pc[1][k] / n[1]) ** 2 / (p[k] * (1 / n[0] + 1 / n[1]))))
        dof += int(k.sum()) - 1
    return max(1.0, chi / dof) if dof else 1.0


def ro_cov(ptrue, sel, A, nro, kappa=1.0):
    """The noise of the RO matrix carried into the baseline s_m . (A p): column j multinomial at nro[j] shots, times the
    RO arm's drift factor."""
    P = np.array([ptrue[b] for b, _ in sel])                     # (R, 8)
    Sm = np.array([SG[MI[m]] for _, m in sel])                   # (R, 8)
    S = np.zeros((len(sel), len(sel)))
    for j in range(8):
        Cj = (np.diag(A[:, j]) - np.outer(A[:, j], A[:, j])) / nro[j]
        G = Sm @ Cj @ Sm.T
        S += np.outer(P[:, j], P[:, j]) * G
    return kappa * S


class DriftCov:
    """kappa Sigma for one arm's readings (module docstring), block by block: Sigma has no covariance between two basis
    circuits, and every PUB runs at its own shuffled position in the job, so a drift moves the readings of one circuit
    half together and those of two circuits independently. In each basis block, Sigma is whitened and S_b is the mean
    outer product of the standing variants' whitened half differences (h1 - h2) / 2; kappa = max(1, sum of trace S_b /
    rank); K_b = S_b with every eigenvalue below kappa raised to kappa; kappa Sigma = Sigma_b^1/2 K_b Sigma_b^1/2 on
    each block and zero between blocks."""

    def __init__(self, Sig, deltas, blocks):
        n = Sig.shape[0]
        self.parts = []
        tr, rank = 0.0, 0
        Ss = []
        eig = [np.linalg.eigh(Sig[np.ix_(bl, bl)]) for bl in blocks]
        top = max(float(w.max()) for w, _ in eig)
        for bl, (w, U) in zip(blocks, eig):
            keep = w > 1e-9 * top
            E = U[:, keep] * np.sqrt(w[keep])
            M = (U[:, keep] / np.sqrt(w[keep])).T
            Z = np.array([M @ d[bl] for d in deltas])
            S = Z.T @ Z / len(Z)
            tr += float(np.trace(S))
            rank += int(keep.sum())
            Ss.append((bl, E, M, S))
        self.rank = rank
        self.kappa = max(1.0, tr / rank)
        self.n = n
        top = 0.0
        for bl, E, M, S in Ss:
            lam, V = np.linalg.eigh(S)
            top = max(top, float(lam.max()) if len(lam) else 0.0)
            self.parts.append((bl, E @ V, V.T @ M, np.maximum(lam, self.kappa)))
        self.largest = top

    def cov(self, scale=1.0):
        C = np.zeros((self.n, self.n))
        for bl, EV, MV, k in self.parts:
            C[np.ix_(bl, bl)] = scale * (EV * k) @ EV.T
        return C

    def winv(self, scale=1.0):
        W = np.zeros((self.n, self.n))
        for bl, EV, MV, k in self.parts:
            W[np.ix_(bl, bl)] = (MV.T / k) @ MV / scale
        return W


def blocks_of(sel):
    return [[i for i, (b, _) in enumerate(sel) if b == bb] for bb in BASES if any(b == bb for b, _ in sel)]


# ── the analysis ──────────────────────────────────────────────────────────

def chi2_points(rank):
    return {q: float(CHI2.ppf(q, rank)) for q in (0.95, 0.99, 0.999)}


def ruling(stat, rank, nulls, full):
    """The chi^2 points, each replaced by the simulated point where that exceeds it by more than 5 % (only when
    every variant of the statistic stands, the configuration the simulated flights drew)."""
    pts = chi2_points(rank)
    notes = []
    if nulls is not None and full and stat in nulls.get("points", {}):
        sim = nulls["points"][stat]
        if sim.get("rank") == rank:
            for q in list(pts):
                s = sim[str(q)]
                if s > SIM_RULE * pts[q]:
                    notes.append(f"{stat}: simulated {q} point {s:.2f} rules over chi^2 {pts[q]:.2f}")
                    pts[q] = s
    return pts, notes


def frame_stat(Ys, W):
    Yb = np.mean(Ys, axis=0)
    return float(sum((y - Yb) @ W @ (y - Yb) for y in Ys))


def class_of(T, pts):
    if T > pts[P_TURNS]:
        return "TURNS"
    if T <= pts[P_BLIND]:
        return "FRAME-BLIND"
    return "INCONCLUSIVE"


def arm_cov(par, par_h, standing, sel, n):
    """Sigma from the measured parities pooled over the standing variants, and the drift matrix from their halves."""
    pooled = {b: np.mean([par[f"{v}_{b[0]}{b[1]}"] for v in standing], axis=0) for b in BASES}
    Sig = shot_cov(pooled, sel, n)
    deltas = [np.array([(par_h[f"{v}_{b[0]}{b[1]}"][0][MI[m]] - par_h[f"{v}_{b[0]}{b[1]}"][1][MI[m]]) / 2
                        for b, m in sel]) for v in standing]
    return Sig, DriftCov(Sig, deltas, blocks_of(sel))


def frame_block(chip, arm, sel, par, standing, dc, nulls, tag, p_ok=True):
    """T5, T4, T6 and E of one arm (A or B) on the readings sel."""
    W = dc.winv()
    rank = dc.rank
    Y = {v: np.array([par[f"{v}_{b[0]}{b[1]}"][MI[m]] for b, m in sel]) for v in chip.variants(arm)}
    out, notes = {}, []
    five = [f"{arm}5{k}" for k in range(4)]
    st5 = [v for v in five if v in standing]
    if len(st5) >= 2:
        T = frame_stat([Y[v] for v in st5], W)
        pts, nt = ruling(f"T5{tag}", (len(st5) - 1) * rank, nulls, len(st5) == 4)
        out["T5"] = {"T": T, "df": (len(st5) - 1) * rank, "points": pts, "class": class_of(T, pts), "variants": st5}
        notes += nt
    else:
        out["T5"] = {"class": "VOID"}
    for role in ("4", "6"):
        v2 = f"{arm}{role}2"
        if v2 in standing and f"{arm}50" in standing:
            T = frame_stat([Y[f"{arm}50"], Y[v2]], W)
            pts, nt = ruling(f"T{role}{tag}", rank, nulls, True)
            out[f"T{role}"] = {"T": T, "df": rank, "points": pts, "class": class_of(T, pts)}
            notes += nt
        else:
            out[f"T{role}"] = {"class": "VOID"}
    if len(st5) == 4:
        E = Y[five[1]] + Y[five[3]] - Y[five[2]] - Y[five[0]]
        t = float(E @ dc.winv(4.0) @ E)
        pts, nt = ruling(f"E{tag}", rank, nulls, True)
        out["E"] = {"T": t, "df": rank, "points": pts, "present": t > pts[0.999] if p_ok else "VOID"}
        notes += nt
    return out, notes, Y


def q_fit(x, y, cov, taus):
    """The VQ fit of one qubit: (delta, theta0, C, g) by weighted least squares with the full covariance."""
    tt = np.array(taus)
    r_obs = np.concatenate([x, y])
    try:
        Lc = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        return {"delta": float("nan"), "theta0": float("nan"), "C": float("nan"), "g": float("nan"),
                "sigma_delta": float("nan"), "chi2": float("nan"), "converged": False, "on_bound": False,
                "failure": "covariance not positive definite"}

    def res(v):
        d, th, C, g = v
        e = C * np.exp(-g * tt)
        return solve_triangular(Lc, np.concatenate([e * np.cos(th + d * tt), e * np.sin(th + d * tt)]) - r_obs, lower=True)
    fit = least_squares(res, [0.0, math.atan2(y[0], x[0]), 1.0, 0.0], bounds=Q_BOUNDS)
    chi = float(fit.fun @ fit.fun)
    J = fit.jac
    try:
        cp = np.linalg.inv(J.T @ J)
        sig = math.sqrt(cp[0, 0]) * max(1.0, math.sqrt(chi / 4))
    except (np.linalg.LinAlgError, ValueError):
        sig = float("inf")
    on_bound = abs(fit.x[0]) >= Q_BOUNDS[1][0] * (1 - 1e-6)
    return {"delta": float(fit.x[0]), "theta0": float(fit.x[1]), "C": float(fit.x[2]), "g": float(fit.x[3]),
            "sigma_delta": sig, "chi2": chi, "converged": bool(fit.success), "on_bound": bool(on_bound)}


def vq(chip, comp, A, nro, q_void, nulls, k_ro=1.0):
    """Arm Q: the readings corrected by the inverse RO matrix, their covariance (shots carried through A^-1 and scaled
    by the qubit's drift weight, plus the RO arm's noise), and the fit per qubit."""
    Ai = np.linalg.inv(A)
    Cjs = [k_ro * (np.diag(A[:, j]) - np.outer(A[:, j], A[:, j])) / nro[j] for j in range(8)]
    masks = {0: (1, 0, 0), 1: (0, 0, 1), 2: (0, 1, 0)}
    keys = [(l, P, t) for l in range(3) for P in "XY" for t in range(len(chip.taus))]
    read, var_tot, xs, halves = {}, {}, {}, {}
    for l, P, t in keys:
        nm = f"Q_{l}{P}{t}"
        pc, _, n = comp[nm]
        s = Ai.T @ SG[MI[masks[l]]]                     # reading = s . p_measured
        pt = pc.sum(0) / n.sum()
        read[nm] = float(s @ pt)
        var_tot[nm] = float(s @ ((np.diag(pt) - np.outer(pt, pt)) / n.sum()) @ s)
        xs[nm] = Ai @ pt
        hh = []
        for h in range(2):
            ph = pc[h] / n[h]
            hh.append((float(s @ ph), float(s @ ((np.diag(ph) - np.outer(ph, ph)) / n[h]) @ s)))
        halves[nm] = hh
    terms = {nm: (h[0][0] - h[1][0]) ** 2 / (h[0][1] + h[1][1]) for nm, h in halves.items() if h[0][1] + h[1][1] > 0}
    kQ = max(1.0, float(np.mean(list(terms.values())))) if terms else 1.0
    out = {"kappa_Q": kQ}
    scale = (nulls or {}).get("delta_scale", {})
    # a phase that wanders from PUB to PUB (a detuning that changes between the shuffled PUBs) moves a Ramsey reading
    # by theta * tau perpendicular to the Bloch vector: d<X> = -<Y> d(theta), d<Y> = <X> d(theta). Per qubit one
    # parameter s^2 (the wander's variance per unit tau^2) is estimated from its six readings with tau > 0: the excess of
    # (half1 - half2)^2 over var1 + var2 is 2 s^2 tau^2 g^2 with g the perpendicular component (<Y> for an X reading,
    # <X> for a Y reading). The reading's variance is then kappa_Q times its shot variance plus s^2 tau^2 g^2 / 2.
    nt = len(chip.taus)
    for l in range(3):
        nms = [f"Q_{l}{P}{t}" for P in "XY" for t in range(nt)]
        num, den = 0.0, 0.0
        g2 = {}
        for P in "XY":
            for t in range(nt):
                nm = f"Q_{l}{P}{t}"
                g2[nm] = read[f"Q_{l}{'Y' if P == 'X' else 'X'}{t}"] ** 2
                if chip.taus[t] > 0:
                    h = halves[nm]
                    num += (h[0][0] - h[1][0]) ** 2 - (h[0][1] + h[1][1])
                    den += 2 * chip.taus[t] ** 2 * g2[nm]
        s2 = max(0.0, num / den) if den > 0 else 0.0
        w = kQ
        sl = SG[MI[masks[l]]] @ Ai
        C = np.diag([kQ * var_tot[n] + s2 * chip.taus[int(n[-1])] ** 2 * g2[n] / 2 for n in nms])
        for a, na in enumerate(nms):
            for b, nb in enumerate(nms):
                C[a, b] += sum(xs[na][j] * xs[nb][j] * (sl @ Cjs[j] @ sl) for j in range(8))
        x = np.array([read[f"Q_{l}X{t}"] for t in range(len(chip.taus))])
        y = np.array([read[f"Q_{l}Y{t}"] for t in range(len(chip.taus))])
        f = q_fit(x, y, C, chip.taus)
        f.update({"weight": w, "phase_wander_s2": s2, "readings_X": x.tolist(), "readings_Y": y.tolist()})
        sc = scale.get(str(l), 1.0)
        f["sigma_delta_ruling"] = f["sigma_delta"] * sc
        f["sigma_scale"] = sc
        void = (not f["converged"]) or f["on_bound"] or any(q_void.get(n) is False for n in nms)
        z = abs(f["delta"]) / f["sigma_delta_ruling"] if f["sigma_delta_ruling"] > 0 else float("inf")
        f["z"] = z
        f["verdict"] = "VOID" if void or not math.isfinite(z) else (
            "DETUNED" if z >= Q_DETUNED_Z else ("NOT DETUNED" if z < 2 else "INCONCLUSIVE"))
        out[str(l)] = f
    return out


def vs(chip, comp, A):
    """VS, the forest, on all A circuits. Under the null the variances are taken at the null's rates (p0 on the A side,
    (ones + 1/2) / (n + 1) on each RO circuit), and both sides are scaled by the forest qubit's drift factor from the
    two halves of every A and RO circuit."""
    a_names = sorted(n for n in comp if n.startswith("A"))
    PA = sum(comp[n][0].sum(0) for n in a_names)
    NA = PA.sum()
    PA = PA / NA
    out = {}
    for i, f in enumerate(chip.forest):
        ones = sum(comp[n][1][:, i].sum() for n in a_names)
        pA = ones / NA
        r = np.zeros(8)
        rt = np.zeros(8)
        nr = np.zeros(8)
        ro_ones = 0.0
        for j in range(8):
            _, fo, n = comp[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"]
            r[j] = fo[:, i].sum() / n.sum()
            rt[j] = (fo[:, i].sum() + 0.5) / (n.sum() + 1)
            nr[j] = n.sum()
            ro_ones += fo[:, i].sum()
        p0 = float(PA @ r)
        p0t = float(PA @ rt)
        terms = []
        for nm in a_names + [f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}" for j in range(8)]:
            _, fo, n = comp[nm]
            q = fo[:, i].sum() / n.sum()
            if 0 < q < 1:
                terms.append((fo[0, i] / n[0] - fo[1, i] / n[1]) ** 2 / (q * (1 - q) * (1 / n[0] + 1 / n[1])))
        kf = max(1.0, float(np.mean(terms))) if terms else 1.0
        var0 = p0t * (1 - p0t) / NA + float(np.sum(PA ** 2 * rt * (1 - rt) / nr))
        sig = math.sqrt(kf * var0)
        if ones == 0 and ro_ones == 0:
            out[str(f)] = {"p_A": 0.0, "p0": 0.0, "sigma": 0.0, "z": float("nan"), "kappa_f": kf,
                           "verdict": "no verdict (no flip on either side)"}
            continue
        z = (pA - p0) / sig
        v = "HEARD" if z >= 5 else ("SILENT" if abs(z) < 3 else "no verdict")
        out[str(f)] = {"p_A": pA, "p0": p0, "sigma": sig, "z": z, "kappa_f": kf, "verdict": v}
    return out


def analyze_compact(chip, comp, nulls=None, turn_ok=None, quiet=True, lines=None, readings=False):
    """Every verdict of both pages for one chip, from compact counts. turn_ok: {variant or Q circuit: True/False/None}
    from the circuits retrieved after the flight (False makes it VOID; None means not checked)."""
    turn_ok = turn_ok or {}
    lines = lines if lines is not None else []
    A, nro = readout_matrix(comp)
    par, par_h, nshots = {}, {}, {}
    for nm, (pc, fo, n) in comp.items():
        if nm.startswith("RO_"):
            continue
        par[nm] = SG @ (pc.sum(0) / n.sum())
        par_h[nm] = [SG @ (pc[h] / n[h]) for h in range(2)]
        nshots[nm] = n.sum()
    k_ro = ro_kappa(comp)
    out = {"readout_matrix": A.tolist(), "notes": [], "kappa_RO": k_ro}
    # the RO arm's own P(1|0), P(0|1) and readout error per path qubit, marginal over the other bits, against section 6
    rog = {}
    for name, bit, q in (("q0", 4, chip.path[0]), ("q2", 2, chip.path[2]), ("q1", 1, chip.path[1])):
        p10 = float(np.mean([sum(A[j, k] for j in range(8) if j & bit) for k in range(8) if not k & bit]))
        p01 = float(np.mean([sum(A[j, k] for j in range(8) if not j & bit) for k in range(8) if k & bit]))
        rog[str(q)] = {"P(1|0)": p10, "P(0|1)": p01,
                       "breaches_guard": p10 > RO_GUARD_P1G0 or (p10 + p01) / 2 > RO_GUARD_ERR}
    out["RO_guard"] = rog
    # ---- the drive-frame page ----
    va = chip.variants("A")
    standing_a = [v for v in va if turn_ok.get(v) is not False]
    out["A_void"] = [v for v in va if v not in standing_a]
    pt_a = {b: chip.p_true(f"A50_{b[0]}{b[1]}") for b in BASES}
    y0a = {b: SG @ (A @ pt_a[b]) for b in BASES}
    Sig_a, dc_a = arm_cov(par, par_h, standing_a, A_IDX, nshots["A50_XX"])
    Ro_a = ro_cov(pt_a, A_IDX, A, nro, k_ro)
    xz, zx = par["A50_XZ"][MI[(1, 1, 0)]], par["A50_ZX"][MI[(1, 1, 0)]]
    p_ok = bool((xz < -0.5) and (-0.75 < zx < -0.3) and (xz < zx - P_GATE_MARGIN))
    out["P"] = {"pass": p_ok, "<X0Z2>": xz, "<Z0X2>": zx,
                "turned variants (reading)": {v: {"<X0Z2>": par[f"{v}_XZ"][MI[(1, 1, 0)]], "<Z0X2>": par[f"{v}_ZX"][MI[(1, 1, 0)]]}
                                              for v in va[1:]}}
    kS = dc_a.cov()
    S_meas = par["A50_YZ"][MI[(1, 1, 0)]]
    S_base = y0a[("Y", "Z")][MI[(1, 1, 0)]]
    sig_S = math.sqrt(kS[S_POS, S_POS] + Ro_a[S_POS, S_POS])
    s_A = S_meas - S_base
    thr = 3.0
    if nulls is not None and "sA_997" in nulls.get("points", {}):
        sim = nulls["points"]["sA_997"]
        if sim > SIM_RULE * 3.0:
            thr = sim
            out["notes"].append(f"V1: simulated 99.7 % point of |s_A|/sigma_S {sim:.2f} rules over 3")
    if not p_ok or "A50" not in standing_a:
        v1 = "VOID"
    elif s_A >= thr * sig_S:
        v1 = "REPLICATED"
    elif s_A <= -thr * sig_S:
        v1 = "REVERSED"
    else:
        v1 = "ABSENT" if abs(s_A) < 2 * sig_S else "INCONCLUSIVE"
    out["V1"] = {"verdict": v1, "s_A": s_A, "sigma_S": sig_S, "S": S_meas, "S_base": S_base, "z": s_A / sig_S}
    fa, nt, Ya = frame_block(chip, "A", A_IDX, par, standing_a, dc_a, nulls, "", p_ok)
    out["notes"] += nt
    out["VF"] = fa
    out["drift_A"] = {"kappa": dc_a.kappa, "largest_direction": dc_a.largest,
                      "kappa_i": [float(x) for x in np.diag(dc_a.cov()) / np.diag(Sig_a)]}
    five = [f"A5{k}" for k in range(4)]
    if all(v in standing_a for v in five):
        Ss = [Ya[v][S_POS] for v in five]
        Sb = float(np.mean(Ss))
        out["share"] = {"S_bar_minus_base": Sb - S_base, "sigma": math.sqrt(kS[S_POS, S_POS] / 4 + Ro_a[S_POS, S_POS]),
                        "S0_minus_S_bar": Ss[0] - Sb, "sigma_t": math.sqrt(kS[S_POS, S_POS]) * math.sqrt(3) / 2}
    out["VS"] = vs(chip, comp, A)
    out["VS_said_with"] = ("the correlation between p_A and P_A(s) from shared shots is neglected, and the forest's "
                           "longer idle in the A circuits is not modelled")
    if not p_ok:                       # a chip that fails the parse gate is VOID on this page; its numbers are kept
        out["notes"].append("P fails: every verdict of the drive-frame page on this chip is VOID")
        for k in ("T5", "T4", "T6"):
            out["VF"][k]["class"] = "VOID"
        for x in out["VS"].values():
            x["verdict"] = "VOID"
    # sightings: single readings of A5(0) beyond 4 sigma against the baseline (a sighting, not a finding)
    za = (np.array([par[f"A50_{b[0]}{b[1]}"][MI[m]] for b, m in A_IDX]) - np.array([y0a[b][MI[m]] for b, m in A_IDX])) \
        / np.sqrt(np.diag(kS) + np.diag(Ro_a))
    out["A50_sightings"] = [(f"{b[0]}{b[1]} {m}", float(z)) for (b, m), z in zip(A_IDX, za) if abs(z) > 4]
    # the even readings of every turned variant against A5(0), with the shot sigma of a difference of two circuits
    ev = {}
    for v in va[1:]:
        worst = 0.0
        for b in BASES:
            for m in rl.MASKS:
                if rl.odd_y(b, m):
                    continue
                pa, pv = par[f"A50_{b[0]}{b[1]}"][MI[m]], par[f"{v}_{b[0]}{b[1]}"][MI[m]]
                sd = math.sqrt(max((1 - pa ** 2) / nshots[f"A50_{b[0]}{b[1]}"] + (1 - pv ** 2) / nshots[f"{v}_{b[0]}{b[1]}"], 1e-12))
                worst = max(worst, abs(pv - pa) / sd)
        ev[v] = worst
    out["turned_even_largest_z"] = ev
    if chip.rider:
        out.update(analyze_b(chip, comp, A, nro, par, par_h, nshots, turn_ok, nulls, k_ro))
        if readings:
            standing_b = [v for v in chip.variants("B") if v not in out["B_void"]]
            out["colouring"] = colouring(chip, comp, A, nro, standing_b, out["VQ"])
    if not quiet:
        report(chip, out, lines)
    return out


def analyze_b(chip, comp, A, nro, par, par_h, nshots, turn_ok, nulls, k_ro=1.0):
    out = {"notes_B": []}
    vb = chip.variants("B")
    standing = [v for v in vb if turn_ok.get(v) is not False]
    out["B_void"] = [v for v in vb if v not in standing]
    pt = {b: chip.p_true(f"B50_{b[0]}{b[1]}") for b in BASES}
    y0 = {b: SG @ (A @ pt[b]) for b in BASES}
    n_b = nshots["B50_XX"]
    Sig, dc = arm_cov(par, par_h, standing, B_IDX, n_b)
    Ro = ro_cov(pt, B_IDX, A, nro, k_ro)
    out["drift_B"] = {"kappa": dc.kappa, "largest_direction": dc.largest}
    xx = par["B50_XX"]
    pb = (xx[MI[(1, 0, 0)]] > PB_GATE[0]) and (xx[MI[(0, 1, 0)]] > PB_GATE[1]) and (xx[MI[(1, 1, 1)]] < PB_GATE[2])
    out["PB"] = {"pass": bool(pb), "<X0>": xx[MI[(1, 0, 0)]], "<X2>": xx[MI[(0, 1, 0)]], "<X0X1X2>": xx[MI[(1, 1, 1)]],
                 "<X0X1X2> baseline": y0[("X", "X")][MI[(1, 1, 1)]]}
    kS = dc.cov()
    y = np.array([par[f"B50_{b[0]}{b[1]}"][MI[m]] for b, m in B_IDX])
    yb = np.array([y0[b][MI[m]] for b, m in B_IDX])
    r = y - yb
    Wt, rk = pinv_rank(kS + Ro)
    if pb:
        Wv = float(r @ Wt @ r)
        pts, nt = ruling("W", rk, nulls, True)
        out["notes_B"] += nt
        wb = "OPEN" if Wv > pts[P_OPEN] else ("SHUT" if Wv < pts[P_SHUT] else "INCONCLUSIVE")
        if "B50" not in standing:
            wb = "VOID"
        kinds = {}
        for kname, ii in B_KINDS.items():
            sub = np.ix_(ii, ii)
            Wk, rkk = pinv_rank((kS + Ro)[sub])
            kinds[kname] = {"W": float(r[ii] @ Wk @ r[ii]), "df": rkk}
        pooled = {}
        for pname, ii in POOLED.items():
            C = (kS + Ro)[np.ix_(ii, ii)]
            pooled[pname] = {"value": float(np.mean(y[ii])), "baseline": float(np.mean(yb[ii])),
                             "sigma": math.sqrt(float(np.sum(C))) / 3}
        # the other 22 forbidden readings, against the baseline, with their own kappa Sigma + Sigma_RO
        _, dco = arm_cov(par, par_h, standing, B_OTHER, n_b)
        Ro_o = ro_cov(pt, B_OTHER, A, nro, k_ro)
        ro_ = np.array([par[f"B50_{b[0]}{b[1]}"][MI[m]] - y0[b][MI[m]] for b, m in B_OTHER])
        Wo, rko = pinv_rank(dco.cov() + Ro_o)
        zo = ro_ / np.sqrt(np.diag(dco.cov() + Ro_o))
        out["B_other22_z"] = [(f"{b[0]}{b[1]} {m}", float(z)) for (b, m), z in zip(B_OTHER, zo)]
        out["B_sightings"] = [(f"{b[0]}{b[1]} {m}", float(z)) for (b, m), z in zip(B_IDX + B_OTHER, list(r / np.sqrt(np.diag(kS + Ro))) + list(zo)) if abs(z) > 4]
        out["WB"] = {"verdict": wb, "W": Wv, "df": rk, "points": pts, "kinds": kinds, "pooled": pooled,
                     "other22": {"W": float(ro_ @ Wo @ ro_), "df": rko}, "residuals": r.tolist(),
                     "sigmas": np.sqrt(np.diag(kS + Ro)).tolist()}
        fb, nt, Yb = frame_block(chip, "B", B_IDX, par, standing, dc, nulls, "B")
        out["notes_B"] += nt
        out["VFB"] = fb
        five = [f"B5{k}" for k in range(4)]
        if all(v in standing for v in five):
            mean = np.mean([Yb[v] for v in five], axis=0) - yb
            Wm, rkm = pinv_rank(dc.cov(0.25) + Ro)
            out["B_turn_mean"] = {"chi2": float(mean @ Wm @ mean), "df": rkm, "largest": float(np.abs(mean).max())}
    else:
        out["WB"] = {"verdict": "no B verdicts (PB fails)"}
        out["VFB"] = {"T5": {"class": "VOID"}, "T4": {"class": "VOID"}, "T6": {"class": "VOID"}}
    q_void = {n: turn_ok.get(n) for n in chip.names("Q_")}
    out["VQ"] = vq(chip, comp, A, nro, q_void, nulls, k_ro)
    # the replay of B5(0) with the three fitted delta as its only letters, a reading beside W
    dl = {("det", ROLE[l]): out["VQ"][str(l)]["delta"] for l in range(3) if math.isfinite(out["VQ"][str(l)]["delta"])}
    if dl:
        yd = np.array([SG[MI[m]] @ (A @ chip.p_true(f"B50_{b[0]}{b[1]}", dl)) for b, m in B_IDX])
        out["B_delta_replay"] = {"W_of_replay_against_y0": float((yd - yb) @ Wt @ (yd - yb)),
                                 "W_of_measurement_against_replay": float((y - yd) @ Wt @ (y - yd)), "df": rk,
                                 "letters": {str(k): v for k, v in dl.items()}}
    fp = {}
    for i, f in enumerate(chip.forest):
        for arm in ("B", "Q"):
            nms = [n for n in comp if n.startswith(arm)]
            fp.setdefault(str(f), {})[arm] = float(sum(comp[n][1][:, i].sum() for n in nms) / sum(comp[n][2].sum() for n in nms))
    out["forest_P1_BQ"] = fp
    return out


XX_MASKS = {"X0": (1, 0, 0), "X1": (0, 0, 1), "X2": (0, 1, 0), "X0X1": (1, 0, 1), "X1X2": (0, 1, 1), "X0X2": (1, 1, 0),
            "X0X1X2": (1, 1, 1)}


def _xx_readings(p):
    """Singles, pairs, triple of the XX basis from outcome probabilities; the connected pairs and the third cumulant."""
    par = SG @ p
    r = {k: float(par[MI[m]]) for k, m in XX_MASKS.items()}
    a, b, c = r["X0"], r["X1"], r["X2"]
    r["c01"], r["c12"], r["c02"] = r["X0X1"] - a * b, r["X1X2"] - b * c, r["X0X2"] - a * c
    r["k3"] = r["X0X1X2"] - r["X0X1"] * c - r["X0X2"] * b - r["X1X2"] * a + 2 * a * b * c
    r["L3"] = -math.log(abs(r["X0X1X2"])) if r["X0X1X2"] else float("inf")
    r["L1"] = -sum(math.log(abs(r[k])) for k in ("X0", "X1", "X2")) if a * b * c else float("inf")
    return r


def colouring(chip, comp, A, nro, standing, vq_out, n_boot=400, seed=20261113):
    """X^3, the colouring of YZ+ZY, read in situ (X^3 window page section 5, a reading without a verdict). The XX-basis
    circuits of B's standing variants, which are logically one circuit, pooled and corrected by the inverse RO matrix;
    beside them the RO-free replay under the freeze snapshot and under arm Q's own coherence rates (T2 = 1/g, held at
    or below 2 T1). In the logical evolution X^3 commutes with H and decays at sum_l 1/T2_l (F158's far element; the
    absorption theorem's 2 gamma per X letter plus 1/(2 T1)), whatever H is; the compiled circuit leaves that frame
    between its cz gates, so the reading is a scale against the replay and not a model-free number. sigma by resampling
    the pooled counts and the RO arm's columns."""
    names = [f"{v}_XX" for v in standing]
    pc = sum(comp[n][0].sum(0) for n in names)
    N = pc.sum()
    pm = pc / N
    meas = _xx_readings(np.linalg.inv(A) @ pm)
    base = _xx_readings(chip.p_true("B50_XX"))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        pb = rng.multinomial(int(N), pm) / N
        Ab = np.stack([rng.multinomial(int(nro[j]), A[:, j] / A[:, j].sum()) / nro[j] for j in range(8)], axis=1)
        try:
            boots.append(_xx_readings(np.linalg.inv(Ab) @ pb))
        except (np.linalg.LinAlgError, ValueError):
            continue
    sd = {k: float(np.std([b[k] for b in boots])) for k in meas}
    out = {"circuits": names, "shots": float(N), "measured": meas, "sigma": sd, "replay_snapshot": base,
           "scale_L3": meas["L3"] / base["L3"], "scale_L3_sigma": sd["L3"] / base["L3"],
           "scale_L1": meas["L1"] / base["L1"], "scale_L1_sigma": sd["L1"] / base["L1"]}
    gs = [vq_out[str(l)] for l in range(3)]
    if all(g["verdict"] != "VOID" and math.isfinite(g["g"]) for g in gs):
        T2q = [1.0 / max(g["g"], 1.0 / (2 * chip.T1[l])) for l, g in enumerate(gs)]
        rho = rl.replay(chip.ops["B50_XX"], ROLE, chip.T1, T2q)
        out["replay_with_Q"] = _xx_readings(np.clip(rl.probs_rho(rho), 0, None))
        out["T2_from_Q_us"] = T2q
    return out


def report(chip, out, lines):
    s = lambda x: say(lines, x)
    s(f"== {chip.name} path {chip.path} forest {chip.forest}; RO drift factor {out['kappa_RO']:.2f}")
    for n in out["notes"] + out.get("notes_B", []):
        s(f"   note: {n}")
    s("   RO arm, per path qubit (reading): " + "; ".join(
        f"q{q} P(1|0) {x['P(1|0)']:.4f} P(0|1) {x['P(0|1)']:.4f}" + (" BREACHES the P(1|0) guard" if x["breaches_guard"] else "")
        for q, x in out["RO_guard"].items()))
    s(f"P  parse gate: {'PASS' if out['P']['pass'] else 'FAIL'}  <X0Z2> {out['P']['<X0Z2>']:+.3f}  <Z0X2> {out['P']['<Z0X2>']:+.3f}")
    s("   turned variants (reading): " + "; ".join(f"{v} {x['<X0Z2>']:+.3f}/{x['<Z0X2>']:+.3f}"
                                                   for v, x in out["P"]["turned variants (reading)"].items()))
    v = out["V1"]
    s(f"V1 replication: {v['verdict']}  s_A {v['s_A']:+.4f} (S {v['S']:+.4f}, baseline {v['S_base']:+.4f}), sigma_S {v['sigma_S']:.4f}, z {v['z']:+.2f}")
    s(f"   drift A: kappa {out['drift_A']['kappa']:.3f}, largest whitened direction {out['drift_A']['largest_direction']:.2f}; "
      f"void variants {out['A_void']}")
    for k in ("T5", "T4", "T6"):
        x = out["VF"][k]
        if "T" not in x:
            s(f"VF {k}: {x['class']}")
        else:
            s(f"VF {k}: {x['class']}  T {x['T']:.1f} on {x['df']}  (99 % {x['points'][0.99]:.1f}, 99.9 % {x['points'][0.999]:.1f})")
    if "E" in out["VF"]:
        e = out["VF"]["E"]
        pres = "VOID" if e["present"] == "VOID" else ("2phi PRESENT" if e["present"] else "not present")
        s(f"   E: {e['T']:.1f} on {e['df']} (99.9 % {e['points'][0.999]:.1f}) -> {pres}")
    if "share" in out:
        sh = out["share"]
        s(f"   share (reading): S_bar - base {sh['S_bar_minus_base']:+.4f} +- {sh['sigma']:.4f}; S(0) - S_bar {sh['S0_minus_S_bar']:+.4f} +- {sh['sigma_t']:.4f}")
    s(f"   A5(0) readings beyond 4 sigma (sightings): {out['A50_sightings']}")
    s("   turned variants' even readings against A5(0), largest |z| (reading): "
      + "; ".join(f"{v} {z:.1f}" for v, z in out["turned_even_largest_z"].items()))
    for f, x in out["VS"].items():
        s(f"VS q{f}: {x['verdict']}  p_A {x['p_A']:.5f}  p0 {x['p0']:.5f}  sigma {x['sigma']:.5f}  z {x['z']:+.2f}  kappa_f {x['kappa_f']:.2f}")
    s(f"   VS said with: {out['VS_said_with']}")
    if chip.rider:
        pb = out["PB"]
        s(f"PB parse gate: {'PASS' if pb['pass'] else 'FAIL'}  <X0> {pb['<X0>']:+.3f} <X2> {pb['<X2>']:+.3f} <X0X1X2> {pb['<X0X1X2>']:+.3f} (baseline {pb['<X0X1X2> baseline']:+.3f})")
        wb = out["WB"]
        if "W" in wb:
            s(f"WB window: {wb['verdict']}  W {wb['W']:.1f} on {wb['df']} (95 % {wb['points'][0.95]:.1f}, 99.9 % {wb['points'][0.999]:.1f})")
            s("   kinds: " + "; ".join(f"{k} W {x['W']:.1f} on {x['df']}" for k, x in wb["kinds"].items())
              + f"; other 22: W {wb['other22']['W']:.1f} on {wb['other22']['df']}")
            s(f"   B readings beyond 4 sigma (sightings): {out.get('B_sightings')}")
            s("   pooled: " + "; ".join(f"<{k}> {x['value']:+.4f} (base {x['baseline']:+.4f}) +- {x['sigma']:.4f}" for k, x in wb["pooled"].items()))
        else:
            s(f"WB: {wb['verdict']}")
        s(f"   drift B: kappa {out['drift_B']['kappa']:.3f}, largest whitened direction {out['drift_B']['largest_direction']:.2f}")
        for k in ("T5", "T4", "T6"):
            x = out["VFB"][k]
            if "T" not in x:
                s(f"VFB {k}: {x['class']}")
            else:
                s(f"VFB {k}: {x['class']}  T {x['T']:.1f} on {x['df']}  (99 % {x['points'][0.99]:.1f}, 99.9 % {x['points'][0.999]:.1f})")
        if "E" in out["VFB"]:
            e = out["VFB"]["E"]
            s(f"   E_B: {e['T']:.1f} on {e['df']} -> {'2phi PRESENT' if e['present'] else 'not present'}")
        if "B_turn_mean" in out:
            m = out["B_turn_mean"]
            s(f"   B turn mean against y0 (reading): chi2 {m['chi2']:.1f} on {m['df']}, largest {m['largest']:.4f}")
        for l in range(3):
            q = out["VQ"][str(l)]
            s(f"VQ q{chip.path[l]}: {q['verdict']}  delta {q['delta']:+.4f} rad/us +- {q['sigma_delta_ruling']:.4f} (z {q['z']:.2f}), "
              f"chi2 {q['chi2']:.1f}, C {q['C']:.3f}, g {q['g']:.3f}/us, kappa_Q {q['weight']:.2f}, phase wander s {math.sqrt(q['phase_wander_s2']):.4f} rad/us")
        if "B_delta_replay" in out:
            d = out["B_delta_replay"]
            s(f"   B5(0) replayed with the three delta (reading): W against y0 {d['W_of_replay_against_y0']:.1f}, "
              f"measurement against it {d['W_of_measurement_against_replay']:.1f} on {d['df']}")
        if "colouring" in out:
            c = out["colouring"]
            m, b, sd = c["measured"], c["replay_snapshot"], c["sigma"]
            s(f"   the colouring X^3, RO-corrected, {len(c['circuits'])} XX circuits pooled (reading): <X0X1X2> {m['X0X1X2']:+.4f} "
              f"+- {sd['X0X1X2']:.4f} (replay {b['X0X1X2']:+.4f}); -ln|XXX| {m['L3']:.4f} +- {sd['L3']:.4f} against "
              f"{b['L3']:.4f}, scale {c['scale_L3']:.2f} +- {c['scale_L3_sigma']:.2f}; singles' sum {m['L1']:.4f} against "
              f"{b['L1']:.4f}, scale {c['scale_L1']:.2f} +- {c['scale_L1_sigma']:.2f}; third cumulant {m['k3']:+.4f} +- {sd['k3']:.4f} "
              f"(replay {b['k3']:+.4f})"
              + (f"; the replay with arm Q's T2 {[round(x) for x in c['T2_from_Q_us']]} us: -ln|XXX| {c['replay_with_Q']['L3']:.4f}"
                 if "replay_with_Q" in c else ""))
        s("   forest P(1) in B and Q (reading): " + "; ".join(f"q{f} B {x['B']:.4f} Q {x['Q']:.4f}" for f, x in out["forest_P1_BQ"].items()))


def run_window(data, chip):
    """Hours from the freeze snapshot's oldest property to the job's start, when the job's metrics hold it."""
    try:
        ts = data.get("job_metrics", {}).get("timestamps", {})
        start = ts.get("running") or ts.get("created")
        if not start:
            return None
        t = datetime.fromisoformat(str(start).replace("Z", "+00:00"))
        old = datetime.fromisoformat(chip.freeze["oldest_property_utc"])
        return (t - old).total_seconds() / 3600
    except Exception:
        return None


def analyze(result_path, reach=False):
    data = json.loads(Path(result_path).read_text(encoding="utf-8"))
    chip = Chip(data["backend"])
    if data.get("qpy_sha256") != sha256(chip.qpy):
        raise SystemExit(f"result was taken with circuits {data.get('qpy_sha256')}, this file is {sha256(chip.qpy)}")
    if data.get("freeze_sha256", chip.freeze_sha) != chip.freeze_sha or \
            ("freeze_snapshot" in data and data["freeze_snapshot"] != chip.freeze["snapshot"]):
        raise SystemExit("the record was taken under another freeze record than the one on disk; check out the commit "
                         "the flight was submitted from")
    lines = []
    nulls = load_nulls(chip.name, strict=False)
    if nulls is None:
        say(lines, "   note: no simulated-flight nulls on disk for this chip; the chi^2 points rule alone")
    comp = compact(data["counts"], chip.nf)
    tc = data.get("turn_ok")
    out = analyze_compact(chip, comp, nulls=nulls, turn_ok=tc if isinstance(tc, dict) else None, quiet=False, lines=lines,
                          readings=True)
    if not isinstance(tc, dict) or not tc:
        say(lines, "   note: the circuits were not checked against the ones retrieved from the job; the variants "
                   "stand on the structure asserted at submission")
    else:
        unchecked = sorted(k for k, v in tc.items() if v is None)
        if unchecked:
            say(lines, f"   note: not checked on retrieved circuits: {unchecked}")
    h = run_window(data, chip)
    out["hours_from_oldest_property_to_run"] = h
    if h is not None:
        say(lines, f"   the job ran {h:.1f} h after the freeze snapshot's oldest property (a reading)")
    if reach:
        out["reach"] = reach_day(chip, comp, out, lines)
    stem = Path(result_path).with_suffix("")
    Path(str(stem) + ".analysis.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    Path(str(stem) + ".analysis.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return out


# ── simulated flights ─────────────────────────────────────────────────────

def true_probs(chip, L=None):
    """Each circuit's outcome distribution after readout under the snapshot's RO matrix, from its own replay."""
    P = {}
    for nm in chip.circs:
        if nm.startswith("RO_"):
            j = 4 * int(nm[3]) + 2 * int(nm[4]) + int(nm[5])
            P[nm] = chip.A_snap[:, j].copy()
        else:
            P[nm] = chip.A_snap @ chip.p_true(nm, L)
    return P


def draw(chip, P, rng, forest_p=None):
    """One simulated job: each circuit's two halves multinomial from P (a distribution, or a pair for the two halves),
    the forest bits drawn independently at forest_p (default the snapshot's P(1|0); a callable of the circuit name
    returning one list per half)."""
    fp = chip.forest_p01 if forest_p is None else forest_p
    comp = {}
    for nm, p in P.items():
        ps = p if isinstance(p, tuple) else (p, p)
        pc = np.array([rng.multinomial(SHOTS_HALF, q / q.sum()) for q in ps], dtype=float)
        fpp = fp(nm) if callable(fp) else [fp, fp]
        fo = np.array([[rng.binomial(SHOTS_HALF, q) for q in fpp[h]] for h in range(2)], dtype=float).reshape(2, chip.nf)
        comp[nm] = (pc, fo, np.array([SHOTS_HALF, SHOTS_HALF], dtype=float))
    return comp


def to_counts(chip, comp, rng):
    """A compact draw as IBM count dicts (the forest bits drawn independently of the path bits)."""
    out = {}
    for nm, (pc, fo, n) in comp.items():
        halves = []
        for h in range(2):
            path_k = np.repeat(np.arange(8), pc[h].astype(int))
            rng.shuffle(path_k)
            bits = np.zeros((len(path_k), chip.nf), dtype=int)
            for i in range(chip.nf):
                col = np.zeros(len(path_k), dtype=int)
                col[rng.choice(len(path_k), int(fo[h, i]), replace=False)] = 1
                bits[:, i] = col
            cnt = {}
            for k, fb in zip(path_k, bits):
                q0, q2, q1 = k >> 2 & 1, k >> 1 & 1, k & 1
                key = "".join(str(x) for x in fb[::-1]) + f"{q1}{q2}{q0}"
                cnt[key] = cnt.get(key, 0) + 1
            halves.append(cnt)
        out[nm] = halves
    return out


NULL_STATS = ("sA", "T5", "T4", "T6", "E", "W", "T5B", "T4B", "T6B", "EB")


def stats_of(out, rider):
    r = {"sA": out["V1"]["s_A"] / out["V1"]["sigma_S"], "T5": out["VF"]["T5"]["T"], "T4": out["VF"]["T4"]["T"],
         "T6": out["VF"]["T6"]["T"], "E": out["VF"]["E"]["T"],
         "ranks": {"T5": out["VF"]["T5"]["df"], "T4": out["VF"]["T4"]["df"], "T6": out["VF"]["T6"]["df"],
                   "E": out["VF"]["E"]["df"]},
         "P": out["P"]["pass"]}
    for f, x in out["VS"].items():
        r[f"VS{f}"] = x["z"]
    if rider and "W" in out["WB"]:
        r.update({"W": out["WB"]["W"], "T5B": out["VFB"]["T5"]["T"], "T4B": out["VFB"]["T4"]["T"],
                  "T6B": out["VFB"]["T6"]["T"], "EB": out["VFB"]["E"]["T"]})
        r["ranks"].update({"W": out["WB"]["df"], "T5B": out["VFB"]["T5"]["df"], "T4B": out["VFB"]["T4"]["df"],
                           "T6B": out["VFB"]["T6"]["df"], "EB": out["VFB"]["E"]["df"]})
    if rider:
        r["PB"] = out["PB"]["pass"]
        for l in range(3):
            q = out["VQ"][str(l)]
            r[f"delta{l}"] = q["delta"]
            r[f"sdelta{l}"] = q["sigma_delta"]
    return r


def _null_worker(args):
    chip_name, seeds = args
    chip = Chip(chip_name)
    P = true_probs(chip)
    res = []
    for sd in seeds:
        out = analyze_compact(chip, draw(chip, P, np.random.default_rng(sd)))
        res.append(stats_of(out, chip.rider))
    return res


def nulls_path(chip_name):
    return HERE / f"drive_frame_{chip_name}_nulls.json"


def load_nulls(chip_name, strict=True):
    """The ruling simulated points: drawn for this frozen file and this freeze record, from N_NULLS flights."""
    p = nulls_path(chip_name)
    if not p.exists():
        if strict:
            raise SystemExit(f"{p.name} is missing")
        return None
    nl = json.loads(p.read_text(encoding="utf-8"))
    qpy, meta = FZ.files(chip_name)
    fr = json.loads(meta.read_text(encoding="utf-8"))
    bad = []
    if nl.get("qpy_sha256") != fr["sha256"]:
        bad.append(f"drawn for circuits {nl.get('qpy_sha256')}, the frozen file is {fr['sha256']}")
    if nl.get("freeze_sha256") != FZ.canonical_sha(fr):
        bad.append("drawn under another freeze record")
    if nl.get("n") != N_NULLS:
        bad.append(f"{nl.get('n')} flights, not {N_NULLS}")
    if bad:
        raise SystemExit(f"{p.name}: " + "; ".join(bad))
    return nl


def nulls(chip_name, n=N_NULLS, workers=16):
    from concurrent.futures import ProcessPoolExecutor
    chip = Chip(chip_name)
    seeds = list(range(1, n + 1))
    chunks = [(chip_name, seeds[i::workers]) for i in range(workers)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = [r for part in ex.map(_null_worker, chunks) for r in part]
    lines = []
    say(lines, f"{chip_name}: {len(res)} simulated flights drawn from the baseline replay ({time.time() - t0:.0f} s); "
               f"path {chip.path}, circuits {chip.freeze['sha256']}")
    say(lines, f"  parse gates: P fails in {1 - np.mean([r['P'] for r in res]):.4f}"
        + (f", PB in {1 - np.mean([r['PB'] for r in res]):.4f}" if chip.rider else ""))
    points = {}
    for st in NULL_STATS:
        if st not in res[0]:
            continue
        x = np.array([r[st] for r in res])
        if st == "sA":
            q = float(np.quantile(np.abs(x), 0.997))
            points["sA_997"] = q
            say(lines, f"  |s_A|/sigma_S: 99.7 % point {q:.3f} (against 3); mean {x.mean():+.3f}, sd {x.std():.3f}")
            continue
        rank = res[0]["ranks"][st]
        cp = chi2_points(rank)
        sim = {str(q): float(np.quantile(x, q)) for q in (0.95, 0.99, 0.999)}
        sim["rank"] = rank
        points[st] = sim
        say(lines, f"  {st:4s} rank {rank:3d}: mean {x.mean():7.2f}; simulated 95/99/99.9 % {sim['0.95']:.2f} / {sim['0.99']:.2f} / "
                   f"{sim['0.999']:.2f} against chi^2 {cp[0.95]:.2f} / {cp[0.99]:.2f} / {cp[0.999]:.2f}"
                   + ("  (rules)" if any(sim[str(q)] > SIM_RULE * cp[q] for q in cp) else ""))
    for f in chip.forest:
        z = np.array([r[f"VS{f}"] for r in res])
        ok = np.isfinite(z)
        if ok.any():
            say(lines, f"  VS q{f}: z defined in {ok.mean():.4f} of the flights; mean {np.mean(z[ok]):+.3f}, sd {np.std(z[ok]):.3f}; "
                       f"HEARD in {np.mean(z[ok] >= 5):.4f}, |z| >= 3 in {np.mean(np.abs(z[ok]) >= 3):.4f}")
        else:
            say(lines, f"  VS q{f}: no flip on either side in any flight (P(1|0) 0 in the snapshot)")
    scale = {}
    if chip.rider:
        for l in range(3):
            d = np.array([r[f"delta{l}"] for r in res])
            sd = np.array([r[f"sdelta{l}"] for r in res])
            ratio = float(np.nanstd(d) / np.nanmedian(sd))
            scale[str(l)] = ratio if ratio > SIM_RULE else 1.0
            say(lines, f"  delta q{chip.path[l]}: mean {np.nanmean(d):+.5f}, spread {np.nanstd(d):.5f} rad/us, median sigma_delta "
                       f"{np.nanmedian(sd):.5f}; ratio {ratio:.3f}" + ("  (sigma_delta scaled)" if ratio > SIM_RULE else ""))
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / f"drive_frame_{chip_name}_nulls_raw_{n}.json").write_text(json.dumps(res), encoding="utf-8")
    if n != N_NULLS:
        print(f"(a trial of {n} flights: the ruling file is written only from {N_NULLS})")
        return
    nulls_path(chip_name).write_text(json.dumps({"chip": chip_name, "qpy_sha256": chip.freeze["sha256"],
                                                 "freeze_sha256": chip.freeze_sha, "n": len(res),
                                                 "points": points, "delta_scale": scale}, indent=1), encoding="utf-8")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / f"{chip_name}_nulls.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


# ── reach: the letter size at which a statistic passes its point in half of 2000 simulated flights ──

def _draw_par(chip, need, L, A_true, rng):
    """{name: (ndraw, 2, 7)} half parities, drawn under the letters L."""
    out = {}
    for v in need:
        for b in BASES:
            nm = f"{v}_{b[0]}{b[1]}"
            p = A_true @ chip.p_true(nm, L)
            c = rng.multinomial(SHOTS_HALF, p / p.sum(), size=(REACH_DRAWS, 2)).astype(float)
            out[nm] = c @ SG.T / SHOTS_HALF
    return out


def _vec(hp, v, sel):
    """(ndraw, R) readings of variant v, both halves together."""
    return np.stack([hp[f"{v}_{b[0]}{b[1]}"][..., MI[m]].mean(1) for b, m in sel], axis=1)


def _drift_per_draw(hp, variants, sel, d):
    par = {f"{v}_{b[0]}{b[1]}": hp[f"{v}_{b[0]}{b[1]}"][d].mean(0) for v in variants for b in BASES}
    par_h = {k: [hp[k][d][0], hp[k][d][1]] for k in par}
    return arm_cov(par, par_h, variants, sel, 2 * SHOTS_HALF)


def reach_frame(chip, arm, sel, A_true, dc_day, pt, letter, seed, lo=0.0, hi=0.02):
    """eps at which T5 of arm (A or B), with the day's kappa Sigma (dc_day) or, if None, each simulated flight's own
    Sigma and drift matrix, passes the point pt in half of REACH_DRAWS flights drawn under letter(eps)."""
    variants = chip.variants(arm)
    five = [f"{arm}5{k}" for k in range(4)]

    def frac(eps):
        rng = np.random.default_rng(seed)
        need = variants if dc_day is None else five
        hp = _draw_par(chip, need, letter(eps), A_true, rng)
        Y = {v: _vec(hp, v, sel) for v in five}
        hits = 0
        W = dc_day.winv() if dc_day is not None else None
        for d in range(REACH_DRAWS):
            if dc_day is None:
                W = _drift_per_draw(hp, variants, sel, d)[1].winv()
            hits += frame_stat([Y[v][d] for v in five], W) > pt
        return hits / REACH_DRAWS
    return bisect(frac, lo, hi)


def reach_window(chip, A_true, nro, pt_base, dc_day, ptv, seed, lo=0.0, hi=0.5, k_ro=1.0):
    """The detuning on c (rad/us) at which W of B5(0) passes the point ptv in half of REACH_DRAWS flights drawn under
    it, the RO matrix redrawn in each flight, Sigma_RO from the true matrix; kappa Sigma the day's (dc_day) or each
    flight's own."""
    Ro = ro_cov(pt_base, B_IDX, A_true, nro, k_ro)
    variants = chip.variants("B")
    P0 = np.array([pt_base[b] for b, _ in B_IDX])
    Sm = np.array([SG[MI[m]] for _, m in B_IDX])

    def frac(det):
        rng = np.random.default_rng(seed)
        L = {("det", 6): det} if det else None
        need = variants if dc_day is None else ["B50"]
        hp = _draw_par(chip, need, L, A_true, rng)
        Ae = np.stack([rng.multinomial(int(nro[j]), A_true[:, j], size=REACH_DRAWS) / nro[j] for j in range(8)], axis=2)
        Y = _vec(hp, "B50", B_IDX)
        hits = 0
        Wfix = pinv_rank(dc_day.cov() + Ro)[0] if dc_day is not None else None
        for d in range(REACH_DRAWS):
            y0 = np.einsum("rk,kj,rj->r", Sm, Ae[d], P0)
            r = Y[d] - y0
            W = Wfix if Wfix is not None else pinv_rank(_drift_per_draw(hp, variants, B_IDX, d)[1].cov() + Ro)[0]
            hits += float(r @ W @ r) > ptv
        return hits / REACH_DRAWS
    return bisect(frac, lo, hi)


def bisect(frac, lo, hi, steps=12, hi_max=None):
    """The size at which frac crosses one half: hi doubled until frac(hi) >= 1/2 (inf past hi_max), then halved."""
    hi_max = hi_max if hi_max is not None else 16 * hi
    while frac(hi) < 0.5:
        lo, hi = hi, 2 * hi
        if hi > hi_max:
            return float("inf")
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if frac(mid) >= 0.5:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def reach_day(chip, comp, out, lines):
    """The reaches the verdict rules report, with the day's kappa Sigma and RO matrix, on the standing variants."""
    A, nro = readout_matrix(comp)
    par, par_h = {}, {}
    for nm, (pc, fo, n) in comp.items():
        if not nm.startswith("RO_"):
            par[nm] = SG @ (pc.sum(0) / n.sum())
            par_h[nm] = [SG @ (pc[h] / n[h]) for h in range(2)]
    r = {}
    n = comp["A50_XX"][2].sum()
    if out["VF"]["T5"]["class"] == "FRAME-BLIND" and len(out["VF"]["T5"]["variants"]) == 4:
        standing = [v for v in chip.variants("A") if v not in out["A_void"]]
        dc = arm_cov(par, par_h, standing, A_IDX, n)[1]
        r["VF_reach_eps"] = reach_frame(chip, "A", A_IDX, A, dc, out["VF"]["T5"]["points"][P_TURNS], xt_b, seed=31)
        say(lines, f"   VF FRAME-BLIND reach (50 % point): eps {100 * r['VF_reach_eps']:.3f} % on (b -> a) and (b -> c)")
    if chip.rider and "T" in out["VFB"]["T5"] and out["VFB"]["T5"]["class"] == "FRAME-BLIND" and len(out["VFB"]["T5"]["variants"]) == 4:
        standing = [v for v in chip.variants("B") if v not in out["B_void"]]
        dc = arm_cov(par, par_h, standing, B_IDX, n)[1]
        r["VFB_reach_eps"] = reach_frame(chip, "B", B_IDX, A, dc, out["VFB"]["T5"]["points"][P_TURNS], xt_b, seed=32)
        say(lines, f"   VFB FRAME-BLIND reach on b: eps {100 * r['VFB_reach_eps']:.3f} %")
    if chip.rider and out.get("WB", {}).get("verdict") == "SHUT":
        standing = [v for v in chip.variants("B") if v not in out["B_void"]]
        dc = arm_cov(par, par_h, standing, B_IDX, n)[1]
        pt = {b: chip.p_true(f"B50_{b[0]}{b[1]}") for b in BASES}
        r["WB_reach_det"] = reach_window(chip, A, nro, pt, dc, out["WB"]["points"][P_SHUT], seed=33, k_ro=out["kappa_RO"])
        say(lines, f"   WB SHUT reach: detuning on c {r['WB_reach_det']:.4f} rad/us ({r['WB_reach_det'] / (2 * math.pi) * 1e3:.1f} kHz)")
    return r


# ── the freeze numbers of the X^3 window page, section 4 ──────────────────

def sx_tilted(ops, qs, t):
    out = []
    for o in ops:
        if o[0] == "sx" and o[1][0] in qs:
            out += [("rz", o[1], [-t], 0.0), o, ("rz", o[1], [t], 0.0)]
        else:
            out.append(o)
    return out


def freeze_numbers_a(chip, nulls_):
    """The drive-frame page's per-chip numbers on a chip without the X^3 arms: A5(0)'s baseline, its leaning on the
    snapshot's T1 and T2, and the FRAME-BLIND reach of T5."""
    lines = []
    s = lambda x: say(lines, x)
    A, n, nro = chip.A_snap, 2.0 * SHOTS_HALF, np.full(8, 2.0 * SHOTS_HALF)
    pa = {b: chip.p_true(f"A50_{b[0]}{b[1]}") for b in BASES}
    ya = {b: SG @ (A @ pa[b]) for b in BASES}
    s(f"{chip.name}: path {chip.path}, forest {chip.forest}, circuits {chip.freeze['sha256'][:16]}; freeze snapshot "
      f"{chip.freeze['snapshot']['calibration']}: T1 {[round(x) for x in chip.T1]} us, T2 {[round(x) for x in chip.T2]} us, "
      f"(P(1|0), P(0|1)) on (a, c, b) {[(round(a, 4), round(b, 4)) for a, b in chip.ro_snap]}; oldest dated property "
      f"{chip.freeze['oldest_property_utc']}")
    s(f"baseline A5(0): parse gate <X0Z2> {ya[('X', 'Z')][MI[(1, 1, 0)]]:+.3f}  <Z0X2> {ya[('Z', 'X')][MI[(1, 1, 0)]]:+.3f}; "
      f"S = <Y0Z2> {ya[('Y', 'Z')][MI[(1, 1, 0)]]:+.4f}")
    sS = math.sqrt(shot_cov(ya, A_IDX, n)[S_POS, S_POS] + ro_cov(pa, A_IDX, A, nro)[S_POS, S_POS])
    T1, T2 = list(chip.T1), list(chip.T2)
    chip.T1, chip.T2, chip._ptrue = [t / 2 for t in T1], [t / 2 for t in T2], {}
    S_h = SG[MI[(1, 1, 0)]] @ (A @ chip.p_true("A50_YZ"))
    chip.T1, chip.T2, chip._ptrue = T1, T2, {}
    s(f"the baseline with the snapshot's T1 and T2 halved: S moves by {(S_h - ya[('Y', 'Z')][MI[(1, 1, 0)]]) / sS:+.3f} sigma_S")
    ra = pinv_rank(shot_cov(ya, A_IDX, n))[1]
    eps_a = reach_frame(chip, "A", A_IDX, A, None, ruling("T5", 3 * ra, nulls_, True)[0][P_TURNS], xt_b, seed=42)
    s(f"FRAME-BLIND reach of T5 on a drift-free chip (each simulated flight with its own Sigma and drift matrix): eps {100 * eps_a:.3f} %")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / f"{chip.name}_freeze_numbers.txt").write_text(chr(10).join(lines) + chr(10), encoding="utf-8", newline=chr(10))


def freeze_numbers(chip_name):
    chip = Chip(chip_name)
    nulls_ = load_nulls(chip_name)
    if not chip.rider:
        return freeze_numbers_a(chip, nulls_)
    lines = []
    s = lambda x: say(lines, x)
    A = chip.A_snap
    nro = np.full(8, 2.0 * SHOTS_HALF)
    n = 2.0 * SHOTS_HALF
    pt = {b: chip.p_true(f"B50_{b[0]}{b[1]}") for b in BASES}
    y0 = {b: SG @ (A @ pt[b]) for b in BASES}
    xx = y0[("X", "X")]
    s(f"{chip_name}: path {chip.path}, forest {chip.forest}, circuits {chip.freeze['sha256'][:16]}; freeze snapshot "
      f"{chip.freeze['snapshot']['calibration']}: T1 {[round(x) for x in chip.T1]} us, T2 {[round(x) for x in chip.T2]} us, "
      f"(P(1|0), P(0|1)) on (a, c, b) {[(round(a, 4), round(b, 4)) for a, b in chip.ro_snap]}")
    s(f"baseline B5(0), snapshot RO matrix: XX basis <X0> {xx[MI[(1, 0, 0)]]:+.3f}  <X2> {xx[MI[(0, 1, 0)]]:+.3f}  "
      f"<X1> {xx[MI[(0, 0, 1)]]:+.3f}  <X0X1X2> {xx[MI[(1, 1, 1)]]:+.3f}")
    s("  pooled single-site readings: " + "  ".join(
        f"<{k}> {np.mean([y0[B_IDX[i][0]][MI[B_IDX[i][1]]] for i in ii]):+.4f}" for k, ii in POOLED.items()))
    pa = {b: chip.p_true(f"A50_{b[0]}{b[1]}") for b in BASES}
    ya = {b: SG @ (A @ pa[b]) for b in BASES}
    s(f"baseline A5(0): parse gate <X0Z2> {ya[('X', 'Z')][MI[(1, 1, 0)]]:+.3f}  <Z0X2> {ya[('Z', 'X')][MI[(1, 1, 0)]]:+.3f}; "
      f"S = <Y0Z2> {ya[('Y', 'Z')][MI[(1, 1, 0)]]:+.4f}")
    Sig = shot_cov(y0, B_IDX, n)
    Ro = ro_cov(pt, B_IDX, A, nro)
    Wt, rk = pinv_rank(Sig + Ro)
    s(f"W on {rk} readings: 95 % point {CHI2.ppf(0.95, rk):.1f}, 99.9 % point {CHI2.ppf(0.999, rk):.1f}; frame-blind yardsticks, "
      f"noncentrality of W (shot + RO covariance, no drift):")
    yb = np.array([y0[b][MI[m]] for b, m in B_IDX])
    BLIND = (("five Z letters (October fit)", FIVE, None),
             ("detuning on a +0.135 rad/us (21.5 kHz)", {("det", 4): 0.135}, None),
             ("detuning on c +0.135 rad/us", {("det", 6): 0.135}, None),
             ("detuning on b +0.1 rad/us", {("det", 5): 0.1}, None),
             ("sx and x amplitude on a +1 %", {("sxa", 4): 0.01}, None),
             ("sx and x amplitude on c +1 %", {("sxa", 6): 0.01}, None),
             ("cz local phase 0.03, all", {("czp", p, q): 0.03 for p in rl.PAIRS for q in p}, None),
             ("sx tilted 0.06 rad against x, on a and c", None, ([0, 2], 0.06)))
    for label, L, tilt in BLIND:
        y = []
        for b, m in B_IDX:
            nm = f"B50_{b[0]}{b[1]}"
            ops = sx_tilted(chip.ops[nm], *tilt) if tilt else None
            p = chip.p_true(nm, L, ops=ops)
            y.append(SG[MI[m]] @ (A @ p))
        d = np.array(y) - yb
        s(f"  {label:42s} {float(d @ Wt @ d):8.1f}   largest shift {np.abs(d).max():.4f}")
    # how much the baseline leans on the snapshot's T1 and T2: S of A5(0) and W of B5(0) with both halved
    sS = math.sqrt(shot_cov(ya, A_IDX, n)[S_POS, S_POS] + ro_cov(pa, A_IDX, A, nro)[S_POS, S_POS])
    T1, T2 = list(chip.T1), list(chip.T2)
    chip.T1, chip.T2, chip._ptrue = [t / 2 for t in T1], [t / 2 for t in T2], {}
    S_h = SG[MI[(1, 1, 0)]] @ (A @ chip.p_true("A50_YZ"))
    yh = np.array([SG[MI[m]] @ (A @ chip.p_true(f"B50_{b[0]}{b[1]}")) for b, m in B_IDX])
    chip.T1, chip.T2, chip._ptrue = T1, T2, {}
    s(f"the baseline with the snapshot's T1 and T2 halved: S moves by {(S_h - ya[('Y', 'Z')][MI[(1, 1, 0)]]) / sS:+.3f} sigma_S, "
      f"W's noncentrality {float((yh - yb) @ Wt @ (yh - yb)):.2f}")
    t0 = time.time()
    eps = reach_frame(chip, "B", B_IDX, A, None, ruling("T5B", 3 * rk, nulls_, True)[0][P_TURNS], xt_b, seed=41)
    s(f"reach of T5_B at its 99.9 % point (each simulated flight with its own Sigma and drift matrix, October phases "
      f"carried by role): eps {100 * eps:.3f} %  ({time.time() - t0:.0f} s)")
    ra = pinv_rank(shot_cov(ya, A_IDX, n))[1]
    eps_a = reach_frame(chip, "A", A_IDX, A, None, ruling("T5", 3 * ra, nulls_, True)[0][P_TURNS], xt_b, seed=42)
    s(f"the same for T5 of the drive-frame page (its FRAME-BLIND reach on a drift-free chip): eps {100 * eps_a:.3f} %")
    det = reach_window(chip, A, nro, pt, None, ruling("W", rk, nulls_, True)[0][P_SHUT], seed=43)
    s(f"reach of a SHUT (W past its 95 % point in half the flights), detuning on c: {det:.4f} rad/us ({det / (2 * math.pi) * 1e3:.1f} kHz)")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / f"{chip_name}_freeze_numbers.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


# ── the gate (offline) ────────────────────────────────────────────────────

def logical_probs(ops):
    """Path probabilities (k = 4 q0 + 2 q2 + q1) of the gate sequence alone, delays dropped, ideal and instantaneous."""
    psi = np.zeros(8, dtype=complex)
    psi[0] = 1
    for name, qs, prm, _ in ops:
        if name == "delay":
            continue
        if name == "rz":
            U = rl.on(qs[0], rl.rz(prm[0]))
        elif name == "sx":
            U = rl.on(qs[0], rl.rx(math.pi / 2))
        elif name == "x":
            U = rl.on(qs[0], rl.rx(math.pi))
        elif name == "cz":
            U = rl.czm(*sorted(qs))
        psi = U @ psi
    return rl.probs_rho(np.outer(psi, psi.conj()))


def gate(chip_name):
    chip = Chip(chip_name)
    lines = []
    s = lambda x: say(lines, x)
    s(f"{chip_name}: frozen file sha256 {chip.freeze['sha256']}, {len(chip.circs)} circuits, path {chip.path}, "
      f"forest {chip.forest}")
    assert sorted(chip.circs) == FZ.expected_names(chip_name)
    from qiskit import qpy
    with open(chip.qpy, "rb") as f:
        FZ.assert_structure(qpy.load(f), chip.freeze)
    s("(1) structure: turns as first and last operations on the turned qubit, otherwise the unturned circuit; Q's idles; "
      "the forest read and not gated, every measurement at one time; RO by x alone: pass")
    with open(FZ.APRIL_QPY, "rb") as f:
        april = {c.name: c for c in qpy.load(f)}
    pmap = dict(zip(FZ.APRIL_PATH, chip.path))
    same = 0
    for arm, cls, tail in (("A", "pi2_odd_pure", chip.freeze["ytail"]), ("B", "pi2_even_nontruly", chip.freeze["xtail"])):
        if arm == "B" and not chip.rider:
            continue
        for b0, b2 in BASES:
            m = FZ.per_qubit(chip.circs[f"{arm}50_{b0}{b2}"])
            t = FZ.per_qubit(april[f"{cls}_{b0}{b2}"])
            strip = lambda ops, mp=None: [(o[0], tuple((mp or {}).get(q, q) for q in o[1]), o[2]) for o in ops
                                          if o[0] not in ("delay", "measure")]
            ok = True
            for lq, ap in enumerate(FZ.APRIL_PATH):
                got = strip(m.get(chip.path[lq], []))
                want = strip(t.get(ap, []), pmap)
                if lq == 1:
                    want = want + [(n_, (chip.path[1],), tuple(p_)) for n_, p_ in tail]
                ok &= got == want
            same += ok
    s(f"(2) A5(0) and B5(0) = April's flown circuits per qubit, carried to the path, + the tail on b: {same}/{18 if chip.rider else 9}")
    assert same == (18 if chip.rider else 9)
    worst = 0.0
    for nm in chip.circs:
        var = nm.split("_")[0]
        if var[0] in "AB" and var[2] != "0":
            ref = f"{var[0]}50_{nm.split('_')[1]}"
            worst = max(worst, float(np.abs(logical_probs(chip.ops[nm]) - logical_probs(chip.ops[ref])).max()))
    s(f"(3) every turned circuit against its unturned one, ideal instantaneous gates: max |dp| {worst:.1e}")
    assert worst < 1e-12
    inf = [1e12] * 3
    za = max(abs(float(SG[MI[m]] @ rl.probs_rho(rl.replay(chip.ops[f"A50_{b[0]}{b[1]}"], ROLE, inf, inf)))) for b, m in A_IDX)
    s(f"(4) the zero: replay without decay, ideal readout, largest odd-Y reading of A5(0) {za:.1e}")
    assert za < 1e-9
    if chip.rider:
        zb = max(abs(float(SG[MI[m]] @ rl.probs_rho(rl.replay(chip.ops[f"B50_{b[0]}{b[1]}"], ROLE, inf, inf))))
                 for b in BASES for m in rl.MASKS if forbidden_b(b, m))
        s(f"    largest forbidden reading of B5(0) (all 36) {zb:.1e}")
        assert zb < 1e-9
    for arm in ("A", "B") if chip.rider else ("A",):
        sel = A_IDX if arm == "A" else B_IDX
        for label, L in (("five Z letters", FIVE), ("crosstalk alone (October fit, by role)", XT_FIT)):
            base = np.array([SG[MI[m]] @ chip.p_true(f"{arm}50_{b[0]}{b[1]}", L) for b, m in sel])
            dmax = 0.0
            for v in chip.variants(arm)[1:]:
                y = np.array([SG[MI[m]] @ chip.p_true(f"{v}_{b[0]}{b[1]}", L) for b, m in sel])
                dmax = max(dmax, float(np.abs(y - base).max()))
            s(f"(5) {arm}: {label}: largest change of a reading under the turns {dmax:.1e}")
            if L is FIVE:
                # the replay's step error where a float sum of durations moves a time boundary by one rounding: 1e-5
                assert dmax < 1e-4
            else:
                assert dmax > 1e-2
    rng = np.random.default_rng(20261113)
    comp = draw(chip, true_probs(chip), rng)
    RESULTS.mkdir(exist_ok=True)
    fake = RESULTS / f"drive_frame_{chip_name}_simulated_flight.json"
    fake.write_text(json.dumps({"job_id": "simulated", "backend": chip_name, "path": chip.path, "forest": chip.forest,
                                "qpy_sha256": chip.freeze["sha256"], "snapshot": chip.freeze["snapshot"],
                                "counts": to_counts(chip, comp, rng),
                                "turn_ok": {k: True for k in chip.variants("A") + (chip.variants("B") + chip.names("Q_") if chip.rider else [])}},
                               default=str), encoding="utf-8")
    t0 = time.time()
    out = analyze(fake)
    s(f"(6) a simulated flight under the baseline through --analyze ({time.time() - t0:.1f} s): P {out['P']['pass']}, "
      f"V1 {out['V1']['verdict']}, VF {out['VF']['T5']['class']}/{out['VF']['T4']['class']}/{out['VF']['T6']['class']}"
      + (f"; PB {out['PB']['pass']}, WB {out['WB']['verdict']}, VFB {out['VFB']['T5']['class']}, "
         f"VQ {[out['VQ'][str(l)]['verdict'] for l in range(3)]}" if chip.rider else ""))
    assert out["P"]["pass"] and (not chip.rider or out["PB"]["pass"])
    comp = draw(chip, true_probs(chip, XT_FIT), np.random.default_rng(7))
    out = analyze_compact(chip, comp)
    s(f"(7) one draw under the crosstalk-alone fit (by role): V1 {out['V1']['verdict']} (z {out['V1']['z']:+.2f}), "
      f"T5 {out['VF']['T5']['T']:.0f} {out['VF']['T5']['class']}, T4 {out['VF']['T4']['class']}, T6 {out['VF']['T6']['class']}"
      + (f"; WB {out['WB']['verdict']}, T5_B {out['VFB']['T5']['T']:.0f} {out['VFB']['T5']['class']}" if chip.rider else ""))
    if chip.rider:
        # (8) the colouring's scale on a chip whose truth is the snapshot, and on one whose T1 and T2 are half of it
        T1, T2 = list(chip.T1), list(chip.T2)
        for label, f in (("the snapshot", 1.0), ("T1 and T2 at half the snapshot's", 0.5)):
            chip.T1, chip.T2, chip._ptrue = [t * f for t in T1], [t * f for t in T2], {}
            P = true_probs(chip)
            chip.T1, chip.T2, chip._ptrue = T1, T2, {}
            for sd in (11, 12):
                c = analyze_compact(chip, draw(chip, P, np.random.default_rng(sd)), readings=True)["colouring"]
                s(f"(8) truth {label}, draw {sd}: colouring scale {c['scale_L3']:.2f} +- {c['scale_L3_sigma']:.2f} "
                  f"(singles {c['scale_L1']:.2f} +- {c['scale_L1_sigma']:.2f}); replay with arm Q's T2: -ln|XXX| "
                  f"{c.get('replay_with_Q', {}).get('L3', float('nan')):.4f} against the snapshot replay's "
                  f"{c['replay_snapshot']['L3']:.4f}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / f"{chip_name}_gate.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


# ── the flight ────────────────────────────────────────────────────────────

def committed_row(chip_name):
    """The row written for this chip in section 6 of the drive-frame page: | chip | path | `circuit sha` | `record sha` |."""
    m = re.findall(r"^\| " + re.escape(chip_name) + r" \| (\[[0-9, ]+\]) \| `([0-9a-f]{64})` \| `([0-9a-f]{64})` \|",
                   PAGE.read_text(encoding="utf-8"), flags=re.M)
    return m[0] if len(m) == 1 else None


def committed_clean(paths):
    """Every file tracked and unmodified against HEAD."""
    bad = []
    for p in paths:
        rel = str(Path(p).resolve().relative_to(REPO)).replace("\\", "/")
        if subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=REPO, capture_output=True).returncode:
            bad.append(f"{rel} is not tracked")
        elif subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=REPO).returncode:
            bad.append(f"{rel} differs from HEAD")
    return bad


def day_of(svc, be, chip):
    """Section 6 of the drive-frame page. Any failure aborts; there is no override."""
    fails = []
    have = sha256(chip.qpy)
    row = committed_row(chip.name)
    if row is None:
        fails.append("no row written for this chip in section 6 of the page")
    else:
        if row[0] != str(chip.path):
            fails.append(f"the page's path {row[0]} is not the freeze's {chip.path}")
        if have != row[1]:
            fails.append(f"frozen file hash {have} != the page's {row[1]}")
        if chip.freeze_sha != row[2]:
            fails.append(f"freeze record {chip.freeze_sha} != the page's {row[2]}")
    for d in (RESULTS, REPO / chip.freeze["record_dir"]):
        for st in d.glob(f"drive_frame_{chip.name}_*.submitted.json") if d.exists() else []:
            try:
                prev = json.loads(st.read_text(encoding="utf-8"))
            except Exception:
                continue
            if prev.get("qpy_sha256") == have:
                fails.append(f"already submitted from this frozen file ({st.name}, job {prev.get('job_id')}): use --recover")
    if chip.freeze["sha256"] != have:
        fails.append("the freeze record does not match the frozen file")
    try:
        load_nulls(chip.name)
    except SystemExit as exc:
        fails.append(f"nulls: {exc}")
    dc_out = REPO / "simulations" / "results" / "drive_frame" / f"{chip.name}_drift_check.txt"
    if not dc_out.exists() or chip.freeze_sha not in dc_out.read_text(encoding="utf-8").splitlines()[0]:
        fails.append("the drift check was not run on this freeze record")
    files = [PAGE, REPO / "experiments" / "X3_WINDOW_HARDWARE_PREDICTION.md", Path(FZ.__file__), Path(__file__),
             chip.qpy, chip.meta, nulls_path(chip.name), Path(rl.__file__),
             REPO / "simulations" / "results" / "drive_frame" / f"{chip.name}_drift_check.txt",
             REPO / "simulations" / "drive_frame_drift_check.py",
             REPO / "simulations" / "reality_window_crosstalk_fit.py"]
    fails += committed_clean(files)
    from qiskit import qpy
    with open(chip.qpy, "rb") as f:
        try:
            FZ.assert_structure(qpy.load(f), chip.freeze)
        except AssertionError as exc:
            fails.append(f"structure: {exc}")
    now = datetime.now(timezone.utc)
    age = (now - datetime.fromisoformat(chip.freeze["oldest_property_utc"])).total_seconds() / 3600
    # the property ages are recorded, not gated (section 6)
    try:
        live = FZ.durations_of(be, chip.path, chip.forest)
        if live != chip.freeze["durations_dt"]:
            fails.append(f"the target's durations {live} differ from the freeze's {chip.freeze['durations_dt']}")
    except Exception as exc:
        live = repr(exc)
        fails.append(f"the target lacks an operation of the frozen circuits (a faulty qubit?): {exc!r}")
    if be.dt != chip.freeze["dt_s"]:
        fails.append(f"dt {be.dt} differs from the freeze's {chip.freeze['dt_s']}")
    p = be.properties()
    recs, bonds, f = FZ.path_check(p, be.target, chip.path, now)
    fails += f
    rem = svc.usage().get("usage_remaining_seconds")
    if rem is None or rem < MIN_REMAINING_S:
        fails.append(f"allocation remaining {rem}")
    st = be.status()
    if not st.operational:
        fails.append("backend not operational")
    snap = {"calibration": str(p.last_update_date), "taken_utc": now.isoformat(), "path": {str(q): r for q, r in recs.items()},
            "bonds": bonds, "remaining_s": rem, "pending_jobs": st.pending_jobs, "operational": st.operational,
            "freeze_age_h": age, "durations_dt": live}
    return fails, snap


def check_retrieved(chip, job, order):
    """The circuits retrieved from the job (the client's submitted circuits, as IBM stores them) against the frozen
    ones, per variant and Q circuit: True/False, or None when they cannot be read."""
    out = {k: None for k in chip.variants("A") + (chip.variants("B") + chip.names("Q_") if chip.rider else [])}
    try:
        pubs = job.inputs["pubs"]
        circs = [pb[0] if isinstance(pb, (tuple, list)) else pb.circuit for pb in pubs]
        if len(circs) != len(order):
            print(f"retrieved {len(circs)} circuits for {len(order)} PUBs: not checked")
            return out
        got = {}
        for (nm, h), c in zip(order, circs):
            got.setdefault(nm, []).append(c)
        for key in out:
            nms = [n for n in chip.circs if n.split("_")[0] == key] if not key.startswith("Q_") else [key]
            if any(len(got.get(nm, [])) != 2 for nm in nms):
                continue
            keep_delay = key.startswith("Q_")
            ok = True
            for nm in nms:
                ref = FZ.per_qubit(chip.circs[nm])
                for c in got[nm]:
                    ret = FZ.per_qubit(c)
                    for q in chip.path:
                        a = [(o[0], o[2]) for o in ret.get(q, []) if keep_delay or o[0] != "delay"]
                        b = [(o[0], o[2]) for o in ref.get(q, []) if keep_delay or o[0] != "delay"]
                        if len(a) != len(b) or any(x[0] != y[0] or not np.allclose(x[1], y[1], atol=1e-9) for x, y in zip(a, b)):
                            ok = False
            out[key] = ok
    except Exception as exc:
        print(f"retrieved circuits not checked: {exc!r}")
    return out


def write_counts(chip, job, order, snap, stamp, qpy_sha):
    """The counts and the order go to disk first; the checks, metrics and usage are added after."""
    res = job.result()
    counts = {}
    for pos, (nm, h) in enumerate(order):
        counts.setdefault(nm, [None, None])[h] = res[pos].data.meas.get_counts()
    rec = REPO / chip.freeze["record_dir"]
    rec.mkdir(parents=True, exist_ok=True)
    out = rec / f"drive_frame_{chip.name}_{stamp}.json"
    record = {"job_id": job.job_id(), "backend": chip.name, "path": chip.path, "forest": chip.forest,
              "shots_per_half": SHOTS_HALF, "qpy_sha256": qpy_sha, "order": order,
              "freeze_snapshot": chip.freeze["snapshot"], "freeze_sha256": chip.freeze_sha, "snapshot": snap,
              "counts": counts}
    out.write_text(json.dumps(record, default=str, indent=1), encoding="utf-8")
    print("counts persisted ->", out)
    for key, fn in (("result_metadata", lambda: res.metadata), ("job_metrics", job.metrics), ("billed_usage", job.usage),
                    ("turn_ok", lambda: check_retrieved(chip, job, order))):
        try:
            record[key] = fn()
        except Exception as exc:
            record[key] = repr(exc)
    out.write_text(json.dumps(record, default=str, indent=1), encoding="utf-8")
    print("billed usage", record["billed_usage"])
    return out


def hardware(chip_name):
    from qiskit_ibm_runtime import SamplerV2
    chip = Chip(chip_name)
    svc = FZ.service()
    be = svc.backend(chip_name)
    fails, snap = day_of(svc, be, chip)
    RESULTS.mkdir(exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    print("day-of snapshot:", json.dumps(snap, default=str))
    if fails:
        (RESULTS / f"drive_frame_{chip_name}_{stamp}.aborted.json").write_text(
            json.dumps({"fails": fails, "snapshot": snap}, default=str, indent=1), encoding="utf-8")
        raise SystemExit("ABORT: " + "; ".join(fails))
    order = [(nm, h) for nm in sorted(chip.circs) for h in (0, 1)]
    random.Random(chip.freeze["shuffle_seed"]).shuffle(order)
    assert len(order) == chip.freeze["n_pubs"]
    stub = RESULTS / f"drive_frame_{chip_name}_{stamp}.submitted.json"
    rec_dir = REPO / chip.freeze["record_dir"]
    rec_dir.mkdir(parents=True, exist_ok=True)
    stubs = (stub, rec_dir / stub.name)
    qpy_sha = sha256(chip.qpy)
    record = {"job_id": None, "backend": chip_name, "stamp": stamp, "snapshot": snap, "order": order, "qpy_sha256": qpy_sha,
              "freeze_sha256": chip.freeze_sha}

    def save():
        for st in stubs:
            st.write_text(json.dumps(record, default=str), encoding="utf-8")
    save()
    sampler = SamplerV2(mode=be)
    sampler.options.default_shots = SHOTS_HALF
    sampler.options.max_execution_time = MAX_EXEC_S
    sampler.options.dynamical_decoupling.enable = False
    sampler.options.twirling.enable_gates = False
    sampler.options.twirling.enable_measure = False
    sampler.options.execution.init_qubits = True
    sampler.options.execution.meas_type = "classified"
    sampler.options.environment.job_tags = ["drive_frame", chip_name, qpy_sha[:12]]
    try:
        job = sampler.run([chip.circs[nm] for nm, _ in order])
    except Exception as exc:
        try:
            recent = [j.job_id() for j in svc.jobs(limit=5, backend_name=chip_name)]
        except Exception as exc2:
            recent = [f"job listing failed: {exc2!r}"]
        record["submit_error"] = repr(exc)
        record["recent_jobs"] = recent
        save()
        raise SystemExit(f"submission raised {exc!r}; the five newest jobs on {chip_name} are {recent}, recorded in "
                         f"{stub}: if one of them is this flight, put its id in the stub and --recover")
    record["job_id"] = job.job_id()
    save()
    print("submitted", record["job_id"], "->", stub)
    analyze(write_counts(chip, job, order, snap, stamp, qpy_sha), reach=True)


def recover(stub_path):
    stub = json.loads(Path(stub_path).read_text(encoding="utf-8"))
    if not stub.get("job_id"):
        raise SystemExit("the stub has no job id: the submission did not return one; look the job up by its time")
    chip = Chip(stub["backend"])
    if stub["qpy_sha256"] != sha256(chip.qpy) or stub.get("freeze_sha256", chip.freeze_sha) != chip.freeze_sha:
        raise SystemExit("the stub was submitted under another frozen file or freeze record than the ones on disk; check "
                         "out the commit it was submitted from")
    job = FZ.service().job(stub["job_id"])
    analyze(write_counts(chip, job, [tuple(x) for x in stub["order"]], stub["snapshot"], stub["stamp"], stub["qpy_sha256"]),
            reach=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        raise SystemExit(__doc__)
    if a[0] == "--gate":
        gate(a[1])
    elif a[0] == "--nulls":
        nulls(a[1], int(a[2]) if len(a) > 2 else N_NULLS)
    elif a[0] == "--freeze-numbers":
        freeze_numbers(a[1])
    elif a[0] == "--hardware":
        hardware(a[1])
    elif a[0] == "--recover":
        recover(a[1])
    elif a[0] == "--analyze":
        analyze(a[1], reach="--reach" in a)
    else:
        raise SystemExit(__doc__)
