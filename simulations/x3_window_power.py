"""The X^3 window on YZ+ZY under a turned drive frame: predictions from below for arms B and Q.

experiments/X3_WINDOW_HARDWARE_PREDICTION.md section 4. Arm B is the 2026-04-30 flown YZ+ZY circuit (pi2_even_nontruly) on
Marrakesh [4, 5, 6], gate for gate, with q1 read in X (the October freezer's from_april with an X tail in place of the Y
tail), scheduled as late as possible. Its forbidden readings are every parity of (q0, q2, q1) whose Pauli string carries
an odd number of Y and Z letters on q0 and q2 (anticommutes with X^3) or an odd number of Y letters (YZ+ZY is imaginary,
|+-+> real). The turn: rz(phi) as the first and rz(-phi) as the last operation on one qubit, a logical identity
(experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md section 1). Arm Q: one path qubit at a time, sx, an idle of tau, read in
X or Y, the other two in |0>.

  (0) April's flown YZ+ZY circuits under the October letter sets: pooled single-site readings and the sum of squared
      residuals over the ten readings the class forbids on (q0, q2), both schedules;
  (1) B's baseline (no letter, October snapshot, October RO matrix): the parse-gate readings and the single-site offsets;
  (2) the noise of the day's RO matrix carried into a baseline comparison: its covariance J C J' over B5(0)'s 14
      readings beside the shot covariance (and, for comparison, the same rates on all 36 forbidden readings), and the false-positive rate of W at its 99 % point without it (with
      it the rate is 1 % by construction, the draws being Gaussian from the same covariance);
  (3) frame-blind letters (Z letters, own-frame letters): invariance under the turn, and the noncentrality they give W;
      and pure dephasing alone, passed through the compiled gates;
  (4) turning letters (October's crosstalk fits): the frame statistics T5, T4, T6 and their power; the turn-averaged
      mean's offset from the baseline (the zeroth harmonic, which a transverse letter carries at second order);
  (5) the bound: power of T5 against a crosstalk eps on b's two pairs at the October phases;
  (6) arm Q, one qubit at a time with the other two in |0>: the replay under a detuning, the four-parameter fit (delta,
      theta0, contrast, decay) recovering it, and the spread of the fitted delta over noisy redraws;
  (7) a drift factor estimated from the data under real drift: the false-positive rates at the 99 % and 99.9 % points.

    python simulations/x3_window_power.py   ->   simulations/results/x3_window/x3_window_power.txt
"""
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as CHI2, ncx2 as NCX2

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "simulations"))
sys.path.insert(0, str(REPO / "simulations" / "flight"))
import run_reality_window as R                       # noqa: E402
import reality_window_letters as rl                  # noqa: E402
import reality_window_crosstalk_fit as CF            # noqa: E402
import drive_frame_power as DP                       # noqa: E402

from qiskit import QuantumCircuit, transpile, qpy   # noqa: E402

SHOTS = 8192          # per circuit, both halves
RO_SHOTS = 8192       # per RO circuit, both halves
KAPPA_OCT = 1.27       # October's drift factor over all its fit arms
KAPPA_A = 1.96         # the same over arm A's 36 odd-Y readings alone
OUT = REPO / "simulations" / "results" / "x3_window" / "x3_window_power.txt"
APRIL_JSON = REPO / "data" / "ibm_f83_signature_april2026" / "f83_signature_ibm_marrakesh_20260430_190035.json"
LINES = []
PHIS = [0.0, math.pi / 2, math.pi, 3 * math.pi / 2]
TAUS = [0.0, 1.0, 2.0, 4.0]
SIGN = {m: np.array([(-1.0) ** (m[0] * (k >> 2 & 1) + m[1] * (k >> 1 & 1) + m[2] * (k & 1)) for k in range(8)]) for m in rl.MASKS}


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def tail(gates):
    q = QuantumCircuit(1)
    for g in gates:
        getattr(q, g)(0)
    return [(i.operation.name, [float(p) for p in i.operation.params])
            for i in transpile(q, basis_gates=["rz", "sx", "x"], optimization_level=1).data]


def dur(nm, qs):
    if nm == "rz":
        return 0.0
    if nm in ("sx", "x"):
        return R.DUR_US[nm][R.PATH[qs[0]]]
    return R.DUR_US["cz"][tuple(sorted(R.PATH[q] for q in qs))]


def alap(ops):
    """As late as possible against a common end, with explicit delays (agrees with qiskit's ALAPScheduleAnalysis +
    PadDelay on the October arm A, operation for operation)."""
    t = [0.0, 0.0, 0.0]
    ends = []
    for nm, qs, prm in reversed(ops):
        end = max(t[q] for q in qs)
        for q in qs:
            t[q] = end + dur(nm, qs)
        ends.append(end)
    T = max(t)
    sched = [(T - e - dur(nm, qs), nm, qs, prm) for (nm, qs, prm), e in zip(reversed(ops), ends)][::-1]
    clock = [0.0, 0.0, 0.0]
    res = []
    for st, nm, qs, prm in sched:
        for q in qs:
            if st > clock[q] + 1e-12:
                res.append(("delay", [q], [], st - clock[q]))
                clock[q] = st
        start = max(clock[q] for q in qs)
        res.append((nm, qs, prm, dur(nm, qs)))
        for q in qs:
            clock[q] = start + dur(nm, qs)
    return res


def b_ops():
    with open(R.APRIL_QPY, "rb") as f:
        april = {c.name: c for c in qpy.load(f)}
    xt = tail(["h"])
    m = {p: i for i, p in enumerate(R.PATH)}
    out = {}
    for b0 in "XYZ":
        for b2 in "XYZ":
            tc = april[f"pi2_even_nontruly_{b0}{b2}"]
            ops = []
            for inst in tc.data:
                qs = [tc.find_bit(q).index for q in inst.qubits]
                nm = inst.operation.name
                if nm in ("measure", "barrier", "delay") or not all(q in m for q in qs):
                    continue
                ops.append([nm, [m[q] for q in qs], [float(p) for p in inst.operation.params]])
            for nm, prm in xt:
                ops.append([nm, [1], prm])
            out[(b0, b2)] = alap(ops)
    return out


def q_ops():
    """Arm Q: one path qubit at a time, the other two left in |0>: sx, an idle of tau, then read in X (h) or in Y
    (sdg, h). 3 qubits x 2 bases x 4 idles = 24 circuits."""
    out = {}
    for q in range(3):
        for basis, gates in (("X", ["h"]), ("Y", ["sdg", "h"])):
            bt = tail(gates)
            for tau in TAUS:
                ops = [("sx", [q], [], dur("sx", [q]))]
                if tau > 0:
                    ops.append(("delay", [q], [], tau))
                ops += [(nm, [q], prm, dur(nm, [q])) for nm, prm in bt]
                out[(q, basis, tau)] = ops
    return out


def sx_tilted(op_list, qs, t):
    """The sx pulses on logical qubits qs with their axis turned by t against the x pulses: rz(-t), sx, rz(t)."""
    out = []
    for o in op_list:
        if o[0] == "sx" and o[1][0] in qs:
            out += [("rz", o[1], [-t], 0.0), o, ("rz", o[1], [t], 0.0)]
        else:
            out.append(o)
    return out


def turned(op_list, lq, phi):
    if phi == 0.0:
        return op_list
    clock = [0.0, 0.0, 0.0]
    for name, qs, prm, d in op_list:
        start = max(clock[q] for q in qs)
        for q in qs:
            clock[q] = start + d
    t_end = max(clock)
    pad = [("delay", [lq], [], t_end - clock[lq])] if t_end > clock[lq] else []
    return [("rz", [lq], [phi], 0.0)] + list(op_list) + pad + [("rz", [lq], [-phi], 0.0)]


def forbidden(basis, m):
    letters = [basis[0] if m[0] else "I", basis[1] if m[1] else "I"]
    return sum(c in "YZ" for c in letters) % 2 == 1 or sum(c == "Y" for c in letters) % 2 == 1


def main():
    data = json.loads(CF.DATA.read_text(encoding="utf-8"))
    tot = {}
    for nm, (h1, h2) in data["counts"].items():
        t = dict(h1)
        for k, v in h2.items():
            t[k] = t.get(k, 0) + v
        tot[nm] = t
    A_oct = R.readout_matrix_measured(tot)
    calib = R.calib_from(data["snapshot"])
    phys = R.PATH
    T1 = [calib["t1_us"][q] for q in phys]
    T2 = [calib["t2_us"][q] for q in phys]
    replay = CF.patched_replay()
    OPS = b_ops()
    bases = sorted(OPS)
    idx36 = [(b, m) for b in bases for m in rl.MASKS if forbidden(b, m)]
    # the statistics read the forbidden parities that carry no X letter: the start is an X eigenstate on every site and
    # YZ+ZY annihilates it, so <X0>, <X2>, <X1> stay near +-1 and every parity with an X partner nearly repeats the one
    # without it; their difference lives on rare outcomes, a near-singular direction of the covariance whose true
    # variance moves with any incoherent error the replay leaves out
    def no_x(b, m):
        return m[2] == 0 and not (m[0] and b[0] == "X") and not (m[1] and b[1] == "X")
    idx18 = [(b, m) for b, m in idx36 if m[2] == 0]
    idx = [(b, m) for b, m in idx36 if no_x(b, m)]

    # (0) April
    import reality_window_april_letters as ap
    april = json.loads(APRIL_JSON.read_text(encoding="utf-8"))["counts"]["pi2_even_nontruly"]
    meas = {(b[0], b[1]): rl.parities_probs(rl.counts_to_probs({"0" + k: v for k, v in cnt.items()})) for b, cnt in april.items()}
    A_apr = rl.assignment_product([ap.RO[0], ap.RO[1], (0.0, 0.0)])
    M = ap.measured()
    F = ap.forbidden("pi2_even_nontruly")

    def pooled(pars):
        out = {}
        for site, mask, pos in (("0", (1, 0, 0), 0), ("2", (0, 1, 0), 1)):
            for P in "XYZ":
                out[P + site] = float(np.mean([pars[b][mask] for b in pars if b[pos] == P]))
        return out
    keys6 = ("Y0", "Z0", "Y2", "Z2", "X0", "X2")
    mp = pooled(meas)
    say("(0) April 2026-04-30 YZ+ZY on [4, 5, 6] (q1 not read): single-site readings pooled over bases, and the sum of")
    say("    squared residuals over the ten readings the class forbids on (q0, q2), April's binomial sigma")
    say("    measured             " + "  ".join(f"{k} {mp[k]:+.3f}" for k in keys6))
    LSETS = (("no letter", None), ("October five Z letters", DP.FIVE), ("October crosstalk alone", DP.XT_FIT), ("October three families", DP.ALL_FIT))
    for sched in ("alap", "asap"):
        aops = ap.flown_ops(sched)
        for label, L in LSETS:
            pars = {b: rl.parities_probs(A_apr @ rl.probs_rho(replay(aops[("pi2_even_nontruly",) + b], ap.PATH, ap.T1, ap.T2, L))) for b in bases}
            acc = {}
            for (b0, b2), p in pars.items():
                acc.setdefault((b0, b2), []).append(p[(1, 1, 0)])
                acc.setdefault((b0, "I"), []).append(p[(1, 0, 0)])
                acc.setdefault(("I", b2), []).append(p[(0, 1, 0)])
            r = np.array([(M[("pi2_even_nontruly",) + k][0] - np.mean(acc[k])) / M[("pi2_even_nontruly",) + k][1] for k in F])
            q = pooled(pars)
            say(f"    {sched} {label:24s} sum {float(r @ r):6.1f}   " + "  ".join(f"{k} {q[k]:+.3f}" for k in keys6))

    # (1) B's baseline
    def parities(L, lq=1, phi=0.0, A=A_oct, tilt=None):
        def ops_of(b):
            o = OPS[b] if tilt is None else sx_tilted(OPS[b], tilt[0], tilt[1])
            return turned(o, lq, phi)
        return {b: rl.parities_probs(A @ rl.probs_rho(replay(ops_of(b), phys, T1, T2, L))) for b in bases}

    base = parities(None)

    def marg(A, bit):
        """P(read 0 | prepared 1) - P(read 1 | prepared 0) on one bit, averaged over the other bits' preparations."""
        p01 = np.mean([sum(A[j, k] for j in range(8) if j & bit) for k in range(8) if not k & bit])
        p10 = np.mean([sum(A[j, k] for j in range(8) if not j & bit) for k in range(8) if k & bit])
        return float(p10 - p01)
    xx = base[("X", "X")]
    say()
    say(f"(1) arm B: {len(bases)} bases, {len(idx36)} forbidden readings, {len(idx18)} without q1's bit, {len(idx)} with no X letter at all (the statistics' readings);")
    say("    baseline (no letter, October snapshot and RO matrix)")
    say(f"    XX basis: <X0> {xx[(1, 0, 0)]:+.3f}  <X2> {xx[(0, 1, 0)]:+.3f}  <X1> {xx[(0, 0, 1)]:+.3f}  <X0 X1 X2> {xx[(1, 1, 1)]:+.3f}")
    say(f"    single-site offsets in the Z bases: <Z0> {np.mean([base[b][(1, 0, 0)] for b in bases if b[0] == 'Z']):+.4f}"
        f"  <Z2> {np.mean([base[b][(0, 1, 0)] for b in bases if b[1] == 'Z']):+.4f}  (October RO asymmetries P(0|1) - P(1|0),"
        f" marginal over the other bits: q4 {marg(A_oct, 4):+.4f}, q6 {marg(A_oct, 2):+.4f})")

    def vec(par):
        return np.array([par[b][m] for b, m in idx])

    def shot_cov(par, n=SHOTS):
        covs = {b: R.cov_from_parities(par[b], n) for b in bases}
        S = np.zeros((len(idx), len(idx)))
        for i, (b1, m1) in enumerate(idx):
            for j, (b2, m2) in enumerate(idx):
                if b1 == b2:
                    S[i, j] = covs[b1][(m1, m2)]
        return S

    # (2) the RO matrix's own noise in y0: y_{c,m} = s_m . (A p_c), column j of A multinomial over RO_SHOTS
    ptrue = {b: rl.probs_rho(replay(OPS[b], phys, T1, T2, None)) for b in bases}
    Cj = [(np.diag(A_oct[:, j]) - np.outer(A_oct[:, j], A_oct[:, j])) / RO_SHOTS for j in range(8)]

    def ro_cov(pt):
        S = np.zeros((len(idx), len(idx)))
        for i, (b1, m1) in enumerate(idx):
            for k, (b2, m2) in enumerate(idx):
                S[i, k] = sum(pt[b1][j] * pt[b2][j] * (SIGN[m1] @ Cj[j] @ SIGN[m2]) for j in range(8))
        return S
    def shot_cov_on(sel, par, n=SHOTS):
        covs = {bb: R.cov_from_parities(par[bb], n) for bb in bases}
        return np.array([[covs[b1][(m1, m2)] if b1 == b2 else 0.0 for (b2, m2) in sel] for (b1, m1) in sel])

    def ro_cov_on(sel, pt):
        return np.array([[sum(pt[b1][j] * pt[b2][j] * (SIGN[m1] @ Cj[j] @ SIGN[m2]) for j in range(8)) for (b2, m2) in sel] for (b1, m1) in sel])
    Ss = shot_cov(base)
    Sr = ro_cov(ptrue)

    def pinv_rank(S):
        w, V = np.linalg.eigh(S)
        keep = w > 1e-9 * w.max()
        return (V[:, keep] / w[keep]) @ V[:, keep].T, int(keep.sum())
    Ws, rk = pinv_rank(Ss)
    Wt, _ = pinv_rank(Ss + Sr)
    say()
    say("(2) the day's RO matrix in a baseline comparison, B5(0)'s forbidden readings")
    say(f"    shot sigma median {np.median(np.sqrt(np.diag(Ss))):.4f}; RO-matrix sigma median {np.median(np.sqrt(np.diag(Sr))):.4f}, max {np.sqrt(np.diag(Sr)).max():.4f}")
    rng = np.random.default_rng(20261113)
    Lc = np.linalg.cholesky(Ss + Sr + 1e-15 * np.eye(len(idx)))
    draws = (Lc @ rng.standard_normal((len(idx), 20000))).T
    pt99 = CHI2.ppf(0.99, rk)
    fp_s = float(np.mean(np.einsum("ij,jk,ik->i", draws, Ws, draws) > pt99))
    # the same on all 36 forbidden readings, as a first draft read them
    S36s = shot_cov_on(idx36, base)
    pt36 = {b_: ptrue[b_] for b_ in bases}
    S36r = ro_cov_on(idx36, pt36)
    W36, rk36 = pinv_rank(S36s)
    L36 = np.linalg.cholesky(S36s + S36r + 1e-15 * np.eye(len(idx36)))
    d36 = (L36 @ rng.standard_normal((len(idx36), 20000))).T
    Lm36 = np.linalg.cholesky(S36s / 4 + S36r + 1e-15 * np.eye(len(idx36)))
    dm36 = (Lm36 @ rng.standard_normal((len(idx36), 20000))).T
    W36m, _ = pinv_rank(S36s / 4)
    say(f"    on all 36 forbidden readings with the shot covariance alone: W {float(np.mean(np.einsum('ij,jk,ik->i', d36, W36, d36) > CHI2.ppf(0.99, rk36))):.3f},"
        f" the turn-averaged mean {float(np.mean(np.einsum('ij,jk,ik->i', dm36, W36m, dm36) > CHI2.ppf(0.99, rk36))):.3f}")
    fp_t = float(np.mean(np.einsum("ij,jk,ik->i", draws, Wt, draws) > pt99))
    say(f"    W at its 99 % point ({pt99:.1f}, rank {rk}) under noise alone: false-positive rate {fp_s:.3f} with the shot covariance alone,"
        f" {fp_t:.3f} with the RO term added")
    Wm4, _ = pinv_rank(Ss / 4 + Sr)
    Wm4s, _ = pinv_rank(Ss / 4)
    Lm = np.linalg.cholesky(Ss / 4 + Sr + 1e-15 * np.eye(len(idx)))
    dm = (Lm @ rng.standard_normal((len(idx), 20000))).T
    say(f"    the turn-averaged mean against the baseline at its 99 % point: false-positive rate {float(np.mean(np.einsum('ij,jk,ik->i', dm, Wm4s, dm) > pt99)):.3f}"
        f" with the shot covariance alone, {float(np.mean(np.einsum('ij,jk,ik->i', dm, Wm4, dm) > pt99)):.3f} with the RO term")

    # (2b) a frame-blind true null the replay leaves out: a global depolarization of size eps scales every parity by
    # (1 - eps); the measured covariance follows the true parities while the statistics weigh with the baseline's
    def false_positive(sel, eps, ndraw=4000):
        ii = [k for k, (bb, mm) in enumerate(idx36) if (bb, mm) in sel]
        true = {bb: {mm: (1 - eps) * v for mm, v in base[bb].items()} for bb in bases}
        def cov_on(par):
            covs = {bb: R.cov_from_parities(par[bb], SHOTS) for bb in bases}
            return np.array([[covs[b1][(m1, m2)] if b1 == b2 else 0.0 for (b2, m2) in sel] for (b1, m1) in sel])
        Sb, St = cov_on(base), cov_on(true)
        Wb, r1 = pinv_rank(Sb)
        pt1 = CHI2.ppf(0.99, r1)
        Wd, r5 = pinv_rank(np.kron(np.array([[2.0, 1, 1], [1, 2, 1], [1, 1, 2]]), Sb))
        pt5 = CHI2.ppf(0.99, r5)
        L1 = np.linalg.cholesky(St + 1e-15 * np.eye(len(sel)))
        Lt5 = np.linalg.cholesky(np.kron(np.array([[2.0, 1, 1], [1, 2, 1], [1, 1, 2]]), St) + 1e-15 * np.eye(3 * len(sel)))
        x1 = (L1 @ rng.standard_normal((len(sel), ndraw))).T
        x5 = (Lt5 @ rng.standard_normal((3 * len(sel), ndraw))).T
        return (float(np.mean(np.einsum("ij,jk,ik->i", x1, Wb, x1) > pt1)), float(np.mean(np.einsum("ij,jk,ik->i", x5, Wd, x5) > pt5)))
    say("    a frame-blind null left out of the replay, global depolarization eps (shot covariance only): false-positive rates of")
    say("    W (the deviation's own fluctuation) and T5 at their 99 % points")
    for eps in (0.01, 0.02, 0.04, 0.1):
        f36 = false_positive(idx36, eps)
        f18 = false_positive(idx18, eps)
        f14 = false_positive(idx, eps)
        say(f"      eps {eps:.2f}: all 36 W {f36[0]:.3f} T5 {f36[1]:.3f};  18 without q1's bit W {f18[0]:.3f} T5 {f18[1]:.3f};"
            f"  {len(idx)} with no X W {f14[0]:.3f} T5 {f14[1]:.3f}")

    def ncW(L):
        d = vec(parities(L)) - vec(base)
        return float(d @ Wt @ d), float(np.abs(d).max())

    # (3) frame-blind letters
    say()
    say("(3) frame-blind letters: largest change under the turn on q5 (four phases) and q4, q6 (pi); W's noncentrality (RO term in)")
    BLIND = (("five Z letters (October fit)", DP.FIVE),
             ("detuning q4 +0.135 rad/us (21.5 kHz)", {("det", 4): 0.135}),
             ("detuning q6 +0.135 rad/us", {("det", 6): 0.135}),
             ("detuning q5 +0.1 rad/us", {("det", 5): 0.1}),
             ("sx and x axis tilt q4, q6 0.06 rad", {("sxt", 4): 0.06, ("sxt", 6): 0.06}),
             ("sx and x amplitude q4 +1 %", {("sxa", 4): 0.01}),
             ("sx and x amplitude q6 +1 %", {("sxa", 6): 0.01}),
             ("cz local phase 0.03, all", {("czp", p, q): 0.03 for p in rl.PAIRS for q in p}))
    for label, L in BLIND:
        p0 = parities(L)
        dmax = max(float(np.abs(vec(parities(L, lq, phi)) - vec(p0)).max())
                   for lq, phis in ((1, PHIS[1:]), (0, [math.pi]), (2, [math.pi])) for phi in phis)
        nc, big = ncW(L)
        say(f"    {label:38s} turn: max |D| {dmax:.1e};  W noncentrality {nc:8.1f}, largest shift {big:.1e}")

    tiltpar = parities(None, tilt=([0, 2], 0.06))
    dmax = max(float(np.abs(vec(parities(None, lq, phi, tilt=([0, 2], 0.06))) - vec(tiltpar)).max())
               for lq, phis in ((1, PHIS[1:]), (0, [math.pi]), (2, [math.pi])) for phi in phis)
    d = vec(tiltpar) - vec(base)
    say(f"    {'sx tilted 0.06 rad against x, q4 and q6':38s} turn: max |D| {dmax:.1e};  W noncentrality {float(d @ Wt @ d):8.1f}, largest shift {np.abs(d).max():.1e}")
    big = max(abs(v) for b in bases for m, v in rl.parities_probs(rl.probs_rho(
        replay(OPS[b], phys, [1e12] * 3, T2, None))).items() if forbidden(b, m))
    say(f"    no letter, T1 off, the snapshot's T2 dephasing interleaved through the compiled gates, ideal readout: largest forbidden reading {big:.1e}")

    # (4) turning letters
    def stacked(par0, pars):
        k, n = len(pars), len(idx)
        D = np.concatenate([vec(p) - vec(par0) for p in pars])
        Sig = np.kron(np.ones((k, k)), shot_cov(par0))
        for i, p in enumerate(pars):
            Sig[i * n:(i + 1) * n, i * n:(i + 1) * n] += shot_cov(p)
        W, r = pinv_rank(Sig)
        return float(D @ W @ D), r

    def power_line(label, T, rank):
        pt = CHI2.ppf(0.99, rank)
        pw = lambda lam: float(NCX2.sf(pt, rank, lam)) if lam > 0 else 0.01
        say(f"    {label}: noncentrality {T:.1f}, rank {rank}, 99 % point {pt:.1f}, power {pw(T):.3f}; with drift 1.27 {T / KAPPA_OCT:.1f},"
            f" power {pw(T / KAPPA_OCT):.3f}; with drift 1.96 {T / KAPPA_A:.1f}, power {pw(T / KAPPA_A):.3f}")
        p3 = CHI2.ppf(0.999, rank)
        pw3 = lambda lam: float(NCX2.sf(p3, rank, lam)) if lam > 0 else 0.001
        say(f"      at the 99.9 % point {p3:.1f}: power {pw3(T):.3f}; with drift 1.27 {pw3(T / KAPPA_OCT):.3f}; with drift 1.96 {pw3(T / KAPPA_A):.3f}")

    say()
    say("(4) turning letters (October's crosstalk fits)")
    for label, L in (("crosstalk alone (8)", DP.XT_FIT), ("all three families (17)", DP.ALL_FIT)):
        say(f"  == {label}")
        pars = [parities(L, 1, phi) for phi in PHIS]
        say(f"    largest change of a forbidden reading under the turn on q5: {max(float(np.abs(vec(p) - vec(pars[0])).max()) for p in pars[1:]):.4f}")
        T, rank = stacked(pars[0], pars[1:])
        power_line("T5, turn on q5 at four phases", T, rank)
        for lq in (0, 2):
            T, rank = stacked(pars[0], [parities(L, lq, math.pi)])
            power_line(f"T{R.PATH[lq]}, turn on q{R.PATH[lq]} at pi", T, rank)
        z2 = [float(np.mean([p[b][(0, 1, 0)] for b in bases if b[1] == "Z"])) for p in pars]
        say("    <Z2> pooled at the four turns " + ", ".join(f"{x:+.4f}" for x in z2) + f"; baseline {np.mean([base[b][(0, 1, 0)] for b in bases if b[1] == 'Z']):+.4f}")
        mean = np.mean([vec(p) for p in pars], axis=0)
        d = mean - vec(base)
        say(f"    turn-averaged mean minus baseline: largest {np.abs(d).max():.4f}; its noncentrality against shot/4 + RO {float(d @ Wm4 @ d):.1f} (99 % point {pt99:.1f})")
        nc, big = ncW(L)
        say(f"    W (B5(0) against the baseline) noncentrality {nc:.1f}")

    # (5) the bound on b's two pairs
    say()
    say("(5) the bound: T5 against eps on (5 -> 4) and (5 -> 6) at the October phases, the other pairs zero")
    for eps in (0.001, 0.002, 0.003, 0.004, 0.005):
        L = {("xtc", (5, 4)): eps * math.cos(-2.24), ("xts", (5, 4)): eps * math.sin(-2.24),
             ("xtc", (5, 6)): eps * math.cos(0.83), ("xts", (5, 6)): eps * math.sin(0.83)}
        pars = [parities(L, 1, phi) for phi in PHIS]
        T, rank = stacked(pars[0], pars[1:])
        power_line(f"eps {eps:.3f}", T, rank)

    # (6) arm Q
    say()
    say("(6) arm Q: one qubit at a time (the other two in |0>), sx, idle tau in {0, 1, 2, 4} us, read in X and Y; replay with a")
    say("    detuning on each qubit, RO-corrected readings,")
    say("    fitted per qubit by (delta, theta0, contrast, decay rate): <X> = C e^(-g tau) cos(theta0 + delta tau), <Y> = C e^(-g tau) sin(...)")
    QO = q_ops()
    Ainv = np.linalg.inv(A_oct)
    for dets in ((0.0, 0.0, 0.0), (0.05, -0.1, 0.2), (0.135, 0.0, -0.05)):
        L = {("det", phys[q]): dets[q] for q in range(3) if dets[q]}
        rd = {}
        for key, ops in QO.items():
            p = Ainv @ (A_oct @ rl.probs_rho(replay(ops, phys, T1, T2, L)))
            rd[key] = rl.parities_probs(p)
        fits = []
        for q, mask in ((0, (1, 0, 0)), (1, (0, 0, 1)), (2, (0, 1, 0))):
            x = np.array([rd[(q, "X", t)][mask] for t in TAUS])
            y = np.array([rd[(q, "Y", t)][mask] for t in TAUS])
            tt = np.array(TAUS)

            def res(v):
                d, th, C, g = v
                e = C * np.exp(-g * tt)
                return np.concatenate([e * np.cos(th + d * tt) - x, e * np.sin(th + d * tt) - y])
            th0 = math.atan2(y[0], x[0])
            r = least_squares(res, [0.0, th0, 1.0, 0.0], bounds=([-1.5, -10, 0, 0], [1.5, 10, 1.2, 1.0]))
            J = r.jac
            sig = np.sqrt(np.diag(np.linalg.inv(J.T @ J))) / math.sqrt(SHOTS) if np.linalg.matrix_rank(J) == 4 else [float("nan")] * 4
            fits.append(f"q{phys[q]} delta {r.x[0]:+.4f} (set {dets[q]:+.3f}, sigma {sig[0]:.4f}) theta0 {r.x[1]:+.3f}")
        say("    " + ";  ".join(fits))
    # noisy redraws: multinomial counts at SHOTS from the replay's distribution, RO-corrected, fitted the same way
    rngq = np.random.default_rng(7)
    L = {("det", 4): 0.05, ("det", 5): -0.1, ("det", 6): 0.2}
    probs = {key: A_oct @ rl.probs_rho(replay(ops, phys, T1, T2, L)) for key, ops in QO.items()}
    est = {0: [], 1: [], 2: []}
    for _ in range(200):
        rd = {key: rl.parities_probs(Ainv @ (rngq.multinomial(SHOTS, np.clip(p, 0, None) / np.clip(p, 0, None).sum()) / SHOTS)) for key, p in probs.items()}
        for q, mask in ((0, (1, 0, 0)), (1, (0, 0, 1)), (2, (0, 1, 0))):
            x = np.array([rd[(q, "X", t)][mask] for t in TAUS]); y = np.array([rd[(q, "Y", t)][mask] for t in TAUS]); tt = np.array(TAUS)
            def res(v):
                d, th, C, g = v
                e = C * np.exp(-g * tt)
                return np.concatenate([e * np.cos(th + d * tt) - x, e * np.sin(th + d * tt) - y])
            est[q].append(least_squares(res, [0.0, math.atan2(y[0], x[0]), 1.0, 0.0], bounds=([-1.5, -10, 0, 0], [1.5, 10, 1.2, 1.0])).x[0])
    say("    200 noisy redraws at 8192 shots, detunings set +0.05, -0.10, +0.20: " + ";  ".join(
        f"q{phys[q]} mean {np.mean(est[q]):+.4f} sd {np.std(est[q]):.4f}" for q in range(3)))
    # (7) the drift factor estimated from the data: T ~ kappa chi2_rank, kappa_hat = max(1, kappa chi2_nu / nu) from the
    # nu half-difference terms (84 for B's 14 readings over six variants), T / kappa_hat against the chi2 points
    say()
    say("(7) a drift factor estimated from the data, under real drift: false-positive rate of T / kappa_hat at the 99 % and 99.9 % points")
    terms_all, terms_a = [], []
    for nm, (h1, h2) in data["counts"].items():
        if R.arm_of(nm) not in R.FIT_ARMS:
            continue
        bb = nm.split("_")[1]
        p1, c1 = rl.parities_counts(h1)
        p2, c2 = rl.parities_counts(h2)
        for mm in rl.MASKS:
            if rl.odd_y((bb[0], bb[1]), mm) and c1[(mm, mm)] + c2[(mm, mm)] > 0:
                t = (p1[mm] - p2[mm]) ** 2 / (c1[(mm, mm)] + c2[(mm, mm)])
                terms_all.append(t)
                if R.arm_of(nm) == "A":
                    terms_a.append(t)
    say(f"    October's drift factor from its two halves: over all fit arms {np.mean(terms_all):.4f} ({len(terms_all)} odd-Y parities),"
        f" over arm A alone {np.mean(terms_a):.4f} ({len(terms_a)})")
    rng7 = np.random.default_rng(1)
    NN = 200000
    for rank in (42, 14):
        for kap in (1.0, 1.3, 2.0):
            T = kap * rng7.chisquare(rank, NN)
            kh = np.maximum(1.0, kap * rng7.chisquare(84, NN) / 84)
            say(f"    rank {rank}, kappa {kap}: {np.mean(T / kh > CHI2.ppf(0.99, rank)):.4f} and {np.mean(T / kh > CHI2.ppf(0.999, rank)):.4f}")
    exc = np.ones(14); exc[:6] = 2.0
    T = (rng7.standard_normal((NN, 42)) ** 2 * np.tile(exc, 3)).sum(1)
    kh = np.maximum(1.0, (rng7.standard_normal((NN, 84)) ** 2 * np.tile(exc, 6)).mean(1))
    say(f"    rank 42, drift doubling the variance of 6 of the 14 readings: {np.mean(T / kh > CHI2.ppf(0.99, 42)):.4f} and {np.mean(T / kh > CHI2.ppf(0.999, 42)):.4f}")
    # the same with drift concentrated on a few readings (one qubit's frequency wandering), and three ways to estimate
    # it: one global factor, one per reading (its six half-difference terms), and the larger of the two per reading
    say("    drift on a few readings; the frame statistic over R readings x 3 differences, six half terms per reading;")
    say("    false-positive rate at the 99 % / 99.9 % points for kappa_hat global, per reading, and max(global, per reading)")
    for R_, label in ((14, "B, 14 readings, rank 42"), (36, "A, 36 readings, rank 108")):
        p99, p999 = CHI2.ppf(0.99, 3 * R_), CHI2.ppf(0.999, 3 * R_)
        for nhot, f in ((0, 1.0), (R_, 2.0), (3, 3.0), (3, 5.0), (3, 10.0), (1, 10.0)):
            k = np.ones(R_); k[:nhot] = f
            terms = k[None, :] * rng7.chisquare(3, (NN, R_))
            half = k[None, :] * rng7.chisquare(6, (NN, R_)) / 6
            kg = np.maximum(1.0, half.mean(1))
            ki = np.maximum(1.0, half)
            kc = np.maximum(ki, kg[:, None])
            row = []
            for T in (terms.sum(1) / kg, (terms / ki).sum(1), (terms / kc).sum(1)):
                row.append(f"{np.mean(T > p99):.4f}/{np.mean(T > p999):.4f}")
            say(f"    {label}, {nhot} readings at {f}x: " + "  ".join(row))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(chr(10).join(LINES) + chr(10), encoding="utf-8", newline=chr(10))


if __name__ == "__main__":
    main()
