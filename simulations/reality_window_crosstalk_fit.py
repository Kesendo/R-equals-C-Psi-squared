"""The reality-window flight: the drive-crosstalk letter and the spectator's Stark letter fitted to the 288 odd-Y readings
through the patched replay (exploratory, not pre-registered; stage 0 of experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md).

The letter: while logical qubit k is pulsed (sx/x, rotation rate `rate` about the X axis of k's circuit frame, whose
physical drive phase is the accumulated circuit rz TH[k]), its neighbour j receives the field
    eps * rate / 2 * (cos phi X_j + sin phi Y_j),   phi = TH[j] - TH[k] + delta,
(the handedness: IBM's rz(theta) deferred past a pulse turns the pulse's axis by -theta, so a field at k's drive
phase sits in the replay's frame for j at TH[j] - TH[k]; gated in simulations/drive_frame_handedness_gate.py)
one (eps, delta) per ordered neighbour pair (k -> j) of the path [4, 5, 6], parametrised as (eps cos delta, eps sin delta)
= letters ('xtc', (k, j)), ('xts', (k, j)) on PHYSICAL qubits; and the spectator's Z letter during the neighbour's pulses
(an AC Stark shift), ('stx', (k, j)) in rad/us. Fitted alone and beside the five letters of the post-flight fits. Gate: with no letter the patched replay equals the frozen replay exactly (baseline path).
Output: simulations/results/reality_window/crosstalk_fit.txt (about half an hour).
"""
import inspect
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as CHI2

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "simulations" / "flight"))
import run_reality_window as R          # noqa: E402
import reality_window_letters as rl     # noqa: E402

DATA = REPO / "data" / "ibm_reality_window_october2026" / "reality_window_ibm_marrakesh_20261029_053203.json"
OUT = REPO / "simulations" / "results" / "reality_window" / "crosstalk_fit.txt"
PAIRS = [(5, 4), (5, 6), (4, 5), (6, 5)]                       # (driver k, spectator j), physical
XT = [("xtc", p) for p in PAIRS] + [("xts", p) for p in PAIRS]
FIVE = [("det", 4), ("det", 5), ("det", 6), ("zz", (4, 5)), ("zz", (5, 6))]
FIVE_OPT = [0.001, 0.201, 0.180, 0.002, -0.111]

EDITS = [
    ("    k = 0\n    for j, t in enumerate(bounds):", "    TH = [0.0, 0.0, 0.0]\n    k = 0\n    for j, t in enumerate(bounds):"),
    ("_, _, kind, data = inst[k]", "_, ordr, kind, data = inst[k]"),
    ("                unitary(on(data[0], rz(data[1])))",
     "                unitary(on(data[0], rz(data[1])))\n"
     "                if ordr == int(ordr):\n"
     "                    TH[data[0]] += data[1]"),
    ("    L = L or {}\n",
     "    L = L or {}\n"
     "    XT_ON = any(k_[0] in ('xtc', 'xts', 'stx') and v_ for k_, v_ in L.items())\n"),
    ("        key = (active, round(h, 12))\n",
     "        key = (active, round(h, 12)) + ((tuple(round(t_, 9) for t_ in TH),) if XT_ON else ())\n"),
    ("                H = H + rate / 2 * on(q, math.cos(tilt[q]) * XP + math.sin(tilt[q]) * YP)\n",
     "                H = H + rate / 2 * on(q, math.cos(tilt[q]) * XP + math.sin(tilt[q]) * YP)\n"
     "                for j_ in range(3):\n"
     "                    if j_ == q:\n"
     "                        continue\n"
     "                    xc_ = L.get(('xtc', (phys[q], phys[j_])), 0.0)\n"
     "                    xs_ = L.get(('xts', (phys[q], phys[j_])), 0.0)\n"
     "                    if xc_ or xs_:\n"
     "                        ph_ = TH[j_] - TH[q] + math.atan2(xs_, xc_)\n"
     "                        H = H + math.hypot(xc_, xs_) * rate / 2 * on(j_, math.cos(ph_) * XP + math.sin(ph_) * YP)\n"
     "                    sz_ = L.get(('stx', (phys[q], phys[j_])), 0.0)\n"
     "                    if sz_:\n"
     "                        H = H + sz_ / 2 * ZOP[j_]\n"),
]


def patched_replay():
    src = inspect.getsource(rl.replay)
    for old, new in EDITS:
        if src.count(old) != 1:
            raise RuntimeError("patch point not found once: " + repr(old[:60]))
        src = src.replace(old, new)
    ns = dict(rl.__dict__)
    ns.update(math=math)
    exec(src, ns)
    return ns["replay"]


def main():
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    data = json.loads(DATA.read_text(encoding="utf-8"))
    tot = {}
    for nm, (h1, h2) in data["counts"].items():
        t = dict(h1)
        for k, v in h2.items():
            t[k] = t.get(k, 0) + v
        tot[nm] = t
    inflate = R.drift_inflation(data["counts"])[0]
    A = R.readout_matrix_measured(tot)
    calib = R.calib_from(data["snapshot"])
    circs = R.load_circuits()
    by = {c.name: c for c in circs}
    par, nsh = {}, {}
    for nm, t in tot.items():
        if not nm.startswith("RO_"):
            par[nm] = rl.parities_counts(t)[0]
            nsh[nm] = sum(t.values())
    names = sorted(n for n in par if R.arm_of(n) in R.FIT_ARMS)
    idx = rl.reading_index(names)
    ops = {nm: R.extract(by[nm]) for nm in names}
    meas = np.array([par[n][m] for n, m in idx])
    replay = patched_replay()

    def model(L=None):
        out = {}
        for nm in names:
            arm = R.arm_of(nm)
            phys = R.phys_of(arm)
            rho = replay(ops[nm], phys, [calib["t1_us"][q] for q in phys], [calib["t2_us"][q] for q in phys], L)
            out[nm] = rl.parities_probs(R.arm_readout(A, arm) @ rl.probs_rho(rho))
        return np.array([out[n][m] for n, m in idx]), out

    base, base_all = model()
    gap = float(np.abs(base - R.model_vector(circs, idx, calib, A)).max())
    say("patched replay with no letter against the frozen replay: max difference " + str(gap))
    assert gap == 0.0
    basepar = {nm: R.model_parities(by[nm], calib, A) for nm in names}
    covs = {nm: R.cov_from_parities(basepar[nm], nsh[nm]) for nm in names}
    W, rank, _ = R.weight_matrix(idx, covs, inflate, None)
    ev, U = np.linalg.eigh(np.linalg.pinv(W))
    keep = ev > ev.max() * 1e-10
    Wa = (U[:, keep] / np.sqrt(ev[keep])).T

    def chi(v):
        return float((meas - v) @ W @ (meas - v))

    def s_of(out, arm):
        nm = [n for n in names if R.arm_of(n) == arm and n.split("_")[1][:2] == "YZ"]
        return float(out[nm[0]][(1, 1, 0)]) if nm else float("nan")

    say("chi2 with no letter %.1f over rank %d (99 %% point %.1f)" % (chi(base), rank, CHI2.ppf(0.99, rank)))
    say("measured <Y0 Z2> raw: A %+.4f  D3 %+.4f  R %+.4f  D6 %+.4f" % tuple(s_of(par, a) for a in ("A", "D3", "R", "D6")))

    def report(label, L, status=None):
        v, out = model(L)
        d = np.abs(v - base)
        say("  %s: chi2 %.1f%s; largest odd-Y shift %.4f (median %.4f); model <Y0 Z2> raw A %+.4f D3 %+.4f R %+.4f D6 %+.4f"
            % (label, chi(v), (" (status %s)" % status) if status is not None else "", float(d.max()),
               float(np.median(d)), s_of(out, "A"), s_of(out, "D3"), s_of(out, "R"), s_of(out, "D6")))

    say("")
    say("1. single-pair scan: eps = 0.01 at delta = 0 and pi/2 for each ordered pair")
    for p in PAIRS:
        report("pair %s eps 0.01 delta 0" % (p,), {("xtc", p): 0.01})
        report("pair %s eps 0.01 delta pi/2" % (p,), {("xts", p): 0.01})

    def fit(sub, starts, scale, nfev):
        def fun(x):
            return Wa @ (meas - model(dict(zip(sub, x)))[0])
        best = None
        for st in starts:
            r = least_squares(fun, st, x_scale=scale, max_nfev=nfev)
            c = float(r.fun @ r.fun)
            say("    start %s -> chi2 %.1f (status %s, nfev %d)" % (np.round(st, 3).tolist(), c, r.status, r.nfev))
            if best is None or c < best[1]:
                best = (r.x, c, r.status)
        return best

    def describe(x, sub):
        parts = []
        for p in PAIRS:
            c = x[sub.index(("xtc", p))]
            s = x[sub.index(("xts", p))]
            parts.append("%s: eps %.4f delta %+.2f" % (p, math.hypot(c, s), math.atan2(s, c)))
        return "  ".join(parts)

    say("")
    say("2. the eight crosstalk parameters alone")
    x, c, st = fit(XT, [[0.0] * 8, [0.01] * 8, [-0.01] * 8], [0.01] * 8, 150)
    say("  best: " + describe(x, XT))
    report("  at the fit", dict(zip(XT, x)), st)
    say("  (99 %% point at rank-8 %.1f)" % CHI2.ppf(0.99, rank - 8))

    say("")
    say("3. the eight beside the five letters")
    sub = XT + FIVE
    x2, c2, st2 = fit(sub, [list(x) + FIVE_OPT, [0.0] * 8 + FIVE_OPT], [0.01] * 8 + [0.1] * 5, 200)
    say("  best: " + describe(x2, sub) + "  five " + str(np.round(x2[8:], 3).tolist()))
    report("  at the fit", dict(zip(sub, x2)), st2)
    say("  (99 %% point at rank-13 %.1f)" % CHI2.ppf(0.99, rank - 13))

    say("")
    say("4. the spectator's Z letter during the neighbour's pulses (AC Stark), rad/us, four ordered pairs alone")
    ST = [("stx", p) for p in PAIRS]
    x3, c3, st3 = fit(ST, [[0.0] * 4, [0.3] * 4, [-0.3] * 4], [0.1] * 4, 150)
    say("  best: " + "  ".join("%s: %+.4f rad/us" % (p, v) for p, v in zip(PAIRS, x3)))
    report("  at the fit", dict(zip(ST, x3)), st3)
    say("  (99 %% point at rank-4 %.1f)" % CHI2.ppf(0.99, rank - 4))
    say("")
    say("5. the spectator's Z letter beside the five letters")
    sub5 = ST + FIVE
    x4, c4, st4 = fit(sub5, [list(x3) + FIVE_OPT, [0.0] * 4 + FIVE_OPT], [0.1] * 4 + [0.1] * 5, 200)
    say("  best: " + "  ".join("%s: %+.4f" % (p, v) for p, v in zip(PAIRS, x4[:4])) + "  five " + str(np.round(x4[4:], 3).tolist()))
    report("  at the fit", dict(zip(sub5, x4)), st4)
    say("  (99 %% point at rank-9 %.1f)" % CHI2.ppf(0.99, rank - 9))
    say("")
    say("6. all three families: eight crosstalk + four Stark + five letters")
    sub6 = XT + ST + FIVE
    x5, c5, st5 = fit(sub6, [list(x2[:8]) + list(x4[:4]) + list(x4[4:]), [0.0] * 12 + FIVE_OPT], [0.01] * 8 + [0.1] * 4 + [0.1] * 5, 250)
    say("  best: " + describe(x5, sub6) + "  stark " + str(np.round(x5[8:12], 4).tolist()) + "  five " + str(np.round(x5[12:], 3).tolist()))
    report("  at the fit", dict(zip(sub6, x5)), st5)
    say("  (99 %% point at rank-17 %.1f)" % CHI2.ppf(0.99, rank - 17))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
