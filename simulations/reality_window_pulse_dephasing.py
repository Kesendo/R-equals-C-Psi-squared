"""The reality-window flight: the cube's incoherent candidate for the odd-Y excess, tested through the replay
(exploratory, not pre-registered; experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md section 7, Read on the cube).

The cube says an incoherent process pays the same rate for its odd-Y arrow (a dephasing in a tilted frame
carries X <-> Y at sin(2 theta) times its rate) and for the contrast it removes from the real sector. The real
sector lost 15 to 22 % of its contrast against the baseline on the day; if that loss is incoherent dephasing
DURING THE DRIVE (a pulse-induced dephasing of the driven qubit, where the frame is tilted), the same rate
could carry an odd-Y weight of the measured size. Tested here: an extra Z dephasing on the driven qubit during
its own sx/x pulses, coherence decay rate G (1/us), global or per physical qubit; fitted to the 288 odd-Y
readings (post-flight weights); the even-sector contrast slope per arm is printed for the baseline and the fit,
as the handshake the cube asks for.  Also, as the control the cube predicts to be weaker: the same extra
dephasing at ALL times (idle included), which spends most of its time in frames where the carried axis is a
Pauli axis.  And the cube's stronger incoherent candidate: an extra amplitude damping (sigma^- at rate D, no
extra dephasing) of the driven qubit during its own pulses, whose arrow I -> (carried Z) is fed by the trace at
the full rate D wherever the frame puts Z on the Y axis.
Gate: with no letter the patched replay must equal the frozen replay exactly (the baseline path, not the new terms).
Output: simulations/results/reality_window/pulse_dephasing.txt (about a quarter of an hour).
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
OUT = REPO / "simulations" / "results" / "reality_window" / "pulse_dephasing.txt"
FIVE = [("det", 4), ("det", 5), ("det", 6), ("zz", (4, 5)), ("zz", (5, 6))]
FIVE_OPT = [0.001, 0.201, 0.180, 0.002, -0.111]

PULSE_OLD = "        for _ in range(n):\n            unitary(U)\n            decay_all(h)\n"
PULSE_NEW = PULSE_OLD + (
    "            for q_, _r in active:\n"
    "                g_ = L.get(('pdp', phys[q_]), 0.0) + L.get(('pdg',), 0.0)\n"
    "                if g_:\n"
    "                    p_ = (1 - math.exp(-g_ * h)) / 2\n"
    "                    rho = (1 - p_) * rho + p_ * ZOP[q_] @ rho @ ZOP[q_]\n"
    "                d_ = L.get(('pdd', phys[q_]), 0.0) + L.get(('pddg',), 0.0)\n"
    "                if d_:\n"
    "                    rho = decay(rho, q_, h, [1 / d_] * 3, [2 / d_] * 3)\n")
ALL_OLD = "    def decay_all(h):\n        nonlocal rho\n        for q in range(3):\n            rho = decay(rho, q, h, T1, T2)\n"
ALL_NEW = ALL_OLD + (
    "            g_ = L.get(('adp', phys[q]), 0.0) + L.get(('adg',), 0.0)\n"
    "            if g_:\n"
    "                p_ = (1 - math.exp(-g_ * h)) / 2\n"
    "                rho = (1 - p_) * rho + p_ * ZOP[q] @ rho @ ZOP[q]\n")


def patched_replay():
    src = inspect.getsource(rl.replay)
    for old, new in ((PULSE_OLD, PULSE_NEW), (ALL_OLD, ALL_NEW)):
        if src.count(old) != 1:
            raise RuntimeError("patch point not found once: " + repr(old[:50]))
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

    def basis(nm):
        b = nm.split("_")[1]
        return (b[0], b[1])

    even_idx = [(nm, m) for nm in names for m in rl.MASKS if not rl.odd_y(basis(nm), m)]
    ops = {nm: R.extract(by[nm]) for nm in names}
    meas = np.array([par[n][m] for n, m in idx])
    meas_even = np.array([par[n][m] for n, m in even_idx])
    replay = patched_replay()

    def model_all(L=None):
        out = {}
        for nm in names:
            arm = R.arm_of(nm)
            phys = R.phys_of(arm)
            rho = replay(ops[nm], phys, [calib["t1_us"][q] for q in phys], [calib["t2_us"][q] for q in phys], L)
            out[nm] = rl.parities_probs(R.arm_readout(A, arm) @ rl.probs_rho(rho))
        return out

    def split(out):
        return np.array([out[n][m] for n, m in idx]), np.array([out[n][m] for n, m in even_idx])

    base_all = model_all()
    base, base_even = split(base_all)
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

    say("chi2 with no letter %.1f over rank %d (99 %% point %.1f); idle 1/T2: " % (chi(base), rank, CHI2.ppf(0.99, rank))
        + "  ".join("q%d %.4f/us" % (q, 1 / calib["t2_us"][q]) for q in (4, 5, 6)))

    def slopes(mod_even):
        """Per arm: slope of measured even readings on the model's, over readings with |model| > 0.1."""
        res = {}
        for arm in R.FIT_ARMS:
            sel = [i for i, (nm, m) in enumerate(even_idx) if R.arm_of(nm) == arm and abs(mod_even[i]) > 0.1]
            x = mod_even[sel]
            y = meas_even[sel]
            res[arm] = float(x @ y / (x @ x))
        return res

    def s_of(out, arm):
        nm = [n for n in names if R.arm_of(n) == arm and basis(n) == ("Y", "Z")]
        return float(out[nm[0]][(1, 1, 0)]) if nm else float("nan")

    def report(label, L, status=None):
        out = model_all(L)
        v, ve = split(out)
        sl = slopes(ve)
        d = np.abs(v - base)
        say("  %s: chi2 %.1f%s; largest shift of an odd-Y reading against the baseline %.4f (median %.4f)"
            % (label, chi(v), (" (status %s)" % status) if status is not None else "", float(d.max()), float(np.median(d))))
        say("    even slope per arm: " + "  ".join("%s %.3f" % (a, sl[a]) for a in R.FIT_ARMS))
        say("    model <Y0 Z2> raw in arm A %+.4f, D3 %+.4f, R %+.4f" % (s_of(out, "A"), s_of(out, "D3"), s_of(out, "R")))

    say("")
    say("baseline even slopes (the contrast the real sector kept against the baseline):")
    sl = slopes(base_even)
    say("    " + "  ".join("%s %.3f" % (a, sl[a]) for a in R.FIT_ARMS))
    say("    measured <Y0 Z2> raw in arm A %+.4f, D3 %+.4f, R %+.4f; baseline A %+.4f"
        % (s_of(par, "A"), s_of(par, "D3"), s_of(par, "R"), s_of(base_all, "A")))

    def fit(sub, starts, scale):
        def fun(x):
            return Wa @ (meas - split(model_all(dict(zip(sub, x))))[0])
        best = None
        for st in starts:
            r = least_squares(fun, st, x_scale=scale, max_nfev=120)
            c = float(r.fun @ r.fun)
            if best is None or c < best[1]:
                best = (r.x, c, r.status)
        return best

    say("")
    say("1. extra dephasing on the driven qubit during its own pulses, coherence rate G (1/us)")
    for g in (0.05, 0.2, 0.5, 1.0):
        report("global G = %s" % g, {("pdg",): g})
    x, c, st = fit([("pdg",)], [[0.1], [0.5]], [0.2])
    say("  fit global G: %+.4f/us" % x[0])
    report("  at the fit", {("pdg",): x[0]}, st)
    sub = [("pdp", 4), ("pdp", 5), ("pdp", 6)]
    x, c, st = fit(sub, [[0.1, 0.1, 0.1], [0.5, 0.5, 0.5]], [0.2] * 3)
    say("  fit per qubit G4, G5, G6: " + str(np.round(x, 4).tolist()))
    report("  at the fit", dict(zip(sub, x)), st)
    sub5 = sub + FIVE
    x, c, st = fit(sub5, [[0.1, 0.1, 0.1] + FIVE_OPT, [0.3, 0.3, 0.3] + [0.0] * 5], [0.2] * 3 + [0.1] * 5)
    say("  fit per qubit G + five letters: G " + str(np.round(x[:3], 4).tolist()) + "  five " + str(np.round(x[3:], 3).tolist()))
    report("  at the fit", dict(zip(sub5, x)), st)
    say("  (99 %% points: rank-1 %.1f, rank-3 %.1f, rank-8 %.1f)" % (CHI2.ppf(0.99, rank - 1), CHI2.ppf(0.99, rank - 3), CHI2.ppf(0.99, rank - 8)))

    say("")
    say("2. control: the same extra dephasing at all times (idle included)")
    for g in (0.05, 0.2, 0.5):
        report("global G_all = %s" % g, {("adg",): g})
    x, c, st = fit([("adg",)], [[0.1], [0.3]], [0.2])
    say("  fit global G_all: %+.4f/us" % x[0])
    report("  at the fit", {("adg",): x[0]}, st)

    say("")
    say("3. extra amplitude damping of the driven qubit during its own pulses, rate D (1/us), no extra dephasing")
    for g in (0.05, 0.2, 0.5):
        report("global D = %s" % g, {("pddg",): g})
    x, c, st = fit([("pddg",)], [[0.1], [0.5]], [0.2])
    say("  fit global D: %+.4f/us" % x[0])
    report("  at the fit", {("pddg",): x[0]}, st)
    subd = [("pdd", 4), ("pdd", 5), ("pdd", 6)]
    x, c, st = fit(subd, [[0.1, 0.1, 0.1], [0.5, 0.5, 0.5]], [0.2] * 3)
    say("  fit per qubit D4, D5, D6: " + str(np.round(x, 4).tolist()))
    report("  at the fit", dict(zip(subd, x)), st)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
