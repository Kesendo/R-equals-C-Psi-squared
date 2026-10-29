"""The reality-window flight: three candidate carriers of the odd-Y excess, tested through the replay
(exploratory, not pre-registered; experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md section 7).

Reads the committed flight record (ibm_marrakesh job dbhcohht9uhs73c9vpj0) through the frozen replay of
simulations/flight/reality_window_letters.py, with the readout matrix the RO arm measured, against the 288 odd-Y
readings of the fit arms (weights as in reality_window_post_flight.py: the drift factor, no bootstrap term), and
prints:
  1. a residual excitation at the start, rho0 = (x) diag(1 - p_q, p_q) on q4, q5, q6, read together with the RO
     arm: its x gates flip the prepared bits from the same start, so A_RO = A_true . M(p) with M the symmetric flip,
     and the model carries A_true = A_RO . M(p)^-1. Each p_q is bounded by the P(1|0) the RO arm measured on that
     qubit in its all-zero preparation, an approximate physical bound (cells of A_RO that read zero make A_true
     slightly negative at every p > 0; the smallest entry is printed). Printed: a scan in p4 with chi2 and the
     excitation's noncentrality d W d (d its shift of the readings: what it would add to chi2 on readings that held
     nothing else), and a fit of all three p within those bounds;
  2. quasi-static detuning: the replay averaged over the eight sign choices of +-sigma on the three qubits (a
     symmetric two-point distribution per qubit), about zero and about the three-detuning optimum of
     reality_window_post_flight.py (THREE_OPT, from its post_flight.txt);
  3. a residual exchange at every cz, U = expm(-i eps G), G = e^{i phi} s+_a s-_b + h.c. (a the lower physical
     qubit), with phi0 fixed either in the hardware frame (phi = phi0 - (theta_a - theta_b), theta the virtual rz
     accumulated so far) or in the replay's frame (phi = phi0); fitted alone and beside the three detunings and
     static ZZ on both pairs, with the post-flight script's budget of 200 evaluations, from five starts (zero and
     +-0.05 in the eps cos phi0 components) and, beside the five letters, a sixth at its five-letter optimum
     (FIVE_OPT); the best start's least-squares status is printed. Left out: the phase
     (omega_a - omega_b) t that a static coupling carries in the rotating frame, which moves with the absolute time
     of each cz (ibm_marrakesh does not publish its qubit frequencies), and any virtual phase the calibrated cz
     carries inside itself, which would move the hardware frame at every cz by an amount not known here.
The frozen modules are not changed: items 1 and 3 run a copy of the replay with the start state and the exchange
patched in. Two gates: with no letter the copy must equal the frozen replay exactly (this guards the patches of the
baseline path, not the new terms), and the exchange generator must equal (XX+YY)/2 at phase 1 and (XY-YX)/2 at
phase i, exactly.
Output: simulations/results/reality_window/carriers.txt (about an hour).
"""
import inspect
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares
from scipy.stats import chi2 as CHI2

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "simulations" / "flight"))
import run_reality_window as R          # noqa: E402
import reality_window_letters as rl     # noqa: E402

DATA = REPO / "data" / "ibm_reality_window_october2026" / "reality_window_ibm_marrakesh_20261029_053203.json"
OUT = REPO / "simulations" / "results" / "reality_window" / "carriers.txt"
PAIRS = [(4, 5), (5, 6)]
FIVE = [("det", 4), ("det", 5), ("det", 6), ("zz", (4, 5)), ("zz", (5, 6))]
EX = [("exc", p) for p in PAIRS] + [("exs", p) for p in PAIRS]     # eps cos phi0, eps sin phi0
SCALE = {"exc": 0.05, "exs": 0.05, "det": 0.1, "zz": 0.1}
THREE_OPT = {("det", 4): -0.012, ("det", 5): 0.204, ("det", 6): 0.176}     # post_flight.txt, the 3-letter fit
FIVE_OPT = [0.001, 0.201, 0.180, 0.002, -0.111]                            # post_flight.txt, the 5-letter fit
SP = np.array([[0, 1], [0, 0]], dtype=complex)                             # sigma+ = |0><1|


def generator(la, lb, phi=None, phase=None):
    """G = z s+_a s-_b + h.c. on logical qubits la, lb of the replay's three, z = e^{i phi} (or z = phase given)."""
    z = phase if phase is not None else np.exp(1j * phi)
    G = z * rl.on(la, SP) @ rl.on(lb, SP.T)
    return G + G.conj().T


def patched_replay():
    """The frozen replay with a start state RHO0[0] and the exchange letter; FRAME[0] in ('nco', 'replay')."""
    src = inspect.getsource(rl.replay)
    edits = [
        ("rho = np.zeros((8, 8), dtype=complex); rho[0, 0] = 1", "rho = RHO0[0].copy()"),
        ("    k = 0\n    for j, t in enumerate(bounds):", "    TH = [0.0, 0.0, 0.0]\n    k = 0\n    for j, t in enumerate(bounds):"),
        ("_, _, kind, data = inst[k]", "_, ordr, kind, data = inst[k]"),
        ("                unitary(on(data[0], rz(data[1])))",
         "                unitary(on(data[0], rz(data[1])))\n"
         "                if ordr == int(ordr):\n"            # a circuit rz is a frame change; an sxp letter is not
         "                    TH[data[0]] += data[1]"),
        ("                U = czm(a, b)\n",
         "                U = czm(a, b)\n"
         "                ec, es = L.get(('exc', (pa, pb)), 0.0), L.get(('exs', (pa, pb)), 0.0)\n"
         "                if ec or es:\n"
         "                    la, lb = (a, b) if phys[a] == pa else (b, a)\n"
         "                    phi = math.atan2(es, ec)\n"
         "                    if FRAME[0] == 'nco':\n"
         "                        phi -= TH[la] - TH[lb]\n"
         "                    U = expm(-1j * math.hypot(ec, es) * GEN(la, lb, phi)) @ U\n"),
    ]
    for old, new in edits:
        if src.count(old) != 1:
            raise RuntimeError(f"replay patch point not found once: {old[:40]!r}")
        src = src.replace(old, new)
    ns = dict(rl.__dict__)
    ns.update(RHO0=[None], FRAME=["nco"], GEN=generator, math=math, expm=expm)
    exec(src, ns)
    return ns["replay"], ns["RHO0"], ns["FRAME"]


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
    replay, RHO0, FRAME = patched_replay()

    def flips(p):                       # measured bit order (q4, q6, q5)
        M = np.array([[1.0]])
        for q in (4, 6, 5):
            M = np.kron(M, np.array([[1 - p[q], p[q]], [p[q], 1 - p[q]]]))
        return M

    def model(L=None, p=None):
        p = p if p is not None else {4: 0.0, 5: 0.0, 6: 0.0}
        At = A @ np.linalg.inv(flips(p))
        out = {}
        for nm in names:
            arm = R.arm_of(nm)
            phys = R.phys_of(arm)
            r0 = np.array([[1.0]])
            for l in range(3):
                r0 = np.kron(r0, np.diag([1 - p[phys[l]], p[phys[l]]]))
            RHO0[0] = r0.astype(complex)
            rho = replay(ops[nm], phys, [calib["t1_us"][q] for q in phys], [calib["t2_us"][q] for q in phys], L)
            out[nm] = rl.parities_probs(R.arm_readout(At, arm) @ rl.probs_rho(rho))
        return np.array([out[n][m] for n, m in idx])

    base = model()
    gap = float(np.abs(base - R.model_vector(circs, idx, calib, A)).max())
    say(f"patched replay with no letter against the frozen replay: max difference {gap}")
    assert gap == 0.0, "the patched replay departs from the frozen one with no letter"
    XX = rl.on(0, rl.XP) @ rl.on(1, rl.XP)
    YY = rl.on(0, rl.YP) @ rl.on(1, rl.YP)
    XY = rl.on(0, rl.XP) @ rl.on(1, rl.YP)
    YX = rl.on(0, rl.YP) @ rl.on(1, rl.XP)
    g0 = float(np.abs(generator(0, 1, phase=1.0) - (XX + YY) / 2).max())
    g1 = float(np.abs(generator(0, 1, phase=1j) - (XY - YX) / 2).max())
    say(f"exchange generator against (XX+YY)/2 at phase 1 and against (XY-YX)/2 at phase i: {g0}, {g1}")
    assert g0 == 0.0 and g1 == 0.0, "exchange generator convention"

    basepar = {nm: R.model_parities(by[nm], calib, A) for nm in names}
    covs = {nm: R.cov_from_parities(basepar[nm], nsh[nm]) for nm in names}
    W, rank, _ = R.weight_matrix(idx, covs, inflate, None)
    ev, U = np.linalg.eigh(np.linalg.pinv(W))
    keep = ev > ev.max() * 1e-10
    Wa = (U[:, keep] / np.sqrt(ev[keep])).T
    chi = lambda v: float((meas - v) @ W @ (meas - v))
    say(f"chi2 with no letter {chi(base):.1f} over rank {rank} (99 % point {CHI2.ppf(0.99, rank):.1f}); "
        f"drift inflation {inflate:.3f}; no bootstrap term")

    ro = {q: float(sum(A[k, 0] for k in range(8) if k >> sh & 1)) for q, sh in ((4, 2), (6, 1), (5, 0))}
    say("\n1. residual excitation at the start, A_RO = A_true . M(p); each p bounded by the RO arm's P(1|0)")
    say("  RO arm P(1|0): " + "  ".join(f"q{q} {ro[q]:.4f}" for q in (4, 5, 6)))
    for p4 in (0.0083, 0.02, ro[4]):
        pp = {4: p4, 5: 0.0, 6: 0.0}
        v = model(p=pp)
        d = v - base
        say(f"  p4 = {p4:.4f}: chi2 {chi(v):.1f}; its own term d W d {float(d @ W @ d):.3f}; largest shift of a "
            f"reading {float(np.abs(d).max()):.4f}; smallest entry of A_true {float((A @ np.linalg.inv(flips(pp))).min()):+.5f}")
    hi = [ro[4], ro[5], ro[6]]
    sol = least_squares(lambda x: Wa @ (meas - model(p={4: x[0], 5: x[1], 6: x[2]})), [h / 2 for h in hi],
                        bounds=([0, 0, 0], hi), x_scale=0.01)
    say(f"  fit of p4, p5, p6 within those bounds: {np.round(sol.x, 4).tolist()}, chi2 {float(sol.fun @ sol.fun):.1f} "
        f"(status {sol.status})")

    say("\n2. quasi-static detuning, +-sigma on each of q4, q5, q6 (eight sign choices averaged)")
    for label, mean in (("about zero", {("det", q): 0.0 for q in (4, 5, 6)}),
                        ("about the three-detuning optimum", THREE_OPT)):
        ref = model(dict(mean))
        say(f"  {label}: chi2 at sigma = 0 {chi(ref):.1f}")
        for s in (0.03, 0.1, 0.2, 0.3):
            v = sum(model({k: mean[k] + sg * s for k, sg in zip(sorted(mean), signs)})
                    for signs in itertools.product((1, -1), repeat=3)) / 8
            say(f"    sigma {s} rad/us: largest shift of an odd-Y reading {float(np.abs(v - ref).max()):.4f}, chi2 {chi(v):.1f}")

    say("\n3. residual exchange at every cz, eps and phi0 per pair; the phase (omega_a - omega_b) t is left out")
    for frame in ("nco", "replay"):
        FRAME[0] = frame
        for sub in (EX, EX + FIVE):
            fun = lambda x: Wa @ (meas - model(dict(zip(sub, x))))
            starts = [[0.0] * len(sub)] + [[s * 0.05 if i == j else 0.0 for i in range(len(sub))]
                                           for j in range(2) for s in (1, -1)]
            if len(sub) > 4:
                starts.append([0.0] * 4 + FIVE_OPT)
            best = None
            for st in starts:
                r = least_squares(fun, st, x_scale=[SCALE[k[0]] for k in sub], max_nfev=200)
                c = float(r.fun @ r.fun)
                if best is None or c < best[1]:
                    best = (r.x, c, r.status)
            x, c, status = best
            ex = "  ".join(f"pair {p}: eps {math.hypot(x[sub.index(('exc', p))], x[sub.index(('exs', p))]):.3f}"
                           f" phi0 {math.atan2(x[sub.index(('exs', p))], x[sub.index(('exc', p))]):+.2f}" for p in PAIRS)
            rest = "  ".join(f"{k[0]}{k[1]} {v:+.3f}" for k, v in zip(sub[4:], x[4:]))
            say(f"  [{frame} frame] {len(sub)} letters: chi2 {c:.1f} (99 % point {CHI2.ppf(0.99, rank - len(sub)):.1f}; "
                f"status {status})  {ex}  {rest}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
