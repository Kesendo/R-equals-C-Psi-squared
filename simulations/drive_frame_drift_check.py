"""The drive-frame and X^3 window verdicts under drift and under a wrong baseline: false-positive rates through the
flight's own analysis (simulations/flight/run_drive_frame.py, analyze_compact), on a chip's frozen circuits.

The simulated flights of the pre-registrations draw from the baseline replay and carry no drift; this script draws
flights the analysis was not built from, all frame-blind, so that every TURNS, 2phi PRESENT, DETUNED or HEARD it
returns is false:
  (a) letters wandering per PUB: every circuit half runs at its own shuffled position in the job, so each PUB draws its
      own detunings on a, b, c, sx/x amplitudes on a, b, c, static ZZ on both bonds and local phases after every cz,
      independent and Gaussian with the sizes below times a factor alpha, carried into the PUB's outcome
      probabilities through the replay's derivatives (each letter's central difference; linear, so the mean of the
      drawn flights is the baseline);
  (b) the readout wandering per PUB: P(1|0) and P(0|1) of every path qubit drawn per PUB, the RO circuits included;
  (c) a wrong baseline: the true state depolarized by 20 % (every parity times 0.8), no drift, which also takes the
      covariance away from the replay's;
  (d) the forest: a common tiny flip rate, and the forest's P(1) wandering per PUB by a relative amount;
  (e) for comparison, the reading-by-reading drift weights the pages first carried, under the letters of (a).
Each scenario: 1000 flights on 16 processes. Printed per scenario: the median kappa of A and B, and the fraction of
flights in which each verdict fires falsely (V1 at 3 sigma, T5, T4, T6 above 99.9 % (and 99 %), E present, W above
its 99.9 % point (OPEN), T5_B, T4_B, T6_B, E_B, VQ DETUNED on any path qubit, VS HEARD and |z| >= 3 on any forest
qubit, P and PB failing).

    python simulations/drive_frame_drift_check.py [CHIP]   ->   simulations/results/drive_frame/<chip>_drift_check.txt
"""
import math
import pickle
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "simulations" / "flight"))
import run_drive_frame as R                          # noqa: E402
import reality_window_letters as rl                  # noqa: E402

N_FLIGHTS = 1000
WORKERS = 16
SIZES = {("det", 4): 0.02, ("det", 5): 0.02, ("det", 6): 0.02,          # rad/us (3.2 kHz)
         ("sxa", 4): 0.002, ("sxa", 5): 0.002, ("sxa", 6): 0.002,       # fraction of the rotation
         ("zz", (4, 5)): 0.01, ("zz", (5, 6)): 0.01,                    # rad/us
         ("czp", (4, 5), 4): 0.005, ("czp", (4, 5), 5): 0.005, ("czp", (5, 6), 5): 0.005, ("czp", (5, 6), 6): 0.005,
         ("sxp", 4): 0.005, ("sxp", 5): 0.005, ("sxp", 6): 0.005,         # phase after each sx/x, rad
         ("czk", (4, 5)): 0.005, ("czk", (5, 6)): 0.005,                  # cz conditional phase, rad
         ("stk", (4, 5), 6): 0.005, ("stk", (5, 6), 4): 0.005}            # the spectator's phase at a cz, rad
FD = {"det": 0.01, "sxa": 0.002, "zz": 0.01, "czp": 0.005, "sxp": 0.005, "czk": 0.005, "stk": 0.005}
P_MARGIN_OCTOBER = 0.15


def fingerprints(chip):
    cache = R.RESULTS / f"drive_frame_{chip.name}_fingerprints_{chip.freeze_sha[:12]}.pkl"
    if cache.exists():
        return pickle.loads(cache.read_bytes())
    F = {}
    for nm in chip.circs:
        if nm.startswith("RO_"):
            continue
        p0 = chip.p_true(nm)
        d = {}
        for k in SIZES:
            e = FD[k[0]]
            d[k] = (chip.p_true(nm, {k: e}) - chip.p_true(nm, {k: -e})) / (2 * e)
        F[nm] = (p0, d)
    R.RESULTS.mkdir(exist_ok=True)
    cache.write_bytes(pickle.dumps(F))
    return F


def flight(chip, F, rng, alpha=0.0, ro_sd=0.0, depol=0.0, forest_rate=None, forest_rel=0.0, letters=None):
    """letters: {letter: sd per PUB} in place of SIZES times alpha (one letter wandering alone)."""
    sizes = letters if letters is not None else {k: alpha * v for k, v in SIZES.items()} if alpha else {}
    comp = {}
    names = sorted(chip.circs)
    for nm in names:
        pc = np.zeros((2, 8))
        fo = np.zeros((2, chip.nf))
        for h in range(2):
            if nm.startswith("RO_"):
                j = 4 * int(nm[3]) + 2 * int(nm[4]) + int(nm[5])
                p = np.zeros(8)
                p[j] = 1.0
            else:
                p0, d = F[nm]
                p = p0.copy()
                for k, sd in sizes.items():
                    p = p + sd * rng.standard_normal() * d[k]
                if depol:
                    p = (1 - depol) * p + depol / 8
                p = np.clip(p, 0, None)
                p /= p.sum()
            ro = chip.ro_snap
            if ro_sd:
                ro = [(min(max(a + ro_sd * rng.standard_normal(), 0.0), 0.5), min(max(b + ro_sd * rng.standard_normal(), 0.0), 0.5))
                      for a, b in ro]
            q = rl.assignment_product(ro) @ p
            pc[h] = rng.multinomial(R.SHOTS_HALF, q / q.sum())
            base = chip.forest_p01 if forest_rate is None else [forest_rate] * chip.nf
            for i, b in enumerate(base):
                b2 = b * (1 + forest_rel * rng.standard_normal()) if forest_rel else b
                fo[h, i] = rng.binomial(R.SHOTS_HALF, min(max(b2, 0.0), 1.0))
        comp[nm] = (pc, fo, np.array([R.SHOTS_HALF, R.SHOTS_HALF], dtype=float))
    return comp


def fires(chip, out):
    xz, zx = out["P"]["<X0Z2>"], out["P"]["<Z0X2>"]
    f = {"P fails": not out["P"]["pass"],
         "P fails at October's margin": not ((xz < -0.5) and (-0.75 < zx < -0.3) and (xz < zx - P_MARGIN_OCTOBER)),
         "V1 REPLICATED": out["V1"]["z"] >= 3,
         "V1 at 3 sigma": abs(out["V1"]["z"]) >= 3,
         "kappa_A": out["drift_A"]["kappa"]}
    for k in ("T5", "T4", "T6"):
        x = out["VF"][k]
        f[f"{k} TURNS"] = x.get("T", 0) > x["points"][0.999] if "T" in x else False
        f[f"{k} > 99 %"] = x.get("T", 0) > x["points"][0.99] if "T" in x else False
    f["E present"] = out["VF"]["E"]["present"] is True
    f["E > 99 %"] = out["VF"]["E"]["T"] > out["VF"]["E"]["points"][0.99]
    zs = [x["z"] for x in out["VS"].values() if math.isfinite(x["z"])]
    f["VS HEARD"] = any(z >= 5 for z in zs)
    f["VS |z| >= 3"] = any(abs(z) >= 3 for z in zs)
    if chip.rider:
        f["PB fails"] = not out["PB"]["pass"]
        f["kappa_B"] = out["drift_B"]["kappa"]
        if "W" in out["WB"]:
            f["W OPEN"] = out["WB"]["W"] > out["WB"]["points"][0.999]
            for k in ("T5", "T4", "T6"):
                x = out["VFB"][k]
                f[f"{k}B TURNS"] = x["T"] > x["points"][0.999]
                f[f"{k}B > 99 %"] = x["T"] > x["points"][0.99]
            f["EB present"] = out["VFB"]["E"]["present"] is True
            f["EB > 99 %"] = out["VFB"]["E"]["T"] > out["VFB"]["E"]["points"][0.99]
            f["W > 95 %"] = out["WB"]["W"] > out["WB"]["points"][0.95]
        f["VQ DETUNED"] = any(out["VQ"][str(l)]["verdict"] == "DETUNED" for l in range(3))
        f["VQ DETUNED per qubit"] = float(np.mean([out["VQ"][str(l)]["verdict"] == "DETUNED" for l in range(3)]))
    return f


class PerReading:
    """The weights the pages first carried, for comparison: kappa_i = mean over the standing variants of the half
    difference's square over its variance, kappa = max(1, mean kappa_i), kappa Sigma = D^1/2 Sigma D^1/2 with
    D = max(kappa, kappa_i), reading by reading."""

    def __init__(self, Sig, deltas, blocks):
        d = np.diag(Sig)
        ki = np.mean([x ** 2 / d for x in deltas], axis=0)
        self.kappa = max(1.0, float(np.mean(ki)))
        D = np.sqrt(np.maximum(self.kappa, ki))
        self.C = Sig * np.outer(D, D)
        self.W, self.rank = R.pinv_rank(self.C)
        self.largest = float(ki.max())

    def cov(self, scale=1.0):
        return scale * self.C

    def winv(self, scale=1.0):
        return self.W / scale


def _worker(args):
    chip_name, seeds, kw = args
    kw = dict(kw)
    if kw.pop("per_reading", False):
        R.DriftCov = PerReading
    chip = R.Chip(chip_name)
    F = fingerprints(chip)
    res = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        res.append(fires(chip, R.analyze_compact(chip, flight(chip, F, rng, **kw))))
    return res


SCENARIOS = (("no drift", {}),
             ("one letter: detuning on a, sd 0.02 rad/us (3.2 kHz)", {"letters": {("det", 4): 0.02}}),
             ("one letter: detuning on b, sd 0.02 rad/us", {"letters": {("det", 5): 0.02}}),
             ("one letter: detuning on a, sd 0.05 rad/us (8 kHz)", {"letters": {("det", 4): 0.05}}),
             ("one letter: detuning on a, sd 0.1 rad/us (16 kHz)", {"letters": {("det", 4): 0.1}}),
             ("one letter: phase on b after cz(b, c), sd 0.03 rad", {"letters": {("czp", (5, 6), 5): 0.03}}),
             ("one letter: phase on b after cz(b, c), sd 0.06 rad", {"letters": {("czp", (5, 6), 5): 0.06}}),
             ("one letter: phase on b after cz(b, c), sd 0.12 rad", {"letters": {("czp", (5, 6), 5): 0.12}}),
             ("readout per PUB, sd 0.02", {"ro_sd": 0.02}),
             ("depolarized 22 % (contrast 0.78), no drift", {"depol": 0.22}),
             ("depolarized 30 % (contrast 0.70), no drift", {"depol": 0.30}),
             ("depolarized 34 % (contrast 0.66), no drift", {"depol": 0.34}),
             ("letters per PUB, alpha 1", {"alpha": 1.0}),
             ("letters per PUB, alpha 2", {"alpha": 2.0}),
             ("letters per PUB, alpha 4", {"alpha": 4.0}),
             ("readout per PUB, sd 0.005", {"ro_sd": 0.005}),
             ("readout per PUB, sd 0.01", {"ro_sd": 0.01}),
             ("depolarized 20 %, no drift", {"depol": 0.2}),
             ("depolarized 20 % and letters alpha 2", {"depol": 0.2, "alpha": 2.0}),
             ("forest rate 2e-5", {"forest_rate": 2e-5}),
             ("forest rate 5e-5", {"forest_rate": 5e-5}),
             ("forest rate 1e-4", {"forest_rate": 1e-4}),
             ("forest rate 0.05, 10 % per PUB", {"forest_rate": 0.05, "forest_rel": 0.1}),
             ("forest rate 0.05, 20 % per PUB", {"forest_rate": 0.05, "forest_rel": 0.2}),
             ("for comparison, the reading-by-reading weights: no drift", {"per_reading": True}),
             ("for comparison, the reading-by-reading weights: letters alpha 1", {"per_reading": True, "alpha": 1.0}),
             ("for comparison, the reading-by-reading weights: letters alpha 2", {"per_reading": True, "alpha": 2.0}))


def main(chip_name):
    from concurrent.futures import ProcessPoolExecutor
    chip = R.Chip(chip_name)
    t0 = time.time()
    fingerprints(chip)
    lines = [f"{chip_name}: path {chip.path}, circuits {chip.freeze['sha256'][:16]}, freeze record {chip.freeze_sha}; {N_FLIGHTS} flights per scenario "
             f"through analyze_compact (fingerprints {time.time() - t0:.0f} s)",
             "letter sizes per PUB at alpha 1: " + ", ".join(f"{k} {v}" for k, v in SIZES.items())]
    print(lines[0], flush=True)
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for label, kw in SCENARIOS:
            seeds = list(range(1, N_FLIGHTS + 1))
            res = [r for part in ex.map(_worker, [(chip_name, seeds[i::WORKERS], kw) for i in range(WORKERS)]) for r in part]
            keys = [k for k in res[0] if not k.startswith("kappa")]
            lines_k = None
            kap = f"kappa_A median {np.median([r['kappa_A'] for r in res]):.2f}" + (
                f", kappa_B {np.median([r['kappa_B'] for r in res]):.2f}" if chip.rider else "")
            lines.append(f"== {label}: {kap}")
            lines.append("   " + "; ".join(f"{k} {np.mean([r.get(k, False) for r in res]):.3f}" for k in keys))
            print(lines[-2] + "\n" + lines[-1], flush=True)
    out = REPO / "simulations" / "results" / "drive_frame" / f"{chip_name}_drift_check.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "ibm_marrakesh")
