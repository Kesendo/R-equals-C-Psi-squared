"""The reality-window flight: XY+YX on N = 3, read through the odd-Y sector.

Pre-registration: experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md.

    python run_reality_window.py --gate                  # offline: the zero, replication, control, a full simulated
                                                         # flight through analyze(), letter power
    python run_reality_window.py --hardware              # day-of conditions, submit, persist, analyze
    python run_reality_window.py --recover STUB.json     # rebuild the counts file from a submitted job id
    python run_reality_window.py --analyze RESULT.json   # the committed verdict rules on measured counts

The 98 circuits are frozen in reality_window_circuits.qpy by reality_window_freeze.py, each with an
explicit as-late-as-possible schedule on the path; the runner never re-transpiles them. Each circuit is
submitted twice at 4096 shots in one job, in an order shuffled with a fixed seed.
Credentials: IBM_QUANTUM_TOKEN from the environment or a .env here (gitignored), else the saved account.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import random
import sys
import time
import warnings
from pathlib import Path

import numpy as np

import reality_window_letters as rl

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
QPY = HERE / "reality_window_circuits.qpy"
FREEZE = HERE / "reality_window_freeze.json"
APRIL_QPY = REPO / "data" / "ibm_f83_signature_april2026" / "flown_circuits_d7pol1e7g7gs73cf7j90.qpy"
RESULTS = HERE / "results"

BACKEND = "ibm_marrakesh"
PATH = [4, 5, 6]
SHOTS_HALF = 4096
SHUFFLE_SEED = 20261028
ARMS = ("A", "R", "D0", "D1", "D2", "D3", "D4", "D6", "Kp", "Km")
FIT_ARMS = ("A", "R", "D0", "D1", "D2", "D3", "D4", "D6")
QPY_SHA256 = "e37a378f97070f6077b00a520993c3364f33847df77d9955645ff0b0479d3f62"

# Section 5 of the pre-registration
SIGMA = 1 / math.sqrt(2 * SHOTS_HALF)          # one-basis two-body reading, before drift inflation
SIGMA_C = SIGMA / math.sqrt(2)
C_LAW = (0.75, 2.0, 3.0)                       # band = [0.75 C_pred - 2 sigma_C, C_pred + 3 sigma_C]
N_BOOT = 200                                   # resamples of the readout counts for the baseline's own noise
CARRY_DCHI2 = 25.0
SEPARATE_DCHI2 = 9.0
P_CRIT = 0.01

# Section 6
MIN_T1, MIN_T2 = 100e-6, 80e-6
MAX_RO, MAX_SX, MAX_CZ = 0.03, 0.001, 0.008
MAX_P10 = 0.01      # P(1|0): the RO arm cannot tell a residual excitation at the start from a readout flip
MAX_CALIB_AGE_H = 24
MIN_REMAINING_S = 600
MAX_EXEC_S = 600

FROZEN = json.loads(FREEZE.read_text(encoding="utf-8"))
DT_US = FROZEN["durations_s"]["dt"] * 1e6
DUR_US = {"sx": {int(q): v * 1e6 for q, v in FROZEN["durations_s"]["sx"].items()},
          "x": {int(q): v * 1e6 for q, v in FROZEN["durations_s"]["x"].items()},
          "cz": {tuple(int(x) for x in k.split("-")): v * 1e6 for k, v in FROZEN["durations_s"]["cz"].items()}}


def calib_from(snapshot):
    return {k: {int(q): float(v) for q, v in snapshot[k].items()} for k in ("t1_us", "t2_us", "p01", "p10")}


CALIB_FROZEN = calib_from(FROZEN["calibration"])
EPS = {"det": 0.05, "sxa": 0.005, "sxp": 0.01, "czp": 0.01, "czk": 0.01, "zz": 0.05, "stk": 0.01}
# Section 5: the largest size each letter can plausibly have on a calibrated Heron (order of magnitude)
PLAUSIBLE = {"det": 0.1, "sxa": 0.01, "sxp": 0.02, "czp": 0.05, "czk": 0.05, "zz": 0.1, "stk": 0.05}


# ── circuits ──────────────────────────────────────────────────────────────

def load_circuits(path=QPY):
    from qiskit import qpy
    with open(path, "rb") as f:
        return qpy.load(f)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def arm_of(name):
    return name.split("_")[0]


def phys_of(arm):
    """phys[l] = physical qubit carrying logical l (0 = measured q0, 1 = middle, 2 = measured q2)."""
    return PATH[::-1] if arm == "R" else PATH


def extract(tc):
    """(name, logical qubits, params, duration_us) on the path, in program order; checks the clbit map."""
    phys = phys_of(arm_of(tc.name))
    m = {p: i for i, p in enumerate(phys)}
    ops, meas = [], {}
    for inst in tc.data:
        name = inst.operation.name
        qs = [tc.find_bit(q).index for q in inst.qubits]
        if name == "barrier":
            continue
        if name == "measure":
            meas[tc.find_bit(inst.clbits[0]).index] = m[qs[0]]
            continue
        if not all(q in m for q in qs):
            raise ValueError(f"{tc.name}: {name} on {qs} leaves the path")
        prm = [float(p) for p in inst.operation.params]
        if name == "delay":
            if inst.operation.unit != "dt":
                raise ValueError(f"{tc.name}: delay in {inst.operation.unit}")
            d = prm[0] * DT_US
        elif name == "rz":
            d = 0.0
        elif name in ("sx", "x"):
            d = DUR_US[name][qs[0]]
        elif name == "cz":
            d = DUR_US["cz"][tuple(sorted(qs))]
        else:
            raise ValueError(name)
        ops.append((name, [m[q] for q in qs], prm, d))
    if meas != {0: 0, 1: 2, 2: 1}:
        raise ValueError(f"{tc.name}: measure map {meas}")
    return ops


def logical(n, b0, b2, b1="Y"):
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import SparsePauliOp
    ps = []
    for (i, j) in ((0, 1), (1, 2)):
        for (a, b) in (("X", "Y"), ("Y", "X")):
            ch = ["I"] * 3; ch[2 - i] = a; ch[2 - j] = b
            ps.append("".join(ch))
    qc = QuantumCircuit(3)
    for q in range(3):
        qc.h(q)
        if q % 2:
            qc.z(q)
    gate = PauliEvolutionGate(SparsePauliOp(ps, coeffs=np.ones(len(ps))), time=0.8 / 3)
    for _ in range(n):
        qc.append(gate, range(3))
    for q, b in ((0, b0), (2, b2), (1, b1)):
        if b == "X":
            qc.h(q)
        elif b == "Y":
            qc.sdg(q); qc.h(q)
    return qc.decompose(gates_to_decompose=["PauliEvolution"])


# ── readout matrix, physical order (q4, q6, q5) = measured order of every arm but R ──

def readout_matrix_snapshot(calib):
    return rl.assignment_product([(calib["p01"][q], calib["p10"][q]) for q in (4, 6, 5)])


def readout_matrix_measured(counts_by_name):
    """Column j = outcome distribution of RO_{b0}{b2}{b1}, j = 4 b0 + 2 b2 + b1 (bits on q4, q6, q5)."""
    A = np.zeros((8, 8))
    for j in range(8):
        A[:, j] = rl.counts_to_probs(counts_by_name[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"])
    return A


SWAP02 = np.zeros((8, 8))
for _k in range(8):
    SWAP02[((_k >> 1) & 1) << 2 | ((_k >> 2) & 1) << 1 | (_k & 1), _k] = 1


def arm_readout(A_phys, arm):
    """In arm R measured bit 0 sits on q6 and bit 1 on q4: the physical matrix with bits 0 and 1 exchanged."""
    return SWAP02 @ A_phys @ SWAP02 if arm == "R" else A_phys


# ── the model of one circuit ──────────────────────────────────────────────

def model_parities(tc, calib, A_phys, L=None, ops=None):
    arm = arm_of(tc.name)
    phys = phys_of(arm)
    T1 = [calib["t1_us"][p] for p in phys]
    T2 = [calib["t2_us"][p] for p in phys]
    rho = rl.replay(ops if ops is not None else extract(tc), phys, T1, T2, L)
    return rl.parities_probs(arm_readout(A_phys, arm) @ rl.probs_rho(rho))


_OPS = {}


def model_vector(circs, idx, calib, A_phys, L=None):
    names = sorted({n for n, _ in idx})
    by = {}
    for c in circs:
        if c.name in names:
            if c.name not in _OPS:
                _OPS[c.name] = extract(c)
            by[c.name] = model_parities(c, calib, A_phys, L, _OPS[c.name])
    return np.array([by[n][m] for n, m in idx])


def fingerprints(circs, idx, calib, A_phys):
    F = {}
    for k in rl.LETTERS:
        e = EPS[k[0]]
        F[k] = (model_vector(circs, idx, calib, A_phys, {k: e}) - model_vector(circs, idx, calib, A_phys, {k: -e})) / (2 * e)
    return F


def cov_from_parities(par, n):
    full = dict(par); full[(0, 0, 0)] = 1.0
    C = {}
    for m1 in rl.MASKS:
        for m2 in rl.MASKS:
            x = tuple(a ^ b for a, b in zip(m1, m2))
            C[(m1, m2)] = (full[x] - par[m1] * par[m2]) / n
    return C


def weight_matrix(idx, covs, inflate=1.0, extra=None):
    n = len(idx)
    S = np.zeros((n, n))
    pos = {}
    for i, (nm, m) in enumerate(idx):
        pos.setdefault(nm, []).append((i, m))
    for nm, lst in pos.items():
        for i, m1 in lst:
            for j, m2 in lst:
                S[i, j] = covs[nm][(m1, m2)]
    S = S * inflate
    if extra is not None:
        S = S + extra
    w, V = np.linalg.eigh(S)
    keep = w > 1e-9 * w.max()          # within one circuit some parities can be exact multiples of each other
    return (V[:, keep] / w[keep]) @ V[:, keep].T, int(keep.sum()), (V[:, keep] / np.sqrt(w[keep])).T


def S_from(par, arm):
    return par[f"{arm}_YZ"][(1, 1, 0)]


# ── the letters, fitted through the full replay ───────────────────────────

_WORKER = {}


def _nl_fit(args):
    """One letter, one amplitude: minimise |Wh (r - (model(a) - base))|^2 from four starts (the linear estimate,
    half and twice it, and its negative); the best converged fit, or the best of all flagged as unconverged."""
    from scipy.optimize import least_squares
    k, a0, r, Wh, idx, calib, A_phys, base = args
    if "circs" not in _WORKER:
        _WORKER["circs"] = load_circuits()
    circs = _WORKER["circs"]
    fun = lambda x: Wh @ (r - (model_vector(circs, idx, calib, A_phys, {k: x[0]}) - base))
    scale = max(abs(a0), EPS[k[0]])
    starts = (a0, 0.5 * a0, 2.0 * a0, -a0) if a0 != 0 else (5 * EPS[k[0]], -5 * EPS[k[0]])
    conv, anyfit = None, None
    for start in starts:
        f = least_squares(fun, [start], x_scale=[scale], diff_step=[1e-3], max_nfev=60)
        c = float(f.fun @ f.fun)
        J = f.jac[:, 0]
        cand = (k, float(f.x[0]), c, 1 / math.sqrt(float(J @ J)) if J @ J > 0 else float("inf"))
        if anyfit is None or c < anyfit[2]:
            anyfit = cand
        if f.success and (conv is None or c < conv[2]):
            conv = cand
    return conv + (True,) if conv is not None else anyfit + (False,)


def letter_fits(r, W, Wh, F, idx, calib, A_phys, base, pool, only=None):
    chi0, lin = rl.one_letter_fits(r, W, F)
    keys = only if only is not None else rl.LETTERS
    jobs = [(k, lin[k][0], r, Wh, idx, calib, A_phys, base) for k in keys]
    fits = {k: (lin[k][0], lin[k][1], lin[k][2], "linear") for k in rl.LETTERS}
    for k, a, c, s, ok in (pool.map(_nl_fit, jobs) if pool is not None else map(_nl_fit, jobs)):
        fits[k] = (a, s, c, "replay" if ok else "replay, unconverged")
    return float(r @ W @ r), fits


def verdict_letters(chi0, fits, rank):
    from scipy.stats import chi2
    crit0, crit1 = chi2.ppf(1 - P_CRIT, rank), chi2.ppf(1 - P_CRIT, rank - 1)
    nl = {k: v for k, v in fits.items() if v[3].startswith("replay")}
    best = min(nl, key=lambda k: nl[k][2])
    cb = nl[best][2]
    near = [k for k in nl if k != best and nl[k][2] - cb < SEPARATE_DCHI2]
    carries = chi0 - cb >= CARRY_DCHI2
    if chi0 <= crit0 and not carries:
        v = "NO SIGNAL"
    elif carries and cb <= crit1:
        v = "LETTER FOUND"
    elif carries:
        v = "SIGNAL, NO SINGLE LETTER"
    else:
        v = "SIGNAL OUTSIDE THE LIBRARY"
    size = "plausible" if abs(nl[best][0]) <= PLAUSIBLE[best[0]] else "implausible"
    return {"verdict": v, "chi0": chi0, "crit0": crit0, "crit1": crit1, "best": best, "chi_best": cb,
            "best_amplitude": nl[best][0], "best_size": size, "not_separated": near,
            "unconverged": [k for k, v in nl.items() if v[3] != "replay"]}


# ── the other verdicts ────────────────────────────────────────────────────

def parse_gate(par):
    xz, zx = par["A_XZ"][(1, 1, 0)], par["A_ZX"][(1, 1, 0)]
    return (xz < -0.5) and (-0.75 < zx < -0.3) and (xz < zx - 0.15), {"A <X0Z2>": xz, "A <Z0X2>": zx}


def control_band(c_pred, sig_c):
    return (C_LAW[0] * c_pred - C_LAW[1] * sig_c, c_pred + C_LAW[2] * sig_c)


def verdict_control(par, c_pred, sig_c):
    C = (S_from(par, "Kp") - S_from(par, "Km")) / 2
    lo, hi = control_band(c_pred, sig_c)
    if C <= -3 * sig_c:
        return "SIGN CHAIN BROKEN", C, (lo, hi)
    if C < 3 * sig_c:
        return "CONTROL DEAD", C, (lo, hi)
    return ("PASS" if lo <= C <= hi else "MAGNITUDE OFF"), C, (lo, hi)


def verdict_replication(s, sig):
    if s >= 3 * sig:
        return "REPLICATED"
    if s <= -3 * sig:
        return "REVERSED"
    return "ABSENT" if abs(s) < 2 * sig else "INCONCLUSIVE"


# ── analysis ──────────────────────────────────────────────────────────────

def drift_inflation(halves):
    """Mean of (half1 - half2)^2 / (var1 + var2) over the odd-Y parities of the fit arms with nonzero variance."""
    terms = []
    for nm, (h1, h2) in halves.items():
        if arm_of(nm) not in FIT_ARMS:
            continue
        b = nm.split("_")[1]
        p1, c1 = rl.parities_counts(h1); p2, c2 = rl.parities_counts(h2)
        for m in rl.MASKS:
            if not rl.odd_y((b[0], b[1]), m):
                continue
            v = c1[(m, m)] + c2[(m, m)]
            if v > 0:
                terms.append((p1[m] - p2[m]) ** 2 / v)
    return max(1.0, float(np.mean(terms))), float(np.mean(terms)), len(terms)


def analyze(result_path, pool=None, quiet=False):
    data = json.loads(Path(result_path).read_text(encoding="utf-8"))
    if data.get("qpy_sha256") != sha256(QPY):
        raise SystemExit(f"result was taken with circuits {data.get('qpy_sha256')}, this file is {sha256(QPY)}")
    circs = {c.name: c for c in load_circuits()}
    calib = calib_from(data["snapshot"])
    halves = data["counts"]
    tot = {}
    for nm, (h1, h2) in halves.items():
        t = dict(h1)
        for k, v in h2.items():
            t[k] = t.get(k, 0) + v
        tot[nm] = t
    inflate, mean_t, n_t = drift_inflation(halves)
    sig = SIGMA * math.sqrt(inflate)
    sig_c = SIGMA_C * math.sqrt(inflate)
    A_meas = readout_matrix_measured(tot)
    par, n_shots = {}, {}
    for nm, t in tot.items():
        if not nm.startswith("RO_"):
            par[nm] = rl.parities_counts(t)[0]
            n_shots[nm] = sum(t.values())
    basepar = {nm: model_parities(circs[nm], calib, A_meas) for nm in par}
    # the multinomial covariance of each circuit's parities, from the baseline's predicted parities: a bit that
    # comes out nearly unanimous in the counts would otherwise carry a variance estimated from a handful of flips
    covs = {nm: cov_from_parities(basepar[nm], n_shots[nm]) for nm in par}
    c_pred = (S_from(basepar, "Kp") - S_from(basepar, "Km")) / 2
    out = {"drift": {"mean_term": mean_t, "n_terms": n_t, "inflation": inflate},
           "readout_matrix_measured": A_meas.tolist(),
           "readout_matrix_snapshot": readout_matrix_snapshot(calib).tolist()}
    ok, pv = parse_gate(par)
    v0 = verdict_control(par, c_pred, sig_c)
    void = (not ok) or v0[0] in ("SIGN CHAIN BROKEN", "CONTROL DEAD")
    s_arm = {a: S_from(par, a) - S_from(basepar, a) for a in ARMS}
    v1 = verdict_replication(s_arm["A"], sig)
    out.update({"P": ("PASS" if ok else "FAIL", pv), "V0": v0[0], "C": v0[1], "C_pred": c_pred, "C_band": v0[2],
                "half_sum_K": (S_from(par, "Kp") + S_from(par, "Km")) / 2 - (S_from(basepar, "Kp") + S_from(basepar, "Km")) / 2,
                "V1": "VOID" if void else v1, "s": s_arm, "sigma_eff": sig})
    fit_names = sorted(n for n in par if arm_of(n) in FIT_ARMS)
    idx = rl.reading_index(fit_names)
    base = model_vector(list(circs.values()), idx, calib, A_meas)
    meas = np.array([par[n][m] for n, m in idx])
    # the baseline carries the readout matrix's own shot noise: resample the RO counts, rebuild the baseline
    rng = np.random.default_rng(SHUFFLE_SEED)
    boots = []
    for _ in range(N_BOOT):
        A_b = np.zeros((8, 8))
        for j in range(8):
            col = rl.counts_to_probs(tot[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"])
            n_ro = sum(tot[f"RO_{j >> 2 & 1}{j >> 1 & 1}{j & 1}"].values())
            A_b[:, j] = rng.multinomial(n_ro, col) / n_ro
        boots.append(model_vector(list(circs.values()), idx, calib, A_b))
    extra = np.cov(np.array(boots).T)
    W, rank, Wh = weight_matrix(idx, covs, inflate, extra)
    F = fingerprints(list(circs.values()), idx, calib, A_meas)
    own = pool is None
    if own:
        from concurrent.futures import ProcessPoolExecutor
        pool = ProcessPoolExecutor(max_workers=17)
    chi0, fits = letter_fits(meas - base, W, Wh, F, idx, calib, A_meas, base, pool)
    if own:
        pool.shutdown()
    v2 = verdict_letters(chi0, fits, rank)
    out.update({"V2": "VOID" if void else v2["verdict"], "V2_detail": v2, "rank": rank, "n_readings": len(idx),
                "fits": {str(k): v for k, v in fits.items()}})
    if not quiet:
        print(f"drift: mean (h1-h2)^2/var = {mean_t:.3f} over {n_t} parities -> inflation {inflate:.3f}; sigma_eff {sig:.4f}")
        print("P  parse gate:", out["P"][0], {k: round(v, 3) for k, v in pv.items()})
        print(f"V0 control: {v0[0]}  C = {v0[1]:+.4f}, predicted {c_pred:+.4f}, band ({v0[2][0]:.3f}, {v0[2][1]:.3f})")
        print(f"V1 replication: {out['V1']}  s_A = {s_arm['A']:+.4f}  (half-sum of K: {out['half_sum_K']:+.4f})")
        print("   s per arm:", {a: round(v, 4) for a, v in s_arm.items()})
        print(f"V2 letters: {out['V2']}  chi2_null {chi0:.1f} vs {v2['crit0']:.1f}; best {v2['best']} chi2 "
              f"{v2['chi_best']:.1f} vs {v2['crit1']:.1f}, amplitude {v2['best_amplitude']:+.4f} "
              f"({v2['best_size']} size, bound {PLAUSIBLE[v2['best'][0]]}); not separated: "
              f"{[str(k) for k in v2['not_separated']]}; unconverged fits: {[str(k) for k in v2['unconverged']]}")
        for k, (a, s, c, how) in sorted(fits.items(), key=lambda kv: kv[1][2])[:8]:
            print(f"   {str(k):24s} a = {a:+.4f} ± {s:.4f}  chi2 {c:.1f}  ({how})")
    return out


# ── the gate (offline) ────────────────────────────────────────────────────

def gate():
    from qiskit.quantum_info import Statevector
    from qiskit import QuantumCircuit
    circs = load_circuits()
    h = sha256(QPY)
    print(f"frozen file sha256 {h}; {len(circs)} circuits")
    assert h == QPY_SHA256 and h == FROZEN["sha256"] and len(circs) == 98
    by = {c.name: c for c in circs}
    # (1) the zero: compiled == logical Trotter (delays are identities)
    worst = 0.0
    for c in circs:
        arm, b = c.name.split("_")
        if arm.startswith("K") or arm == "RO":
            continue
        n = 3 if arm in ("A", "R") else int(arm[1:])
        qc = QuantumCircuit(3)
        for name, q, p, _ in extract(c):
            if name != "delay":
                getattr(qc, name)(*p, *q) if p else getattr(qc, name)(*q)
        worst = max(worst, float(np.max(np.abs(Statevector(qc).probabilities([0, 2, 1])
                                                - Statevector(logical(n, b[0], b[1])).probabilities([0, 2, 1])))))
    print(f"(1) compiled vs logical Trotter, max |dp| over 72 circuits: {worst:.1e}")
    assert worst < 1e-12
    second_route(by)
    # (2) arm A = April flown + a q1-only tail: the same gate sequence on every qubit (scheduling may
    # interleave independent qubits differently, which changes neither the unitary nor the explicit timing)
    april = {c.name: c for c in load_circuits(APRIL_QPY)}

    def per_qubit(tc):
        d = {}
        for i in tc.data:
            if i.operation.name in ("measure", "delay", "barrier"):
                continue
            qs = tuple(tc.find_bit(q).index for q in i.qubits)
            for q in qs:
                d.setdefault(q, []).append((i.operation.name, qs, tuple(float(p) for p in i.operation.params)))
        return d
    same = 0
    for b0 in "XYZ":
        for b2 in "XYZ":
            m, t = per_qubit(by[f"A_{b0}{b2}"]), per_qubit(april[f"pi2_odd_pure_{b0}{b2}"])
            same += (m[4] == t[4] and m[6] == t[6] and m[5][:len(t[5])] == t[5]
                     and all(x[1] == (5,) for x in m[5][len(t[5]):]))
    print(f"(2) arm A = April flown circuit, per qubit, + a q1-only tail: {same}/9")
    assert same == 9
    # (3) predictions, frozen calibration, ideal gates
    A0 = np.eye(8)
    A1 = readout_matrix_snapshot(CALIB_FROZEN)
    par0 = {c.name: model_parities(c, CALIB_FROZEN, A0) for c in circs if arm_of(c.name) != "RO"}
    par1 = {c.name: model_parities(c, CALIB_FROZEN, A1) for c in circs if arm_of(c.name) != "RO"}
    floor = 0.0
    for arm in ARMS:
        odd = lambda P: max(abs(P[f"{arm}_{b0}{b2}"][m]) for b0 in "XYZ" for b2 in "XYZ" for m in rl.MASKS
                            if rl.odd_y((b0, b2), m))
        print(f"(3) {arm:3s} no readout: S {S_from(par0, arm):+.4f}, largest odd-Y {odd(par0):.4f} | "
              f"snapshot readout: S {S_from(par1, arm):+.4f}, largest odd-Y {odd(par1):.4f}")
        if arm in FIT_ARMS:
            floor = max(floor, odd(par0))
            # the floor is T1 acting inside the compiled frames, where the cooling jump is not the logical one;
            # it has to stay well under the shot noise of one reading, or the zero of section 1 is not a zero
            assert odd(par0) < SIGMA / 4
    print(f"    floor of the zero over the fit arms: {floor:.4f} = {floor / SIGMA:.2f} sigma")
    c_pred = (S_from(par1, "Kp") - S_from(par1, "Km")) / 2
    v0 = verdict_control(par1, c_pred, SIGMA_C)
    ok, pv = parse_gate(par1)
    print(f"(4) control C {c_pred:+.4f}, band ({v0[2][0]:.3f}, {v0[2][1]:.3f}); parse gate {ok} "
          f"{ {k: round(v, 3) for k, v in pv.items()} }")
    assert v0[0] == "PASS" and ok
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=17) as pool:
        simulated_flight(circs, pool)
        power(circs, par1, pool)


def _independent(ops, phys, L):
    """A second route for the replay without decay: pure state, piecewise-constant Hamiltonian, the schedule's
    segments; sx/x as finite pulses with detuning and static ZZ on during them, cz at its midpoint, rz at its
    start. Written apart from reality_window_letters.replay; shares only the elementary matrices."""
    from scipy.linalg import expm
    det = [L.get(("det", phys[l]), 0.0) for l in range(3)]
    lg = {p: l for l, p in enumerate(phys)}
    zz = {tuple(sorted((lg[a], lg[b]))): L.get(("zz", (a, b)), 0.0) for (a, b) in rl.PAIRS}
    Hb = np.zeros((8, 8), dtype=complex)
    for q in range(3):
        Hb += det[q] / 2 * rl.ZOP[q]
    for (a, b), z in zz.items():
        Hb += z / 2 * rl.ZOP[a] @ rl.ZOP[b]
    clock, segs = [0.0, 0.0, 0.0], []
    for i, (name, qs, prm, d) in enumerate(ops):
        st = max(clock[q] for q in qs)
        for q in qs:
            clock[q] = st + d
        if name != "delay":
            segs.append((round(st, 9), round(st + d, 9), name, qs, prm, i))
    inst = [(s_[0] if s_[2] == "rz" else round((s_[0] + s_[1]) / 2, 9), s_[5], s_) for s_ in segs if s_[2] in ("rz", "cz")]
    times = sorted({0.0, round(max(clock), 9)} | {x for s_ in segs for x in s_[:2]} | {t for t, _, _ in inst})
    psi = np.zeros(8, dtype=complex); psi[0] = 1
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    for k, t in enumerate(times):
        for _, _, s_ in sorted([x for x in inst if x[0] == t], key=lambda x: x[1]):
            psi = (rl.on(s_[3][0], rl.rz(s_[4][0])) if s_[2] == "rz" else rl.czm(*sorted(s_[3]))) @ psi
        if k + 1 < len(times):
            t1, H = times[k + 1], Hb.copy()
            for s_ in segs:
                if s_[2] in ("sx", "x") and s_[0] <= t + 1e-12 and s_[1] >= t1 - 1e-12:
                    H = H + (math.pi / 2 if s_[2] == "sx" else math.pi) / (s_[1] - s_[0]) / 2 * rl.on(s_[3][0], X)
            psi = expm(-1j * H * (t1 - t)) @ psi
    return np.outer(psi, psi.conj())


def second_route(by):
    """(1b) the replay against the independent evaluation, without decay, at sizes up to those the April
    reading needs from a single letter."""
    inf = [1e12] * 3
    worst = 0.0
    for L in ({}, {("det", 5): 0.5}, {("det", 4): -0.46}, {("zz", (4, 5)): -1.3}, {("zz", (5, 6)): -1.08},
              {("det", 6): -0.72, ("zz", (5, 6)): -0.5}):
        for nm in ("A_YZ", "R_YZ", "D3_YZ", "D6_ZY", "Kp_YZ", "D1_XY"):
            ops = extract(by[nm]); phys = phys_of(arm_of(nm))
            a = rl.parities_rho(rl.replay(ops, phys, inf, inf, L))
            b = rl.parities_rho(_independent(ops, phys, L))
            worst = max(worst, max(abs(a[m] - b[m]) for m in rl.MASKS))
    print(f"(1b) replay against an independent finite-pulse evaluation, 36 cases up to 1.3 rad/us: max |diff| {worst:.1e}")
    assert worst < 1e-9


def three_qubit(tc):
    from qiskit import QuantumCircuit, ClassicalRegister
    m = {p: i for i, p in enumerate(phys_of(arm_of(tc.name)))}
    qc = QuantumCircuit(3); qc.add_register(ClassicalRegister(3, "meas"))
    for inst in tc.data:
        if inst.operation.name == "barrier":
            continue
        qc.append(inst.operation, [m[tc.find_bit(q).index] for q in inst.qubits],
                  [tc.find_bit(c).index for c in inst.clbits])
    return qc


def simulated_flight(circs, pool):
    """(5) Aer runs the 98 frozen circuits (ideal gates, the snapshot's readout errors on the right physical
    qubits) in the hardware's two-half shuffled format; the file goes through analyze() unchanged."""
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, ReadoutError
    sims = {}
    for arm in ("plain", "R"):
        nm = NoiseModel()
        for l, p in enumerate(phys_of("R" if arm == "R" else "A")):
            p01, p10 = CALIB_FROZEN["p01"][p], CALIB_FROZEN["p10"][p]
            nm.add_readout_error(ReadoutError([[1 - p01, p01], [p10, 1 - p10]]), [l])
        sims[arm] = AerSimulator(method="density_matrix", noise_model=nm)
    by = {c.name: c for c in circs}
    order = [(nm, h) for nm in sorted(by) for h in (0, 1)]
    random.Random(SHUFFLE_SEED).shuffle(order)
    counts = {c.name: [None, None] for c in circs}
    for pos, (nm, h) in enumerate(order):
        sim = sims["R" if arm_of(nm) == "R" else "plain"]
        counts[nm][h] = sim.run(three_qubit(by[nm]), shots=SHOTS_HALF, seed_simulator=pos + 1).result().get_counts()
    RESULTS.mkdir(exist_ok=True)
    fake = RESULTS / "reality_window_simulated_flight.json"
    snap = dict(FROZEN["calibration"])
    fake.write_text(json.dumps({"job_id": "simulated", "qpy_sha256": sha256(QPY), "order": order,
                                "snapshot": snap, "counts": counts}, default=str), encoding="utf-8")
    t0 = time.time()
    out = analyze(fake, pool=pool)
    print(f"(5) simulated flight through analyze() in {time.time() - t0:.0f} s: P {out['P'][0]}, V0 {out['V0']}, "
          f"V1 {out['V1']}, V2 {out['V2']} (rank {out['rank']} of {out['n_readings']})")
    worst = max(abs(v) for a, v in out["s"].items() if a in FIT_ARMS) / out["sigma_eff"]
    print(f"    largest |s| over the fit arms: {worst:.2f} sigma_eff (one draw; ABSENT needs < 2, so INCONCLUSIVE "
          f"is expected in about one draw of twenty)")
    assert out["P"][0] == "PASS" and out["V0"] == "PASS" and out["V1"] in ("ABSENT", "INCONCLUSIVE")
    assert worst < 4 and out["V2"] == "NO SIGNAL"


def _scenario(args):
    """One letter at the size that makes s_A = target through the full replay, one noise draw (counts noise on
    the readings and on the readout-calibration circuits), then the verdict rule with nonlinear refits."""
    k, target, seed = args          # k: a library letter, or ("outside", name, {letter: weight}) for a mixture
    rng = np.random.default_rng(seed)
    if "circs" not in _WORKER:
        _WORKER["circs"] = load_circuits()
    circs = _WORKER["circs"]
    fit = [c for c in circs if arm_of(c.name) in FIT_ARMS]
    idx = rl.reading_index([c.name for c in fit])
    A_true = readout_matrix_snapshot(CALIB_FROZEN)
    base_true = model_vector(circs, idx, CALIB_FROZEN, A_true)
    iA = idx.index(("A_YZ", (1, 1, 0)))
    F = fingerprints(circs, idx, CALIB_FROZEN, A_true)
    if k[0] == "outside":
        mix = k[2]
        Lof = lambda a: {q: a * w for q, w in mix.items()}
        slope = (model_vector(circs, idx, CALIB_FROZEN, A_true, Lof(1e-3)) - base_true)[iA] / 1e-3
    else:
        Lof = lambda a: {k: a}
        slope = F[k][iA]
    g = lambda a: (model_vector(circs, idx, CALIB_FROZEN, A_true, Lof(a)) - base_true)[iA]
    a = target / slope if abs(slope) > 1e-6 else 0.0
    grid = sorted({f * a for f in (0.25, 0.5, 1, 1.5, 2, 3, 4)} | {-f * a for f in (0.25, 0.5, 1, 1.5, 2, 3, 4)},
                  key=abs) if a != 0 else []
    reached, hit = 0.0, None
    vals = {x: g(x) for x in grid}
    for sign in (1, -1):                               # the smallest size on each side that reaches the target
        side = [x for x in grid if x * sign > 0]
        prev = 0.0
        for x in side:
            if vals[x] * target > 0 and abs(vals[x]) >= abs(target):
                lo, hi = prev, x
                for _ in range(25):
                    mid = 0.5 * (lo + hi)
                    if g(mid) * target > 0 and abs(g(mid)) >= abs(target):
                        hi = mid
                    else:
                        lo = mid
                if hit is None or abs(hi) < abs(hit):
                    hit = hi
                break
            prev = x
    if hit is not None:
        a = hit; reached = g(a)
    elif grid:
        a = max(grid, key=lambda x: vals[x] * np.sign(target)); reached = vals[a]
    Ltrue = Lof(a)
    # sample readings: multinomial per circuit from the true model, two halves; readout circuits likewise
    A_est = np.zeros((8, 8))
    for j in range(8):
        A_est[:, j] = rng.multinomial(2 * SHOTS_HALF, A_true[:, j]) / (2 * SHOTS_HALF)
    covs, meas = {}, {}
    for c in fit:
        arm = arm_of(c.name)
        phys = phys_of(arm)
        rho = rl.replay(_OPS.setdefault(c.name, extract(c)), phys, [CALIB_FROZEN["t1_us"][p] for p in phys],
                        [CALIB_FROZEN["t2_us"][p] for p in phys], Ltrue)
        p = np.clip(arm_readout(A_true, arm) @ rl.probs_rho(rho), 0, None); p /= p.sum()
        draw = rng.multinomial(2 * SHOTS_HALF, p) / (2 * SHOTS_HALF)
        meas[c.name] = rl.parities_probs(draw)
        covs[c.name] = cov_from_parities(meas[c.name], 2 * SHOTS_HALF)
    r_meas = np.array([meas[n][m] for n, m in idx])
    base = model_vector(circs, idx, CALIB_FROZEN, A_est)
    W, rank, Wh = weight_matrix(idx, covs)
    F_est = fingerprints(circs, idx, CALIB_FROZEN, A_est)
    chi0, fits = letter_fits(r_meas - base, W, Wh, F_est, idx, CALIB_FROZEN, A_est, base, None)
    v = verdict_letters(chi0, fits, rank)
    chik = fits[k][2] if k[0] != "outside" else float("nan")
    return k, target, a, reached, v["verdict"], v["best"], v["not_separated"], chik, v["crit1"], v["best_size"]


OUTSIDE = [("outside", "sx axis tilt on q5", {("sxt", 5): 1.0}),
           ("outside", "sx axis tilt on q4", {("sxt", 4): 1.0}),
           ("outside", "detuning q5 + static ZZ 4-5, equal parts", {("det", 5): 1.0, ("zz", (4, 5)): 1.0}),
           ("outside", "CZ phase on q5 from both pairs, equal", {("czp", (4, 5), 5): 1.0, ("czp", (5, 6), 5): 1.0}),
           ("outside", "sx amplitude on all three, equal", {("sxa", 4): 1.0, ("sxa", 5): 1.0, ("sxa", 6): 1.0})]


def power(circs, par, pool, draws_seed=11):
    """(6) An upper bound on what V2 can do: data drawn from the same model family the fit uses, one draw per
    scenario, every letter refitted through the full replay; (7) errors outside the library, and what the rule
    calls them."""
    jobs = [(k, t, draws_seed + 17 * i + (0 if t > 0.05 else 1)) for i, k in enumerate(rl.LETTERS) for t in (0.067, 0.033)]
    print("(6) letter power (upper bound; one draw per scenario, full replay, readout matrix re-estimated):")
    for k, t, a, reached, v, best, near, chik, crit1, size in pool.map(_scenario, jobs):
        tag = "itself best" if best == k else f"best {best}"
        print(f"    {str(k):24s} target {t:.3f}: amp {a:+.4f} reaches {reached:+.3f} -> {v}; {tag} ({size}); "
              f"not separated {[str(q) for q in near]}; chi2 itself {chik:.0f} vs {crit1:.0f}")
    jobs = [(k, 0.067, 1000 + i) for i, k in enumerate(OUTSIDE)]
    print("(7) errors outside the library, at s_A = 0.067 where reachable:")
    for k, t, a, reached, v, best, near, chik, crit1, size in pool.map(_scenario, jobs):
        print(f"    {k[1]:44s} weight {a:+.4f} reaches {reached:+.3f} -> {v}; best {best} ({size}); "
              f"not separated {[str(q) for q in near]}")


# ── the flight ────────────────────────────────────────────────────────────

def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    pipe = Path(r"D:/Entwicklung/Projekte/.NET Projekte/AIEvolution/AIEvolution.UI")
    for p in (HERE / ".env", pipe / ".env", pipe / "experiments" / "ibm_quantum_tomography" / ".env"):
        if p.exists():
            from dotenv import load_dotenv
            load_dotenv(p)
    tok = os.environ.get("IBM_QUANTUM_TOKEN")
    return (QiskitRuntimeService(channel="ibm_quantum_platform", token=tok) if tok
            else QiskitRuntimeService(channel="ibm_quantum_platform"))


def day_of(svc, be):
    """Section 6 of the pre-registration. Any failure aborts; there is no override."""
    from datetime import datetime, timezone
    fails = []
    if sha256(QPY) != QPY_SHA256:
        fails.append(f"frozen file hash {sha256(QPY)} != committed {QPY_SHA256}")
    p = be.properties()
    now = datetime.now(timezone.utc)
    dates = []
    snap = {"calibration": str(p.last_update_date), "t1_us": {}, "t2_us": {}, "p01": {}, "p10": {},
            "readout": {}, "sx": {}, "cz": {}}
    for q in PATH:
        props = p.qubit_property(q)
        missing = False
        for key in ("T1", "T2", "readout_error", "prob_meas1_prep0", "prob_meas0_prep1"):
            if key not in props or props[key][0] is None:
                fails.append(f"q{q} has no {key}")
                missing = True
                continue
            dates.append(props[key][1])
        if missing:
            continue
        snap["t1_us"][q] = props["T1"][0] * 1e6
        snap["t2_us"][q] = props["T2"][0] * 1e6
        snap["p01"][q] = props["prob_meas1_prep0"][0]
        snap["p10"][q] = props["prob_meas0_prep1"][0]
        snap["readout"][q] = props["readout_error"][0]
        snap["sx"][q] = p.gate_error("sx", [q])
        if snap["t1_us"][q] < MIN_T1 * 1e6: fails.append(f"T1 q{q} {snap['t1_us'][q]:.0f} us")
        if snap["t2_us"][q] < MIN_T2 * 1e6: fails.append(f"T2 q{q} {snap['t2_us'][q]:.0f} us")
        if snap["readout"][q] > MAX_RO: fails.append(f"readout q{q} {snap['readout'][q]:.4f}")
        if snap["p01"][q] > MAX_P10: fails.append(f"P(1|0) q{q} {snap['p01'][q]:.4f}")
        if snap["sx"][q] > MAX_SX: fails.append(f"sx q{q} {snap['sx'][q]:.5f}")
    for a, b in ((4, 5), (5, 6)):
        e, d = p.gate_property("cz", [a, b], "gate_error")
        dates.append(d)
        snap["cz"][f"{a}-{b}"] = e
        if e > MAX_CZ: fails.append(f"cz {a}-{b} {e:.5f}")
    if dates:
        age_h = max((now - d).total_seconds() / 3600 for d in dates)
        snap["oldest_property_age_h"] = age_h
        if age_h >= MAX_CALIB_AGE_H:
            fails.append(f"oldest calibration property {age_h:.1f} h old")
    rem = svc.usage().get("usage_remaining_seconds")
    snap["remaining_s"] = rem
    if rem is None or rem < MIN_REMAINING_S:
        fails.append(f"allocation remaining {rem}")
    return fails, snap


def write_counts(job, order, snap, stamp):
    res = job.result()
    counts = {}
    for pos, (nm, h) in enumerate(order):
        counts.setdefault(nm, [None, None])[h] = res[pos].data.meas.get_counts()
    try:
        timestamps = job.metrics().get("timestamps")
    except Exception as exc:
        timestamps = repr(exc)
    out = RESULTS / f"reality_window_{BACKEND}_{stamp}.json"
    out.write_text(json.dumps({"job_id": job.job_id(), "backend": BACKEND, "path": PATH, "shots_per_half": SHOTS_HALF,
                               "qpy_sha256": sha256(QPY), "order": order, "snapshot": snap,
                               "job_timestamps": timestamps, "counts": counts},
                              default=str, indent=1), encoding="utf-8")
    print("counts persisted ->", out)
    try:
        print("billed usage:", job.usage())
    except Exception:
        pass
    return out


def hardware():
    from qiskit_ibm_runtime import SamplerV2
    svc = service()
    be = svc.backend(BACKEND)
    fails, snap = day_of(svc, be)
    print("day-of snapshot:", json.dumps(snap, default=str))
    if fails:
        raise SystemExit("ABORT: " + "; ".join(fails))
    circs = {c.name: c for c in load_circuits()}
    order = [(nm, h) for nm in sorted(circs) for h in (0, 1)]
    random.Random(SHUFFLE_SEED).shuffle(order)
    RESULTS.mkdir(exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    stub = RESULTS / f"reality_window_{stamp}.submitted.json"
    record = {"job_id": None, "stamp": stamp, "snapshot": snap, "order": order, "qpy_sha256": sha256(QPY)}
    stub.write_text(json.dumps(record, default=str), encoding="utf-8")
    sampler = SamplerV2(mode=be)
    sampler.options.default_shots = SHOTS_HALF
    sampler.options.max_execution_time = MAX_EXEC_S
    sampler.options.dynamical_decoupling.enable = False
    sampler.options.twirling.enable_gates = False
    sampler.options.twirling.enable_measure = False
    sampler.options.execution.init_qubits = True
    sampler.options.execution.meas_type = "classified"
    t_submit = time.time()
    try:
        job = sampler.run([circs[nm] for nm, _ in order])
    except Exception as exc:
        try:
            recent = [j.job_id() for j in svc.jobs(limit=5, backend_name=BACKEND)]
        except Exception as exc2:
            recent = [f"job listing failed: {exc2!r}"]
        record["submit_error"] = repr(exc); record["recent_jobs"] = recent
        stub.write_text(json.dumps(record, default=str), encoding="utf-8")
        raise SystemExit(f"submission raised {exc!r}; the five newest jobs on {BACKEND} are {recent}, "
                         f"recorded in {stub}: if one of them is this flight, put its id in the stub and --recover")
    record["job_id"] = job.job_id()
    stub.write_text(json.dumps(record, default=str), encoding="utf-8")
    print("submitted", record["job_id"], "->", stub)
    analyze(write_counts(job, order, snap, stamp))


def recover(stub_path):
    stub = json.loads(Path(stub_path).read_text(encoding="utf-8"))
    if not stub.get("job_id"):
        raise SystemExit("the stub has no job id: the submission did not return one; look the job up by its time")
    job = service().job(stub["job_id"])
    analyze(write_counts(job, [tuple(x) for x in stub["order"]], stub["snapshot"], stub["stamp"]))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    if sys.argv[1] == "--gate":
        gate()
    elif sys.argv[1] == "--hardware":
        hardware()
    elif sys.argv[1] == "--recover":
        recover(sys.argv[2])
    elif sys.argv[1] == "--analyze":
        analyze(sys.argv[2])
    else:
        raise SystemExit(__doc__)
