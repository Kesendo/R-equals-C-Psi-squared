"""Freeze the drive-frame circuits of one chip (needs the chip's target and properties, read only; no job is submitted).

    python drive_frame_freeze.py ibm_marrakesh | ibm_kingston | ibm_fez

Pre-registrations: experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md (arms A, RO, the forest; sections 3 and 6) and, on
Marrakesh and Kingston, experiments/X3_WINDOW_HARDWARE_PREDICTION.md (arms B and Q; sections 3 and 6).

The path is chosen by the rule of section 3 of the drive-frame page: the preferred path (Marrakesh [4, 5, 6], Kingston
[43, 56, 63]) whenever its three qubits and two bonds pass the guards of section 6 (the properties' ages recorded and
not gated), else
among the linear triples a-b-c (a < c) whose middle b has exactly two neighbours and which pass those guards, the one
with the largest smallest T2. Logical q0, q1, q2 sit on a, b, c. The forest is every coupling-map neighbour of a, b, c
off the path, in ascending physical index.

Circuit names, VARIANT_BASIS:
  A5k_{b0}{b2}  the 2026-04-30 flown XY+YX circuits (pi2_odd_pure), gate for gate with April's physical 4, 5, 6 carried
                to a, b, c, then q1 read in Y; turned by phi = k pi/2 on b (k = 0 is the October arm A)
  A42, A62      the same turned by pi on a and on c
  B5k, B42, B62 the 2026-04-30 flown YZ+ZY circuits (pi2_even_nontruly), q1 read in X, turned the same way
                (Marrakesh and Kingston only)
  Q_{l}{P}{t}   one path qubit l at a time, the other two in |0>: sx, an idle of TAUS_US[t], read in P = X or Y
                (Marrakesh and Kingston only)
  RO_{b0}{b2}{b1} the eight computational states of (q0, q2, q1), prepared with x
The turn is rz(phi) as the first operation on the turned qubit and rz(-phi) as its last before its measurement. Every
circuit measures q0 -> clbit 0, q2 -> clbit 1, q1 -> clbit 2 and forest[i] -> clbit 3 + i, and is scheduled as late as
possible with explicit delays on the path and the forest (the padding on every other qubit is removed), so that the
backend has no timing left to decide and the forest is measured with the path.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

warnings.filterwarnings("ignore")
from qiskit import QuantumCircuit, ClassicalRegister, transpile, qpy
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import ALAPScheduleAnalysis, PadDelay

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
APRIL_QPY = REPO / "data" / "ibm_f83_signature_april2026" / "flown_circuits_d7pol1e7g7gs73cf7j90.qpy"
APRIL_PATH = [4, 5, 6]
CHIPS = ("ibm_marrakesh", "ibm_kingston", "ibm_fez")
PREFERRED = {"ibm_marrakesh": [4, 5, 6], "ibm_kingston": [43, 56, 63]}
RIDERS = ("ibm_marrakesh", "ibm_kingston")          # chips that also carry the X^3 window arms B and Q
SHUFFLE_SEED = {"ibm_fez": 20261030, "ibm_marrakesh": 20261113, "ibm_kingston": 20261113}
TAUS_US = (0.0, 1.0, 2.0, 4.0)
TURNS = (("5", 0), ("5", 1), ("5", 2), ("5", 3), ("4", 2), ("6", 2))   # (role of the turned qubit, k), phi = k pi/2
ROLE_LOGICAL = {"4": 0, "5": 1, "6": 2}             # the role digit names a, b, c by April's Marrakesh qubits

# Section 6 of the drive-frame page
MIN_T1_US, MIN_T2_US = 100.0, 80.0
MAX_RO, MAX_P1G0, MAX_SX, MAX_CZ = 0.03, 0.01, 0.001, 0.008
# The vendor's property ages are recorded and not gated (section 6 of the drive-frame page): the frame verdicts compare
# the turns among themselves, the RO arm measures the day's readout, and the colouring and arm Q read the day's coherence.
MAX_CALIB_AGE_H = float("inf")
FREEZE_MAX_AGE_H = float("inf")


def canonical_sha(obj):
    """sha256 of a record's canonical JSON, independent of line endings and key order."""
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")).hexdigest()


def files(chip):
    return HERE / f"drive_frame_{chip}_circuits.qpy", HERE / f"drive_frame_{chip}_freeze.json"


def service():
    from dotenv import load_dotenv
    from qiskit_ibm_runtime import QiskitRuntimeService
    pipe = Path(r"D:/Entwicklung/Projekte/.NET Projekte/AIEvolution/AIEvolution.UI")
    for p in (HERE / ".env", pipe / ".env", pipe / "experiments" / "ibm_quantum_tomography" / ".env"):
        if p.exists():
            load_dotenv(p)
    tok = os.environ.get("IBM_QUANTUM_TOKEN")
    return (QiskitRuntimeService(channel="ibm_quantum_platform", token=tok) if tok
            else QiskitRuntimeService(channel="ibm_quantum_platform"))


# ── the guards of section 6, on the chip's properties ─────────────────────

def bond_key(target, a, b):
    return (a, b) if (a, b) in target["cz"] else (b, a)


def qubit_record(props, q, now, max_age=MAX_CALIB_AGE_H):
    """(record, failures) for one qubit against the guards; a missing property is a failure."""
    pr = props.qubit_property(q)
    rec, fails, ages = {}, [], []
    for key, name in (("T1", "t1_us"), ("T2", "t2_us"), ("readout_error", "readout"),
                      ("prob_meas1_prep0", "p01"), ("prob_meas0_prep1", "p10")):
        if key not in pr or pr[key][0] is None:
            fails.append(f"q{q} has no {key}")
            continue
        v = float(pr[key][0])
        rec[name] = v * 1e6 if key in ("T1", "T2") else v
        ages.append((now - pr[key][1]).total_seconds() / 3600)
    try:
        # the sx error's date is recorded beside the others (IBM refreshes it rarely on some chips: 209 h on Kingston on
        # 2026-11-13)
        e, d = props.gate_property("sx", [q], "gate_error")
        rec["sx"] = float(e)
        rec["sx_age_h"] = (now - d).total_seconds() / 3600
    except Exception:
        fails.append(f"q{q} has no sx error")
    if fails:
        return rec, fails
    rec["oldest_age_h"] = max(ages)
    if rec["t1_us"] < MIN_T1_US: fails.append(f"T1 q{q} {rec['t1_us']:.0f} us")
    if rec["t2_us"] < MIN_T2_US: fails.append(f"T2 q{q} {rec['t2_us']:.0f} us")
    if rec["readout"] > MAX_RO: fails.append(f"readout q{q} {rec['readout']:.4f}")
    if rec["p01"] > MAX_P1G0: fails.append(f"P(1|0) q{q} {rec['p01']:.4f}")
    if rec["sx"] > MAX_SX: fails.append(f"sx q{q} {rec['sx']:.5f}")
    if rec["oldest_age_h"] >= max_age: fails.append(f"q{q} property {rec['oldest_age_h']:.1f} h old")
    return rec, fails


def bond_record(props, target, a, b, now, max_age=MAX_CALIB_AGE_H):
    k = bond_key(target, a, b)
    try:
        e, d = props.gate_property("cz", list(k), "gate_error")
    except Exception:
        return {}, [f"cz {a}-{b} has no gate_error"]
    age = (now - d).total_seconds() / 3600
    fails = []
    if e > MAX_CZ: fails.append(f"cz {a}-{b} {e:.5f}")
    if age >= max_age: fails.append(f"cz {a}-{b} property {age:.1f} h old")
    return {"cz": float(e), "age_h": age, "key": list(k)}, fails


def path_check(props, target, path, now, max_age=MAX_CALIB_AGE_H):
    recs, fails = {}, []
    for q in path:
        try:
            r, f = qubit_record(props, q, now, max_age)
        except Exception as exc:
            r, f = {}, [f"q{q}: {exc!r}"]
        recs[q] = r
        fails += f
    bonds = {}
    for a, b in zip(path, path[1:]):
        try:
            r, f = bond_record(props, target, a, b, now, max_age)
        except Exception as exc:
            r, f = {}, [f"bond {a}-{b}: {exc!r}"]
        bonds[f"{a}-{b}"] = r
        fails += f
    return recs, bonds, fails


def choose_path(chip, be, props, now):
    """The rule of section 3. Returns (path, record of the choice)."""
    target = be.target
    cm = be.coupling_map
    nb = {q: sorted(set(cm.neighbors(q))) for q in range(be.num_qubits)}
    record = {}
    pref = PREFERRED.get(chip)
    if pref is not None:
        _, _, f = path_check(props, target, pref, now, FREEZE_MAX_AGE_H)
        record["preferred"] = pref
        record["preferred_failures"] = f
        if not f:
            record["rule"] = "preferred path passes"
            return pref, record
    cands, kinds, n_triples = [], {}, 0
    for b in range(be.num_qubits):
        if len(nb[b]) != 2:
            continue
        a, c = nb[b]
        path = [a, b, c]
        n_triples += 1
        recs, _, f = path_check(props, target, path, now, FREEZE_MAX_AGE_H)
        if not f:
            cands.append((min(recs[q]["t2_us"] for q in path), path))
        for x in f:
            k = re.sub(r"[0-9.]+", "#", x)
            kinds[k] = kinds.get(k, 0) + 1
    if not cands:
        # the refusal is recorded, so that a page can say why a chip was not frozen
        out = REPO / "simulations" / "results" / "drive_frame" / f"{chip}_freeze_refused.txt"
        out.parent.mkdir(parents=True, exist_ok=True)
        if pref is not None:
            kinds_pref = record.get("preferred_failures")
        else:
            kinds_pref = None
        t1_ages = sorted((now - props.qubit_property(q)["T1"][1]).total_seconds() / 3600 for q in range(be.num_qubits)
                         if "T1" in props.qubit_property(q) and props.qubit_property(q)["T1"][1] is not None)
        lines = [f"{chip}, {now.isoformat()}: none of {n_triples} linear triples with a middle of degree two passes "
                 f"the guards (property ages recorded, not gated); T1 ages over the chip's qubits: "
                 f"youngest {t1_ages[0]:.1f} h, median {t1_ages[len(t1_ages) // 2]:.1f} h; failure kinds over all "
                 "triples (# a number):"]
        lines += [f"  {v:4d}  {k}" for k, v in sorted(kinds.items(), key=lambda kv: -kv[1])]
        if pref is not None:
            lines.append(f"preferred path {pref}: {kinds_pref}")
        out.write_text(chr(10).join(lines) + chr(10), encoding="utf-8", newline=chr(10))
        raise SystemExit(f"{chip}: no linear triple passes the guards; not frozen ({out.name})")
    cands.sort(key=lambda x: (-x[0], x[1]))
    record["rule"] = "largest smallest T2 among the passing triples with a middle of degree two"
    record["n_candidates"] = len(cands)
    record["top_candidates"] = [{"path": p, "min_t2_us": t} for t, p in cands[:5]]
    return cands[0][1], record


# ── circuits ──────────────────────────────────────────────────────────────

def tail_of(gates):
    q = QuantumCircuit(1)
    for g in gates:
        getattr(q, g)(0)
    return [(i.operation.name, [float(p) for p in i.operation.params])
            for i in transpile(q, basis_gates=["rz", "sx", "x"], optimization_level=1).data]


def new_circuit(nq, nclb):
    out = QuantumCircuit(nq)
    out.add_register(ClassicalRegister(nclb, "meas"))
    return out


def measure_all(out, path, forest):
    out.measure(path[0], 0)
    out.measure(path[2], 1)
    out.measure(path[1], 2)
    for i, f in enumerate(forest):
        out.measure(f, 3 + i)


def from_april(tc, nq, path, forest, tail, turn):
    """April's flown gates with physical 4, 5, 6 carried to path a, b, c; the q1 tail; the turn (lq, phi) or None."""
    pmap = dict(zip(APRIL_PATH, path))
    out = new_circuit(nq, 3 + len(forest))
    if turn is not None:
        out.rz(turn[1], path[turn[0]])
    for inst in tc.data:
        name = inst.operation.name
        if name in ("measure", "barrier", "delay"):
            continue
        qs = [tc.find_bit(q).index for q in inst.qubits]
        if not all(q in pmap for q in qs):
            raise ValueError(f"{tc.name}: {name} on {qs} leaves April's path")
        out.append(inst.operation, [pmap[q] for q in qs])
    for nm, prm in tail:
        getattr(out, nm)(*prm, path[1]) if prm else getattr(out, nm)(path[1])
    if turn is not None:
        out.rz(-turn[1], path[turn[0]])
    measure_all(out, path, forest)
    return out


def q_circuit(nq, path, forest, lq, tail, tau_dt):
    out = new_circuit(nq, 3 + len(forest))
    q = path[lq]
    out.sx(q)
    if tau_dt > 0:
        out.delay(tau_dt, q, unit="dt")
    for nm, prm in tail:
        getattr(out, nm)(*prm, q) if prm else getattr(out, nm)(q)
    measure_all(out, path, forest)
    return out


def ro_circuit(nq, path, forest, b0, b2, b1):
    out = new_circuit(nq, 3 + len(forest))
    for q, b in ((path[0], b0), (path[2], b2), (path[1], b1)):
        if b:
            out.x(q)
    measure_all(out, path, forest)
    return out


def keep_delays(tc, keep):
    """The explicit schedule on the path and the forest; the padding on every other qubit dropped."""
    out = tc.copy_empty_like()
    for inst in tc.data:
        qs = [tc.find_bit(q).index for q in inst.qubits]
        if inst.operation.name == "delay" and qs[0] not in keep:
            continue
        out.append(inst)
    return out


def build(chip, be, path, forest):
    nq = be.num_qubits
    ytail, xtail = tail_of(["sdg", "h"]), tail_of(["h"])
    dt = be.dt
    taus_dt = [int(round(t * 1e-6 / dt)) for t in TAUS_US]
    gran = be.target.granularity
    if any(t % gran for t in taus_dt):
        raise SystemExit(f"idle lengths {taus_dt} dt are not multiples of the granularity {gran}")
    with open(APRIL_QPY, "rb") as f:
        april = {c.name: c for c in qpy.load(f)}
    arms = {}
    rider = chip in RIDERS
    for role, k in TURNS:
        turn = None if k == 0 else (ROLE_LOGICAL[role], k * math.pi / 2)
        for b0 in "XYZ":
            for b2 in "XYZ":
                arms[f"A{role}{k}_{b0}{b2}"] = from_april(april[f"pi2_odd_pure_{b0}{b2}"], nq, path, forest, ytail, turn)
                if rider:
                    arms[f"B{role}{k}_{b0}{b2}"] = from_april(april[f"pi2_even_nontruly_{b0}{b2}"], nq, path, forest,
                                                              xtail, turn)
    if rider:
        for lq in range(3):
            for P, tl in (("X", xtail), ("Y", ytail)):
                for t, tdt in enumerate(taus_dt):
                    arms[f"Q_{lq}{P}{t}"] = q_circuit(nq, path, forest, lq, tl, tdt)
    for b0 in (0, 1):
        for b2 in (0, 1):
            for b1 in (0, 1):
                arms[f"RO_{b0}{b2}{b1}"] = ro_circuit(nq, path, forest, b0, b2, b1)
    pm = PassManager([ALAPScheduleAnalysis(target=be.target), PadDelay(target=be.target)])
    keep = set(path) | set(forest)
    frozen = []
    for nm in sorted(arms):
        c = keep_delays(pm.run(arms[nm]), keep)
        c.name = nm
        frozen.append(c)
    return frozen, {"ytail": ytail, "xtail": xtail, "taus_dt": taus_dt, "taus_us": [t * dt * 1e6 for t in taus_dt]}


# ── the structural assertions, at the freeze and again at submission ──────

def per_qubit(tc):
    """{physical qubit: [(name, qubits, params, duration_dt or None)]} in program order, measures included."""
    d = {}
    for inst in tc.data:
        name = inst.operation.name
        if name == "barrier":
            continue
        qs = tuple(tc.find_bit(q).index for q in inst.qubits)
        prm = tuple(float(p) for p in inst.operation.params)
        for q in qs:
            d.setdefault(q, []).append((name, qs, prm))
    return d


def assert_structure(circs, freeze):
    """Every turned variant carries rz(phi) as the first operation on its turned qubit and rz(-phi) as its last before
    the measurement, and is otherwise the k = 0 circuit on every qubit; every Q circuit carries its idle; the forest
    and, in Q, the two idle path qubits carry nothing but a delay and a measurement; RO prepares with x alone; the
    clbit map is (q0, q2, q1, forest). Raises on the first violation."""
    path, forest = freeze["path"], freeze["forest"]
    clmap = {path[0]: 0, path[2]: 1, path[1]: 2}
    clmap.update({f: 3 + i for i, f in enumerate(forest)})
    by = {c.name: c for c in circs}
    for nm, c in by.items():
        meas = {}
        for inst in c.data:
            if inst.operation.name == "measure":
                meas[c.find_bit(inst.qubits[0]).index] = c.find_bit(inst.clbits[0]).index
        if meas != clmap:
            raise AssertionError(f"{nm}: measure map {meas}")
        pq = per_qubit(c)
        for f in forest:
            if [o[0] for o in pq.get(f, []) if o[0] != "delay"] != ["measure"]:
                raise AssertionError(f"{nm}: forest qubit {f} carries {pq.get(f)}")
        used = set(pq) - set(path) - set(forest)
        if used:
            raise AssertionError(f"{nm}: operations on {sorted(used)}")
        # every measurement starts at the same time (the forest is read with the path)
        starts = set()
        for q, ops in pq.items():
            t = 0
            for name, qs, prm in ops:
                if name == "measure":
                    starts.add(t)
                    break
                t += duration_dt(c, name, qs, prm, freeze)
        if len(starts) != 1:
            raise AssertionError(f"{nm}: measurements start at {sorted(starts)} dt")
    for nm, c in by.items():
        var, _, rest = nm.partition("_")
        if var[0] in "AB":
            role, k = var[1], int(var[2])
            if k == 0:
                continue
            ref = per_qubit(by[f"{var[0]}50_{rest}"])
            pq = per_qubit(c)
            tq = path[ROLE_LOGICAL[role]]
            phi = k * math.pi / 2
            for q in set(ref) | set(pq):
                a = [o for o in pq.get(q, []) if o[0] != "delay"]
                b = [o for o in ref.get(q, []) if o[0] != "delay"]
                if q != tq:
                    if a != b:
                        raise AssertionError(f"{nm}: qubit {q} differs from the unturned circuit")
                    continue
                if not (a[0][0] == "rz" and abs(a[0][2][0] - phi) < 1e-12):
                    raise AssertionError(f"{nm}: first operation on {q} is {a[0]}, not rz({phi})")
                if not (a[-2][0] == "rz" and abs(a[-2][2][0] + phi) < 1e-12 and a[-1][0] == "measure"):
                    raise AssertionError(f"{nm}: last operations on {q} are {a[-2:]}, not rz({-phi}), measure")
                if a[1:-2] != b[:-1] or a[-1] != b[-1]:
                    raise AssertionError(f"{nm}: qubit {q} between the turns differs from the unturned circuit")
        elif var == "Q":
            lq, P, t = int(rest[0]), rest[1], int(rest[2])
            tail = freeze["xtail"] if P == "X" else freeze["ytail"]
            pq = per_qubit(c)
            ops = pq[path[lq]]
            want = [("sx", ())] + ([("delay", (float(freeze["taus_dt"][t]),))] if freeze["taus_dt"][t] else []) \
                + [(n_, tuple(p_)) for n_, p_ in tail] + [("measure", ())]
            got = [(o[0], o[2]) for o in ops]
            if got != want:
                raise AssertionError(f"{nm}: qubit {path[lq]} carries {got}, wanted {want}")
            for l2 in range(3):
                if l2 != lq and [o[0] for o in pq[path[l2]] if o[0] != "delay"] != ["measure"]:
                    raise AssertionError(f"{nm}: idle path qubit {path[l2]} carries {pq[path[l2]]}")
        elif var == "RO":
            pq = per_qubit(c)
            for l, bit in ((0, rest[0]), (2, rest[1]), (1, rest[2])):
                want = (["x"] if bit == "1" else []) + ["measure"]
                if [o[0] for o in pq[path[l]] if o[0] != "delay"] != want:
                    raise AssertionError(f"{nm}: path qubit {path[l]} carries {pq[path[l]]}")
        else:
            raise AssertionError(f"unknown circuit {nm}")


def duration_dt(c, name, qs, prm, freeze):
    if name == "delay":
        return int(round(prm[0]))
    if name == "rz":
        return 0
    if name in ("sx", "x"):
        return freeze["durations_dt"][name][str(qs[0])]
    if name == "cz":
        a, b = sorted(qs)
        return freeze["durations_dt"]["cz"][f"{a}-{b}"]
    if name == "measure":
        return freeze["durations_dt"]["measure"][str(qs[0])]
    raise ValueError(name)


def durations_of(be, path, forest):
    """The target's durations in dt of every operation the frozen circuits carry."""
    T, dt = be.target, be.dt
    return {"sx": {str(q): int(round(T["sx"][(q,)].duration / dt)) for q in path},
            "x": {str(q): int(round(T["x"][(q,)].duration / dt)) for q in path},
            "cz": {f"{min(a, b)}-{max(a, b)}": int(round(T["cz"][bond_key(T, a, b)].duration / dt))
                   for a, b in zip(path, path[1:])},
            "measure": {str(q): int(round(T["measure"][(q,)].duration / dt)) for q in path + forest}}


def oldest_property_utc(freeze):
    """The date of the oldest guarded property in the freeze snapshot (qubits, sx and bonds)."""
    snap = freeze["snapshot"]
    taken = datetime.fromisoformat(snap["taken_utc"])
    ages = [r["oldest_age_h"] for r in snap["path"].values()] + [b["age_h"] for b in snap["bonds"].values()]
    from datetime import timedelta
    return taken - timedelta(hours=max(ages))


def expected_names(chip):
    names = []
    for role, k in TURNS:
        for b0 in "XYZ":
            for b2 in "XYZ":
                names.append(f"A{role}{k}_{b0}{b2}")
                if chip in RIDERS:
                    names.append(f"B{role}{k}_{b0}{b2}")
    if chip in RIDERS:
        names += [f"Q_{l}{P}{t}" for l in range(3) for P in "XY" for t in range(len(TAUS_US))]
    names += [f"RO_{a}{b}{c}" for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    return sorted(names)


def main(chip):
    if chip not in CHIPS:
        raise SystemExit(__doc__)
    svc = service()
    be = svc.backend(chip)
    props = be.properties()
    now = datetime.now(timezone.utc)
    path, choice = choose_path(chip, be, props, now)
    cm = be.coupling_map
    forest = sorted({n for q in path for n in cm.neighbors(q)} - set(path))
    T = be.target
    dt = be.dt
    dur = durations_of(be, path, forest)
    frozen, extra = build(chip, be, path, forest)
    names = [c.name for c in frozen]
    if names != expected_names(chip):
        raise SystemExit("the frozen names are not the expected set")
    qpy_path, meta_path = files(chip)
    recs, bonds, fails = path_check(props, T, path, now, FREEZE_MAX_AGE_H)
    forest_recs = {}
    for f in forest:
        forest_recs[f], _ = qubit_record(props, f, now)
    freeze = {"chip": chip, "path": path, "forest": forest, "num_qubits": be.num_qubits, "dt_s": dt,
              "durations_dt": dur, **extra, "path_choice": choice,
              "snapshot": {"calibration": str(props.last_update_date), "taken_utc": now.isoformat(),
                           "path": {str(q): r for q, r in recs.items()}, "bonds": bonds,
                           "forest": {str(q): r for q, r in forest_recs.items()}},
              "clbits": {"q0": path[0], "q2": path[2], "q1": path[1],
                         **{f"forest{i}": f for i, f in enumerate(forest)}},
              "shuffle_seed": SHUFFLE_SEED[chip], "shots_per_half": 4096, "n_circuits": len(frozen),
              "n_pubs": 2 * len(frozen), "turns": [list(t) for t in TURNS],
              "record_dir": f"data/ibm_drive_frame_{now.strftime('%B').lower()}{now.year}"}
    if fails:
        raise SystemExit(f"{chip}: the chosen path fails {fails}")
    assert_structure(frozen, freeze)
    with open(qpy_path, "wb") as f:
        qpy.dump(frozen, f)
    freeze["sha256"] = hashlib.sha256(qpy_path.read_bytes()).hexdigest()
    freeze["oldest_property_utc"] = oldest_property_utc(freeze).isoformat()
    with open(qpy_path, "rb") as f:
        assert_structure(qpy.load(f), freeze)
    meta_path.write_text(json.dumps(freeze, indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"freeze record canonical sha256 {canonical_sha(json.loads(meta_path.read_text(encoding='utf-8')))}; "
          f"oldest dated property {freeze['oldest_property_utc']}")
    print(f"{chip}: path {path} ({choice['rule']}), forest {forest}; froze {len(frozen)} circuits, "
          f"{2 * len(frozen)} PUBs, sha256 {freeze['sha256']}")
    print("durations (dt)", dur, "idles (dt)", extra["taus_dt"])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
