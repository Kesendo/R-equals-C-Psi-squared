"""Freeze the reality-window circuits (needs the ibm_marrakesh target, read only; no job is submitted).

    python reality_window_freeze.py

Arms (section 3 of experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md):
  A      the 2026-04-30 flown XY+YX circuits verbatim, then q1 read in Y
  R      the same gates with physical 4 and 6 exchanged
  D0..D6 the logical circuit at n_trotter 0, 1, 2, 3, 4, 6, compiled today (level 3, seed 11)
  Kp, Km arm A with a virtual rz(+-0.2) on q1 right after its sixth cz
  RO     the eight computational states of (q0, q2, q1), prepared with x, for the readout matrix
Every circuit is then scheduled as late as possible with explicit delays on the three path qubits,
so the backend has no timing left to decide.
"""
import hashlib
import json
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile, qpy
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import ALAPScheduleAnalysis, PadDelay

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
APRIL_QPY = REPO / "data" / "ibm_f83_signature_april2026" / "flown_circuits_d7pol1e7g7gs73cf7j90.qpy"
OUT = HERE / "reality_window_circuits.qpy"
META = HERE / "reality_window_freeze.json"
PATH = [4, 5, 6]
SEED = 11
THETA, M_CZ = 0.2, 6


def backend():
    from dotenv import load_dotenv
    from qiskit_ibm_runtime import QiskitRuntimeService
    pipe = Path(r"D:/Entwicklung/Projekte/.NET Projekte/AIEvolution/AIEvolution.UI")
    for p in (HERE / ".env", pipe / ".env", pipe / "experiments" / "ibm_quantum_tomography" / ".env"):
        if p.exists():
            load_dotenv(p)
    tok = os.environ.get("IBM_QUANTUM_TOKEN")
    s = (QiskitRuntimeService(channel="ibm_quantum_platform", token=tok) if tok
         else QiskitRuntimeService(channel="ibm_quantum_platform"))
    return s.backend("ibm_marrakesh")


def logical(n, b0, b2, b1="Y"):
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import SparsePauliOp
    ps = []
    for (i, j) in ((0, 1), (1, 2)):
        for (a, b) in (("X", "Y"), ("Y", "X")):
            ch = ["I"] * 3; ch[2 - i] = a; ch[2 - j] = b
            ps.append("".join(ch))
    qc = QuantumCircuit(3, name=f"n{n}")
    qc.add_register(ClassicalRegister(3, "meas"))
    for q in range(3):
        qc.h(q)
        if q % 2:
            qc.z(q)
    g = PauliEvolutionGate(SparsePauliOp(ps, coeffs=np.ones(len(ps))), time=0.8 / 3)
    for _ in range(n):
        qc.append(g, range(3))
    for q, b in ((0, b0), (2, b2), (1, b1)):
        if b == "X":
            qc.h(q)
        elif b == "Y":
            qc.sdg(q); qc.h(q)
    qc.measure(0, 0); qc.measure(2, 1); qc.measure(1, 2)
    return qc


def from_april(tc, nq, phys_map, ytail, theta=None):
    out = QuantumCircuit(nq); out.add_register(ClassicalRegister(3, "meas"))
    meas, ncz, done = [], 0, theta is None
    p1 = phys_map[PATH[1]]
    for inst in tc.data:
        qs = [phys_map.get(tc.find_bit(q).index, tc.find_bit(q).index) for q in inst.qubits]
        if inst.operation.name == "measure":
            meas.append((qs[0], tc.find_bit(inst.clbits[0]).index)); continue
        if inst.operation.name in ("barrier", "delay"):
            continue
        out.append(inst.operation, qs)
        if inst.operation.name == "cz":
            ncz += 1
            if not done and ncz == M_CZ:
                out.rz(theta, p1); done = True
    assert done
    for n, p in ytail:
        getattr(out, n)(*p, p1) if p else getattr(out, n)(p1)
    for q, c in meas:
        out.measure(q, c)
    out.measure(p1, 2)
    return out


def readout_circuit(nq, b0, b2, b1):
    out = QuantumCircuit(nq); out.add_register(ClassicalRegister(3, "meas"))
    for q, b in ((PATH[0], b0), (PATH[2], b2), (PATH[1], b1)):
        if b:
            out.x(q)
    out.measure(PATH[0], 0); out.measure(PATH[2], 1); out.measure(PATH[1], 2)
    return out


def strip_idle_delays(tc):
    """Keep the explicit schedule on the path, drop the padding on the 153 untouched qubits."""
    out = tc.copy_empty_like()
    for inst in tc.data:
        qs = [tc.find_bit(q).index for q in inst.qubits]
        if inst.operation.name == "delay" and qs[0] not in PATH:
            continue
        out.append(inst)
    return out


def main():
    be = backend()
    nq = be.num_qubits
    T = be.target
    dur = {"sx": {q: T["sx"][(q,)].duration for q in PATH}, "x": {q: T["x"][(q,)].duration for q in PATH},
           "cz": {f"{a}-{b}": T["cz"][(a, b)].duration for a, b in ((4, 5), (5, 6))}, "dt": be.dt}
    yq = QuantumCircuit(1); yq.sdg(0); yq.h(0)
    ytail = [(i.operation.name, [float(p) for p in i.operation.params])
             for i in transpile(yq, basis_gates=["rz", "sx", "x"], optimization_level=1).data]
    with open(APRIL_QPY, "rb") as f:
        april = {c.name: c for c in qpy.load(f)}
    arms = {}
    ident, mirror = {q: q for q in PATH}, {4: 6, 5: 5, 6: 4}
    for b0 in "XYZ":
        for b2 in "XYZ":
            tc = april[f"pi2_odd_pure_{b0}{b2}"]
            arms[f"A_{b0}{b2}"] = from_april(tc, nq, ident, ytail)
            arms[f"R_{b0}{b2}"] = from_april(tc, nq, mirror, ytail)
            arms[f"Kp_{b0}{b2}"] = from_april(tc, nq, ident, ytail, +THETA)
            arms[f"Km_{b0}{b2}"] = from_april(tc, nq, ident, ytail, -THETA)
            for n in (0, 1, 2, 3, 4, 6):
                arms[f"D{n}_{b0}{b2}"] = transpile(logical(n, b0, b2), be, optimization_level=3,
                                                   initial_layout=PATH, seed_transpiler=SEED)
    for b0 in (0, 1):
        for b2 in (0, 1):
            for b1 in (0, 1):
                arms[f"RO_{b0}{b2}{b1}"] = readout_circuit(nq, b0, b2, b1)
    pm = PassManager([ALAPScheduleAnalysis(target=T), PadDelay(target=T)])
    names = sorted(arms)
    frozen = []
    for nm in names:
        c = strip_idle_delays(pm.run(arms[nm]))
        c.name = nm
        frozen.append(c)
    with open(OUT, "wb") as f:
        qpy.dump(frozen, f)
    h = hashlib.sha256(OUT.read_bytes()).hexdigest()
    p = be.properties()
    calib = {"calibration": str(p.last_update_date),
             "t1_us": {q: p.t1(q) * 1e6 for q in PATH}, "t2_us": {q: p.t2(q) * 1e6 for q in PATH},
             "p01": {q: p.qubit_property(q)["prob_meas1_prep0"][0] for q in PATH},
             "p10": {q: p.qubit_property(q)["prob_meas0_prep1"][0] for q in PATH}}
    META.write_text(json.dumps({"sha256": h, "n_circuits": len(frozen), "durations_s": dur, "ytail": ytail,
                                "theta": THETA, "m_cz": M_CZ, "seed": SEED, "calibration": calib},
                               indent=1, default=str), encoding="utf-8")
    print("froze", len(frozen), "circuits, sha256", h)
    print("durations", dur)


if __name__ == "__main__":
    main()
