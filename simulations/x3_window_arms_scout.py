"""Scout for experiments/X3_WINDOW_HARDWARE_PREDICTION.md: a sign twin and a cz twirl through the letter replay.

YZ+ZY (or truly) on Marrakesh [4, 5, 6], compiled today against FakeMarrakesh (level 3, seed 11) and scheduled as
late as possible: from |+,-,+> (anchor) and from |-,+,-> (twin); each compiled circuit also with a random Pauli pair
before every cz and its conjugate after, single-qubit runs merged (twirl, K instances averaged); and with an idle of
4 us on the end qubits before the basis change (delay). Replayed with April's T1, T2 and one letter at a time; printed
are the forbidden single-site readings <Y0>, <Z0>, <Y2>, <Z2> pooled over bases (and <X0>, <X2>). The twin changes the
sign of some single-site readings and keeps others (under a cz phase of 0.05, <Z2> -0.34 -> +0.11 while <Z0> stays near
-0.15), and a twirl suppresses single-qubit letters as well as cz letters, so neither sorts the letter classes by a sign
or a suppression; neither is flown.

    python simulations/x3_window_arms_scout.py yzzy 4  >  simulations/results/x3_window/x3_window_arms_scout.txt
"""
import math
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "flight"))
import reality_window_april_letters as ap  # noqa: E402
import reality_window_letters as rl  # noqa: E402

from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import ALAPScheduleAnalysis, PadDelay, Optimize1qGatesDecomposition
from qiskit_ibm_runtime.fake_provider import FakeMarrakesh

BACK = FakeMarrakesh()
PATH = [4, 5, 6]
TERMS = {"yzzy": [("Y", "Z"), ("Z", "Y")], "truly": [("X", "X"), ("Y", "Y")]}


def ham(cls):
    ps = []
    for (i, j) in ((0, 1), (1, 2)):
        for (a, b) in TERMS[cls]:
            ch = ["I"] * 3; ch[2 - i] = a; ch[2 - j] = b
            ps.append("".join(ch))
    return SparsePauliOp(ps, coeffs=np.ones(len(ps)))


def logical(cls, state, b0, b2, delay_us=0.0):
    cr = ClassicalRegister(2, "meas")
    qc = QuantumCircuit(3); qc.add_register(cr)
    signs = {"pmp": (0, 1, 0), "mpm": (1, 0, 1)}[state]
    for q in range(3):
        qc.h(q)
        if signs[q]:
            qc.z(q)
    g = PauliEvolutionGate(ham(cls), time=0.8 / 3)
    for _ in range(3):
        qc.append(g, range(3))
    if delay_us:
        qc.barrier()
        qc.delay(delay_us * 1000, [0, 2], unit="ns")
        qc.barrier()
    for q, bb in ((0, b0), (2, b2)):
        if bb == "X":
            qc.h(q)
        elif bb == "Y":
            qc.sdg(q); qc.h(q)
    qc.measure(0, cr[0]); qc.measure(2, cr[1])
    return qc


CZ_CONJ = {}  # (Pa, Pb) before cz -> (Pa', Pb') after, up to sign: cz X_a cz = X_a Z_b


def _conj(pa, pb):
    # Paulis as (x, z) bits; cz maps x_a -> x_a, z_b ^= x_a; x_b -> x_b, z_a ^= x_b
    xa, za = pa; xb, zb = pb
    return (xa, za ^ xb), (xb, zb ^ xa)


def twirl(tc, rng):
    """Insert a random Pauli pair before each cz and its conjugate after; merge 1q runs."""
    out = tc.copy_empty_like()
    names = {(0, 0): None, (1, 0): "x", (0, 1): "z", (1, 1): "y"}
    for inst in tc.data:
        if inst.operation.name == "cz":
            a, b = inst.qubits
            pa = (int(rng.integers(2)), int(rng.integers(2))); pb = (int(rng.integers(2)), int(rng.integers(2)))
            qa, qb = _conj(pa, pb)
            for p, q in ((pa, a), (pb, b)):
                if names[p]:
                    getattr(out, names[p])(q)
            out.append(inst)
            for p, q in ((qa, a), (qb, b)):
                if names[p]:
                    getattr(out, names[p])(q)
        else:
            out.append(inst)
    pm = PassManager([Optimize1qGatesDecomposition(basis=["rz", "sx", "x"])])
    return pm.run(out)


def to_ops(tc):
    pm = PassManager([ALAPScheduleAnalysis(target=BACK.target), PadDelay(target=BACK.target)])
    tc = pm.run(tc)
    m = {p: i for i, p in enumerate(PATH)}
    ops = []
    for inst in tc.data:
        name = inst.operation.name
        qs = [tc.find_bit(q).index for q in inst.qubits]
        if name in ("measure", "barrier") or not all(q in m for q in qs):
            continue
        prm = [float(p) for p in inst.operation.params]
        d = prm[0] * ap.DUR["dt"] if name == "delay" else ap.DUR.get(name, 0.0)
        if name == "delay" and inst.operation.unit != "dt":
            d = prm[0] * {"ns": 1e-3, "us": 1.0, "s": 1e6}[inst.operation.unit]
        ops.append((name, [m[q] for q in qs], prm, d))
    last = max(i for i, o in enumerate(ops) if o[0] == "cz" and 1 in o[1])
    ops = [o for i, o in enumerate(ops) if not (i > last and o[1] == [1])]
    return ops


def compile_(qc, seed=11):
    return transpile(qc, BACK, optimization_level=3, initial_layout=PATH, seed_transpiler=seed)


def readings(opslist, L):
    """opslist: {(b0, b2): [ops, ...]} (several = twirl instances, averaged). Pooled single-site readings."""
    acc = {}
    for (b0, b2), variants in opslist.items():
        p = np.mean([rl.probs_rho(rl.replay(o, PATH, ap.T1, ap.T2, L)) for o in variants], axis=0)
        par = rl.parities_probs(p)
        acc.setdefault(("0", b0), []).append(par[(1, 0, 0)])
        acc.setdefault(("2", b2), []).append(par[(0, 1, 0)])
    return {k: float(np.mean(v)) for k, v in acc.items()}


MECH = {
    "ideal(T1,T2)": {},
    "det end +0.04": {("det", 4): 0.04, ("det", 6): 0.04},
    "sxt end +0.06": {("sxt", 4): 0.06, ("sxt", 6): 0.06},
    "sxp end +0.02": {("sxp", 4): 0.02, ("sxp", 6): 0.02},
    "sxa end +0.01": {("sxa", 4): 0.01, ("sxa", 6): 0.01},
    "czp all +0.05": {("czp", p, q): 0.05 for p in rl.PAIRS for q in p},
    "czk all +0.05": {("czk", p): 0.05 for p in rl.PAIRS},
    "stk +0.05": {("stk", (4, 5), 6): 0.05, ("stk", (5, 6), 4): 0.05},
}

if __name__ == "__main__":
    cls = sys.argv[1] if len(sys.argv) > 1 else "yzzy"
    rng = np.random.default_rng(7)
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    bases = [(a, b) for a in "XYZ" for b in "XYZ"]
    diag = [("X", "X"), ("Y", "Y"), ("Z", "Z")]
    arms = {}
    for state in ("pmp", "mpm"):
        comp = {bb: compile_(logical(cls, state, *bb)) for bb in bases}
        arms[f"anchor {state}"] = {bb: [to_ops(c)] for bb, c in comp.items()}
        arms[f"twirl  {state}"] = {bb: [to_ops(twirl(c, rng)) for _ in range(K)] for bb, c in comp.items()}
        arms[f"delay4 {state}"] = {bb: [to_ops(compile_(logical(cls, state, *bb, delay_us=4.0)))] for bb in diag}
    print(f"class {cls}; readings <Y0> <Z0> <Y2> <Z2> pooled over bases (twirl: {K} instances)")
    for mname, L in MECH.items():
        print(f"\n  {mname}")
        for aname, ol in arms.items():
            r = readings(ol, L)
            print(f"    {aname:12s}  Y0 {r[('0','Y')]:+.4f}  Z0 {r[('0','Z')]:+.4f}  Y2 {r[('2','Y')]:+.4f}  Z2 {r[('2','Z')]:+.4f}"
                  f"   X0 {r[('0','X')]:+.3f} X2 {r[('2','X')]:+.3f}")
