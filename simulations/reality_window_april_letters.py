"""The 2026-04-30 Marrakesh F83 run read through the reality window, on the circuits actually flown.

Reproduces section 2 of experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md:
  (a) the 36 flown circuits (retrieved from IBM, job d7pol1e7g7gs73cf7j90) implement the runner's
      first-order Trotter circuits exactly;
  (b) with ideal gates and the day's T1/T2, every reading the class forbids is near zero, while the
      hardware reads up to 6σ there;
  (c) one set of physical error letters, shared by the four classes (the library of
      simulations/flight/reality_window_letters.py), fitted to the forbidden readings of all four:
      the YZ+ZY window closes, the XY+YX ⟨Y₀Z₂⟩ does not;
  (d) static ZZ alone, fitted to XY+YX alone, and the size it needs.
The April job carried no delays, so where IBM placed the gates in time is not recorded; every fit is
run under both schedules, as late as possible (what the 2026-10-28 flight pins explicitly) and as soon
as possible, because letters on an end qubit and static ZZ read differently under the two.

    python simulations/reality_window_april_letters.py
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "flight"))
import reality_window_letters as rl  # noqa: E402

DATA = HERE.parent / "data" / "ibm_f83_signature_april2026"
PATH = [4, 5, 6]
T1 = [206.0, 144.0, 351.0]          # us, 2026-04-30 calibration (data README)
T2 = [184.0, 121.0, 151.0]
# readout of 2026-04-30 (data README): error e and asymmetry d = P(0|1) - P(1|0) on q4 and q6; q5 not read
RO = [(0.0042 - 0.0059 / 2, 0.0042 + 0.0059 / 2), (0.0037 - 0.0039 / 2, 0.0037 + 0.0039 / 2)]   # (p01, p10)
DUR = {"sx": 0.036, "x": 0.036, "cz": 0.068, "dt": 0.004}   # us, ibm_marrakesh target
TERMS = {"truly_unbroken": [("X", "X"), ("Y", "Y")], "pi2_odd_pure": [("X", "Y"), ("Y", "X")],
         "pi2_even_nontruly": [("Y", "Z"), ("Z", "Y")], "mixed_anti_one_sixth": [("X", "Y"), ("Y", "Z")]}
# sector steps (k_X, n_Y) of each class's Hamiltonian letters; a reading outside the subgroup they generate
# from (0, 0) is forbidden
STEPS = {"truly_unbroken": [(0, 1)], "pi2_odd_pure": [(1, 0)], "pi2_even_nontruly": [], "mixed_anti_one_sixth": [(1, 0)]}
KEYS = [(a, b) for a in "IXYZ" for b in "IXYZ" if (a, b) != ("I", "I")]


def sector(a, b):
    return (sum(x in "YZ" for x in (a, b)) % 2, sum(x == "Y" for x in (a, b)) % 2)


def forbidden(cls):
    S = {(0, 0)}
    for s in STEPS[cls]:
        S |= {((x + s[0]) % 2, (y + s[1]) % 2) for (x, y) in S}
    return [k for k in KEYS if sector(*k) not in S]


def measured():
    """Pauli readings with the basis-averaged estimator for singles; σ from the binomial count."""
    d = json.loads(next(DATA.glob("f83_signature_*.json")).read_text(encoding="utf-8"))
    out = {}
    for cls, bases in d["counts"].items():
        for (a, b) in KEYS:
            s = n = 0
            for basis, cnt in bases.items():
                if (a != "I" and basis[0] != a) or (b != "I" and basis[1] != b):
                    continue
                for k, c in cnt.items():
                    q2, q0 = int(k[0]), int(k[1])
                    s += c * (-1) ** (q0 * (a != "I") + q2 * (b != "I")); n += c
            e = s / n
            out[(cls, a, b)] = (e, math.sqrt(max(1 - e * e, 1e-9) / n))
    return out


def flown_ops(schedule):
    """The flown circuits as (name, logical qubits, params, duration_us); schedule 'alap' pads explicit
    delays (the replay then follows them), 'asap' leaves the gaps for the replay's as-soon-as-possible clock."""
    from qiskit import qpy
    from qiskit.transpiler import PassManager
    from qiskit.transpiler.passes import ALAPScheduleAnalysis, PadDelay
    from qiskit_ibm_runtime.fake_provider import FakeMarrakesh
    with open(DATA / "flown_circuits_d7pol1e7g7gs73cf7j90.qpy", "rb") as f:
        circs = qpy.load(f)
    pm = PassManager([ALAPScheduleAnalysis(target=FakeMarrakesh().target), PadDelay(target=FakeMarrakesh().target)])
    out = {}
    for tc in circs:
        cls, b = tc.name.rsplit("_", 1)
        if schedule == "alap":
            tc = pm.run(tc)
        m = {p: i for i, p in enumerate(PATH)}
        ops = []
        for inst in tc.data:
            name = inst.operation.name
            qs = [tc.find_bit(q).index for q in inst.qubits]
            if name in ("measure", "barrier") or not all(q in m for q in qs):
                continue
            prm = [float(p) for p in inst.operation.params]
            d = prm[0] * DUR["dt"] if name == "delay" else DUR.get(name, 0.0)
            ops.append((name, [m[q] for q in qs], prm, d))
        # April did not read q1, so the as-late-as-possible schedule pushes q1's gates after its last cz past the
        # 2.6 us measurement of q0 and q2; those gates cannot change q0's and q2's statistics, and keeping them
        # would run the replay's global clock, and q0's and q2's decay, through their own readout. Dropped.
        last = max(i for i, o in enumerate(ops) if o[0] == "cz" and 1 in o[1])
        ops = [o for i, o in enumerate(ops) if not (i > last and o[1] == [1])]
        out[(cls, b[0], b[1])] = ops
    return out


def logical_probs(cls, b0, b2):
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import SparsePauliOp, Statevector
    ps = []
    for (i, j) in ((0, 1), (1, 2)):
        for (a, b) in TERMS[cls]:
            ch = ["I"] * 3; ch[2 - i] = a; ch[2 - j] = b
            ps.append("".join(ch))
    qc = QuantumCircuit(3)
    for q in range(3):
        qc.h(q)
        if q % 2:
            qc.z(q)
    g = PauliEvolutionGate(SparsePauliOp(ps, coeffs=np.ones(len(ps))), time=0.8 / 3)
    for _ in range(3):
        qc.append(g, range(3))
    for q, bb in ((0, b0), (2, b2)):
        if bb == "X":
            qc.h(q)
        elif bb == "Y":
            qc.sdg(q); qc.h(q)
    return Statevector(qc.decompose(gates_to_decompose=["PauliEvolution"])).probabilities([0, 2])


def ops_probs(ops):
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(3)
    for n, q, p, _ in ops:
        if n != "delay":
            getattr(qc, n)(*p, *q) if p else getattr(qc, n)(*q)
    return Statevector(qc).probabilities([0, 2])


def readings(OPS, L):
    """Pauli readings on (q0, q2) per class, singles averaged over the three bases that carry them."""
    out = {}
    for cls in TERMS:
        acc = {}
        for b0 in "XYZ":
            for b2 in "XYZ":
                p8 = rl.assignment_product([RO[0], RO[1], (0.0, 0.0)]) @ rl.probs_rho(rl.replay(OPS[(cls, b0, b2)], PATH, T1, T2, L))
                par = rl.parities_probs(p8)
                acc.setdefault((b0, b2), []).append(par[(1, 1, 0)])
                acc.setdefault((b0, "I"), []).append(par[(1, 0, 0)])
                acc.setdefault(("I", b2), []).append(par[(0, 1, 0)])
        for k, v in acc.items():
            out[(cls,) + k] = float(np.mean(v))
    return out


def main():
    for schedule in ("alap", "asap"):
        print(f"== schedule {schedule}")
        run(flown_ops(schedule))


def run(OPS):
    worst = max(float(np.max(np.abs(ops_probs(OPS[k]) - logical_probs(*k)))) for k in OPS)
    print(f"(a) flown vs logical Trotter, max |dp| over 36 circuits: {worst:.1e}")
    M = measured()
    FORB = {c: forbidden(c) for c in TERMS}

    def resid(L, classes=tuple(TERMS)):
        R = readings(OPS, L)
        return np.array([(M[(c,) + k][0] - R[(c,) + k]) / M[(c,) + k][1] for c in classes for k in FORB[c]])

    def per_class(L):
        R = readings(OPS, L)
        return {c: sum(((M[(c,) + k][0] - R[(c,) + k]) / M[(c,) + k][1]) ** 2 for k in FORB[c]) for c in TERMS}, R

    pc, R0 = per_class({})
    print("(b) ideal gates, forbidden-reading chi2 per class:", {c: round(v, 1) for c, v in pc.items()},
          f"| XY+YX <Y0Z2>: model {R0[('pi2_odd_pure', 'Y', 'Z')]:+.4f}, hardware {M[('pi2_odd_pure', 'Y', 'Z')][0]:+.4f}")
    letters = rl.LETTERS
    f = least_squares(lambda x: resid(dict(zip(letters, x))), np.zeros(len(letters)), x_scale=0.01)
    L = dict(zip(letters, f.x))
    pc, R = per_class(L)
    print(f"(c) all {len(letters)} letters jointly, forbidden chi2 per class:", {c: round(v, 1) for c, v in pc.items()},
          f"| XY+YX <Y0Z2>: model {R[('pi2_odd_pure', 'Y', 'Z')]:+.4f}")
    print("    letters:", {str(k): round(v, 4) for k, v in L.items()})
    zz = [("zz", p) for p in rl.PAIRS]
    g = least_squares(lambda x: resid(dict(zip(zz, x)), ("pi2_odd_pure",)), np.full(2, 0.05))
    pc, R = per_class(dict(zip(zz, g.x)))
    print(f"(d) static ZZ alone on XY+YX: zeta = {np.round(g.x, 3)} rad/us, XY+YX chi2 {pc['pi2_odd_pure']:.1f}, "
          f"<Y0Z2> model {R[('pi2_odd_pure', 'Y', 'Z')]:+.4f}")


if __name__ == "__main__":
    main()
