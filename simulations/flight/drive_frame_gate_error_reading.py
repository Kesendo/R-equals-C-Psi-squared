"""A post-flight reading: the gate errors the flight's replay leaves out, and the baseline at the day's coherence.

A reading without a verdict, written after the flight; none of its numbers is a pre-registered quantity of
experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md or experiments/X3_WINDOW_HARDWARE_PREDICTION.md. The flight's replay
(reality_window_letters.replay through run_drive_frame) carries T1 and T2 and the letters of its library, and no
incoherent gate error: a cz is an ideal unitary. This script asks three things of that gap, read-only, printing only.

1. The colouring X^⊗3 (the X^⊗3 window page's section 5). B50_XX is run on the three path qubits with every sx and x
   followed by a one-qubit depolarizing channel and every cz by a two-qubit one, at the error the freeze snapshot
   reports for that gate (for x, which it does not report, the sx error; depolarizing parameter d·e/(d − 1), d = 2^n, as Aer's device model converts a gate error
   with no relaxation), with ideal readout and no relaxation. Its moments, in the normalised signs of
   drive_frame_colouring_reading.py (U0 = X0, U1 = -X1, U2 = X2), stand beside the measured ones and beside the
   replay under the snapshot and under arm Q's T2; the two are combined to first order (the replay's moments plus the
   gates' departures from the ideal), which double-counts the part of a gate's error that is decoherence during it.
2. The pre-registered zeros under the same channels. S (the drive-frame page's <Y0Z2> of A50_YZ) and the 14 forbidden
   readings of B5(0) are read on the depolarized state at one and at seven times the reported errors.
3. The baseline at another coherence. The analysis is rerun unchanged with the chip's T2 replaced by arm Q's
   (T2 = 1/g held at or below 2 T1, as the runner's colouring forms it) and with T1 and T2 divided by 7, and s_A,
   its z and W are printed beside the record's.
4. Arm Q's g with its fit uncertainty (the runner's fit, its Jacobian's g entry inflated as it inflates sigma_delta),
   and the colouring's replay at g - sigma_g and g + sigma_g.
5. A population p left in |1> at the start on one path qubit, which the RO arm books as readout (A_RO = A_true M(p),
   M flipping that qubit's bit with probability p): the corrected distribution is then M(p)^-1 applied to the state's,
   whatever A_true is. B50_XX with the gate errors, for p up to the RO arm's P(0|1) on that qubit.
6. A forward check with nothing inverted: two models, the replay with T1 and T2 divided by 7, and the replay with arm
   Q's T2 plus the depolarizing gate errors (to first order), each mapped through the day's RO matrix, against the raw counts' covariance of q0 and q2, third cumulant and rate of all
   three failing, with sigma from resampling the raw counts.
7. The raw counts of all three failing, per circuit half of B's XX circuits, against one common rate (chi^2): a burst
   that struck one PUB would gather them there.

    python drive_frame_gate_error_reading.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, Pauli, SuperOp
from qiskit_aer.noise import depolarizing_error

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
import run_drive_frame as R  # noqa: E402  read-only use
import drive_frame_colouring_reading as CR  # noqa: E402  read-only use: readings, replay

DATA = REPO / "data" / "ibm_drive_frame_november2026"
RESULTS = {"ibm_fez": "drive_frame_ibm_fez_20261113_203626.json",
           "ibm_marrakesh": "drive_frame_ibm_marrakesh_20261113_203627.json",
           "ibm_kingston": "drive_frame_ibm_kingston_20261113_203653.json"}
COLS = ["U0", "U1", "U2", "cov01", "cov12", "cov02", "k3", "allfail", "L3", "L1"]


def depolarized(chip, name, scale, start="000"):
    """The three path qubits' state after the frozen circuit `name`, each sx and cz followed by its depolarizing
    channel at `scale` times the snapshot's error, each x at the sx error (the snapshot reports none for x); delays and
    measurements dropped. Index i = path position; `start` is a Qiskit label over the path positions (q2 q1 q0)."""
    c = chip.circs[name]
    snap = chip.freeze["snapshot"]
    idx = {c.find_bit(c.qubits[q]).index: i for i, q in enumerate(chip.path)}
    e1 = {i: snap["path"][str(q)]["sx"] for i, q in enumerate(chip.path)}
    e2 = {frozenset(chip.path.index(int(x)) for x in k.split("-")): v["cz"] for k, v in snap["bonds"].items()}
    rho = DensityMatrix.from_label(start)
    for ins in c.data:
        n = ins.operation.name
        if n in ("delay", "measure", "barrier"):
            continue
        qs = [idx[c.find_bit(q).index] for q in ins.qubits]
        sub = QuantumCircuit(len(qs))
        sub.append(ins.operation, list(range(len(qs))))
        rho = rho.evolve(sub, qargs=qs)
        e = e2[frozenset(qs)] if n == "cz" else (e1[qs[0]] if n in ("sx", "x") else 0.0)
        if e and scale:
            d = 2 ** len(qs)
            rho = rho.evolve(SuperOp(depolarizing_error(scale * d * e / (d - 1), len(qs))), qargs=qs)
    return rho


def outcome_probs(rho):
    """Outcome probabilities in the runner's index k = 4 q0 + 2 q2 + q1 (Qiskit's index is q0 + 2 q1 + 4 q2)."""
    p = np.zeros(8)
    for i, v in enumerate(rho.probabilities()):
        p[4 * (i & 1) + 2 * (i >> 2 & 1) + (i >> 1 & 1)] += v
    return p


def reading(rho, mask):
    """A parity (mask over q0, q2, q1) read after the circuit's own basis change: Z on each masked qubit."""
    lab = ["I", "I", "I"]                      # Qiskit labels run q2 q1 q0
    for bit, pos in zip(mask, (0, 2, 1)):
        if bit:
            lab[2 - pos] = "Z"
    return float(rho.expectation_value(Pauli("".join(lab))).real)


def q_fit_with_g(x, y, cov, taus):
    """The runner's q_fit, returning also sigma_g from the same Jacobian, inflated as sigma_delta is."""
    import math
    from scipy.linalg import solve_triangular
    from scipy.optimize import least_squares
    tt = np.array(taus, dtype=float)
    r_obs = np.concatenate([x, y])
    Lc = np.linalg.cholesky(cov)

    def res(v):
        d, th, C, g = v
        e = C * np.exp(-g * tt)
        return solve_triangular(Lc, np.concatenate([e * np.cos(th + d * tt), e * np.sin(th + d * tt)]) - r_obs, lower=True)
    fit = least_squares(res, [0.0, math.atan2(y[0], x[0]), 1.0, 0.0], bounds=R.Q_BOUNDS)
    chi = float(fit.fun @ fit.fun)
    cp = np.linalg.inv(fit.jac.T @ fit.jac)
    infl = max(1.0, math.sqrt(chi / 4))
    return {"delta": float(fit.x[0]), "theta0": float(fit.x[1]), "C": float(fit.x[2]), "g": float(fit.x[3]),
            "sigma_delta": math.sqrt(cp[0, 0]) * infl, "sigma_g": math.sqrt(cp[3, 3]) * infl, "chi2": chi,
            "converged": bool(fit.success), "on_bound": False}


def raw_moments(p):
    """Moments in the normalised signs of probabilities p, without any readout correction."""
    return CR.readings(p)


def row(label, r):
    print(f"  {label:34s}" + " ".join(f"{r[k]:+9.5f}" for k in COLS))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    for name, fname in RESULTS.items():
        data = json.loads((DATA / fname).read_text(encoding="utf-8"))
        chip = R.Chip(name)
        if data.get("qpy_sha256") != R.sha256(chip.qpy):
            raise SystemExit(f"{fname} was taken with other circuits than {chip.qpy}")
        snap = chip.freeze["snapshot"]
        print(f"== {name} path {chip.path}; snapshot errors: sx/x {[round(snap['path'][str(q)]['sx'], 5) for q in chip.path]}, "
              f"cz {[round(v['cz'], 5) for v in snap['bonds'].values()]}")
        comp = R.compact(data["counts"], chip.nf)
        A, nro = R.readout_matrix(comp)
        k_ro = R.ro_kappa(comp)

        print("-- 2. the pre-registered zeros on the depolarized state (ideal readout)")
        for scale in (1, 7):
            line = f"   x{scale}: S {reading(depolarized(chip, 'A50_YZ', scale), (1, 1, 0)):+.2e}"
            if chip.rider:
                vals = [reading(depolarized(chip, f"B50_{b[0]}{b[1]}", scale), m) for b, m in R.B_IDX]
                line += f"; B5(0)'s 14 forbidden readings, largest |.| {max(abs(v) for v in vals):.2e}"
            print(line)

        nulls = R.load_nulls(chip.name, strict=False)
        tc = data.get("turn_ok")
        tc = tc if isinstance(tc, dict) else None
        variants = []
        if chip.rider:
            vq_out = R.vq(chip, comp, A, nro, {}, None, k_ro)
            T2q = [1.0 / max(vq_out[str(l)]["g"], 1.0 / (2 * chip.T1[l])) for l in range(3)]
            variants.append((f"T2 from arm Q {[round(t, 1) for t in T2q]} us", list(chip.T1), T2q))
        variants.append(("T1 and T2 / 7", [t / 7 for t in chip.T1], [t / 7 for t in chip.T2]))
        print("-- 3. the analysis rerun with the baseline at another coherence (readings, not verdicts)")

        def say(lab, o):
            v1 = o["V1"]
            s = f"   {lab:48s} s_A {v1['s_A']:+.4f}  sigma_S {v1['sigma_S']:.4f}  z {v1['z']:+.2f}"
            if chip.rider and "WB" in o:
                s += f"  W {o['WB']['W']:.1f}"
            print(s)
        say("the record's (freeze snapshot)", R.analyze_compact(chip, comp, nulls=nulls, turn_ok=tc))
        T1s, T2s = list(chip.T1), list(chip.T2)
        for lab, T1, T2 in variants:
            chip.T1, chip.T2, chip._ptrue = T1, T2, {}
            say(lab, R.analyze_compact(chip, comp, nulls=nulls, turn_ok=tc))
        chip.T1, chip.T2, chip._ptrue = T1s, T2s, {}

        if chip.rider:
            print("-- 1. the colouring: moments of B's pooled XX circuits (RO-corrected) and of the models")
            print(f"  {'':34s}" + " ".join(f"{k:>9s}" for k in COLS))
            names = [f"{v}_XX" for v in chip.variants("B")]
            pc = sum(comp[n][0].sum(0) for n in names)
            meas = CR.readings(np.linalg.inv(A) @ (pc / pc.sum()))
            g = CR.readings(outcome_probs(depolarized(chip, "B50_XX", 1)))
            snapr = CR.readings(CR.replay(chip, chip.T1, chip.T2))
            qr = CR.readings(CR.replay(chip, chip.T1, variants[0][2]))
            row("measured", meas)
            row("gate errors alone (depolarizing)", g)
            row("replay, freeze snapshot", snapr)
            row("replay, T2 from arm Q", qr)
            for lab, b in (("snapshot", snapr), ("arm Q's T2", qr)):
                comb = {k: b[k] + g[k] - (1.0 if k in ("U0", "U1", "U2") else 0.0) for k in COLS}
                row(f"{lab} + gate errors (1st order)", comb)
                row("  measured minus that", {k: meas[k] - comb[k] for k in COLS})
            print(f"   the corrected distribution's cells (k = 4 q0 + 2 q2 + q1): "
                  f"{[round(float(v), 4) for v in np.linalg.inv(A) @ (pc / pc.sum())]}")

            print("-- 4. arm Q's g with its fit uncertainty, and the colouring's replay at g -/+ sigma_g")
            orig = R.q_fit
            R.q_fit = q_fit_with_g
            try:
                vq_g = R.vq(chip, comp, A, nro, {}, None, k_ro)
            finally:
                R.q_fit = orig
            gs = [vq_g[str(l)]["g"] for l in range(3)]
            sg = [vq_g[str(l)]["sigma_g"] for l in range(3)]
            print(f"   g {[round(v, 4) for v in gs]} +- {[round(v, 4) for v in sg]} per us (a Ramsey decay, no echo)")
            for lab, sgn in (("g - sigma_g", -1), ("g + sigma_g", +1)):
                T2v = [1.0 / max(gs[l] + sgn * sg[l], 1.0 / (2 * chip.T1[l]), 1e-6) for l in range(3)]
                rv = CR.readings(CR.replay(chip, chip.T1, T2v))
                print(f"   replay at {lab}: T2 {[round(t, 1) for t in T2v]} us, -ln|<XXX>| {rv['L3']:.4f}, "
                      f"with the gate errors (1st order) {rv['L3'] + g['L3']:.4f}, cov02 {rv['cov02']:+.4f}")

            print("-- 5. a population p left in |1> at the start, booked as readout (B50_XX with the gate errors)")
            print(f"  {'':34s}" + " ".join(f"{k:>9s}" for k in COLS))
            row("measured", meas)
            ro = {i: (float(np.mean([sum(A[j, k] for j in range(8) if j & b) for k in range(8) if not k & b])),
                      float(np.mean([sum(A[j, k] for j in range(8) if not j & b) for k in range(8) if k & b])))
                  for i, b in ((0, 4), (1, 1), (2, 2))}
            for i, b in ((0, 4), (1, 1), (2, 2)):
                lab_start = ["0", "0", "0"]
                lab_start[2 - i] = "1"
                p1 = outcome_probs(depolarized(chip, "B50_XX", 1, start="".join(lab_start)))
                p0s = outcome_probs(depolarized(chip, "B50_XX", 1))
                for frac in (0.5, 1.0):
                    q = frac * ro[i][1]
                    M = np.zeros((8, 8))
                    for k in range(8):
                        M[k, k] += 1 - q
                        M[k ^ b, k] += q
                    r = CR.readings(np.linalg.solve(M, (1 - q) * p0s + q * p1))
                    row(f"q{chip.path[i]} start p {q:.4f}", r)
            print(f"   (the RO arm's P(1|0), P(0|1) per path qubit: "
                  f"{[(chip.path[i], round(ro[i][0], 4), round(ro[i][1], 4)) for i in range(3)]})")
            # the ideal B outcome is q0 = 0, q2 = 0, q1 = 1 (k = 1): each qubit's false-flip rate with the other two
            # prepared that way, against the B circuits' raw rate of that qubit reading its failed value
            raw = pc / pc.sum()
            for i, b in ((0, 4), (1, 1), (2, 2)):
                fail = 0 if b == 1 else b                  # q1 fails by reading 0, q0 and q2 by reading 1
                ro_fail = sum(A[j, 1] for j in range(8) if (j & b) == fail)
                raw_fail = sum(raw[j] for j in range(8) if (j & b) == fail)
                print(f"   q{chip.path[i]}: RO arm's false-flip rate with the others prepared as B ideally leaves them "
                      f"{ro_fail:.4f}; B's raw failure rate {raw_fail:.4f}")

            print("-- 6. forward check, nothing inverted: models mapped through the day's RO matrix against the raw counts")
            raw = pc / pc.sum()
            rng = np.random.default_rng(20261116)
            boots = [raw_moments(rng.multinomial(int(pc.sum()), raw) / pc.sum()) for _ in range(400)]
            sraw = {k: float(np.std([bb[k] for bb in boots])) for k in ("cov02", "k3", "allfail")}
            mr = raw_moments(raw)
            ideal = outcome_probs(depolarized(chip, "B50_XX", 0))
            gp = outcome_probs(depolarized(chip, "B50_XX", 1))
            p7 = CR.replay(chip, [t / 7 for t in chip.T1], [t / 7 for t in chip.T2])
            pq = CR.replay(chip, chip.T1, variants[0][2])
            for lab, pmod in (("T1, T2 / 7", p7), ("arm Q's T2 + gate errors", pq + gp - ideal)):
                m = raw_moments(A @ pmod)
                print(f"   {lab:26s} raw minus model: " + "  ".join(
                    f"{k} {mr[k] - m[k]:+.4f} ({(mr[k] - m[k]) / sraw[k]:+.1f} sigma)" for k in ("cov02", "k3", "allfail")))

            print("-- 7. raw counts of all three failing, per circuit half of B's XX circuits")
            af, nn = [], []
            for n in names:
                c, _, ns = comp[n]
                for h in range(2):
                    af.append(int(c[h][6]))
                    nn.append(int(c[h].sum()))
            af, nn = np.array(af), np.array(nn)
            ex = af.sum() / nn.sum() * nn
            print(f"   {af.tolist()} of {nn[0]} shots each; chi^2 against one rate {float(((af - ex) ** 2 / ex).sum()):.1f} "
                  f"on {len(af) - 1}")
        print()


if __name__ == "__main__":
    main()
