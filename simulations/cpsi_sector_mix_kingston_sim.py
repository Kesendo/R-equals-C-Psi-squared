#!/usr/bin/env python3
"""
Kingston-pipeline simulation of the sector-mix test.

Reproduces the proposed hardware run end-to-end in Qiskit Aer:
  1. Choose a 3-qubit linear chain on Kingston (from real calibration).
  2. Build a Kingston-noise model (per-qubit T1/T2).
  3. Prepare each of four test states via qc.initialize():
       |+⟩^3, GHZ_3, W_3, OPT = 0.6125 |GHZ⟩ + 0.7905 |W⟩
  4. Sweep delay via id-chunks with the thermal_relaxation_error noise model.
  5. save_density_matrix → partial_trace to pair (q_0, q_1) → CΨ.
  6. Compare to analytical Lindblad evolution (no gate errors).

Purpose: verify the hardware-visible signal before burning QPU time.
If the OPT state still crosses ¼ well above the state-prep noise floor,
the hardware run will work.

Date: 2026-04-16
"""

from __future__ import annotations

import csv
import itertools
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RCPSI_ROOT = Path(__file__).parent.parent
CALIB_DIR = RCPSI_ROOT / "ClaudeTasks" / "IBM_R2_calibrations"
RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)
OUT = RESULTS / "cpsi_sector_mix_kingston_sim.txt"
CHUNK_US = 1.0
SHOTS = 2048  # matches the hardware plan


# ── state vectors (big-endian to match qiskit initialize semantics) ──
def sv_ghz3() -> np.ndarray:
    v = np.zeros(8, dtype=complex)
    v[0] = 1 / np.sqrt(2)
    v[7] = 1 / np.sqrt(2)
    return v


def sv_w3() -> np.ndarray:
    v = np.zeros(8, dtype=complex)
    for i in range(3):
        v[1 << i] = 1 / np.sqrt(3)
    return v


def sv_opt() -> np.ndarray:
    alpha = 0.6125
    beta = math.sqrt(1 - alpha * alpha)  # ≈ 0.7905
    return alpha * sv_ghz3() + beta * sv_w3()


def sv_plus3() -> np.ndarray:
    return np.ones(8, dtype=complex) / math.sqrt(8)


# ── calibration parsing (subset of run_bonding_mode.py) ──
def _safe_float(s, default=0.0):
    if s is None:
        return default
    s = str(s).strip()
    if not s:
        return default
    try:
        return float(s)
    except ValueError:
        return default


def _parse_neighbors(s: str):
    out = {}
    if not s:
        return out
    for item in s.split(";"):
        m = re.match(r"(\d+)\s*:\s*([\d.eE+\-]+)", item.strip())
        if m:
            out[int(m.group(1))] = float(m.group(2))
    return out


def load_calibration():
    csvs = sorted(CALIB_DIR.glob("ibm_kingston_calibrations_*.csv"))
    if not csvs:
        raise FileNotFoundError("No Kingston calibration CSV")
    latest = csvs[-1]
    qubits = {}
    edges = {}
    with open(latest, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            q = int(row["Qubit"])
            qubits[q] = {
                "T1": _safe_float(row["T1 (us)"]) * 1e-6,
                "T2": _safe_float(row["T2 (us)"]) * 1e-6,
                "sq_len": _safe_float(row["Single-qubit gate length (ns)"]) * 1e-9,
            }
            for nb, err in _parse_neighbors(row.get("CZ error", "")).items():
                edges[(min(q, nb), max(q, nb))] = {"cz_err": err}
            for nb, length in _parse_neighbors(row.get("Gate length (ns)", "")).items():
                key = (min(q, nb), max(q, nb))
                if key not in edges:
                    edges[key] = {}
                edges[key]["gate_len"] = length * 1e-9
    # adjacency
    adj = {q: [] for q in qubits}
    for (a, b) in edges:
        adj[a].append(b)
        adj[b].append(a)
    return {"src": str(latest), "qubits": qubits, "edges": edges, "adj": adj}


def best_3_chain(cal):
    """Pick a length-3 linear path maximizing min(T2) × product(1-cz_err)."""
    qubits = [q for q, d in cal["qubits"].items() if d["T2"] > 50e-6]
    best = None
    for q0 in qubits:
        for q1 in cal["adj"].get(q0, []):
            if q1 not in cal["qubits"]:
                continue
            for q2 in cal["adj"].get(q1, []):
                if q2 == q0 or q2 not in cal["qubits"]:
                    continue
                path = [q0, q1, q2]
                t2s = [cal["qubits"][q]["T2"] for q in path]
                t1s = [cal["qubits"][q]["T1"] for q in path]
                e01 = cal["edges"].get((min(q0, q1), max(q0, q1)), {})
                e12 = cal["edges"].get((min(q1, q2), max(q1, q2)), {})
                cz0 = e01.get("cz_err", 0.01)
                cz1 = e12.get("cz_err", 0.01)
                score = min(t2s) * (1 - cz0) * (1 - cz1)
                if best is None or score > best["score"]:
                    best = {"path": path, "t2s_us": [t * 1e6 for t in t2s],
                            "t1s_us": [t * 1e6 for t in t1s],
                            "cz_errs": [cz0, cz1],
                            "gate_lens": [e01.get("gate_len", 68e-9),
                                          e12.get("gate_len", 68e-9)],
                            "score": score}
    return best


# ── analytical Lindblad evolution ──
I2 = np.eye(2, dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site_z(site: int, n: int) -> np.ndarray:
    out = I2 if site != 0 else Z
    for s in range(1, n):
        out = np.kron(out, Z if s == site else I2)
    return out


def lindblad_op(gammas, n: int) -> np.ndarray:
    d = 2 ** n
    Id = np.eye(d, dtype=complex)
    L = np.zeros((d * d, d * d), dtype=complex)
    for site, g in enumerate(gammas):
        if g <= 0:
            continue
        Lk = math.sqrt(g) * site_z(site, n)
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk) - 0.5 * np.kron(Id, LdL) - 0.5 * np.kron(LdL.T, Id))
    return L


def partial_trace_01(rho, n):
    """Trace all qubits except 0 and 1 (big-endian: MSB = q0)."""
    rho_r = rho.reshape([2] * (2 * n))
    current = list(range(n))
    traced = rho_r
    while len(current) > 2:
        rem = next(q for q in current if q not in (0, 1))
        idx = current.index(rem)
        row_ax = len(current) - 1 - idx
        col_ax = 2 * len(current) - 1 - idx
        traced = np.trace(traced, axis1=row_ax, axis2=col_ax)
        current.pop(idx)
    return traced.reshape(4, 4)


def cpsi_pair(rho2):
    C = float(np.real(np.trace(rho2 @ rho2)))
    diag = np.diag(np.diag(rho2))
    L1 = float(np.sum(np.abs(rho2 - diag)))
    return C, L1 / 3, C * L1 / 3


# ── Aer pipeline ──
def run_aer_pipeline(state_vec, delays_us, chain_info, include_gate_noise=True):
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, thermal_relaxation_error
    from qiskit.quantum_info import partial_trace, DensityMatrix

    n = 3
    T1s = [t * 1e-6 for t in chain_info["t1s_us"]]
    T2s = [min(t * 1e-6, 2.0 * T1) for t, T1 in zip(chain_info["t2s_us"], T1s)]

    nm = NoiseModel()
    if include_gate_noise:
        # single-qubit gate noise (H, RY, RZ, SX, X on each virtual qubit)
        for vq in range(n):
            err = thermal_relaxation_error(T1s[vq], T2s[vq], 32e-9)
            nm.add_quantum_error(err, ["sx", "rz", "x", "h", "ry", "u", "u1", "u2", "u3"], [vq])
            idle = thermal_relaxation_error(T1s[vq], T2s[vq], CHUNK_US * 1e-6)
            nm.add_quantum_error(idle, ["id"], [vq])
        # 2-qubit cx / cz on (0,1) and (1,2)
        for (a, b), glen in zip([(0, 1), (1, 2)], chain_info["gate_lens"]):
            err_2q = thermal_relaxation_error(T1s[a], T2s[a], glen).tensor(
                thermal_relaxation_error(T1s[b], T2s[b], glen))
            for name in ("cz", "cx", "rzz", "swap", "ecr"):
                nm.add_quantum_error(err_2q, [name], [a, b])
                nm.add_quantum_error(err_2q, [name], [b, a])
    else:
        # Pure dephasing only, no T1. Pure dephasing time ≈ T2 (worst case).
        for vq in range(n):
            idle = thermal_relaxation_error(
                1e6, T2s[vq], CHUNK_US * 1e-6)  # huge T1 → pure dephasing
            nm.add_quantum_error(idle, ["id"], [vq])

    backend = AerSimulator(noise_model=nm, method="density_matrix")

    def build_circ(delay_us):
        qc = QuantumCircuit(n)
        qc.initialize(state_vec, list(range(n)))
        if delay_us > 0:
            n_chunks = max(1, int(round(delay_us / CHUNK_US)))
            for _ in range(n_chunks):
                for q in range(n):
                    qc.id(q)
        qc.save_density_matrix()
        return qc

    results = []
    for t in delays_us:
        qc = build_circ(float(t))
        qc_t = transpile(qc, backend=backend, optimization_level=0)
        rho = np.asarray(backend.run(qc_t, shots=SHOTS).result().data(0)["density_matrix"])
        # trace to (0, 1)
        rho2 = np.asarray(partial_trace(DensityMatrix(rho), qargs=[2]).data)
        C, P, CP = cpsi_pair(rho2)
        results.append({"t_us": float(t), "C": C, "Psi": P, "CPsi": CP})
    return results


# ── main ──
if __name__ == "__main__":
    with open(OUT, "w", encoding="utf-8") as f:
        def log(msg=""):
            print(msg, flush=True)
            f.write(msg + "\n")

        log("=" * 76)
        log("  SECTOR-MIX TEST: KINGSTON PIPELINE SIMULATION")
        log("=" * 76)

        cal = load_calibration()
        log(f"\n  calibration: {Path(cal['src']).name}")
        chain = best_3_chain(cal)
        log(f"  best 3-qubit chain: {chain['path']}  "
            f"T2=[{', '.join(f'{x:.0f}' for x in chain['t2s_us'])}] μs, "
            f"T1=[{', '.join(f'{x:.0f}' for x in chain['t1s_us'])}] μs")
        log(f"  CZ errors: {['%.4f'%e for e in chain['cz_errs']]}")

        states = {
            "|+>^3 (product)": sv_plus3(),
            "GHZ_3": sv_ghz3(),
            "W_3": sv_w3(),
            "OPT (0.6125·GHZ + 0.7905·W)": sv_opt(),
        }

        delays = np.array([0.0, 3.0, 7.0, 11.0, 15.0, 20.0, 30.0, 50.0])
        log()
        log(f"  delay sweep: {list(delays)} μs")
        log(f"  shots per circuit: {SHOTS}")
        log()

        # Analytical reference (Lindblad, no gate errors)
        gammas = [1.0 / (t * 1e-6) for t in chain["t2s_us"]]  # rad/s → 1/s
        # convert to 1/μs for convenience:
        gammas_per_us = [1.0 / t for t in chain["t2s_us"]]
        L_us = lindblad_op(gammas_per_us, 3)
        log("  Column legend for each state:")
        log("    t        CΨ_analytical   CΨ_gates+T2     CΨ_pure-T2")
        log()

        all_data = {}
        for name, sv in states.items():
            rho0 = np.outer(sv, sv.conj())
            log(f"  ── {name} ──")
            aer_full = run_aer_pipeline(sv, delays, chain, include_gate_noise=True)
            aer_pure = run_aer_pipeline(sv, delays, chain, include_gate_noise=False)
            analyt = []
            for t in delays:
                rho_vec = expm(L_us * float(t)) @ rho0.flatten(order="F")
                rho_t = rho_vec.reshape(8, 8, order="F")
                rho2 = partial_trace_01(rho_t, 3)
                _, _, cp = cpsi_pair(rho2)
                analyt.append(cp)
            for i, t in enumerate(delays):
                mark_analyt = "  *" if analyt[i] < 0.25 and (i == 0 or analyt[i - 1] >= 0.25) else ""
                mark_aer = "  *" if aer_full[i]["CPsi"] < 0.25 and (i == 0 or aer_full[i - 1]["CPsi"] >= 0.25) else ""
                log(f"    t={t:6.1f} μs   analyt={analyt[i]:.4f}{mark_analyt:<4s}   "
                    f"full-noise={aer_full[i]['CPsi']:.4f}{mark_aer:<4s}   "
                    f"pure-T2={aer_pure[i]['CPsi']:.4f}")
            all_data[name] = {"analyt": analyt, "aer_full": aer_full, "aer_pure": aer_pure,
                              "delays_us": delays.tolist()}
            log()

        # Summary: where does the OPT state cross 1/4 in each model?
        log("=" * 76)
        log("  OPT state crossing summary")
        log("=" * 76)
        opt = all_data["OPT (0.6125·GHZ + 0.7905·W)"]
        for label, series in [("analyt", opt["analyt"]),
                               ("aer_full (gates+T2)", [r["CPsi"] for r in opt["aer_full"]]),
                               ("aer_pure-T2", [r["CPsi"] for r in opt["aer_pure"]])]:
            t_cross = None
            for i in range(1, len(delays)):
                if series[i - 1] >= 0.25 > series[i]:
                    t_cross = delays[i - 1] + (0.25 - series[i - 1]) * (
                        delays[i] - delays[i - 1]) / (series[i] - series[i - 1])
                    break
            log(f"  {label:<22s}  CΨ(0) = {series[0]:.4f}   "
                f"t*(CΨ=1/4) = {f'{t_cross:.2f} μs' if t_cross else '(not crossed in range)'}")

        log()
        log("  INTERPRETATION")
        log("  ──────────────")
        log("  • analyt = pure Lindblad with site-dependent Z-dephasing (no state-prep errors)")
        log("  • aer_full = same chain, Kingston gate errors on initialize() + delay")
        log("  • aer_pure-T2 = pure dephasing only, infinite T1 (isolates dephasing contribution)")
        log()
        log("  The three OPT crossings should bracket the expected hardware crossing.")
        log("  If aer_full and analyt agree within ~0.02 on CΨ(0), the state-prep circuit")
        log("  is faithful enough for the hardware run.")
        log()
        log(f"  Output: {OUT}")
