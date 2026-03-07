"""
star_topology_v3.py

Extended star-topology simulator.
Generalizes star_topology_v2.py from 3 to N+1 qubits.

What this adds over v2:
- arbitrary N+1 qubits: S + N observers
- Bell_SA ⊗ |+>^(N-1) initial state
- equal S<->observer couplings for observers 2..N
- optional direct A<->B coupling J_AB in the 3-qubit case
- threshold sweeps for AB crossing of 1/4 in C·Psi
- shadow-effect analysis after projective Z measurement on A
- analytic threshold fitting helpers

Notes
-----
- Qubit 0 is always S
- Qubit 1 is always A
- Qubit 2 is the primary B used for AB metrics
- For N-observer runs, qubits 2..N are fresh observers with equal coupling to S
- Pair metrics:
    C·Psi = concurrence(pair) * normalized_l1_coherence(pair)
    R     = concurrence(pair) * normalized_l1_coherence(pair)^2

See: experiments/STAR_TOPOLOGY_OBSERVERS.md
"""

from __future__ import annotations

import math
import json
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Sequence, Tuple, Any

import numpy as np


# ---------------------------------------------------------------------
# Basic operators
# ---------------------------------------------------------------------

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
P0 = np.array([[1, 0], [0, 0]], dtype=complex)
P1 = np.array([[0, 0], [0, 1]], dtype=complex)


# ---------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------

def kron_all(ops: Sequence[np.ndarray]) -> np.ndarray:
    """Kronecker product of a sequence of operators."""
    out = np.array([[1.0 + 0.0j]])
    for op in ops:
        out = np.kron(out, op)
    return out


def op_on_qubit(op: np.ndarray, qubit: int, n_qubits: int) -> np.ndarray:
    """Place a single-qubit operator on the chosen qubit."""
    ops = [I2] * n_qubits
    ops[qubit] = op
    return kron_all(ops)


def two_qubit_heisenberg_term(q1: int, q2: int, n_qubits: int, J: float) -> np.ndarray:
    """J * (σxσx + σyσy + σzσz) between q1 and q2."""
    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for sigma in (sx, sy, sz):
        ops = [I2] * n_qubits
        ops[q1] = sigma
        ops[q2] = sigma
        H += J * kron_all(ops)
    return H


def basis_state(bits: Sequence[int]) -> np.ndarray:
    """Computational basis ket |bits>."""
    n = len(bits)
    idx = 0
    for b in bits:
        idx = (idx << 1) | int(b)
    ket = np.zeros(2**n, dtype=complex)
    ket[idx] = 1.0
    return ket


def plus_state() -> np.ndarray:
    """|+> = (|0> + |1>)/sqrt(2)."""
    return np.array([1.0, 1.0], dtype=complex) / math.sqrt(2.0)


def bell_phi_plus() -> np.ndarray:
    """|Φ+> = (|00> + |11>)/sqrt(2)."""
    return np.array([1.0, 0.0, 0.0, 1.0], dtype=complex) / math.sqrt(2.0)


def density_from_statevector(psi: np.ndarray) -> np.ndarray:
    """ρ = |ψ><ψ|."""
    return np.outer(psi, psi.conj())


def round4(x: float) -> float:
    return round(float(x), 4)


def round6(x: float) -> float:
    return round(float(x), 6)


# ---------------------------------------------------------------------
# Partial trace and metrics
# ---------------------------------------------------------------------

def partial_trace_keep(rho: np.ndarray, keep: Sequence[int], n_qubits: int) -> np.ndarray:
    """
    Partial trace keeping the qubits listed in 'keep'.

    Returns the reduced density matrix in the qubit order given by 'keep'.
    """
    keep = list(keep)
    trace_out = [q for q in range(n_qubits) if q not in keep]

    dims = [2] * n_qubits
    reshaped = rho.reshape(dims + dims)

    current_n = n_qubits
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1

    d_keep = 2 ** len(keep)
    return reshaped.reshape((d_keep, d_keep))


def l1_coherence(rho: np.ndarray) -> float:
    d = rho.shape[0]
    total = 0.0
    for i in range(d):
        for j in range(d):
            if i != j:
                total += abs(rho[i, j])
    return float(total)


def psi_norm(rho: np.ndarray) -> float:
    d = rho.shape[0]
    if d <= 1:
        return 0.0
    return l1_coherence(rho) / (d - 1)


def purity(rho: np.ndarray) -> float:
    return float(np.trace(rho @ rho).real)


def concurrence_two_qubit(rho4: np.ndarray) -> float:
    """Wootters concurrence for a 4x4 two-qubit density matrix."""
    sy2 = np.kron(sy, sy)
    rho_tilde = sy2 @ rho4.conj() @ sy2
    R = rho4 @ rho_tilde
    evals = np.linalg.eigvals(R)
    evals = np.real(np.sqrt(np.maximum(evals, 0.0)))
    evals = np.sort(evals)[::-1]
    c = evals[0] - evals[1] - evals[2] - evals[3]
    return float(max(0.0, c))


def pair_metrics(rho_pair: np.ndarray) -> Dict[str, float]:
    p = psi_norm(rho_pair)
    c = concurrence_two_qubit(rho_pair)
    cpsi = c * p
    r = c * (p ** 2)
    return {
        "concurrence": c,
        "psi": p,
        "cpsi": cpsi,
        "R": r,
        "purity": purity(rho_pair),
    }


# ---------------------------------------------------------------------
# Dynamics
# ---------------------------------------------------------------------

def star_hamiltonian_n(
    n_observers: int,
    J_SA: float = 1.0,
    J_SB: float = 1.0,
    J_AB: float = 0.0,
) -> np.ndarray:
    """
    Build the N-observer star topology Hamiltonian.

    Total qubits = 1 + n_observers
    qubit 0: S
    qubit 1: A
    qubits 2..n_observers: B_i

    H = J_SA * (S·A) + J_SB * sum_i (S·B_i)
    plus optional J_AB * (A·B) for the 3-qubit case or whenever qubit 2 exists.
    """
    n_qubits = 1 + n_observers
    if n_qubits < 3:
        raise ValueError("Need at least S + A + B, i.e. n_observers >= 2.")

    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)

    # S <-> A
    H += two_qubit_heisenberg_term(0, 1, n_qubits, J_SA)

    # S <-> B_i for i >= 2
    for q in range(2, n_qubits):
        H += two_qubit_heisenberg_term(0, q, n_qubits, J_SB)

    # Optional direct A <-> B (use the first B, i.e. qubit 2)
    if abs(J_AB) > 0.0 and n_qubits >= 3:
        H += two_qubit_heisenberg_term(1, 2, n_qubits, J_AB)

    return H


def dephasing_ops_n(gammas: Sequence[float]) -> List[np.ndarray]:
    """
    Local sigma_z dephasing Lindblad operators for each qubit.

    Convention matches star_topology_v2.py:
        L_q = sqrt(gamma_q) * sigma_z(q)
    """
    n_qubits = len(gammas)
    ops: List[np.ndarray] = []
    for q, gamma in enumerate(gammas):
        if gamma > 0:
            ops.append(math.sqrt(gamma) * op_on_qubit(sz, q, n_qubits))
    return ops


def lindblad_rhs(rho: np.ndarray, H: np.ndarray, L_ops: Sequence[np.ndarray]) -> np.ndarray:
    drho = -1j * (H @ rho - rho @ H)
    for L in L_ops:
        Ld = L.conj().T
        drho += L @ rho @ Ld - 0.5 * (Ld @ L @ rho + rho @ Ld @ L)
    return drho


def rk4_step(rho: np.ndarray, H: np.ndarray, L_ops: Sequence[np.ndarray], dt: float) -> np.ndarray:
    k1 = lindblad_rhs(rho, H, L_ops)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, L_ops)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, L_ops)
    k4 = lindblad_rhs(rho + dt * k3, H, L_ops)

    rho_new = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    # Numerical hygiene
    rho_new = 0.5 * (rho_new + rho_new.conj().T)
    tr = np.trace(rho_new).real
    if abs(tr) < 1e-14:
        raise FloatingPointError("Trace collapsed to ~0 during RK4 step.")
    rho_new /= tr
    return rho_new


# ---------------------------------------------------------------------
# State preparation and measurement
# ---------------------------------------------------------------------

def bell_sa_plus_rest(n_observers: int) -> np.ndarray:
    """
    Initial state:
        Bell_SA ⊗ |+>^(N-1)
    where total qubits are:
        S, A, B1, B2, ..., B_{N-1}
    """
    if n_observers < 2:
        raise ValueError("Need at least 2 observers: A and B.")

    psi = bell_phi_plus()
    for _ in range(n_observers - 1):
        psi = np.kron(psi, plus_state())
    return density_from_statevector(psi)


def projective_measure_z_on_qubit(rho: np.ndarray, qubit: int, n_qubits: int) -> np.ndarray:
    """Unread projective Z measurement on the chosen qubit."""
    P0_full = op_on_qubit(P0, qubit, n_qubits)
    P1_full = op_on_qubit(P1, qubit, n_qubits)
    measured = P0_full @ rho @ P0_full + P1_full @ rho @ P1_full
    tr = np.trace(measured).real
    if tr > 1e-14:
        measured /= tr
    measured = 0.5 * (measured + measured.conj().T)
    return measured


# ---------------------------------------------------------------------
# Simulation record
# ---------------------------------------------------------------------

@dataclass
class PairRecord:
    concurrence: List[float]
    psi: List[float]
    cpsi: List[float]
    R: List[float]
    purity: List[float]


@dataclass
class SimulationRecord:
    metadata: Dict[str, Any]
    t: List[float]
    purity_full: List[float]
    pairs: Dict[str, PairRecord]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "metadata": self.metadata,
            "t": self.t,
            "purity_full": self.purity_full,
            "pairs": {k: asdict(v) for k, v in self.pairs.items()},
        }


def _empty_pair_record() -> PairRecord:
    return PairRecord(concurrence=[], psi=[], cpsi=[], R=[], purity=[])


def run_star_topology(
    n_observers: int,
    J_SA: float = 1.0,
    J_SB: float = 1.0,
    J_AB: float = 0.0,
    gamma: float = 0.05,
    gammas: Optional[Sequence[float]] = None,
    dt: float = 0.005,
    t_max: float = 5.0,
    sample_every: int = 20,
    measure_a_at: Optional[float] = None,
) -> SimulationRecord:
    """
    Main simulator.

    Parameters
    ----------
    n_observers:
        Number of observers including A.
        Example:
            n_observers = 2 means total qubits = 3 => S + A + B
            n_observers = 5 means total qubits = 6 => S + A + B1 + B2 + B3 + B4
    gamma:
        Used for all qubits if 'gammas' is not provided.
    gammas:
        Optional per-qubit dephasing rates of length 1 + n_observers.
    sample_every:
        Record every 'sample_every' RK4 steps.
    measure_a_at:
        If given, apply unread projective Z measurement on qubit A (qubit 1)
        at the first time t >= measure_a_at.
    """
    n_qubits = 1 + n_observers
    if n_qubits < 3:
        raise ValueError("Need at least 3 total qubits (S+A+B).")

    if gammas is None:
        gammas = [gamma] * n_qubits
    else:
        gammas = list(gammas)
        if len(gammas) != n_qubits:
            raise ValueError(f"Expected {n_qubits} gamma values, got {len(gammas)}.")

    H = star_hamiltonian_n(n_observers=n_observers, J_SA=J_SA, J_SB=J_SB, J_AB=J_AB)
    L_ops = dephasing_ops_n(gammas)
    rho = bell_sa_plus_rest(n_observers)

    steps = int(round(t_max / dt))
    measured = False

    pair_keys: Dict[str, Tuple[int, int]] = {
        "SA": (0, 1),
        "SB": (0, 2),
        "AB": (1, 2),
    }

    # For N > 2, also record A-B_i and S-B_i for all fresh observers
    for q in range(3, n_qubits):
        pair_keys[f"S{q}"] = (0, q)
        pair_keys[f"A{q}"] = (1, q)

    pairs = {name: _empty_pair_record() for name in pair_keys}

    rec_t: List[float] = []
    rec_purity_full: List[float] = []

    for step in range(steps + 1):
        t = step * dt

        if measure_a_at is not None and (not measured) and t >= measure_a_at:
            rho = projective_measure_z_on_qubit(rho, qubit=1, n_qubits=n_qubits)
            measured = True

        if step % sample_every == 0:
            rec_t.append(round6(t))
            rec_purity_full.append(round6(purity(rho)))

            for name, (q1, q2) in pair_keys.items():
                rho_pair = partial_trace_keep(rho, keep=[q1, q2], n_qubits=n_qubits)
                m = pair_metrics(rho_pair)
                pairs[name].concurrence.append(round6(m["concurrence"]))
                pairs[name].psi.append(round6(m["psi"]))
                pairs[name].cpsi.append(round6(m["cpsi"]))
                pairs[name].R.append(round6(m["R"]))
                pairs[name].purity.append(round6(m["purity"]))

        if step < steps:
            rho = rk4_step(rho, H, L_ops, dt)

    metadata = {
        "n_observers": n_observers,
        "n_qubits": n_qubits,
        "J_SA": J_SA,
        "J_SB": J_SB,
        "J_AB": J_AB,
        "gamma": gamma if gammas is None else None,
        "gammas": list(gammas),
        "dt": dt,
        "t_max": t_max,
        "sample_every": sample_every,
        "measure_a_at": measure_a_at,
        "initial_state": "Bell_SA ⊗ |+>^(N-1)",
    }

    return SimulationRecord(metadata=metadata, t=rec_t, purity_full=rec_purity_full, pairs=pairs)


# ---------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------

def first_upward_crossing(times: Sequence[float], vals: Sequence[float], threshold: float = 0.25) -> Optional[float]:
    """Linear interpolation for first upward crossing of threshold."""
    for i in range(1, len(vals)):
        v0, v1 = vals[i - 1], vals[i]
        if v0 < threshold <= v1 and abs(v1 - v0) > 1e-14:
            t0, t1 = times[i - 1], times[i]
            frac = (threshold - v0) / (v1 - v0)
            return float(t0 + frac * (t1 - t0))
    return None


def max_after_time(times: Sequence[float], vals: Sequence[float], t_min: float) -> float:
    arr = [v for t, v in zip(times, vals) if t >= t_min]
    if not arr:
        return float("nan")
    return float(max(arr))


def value_at_time_linear(times: Sequence[float], vals: Sequence[float], target_t: float) -> float:
    if target_t <= times[0]:
        return float(vals[0])
    if target_t >= times[-1]:
        return float(vals[-1])
    for i in range(1, len(times)):
        if times[i] >= target_t:
            t0, t1 = times[i - 1], times[i]
            v0, v1 = vals[i - 1], vals[i]
            if abs(t1 - t0) < 1e-14:
                return float(v1)
            frac = (target_t - t0) / (t1 - t0)
            return float(v0 + frac * (v1 - v0))
    return float(vals[-1])


def peak_pair_R(rec: SimulationRecord, pair_key: str = "AB") -> float:
    return float(max(rec.pairs[pair_key].R))


def peak_pair_cpsi(rec: SimulationRecord, pair_key: str = "AB") -> float:
    return float(max(rec.pairs[pair_key].cpsi))


def compute_shadow_suppression(
    n_observers: int,
    J_SA: float,
    J_SB: float,
    gamma: float,
    J_AB: float = 0.0,
    measure_at: float = 1.0,
    compare_pair: str = "SB",
    dt: float = 0.005,
    t_max: float = 5.0,
) -> Dict[str, float]:
    """
    Compare peak post-measurement R for a pair with and without projective measurement on A.
    """
    base = run_star_topology(
        n_observers=n_observers,
        J_SA=J_SA,
        J_SB=J_SB,
        J_AB=J_AB,
        gamma=gamma,
        dt=dt,
        t_max=t_max,
        measure_a_at=None,
    )
    meas = run_star_topology(
        n_observers=n_observers,
        J_SA=J_SA,
        J_SB=J_SB,
        J_AB=J_AB,
        gamma=gamma,
        dt=dt,
        t_max=t_max,
        measure_a_at=measure_at,
    )

    peak_base = max_after_time(base.t, base.pairs[compare_pair].R, measure_at)
    peak_meas = max_after_time(meas.t, meas.pairs[compare_pair].R, measure_at)

    if abs(peak_base) < 1e-14:
        rel = float("nan")
    else:
        rel = 100.0 * (peak_base - peak_meas) / peak_base

    return {
        "peak_without_measurement": round6(peak_base),
        "peak_with_measurement": round6(peak_meas),
        "relative_suppression_percent": round6(rel),
    }


# ---------------------------------------------------------------------
# Threshold sweeps
# ---------------------------------------------------------------------

def _scan_range(start: float, stop: float, step: float) -> List[float]:
    values = []
    x = start
    while x <= stop + 1e-12:
        values.append(round(x, 10))
        x += step
    return values


def find_jsb_threshold(
    n_observers: int,
    J_SA: float = 1.0,
    gamma: float = 0.05,
    J_AB: float = 0.0,
    cpsi_threshold: float = 0.25,
    J_start: float = 1.0,
    J_stop: float = 5.0,
    J_step: float = 0.1,
    fine_half_width: float = 0.05,
    fine_step: float = 0.0025,
    dt: float = 0.005,
    t_max: float = 5.0,
) -> Dict[str, Any]:
    """
    Find the minimal J_SB for which AB crosses cpsi_threshold.
    """
    coarse_hits: List[Tuple[float, float, Optional[float]]] = []

    for J in _scan_range(J_start, J_stop, J_step):
        rec = run_star_topology(
            n_observers=n_observers,
            J_SA=J_SA,
            J_SB=J,
            J_AB=J_AB,
            gamma=gamma,
            dt=dt,
            t_max=t_max,
        )
        peak_cpsi = peak_pair_cpsi(rec, "AB")
        tc = first_upward_crossing(rec.t, rec.pairs["AB"].cpsi, cpsi_threshold)
        coarse_hits.append((J, peak_cpsi, tc))

    passing = [item for item in coarse_hits if item[1] >= cpsi_threshold]
    if not passing:
        return {
            "crosses": False,
            "threshold_J_SB": None,
            "coarse_scan": coarse_hits,
            "fine_scan": [],
            "message": f"No crossing found up to J_SB = {J_stop:.4f}",
        }

    J0 = passing[0][0]
    fine_start = max(J_start, J0 - fine_half_width)
    fine_stop = min(J_stop, J0 + fine_half_width)

    fine_hits: List[Tuple[float, float, Optional[float]]] = []
    for J in _scan_range(fine_start, fine_stop, fine_step):
        rec = run_star_topology(
            n_observers=n_observers,
            J_SA=J_SA,
            J_SB=J,
            J_AB=J_AB,
            gamma=gamma,
            dt=dt,
            t_max=t_max,
        )
        peak_cpsi = peak_pair_cpsi(rec, "AB")
        tc = first_upward_crossing(rec.t, rec.pairs["AB"].cpsi, cpsi_threshold)
        fine_hits.append((J, peak_cpsi, tc))

    fine_passing = [item for item in fine_hits if item[1] >= cpsi_threshold]
    best = fine_passing[0] if fine_passing else passing[0]

    return {
        "crosses": True,
        "threshold_J_SB": round6(best[0]),
        "crossing_time": None if best[2] is None else round6(best[2]),
        "peak_AB_cpsi": round6(best[1]),
        "coarse_scan": coarse_hits,
        "fine_scan": fine_hits,
    }


def sweep_n_observers(
    observer_counts: Sequence[int] = (2, 3, 4, 5),
    J_SA: float = 1.0,
    gamma: float = 0.05,
    J_start: float = 1.0,
    J_stop: float = 5.0,
    J_step: float = 0.1,
    fine_half_width: float = 0.05,
    fine_step: float = 0.0025,
    dt: float = 0.005,
    t_max: float = 5.0,
) -> List[Dict[str, Any]]:
    rows = []
    for n_obs in observer_counts:
        result = find_jsb_threshold(
            n_observers=n_obs,
            J_SA=J_SA,
            gamma=gamma,
            J_start=J_start,
            J_stop=J_stop,
            J_step=J_step,
            fine_half_width=fine_half_width,
            fine_step=fine_step,
            dt=dt,
            t_max=t_max,
        )

        if result["crosses"]:
            J_use = result["threshold_J_SB"]
        else:
            J_use = J_stop

        rec = run_star_topology(
            n_observers=n_obs,
            J_SA=J_SA,
            J_SB=J_use,
            gamma=gamma,
            dt=dt,
            t_max=t_max,
        )
        shadow = compute_shadow_suppression(
            n_observers=n_obs,
            J_SA=J_SA,
            J_SB=J_use,
            gamma=gamma,
            measure_at=1.0,
            compare_pair="SB",
            dt=dt,
            t_max=t_max,
        )

        rows.append({
            "n_observers": n_obs,
            "n_qubits": n_obs + 1,
            "crosses": result["crosses"],
            "J_SB_threshold": result["threshold_J_SB"],
            "crossing_time": result.get("crossing_time"),
            "J_SB_used_for_metrics": round6(J_use),
            "peak_R_AB": round6(peak_pair_R(rec, "AB")),
            "peak_R_SB": round6(peak_pair_R(rec, "SB")),
            "shadow_peak_R_SB_without": shadow["peak_without_measurement"],
            "shadow_peak_R_SB_with": shadow["peak_with_measurement"],
            "shadow_suppression_percent": shadow["relative_suppression_percent"],
        })
    return rows


def sweep_direct_coupling(
    J_AB_values: Sequence[float],
    J_SA: float = 1.0,
    J_SB: float = 1.466,
    gamma: float = 0.05,
    dt: float = 0.005,
    t_max: float = 5.0,
) -> List[Dict[str, Any]]:
    rows = []
    for J_AB in J_AB_values:
        rec = run_star_topology(
            n_observers=2,
            J_SA=J_SA,
            J_SB=J_SB,
            J_AB=J_AB,
            gamma=gamma,
            dt=dt,
            t_max=t_max,
        )
        shadow = compute_shadow_suppression(
            n_observers=2,
            J_SA=J_SA,
            J_SB=J_SB,
            J_AB=J_AB,
            gamma=gamma,
            measure_at=1.0,
            compare_pair="SB",
            dt=dt,
            t_max=t_max,
        )
        rows.append({
            "J_AB": round6(J_AB),
            "max_CPsi_AB": round6(max(rec.pairs["AB"].cpsi)),
            "max_R_AB": round6(max(rec.pairs["AB"].R)),
            "max_R_SB": round6(max(rec.pairs["SB"].R)),
            "first_AB_crossing": None if first_upward_crossing(rec.t, rec.pairs["AB"].cpsi, 0.25) is None
                                 else round6(first_upward_crossing(rec.t, rec.pairs["AB"].cpsi, 0.25)),
            "shadow_peak_R_SB_without": shadow["peak_without_measurement"],
            "shadow_peak_R_SB_with": shadow["peak_with_measurement"],
            "shadow_relative_change_percent": shadow["relative_suppression_percent"],
        })
    return rows


# ---------------------------------------------------------------------
# Simple curve fitting helpers
# ---------------------------------------------------------------------

def _r2_score(y: np.ndarray, yhat: np.ndarray) -> float:
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot == 0:
        return 1.0
    return 1.0 - ss_res / ss_tot


def fit_linear(gamma_vals: Sequence[float], threshold_vals: Sequence[float]) -> Dict[str, Any]:
    x = np.asarray(gamma_vals, dtype=float)
    y = np.asarray(threshold_vals, dtype=float)

    A = np.column_stack([x, np.ones_like(x)])
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    a, b = coeffs
    yhat = a * x + b
    return {
        "model": "linear",
        "formula": "J_th(gamma) = a*gamma + b",
        "params": {"a": float(a), "b": float(b)},
        "r2": _r2_score(y, yhat),
        "yhat": yhat.tolist(),
    }


def fit_power_plus_constant_grid(
    gamma_vals: Sequence[float],
    threshold_vals: Sequence[float],
    b_grid: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Fit J = a * gamma^b + c using a 1D grid over b.
    For each b, solve linear least squares in (a, c).
    """
    x = np.asarray(gamma_vals, dtype=float)
    y = np.asarray(threshold_vals, dtype=float)

    if b_grid is None:
        b_grid = np.linspace(0.2, 3.0, 4000)

    best = None
    for b in b_grid:
        xb = x ** b
        A = np.column_stack([xb, np.ones_like(x)])
        coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        a, c = coeffs
        yhat = a * xb + c
        r2 = _r2_score(y, yhat)
        cand = {
            "model": "power_plus_constant",
            "formula": "J_th(gamma) = a*gamma^b + c",
            "params": {"a": float(a), "b": float(b), "c": float(c)},
            "r2": float(r2),
            "yhat": yhat.tolist(),
        }
        if best is None or cand["r2"] > best["r2"]:
            best = cand
    return best


def fit_exponential_plus_constant_grid(
    gamma_vals: Sequence[float],
    threshold_vals: Sequence[float],
    b_grid: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Fit J = a * exp(b*gamma) + c using a 1D grid over b.
    For each b, solve linear least squares in (a, c).
    """
    x = np.asarray(gamma_vals, dtype=float)
    y = np.asarray(threshold_vals, dtype=float)

    if b_grid is None:
        b_grid = np.linspace(0.1, 10.0, 4000)

    best = None
    for b in b_grid:
        eb = np.exp(b * x)
        A = np.column_stack([eb, np.ones_like(x)])
        coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        a, c = coeffs
        yhat = a * eb + c
        r2 = _r2_score(y, yhat)
        cand = {
            "model": "exp_plus_constant",
            "formula": "J_th(gamma) = a*exp(b*gamma) + c",
            "params": {"a": float(a), "b": float(b), "c": float(c)},
            "r2": float(r2),
            "yhat": yhat.tolist(),
        }
        if best is None or cand["r2"] > best["r2"]:
            best = cand
    return best


def fit_rational_fixed_pole(
    gamma_vals: Sequence[float],
    threshold_vals: Sequence[float],
    pole: float = 0.2,
) -> Dict[str, Any]:
    """
    Fit J = a / (pole - gamma) + c by linear least squares.
    """
    x = np.asarray(gamma_vals, dtype=float)
    y = np.asarray(threshold_vals, dtype=float)

    basis = 1.0 / (pole - x)
    A = np.column_stack([basis, np.ones_like(x)])
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    a, c = coeffs
    yhat = a * basis + c
    return {
        "model": "rational_fixed_pole",
        "formula": "J_th(gamma) = a/(pole-gamma) + c",
        "params": {"a": float(a), "pole": float(pole), "c": float(c)},
        "r2": _r2_score(y, yhat),
        "yhat": yhat.tolist(),
    }


def fit_threshold_models(gamma_vals: Sequence[float], threshold_vals: Sequence[float]) -> List[Dict[str, Any]]:
    fits = [
        fit_linear(gamma_vals, threshold_vals),
        fit_power_plus_constant_grid(gamma_vals, threshold_vals),
        fit_exponential_plus_constant_grid(gamma_vals, threshold_vals),
        fit_rational_fixed_pole(gamma_vals, threshold_vals, pole=0.2),
    ]
    return sorted(fits, key=lambda x: x["r2"], reverse=True)


# ---------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------

def save_json(obj: Any, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


# ---------------------------------------------------------------------
# Example CLI-style usage
# ---------------------------------------------------------------------

def _print_table(rows: List[Dict[str, Any]], columns: Sequence[str]) -> None:
    widths = {}
    for col in columns:
        widths[col] = max(len(col), *(len(str(r.get(col, ""))) for r in rows)) if rows else len(col)

    header = " | ".join(col.ljust(widths[col]) for col in columns)
    sep = "-+-".join("-" * widths[col] for col in columns)
    print(header)
    print(sep)

    for row in rows:
        print(" | ".join(str(row.get(col, "")).ljust(widths[col]) for col in columns))


if __name__ == "__main__":
    print("Running example sweeps for extended star topology...\n")

    # Task 1 style sweep
    rows_n = sweep_n_observers(observer_counts=(2, 3, 4, 5))
    print("N-observer sweep:")
    _print_table(
        rows_n,
        [
            "n_observers",
            "n_qubits",
            "crosses",
            "J_SB_threshold",
            "crossing_time",
            "peak_R_AB",
            "peak_R_SB",
            "shadow_suppression_percent",
        ],
    )

    print("\nTask-2 style fits:")
    gamma_vals = [0.001, 0.010, 0.050, 0.100, 0.150]
    threshold_vals = [1.179, 1.237, 1.465, 1.776, 2.145]
    fits = fit_threshold_models(gamma_vals, threshold_vals)
    for fit in fits:
        print(f"{fit['model']}: R^2 = {fit['r2']:.6f}, params = {fit['params']}")

    print("\nTask-3 style direct coupling sweep:")
    rows_jab = sweep_direct_coupling(J_AB_values=[0.0, 0.1, 0.3, 0.5, 1.0])
    _print_table(
        rows_jab,
        [
            "J_AB",
            "max_CPsi_AB",
            "max_R_AB",
            "max_R_SB",
            "first_AB_crossing",
            "shadow_relative_change_percent",
        ],
    )

    save_json(rows_n, "star_topology_n_observers_sweep.json")
    save_json(rows_jab, "star_topology_direct_coupling_sweep.json")
    print("\nWrote example JSON outputs:")
    print("- star_topology_n_observers_sweep.json")
    print("- star_topology_direct_coupling_sweep.json")
