#!/usr/bin/env python3
"""Finale Aer design: 150-qubit bonding:1 receiver state preparation on Marrakesh.

The plan: prepare F67 bonding:1 single-excitation eigenmode |ψ_1⟩ on N qubits
via Givens-rotation cascade, measure Z-populations on every qubit. The
populations should follow F65: |c_j|² = (2/(N+1)) sin²(π(j+1)/(N+1)).

Why bonding:1 and not Trotter evolution:
- State prep alone: ~N-1 Givens rotations × 2 CZ each = ~2(N-1) CZ gates
- No Trotter: zero further evolution, no compounding gate noise
- Measurement: single Z-basis tomography (1 circuit), 4096 shots

At N=150 the gate budget is ~300 CZ on Marrakesh. With CZ error ~0.5%:
fidelity ≈ 0.995^300 ≈ 0.22. Population signal:
  P(j) = fidelity · |c_j|² + (1 − fidelity) · 0.5
  P_peak ≈ 0.22 · 0.013 + 0.78 · 0.5 = 0.39 (signal swamped by noise floor)

So pure bonding:1 at N=150 will be marginal. Better strategy:
  - Compare bonding:1 vs bonding:2 vs bonding:3 (different node patterns)
  - The QUALITATIVE shape difference (1 arch vs 2 arches vs 3 arches) is
    robust against depolarization

This script:
  1. Builds Givens-cascade preparation for bonding:k at multiple N
  2. Validates on Aer for N = 10, 30, 50 (full state vector)
  3. Uses analytical noise model for N = 100, 150 (state vector too large)
  4. Decides whether to commit to N=150 hardware or pull back

Tom & Claude — 2026-04-26 (finale)
"""
import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ──────────────────────────────────────────────────────────────────────
# Bonding-mode amplitudes (F65)
# ──────────────────────────────────────────────────────────────────────

def bonding_k_amplitudes(N: int, k: int) -> np.ndarray:
    """F65: c_j = √(2/(N+1)) sin(πk(j+1)/(N+1)) for j = 0..N-1."""
    norm = math.sqrt(2.0 / (N + 1))
    return np.array([
        norm * math.sin(math.pi * k * (j + 1) / (N + 1))
        for j in range(N)
    ])


def bonding_k_populations(N: int, k: int) -> np.ndarray:
    """|c_j|² for bonding:k: (2/(N+1)) sin²(πk(j+1)/(N+1))."""
    return bonding_k_amplitudes(N, k) ** 2


# ──────────────────────────────────────────────────────────────────────
# Givens-cascade preparation
# ──────────────────────────────────────────────────────────────────────

def givens_angles(amplitudes: np.ndarray) -> list[float]:
    """Compute Givens rotation angles for single-excitation state preparation.

    Starts with X|0⟩ at site 0. After k Givens rotations on (k, k+1),
    the amplitude at site k is c_k. The remaining amplitude r_k cascades
    to higher sites.

    Returns angles θ_0, θ_1, ..., θ_{N-2}.
    """
    N = len(amplitudes)
    c = np.array(amplitudes, dtype=float)
    angles = []
    remaining = 1.0
    for k in range(N - 1):
        if abs(remaining) < 1e-15:
            angles.append(0.0)
            continue
        ratio = c[k] / remaining
        # Clamp to valid arccos range
        ratio = max(-1.0, min(1.0, ratio))
        theta_k = math.acos(ratio)
        angles.append(theta_k)
        # Update remaining: r_{k+1}² = r_k² − c_k²
        rem_sq = remaining ** 2 - c[k] ** 2
        remaining = math.sqrt(max(0.0, rem_sq))
    return angles


def verify_givens_amplitudes(N: int, angles: list[float]) -> np.ndarray:
    """Apply Givens cascade analytically to |1_0⟩ and return resulting amplitudes."""
    amp = np.zeros(N)
    amp[0] = 1.0  # start with |1_0⟩
    for k in range(N - 1):
        # Givens on (k, k+1) with angle theta_k:
        # amp[k] → cos(theta_k) · amp[k]
        # amp[k+1] → sin(theta_k) · amp[k] + amp[k+1] (latter is 0 here)
        c, s = math.cos(angles[k]), math.sin(angles[k])
        new_k = c * amp[k]
        new_k1 = s * amp[k]
        amp[k] = new_k
        amp[k + 1] = new_k1
    return amp


# ──────────────────────────────────────────────────────────────────────
# Aer simulation (small N, full state vector)
# ──────────────────────────────────────────────────────────────────────

def simulate_aer(N: int, k: int, n_shots: int = 4096, with_noise: bool = True):
    """Build Givens cascade as Qiskit circuit, run on Aer, return Z-populations."""
    from qiskit import QuantumCircuit, ClassicalRegister, transpile
    from qiskit.circuit.library import XXPlusYYGate
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, thermal_relaxation_error, ReadoutError

    amps = bonding_k_amplitudes(N, k)
    angles = givens_angles(amps)

    qc = QuantumCircuit(N, name=f"bonding_{k}_N{N}")
    cr = ClassicalRegister(N, name='meas')
    qc.add_register(cr)
    qc.x(0)  # start with |1_0⟩
    for q, theta_q in enumerate(angles):
        # XXPlusYYGate(theta=2θ, beta=0) implements Givens rotation in single-exc subspace
        # Specifically: |10⟩ → cos(θ) |10⟩ − i sin(θ) |01⟩  for theta=2θ
        # We want |10⟩ → cos(θ) |10⟩ + sin(θ) |01⟩, so adjust phase via β=−π/2
        qc.append(XXPlusYYGate(2 * theta_q, beta=-math.pi / 2), [q, q + 1])
    qc.measure(range(N), cr)

    if with_noise:
        nm = NoiseModel()
        T1, T2 = 250e-6, 240e-6
        therm_1q = thermal_relaxation_error(T1, T2, 50e-9)
        depol_1q = depolarizing_error(0.0005, 1)
        err_1q = therm_1q.compose(depol_1q)
        therm_2q_a = thermal_relaxation_error(T1, T2, 200e-9)
        therm_2q_b = thermal_relaxation_error(T1, T2, 200e-9)
        err_2q = therm_2q_a.tensor(therm_2q_b).compose(depolarizing_error(0.005, 2))
        nm.add_all_qubit_quantum_error(err_1q, ['x', 'sx', 'rz', 'h', 's', 'sdg', 'z', 'y'])
        nm.add_all_qubit_quantum_error(err_2q, ['cx', 'cz'])
        ro_err = ReadoutError([[0.97, 0.03], [0.03, 0.97]])
        for q in range(N):
            nm.add_readout_error(ro_err, [q])
        sim = AerSimulator(noise_model=nm)
    else:
        sim = AerSimulator()

    qc_t = transpile(qc, sim, optimization_level=1)
    # Count gates
    gates = qc_t.count_ops()
    n_2q = gates.get('cx', 0) + gates.get('cz', 0) + gates.get('ecr', 0)

    result = sim.run(qc_t, shots=n_shots).result()
    counts = result.get_counts()
    total = sum(counts.values())

    # Per-qubit population P(qubit j = 1)
    pop = np.zeros(N)
    for bitstring, c in counts.items():
        # Qiskit bitstring is little-endian: rightmost char is qubit 0
        bits = bitstring.replace(' ', '')[::-1]  # reverse to make leftmost = qubit 0
        for j, b in enumerate(bits):
            if b == '1':
                pop[j] += c
    pop /= total

    return {
        'N': N, 'k': k,
        'shots': n_shots,
        'gates_2q': n_2q,
        'depth': qc_t.depth(),
        'populations': pop,
        'predicted': bonding_k_populations(N, k),
    }


# ──────────────────────────────────────────────────────────────────────
# Analytical noise model (large N, when state vector too big)
# ──────────────────────────────────────────────────────────────────────

def predict_with_depolarization(N: int, k: int, gate_error: float = 0.005):
    """Predict measured populations under uniform depolarization.

    Approximate: state ρ = p |ψ⟩⟨ψ| + (1−p) I/2^N
    where p = (1 − gate_error)^(n_gates).
    Per-qubit population: P(j) = p · |c_j|² + (1−p) · 0.5.
    """
    n_gates = 2 * (N - 1)  # ~2 CZ per Givens, N-1 Givens
    p = (1 - gate_error) ** n_gates
    populations_ideal = bonding_k_populations(N, k)
    populations_predicted = p * populations_ideal + (1 - p) * 0.5
    return {
        'N': N, 'k': k,
        'fidelity_estimate': p,
        'gates_2q_estimate': n_gates,
        'predicted_populations_with_noise': populations_predicted,
        'predicted_populations_ideal': populations_ideal,
    }


# ──────────────────────────────────────────────────────────────────────
# Reporting
# ──────────────────────────────────────────────────────────────────────

def report_aer(result):
    N = result['N']
    k = result['k']
    pred = result['predicted']
    meas = result['populations']
    print(f"\n{'─' * 80}")
    print(f"N = {N}, bonding:{k}, shots = {result['shots']}, "
          f"2q gates = {result['gates_2q']}, depth = {result['depth']}")
    print(f"{'─' * 80}")
    print(f"\n{'site j':>8s} {'predicted P_ideal':>20s} {'measured P_aer':>18s} {'diff':>10s}")
    # Show first few, middle, last few
    indices = list(range(min(5, N))) + ([] if N <= 10 else [N // 4, N // 2, 3 * N // 4]) + list(range(max(0, N - 5), N))
    indices = sorted(set(indices))
    for j in indices:
        print(f"{j:>8d} {pred[j]:>20.5f} {meas[j]:>18.5f} {meas[j] - pred[j]:>10.5f}")

    # RMS deviation
    rms = math.sqrt(np.mean((meas - pred) ** 2))
    peak_pred = pred.max()
    print(f"\nRMS(meas − ideal) = {rms:.5f}")
    print(f"Peak predicted    = {peak_pred:.5f}")
    print(f"Peak measured     = {meas.max():.5f}")
    print(f"Signal-to-noise ratio (peak / RMS) = {peak_pred / rms if rms > 0 else float('inf'):.2f}")


def report_analytical(pred):
    N = pred['N']
    k = pred['k']
    p = pred['fidelity_estimate']
    print(f"\n{'─' * 80}")
    print(f"N = {N}, bonding:{k}: ANALYTICAL prediction (state vector too big)")
    print(f"{'─' * 80}")
    print(f"Estimated 2q gates: {pred['gates_2q_estimate']}")
    print(f"Estimated fidelity p = (1 − 0.005)^{pred['gates_2q_estimate']} = {p:.4f}")
    print(f"Population formula: P(j) = p · |c_j|² + (1−p) · 0.5")

    p_id = pred['predicted_populations_ideal']
    p_pred = pred['predicted_populations_with_noise']
    peak = p_id.argmax()
    trough_zero = (np.abs(p_id) < 1e-6).any()

    print(f"\nIdeal P_peak (j={peak})      = {p_id[peak]:.5f}")
    print(f"Predicted P_peak (with noise) = {p_pred[peak]:.5f}")
    print(f"Predicted P_floor (j with |c_j|≈0) = {p_pred.min():.5f}")
    print(f"Contrast (P_peak − P_floor)   = {p_pred[peak] - p_pred.min():.5f}")
    stat_err = math.sqrt(0.25 / 4096)
    print(f"Statistical error per site at 4096 shots ≈ {stat_err:.5f}")
    print(f"Signal/stat-noise ≈ {(p_pred[peak] - p_pred.min()) / stat_err:.2f}σ")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 80)
    print("Finale design: 150-qubit bonding:1 on Marrakesh")
    print("Aer validation at small N + analytical extrapolation to N=150")
    print("=" * 80)

    # Aer-tractable N
    for N in [10, 20, 30, 50]:
        for k in [1, 2]:
            result = simulate_aer(N, k, n_shots=4096, with_noise=True)
            report_aer(result)

    # Analytical extrapolation
    for N in [100, 150]:
        for k in [1, 2]:
            pred = predict_with_depolarization(N, k, gate_error=0.005)
            report_analytical(pred)

    print()
    print("=" * 80)
    print("Decision criteria:")
    print("  - At N=150, signal/stat-noise ≥ 5σ in the contrast P_peak − P_floor")
    print("    means hardware run is justified.")
    print("  - Below 3σ: drop to smaller N.")
    print("=" * 80)


if __name__ == "__main__":
    main()
