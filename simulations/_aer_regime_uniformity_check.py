"""Aer regime-uniformity check with per-qubit T1/T2 profiles (no QPU).

Follow-up to `_lindblad_regime_uniformity_check.py` (commit 02ce722) which
ruled out pure Z-dephasing + Trotter as the source of the 22.8× truly-baseline
mixing penalty. This script adds amplitude damping (T1) and gate errors via
Qiskit Aer's NoiseModel, with per-qubit T1/T2 to test whether non-uniform
T1 or T2 across the chain reproduces the hardware effect.

Profiles tested (4 Hamiltonian categories × N profiles, all at N=3, t=0.8,
n_trotter=3, 4096 shots, |+−+⟩ initial state):

  uniform-classical:  T1 = 300μs, T2 = 240μs uniform (matches r ≈ 0.4)
  uniform-quantum:    T1 = 300μs, T2 = 50μs uniform  (matches r ≈ 0.083)
  mixed-T2-q0:        T1 = 300μs uniform, T2 = [50, 240, 240]μs
  mixed-T2-q1:        T1 = 300μs uniform, T2 = [240, 50, 240]μs
  mixed-T1T2-q0:      T1 = [50, 300, 300]μs, T2 = [50, 240, 240]μs (degraded q0)

Reading: which profile (if any) reproduces the hardware mixing penalty?
The Lindblad-only check showed T2-only variation can't do it. If
mixed-T1T2-q0 reproduces it, the F82 amplitude-damping mechanism is
identified. If neither does, the source is gate errors / crosstalk /
readout / Trotter+T1 interaction in some way the simple per-qubit
thermal_relaxation_error doesn't capture.

Hardware anchors (Marrakesh + Kingston, 2026-04-26 / 05-05):
  uniform-classical [48, 49, 50] truly = 0.0013
  uniform-quantum   [43, 56, 63] truly = 0.0022
  regime-mixed      [0, 1, 2]    truly = 0.0297  (~22.8× uniform-classical)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

# Wire to AIEvolution's run_soft_break helpers
AIE_DIR = Path("D:/Entwicklung/Projekte/.NET Projekte/AIEvolution/AIEvolution.UI/experiments/ibm_quantum_tomography")
sys.path.insert(0, str(AIE_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_soft_break import (
    N, J, T_EVAL, N_TROTTER, CATEGORIES,
    hamiltonian_to_sparse_pauli, prepare_xneel_circuit, build_circuit,
    expectations_from_counts,
)
from _f88_lens_ibm_framework_snapshots import reconstruct_2qubit_rho, f88_lens_2qubit

from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (
    NoiseModel, depolarizing_error, thermal_relaxation_error, ReadoutError,
)


SHOTS = 4096


def build_per_qubit_simulator(t1_per_qubit_us: List[float],
                               t2_per_qubit_us: List[float],
                               n: int = N) -> AerSimulator:
    """AerSimulator with per-qubit thermal_relaxation_error. Two-qubit gates
    use the average of the two participating qubits' T1/T2 (rough but
    consistent across profiles)."""
    nm = NoiseModel()

    GATE_1Q_S = 50e-9
    GATE_2Q_S = 200e-9

    for q in range(n):
        t1_s = t1_per_qubit_us[q] * 1e-6
        t2_s = min(t2_per_qubit_us[q] * 1e-6, 2 * t1_s)
        therm = thermal_relaxation_error(t1_s, t2_s, GATE_1Q_S)
        depol = depolarizing_error(0.0005, 1)
        err_1q = therm.compose(depol)
        nm.add_quantum_error(
            err_1q, ['x', 'sx', 'rz', 'h', 's', 'sdg', 'z', 'y'], qubits=[q],
        )

    for qa in range(n):
        for qb in range(n):
            if qa == qb:
                continue
            t1_a = t1_per_qubit_us[qa] * 1e-6
            t2_a = min(t2_per_qubit_us[qa] * 1e-6, 2 * t1_a)
            t1_b = t1_per_qubit_us[qb] * 1e-6
            t2_b = min(t2_per_qubit_us[qb] * 1e-6, 2 * t1_b)
            therm_a = thermal_relaxation_error(t1_a, t2_a, GATE_2Q_S)
            therm_b = thermal_relaxation_error(t1_b, t2_b, GATE_2Q_S)
            therm = therm_a.tensor(therm_b)
            depol = depolarizing_error(0.005, 2)
            err_2q = therm.compose(depol)
            nm.add_quantum_error(err_2q, ['cx', 'cz'], qubits=[qa, qb])

    ro_matrix = [[0.97, 0.03], [0.03, 0.97]]
    ro_err = ReadoutError(ro_matrix)
    for q in range(n):
        nm.add_readout_error(ro_err, [q])

    return AerSimulator(noise_model=nm, method='density_matrix')


def run_one_profile(t1_us: List[float], t2_us: List[float], shots: int = SHOTS) -> Dict[str, Dict]:
    """Run all 4 categories × 9 Pauli bases under the given T1/T2 profile.
    Returns expectations per category (16-Pauli reduced ρ on q0, q2)."""
    sim = build_per_qubit_simulator(t1_us, t2_us)
    bonds = [(0, 1), (1, 2)]
    prep_circ = prepare_xneel_circuit(N)
    bases = ['X', 'Y', 'Z']

    expectations_per_category = {}
    for category, terms in CATEGORIES:
        ham_op = hamiltonian_to_sparse_pauli(N, bonds, terms, J=J)
        counts_by_basis = {}
        for b0 in bases:
            for b2 in bases:
                qc = build_circuit(N, prep_circ, ham_op, T_EVAL, N_TROTTER, b0, b2)
                qc_t = transpile(qc, sim, optimization_level=1)
                result = sim.run(qc_t, shots=shots).result()
                counts_by_basis[(b0, b2)] = result.get_counts()
        expectations_per_category[category] = expectations_from_counts(counts_by_basis)
    return expectations_per_category


def f88_lens_per_category(expectations_per_category: Dict[str, Dict]) -> Dict[str, float]:
    """Convert per-category Pauli-pair expectations to F88-Lens Π²-odd-memory."""
    out = {}
    for cat, exps in expectations_per_category.items():
        exp_str = {f"{k[0]},{k[1]}": v for k, v in exps.items()}
        rho = reconstruct_2qubit_rho(exp_str)
        out[cat] = f88_lens_2qubit(rho)["pi2_odd_in_memory"]
    return out


PROFILES = {
    "uniform-classical":     ([300, 300, 300], [240, 240, 240]),
    "uniform-quantum":       ([300, 300, 300], [ 50,  50,  50]),
    "mixed-T2-q0":           ([300, 300, 300], [ 50, 240, 240]),
    "mixed-T2-q1":           ([300, 300, 300], [240,  50, 240]),
    "mixed-T2-q2":           ([300, 300, 300], [240, 240,  50]),
    "mixed-T1T2-q0":         ([ 50, 300, 300], [ 50, 240, 240]),
    # Actual Marrakesh Apr-25 calibration values (in μs):
    "marrakesh-[0,1,2]":     ([277.8, 177.5, 320.8], [ 41.4, 160.8, 406.4]),
    "marrakesh-[48,49,50]":  ([225.8, 290.5, 446.3], [231.7, 298.7, 188.6]),
    # Approximate Kingston May-5 calibration values for [43,56,63]:
    "kingston-[43,56,63]":   ([310.8, 251.7, 396.3], [ 44.5,  33.1,  57.6]),
}


def main():
    print("Aer per-qubit T1/T2 regime-uniformity scan")
    print(f"N = {N}, J = {J}, t = {T_EVAL}, n_trotter = {N_TROTTER}, shots = {SHOTS}")
    print("=" * 78)
    print()
    print("Hardware anchors:")
    print("  uniform-classical [48, 49, 50] truly = 0.0013")
    print("  uniform-quantum   [43, 56, 63] truly = 0.0022")
    print("  regime-mixed      [0, 1, 2]    truly = 0.0297  (22.8× uniform-classical)")
    print()

    results = {}
    for profile_name, (t1, t2) in PROFILES.items():
        print(f"Running profile {profile_name}: T1 = {t1} μs, T2 = {t2} μs ...")
        exps = run_one_profile(t1, t2)
        lens = f88_lens_per_category(exps)
        results[profile_name] = lens

    print()
    print("=" * 78)
    print("F88-Lens Π²-odd-memory readings per profile + category")
    print("=" * 78)
    print()
    cats = ["truly_unbroken", "pi2_odd_pure", "pi2_even_nontruly", "mixed_anti_one_sixth"]
    print(f"  {'profile':<22}", *(f"{c[:12]:>14}" for c in cats))
    print("  " + "-" * 80)
    for profile_name, lens in results.items():
        row = [f"{lens[c]:>14.4f}" for c in cats]
        print(f"  {profile_name:<22}", *row)

    print()
    print("=" * 78)
    print("READING")
    print()

    truly_uc = results["uniform-classical"]["truly_unbroken"]
    truly_uq = results["uniform-quantum"]["truly_unbroken"]
    print(f"truly-baseline:")
    for prof in PROFILES:
        v = results[prof]["truly_unbroken"]
        ratio_str = f"({v / truly_uc:.2f}× uniform-classical)" if truly_uc > 1e-9 else ""
        print(f"  {prof:<22} = {v:.4f}  {ratio_str}")

    print()
    max_mixed = max(
        results["mixed-T2-q0"]["truly_unbroken"],
        results["mixed-T2-q1"]["truly_unbroken"],
        results["mixed-T2-q2"]["truly_unbroken"],
        results["mixed-T1T2-q0"]["truly_unbroken"],
    )
    if truly_uc > 1e-9 and max_mixed > 5 * truly_uc:
        print(f"  → A mixed profile (max ratio {max_mixed / truly_uc:.1f}×) reproduces a substantial")
        print(f"    fraction of the hardware 22.8× mixing penalty. The mechanism is per-qubit")
        print(f"    T1/T2 variation captured by Aer's thermal_relaxation_error.")
    else:
        print(f"  → No mixed profile reproduces the 22.8× hardware penalty in Aer (max")
        print(f"    {max_mixed / max(truly_uc, 1e-9):.1f}×). The mechanism is not pure per-qubit")
        print(f"    thermal_relaxation_error; candidates remaining: gate-error variation,")
        print(f"    crosstalk, readout, or the interaction of multiple noise channels.")


if __name__ == "__main__":
    main()
