"""Run-planner for the C_block / Theorem-2-saturation IBM tomography run.

We want to probe the canonical Dicke superposition (|D_0⟩+|D_1⟩)/√2 at N=2
and watch C_block(t) decay under Z-dephasing-equivalent dynamics. The
trajectory under pure dephasing is closed-form:

    C_block(t) = (1/4) · exp(−4γ_eff · t)        (N=2, n=0; only HD=1 channel)

so a fit of the measured C_block(t) curve extracts γ_eff and validates the
Tier-1 bound (C_block(0) = 1/4 exactly per Theorem 1).

This script reads the 2026-05-08 morning IBM calibration CSVs (Marrakesh,
Kingston, Fez) and recommends per-backend the best 2-qubit pair for the run:
maximises min coherence, minimises CZ + single-qubit gate errors. Outputs an
estimated γ_eff = 1/T2_min and a t-grid spanning the decoherence window so
the trajectory has both the ceiling-saturation point (t ≈ 0) and a clear
decay tail.
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ibm_calibration import (  # noqa: E402
    QubitData,
    best_chain,
    chain_score,
    load_calibration,
    score_qubit,
)


CALIB_DIR = Path(
    r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI"
    r"\experiments\ibm_quantum_tomography\Calibrations"
)
CALIB_FILES = {
    "ibm_marrakesh": CALIB_DIR / "ibm_marrakesh_calibrations_2026-05-08T01_54_54Z.csv",
    "ibm_kingston":  CALIB_DIR / "ibm_kingston_calibrations_2026-05-08T02_37_30Z.csv",
    "ibm_fez":       CALIB_DIR / "ibm_fez_calibrations_2026-05-08T02_21_29Z.csv",
}

RESULTS_FILE = Path(__file__).resolve().parent / "results" / "_block_cpsi_run_planner_2026_05_08.txt"


def t_grid_for_t2(t2_min_us: float, n_points: int = 5) -> list[float]:
    """Recommend evolution times spanning the decoherence window. We want the
    full range from "fresh" (t ≈ 0) to "well-decayed" (t ≈ T2). Uniform
    spacing on [0, T2_min] with n_points samples gives a fit-friendly grid;
    a linear spacing here is sufficient because exp(−4γt) over [0, T2] covers
    [1.00 → e^{−4} = 0.018], so the curve is sampled across more than two
    orders of magnitude on a linear t-grid of 5 points."""
    return [k * t2_min_us / (n_points - 1) for k in range(n_points)]


def predict_c_block(gamma_eff_per_us: float, t_us: float) -> float:
    """Closed form for N=2, n=0 (c=1, single HD=1 channel): C_block(t) =
    (1/4) · exp(−4γ_eff · t)."""
    import math
    return 0.25 * math.exp(-4.0 * gamma_eff_per_us * t_us)


def best_pair_for_backend(qubits: list[QubitData]) -> tuple[list[int], float, QubitData, QubitData, float]:
    """Best 2-qubit chain on this backend. Returns
    (path = [q1, q2], chain_score, qubit1_data, qubit2_data, cz_error)."""
    path, score = best_chain(qubits, length=2)
    by_id = {q.qubit: q for q in qubits}
    q1 = by_id[path[0]]
    q2 = by_id[path[1]]
    cz_err = q1.cz_neighbours.get(q2.qubit, float("nan"))
    return path, score, q1, q2, cz_err


def main() -> None:
    out: list[str] = []

    def emit(line: str = "") -> None:
        out.append(line)
        print(line)

    emit("C_block / Theorem-2-saturation run planner")
    emit("=" * 80)
    emit(f"Calibration source: {CALIB_DIR}")
    emit("Run target: |ψ_init⟩ = (|D_0⟩+|D_1⟩)/√2 = (|00⟩ + (|01⟩+|10⟩)/√2)/√2 at N=2.")
    emit("Theorem 1: C_block(t=0) = 1/4 exactly (Mandelbrot-cardioid saturation).")
    emit("Theorem 2 trajectory under Z-dephasing: C_block(t) = (1/4) · exp(−4γ_eff · t).")
    emit()

    overall_recommendation: tuple[float, str, list[int], float, float, float] | None = None

    for backend, csv_path in CALIB_FILES.items():
        qubits = load_calibration(csv_path)
        n_op = sum(1 for q in qubits if q.operational)
        path, score, q1, q2, cz_err = best_pair_for_backend(qubits)
        t2_min = min(q1.t2_us, q2.t2_us)
        t1_min = min(q1.t1_us, q2.t1_us)
        gamma_eff = 1.0 / t2_min  # per microsecond
        sx_err_sum = q1.sx_error + q2.sx_error
        readout_err_max = max(q1.readout_error, q2.readout_error)

        t_grid = t_grid_for_t2(t2_min, n_points=5)

        emit(f"=== {backend} ===")
        emit(f"  source: {csv_path.name}  ({n_op} operational qubits)")
        emit(f"  best pair: q{path[0]}–q{path[1]}  (chain score {score:.3f})")
        emit(f"    q{q1.qubit}: T1 = {q1.t1_us:6.1f} μs, T2 = {q1.t2_us:6.1f} μs, "
             f"sx err = {q1.sx_error:.2e}, readout err = {q1.readout_error:.3f}")
        emit(f"    q{q2.qubit}: T1 = {q2.t1_us:6.1f} μs, T2 = {q2.t2_us:6.1f} μs, "
             f"sx err = {q2.sx_error:.2e}, readout err = {q2.readout_error:.3f}")
        emit(f"  CZ({path[0]},{path[1]}) error: {cz_err:.2e}")
        emit(f"  pair coherence floor: T1_min = {t1_min:6.1f} μs, T2_min = {t2_min:6.1f} μs")
        emit(f"  effective Z-dephasing rate: γ_eff = 1/T2_min = {gamma_eff:.4e} μs⁻¹")
        emit(f"  recommended t-grid (μs): {[f'{t:.1f}' for t in t_grid]}")
        emit(f"  predicted C_block(t):    {[f'{predict_c_block(gamma_eff, t):.4f}' for t in t_grid]}")
        emit()

        # Overall recommendation: pick the backend with highest chain score.
        if overall_recommendation is None or score > overall_recommendation[0]:
            overall_recommendation = (score, backend, path, t2_min, gamma_eff, cz_err)

    emit("=" * 80)
    emit("Overall recommendation (highest 2-pair chain score):")
    if overall_recommendation is not None:
        score, backend, path, t2_min, gamma_eff, cz_err = overall_recommendation
        emit(f"  → {backend} on q{path[0]}–q{path[1]} "
             f"(score {score:.3f}, T2_min = {t2_min:.1f} μs, "
             f"γ_eff = {gamma_eff:.4e} μs⁻¹, CZ err = {cz_err:.2e})")
    emit()

    emit("Run-cost outline (rough):")
    emit("  - state prep: 1 RY + 1 CX (≈ 100 ns wall-time, gate-error ≈ sx_err + cz_err)")
    emit("  - 16-Pauli tomography (2 qubits): 9 unique measurement settings (II is trivial,")
    emit("    other 15 split into 9 mutually-commuting groups for simultaneous readout)")
    emit("  - per measurement setting: ~2048 shots is standard for ~1% per-Pauli precision")
    emit("  - full snapshot at one t-point: 9 × ~5–6 s/job (Heron r2 session overhead per")
    emit("    reference_ibm_qpu_billing) ≈ 50 s")
    emit("  - 5 t-points × 1 backend = 250 s ≈ 4 min QPU")
    emit("  - 5 t-points × 3 backends = 12 min (within 15 min/month budget)")
    emit()
    emit("Decisive output of the run:")
    emit("  1. Theorem 1 saturation: at t=0, hardware ρ_q0q1 must satisfy C_block ≈ 1/4")
    emit("     (state-prep + readout fidelity bound; no decoherence).")
    emit("  2. Theorem 2 trajectory fit: log(C_block(t)) vs t is linear with slope −4γ_fit;")
    emit("     compare γ_fit against γ_eff = 1/T2_min from calibration.")
    emit("     A clean fit confirms the universal-shape closed form (Tier-1-grade evidence).")
    emit("  3. Per-backend γ_fit comparison: differences across Marrakesh/Kingston/Fez")
    emit("     diagnose whether T2 alone explains decoherence or extra noise channels.")

    RESULTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    RESULTS_FILE.write_text("\n".join(out), encoding="utf-8")
    print(f"\nResults written to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
