#!/usr/bin/env python3
"""EQ-030 fine question on hardware (Option A, zero QPU): reconstruct CΨ
and θ on the (q0, q2) reduced 2-qubit state from existing 16-Pauli
tomography snapshots.

Each framework_snapshots_ibm_*.json contains 'expectations_per_category'
with all 16 Pauli expectations ⟨P_q0 ⊗ P_q2⟩ for P ∈ {I, X, Y, Z}.
That is exactly what's needed to reconstruct the reduced 2-qubit density
matrix:
    ρ_02 = (1/4) Σ_{P_a, P_b} ⟨P_a, P_b⟩ · (P_a ⊗ P_b)

From ρ_02, compute:
    Purity = Tr(ρ_02²)
    Ψ-norm = L₁(ρ_02) / (d−1) = L₁ / 3
    CΨ = Purity × Ψ-norm
    θ = arctan(√(4·CΨ − 1)) if CΨ > 1/4, else 0

Compare to the simulator prediction at t=0.8 under (a) pure-Z, (b) +T1.
The discrepancy is the hardware noise floor; the closeness tells us
whether the cusp structure survives at the measured time point.

Important caveat: these snapshots are at t=0.8 (early), BEFORE all
boundary crossings (t_cross = 1.5-3.1). At t=0.8 we expect θ in the
20-50° range for all categories. The fragile-tail discrimination (θ
in the sub-degree regime, near t_cross) is NOT captured here. To
probe the tail we'd need NEW time-resolved tomography near t_cross.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw

RESULTS_DIR = Path(
    r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI"
    r"\experiments\ibm_quantum_tomography\results"
)

PAULI_2x2 = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def reconstruct_rho2q(expectations):
    """Reconstruct 2-qubit ρ from 16 ⟨P_a, P_b⟩ values keyed as 'A,B' strings.

    qiskit convention: 'A,B' means ⟨A_q0 ⊗ B_q2⟩ where q0 is the first
    measured qubit (cr[0]) and q2 the second (cr[1]). For ρ_{q0,q2}
    in standard tensor order |q0⟩|q2⟩, the corresponding Pauli operator
    is A_q0 ⊗ B_q2.
    """
    rho = np.zeros((4, 4), dtype=complex)
    for a in 'IXYZ':
        for b in 'IXYZ':
            key = f"{a},{b}"
            val = expectations.get(key, 0.0)
            P = np.kron(PAULI_2x2[a], PAULI_2x2[b])
            rho += val * P
    return rho / 4.0


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


def psi_norm(rho):
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1)


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def theta_from_cpsi(c):
    if c <= 0.25:
        return 0.0
    return math.degrees(math.atan(math.sqrt(4 * c - 1)))


def evolve_full3q(L, rho_0, t):
    d = rho_0.shape[0]
    rho_vec = rho_0.T.reshape(-1).copy()
    rho_t_vec = expm(L * t) @ rho_vec
    return rho_t_vec.reshape(d, d).T


def reduce_to_q0_q2(rho_3q):
    """Trace out the middle qubit (q1) to get ρ_{q0, q2}."""
    rho = rho_3q.reshape(2, 2, 2, 2, 2, 2)
    return np.einsum('ijkijk->', rho.reshape(8, 8)) * 0  # placeholder; redo below

def partial_trace_middle(rho, n_qubits=3):
    """For 3-qubit ρ in tensor order |q0 q1 q2⟩, trace out q1."""
    # Reshape ρ as (q0, q1, q2, q0', q1', q2') then sum over q1=q1'
    R = rho.reshape(2, 2, 2, 2, 2, 2)
    # Sum over q1 == q1' (axes 1 and 4)
    out = np.zeros((2, 2, 2, 2), dtype=complex)
    for q1 in range(2):
        out += R[:, q1, :, :, q1, :]
    return out.reshape(4, 4)


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0
    T_EVAL = 0.8
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())
    rho_0_q02 = partial_trace_middle(rho_0, 3)

    # Sanity: at t=0, ρ_02 = |+⟩⟨+| ⊗ |+⟩⟨+|
    p0 = purity(rho_0_q02)
    psi0 = psi_norm(rho_0_q02)
    c0 = cpsi(rho_0_q02)
    th0 = theta_from_cpsi(c0)
    print(f"Initial check: ρ_q0q2 (|+⟩⟨+|⊗|+⟩⟨+|)")
    print(f"  Purity={p0:.4f}, Ψ={psi0:.4f}, CΨ={c0:.4f}, θ={th0:.2f}°")
    print()

    # Build the 3 Hamiltonians
    cases = {
        'truly': [('X', 'X', J), ('Y', 'Y', J)],
        'soft':  [('X', 'Y', J), ('Y', 'X', J)],
        'hard':  [('X', 'X', J), ('X', 'Y', J)],
    }
    Hs = {k: fw._build_bilinear(N, bonds, v) for k, v in cases.items()}

    # Simulator predictions at t=0.8 under (a) pure-Z, (b) +T1
    print(f"Simulator predictions at t={T_EVAL}, |+−+⟩, "
          f"γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}:")
    print(f"  Each case: ρ_02 reduced to (q0,q2), then CΨ, θ.")
    print()
    sim_predictions = {}
    for cat, H in Hs.items():
        L_pure = fw.lindbladian_z_dephasing(H, [GAMMA_DEPH] * N)
        L_t1 = fw.lindbladian_z_plus_t1(H, [GAMMA_DEPH] * N, [GAMMA_T1] * N)
        rho_pure = evolve_full3q(L_pure, rho_0, T_EVAL)
        rho_t1 = evolve_full3q(L_t1, rho_0, T_EVAL)
        rho_pure_q02 = partial_trace_middle(rho_pure, 3)
        rho_t1_q02 = partial_trace_middle(rho_t1, 3)

        p_pure = purity(rho_pure_q02)
        c_pure = cpsi(rho_pure_q02)
        th_pure = theta_from_cpsi(c_pure)
        p_t1 = purity(rho_t1_q02)
        c_t1 = cpsi(rho_t1_q02)
        th_t1 = theta_from_cpsi(c_t1)

        sim_predictions[cat] = {
            'pure': {'purity': p_pure, 'cpsi': c_pure, 'theta': th_pure},
            't1':   {'purity': p_t1, 'cpsi': c_t1, 'theta': th_t1},
        }
        print(f"  {cat:>6s}  pure-Z: Pur={p_pure:.4f}, CΨ={c_pure:.4f}, θ={th_pure:.2f}°  "
              f"|  +T1: Pur={p_t1:.4f}, CΨ={c_t1:.4f}, θ={th_t1:.2f}°")
    print()

    # Hardware extraction
    print("Hardware: reconstruct ρ_q0q2 from 16-Pauli tomography per backend.")
    print()
    print(f"  {'backend':<22s}  {'cat':<6s}  "
          f"{'Pur_HW':>7s}  {'CΨ_HW':>7s}  {'θ_HW':>6s}  "
          f"{'CΨ_sim+T1':>10s}  {'θ_sim+T1':>9s}  {'Δθ':>7s}")
    print('-' * 100)

    rows = []
    for json_path in sorted(RESULTS_DIR.glob("framework_snapshots_ibm_*.json")):
        with open(json_path, encoding='utf-8') as f:
            data = json.load(f)
        if 'snapshot_d_softbreak_trichotomy' not in data:
            continue
        backend = data.get('backend', json_path.stem)
        ts = json_path.stem.split('_')[-1]
        key = f"{backend}_{ts}"
        d_data = data['snapshot_d_softbreak_trichotomy']['expectations_per_category']
        for cat in ['truly', 'soft', 'hard']:
            sub = d_data.get(cat)
            if sub is None:
                continue
            rho = reconstruct_rho2q(sub)
            p = purity(rho)
            cpsi_hw = cpsi(rho)
            th_hw = theta_from_cpsi(cpsi_hw)

            sim_t1 = sim_predictions[cat]['t1']
            dtheta = th_hw - sim_t1['theta']
            rows.append({
                'backend_key': key,
                'cat': cat,
                'purity_hw': p,
                'cpsi_hw': cpsi_hw,
                'theta_hw': th_hw,
                'theta_sim_t1': sim_t1['theta'],
            })
            print(f"  {key:<22s}  {cat:<6s}  "
                  f"{p:>7.4f}  {cpsi_hw:>7.4f}  {th_hw:>5.2f}°  "
                  f"{sim_t1['cpsi']:>10.4f}  {sim_t1['theta']:>8.2f}°  "
                  f"{dtheta:>+6.2f}°")
        print()

    # Aggregate per-category averages
    print()
    print("Per-category aggregate (mean over 7 backends):")
    print(f"  {'cat':<6s}  {'mean θ_HW':>10s}  {'std θ_HW':>9s}  "
          f"{'sim+T1 θ':>9s}  {'sim pure-Z θ':>13s}  {'reading':<28s}")
    for cat in ['truly', 'soft', 'hard']:
        cat_rows = [r for r in rows if r['cat'] == cat]
        if not cat_rows:
            continue
        mean_th = float(np.mean([r['theta_hw'] for r in cat_rows]))
        std_th = float(np.std([r['theta_hw'] for r in cat_rows]))
        sim_t1_th = sim_predictions[cat]['t1']['theta']
        sim_pure_th = sim_predictions[cat]['pure']['theta']
        if abs(mean_th - sim_t1_th) < abs(mean_th - sim_pure_th):
            reading = "closer to +T1 prediction"
        else:
            reading = "closer to pure-Z prediction"
        print(f"  {cat:<6s}  {mean_th:>9.2f}°  {std_th:>8.2f}°  "
              f"{sim_t1_th:>8.2f}°  {sim_pure_th:>12.2f}°  {reading:<28s}")

    print()
    print("Reading guide:")
    print("  At t=0.8 (before any crossing), θ should be 30-60° for all")
    print("  categories in simulator. Hardware lower than simulator means")
    print("  decoherence has eaten some θ. Hardware close to simulator")
    print("  means the cusp structure is intact at this early time.")
    print("  This snapshot does NOT probe the fragile tail (which lives")
    print("  near t_cross). It probes the bulk of the trajectory.")


if __name__ == "__main__":
    main()
