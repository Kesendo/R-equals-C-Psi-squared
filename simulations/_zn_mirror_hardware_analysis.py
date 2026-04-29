#!/usr/bin/env python3
"""Z⊗N-Mirror Hardware-Test: post-run analysis on JSON tomography data.

Usage:
  python simulations/_zn_mirror_hardware_analysis.py path/to/zn_mirror_data.json

Expected JSON structure (mirroring data/ibm_soft_break_april2026 format):
  {
    "mode": "hardware",
    "backend": "ibm_marrakesh",
    "path": [48, 49, 50],
    "job_id": "...",
    "parameters": {"N": 3, "J": 1.0, "t": 0.8, "n_trotter": 3, "shots": 4096},
    "expectations": {
      "state_a": {"I,I": 1.0, "X,X": ..., ..., "Z,Z": ...},
      "state_b": {"I,I": 1.0, "X,X": ..., ..., "Z,Z": ...}
    }
  }
"""
import sys
import json
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

import warnings
warnings.simplefilter('ignore')
import framework as fw

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def reconstruct_rho_2qubit(expect_dict):
    """ρ = (1/4) Σ_{a,b} <P_a P_b> · (P_a ⊗ P_b)."""
    rho = np.zeros((4, 4), dtype=complex)
    for a in 'IXYZ':
        for b in 'IXYZ':
            key = f"{a},{b}"
            if key in expect_dict:
                rho = rho + expect_dict[key] * np.kron(PAULI[a], PAULI[b])
    return rho / 4.0


def per_pauli_zn_test(expect_a, expect_b):
    """For each Pauli string, compare HW <P>_b vs (-1)^n_XY · <P>_a."""
    rows = []
    for a in 'IXYZ':
        for b in 'IXYZ':
            if (a, b) == ('I', 'I'):
                continue
            key = f"{a},{b}"
            if key not in expect_a or key not in expect_b:
                continue
            n_xy = (1 if a in 'XY' else 0) + (1 if b in 'XY' else 0)
            sign = (-1) ** n_xy
            e_a = expect_a[key]
            e_b = expect_b[key]
            expected_b = sign * e_a
            violation = abs(e_b - expected_b)
            rows.append((key, n_xy, sign, e_a, e_b, expected_b, violation))
    return rows


def main(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)

    print(f"=== Z⊗N-Mirror Hardware Analysis ===\n")
    print(f"  Backend: {data.get('backend', '?')}")
    print(f"  Job ID:  {data.get('job_id', '?')}")
    print(f"  Path:    {data.get('path', '?')}")
    print(f"  Params:  {data.get('parameters', '?')}")
    print()

    if 'expectations' not in data:
        print("ERROR: 'expectations' key missing in JSON")
        return
    if 'state_a' not in data['expectations'] or 'state_b' not in data['expectations']:
        print("ERROR: expectations must have keys 'state_a' and 'state_b'")
        return

    expect_a = data['expectations']['state_a']
    expect_b = data['expectations']['state_b']

    # Reduced 2-qubit density matrices on (q_first, q_last)
    rho_a_2q = reconstruct_rho_2qubit(expect_a)
    rho_b_2q = reconstruct_rho_2qubit(expect_b)

    print("=== Per-Pauli Z⊗N test on (q0, q2) tomography ===\n")
    print(f"{'Pauli':<8s} {'n_XY':>5s} {'sign':>5s} {'<P>_a':>10s} {'<P>_b':>10s} "
          f"{'expect b':>10s} {'violation':>10s}")
    print("-" * 70)
    rows = per_pauli_zn_test(expect_a, expect_b)
    rows_sorted = sorted(rows, key=lambda r: -r[6])
    max_viol = 0.0
    sum_viol_sq = 0.0
    for key, n_xy, sign, e_a, e_b, e_exp, viol in rows_sorted:
        marker = '←' if viol > 0.05 else ''
        print(f"  {key:<6s} {n_xy:>5d} {sign:>+5d} {e_a:>10.4f} {e_b:>10.4f} "
              f"{e_exp:>10.4f} {viol:>10.4f}  {marker}")
        max_viol = max(max_viol, viol)
        sum_viol_sq += viol**2

    rms_viol = np.sqrt(sum_viol_sq / len(rows))
    print()
    print(f"  Max violation (over all 15 Pauli strings):  {max_viol:.6f}")
    print(f"  RMS violation (over all 15 Pauli strings):  {rms_viol:.6f}")

    # Effective transverse-field estimate (linear scaling, both hypotheses)
    print()
    print("=== Effective transverse-field estimate ===\n")
    print("Linear scalings from framework (see _zn_mirror_hardware_prediction.py):")
    print("  max_viol(h_x) ≈ 0.085 · h_x  (single-X bit_b-even, weak coupling)")
    print("  max_viol(h_y) ≈ 3.5 · h_y    (single-Y bit_b-odd, ~40× stronger)")
    print()
    h_x_est = max_viol / 0.085
    h_y_est = max_viol / 3.5
    print(f"  Hypothesis A (purely h_x): |h_x_eff| ≈ {h_x_est:.4f}")
    print(f"  Hypothesis B (purely h_y): |h_y_eff| ≈ {h_y_est:.4f}")
    print()
    print("  Which is more plausible? Y-Z-mixing structure of the worst-violating")
    print("  Pauli strings (e.g. Y,Z and Z,Y both broken implies Y-axis rotation)")
    print("  usually points to h_y rather than h_x.")
    print()
    if max_viol < 0.02:
        print("  Verdict: max_viol < 0.02 — consistent with NO transverse-field "
              "(or below readout-noise floor at 4096 shots).")
    elif max_viol < 0.05:
        print(f"  Verdict: max_viol ∈ [0.02, 0.05] — weak signal. Either small "
              f"h_y_eff ≈ {h_y_est:.4f}, or shot-noise/readout asymmetry. "
              f"Cross-check with more shots if budget allows.")
    elif h_y_est < 0.5 and h_x_est > 0.5:
        print(f"  Verdict: max_viol = {max_viol:.4f} → likely effective h_y_eff ≈ "
              f"{h_y_est:.4f} (h_x_eff ≈ {h_x_est:.4f} would be unrealistically large).")
        worst = rows_sorted[0][0]
        print(f"           Worst Pauli string: {worst}.")
    else:
        print(f"  Verdict: max_viol = {max_viol:.4f} → significant transverse-field. "
              f"h_x_eff ≈ {h_x_est:.4f} OR h_y_eff ≈ {h_y_est:.4f} "
              f"(inspect worst-Pauli structure to disambiguate).")
        worst = rows_sorted[0][0]
        print(f"           Worst Pauli string: {worst}.")

    # Diagnostic via the framework cockpit method (3-qubit version requires full ρ)
    print()
    print("=== Cockpit zn_mirror_diagnostic on reduced 2-qubit ρ ===\n")
    print("(Note: this runs on the (q0,q2) reduced state, not the full 3-qubit state.")
    print(" Full-3-qubit Z⊗N test would require full Pauli tomography.)")
    chain_2q = fw.ChainSystem(N=2, gamma_0=0.05)
    diag = chain_2q.zn_mirror_diagnostic(rho_a_2q, rho_b_2q, tol=0.02)
    print(f"  max_violation: {diag['max_violation']:.6f}")
    print(f"  verdict:       {diag['verdict']}")
    print(f"  worst_string:  {diag['worst_string']}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
