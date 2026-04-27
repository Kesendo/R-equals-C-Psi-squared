#!/usr/bin/env python3
"""Domain 2: re-read all existing IBM hardware data through the 4-panel
cockpit lens. Compute per-(backend, H-case) saturation ratios — how close
each hardware run came to the theoretical Π-protection ceiling.

Source data (all in the IBM tomography results dir):
  framework_snapshots_ibm_*.json   (7 Snapshot-D runs across 3 Heron r2)
    Each contains 9-Pauli tomography on (q0, q2) for truly/soft/hard
    Hamiltonians at t=0.8, |+−+⟩.
  lebensader_cusp_ibm_marrakesh_*.json   (1 Lebensader run, 2H × 4t)
  iy_yi_completion_ibm_marrakesh_*.json  (1 IY+YI run, 1H × 4t)

For each (backend, H-case, time):
  - Reconstruct ρ_{q0,q2} from 16 Pauli expectations
  - Compute simulator prediction (continuous Lindblad +T1)
  - Identify cells where |sim| < 0.05 (predicted-zero) and
    cells where |hw| < 0.20 (observed-zero)
  - Saturation ratio = predicted-zero ∩ observed-zero / predicted-zero

Cross-tabulate by (H-case, t, backend).
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

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}

PRED_ZERO_THRESHOLD = 0.05    # simulator value below this → predicted zero
NOISE_THRESHOLD = 0.20         # hardware value below this → observed zero


def simulate_q0q2_at_t(H, gamma_l, gamma_t1_l, rho_0_full, t, N=3):
    """Compute simulator prediction for all 16 Pauli expectations on (q0, q2)."""
    L = fw.lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    rho_vec0 = rho_0_full.T.reshape(-1).copy()
    rho_t_vec = expm(L * t) @ rho_vec0
    rho_t = rho_t_vec.reshape(2 ** N, 2 ** N).T

    out = {}
    for a in 'IXYZ':
        for b in 'IXYZ':
            P = np.kron(PAULI[a], np.kron(PAULI['I'], PAULI[b]))
            out[f"{a},{b}"] = float(np.real(np.trace(rho_t @ P)))
    return out


def load_snapshot_d_files():
    """Load Snapshot-D from all framework_snapshots_ibm_*.json."""
    out = []
    for json_path in sorted(RESULTS_DIR.glob("framework_snapshots_ibm_*.json")):
        with open(json_path, encoding='utf-8') as f:
            d = json.load(f)
        if 'snapshot_d_softbreak_trichotomy' not in d:
            continue
        backend = d.get('backend', json_path.stem)
        ts = json_path.stem.split('_')[-1]
        params = d.get('parameters', {})
        out.append({
            'file': json_path.name,
            'backend': f"{backend}_{ts}",
            'param_t': params.get('t_eval', 0.8),
            'expectations': d['snapshot_d_softbreak_trichotomy']['expectations_per_category'],
            'source': 'snapshot_d',
        })
    return out


def load_lebensader_cusp():
    out = []
    for json_path in sorted(RESULTS_DIR.glob("lebensader_cusp_ibm_*.json")):
        with open(json_path, encoding='utf-8') as f:
            d = json.load(f)
        out.append({
            'file': json_path.name,
            'backend': f"{d.get('backend')}_{json_path.stem.split('_')[-1]}",
            'expectations_per_h_t': d.get('hardware_expectations', {}),
            'cases': d.get('cases', {}),
            'source': 'lebensader_cusp',
        })
    return out


def load_iy_yi_completion():
    out = []
    for json_path in sorted(RESULTS_DIR.glob("iy_yi_completion_ibm_*.json")):
        with open(json_path, encoding='utf-8') as f:
            d = json.load(f)
        out.append({
            'file': json_path.name,
            'backend': f"{d.get('backend')}_{json_path.stem.split('_')[-1]}",
            'expectations_per_t': d.get('hardware_expectations', {}),
            'source': 'iy_yi',
        })
    return out


def saturation_for(sim_dict, hw_dict):
    """Compute saturation ratio for a (sim, hw) pair on 16 cells."""
    pred_zero_cells = set()
    pred_active_cells = set()
    obs_zero_cells = set()
    obs_active_cells = set()
    cells_with_data = []

    for k in sim_dict:
        if k == 'I,I':
            continue  # always 1
        sim_val = sim_dict[k]
        hw_val = hw_dict.get(k)
        if hw_val is None:
            continue
        cells_with_data.append((k, sim_val, hw_val))
        if abs(sim_val) < PRED_ZERO_THRESHOLD:
            pred_zero_cells.add(k)
        else:
            pred_active_cells.add(k)
        if abs(hw_val) < NOISE_THRESHOLD:
            obs_zero_cells.add(k)
        else:
            obs_active_cells.add(k)

    n_pred_zero = len(pred_zero_cells)
    n_obs_zero = len(obs_zero_cells)
    n_pred_zero_obs_zero = len(pred_zero_cells & obs_zero_cells)
    sat_ratio = n_pred_zero_obs_zero / max(n_pred_zero, 1)

    n_pred_active = len(pred_active_cells)
    n_pred_active_obs_active = len(pred_active_cells & obs_active_cells)
    active_recall = n_pred_active_obs_active / max(n_pred_active, 1)

    return {
        'n_pred_zero': n_pred_zero,
        'n_obs_zero': n_obs_zero,
        'n_pred_zero_obs_zero': n_pred_zero_obs_zero,
        'saturation': sat_ratio,
        'n_pred_active': n_pred_active,
        'n_pred_active_obs_active': n_pred_active_obs_active,
        'active_recall': active_recall,
        'cells_with_data': cells_with_data,
    }


def make_rho_xneel(N=3):
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    return np.outer(psi, psi.conj())


def main():
    N = 3
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]
    rho_xneel = make_rho_xneel(N)

    print(f"Domain 2: existing IBM data through 4-panel cockpit lens")
    print(f"  Predicted-zero threshold (sim): |sim| < {PRED_ZERO_THRESHOLD}")
    print(f"  Hardware noise threshold:       |hw|  < {NOISE_THRESHOLD}")
    print()

    # ========== Snapshot D ==========
    print(f"=" * 110)
    print(f"Snapshot D: 7 backends × 3 H-cases (truly XX+YY, soft XY+YX, hard XX+XY)")
    print(f"  Initial: |+−+⟩, t = 0.8")
    print(f"=" * 110)
    print()

    snap_d_data = load_snapshot_d_files()
    H_cases_d = {
        'truly': [('X', 'X', J), ('Y', 'Y', J)],
        'soft':  [('X', 'Y', J), ('Y', 'X', J)],
        'hard':  [('X', 'X', J), ('X', 'Y', J)],
    }

    # Pre-compute simulator predictions
    sim_pred = {}
    for cat, terms in H_cases_d.items():
        H = fw._build_bilinear(N, bonds, terms)
        sim_pred[cat] = simulate_q0q2_at_t(
            H, [GAMMA] * N, [GAMMA_T1] * N, rho_xneel, 0.8, N=N
        )

    print(f"  {'backend':<22s}  {'cat':<6s}  "
          f"{'pred|0':>6s}  {'obs|0':>5s}  {'sat ratio':>10s}  "
          f"{'pred|act':>8s}  {'recall':>7s}")
    print('-' * 95)

    overall_results = []
    for snap in snap_d_data:
        backend = snap['backend']
        for cat in ['truly', 'soft', 'hard']:
            hw_dict = snap['expectations'].get(cat, {})
            if not hw_dict:
                continue
            sat = saturation_for(sim_pred[cat], hw_dict)
            print(f"  {backend:<22s}  {cat:<6s}  "
                  f"{sat['n_pred_zero']:>6d}  {sat['n_obs_zero']:>5d}  "
                  f"{sat['saturation']:>9.0%}  "
                  f"{sat['n_pred_active']:>8d}  {sat['active_recall']:>6.0%}")
            overall_results.append({
                'source': 'snap_d', 'backend': backend, 'cat': cat,
                **sat,
            })
        print()

    # Aggregate per category
    print()
    print(f"Aggregate per category over 7 Snapshot-D backends:")
    print(f"  {'cat':<6s}  {'mean sat':>9s}  {'std sat':>8s}  {'mean recall':>11s}")
    for cat in ['truly', 'soft', 'hard']:
        cat_rows = [r for r in overall_results if r['cat'] == cat]
        if not cat_rows:
            continue
        mean_sat = float(np.mean([r['saturation'] for r in cat_rows]))
        std_sat = float(np.std([r['saturation'] for r in cat_rows]))
        mean_recall = float(np.mean([r['active_recall'] for r in cat_rows]))
        print(f"  {cat:<6s}  {mean_sat:>8.0%}  {std_sat:>7.0%}  "
              f"{mean_recall:>10.0%}")
    print()

    # ========== Lebensader cusp ==========
    print(f"=" * 110)
    print(f"Lebensader cusp: Marrakesh, 2 H (XY+YX, YZ+ZY) × 4 t")
    print(f"=" * 110)
    print()

    leb_data = load_lebensader_cusp()
    H_cases_leb = {
        'XY+YX': [('X', 'Y', J), ('Y', 'X', J)],
        'YZ+ZY': [('Y', 'Z', J), ('Z', 'Y', J)],
    }
    print(f"  {'backend':<22s}  {'H':<6s}  {'t':>4s}  "
          f"{'pred|0':>6s}  {'obs|0':>5s}  {'sat':>5s}  {'recall':>7s}")
    print('-' * 80)
    for snap in leb_data:
        for h_label, terms in H_cases_leb.items():
            H = fw._build_bilinear(N, bonds, terms)
            for t_str in ['0.8', '1.4', '1.7', '2.2']:
                hw_dict = snap['expectations_per_h_t'].get(h_label, {}).get(t_str, {})
                if not hw_dict:
                    continue
                t = float(t_str)
                sim = simulate_q0q2_at_t(
                    H, [GAMMA] * N, [GAMMA_T1] * N, rho_xneel, t, N=N
                )
                sat = saturation_for(sim, hw_dict)
                print(f"  {snap['backend']:<22s}  {h_label:<6s}  {t:>4.1f}  "
                      f"{sat['n_pred_zero']:>6d}  {sat['n_obs_zero']:>5d}  "
                      f"{sat['saturation']:>4.0%}  "
                      f"{sat['active_recall']:>6.0%}")
    print()

    # ========== IY+YI completion ==========
    print(f"=" * 110)
    print(f"IY+YI completion: Marrakesh, 1 H × 4 t")
    print(f"=" * 110)
    print()

    iy_data = load_iy_yi_completion()
    print(f"  {'backend':<22s}  {'H':<6s}  {'t':>4s}  "
          f"{'pred|0':>6s}  {'obs|0':>5s}  {'sat':>5s}  {'recall':>7s}")
    print('-' * 80)
    H_iyyi = fw._build_bilinear(N, bonds, [('I', 'Y', J), ('Y', 'I', J)])
    for snap in iy_data:
        for t_str in ['0.8', '1.4', '1.7', '2.2']:
            hw_dict = snap['expectations_per_t'].get(t_str, {})
            if not hw_dict:
                continue
            t = float(t_str)
            sim = simulate_q0q2_at_t(
                H_iyyi, [GAMMA] * N, [GAMMA_T1] * N, rho_xneel, t, N=N
            )
            sat = saturation_for(sim, hw_dict)
            print(f"  {snap['backend']:<22s}  {'IY+YI':<6s}  {t:>4.1f}  "
                  f"{sat['n_pred_zero']:>6d}  {sat['n_obs_zero']:>5d}  "
                  f"{sat['saturation']:>4.0%}  "
                  f"{sat['active_recall']:>6.0%}")

    print()
    print("=" * 110)
    print("Summary:")
    print("  saturation = predicted-zero cells that hardware actually shows in noise band")
    print("  recall     = predicted-active cells that hardware shows above noise band")
    print("  100% in both columns = framework lens fully matches hardware reality")


if __name__ == "__main__":
    main()
