#!/usr/bin/env python3
"""Post-process the Lebensader hardware run on Marrakesh (Apr 26 evening).

Reads `lebensader_cusp_ibm_marrakesh_*.json`, reconstructs ρ_{q0,q2} from
the 16-Pauli expectations at each (H, t), and computes:
  - Purity, Ψ-norm, CΨ, θ on the reduced 2-qubit state
  - Compares hardware values to simulator predictions saved in the JSON

The Lebensader prediction at full N=3:
  XY+YX (drop=1)  → long fragile θ-tail
  YZ+ZY (drop=28) → no fragile tail, generic cusp

On the reduced (q0, q2) the picture is different (full N=3 ↛ reduced 2q),
but the simulator's predictions for the reduced state are saved alongside
hardware so we can do a direct apples-to-apples comparison.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

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
    rho = np.zeros((4, 4), dtype=complex)
    for a in 'IXYZ':
        for b in 'IXYZ':
            key = f"{a},{b}"
            val = expectations.get(key, 0.0)
            P = np.kron(PAULI_2x2[a], PAULI_2x2[b])
            rho += val * P
    return rho / 4.0


def metrics(rho):
    p = float(np.real(np.trace(rho @ rho)))
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi_n = l1 / (d - 1)
    cp = p * psi_n
    if cp <= 0.25:
        th = 0.0
    else:
        th = math.degrees(math.atan(math.sqrt(4 * cp - 1)))
    return p, psi_n, cp, th


def main():
    files = sorted(RESULTS_DIR.glob("lebensader_cusp_ibm_*.json"))
    if not files:
        print("No lebensader_cusp_ibm_*.json found.")
        return
    latest = files[-1]
    with open(latest, encoding='utf-8') as f:
        d = json.load(f)
    print(f"File: {latest.name}")
    print(f"Backend: {d.get('backend')}, Job: {d.get('job_id')}, Path: {d.get('path')}")
    print()

    sim_pred = d.get('simulator_predictions', {})
    hw_exp = d.get('hardware_expectations', {})

    print(f"{'H':<8s}  {'t':>4s}  "
          f"{'CΨ_HW':>7s}  {'θ_HW':>7s}  "
          f"{'CΨ_sim':>7s}  {'θ_sim':>7s}  "
          f"{'Δθ':>7s}  {'reading':<22s}")
    print('-' * 88)

    for h_label in ['XY+YX', 'YZ+ZY']:
        for t_str in ['0.8', '1.4', '1.7', '2.2']:
            t = float(t_str)
            sim = sim_pred.get(h_label, {}).get(t_str, {})
            exps = hw_exp.get(h_label, {}).get(t_str, {})
            if not exps:
                continue
            rho_hw = reconstruct_rho2q(exps)
            p_hw, psi_hw, cp_hw, th_hw = metrics(rho_hw)
            cp_sim = sim.get('cpsi', None)
            th_sim = sim.get('theta', None)
            if cp_sim is None or th_sim is None:
                continue
            dtheta = th_hw - th_sim
            reading = ""
            if cp_hw > 0.25 and cp_sim > 0.25:
                reading = "both above cusp"
            elif cp_hw <= 0.25 and cp_sim <= 0.25:
                reading = "both below cusp"
            elif cp_hw > 0.25 and cp_sim <= 0.25:
                reading = "HW above, sim below"
            else:
                reading = "HW below, sim above"
            print(f"{h_label:<8s}  {t:>4.1f}  "
                  f"{cp_hw:>7.4f}  {th_hw:>6.2f}°  "
                  f"{cp_sim:>7.4f}  {th_sim:>6.2f}°  "
                  f"{dtheta:>+6.2f}°  {reading:<22s}")
        print()

    print()
    print("Lebensader reading on the reduced (q0, q2) state:")
    print()
    for h_label in ['XY+YX', 'YZ+ZY']:
        # Extract θ_HW(t) trajectory
        thetas = []
        cpsis = []
        ts = []
        for t_str in ['0.8', '1.4', '1.7', '2.2']:
            t = float(t_str)
            exps = hw_exp.get(h_label, {}).get(t_str)
            if not exps:
                continue
            rho_hw = reconstruct_rho2q(exps)
            _, _, cp, th = metrics(rho_hw)
            ts.append(t)
            thetas.append(th)
            cpsis.append(cp)

        crosses = sum(1 for i in range(len(cpsis) - 1)
                      if (cpsis[i] - 0.25) * (cpsis[i + 1] - 0.25) < 0)
        n_above = sum(1 for c in cpsis if c > 0.25)
        n_below = sum(1 for c in cpsis if c <= 0.25)
        print(f"  {h_label:<8s}: θ trajectory = "
              f"{', '.join(f'{th:.1f}°' for th in thetas)}")
        print(f"           CΨ = "
              f"{', '.join(f'{c:.3f}' for c in cpsis)}")
        print(f"           crossings of CΨ=1/4: {crosses}; "
              f"above: {n_above}/{len(ts)}, below: {n_below}/{len(ts)}")
        print()

    # Also show the soft-break ⟨X₀Z₂⟩ signature
    print()
    print("Soft-break signature ⟨X₀Z₂⟩(t) on hardware vs simulator:")
    print()
    print(f"  {'H':<8s}  ", end='')
    for t_str in ['0.8', '1.4', '1.7', '2.2']:
        print(f"  t={t_str}: HW vs sim         ", end='')
    print()
    for h_label in ['XY+YX', 'YZ+ZY']:
        print(f"  {h_label:<8s}  ", end='')
        for t_str in ['0.8', '1.4', '1.7', '2.2']:
            xz_hw = hw_exp.get(h_label, {}).get(t_str, {}).get('X,Z', None)
            xz_sim = sim_pred.get(h_label, {}).get(t_str, {}).get('xz_full', None)
            if xz_hw is None or xz_sim is None:
                cell = "    —     "
            else:
                cell = f"{xz_hw:>+7.4f}/{xz_sim:>+7.4f}  "
            print(cell, end='')
        print()


if __name__ == "__main__":
    main()
