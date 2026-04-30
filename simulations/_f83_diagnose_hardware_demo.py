"""F-toolkit lens-reading demo: re-analyze the F83 Marrakesh hardware run.

Codifies in one cockpit call what previously took four hand-rolled scripts
(_f83_signature_predictions, _f83_aer_preflight, _f83_hy_field_check,
_f83_gamma_z_sweep) plus manual analysis. The hypothesis from
THE_POLARITY_LAYER.md says the F-chain F77 → F80 → F81 → F82 → F83 → F84
should be readable as one pass through the hardware data; this demo runs
that pass via the diagnose_hardware workflow.

Reproduces the lens-readings documented in
data/ibm_f83_signature_april2026/README.md without any manual analysis:
  - F77 trichotomy classification per category
  - F83 anti-fraction (closed form, structural)
  - F82/F84 amplitude-damping signature in truly's ⟨Z,Z⟩ damping
  - Pure-class quantitative match (RMS in σ-units)
  - Y/Z asymmetry on truly attributed to per-qubit T2 inhomogeneity
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw


REPO_ROOT = Path(__file__).resolve().parents[1]
F83_HARDWARE_JSON = REPO_ROOT / 'data' / 'ibm_f83_signature_april2026' / 'f83_signature_ibm_marrakesh_20260430_190035.json'

F83_TERMS_PER_CATEGORY = {
    'truly_unbroken':       [('X', 'X'), ('Y', 'Y')],
    'pi2_odd_pure':         [('X', 'Y'), ('Y', 'X')],
    'pi2_even_nontruly':    [('Y', 'Z'), ('Z', 'Y')],
    'mixed_anti_one_sixth': [('X', 'Y'), ('Y', 'Z')],
}

# Calibration values from 2026-04-30T16:25Z, ibm_marrakesh path [4, 5, 6]
F83_CALIBRATION = {
    'T1': [206.115, 144.477, 351.240],
    'T2': [183.978, 120.938, 151.412],
}


def load_measured(json_path):
    with open(json_path) as f:
        data = json.load(f)
    measured = {}
    for cat, exp_dict in data['expectations'].items():
        cat_dict = {}
        for key, val in exp_dict.items():
            if ',' in key:
                p0, p2 = key.split(',')
                cat_dict[(p0, p2)] = val
        measured[cat] = cat_dict
    return measured, data.get('backend'), data.get('path'), data.get('job_id')


def main():
    measured, backend, path, job_id = load_measured(F83_HARDWARE_JSON)

    print("=" * 78)
    print("F83 hardware: F-toolkit lens-reading via diagnose_hardware workflow")
    print("=" * 78)
    print(f"  Backend: {backend}")
    print(f"  Path:    {path}")
    print(f"  Job ID:  {job_id}")
    print(f"  Calibration: T1={F83_CALIBRATION['T1']} μs, T2={F83_CALIBRATION['T2']} μs")
    print()

    chain = fw.ChainSystem(N=3)
    result = fw.diagnose_hardware(
        chain, measured, F83_TERMS_PER_CATEGORY,
        calibration=F83_CALIBRATION,
        t=0.8, n_trotter=3, gamma_z=0.1,
        shots=4096,
    )

    # Per-category lens-readings
    for cat in F83_TERMS_PER_CATEGORY:
        pc = result['per_category'][cat]
        print(f"--- {cat} ---")
        print(f"  F77 class:        {pc['F77_class']}")
        if pc['F83_anti_fraction'] is not None:
            r_str = f"{pc['F83_r']:.3f}" if pc['F83_r'] is not None else "∞"
            print(f"  F83 anti-fraction: {pc['F83_anti_fraction']:.4f}  (r = {r_str})")
        else:
            print(f"  F83 anti-fraction: n/a (M = 0 idealized)")
        print(f"  RMS residual:     {pc['rms_residual']:.4f}")
        for r in pc['lens_readings']:
            lens = r.get('lens', '?')
            reading = r.get('reading', '')
            print(f"    [{lens:>8}] {reading}")
        print()

    # Cross-category readings
    print("--- Cross-category ---")
    cc = result['cross_category']
    print(f"  Shot-noise σ (per ⟨P⟩ at 4096 shots): ±{cc['shot_noise_sigma']:.4f}")
    yz = cc.get('Y_Z_asymmetry_on_truly')
    if yz is not None:
        print(f"  Y/Z asymmetry on truly:")
        print(f"    ⟨Y₀ Z₂⟩ = {yz['YZ']:+.4f},  ⟨Z₀ Y₂⟩ = {yz['ZY']:+.4f}")
        print(f"    |asymmetry| = {yz['asymmetry_magnitude']:.4f}  ({yz['sigma_units']:.1f}σ)")
        print(f"    attribution: {yz['attribution']}")
    print()
    print("  Top discriminating Pauli observables (largest pairwise spreads):")
    for d in cc['discrimination']['top_observables']:
        p0, p2 = d['pauli']
        print(f"    {p0},{p2}: spread = {d['spread']:.4f}  ({d['sigma_units']:.0f}σ)")
    print()

    # Summary
    print("=" * 78)
    print("Hypothesis test summary")
    print("=" * 78)
    truly_f82 = next((r for r in result['per_category']['truly_unbroken']['lens_readings']
                      if r.get('lens') == 'F82/F84'), None)
    pi2_odd_quant = next((r for r in result['per_category']['pi2_odd_pure']['lens_readings']
                          if r.get('lens') == 'F83-quantitative'), None)
    pi2_even_quant = next((r for r in result['per_category']['pi2_even_nontruly']['lens_readings']
                           if r.get('lens') == 'F83-quantitative'), None)

    print(f"  F82/F84 truly signature detected: "
          f"{'YES' if truly_f82 and truly_f82.get('significant') else 'no'}"
          f" (⟨Z,Z⟩ damping {truly_f82['damping_fraction']*100:.1f}%)" if truly_f82 else "")
    if pi2_odd_quant:
        print(f"  pi2_odd_pure F83 match:          {pi2_odd_quant['quality']}")
    if pi2_even_quant:
        print(f"  pi2_even_nontruly F83 match:     {pi2_even_quant['quality']}")
    if yz is not None:
        print(f"  Y/Z asymmetry on truly attributed to T2: "
              f"{'YES' if yz['attributed_to_calibration'] else 'no calibration provided'}")
    print()
    print("This output codifies in one cockpit call the manual analysis")
    print("documented in data/ibm_f83_signature_april2026/README.md.")


if __name__ == "__main__":
    main()
