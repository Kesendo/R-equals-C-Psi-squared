"""Hardware-data diagnostic: apply the F-toolkit lens to measured 2-qubit Pauli expectations."""
from __future__ import annotations

from typing import Optional

import numpy as np

from ..diagnostics.f77_trichotomy import classify_pauli_pair
from ..diagnostics.f83_anti_fraction import predict_pi_decomposition
from .predict_signature_table import predict_signature_table


def diagnose_hardware(
    chain,
    measured,
    terms_per_category,
    calibration: Optional[dict] = None,
    predictions: Optional[dict] = None,
    t: float = 0.8,
    n_trotter: int = 3,
    gamma_z: Optional[float] = None,
    shots: int = 4096,
):
    """F-toolkit lens-reading workflow for hardware data.

    Operationalizes the lens-first discipline (memory feedback_framework_lens_first):
    given measured 2-qubit Pauli expectations per category and the corresponding
    Hamiltonians, walk through the F-chain F77 → F83 → F82/F84, label each
    structural feature of the residuals with the F-theorem it belongs to,
    and return a typed diagnostic.

    The diagnostic codifies the lens-readings developed in the 2026-04-30
    brainstorming arc:
      - F77 trichotomy per category (truly / soft / hard).
      - F83 anti-fraction prediction (closed form on H letters).
      - F82/F84 amplitude-damping signature: truly's ⟨Z,Z⟩ damping IS the
        operational F82 prediction (M=0 lens; environment shines through).
      - Pure-class quantitative match: pi2_odd_pure and pi2_even_nontruly
        should match Trotter prediction at small RMS; deviation indicates
        non-Z noise or path-dependent γ_Z_eff.
      - Y/Z asymmetry on truly: per-qubit T2 inhomogeneity reads as
        ⟨Y₀ Z_{N-1}⟩ ≠ ⟨Z₀ Y_{N-1}⟩ when Hamiltonian symmetry guarantees
        equality.

    Args:
        chain: ChainSystem providing N, J, bonds, gamma_0.
        measured: dict {category_name: {(p0, p2): expectation_value, ...}}.
        terms_per_category: dict {category_name: [(P, Q), ...]} — Hamiltonians.
        calibration: optional dict with 'T1', 'T2' lists per qubit (used to
            attribute Y/Z asymmetry to T2 inhomogeneity when present).
        predictions: optional dict in same format as `measured`. If None,
            compute Trotter+γ_Z predictions internally via
            predict_signature_table.
        t, n_trotter, gamma_z: passed to predict_signature_table.
        shots: per-Pauli-basis shot count, for σ-error estimate (1/√shots).

    Returns:
        dict with keys:
          'per_category':
              {category: {
                  'F77_class': str,
                  'F83_anti_fraction': float or None,
                  'F83_r': float or None,
                  'predictions': dict,
                  'measurements': dict,
                  'residuals': dict,
                  'rms_residual': float,
                  'lens_readings': list of dicts (typed F-theorem labels),
              }}
          'cross_category':
              {'Y_Z_asymmetry_on_truly': dict or None,
               'discrimination': dict (pairwise σ-units),
               'shot_noise_sigma': float}
    """
    N = chain.N
    if predictions is None:
        predictions = predict_signature_table(
            chain, terms_per_category,
            t=t, n_trotter=n_trotter, gamma_z=gamma_z,
        )

    sigma = 1.0 / np.sqrt(shots)

    per_category = {}
    for category, terms in terms_per_category.items():
        f77 = classify_pauli_pair(chain, terms)
        f83 = predict_pi_decomposition(chain, terms)

        meas = measured[category]
        pred = predictions[category]
        common_keys = [k for k in pred if k in meas]
        residuals = {k: meas[k] - pred[k] for k in common_keys}
        rms = float(np.sqrt(np.mean([r ** 2 for r in residuals.values()])))

        readings = []

        # F77 lens
        readings.append({
            'lens': 'F77',
            'class': f77,
            'reading': f"{category} is {f77} per F77 trichotomy classifier",
        })

        # F83 lens
        if f83['M_sq'] < 1e-12:
            readings.append({
                'lens': 'F83',
                'reading': "M = 0 idealized: avoidance defense (Mother lens). "
                           "Any HW residual is environmental signature.",
                'M_sq': float(f83['M_sq']),
            })
        else:
            r_val = f83['r']
            anti = f83['anti_fraction']
            r_str = f'{r_val:.3f}' if r_val != float('inf') else '∞'
            if abs(r_val) < 1e-10:
                role = 'pure Π²-odd (Father lens, active recirculation)'
            elif r_val == float('inf'):
                role = 'pure Π²-even non-truly (Child lens, passive reflection)'
            else:
                role = f'mixed (drive/echo ratio r={r_str})'
            readings.append({
                'lens': 'F83',
                'reading': f"anti-fraction = {anti:.4f}, r = {r_str}; {role}",
                'anti_fraction': float(anti),
                'r': float(r_val) if r_val != float('inf') else None,
            })

        # F82/F84 amplitude-damping signature on truly
        if f77 == 'truly':
            zz_pred = pred.get(('Z', 'Z'), 0.0)
            zz_meas = meas.get(('Z', 'Z'))
            if zz_meas is not None and abs(zz_pred) > 0.05:
                damping_frac = (zz_pred - zz_meas) / zz_pred
                # σ⁻ amplitude damping pulls Z toward +1, breaks ⟨Z⟩-conservation
                # ⟨Z,Z⟩ correlation drops if Z population redistributes via T1.
                significant = abs(damping_frac) > 0.20
                interpretation = (
                    "F82/F84 amplitude-damping signature: σ⁻ destroys |1⟩ population, "
                    "breaks ⟨Z⟩-conservation that the truly Hamiltonian alone protects. "
                    "Truly is the M=0 lens that lets this signature shine through cleanly."
                ) if significant else (
                    "No significant T1 contamination; truly's avoidance defense holds."
                )
                readings.append({
                    'lens': 'F82/F84',
                    'reading': interpretation,
                    'observable': '⟨Z,Z⟩',
                    'predicted': float(zz_pred),
                    'measured': float(zz_meas),
                    'damping_fraction': float(damping_frac),
                    'significant': bool(significant),
                })

        # F83 quantitative match for pure classes
        if f83['M_sq'] > 1e-12:
            r_val = f83['r']
            is_pure = abs(r_val) < 1e-10 or r_val == float('inf')
            if is_pure:
                # RMS in σ units
                rms_sigma = rms / sigma if sigma > 0 else float('inf')
                quality = (
                    "EXCELLENT match (within shot noise)" if rms_sigma < 5
                    else "GOOD match (small structural drift)" if rms_sigma < 15
                    else "DRIFT present (path-specific γ_Z_eff or ordering issue)"
                )
                readings.append({
                    'lens': 'F83-quantitative',
                    'reading': f"pure-class anti-fraction match: RMS={rms:.4f} ({rms_sigma:.1f}σ at {shots} shots) — {quality}",
                    'rms': float(rms),
                    'rms_sigma_units': float(rms_sigma),
                    'quality': quality,
                })

        per_category[category] = {
            'F77_class': f77,
            'F83_anti_fraction': float(f83['anti_fraction']) if f83['M_sq'] > 1e-12 else None,
            'F83_r': float(f83['r']) if (f83['M_sq'] > 1e-12 and f83['r'] != float('inf')) else None,
            'predictions': pred,
            'measurements': meas,
            'residuals': residuals,
            'rms_residual': rms,
            'lens_readings': readings,
        }

    # Cross-category: Y/Z asymmetry on truly
    cross_category = {'shot_noise_sigma': float(sigma)}

    truly_cats = [c for c, t_ in terms_per_category.items()
                  if classify_pauli_pair(chain, t_) == 'truly']
    if truly_cats:
        truly_cat = truly_cats[0]
        m = measured[truly_cat]
        yz = m.get(('Y', 'Z'))
        zy = m.get(('Z', 'Y'))
        if yz is not None and zy is not None:
            asymmetry = abs(yz - zy)
            t2 = (calibration or {}).get('T2')
            if t2 is not None and len(set(t2)) > 1:
                t2_min = min(t2)
                t2_max = max(t2)
                attribution = (
                    f"per-qubit T2 inhomogeneity (T2 range {t2_min:.0f}-{t2_max:.0f} μs); "
                    "Y measurement on shorter-T2 qubit damps more"
                )
                attributed = True
            else:
                attribution = "asymmetry detected; T2 calibration not provided to attribute"
                attributed = False
            cross_category['Y_Z_asymmetry_on_truly'] = {
                'category': truly_cat,
                'YZ': float(yz),
                'ZY': float(zy),
                'asymmetry_magnitude': float(asymmetry),
                'sigma_units': float(asymmetry / sigma) if sigma > 0 else float('inf'),
                'expected_symmetric': True,
                'attribution': attribution,
                'attributed_to_calibration': attributed,
            }
        else:
            cross_category['Y_Z_asymmetry_on_truly'] = None

    # Pairwise discrimination at top observables
    cats = list(terms_per_category.keys())
    discrimination = {}
    if len(cats) >= 2 and cats:
        all_observables = set()
        for c in cats:
            all_observables.update(measured[c].keys())
        # Top discriminators: largest pairwise spread in measured values
        spreads = []
        for obs in all_observables:
            vals = [measured[c].get(obs) for c in cats]
            vals_present = [v for v in vals if v is not None]
            if len(vals_present) < 2:
                continue
            spread = max(vals_present) - min(vals_present)
            spreads.append((obs, spread))
        spreads.sort(key=lambda x: -x[1])
        top_discriminators = []
        for obs, spread in spreads[:6]:
            top_discriminators.append({
                'pauli': obs,
                'spread': float(spread),
                'sigma_units': float(spread / sigma) if sigma > 0 else float('inf'),
                'values': {c: float(measured[c][obs]) for c in cats if obs in measured[c]},
            })
        discrimination['top_observables'] = top_discriminators
    cross_category['discrimination'] = discrimination

    return {
        'per_category': per_category,
        'cross_category': cross_category,
    }
