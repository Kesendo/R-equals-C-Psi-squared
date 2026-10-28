"""Hardware-data diagnostic: apply the F-toolkit lens to measured 2-qubit Pauli expectations.

Exposes raw structural properties of the Hamiltonian via PauliHamiltonian
(klein_set, y_parity_set, full Z₂³ signature, per-term properties) without
the Klein-4 letter labels or Trinity Mother/Father/Child collapse.
Application-layer interpretations live elsewhere; the diagnostic only
exposes structural facts.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

from ..diagnostics.f77_trichotomy import classify_pauli_pair
from ..diagnostics.f83_anti_fraction import predict_pi_decomposition
from ..diagnostics.polarity_coordinates import (
    polarity_coordinates,
    is_structurally_degenerate,
)
from ..pauli_hamiltonian import PauliHamiltonian, PauliTerm
from .predict_signature_table import predict_signature_table


def _structural_summary(terms, chain_length):
    """Build a PauliHamiltonian and extract its raw structural properties.

    Returns a dict with the Hamiltonian object plus its derived properties
    exposed as plain Python values for serialization. No labels, no Trinity
    collapse — just the structure itself.

    Three independent Z₂ axes per term:
      bit_a: (#X + #Y) mod 2 — Z⊗N parity break
      bit_b: (#Y + #Z) mod 2 — Π² parity
      Y-par: #Y mod 2          : fixed by Klein at one non-identity-count parity

    Klein-homogeneity (all terms share Klein index) forces soft-or-truly only
    among identity-free two-letter pairs, and even there it is sufficient, never
    necessary: plenty of Klein-INhomogeneous pairs are soft too. It does not
    generalise: identity-bearing two-letter pairs are 5/15 hard (e.g. IZ+ZI,
    Klein (0,1), hard at N=3 and N=4), and at k=3 the 294-pair Z₂³-homogeneous
    sweep finds 50 hard. The reading below therefore reports the Klein structure,
    not a verdict; call `classify_pauli_pair` for the class.
    """
    H = PauliHamiltonian.from_letter_tuples(terms, chain_length=chain_length)

    # Canonical ordering for display: lex-sorted tuples
    klein_set_sorted = sorted(H.klein_set)
    y_parity_set_sorted = sorted(H.y_parity_set)
    z2_signature_set_sorted = sorted(H.full_z2_signature_set)

    # Build a structural reading text
    if not H.terms:
        reading = "empty Hamiltonian"
    elif H.is_klein_homogeneous:
        klein_str = str(klein_set_sorted[0])
        reading = (
            f"Klein-homogeneous (all terms share Klein index {klein_str}); "
            f"Klein index alone does NOT settle the F87 class here: it forces "
            f"soft-or-truly only for identity-free two-letter pairs, and even "
            f"there only one way round; call classify_pauli_pair for the verdict."
        )
        if not H.is_z2_homogeneous:
            reading += (
                f" Klein-homogeneous but Z₂³-inhomogeneous: Y-parities differ "
                f"({y_parity_set_sorted}), full Z₂³ signatures differ "
                f"({z2_signature_set_sorted}); this already separates two-letter terms, "
                f"e.g. XY against ZI."
            )
    else:
        klein_str = ', '.join(str(k) for k in klein_set_sorted)
        reading = (
            f"Klein-inhomogeneous (terms span Klein indices {{{klein_str}}}); "
            f"F77 hard is possible but not guaranteed (depends on bond "
            f"combinatorics within the inhomogeneous span)."
        )

    return {
        'hamiltonian': H,
        'klein_set': klein_set_sorted,
        'is_klein_homogeneous': H.is_klein_homogeneous,
        'y_parity_set': y_parity_set_sorted,
        'is_y_parity_homogeneous': H.is_y_parity_homogeneous,
        'full_z2_signature_set': z2_signature_set_sorted,
        'is_z2_homogeneous': H.is_z2_homogeneous,
        'bit_b_set': sorted(H.bit_b_set),
        'bit_b_homogeneous_h': H.is_bit_b_homogeneous,
        'per_term_klein_indices': H.per_term_klein_indices,
        'per_term_y_parities': H.per_term_y_parities,
        'per_term_pi2_classes': H.per_term_pi2_classes,
        'per_term_full_z2_signatures': [t.full_z2_signature for t in H.terms],
        'k_body_set': sorted(H.k_body_set),
        'is_uniform_body': H.is_uniform_body,
        'has_truly_term': H.has_truly_term,
        'reading': reading,
    }


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
      - F87 trichotomy per category (truly / soft / hard).
      - F83 anti-fraction prediction (closed form on H letters).
      - F82/F84 window: for classes that commute with X^⊗N (truly XX+YY,
        also YZ+ZY) the single-site ⟨Z⟩, ⟨Y⟩ at the ends, held at zero by
        X^⊗N in any Z-only model, in shot-noise units; beyond 3σ a candidate
        for X^⊗N-breaking noise (σ⁻ among others), never an attribution.
      - Pure-class quantitative match: pi2_odd_pure and pi2_even_nontruly
        should match the Trotter prediction at small RMS.
      - Y/Z asymmetry on truly: compared with the PREDICTED pair, which the
        Trotter order already makes asymmetric; an excess beyond it over 3σ
        (σ = √2/√shots for a difference) names per-qubit T2 inhomogeneity as
        a candidate when the calibration shows it.

    The prediction's Trotter step is the hardware runners' (bonds outside,
    terms inside; gated against Qiskit's synthesis in
    tests/workflows/test_predict_signature_table.py).

    Args:
        chain: ChainSystem providing N, J, bonds, gamma_0.
        measured: dict {category_name: {(p0, p2): expectation_value, ...}}.
        terms_per_category: dict {category_name: [(P, Q), ...]} — Hamiltonians.
        calibration: optional dict with 'T1', 'T2' lists per qubit (used to
            name T2 inhomogeneity as a candidate for a Y/Z excess;
            'attributed_to_calibration' is True when such a candidate is named).
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
            'reading': f"{category} is {f77} per F87 trichotomy classifier",
        })

        # Structure lens: raw Z₂³ axes (bit_a, bit_b, Y-parity) per term and
        # aggregate homogeneity flags from PauliHamiltonian.
        struct = _structural_summary(terms, chain.N)
        readings.append({
            'lens': 'Structure',
            'reading': struct['reading'],
            'klein_set': struct['klein_set'],
            'is_klein_homogeneous': struct['is_klein_homogeneous'],
            'y_parity_set': struct['y_parity_set'],
            'is_y_parity_homogeneous': struct['is_y_parity_homogeneous'],
            'is_z2_homogeneous': struct['is_z2_homogeneous'],
            'per_term_klein_indices': struct['per_term_klein_indices'],
            'per_term_y_parities': struct['per_term_y_parities'],
        })

        # F83 lens
        if f83['M_sq'] < 1e-12:
            readings.append({
                'lens': 'F83',
                'reading': "M = 0 idealized: avoidance defense (Mother lens). "
                           "A residual against a correctly modelled circuit is "
                           "where the environment would show.",
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

        # F112 lens: polarity Π-eigenvalue +i/−i balance on the framework-predicted L.
        # Typed Tier1Derived scope: Hermitian H + bit_b-homogeneous c_k implies
        # asymmetry = 0 exactly AS A THEOREM (the computed value is a float difference
        # and is a different question; see the branch below). chain.L's single-Pauli Z-dephasing is trivially
        # bit_b-homogeneous, so for pure Z-deph the lens predicts BALANCED. Adding
        # T1 / T2 dissipators leaves F112's typed scope (σ⁻ = (X + iY)/2 is not
        # bit_b-homogeneous), and there the balance is NOT guaranteed: F113 gives
        # the exact gap once H carries a single-site Z moment.
        # Note this call cannot see that, but for the dissipator's sake, not the
        # Hamiltonian's: it passes no gamma_t1, so F112's typed scope holds and
        # BALANCED is forced whatever the terms carry. (A Z moment IS reachable
        # through this entry point, as a bond term with an identity leg; it just
        # needs an amplitude-damping channel beside it to show.) Treat it as a wiring check on the terms given, not
        # as evidence about the polarity axis; the L-bound entry point can break.
        pol = polarity_coordinates(chain, terms)
        m_sq = pol['norm_sq']['M']  # reported below; not the denominator any more
        asym = pol['asymmetry']
        # The denominator is the POLARITY CONTENT ||M_anti||^2 = ||M_+||^2 + ||M_-||^2, the
        # whole of what asym is a difference of, so rel_asym is a contrast ratio in [0, 1]
        # with no tolerance constant in it. It replaces max(||M||^2, 1e-15), which was
        # saturated: for the Heisenberg + pure Z-dephasing case this workflow runs most, M is
        # the F81 residual and 'truly' H makes it vanish by the F83 closed form, so ||M||^2 is
        # exactly 0.0 at N=2 and ~1e-31 at N=3..5 (do NOT write 'exactly 0.0' unqualified; an
        # earlier draft of this comment did and contradicted OpenArcsRegistry's own 1e-31).
        # Either way the 1e-15 FLOOR was the divisor and the 1e-10 cut admitted any |asym|
        # below 1e-25. The two denominators convert by F83's anti-fraction, new = old*(2+4r)
        # with r = ||H_even_nontruly||^2/||H_odd||^2. Mirrors
        # PolarityCoordinatesResult.RelativeAsymmetry in C#; the same denominator has been in
        # simulations/polarity_step5_stress.py line 86 since it was written.
        m_anti_sq = pol['norm_sq']['M_plus_half'] + pol['norm_sq']['M_minus_half']
        # STRUCTURAL first, float second. The exact test reads the Pauli letters and holds
        # at every N; the float `== 0.0` fallback only fires at N = 2, because the
        # structurally-zero M_anti arrives as ~1e-33 at N=3 and ~1e-31 at N=5. Keeping only
        # the float test would let a silent row report a ratio of two noise quantities.
        # This entry point passes no amplitude damping (the call above is
        # polarity_coordinates(chain, terms)), so c is pure Z-dephasing and always
        # bit_b-homogeneous; the H parity alone decides.
        structurally_silent = is_structurally_degenerate(terms)
        degenerate = structurally_silent or (m_anti_sq == 0.0)
        # NaN, not 0.0, and for the same reason as in polarity_fingerprint: where there is
        # no polarity content the quotient is 0/0, and a 0.0 in a contrast-ratio field reads
        # as "perfectly balanced", the one reading DEGENERATE exists to refuse. NaN says
        # "no number" and cannot be silently aggregated away (it loses every comparison, so
        # a max() or a mean() over mixed rows shows it instead of hiding it).
        rel_asym = math.nan if degenerate else abs(asym) / m_anti_sq
        if degenerate:
            # A third verdict word, matching the C# witnesses, which gained DEGENERATE on
            # 2026-08-07. The earlier comment here said the opposite and justified staying
            # two-valued by pointing at C#; that stopped being true in the same change.
            # It matters because the stored 'verdict' below is verdict.split(' ')[0], so a
            # two-valued word made a silent row report 'BALANCED' with the note truncated off.
            # Decided BEFORE the thresholds, not patched on after: rel_asym is NaN here and
            # NaN fails every comparison, so running the ladder first would route a silent
            # row through the BROKEN branch to reach the same word by the wrong road.
            # The REASON is stated only when the structural test is what fired. Reached
            # through the float arm instead, all that is known is the measured zero, and
            # naming a structural cause there would assert what the flag does not carry.
            why = ('H is Pi^2-even and c is bit_b-homogeneous, so both halves vanish as a '
                   'theorem' if structurally_silent else
                   'both halves measured exactly 0.0, though the structural test did not fire')
            verdict = ('DEGENERATE (no ratio) [vacuous: no polarity content, ||M_anti||^2 = '
                       f'{m_anti_sq:.2e}; {why}, and there is nothing for a threshold '
                       'to decide]')
        # Report the measured number, never the word. F112's THEOREM is exact, but the ladder
        # below is a threshold on a float difference of two Frobenius norms and cannot witness
        # exactness; GLOSSARY.md:304 is the governing rule. Measured 2026-08-06, and what the
        # old 'bit-exact 0' literal hid: the asymmetry is exactly 0.0 on the bilinear-Pauli
        # Hamiltonians the repo tests, but NOT in general inside F112's typed scope: random
        # Hermitian H with per-site Z c-ops give a nonzero asymmetry in a few percent of draws
        # (the COUNT is seed noise and must not be quoted; the ~1e-17 RELATIVE size is not).
        # Exactly-zero here is a property of those inputs, not of the float route. Recorded in
        # the bit_exact_vocabulary arc.
        elif rel_asym < 1e-10:
            verdict = f'BALANCED (rel asymmetry {rel_asym:.2e})'
        elif rel_asym < 1e-6:
            verdict = f'near-BALANCED (rel {rel_asym:.2e})'
        else:
            verdict = f'BROKEN (rel asymmetry {rel_asym:.2e})'
        bit_b_homog = struct['bit_b_homogeneous_h']
        scope_note = (
            'H is bit_b-homogeneous (terms share bit_b parity)'
            if bit_b_homog
            else 'H is bit_b-inhomogeneous (terms span both bit_b parities)'
        )
        readings.append({
            'lens': 'F112',
            'reading': f"polarity asymmetry = {verdict}; {scope_note}",
            'asymmetry': float(asym),
            'rel_asymmetry': float(rel_asym),
            'M_norm_sq': float(m_sq),
            'M_anti_norm_sq': float(m_anti_sq),
            # Read WITH 'verdict'. Since 2026-08-07 'rel_asymmetry' also tells the two apart
            # on its own (a vacuous row carries NaN, a genuine balance a real 0.0), but this
            # flag is the one to branch on: it is a bool, and NaN needs math.isnan.
            'polarity_degenerate': bool(degenerate),
            'verdict': verdict.split(' ', 1)[0],
            'h_bit_b_homogeneous': bool(bit_b_homog),
        })

        # F82/F84 window: the X^⊗N structural zeros. A Hamiltonian whose every
        # term carries an even number of Y/Z letters (truly XX+YY, also YZ+ZY)
        # commutes with X^⊗N; Z-dephasing and the X-Néel state respect it too, so
        # a Z-only model predicts the single-site ⟨Z_i⟩, ⟨Y_i⟩ (which anticommute
        # with X^⊗N) to be exactly 0. σ⁻ breaks X^⊗N and pulls ⟨Z_i⟩ positive.
        # Two-body ⟨Z₀Z₂⟩ commutes with X^⊗N and is nearly blind to σ⁻ on truly
        # (0.231 → 0.221 at γ_T1 = 0.046), so it is not the window. Error model:
        # shot noise only (σ = 1/√shots bounds one expectation); beyond 3σ marks a
        # candidate for anything that breaks X^⊗N (σ⁻, coherent single-site
        # errors, leakage), never an attribution.
        if all(sum(letter in 'YZ' for letter in term) % 2 == 0 for term in terms):
            sites = []
            for key, label in ((('Z', 'I'), '⟨Z₀⟩'), (('I', 'Z'), f'⟨Z_{N - 1}⟩'),
                               (('Y', 'I'), '⟨Y₀⟩'), (('I', 'Y'), f'⟨Y_{N - 1}⟩')):
                if key in meas and key in pred and abs(pred[key]) < 1e-12:
                    dev = meas[key] - pred[key]
                    sites.append({'observable': label, 'predicted': float(pred[key]),
                                  'measured': float(meas[key]),
                                  'sigma_units': float(abs(dev) / sigma)})
            if sites:
                worst = max(sites, key=lambda r: r['sigma_units'])
                significant = worst['sigma_units'] > 3.0
                interpretation = (
                    f"{worst['observable']} reads {worst['sigma_units']:.1f}σ from the X^⊗N "
                    "structural zero: a candidate for noise that breaks X^⊗N (F82/F84 σ⁻, "
                    "coherent single-site errors, leakage), to be checked against a "
                    "T1-carrying prediction before it is read as any one of them."
                ) if significant else (
                    f"X^⊗N structural zeros within {worst['sigma_units']:.1f}σ; "
                    "no X^⊗N-breaking signature above 3σ."
                )
                readings.append({
                    'lens': 'F82/F84',
                    'reading': interpretation,
                    'observables': sites,
                    'sigma_units': worst['sigma_units'],
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
                    "EXCELLENT match (RMS under 5σ)" if rms_sigma < 5
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
            'structure': struct,
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
            # The predicted pair carries the circuit's own asymmetry: the Trotter
            # order (bonds outside) breaks the q0 ↔ q_{N-1} reflection at O(dt²).
            # Only the excess over it is a residual. The excess is a difference of
            # two independently measured expectations, so its shot-noise bound is
            # √2/√shots. A significant excess names T2 inhomogeneity as a
            # candidate, never as an attribution: the prediction holds no T2.
            p = predictions[truly_cat]
            pred_yz = p.get(('Y', 'Z'))
            pred_zy = p.get(('Z', 'Y'))
            have_pred = pred_yz is not None and pred_zy is not None
            # A reflection-symmetric circuit gives the pair equal up to the roundoff
            # of one 2^N propagation, O(100·eps); the circuit order gives O(0.01-0.5).
            expected_symmetric = (abs(pred_yz - pred_zy) < 1e-12) if have_pred else None
            excess = ((yz - zy) - (pred_yz - pred_zy)) if have_pred else None
            sigma_excess = np.sqrt(2.0) * sigma
            excess_sigma_units = abs(excess) / sigma_excess if have_pred else None
            excess_significant = have_pred and excess_sigma_units > 3.0
            t2 = (calibration or {}).get('T2')
            if not have_pred:
                attribution = "prediction lacks the ⟨Y,Z⟩/⟨Z,Y⟩ pair; asymmetry not compared"
                attributed = False
            elif not excess_significant:
                attribution = (
                    f"measured asymmetry matches the predicted pair "
                    f"(excess {excess:+.3f}, {excess_sigma_units:.1f}σ)"
                )
                attributed = False
            elif t2 is not None and len(set(t2)) > 1:
                attribution = (
                    f"excess {excess:+.3f} ({excess_sigma_units:.1f}σ) beyond the predicted pair; "
                    f"candidate: per-qubit T2 inhomogeneity (T2 range {min(t2):.0f}-{max(t2):.0f} μs), "
                    "to be checked against a T2-carrying prediction"
                )
                attributed = True
            elif t2 is not None:
                attribution = (
                    f"excess {excess:+.3f} ({excess_sigma_units:.1f}σ) beyond the predicted pair; "
                    "T2 is uniform on the path, so not a T2 difference"
                )
                attributed = False
            else:
                attribution = (
                    f"excess {excess:+.3f} ({excess_sigma_units:.1f}σ) beyond the predicted pair; "
                    "no T2 calibration provided"
                )
                attributed = False
            cross_category['Y_Z_asymmetry_on_truly'] = {
                'category': truly_cat,
                'YZ': float(yz),
                'ZY': float(zy),
                'asymmetry_magnitude': float(asymmetry),
                'sigma_units': float(asymmetry / sigma) if sigma > 0 else float('inf'),
                'predicted_YZ': None if pred_yz is None else float(pred_yz),
                'predicted_ZY': None if pred_zy is None else float(pred_zy),
                'excess_asymmetry': None if excess is None else float(excess),
                'excess_sigma_units': None if excess_sigma_units is None else float(excess_sigma_units),
                'excess_significant': bool(excess_significant),
                'expected_symmetric': expected_symmetric if expected_symmetric is None else bool(expected_symmetric),
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
