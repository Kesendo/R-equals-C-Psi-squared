"""Tests for the Confirmations registry of hardware-confirmed predictions."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_price_pair_framework_primitive_names_its_rate_book():
    """The total hardware coherence rate cannot become F1's pure-Z price.

    Local T1 contributes to the former; F1's exact spectral identity has a
    pure-Z domain. Mirrored in the C# registry test.
    """
    primitive = fw.Confirmations.lookup(
        'price_pair_locality_marrakesh_july2026')['framework_primitive']
    assert 'pure-Z structural context' in primitive
    assert '2·Σγ_Z,j' in primitive
    assert 'total fitted transverse rates' in primitive
    assert 'including local T1' in primitive
    assert 'not an exact F1 center measurement' in primitive
    for wrong in ('γ_F1,j = 1/(2·T2*_j)', 'P = 2·Σγ_F1'):
        assert wrong not in primitive


def test_confirmations_has_twenty_four_entries():
    names = fw.Confirmations.list_names()
    assert len(names) == 24
    assert 'palindrome_trichotomy' in names
    assert 'lebensader_skeleton_trace_decoupling' in names
    assert 'gamma_0_marrakesh_calibration' in names
    assert 'marrakesh_transverse_y_field_detection' in names
    assert 'f83_pi2_class_signature_marrakesh' in names
    assert 'd_zero_sector_trichotomy_marrakesh' in names
    assert 'block_cpsi_saturation_kingston_may2026' in names
    assert 'f95_angle_steering_kingston_may2026' in names
    # 2026-06-08 reconciliation with the C# ConfirmationsRegistry (both now hold the
    # union; these two were previously C#-only).
    assert 'regime_uniformity_kingston_uniform_quantum' in names
    assert 'gamma0_off_the_lever_kingston_may2026' in names
    # 2026-06-10: Kingston EP-onset run added to BOTH registries (union of 16).
    assert 'ibm_ep_onset_may2026' in names
    # 2026-06-11: F120 moment-tower pump channel, first hardware reading (union of 17).
    assert 'f120_moment_tower_kingston_june2026' in names
    # 2026-06-18: the three Torino calibration-era runs (Feb-Mar 2026) registered at last
    # (front_matter_truth arc); they predate the systematic April-June campaign (union of 20).
    assert 'cpsi_quarter_crossing_torino_feb2026' in names
    assert 'absorption_theorem_ratio_torino' in names
    assert 'cpsi_quarter_crossing_torino_q80_mar2026' in names
    # 2026-07-05: the two-leg cold-bath attribution on the f95/F113 qubits (union of 22).
    assert 'f84_heating_leg_attribution_kingston_july2026' in names
    # 2026-07-12: the concentrator reload's site-resolved contrast, A-sign only per the
    # pre-registered verdict split (union of 23).
    assert 'concentrator_site_contrast_kingston_july2026' in names


def test_confirmations_lookup_concentrator_site_contrast():
    e = fw.Confirmations.lookup('concentrator_site_contrast_kingston_july2026')
    assert e['date'] == '2026-07-11'
    assert e['machine'] == 'ibm_kingston'
    assert e['job_id'] == 'd99a970tcv6s73dn2atg'
    assert 'slope(MP) − slope(E) = −0.05337/step' in e['measured_value']
    assert 'A-sign CONFIRMED' in e['measured_value']
    # the honest band story is part of the entry, not a hidden footnote
    assert 'magnitude off-prediction' in e['measured_value']
    assert 'A-SIGN ONLY' in e['predicted_value']
    assert '−0.07325' in e['predicted_value']          # the 7a dressed central prediction
    assert 'shot-only' in e['description']             # the binding-pooling instrument deviation
    assert 'concentrator_reloaded_7a.py' in e['framework_primitive']
    assert e['experiment_doc'] == 'experiments/IBM_CONCENTRATOR_RELOADED.md'


def test_confirmations_lookup_f84_heating_leg():
    e = fw.Confirmations.lookup('f84_heating_leg_attribution_kingston_july2026')
    assert e['date'] == '2026-07-05'
    assert e['machine'] == 'ibm_kingston + ibm_marrakesh'
    assert 'd951mhkql68s73ca3u0g' in e['job_id']       # run 1, Kingston [82,83,13]
    assert 'd953ti5gc6cc73ffomig' in e['job_id']       # run 2, Marrakesh [93,94,95]
    assert 'SPLIT on all three qubits' in e['measured_value']
    assert 'the bath is COLD' in e['measured_value']
    assert 'V_σ± = 0.02491' in e['measured_value']     # run 2's first valid corrected-recipe V
    assert 'run_heating_leg.py' in e['framework_primitive']
    assert e['experiment_doc'] == 'experiments/F81_VIOLATION_HARDWARE_BRIDGE.md'


def test_confirmations_lookup_f120_moment_tower():
    e = fw.Confirmations.lookup('f120_moment_tower_kingston_june2026')
    assert e['date'] == '2026-06-11'
    assert e['machine'] == 'ibm_kingston'
    assert 'd8l6c7rqv2lc73863acg' in e['job_id']
    assert 'd8l6c832d42s73cb16a0' in e['job_id']
    assert 'd8l6h03nn5bs738rmrug' in e['job_id']
    assert e['measured_value']['per_qubit_pump_slopes_per_us']['q9'] == [5.794e-3, 5.779e-3]
    assert e['measured_value']['arbiter_T1_us']['q13'] == 430.3
    assert 'HOLDS everywhere in-situ' in e['measured_value']['model_test_in_situ']
    assert 'EPOCH ARTIFACT' in e['measured_value']['corrected_same_day']
    assert 'f120_prep_split_reanalysis' in e['measured_value']['corrected_same_day']
    assert 'moment_tower' in e['framework_primitive']
    assert e['experiment_doc'] == 'experiments/F120_MOMENT_TOWER_KINGSTON.md'


def test_confirmations_lookup_f95_angle_steering():
    e = fw.Confirmations.lookup('f95_angle_steering_kingston_may2026')
    assert e['date'] == '2026-05-16'
    assert e['machine'] == 'ibm_kingston'
    assert 'bxyj5yd4j' in e['job_id']
    assert 'bzklqwt7f' in e['job_id']
    assert e['measured_value']['pair_A_mid_omega_0.13']['t_cross_us'] == 1.395
    assert e['measured_value']['pair_A_mid_omega_0.25']['t_cross_us'] == 1.242
    assert e['measured_value']['pair_B_high_omega_0.25']['t_cross_us'] == 2.814
    assert 'F95AngleAtQuadraticZeroPi2Inheritance' in e['framework_primitive']


def test_confirmations_lookup_palindrome_trichotomy():
    e = fw.Confirmations.lookup('palindrome_trichotomy')
    assert e['date'] == '2026-04-26'
    assert e['machine'] == 'ibm_marrakesh'
    assert e['job_id'] == 'd7mjnjjaq2pc73a1pk4g'
    assert e['measured_value']['delta_soft_minus_truly'] == -0.722


def _slowest_se_block_mode(q, n=3):
    """Slowest non-kernel mode of the single-excitation (1,1) block: -iq[A, rho] - 4(rho - diag rho)."""
    import numpy as np
    a = np.diag(np.ones(n - 1), 1)
    a = a + a.T
    ident = np.eye(n)
    lmat = -1j * q * (np.kron(a, ident) - np.kron(ident, a.T))
    lmat += np.diag(-4.0 * (np.ones((n, n)) - ident).flatten())
    w = np.linalg.eigvals(lmat)
    w = w[np.abs(w) > 1e-7]
    top = w.real.max()
    ties = w[np.abs(w.real - top) < 1e-7]
    return ties[np.argmax(np.abs(ties.imag))]


def test_confirmations_lookup_ibm_ep_onset():
    # The historical slug is retained, but the result is a population handover in
    # the runner's coherence-rate book, not a critical-damping or EP certificate.
    e = fw.Confirmations.lookup('ibm_ep_onset_may2026')
    assert e['date'] == '2026-05-31'
    assert e['machine'] == 'ibm_kingston'
    assert 'd8dr7dfd0j8c73f4man0' in e['job_id']
    assert 'd8drjbfd0j8c73f4mobg' in e['job_id']
    assert e['measured_value']['Q_label_grid'] == [0.5, 1.0, 1.5, 2.5, 5.0, 20.0]
    assert e['measured_value']['Q_lindblad_grid'] == [1.0, 2.0, 3.0, 5.0, 10.0, 40.0]
    assert e['measured_value']['revival'] == [0.2978515625, 0.3623809814453125,
                                              0.34356689453125, 0.4898834228515625,
                                              0.560699462890625, 0.70330810546875]
    assert '1/N = 1/3 reference level' in e['predicted_value']['se_walk_handover_reading']
    assert '0.28 → 0.84' in e['predicted_value']['twirl_simulate_revival']
    assert '0.8417853730254796' in e['measured_value']['high_Q_discrepancy']
    assert 'Q_Lindblad = 2 Q_label' in e['measured_value']['rate_book']
    # The walk's spectral transition is placed, not left open: the (1,1)-block EP Q*(3) = √2
    # (Q_label ≈ 0.71) lies below the handover, and the handover moves with the probe time.
    assert 'Q*(3) = √2' in e['description']
    assert 'Q_label ≈ 0.71' in e['description']
    assert 'probe-time crossover' in e['description']
    assert 'CoherenceHorizonClaim' in e['framework_primitive']
    assert 'EpCharacterWitness' in e['framework_primitive']
    assert 'K=16 twirl simulation' in e['measured_value']['high_Q_discrepancy']
    assert 'spectral character remains open' not in e['description']
    assert 'exact-twirl' not in e['measured_value']['high_Q_discrepancy']
    # The placement, computed on the walk's own (1,1) block from the registry's own grid: the slowest
    # non-kernel mode is real at the first flown Q_Lindblad and a rotating pair at the second, so the
    # EP lies between them. Both points are far from Q* = sqrt 2 (well-conditioned eigenvalues: |Im|
    # ~1e-15 on the real mode, exactly 2 on the pair at Q = 2), so the 1e-7 cut has decades of room.
    q0, q1 = e['measured_value']['Q_lindblad_grid'][:2]
    assert abs(_slowest_se_block_mode(q0).imag) < 1e-7
    assert abs(_slowest_se_block_mode(q1).imag) > 1.0
    assert 'do not certify convergence' in e['description']
    assert 'no calibrated error model' in e['description']
    assert 'floor and the onset are clean' not in e['description']
    assert 'critical-damping' not in e['predicted_value']['se_walk_handover_reading']
    assert 'ExceptionalPointClock' not in e['framework_primitive']
    assert e['experiment_doc'] == 'experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md'


def test_gamma0_lever_entry_is_a_finite_grid_bracket_not_a_carrier_calibration():
    e = fw.Confirmations.lookup('gamma0_off_the_lever_kingston_may2026')
    assert 'finite-time transfer maximum' in e['predicted_value']
    assert 'brackets an effective response change' in e['measured_value']
    assert 'does not isolate' in e['description']
    assert 'Confirms the typed' not in e['description']
    assert 'hardware-anchored' not in e['description']


def test_confirmations_unknown_raises():
    with pytest.raises(KeyError):
        fw.Confirmations.lookup('does_not_exist')


def test_confirmations_by_machine():
    marrakesh = fw.Confirmations.by_machine('ibm_marrakesh')
    kingston = fw.Confirmations.by_machine('ibm_kingston')
    assert len(marrakesh) >= 6
    assert len(kingston) >= 4
