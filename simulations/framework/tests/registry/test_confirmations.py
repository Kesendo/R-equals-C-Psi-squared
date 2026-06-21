"""Tests for the Confirmations registry of hardware-confirmed predictions."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_confirmations_has_twenty_entries():
    names = fw.Confirmations.list_names()
    assert len(names) == 20
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


def test_confirmations_lookup_ibm_ep_onset():
    # The Kingston EP-onset run: the revival lifts off the 1/N floor as Q crosses
    # Q_EP ≈ 1.5. Registry anchor for the hard-coded hardware table in
    # compute/RCPsiSquared.Diagnostics/Foundation/EpField.cs (node 5, "the hardware").
    e = fw.Confirmations.lookup('ibm_ep_onset_may2026')
    assert e['date'] == '2026-05-31'
    assert e['machine'] == 'ibm_kingston'
    assert 'd8dr7dfd0j8c73f4man0' in e['job_id']
    assert 'd8drjbfd0j8c73f4mobg' in e['job_id']
    assert e['measured_value']['Q_grid'] == [0.5, 1.0, 1.5, 2.5, 5.0, 20.0]
    assert e['measured_value']['revival'] == [0.30, 0.36, 0.34, 0.49, 0.56, 0.70]
    assert '1/N = 1/3 equipartition floor' in e['predicted_value']['se_walk_handover_reading']
    assert 'ExceptionalPointClock' in e['framework_primitive']
    assert e['experiment_doc'] == 'experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md'


def test_confirmations_unknown_raises():
    with pytest.raises(KeyError):
        fw.Confirmations.lookup('does_not_exist')


def test_confirmations_by_machine():
    marrakesh = fw.Confirmations.by_machine('ibm_marrakesh')
    kingston = fw.Confirmations.by_machine('ibm_kingston')
    assert len(marrakesh) >= 6
    assert len(kingston) >= 3
