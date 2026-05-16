"""Tests for the Confirmations registry of hardware-confirmed predictions."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_confirmations_has_thirteen_entries():
    names = fw.Confirmations.list_names()
    assert len(names) == 13
    assert 'palindrome_trichotomy' in names
    assert 'lebensader_skeleton_trace_decoupling' in names
    assert 'gamma_0_marrakesh_calibration' in names
    assert 'marrakesh_transverse_y_field_detection' in names
    assert 'f83_pi2_class_signature_marrakesh' in names
    assert 'd_zero_sector_trichotomy_marrakesh' in names
    assert 'block_cpsi_saturation_kingston_may2026' in names
    assert 'f95_angle_steering_kingston_may2026' in names


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


def test_confirmations_unknown_raises():
    with pytest.raises(KeyError):
        fw.Confirmations.lookup('does_not_exist')


def test_confirmations_by_machine():
    marrakesh = fw.Confirmations.by_machine('ibm_marrakesh')
    kingston = fw.Confirmations.by_machine('ibm_kingston')
    assert len(marrakesh) >= 6
    assert len(kingston) >= 3
