"""Production gates must reject damaged surveys before publishing an artifact."""
import copy
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path
import sys

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'simulations'))

import numpy as np
import pytest
import route_b_other_n_unfolding as unfolding


def required_api(name):
    function = getattr(unfolding, name, None)
    assert callable(function), f'missing production API: {name}'
    return function


def inventory():
    return json.loads((ROOT/'simulations/results/route_b_a2_n5.json').read_text())


@lru_cache(maxsize=1)
def fresh_survey():
    return required_api('build_coefficient_survey')(inventory())


def survey():
    return copy.deepcopy(fresh_survey())


@lru_cache(maxsize=1)
def fresh_continuations():
    selected = {'N4-real-lambda-plus', 'N5-O-A2-W-006-Q-plus'}
    return {row['id']: unfolding.follow(row['n'], unfolding.z(row['q']),
            unfolding.z(row['lambda0']), row)
            for row in fresh_survey() if row['id'] in selected}


def test_survey_uses_every_inventory_locus_once():
    rows = survey()
    required_api('validate_coefficient_survey')(rows, inventory())
    assert len(rows) == 60
    assert [row['n'] for row in rows].count(4) == 2
    assert [row['n'] for row in rows].count(5) == 58
    expected = {loc['id'] for sector in inventory()['sectors']
                for root in sector['a2Roots'] for loc in root['qLoci']}
    assert {row['id'] for row in rows if row['n'] == 5} == expected


@pytest.mark.parametrize('damage', ['delete_locus', 'duplicate_locus', 'wrong_q',
                                   'wrong_lambda', 'wrong_dimension', 'odd_order',
                                   'alpha_zero', 'omega_zero', 'nonfinite',
                                   'identity_above_law', 'scalar_above_law',
                                   'missing_profile'])
def test_coefficient_validator_rejects_mutated_data(damage):
    validate = required_api('validate_coefficient_survey')
    rows = survey()
    if damage == 'delete_locus':
        rows.pop()
    elif damage == 'duplicate_locus':
        rows[-1] = copy.deepcopy(rows[-2])
    elif damage == 'wrong_q':
        rows[-1]['q']['real'] += 1
    elif damage == 'wrong_lambda':
        rows[-1]['lambda0']['real'] += 1
    elif damage == 'wrong_dimension':
        rows[-1]['dimension'] = 49
    elif damage == 'odd_order':
        rows[0]['profiles']['odd']['order'] = 1
    elif damage in ('alpha_zero', 'omega_zero'):
        rows[0]['profiles']['odd'][damage[:-5]] = {'real': 0., 'imag': 0.}
    elif damage == 'nonfinite':
        rows[0]['profiles']['one']['slopes']['real'][0] = float('nan')
    elif damage == 'identity_above_law':
        rows[0]['identity_relative_residual'] = 2*rows[0]['roundoff']['identity_relative_bound']
        scale = max(1., max(abs(np.asarray(rows[0]['profiles']['odd']['slopes']['real'])
                            + 1j*np.asarray(rows[0]['profiles']['odd']['slopes']['imag']))))
        rows[0]['identity_absolute_residual'] = rows[0]['identity_relative_residual']*scale
    elif damage == 'scalar_above_law':
        rows[0]['scalar_plane_residual'] = 1.
    elif damage == 'missing_profile':
        del rows[0]['profiles']['even']
    with pytest.raises(ValueError):
        validate(rows, inventory())


def test_corrupted_inventory_is_not_a_smaller_successful_survey():
    build = required_api('build_coefficient_survey')
    damaged = inventory()
    damaged['sectors'][0]['a2Roots'][0]['qLoci'].pop()
    with pytest.raises(ValueError, match='58|coverage'):
        build(damaged)


def test_continuation_validator_accepts_fresh_selected_pairs():
    validate = required_api('validate_continuations')
    validate(copy.deepcopy(fresh_continuations()))


@pytest.mark.parametrize('replacement_kind', ['opposite_lambda', 'opposite_q'])
def test_other_valid_spectral_cluster_cannot_replace_selected_seed(replacement_kind):
    data = copy.deepcopy(fresh_continuations())
    if replacement_kind == 'opposite_lambda':
        minus = next(row for row in fresh_survey() if row['id'] == 'N4-real-lambda-minus')
    else:
        plus = next(row for row in fresh_survey() if row['id'] == 'N4-real-lambda-plus')
        minus = unfolding.coefficients(4, -unfolding.z(plus['q']), unfolding.z(plus['lambda0']))
    replacement = unfolding.follow(4, unfolding.z(minus['q']),
                                   unfolding.z(minus['lambda0']), minus)
    data['N4-real-lambda-plus'] = copy.deepcopy(replacement)
    with pytest.raises(ValueError, match='selected.*seed|selected.*cluster'):
        unfolding.validate_continuations(data)


@pytest.mark.parametrize('field', ['q', 'lambda0', 'reflection_on_plane',
                                  'alpha', 'omega', 'slopes', 'next_coefficients'])
def test_every_encoded_field_rejects_extra_metadata(field):
    rows = survey()
    target = rows[0] if field in ('q', 'lambda0', 'reflection_on_plane') else rows[0]['profiles']['one']
    target[field]['metadata'] = 'not part of a complex encoding'
    with pytest.raises(ValueError, match='complex.*schema|complex.*keys'):
        unfolding.validate_coefficient_survey(rows, inventory())


def test_reflection_extra_key_cannot_hide_nested_nan():
    rows = survey()
    rows[0]['reflection_on_plane']['metadata'] = 'extra'
    rows[0]['reflection_on_plane']['real'][0][0] = 'nan'
    with pytest.raises(ValueError):
        unfolding.validate_coefficient_survey(rows, inventory())


@pytest.mark.parametrize('shape', ['scalar', 'vector', 'wrong_matrix'])
def test_reflection_requires_a_two_by_two_complex_matrix(shape):
    rows = survey()
    component = {'scalar': 0., 'vector': [0., 0.], 'wrong_matrix': [[0.], [0.]]}[shape]
    rows[0]['reflection_on_plane'] = {'real': component, 'imag': copy.deepcopy(component)}
    with pytest.raises(ValueError, match='complex.*shape'):
        unfolding.validate_coefficient_survey(rows, inventory())


@pytest.mark.parametrize('field', ['eigenvalue', 'reflection_character', 'odd_slopes'])
def test_optional_residue_complex_fields_reject_extra_metadata(field):
    rows = survey()
    control = unfolding.residue_control(rows[0])
    target = (control['responses']['full'] if field == 'odd_slopes'
              else control['nearest_modes'][0])
    target[field]['metadata'] = 'extra'
    rows[0]['residue_control'] = control
    with pytest.raises(ValueError, match='complex.*schema|complex.*keys'):
        unfolding.validate_coefficient_survey(rows, inventory())


@pytest.mark.parametrize('damage', ['shift_lambda', 'merged_singular_cluster',
                                   'rescaled_singular_values', 'rescaled_traceless_norm'])
def test_independent_continuation_diagnostics_reject_corruption(damage):
    data = copy.deepcopy(fresh_continuations())
    row = data['N4-real-lambda-plus']['odd'][0]
    if damage == 'shift_lambda':
        row['lambda0']['real'] += 100
    elif damage == 'merged_singular_cluster':
        row['smallest_singular_values'][-2] = row['smallest_singular_values'][-3]
    elif damage == 'rescaled_singular_values':
        row['smallest_singular_values'] = [2*s for s in row['smallest_singular_values']]
    elif damage == 'rescaled_traceless_norm':
        row['traceless_norm'] *= 2
    with pytest.raises(ValueError):
        unfolding.validate_continuations(data)


@pytest.mark.parametrize('field', ['alpha', 'omega'])
@pytest.mark.parametrize('encoded', ['nan', 'inf', '-inf', '1e9999'])
def test_encoded_nonfinite_coefficients_are_rejected(field, encoded):
    rows = survey()
    rows[0]['profiles']['odd'][field]['real'] = encoded
    with pytest.raises(ValueError, match='finite|numeric'):
        unfolding.validate_coefficient_survey(rows, inventory())


@pytest.mark.parametrize('encoded', [' 1 ', '+1', True])
def test_complex_decoder_rejects_noncanonical_numeric_components(encoded):
    with pytest.raises(ValueError, match='finite|numeric'):
        unfolding.z({'real': encoded, 'imag': 0})


@pytest.mark.parametrize('damage', ['missing_key', 'extra_key', 'missing_branch',
                                   'duplicate_branch', 'coalesced_q',
                                   'zero_second_singular', 'large_smallest_singular',
                                   'large_discriminant', 'nonfinite', 'wrong_separation'])
def test_continuation_validator_rejects_mutated_data(damage):
    validate = required_api('validate_continuations')
    data = copy.deepcopy(fresh_continuations())
    rows = data['N4-real-lambda-plus']['odd']
    if damage == 'missing_key':
        del data['N5-O-A2-W-006-Q-plus']
    elif damage == 'extra_key':
        data['unselected'] = copy.deepcopy(data['N4-real-lambda-plus'])
    elif damage == 'missing_branch':
        rows.pop()
    elif damage == 'duplicate_branch':
        rows[3]['branch'] = rows[0]['branch']
    elif damage == 'coalesced_q':
        rows[3]['q'] = copy.deepcopy(rows[0]['q'])
        rows[3]['generator_norm'] = rows[0]['generator_norm']
    elif damage == 'zero_second_singular':
        rows[0]['smallest_singular_values'][-2:] = [0., 0.]
    elif damage == 'large_smallest_singular':
        rows[0]['smallest_singular_values'][-1] = 1.
    elif damage == 'large_discriminant':
        rows[0]['scaled_discriminant_residual'] = 1.
    elif damage == 'nonfinite':
        rows[0]['lambda0']['imag'] = float('inf')
    elif damage == 'wrong_separation':
        rows[0]['branch_separation'] *= 2
    reason = {'coalesced_q': 'branches coalesced',
              'zero_second_singular': 'singular-value clearance',
              'large_discriminant': 'discriminant exceeds roundoff'}.get(damage)
    with pytest.raises(ValueError, match=reason):
        validate(data)


def test_provenance_primitives_have_a_neutral_module():
    assert importlib.util.find_spec('route_b_artifact_provenance') is not None
    import route_b_artifact_provenance as provenance
    import route_b_n4_virtual_readout as virtual
    assert virtual._canonical_text_sha256 is provenance.canonical_text_sha256
    assert provenance.canonical_text_sha256.__module__ == 'route_b_artifact_provenance'


@pytest.mark.parametrize('field', ['source_sha256', 'script_sha256', 'dependencies'])
def test_provenance_validator_rejects_wrong_hash(field):
    validate = required_api('validate_provenance')
    payload = unfolding.artifact_provenance(Path(unfolding.__file__),
                ROOT/'simulations/results/route_b_a2_n5.json')
    validate(copy.deepcopy(payload))
    damaged = copy.deepcopy(payload)
    if field == 'dependencies':
        damaged[field][next(iter(damaged[field]))] = '0'*64
    else:
        damaged[field] = '0'*64
    with pytest.raises(ValueError, match='provenance'):
        validate(damaged)


def test_scale_law_is_measured_across_more_than_seven_decades():
    required_api('build_coefficient_survey')
    q = np.sqrt((-1+np.sqrt(13))/6)
    ratios = []
    for scale in (2.**-12, 1., 2.**12):
        row = unfolding.coefficients(4, q, -4+2j*q, generator_scale=scale)
        error = row['roundoff']
        assert error['generator_norm'] > 0 and error['reduced_resolvent_norm'] > 0
        unit = (np.finfo(float).eps*row['dimension']*error['generator_norm']
                *error['reduced_resolvent_norm']*error['dual_norm']**2)
        assert error['identity_relative_bound'] == 16*unit
        ratios.append(row['identity_relative_residual']/unit)
    assert max(ratios) < 16


def test_every_locus_records_independent_scale_calibration():
    rows = survey()
    for row in rows:
        readings = row['roundoff'].get('scale_checks')
        assert readings is not None, 'missing measured calibration across scales'
        assert [r['generator_scale'] for r in readings] == [2.**-12, 1., 2.**12]
        for reading in readings:
            assert 0 <= reading['identity_ratio'] < 16
            assert reading['alpha_min_clearance'] > 1
            assert reading['omega_min_clearance'] > 1


def test_scale_law_mutation_reaches_production_validator():
    rows = survey()
    readings = rows[0]['roundoff'].get('scale_checks')
    assert readings is not None, 'missing measured calibration across scales'
    readings[-1]['identity_relative_residual'] = 2*readings[-1]['identity_relative_bound']
    with pytest.raises(ValueError, match='calibration|roundoff'):
        unfolding.validate_coefficient_survey(rows, inventory())


@pytest.mark.parametrize('stage', ['coefficients', 'continuations', 'provenance'])
def test_main_rejects_bad_data_before_any_write(stage, monkeypatch):
    required_api('validate_coefficient_survey')
    required_api('validate_continuations')
    required_api('validate_provenance')
    rows = survey()
    data = copy.deepcopy(fresh_continuations())
    if stage == 'coefficients':
        rows.pop()
    if stage == 'continuations':
        pair = data['N4-real-lambda-plus']['odd']
        pair[3]['q'] = copy.deepcopy(pair[0]['q'])
        pair[3]['generator_norm'] = pair[0]['generator_norm']
    monkeypatch.setattr(unfolding, 'build_coefficient_survey', lambda _: rows)
    monkeypatch.setattr(unfolding, 'follow', lambda n,q,lam,r: data[r['id']])
    if stage == 'provenance':
        original = unfolding.artifact_provenance
        def damaged_provenance(*args):
            result = original(*args)
            result['dependencies'][next(iter(result['dependencies']))] = '0'*64
            return result
        monkeypatch.setattr(unfolding, 'artifact_provenance', damaged_provenance)
    def no_write(*args, **kwargs):
        pytest.fail('invalid output reached a file write')
    monkeypatch.setattr(Path, 'write_text', no_write)
    reason = {'coefficients': 'coverage', 'continuations': 'branches coalesced',
              'provenance': 'provenance'}[stage]
    with pytest.raises(ValueError, match=reason):
        unfolding.main()


def test_saved_artifact_passes_the_same_production_validators():
    payload = json.loads((ROOT/'simulations/results/route_b_other_n_unfolding.json').read_text())
    unfolding.validate_coefficient_survey(payload['loci'], inventory())
    unfolding.validate_continuations(payload['continuations'])
    unfolding.validate_provenance(payload)
