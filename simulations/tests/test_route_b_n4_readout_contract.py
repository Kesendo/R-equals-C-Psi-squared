"""The EP anchor alone cannot identify the end-bias convention."""
import importlib
import os
from pathlib import Path
import sys

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'simulations'))

import numpy as np
import pytest

import route_b_n4_virtual_readout as virtual
from route_b_other_n_unfolding import parts

PRODUCERS = (
    'route_b_n4_virtual_readout',
    'route_b_n4_readout_search',
    'route_b_n4_histogram_filter',
)


def required_api(name):
    function = getattr(virtual, name, None)
    assert callable(function), f'missing production validator: {name}'
    return function


def doubled_end_bias(epsilon, q=2):
    return virtual.generator(virtual.EP+2*(epsilon-virtual.EP), q=q)


def test_generator_contract_accepts_independent_block_at_all_three_biases():
    checked = []

    def recording_builder(epsilon, q=2):
        checked.append((epsilon, q))
        return virtual.generator(epsilon, q=q)

    verify = required_api('verify_generator_contract')
    verify()  # The public default is also the production generator.
    verify(builder=recording_builder)
    assert checked == [(virtual.EP+shift, 2) for shift in (0., -.05, .05)]


def test_factor_two_mutation_passes_old_ep_anchor_but_fails_new_gate():
    d, c, left, right, _ = parts(4)
    indices = [16*a+b for a in (1, 2, 4, 8) for b in (3, 5, 6, 9, 10, 12)]
    expected = d+2*(c+virtual.EP*(left+right)/2)
    # This is the old production anchor, reached with the wrong builder.
    residual = np.linalg.norm(doubled_end_bias(virtual.EP)[np.ix_(indices, indices)]-expected)
    assert residual < 1e-13
    assert residual == 0.0
    verify = required_api('verify_generator_contract')
    with pytest.raises(ValueError, match='generator contract.*epsilon'):
        verify(builder=doubled_end_bias)


@pytest.mark.parametrize('producer', PRODUCERS)
def test_every_entrypoint_rejects_mutated_generator_before_writing(producer, monkeypatch):
    required_api('verify_generator_contract')
    module = importlib.import_module(producer)
    generator = module.generator

    def wrong(epsilon, q=2):
        return generator(virtual.EP+2*(epsilon-virtual.EP), q=q)

    monkeypatch.setattr(module, 'generator', wrong)
    writes = []

    def forbidden_write(*args, **kwargs):
        writes.append(args)
        pytest.fail('output write reached before rejecting the wrong generator')

    monkeypatch.setattr(Path, 'write_text', forbidden_write)
    with pytest.raises(ValueError, match='generator contract.*epsilon'):
        module.main()
    assert writes == []
