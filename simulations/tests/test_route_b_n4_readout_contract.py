"""The EP anchor alone cannot identify the end-bias convention."""
import copy
import hashlib
import importlib
import json
import os
from pathlib import Path
import subprocess
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
SOURCE = ROOT/'simulations/results/route_b_n4_self_fold.json'


def imported_local_paths(producer):
    """Observe actual imports in a fresh process, independently of the manifest."""
    result = subprocess.run([sys.executable, '-c', '''
import importlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root/'simulations'))
producer = importlib.import_module(sys.argv[2])
paths = {Path(module.__file__).resolve() for module in sys.modules.copy().values()
         if getattr(module, '__file__', None)}
print(json.dumps(sorted(path.relative_to(root).as_posix() for path in paths
                        if path.is_relative_to(root) and path.suffix == '.py'
                        and path != Path(producer.__file__).resolve())))
''', str(ROOT), producer], check=True, capture_output=True, text=True, cwd=ROOT)
    return set(json.loads(result.stdout))


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


@pytest.mark.parametrize('producer', PRODUCERS)
def test_provenance_recomputes_source_script_and_all_producer_dependencies(producer):
    # Text-mode universal newlines are independent of the production hash route.
    def digest(path):
        return hashlib.sha256(path.read_text(encoding='utf-8').encode('utf-8')).hexdigest()

    provenance = required_api('artifact_provenance')(ROOT/f'simulations/{producer}.py')
    assert provenance['source_sha256'] == digest(SOURCE)
    assert provenance['script_sha256'] == digest(ROOT/f'simulations/{producer}.py')
    observed = imported_local_paths(producer)
    assert 'simulations/framework/__init__.py' in observed
    assert 'simulations/framework/workflows/_propagation.py' in observed
    assert observed <= set(provenance['dependencies'])
    for path, recorded_digest in provenance['dependencies'].items():
        assert recorded_digest == digest(ROOT/path)
    required_api('verify_artifact_provenance')(provenance, ROOT/f'simulations/{producer}.py')


@pytest.mark.parametrize('producer', PRODUCERS)
@pytest.mark.parametrize('field', ['script_sha256', 'source_sha256', 'dependencies'])
def test_provenance_survives_checkout_line_endings(producer, field, monkeypatch):
    script = ROOT/f'simulations/{producer}.py'
    provenance = virtual.artifact_provenance(script)
    target = {'script_sha256': script, 'source_sha256': SOURCE,
              'dependencies': ROOT/'simulations/framework/__init__.py'}[field].resolve()
    read_bytes = Path.read_bytes
    original = read_bytes(target)
    lf = original.replace(b'\r\n', b'\n').replace(b'\r', b'\n')
    assert b'\n' in lf
    for ending in (b'\n', b'\r\n', b'\r'):
        checkout = lf.replace(b'\n', ending)

        def checkout_bytes(path):
            return checkout if path.resolve() == target else read_bytes(path)

        with monkeypatch.context() as context:
            context.setattr(Path, 'read_bytes', checkout_bytes)
            virtual.verify_artifact_provenance(provenance, script)
            assert virtual.artifact_provenance(script) == provenance


def test_provenance_rejects_binary_source_extension(monkeypatch):
    monkeypatch.setattr(virtual, 'SOURCE', ROOT/'visualizations/route_b_n4_virtual_readout.png')
    with pytest.raises(ValueError, match='provenance.*text'):
        virtual.artifact_provenance(virtual.__file__)


@pytest.mark.parametrize('binary_content', [b'\xff\r\nbinary', b'\x00\r\nbinary'])
def test_provenance_rejects_binary_bytes_in_json_source(binary_content, monkeypatch):
    read_bytes = Path.read_bytes

    def binary_source(path):
        return binary_content if path.resolve() == SOURCE.resolve() else read_bytes(path)

    monkeypatch.setattr(Path, 'read_bytes', binary_source)
    with pytest.raises(ValueError, match='provenance.*text'):
        virtual.artifact_provenance(virtual.__file__)


@pytest.mark.parametrize('producer', PRODUCERS)
def test_framework_initializer_byte_mutation_invalidates_provenance(producer, monkeypatch):
    script = ROOT/f'simulations/{producer}.py'
    provenance = virtual.artifact_provenance(script)
    initializer = (ROOT/'simulations/framework/__init__.py').resolve()
    read_bytes = Path.read_bytes

    def mutated_bytes(path):
        content = read_bytes(path)
        if path.resolve() == initializer:
            return content+b'\n# framework initializer provenance mutation\n'
        return content

    monkeypatch.setattr(Path, 'read_bytes', mutated_bytes)
    with pytest.raises(ValueError, match='provenance.*dependencies'):
        virtual.verify_artifact_provenance(provenance, script)


@pytest.mark.parametrize('producer', PRODUCERS)
@pytest.mark.parametrize('field', ['script_sha256', 'source_sha256', 'dependencies'])
def test_wrong_hash_is_rejected_through_production_validator(producer, field):
    script = ROOT/f'simulations/{producer}.py'
    provenance = required_api('artifact_provenance')(script)
    wrong = copy.deepcopy(provenance)
    if field == 'dependencies':
        for dependency in provenance['dependencies']:
            wrong = copy.deepcopy(provenance)
            wrong[field][dependency] = '0'*64
            with pytest.raises(ValueError, match='provenance.*dependencies'):
                required_api('verify_artifact_provenance')(wrong, script)
    else:
        wrong[field] = '0'*64
        with pytest.raises(ValueError, match=f'provenance.*{field}'):
            required_api('verify_artifact_provenance')(wrong, script)


@pytest.mark.parametrize('producer', PRODUCERS)
def test_committed_artifact_hashes_match_current_files(producer):
    output = json.loads((ROOT/f'simulations/results/{producer}.json').read_text(encoding='utf-8'))
    assert 'script_sha256' in output, 'artifact lacks script provenance'
    assert 'source_sha256' in output, 'artifact lacks source provenance'
    required_api('verify_artifact_provenance')(output, ROOT/f'simulations/{producer}.py')
