from pathlib import Path

import numpy as np

from simulations import rmt_analysis
import pytest
import shutil
import subprocess
import sys


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "rmt_analysis.py"
RESULT = ROOT / "results" / "rmt_analysis.txt"


def test_gap_ratios_are_affine_invariant_but_not_non_affine_invariant():
    levels = np.array([0.0, 1.0, 2.0])
    affine = 7.0 * levels + 3.0
    squared = levels**2

    base, _ = rmt_analysis.spacing_ratios(levels)
    affine_value, _ = rmt_analysis.spacing_ratios(affine)
    squared_value, _ = rmt_analysis.spacing_ratios(squared)

    assert np.isclose(base, affine_value)
    assert not np.isclose(base, squared_value)


@pytest.mark.parametrize("levels, expected", [
    ([0., 0., 1., 2.], [0., 1.]),  # Deleting the zero would invent a different adjacent pair.
    ([0., 0., 0., 1.], [np.nan, 0.]),  # 0/0 is undefined, not zero or silently reconnected.
    ([1., 1., 1., 1.], [np.nan, np.nan]),
    ([0., 1e-14, 3e-14], [.5]),  # No absolute cutoff may delete small positive gaps.
])
def test_direct_multiset_gap_ratios_preserve_each_adjacent_gap_pair(levels, expected):
    mean, ratios = rmt_analysis.spacing_ratios(np.array(levels))
    np.testing.assert_allclose(ratios, expected, rtol=1e-15, atol=0, equal_nan=True)
    defined = np.array(expected)[~np.isnan(expected)]
    if len(defined):
        assert mean == pytest.approx(np.mean(defined), rel=1e-15, abs=0)
    else:
        assert np.isnan(mean)


def test_population_counts_distinguish_zero_gaps_from_undefined_ratios():
    reading = rmt_analysis.spacing_population(np.array([0., 0., 0., 1., 2.]))
    assert reading["levels"] == 5
    assert reading["gaps"] == 4
    assert reading["zero_gaps"] == 2
    assert reading["adjacent_pairs"] == 3
    assert reading["defined_ratios"] == 2
    assert reading["undefined_zero_zero"] == 1
    assert reading["mean"] == .5
    control = rmt_analysis.spacing_population(np.array([0., 1., 2., 3., 4.]))
    assert control["zero_gaps"] == control["undefined_zero_zero"] == 0
    assert control["mean"] == 1.


def test_rank_equalizing_map_is_rejected_because_it_forces_unit_spacings():
    levels = np.array([0.0, 0.2, 1.7, 9.0, 20.0])
    rank_equalized = np.arange(1, len(levels) + 1) / len(levels)

    assert not np.allclose(np.diff(levels), np.diff(levels)[0])
    assert np.allclose(np.diff(rank_equalized), np.diff(rank_equalized)[0])
    assert not hasattr(rmt_analysis, "unfold")


def test_source_and_artifact_report_diagnostics_without_class_verdicts():
    for path in (SCRIPT, RESULT):
        source = path.read_text(encoding="utf-8")
        for stale in ("NNSD", "Nearest:", "Class (all)", "Class (half)", "classify_r", "def unfold", "Poisson=0.386"):
            assert stale not in source
        for required in ("affine-invariant", "no irreducible-class inference",
                         "no cutoff or deduplication", "0/0 = undefined",
                         "No standard Poisson/GOE/GUE calibration"):
            assert required in source


def test_independent_artifact_gate_preserves_zero_gap_population_and_detects_mutations(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT))
    from simulations.docs_verify import rmt_artifact_row_errors
    levels = np.array([0., 0., 0., 1., 2.])
    artifact = '''N=2: 5 eigenvalues
  (a) Raw full multiset:
      levels=5, gaps=4, zero_gaps=2, adjacent_pairs=3, defined_ratios=2, undefined_zero_zero=1
      <r> defined = 0.5
  (b) Raw rates below N*gamma (not an irreducible sector or a decorrelation certificate):
      levels=3, gaps=2, zero_gaps=2, adjacent_pairs=1, defined_ratios=0, undefined_zero_zero=1
      <r> defined = nan
'''
    inputs = {2: -levels.astype(complex)}
    assert not rmt_artifact_row_errors(artifact, inputs)
    for old, new in (("zero_gaps=2", "zero_gaps=0"), ("defined_ratios=2", "defined_ratios=3"),
                     ("undefined_zero_zero=1", "undefined_zero_zero=0"),
                     ("<r> defined = 0.5", "<r> defined = 1"), ("levels=3", "levels=0")):
        assert rmt_artifact_row_errors(artifact.replace(old, new, 1), inputs)


def test_band_population_is_a_multiset_not_a_rounded_unique_xy_sector(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT))
    from simulations import rmt_goe_hint_verdict as bands
    # gamma=.05, w=2: rate window (.185,.215); include zero, duplicates,
    # and tiny distinct frequencies, excluding only values outside the band.
    eigenvalues = np.array([-.2+0j, -.2+1j, -.2-1j, -.2+1.000000001j, -.3+2j])
    actual = bands.band_frequency_multiset(eigenvalues, .05, 2)
    np.testing.assert_array_equal(actual, [0., 1., 1., 1.000000001])
    np.testing.assert_array_equal(bands.band_frequency_multiset(eigenvalues, .05, 3), [2.])


def test_band_producer_has_no_ensemble_calibration_or_deleted_api():
    source = (ROOT / 'rmt_goe_hint_verdict.py').read_text(encoding='utf-8')
    for stale in ('classify_r', 'np.unique', 'np.round', 'poisson_band', 'Poisson <r>=', 'XY-weight-w sector'):
        assert stale not in source
    for boundary in ('average-light bands', 'zero frequencies and multiplicities retained',
                     'no standard-ensemble calibration'):
        assert boundary in source


def test_independent_band_artifact_gate_rejects_deduplication_and_wrong_window(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT))
    from simulations.docs_verify import rmt_band_artifact_errors
    inputs = {2: np.array([-.1+0j, -.1+1j, -.1-1j, -.2+2j])}
    artifact = '''N=2: 4 eigenvalues
  w=0:
      levels=0, gaps=0, zero_gaps=0, adjacent_pairs=0, defined_ratios=0, undefined_zero_zero=0
      <r> defined = nan
  w=1:
      levels=3, gaps=2, zero_gaps=1, adjacent_pairs=1, defined_ratios=1, undefined_zero_zero=0
      <r> defined = 0
  w=2:
      levels=1, gaps=0, zero_gaps=0, adjacent_pairs=0, defined_ratios=0, undefined_zero_zero=0
      <r> defined = nan
'''
    assert not rmt_band_artifact_errors(artifact, inputs)
    for old, new in (("levels=3", "levels=2"), ("zero_gaps=1", "zero_gaps=0"),
                     ("<r> defined = 0", "<r> defined = 1"), ("w=2:", "w=3:")):
        assert rmt_band_artifact_errors(artifact.replace(old, new, 1), inputs)


@pytest.mark.parametrize('producer, output', [('rmt_analysis.py', 'rmt_analysis.txt'),
                                           ('rmt_goe_hint_verdict.py', 'rmt_band_multiset.txt')])
def test_csv_producer_regenerates_artifact_in_isolated_directory(tmp_path, producer, output):
    # Running copies is a safe output seam: no tracked file can be mutated by this test.
    for script in ('rmt_analysis.py', 'rmt_goe_hint_verdict.py'):
        shutil.copy2(ROOT / script, tmp_path / script)
    result_dir = tmp_path / 'results'
    result_dir.mkdir()
    for n in range(2, 8):
        shutil.copy2(ROOT / 'results' / f'rmt_eigenvalues_N{n}.csv', result_dir)
    tracked = (ROOT / 'results' / output).read_bytes()
    subprocess.run([sys.executable, str(tmp_path / producer)], check=True, capture_output=True, timeout=30)
    assert (result_dir / output).read_bytes() == tracked
    assert (ROOT / 'results' / output).read_bytes() == tracked
    assert b'undefined_zero_zero=' in tracked  # Wrong-output control: an old summary lacks this contract.


def test_document_tables_are_tied_to_executed_artifacts(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT))
    from simulations.docs_verify import rmt_document_table_errors
    doc = (ROOT.parent / 'experiments/RANDOM_MATRIX_THEORY.md').read_text(encoding='utf-8')
    artifacts = [(ROOT / 'results' / name).read_text(encoding='utf-8')
                 for name in ('rmt_analysis.txt', 'rmt_band_multiset.txt')]
    assert not rmt_document_table_errors(doc, *artifacts)
    for old, new in (('| 16384 | 136 | 16366 | 16 |', '| 16384 | 136 | 16382 | 0 |'),
                     ('0.2021456120489688 |', '0.383 |'),
                     ('0.2021456120489688 |', 'nan |'),
                     ('| 7 | 3 | 1710 |', '| 7 | 3 | 414 |')):
        assert old in doc
        assert rmt_document_table_errors(doc.replace(old, new, 1), *artifacts)


def test_document_does_not_infer_cptp_from_eigenvalue_real_parts():
    doc = ' '.join((ROOT.parent / 'experiments/RANDOM_MATRIX_THEORY.md')
                   .read_text(encoding='utf-8').split())
    assert 'does not establish complete positivity or trace preservation' in doc
    assert 'The Lindblad construction supplies those properties independently' in doc
    assert 'confirming that the Liouvillian is a proper generator' not in doc


def test_document_names_f1_palindrome_without_global_chiral_shortcut():
    doc = ' '.join((ROOT.parent / 'experiments/RANDOM_MATRIX_THEORY.md')
                   .read_text(encoding='utf-8').split())
    assert 'The F1 palindrome (the centered spectrum has exact ± pairing)' in doc
    assert 'The chiral symmetry (centered spectrum has exact ± pairing)' not in doc
    assert 'chiral symmetry Liouvillian eigenvalues' not in doc


def test_result2_labels_the_executed_real_part_projection_not_full_complex_matching():
    source = SCRIPT.read_text(encoding='utf-8')
    result = RESULT.read_text(encoding='utf-8')
    doc = (ROOT.parent / 'experiments/RANDOM_MATRIX_THEORY.md').read_text(encoding='utf-8')
    combined = source + result + doc
    assert 'Centered Decay-Rate Reflection Check' in doc
    assert 'multiplicity-preserving perfect matching of the full real-part projection' in result
    assert 'does not match the full complex multiset' in combined
    assert 'Result 2: Exact Centered Spectral Pairing' not in doc
    assert '| ± pairs | Mean pairing error |' not in doc


def test_centered_rate_projection_rejects_an_unmatched_reflected_multiplicity():
    # Centered rates {-3,-1,+1} have one tempting +/- pair but an unmatched -3.
    # Truncating the two sign lists would report zero; a bijective matcher must not.
    eigenvalues = -np.array([-3.0, -1.0, 1.0], dtype=complex)
    reading = rmt_analysis.centered_rate_projection_reading(eigenvalues, sigma_gamma=0.0)
    assert reading['positive'] == 1
    assert reading['negative'] == 2
    assert reading['central'] == 0
    assert reading['bottleneck_error'] == pytest.approx(2.0)


def test_framework_connection_does_not_claim_unmeasured_cause_of_pooled_ratios():
    doc = (ROOT.parent / 'experiments/RANDOM_MATRIX_THEORY.md').read_text(encoding='utf-8')
    normalized = ' '.join(doc.split())
    assert 'does not isolate block mixing as the cause' in normalized
    assert 'does not measure whether that pairing changes the pooled statistic' in normalized
    assert 'What does the work is the conserved structure' not in doc
    assert 'it is not the one that suppresses repulsion here' not in doc
