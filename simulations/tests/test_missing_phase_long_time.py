"""Exact output-map controls for the N=7 doubled-copy experiment."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from missing_phase_long_time import Run, decoded, pair_distances


def test_two_site_map_separates_real_from_imaginary_coherence():
    a = np.zeros((1, 7, 7), complex)
    a[0, 3, 0] = .5j
    a[0, 0, 3] = -.5j
    assert np.max(pair_distances(a)) == 0
    a[0, 3, 0] = a[0, 0, 3] = .5
    assert np.max(pair_distances(a)) == .5


def test_decoder_reads_B_not_A_on_the_centre_row():
    a = np.zeros((1, 7, 7), complex)
    b = a.copy()
    a[0, 3, 0] = .5
    a[0, 0, 3] = .5
    assert np.count_nonzero(decoded(a, b)) == 0
    b[0, 3, 0] = .5j
    b[0, 0, 3] = -.5j
    out = decoded(a, b)[0]
    assert out[6, 0] == .5j
    assert out[0, 6] == -.5j
    assert np.trace(out @ out) == .5


def test_decoder_retains_population_and_noncentre_coherence():
    a = np.zeros((1, 7, 7), complex)
    a[0, 3, 3] = .5
    a[0, 0, 0] = .5
    a[0, 0, 6] = a[0, 6, 0] = .25
    out = decoded(a, np.zeros_like(a))[0]
    assert out[6, 6] == .5
    assert out[0, 0] == .5
    assert out[0, 5] == out[5, 0] == .25
    assert np.trace(out) == 1


def test_blocks_against_full_doubled_lindblad_generator():
    from framework.lindblad import lindbladian_general
    run = Run(.1, .25)
    z = np.eye(7); z[3, 3] = -1
    h_full = np.kron(np.eye(2), run.h)
    z_full = np.kron(np.diag([1, -1]), z)
    full = lindbladian_general(h_full, [.5*z_full])
    a_indices = [i*14+j for i in range(7) for j in range(7)]
    b_indices = [i*14+7+j for i in range(7) for j in range(7)]
    assert np.array_equal(full[np.ix_(a_indices, a_indices)], run.la)
    assert np.array_equal(full[np.ix_(b_indices, b_indices)], run.lb)
    assert not np.array_equal(full[np.ix_(a_indices, a_indices)], run.lb)


def test_partner_rate_uses_gamma_without_a_square_root_round_trip():
    run = Run(.1, .3)
    assert np.array_equal(run.lb + run.la.conj().T, -.6*np.eye(49))
