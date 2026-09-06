import numpy as np
import pytest

import qudit_ti_intermediate as qti
from qudit_ti_intermediate import validate_cited_finite_rows
from qutrit_partial_palindrome import palindrome_pairs
from su3_heisenberg_rep_theory import main as su3_main, right_hs_rayleigh


def test_qudit_ti_cited_finite_rows_are_literal_and_mutation_sensitive():
    rows = [(3, 2, 36, 54, 54), (3, 3, 216, 378, 378), (4, 2, 64, 128, 128)]
    validate_cited_finite_rows(rows)

    mutated = list(rows)
    mutated[0] = (3, 2, 36, 36, 54)
    with pytest.raises(AssertionError, match="cited finite TI construction"):
        validate_cited_finite_rows(mutated)


@pytest.mark.parametrize("d,N,expected", [(3, 2, 54), (3, 3, 378), (4, 2, 128)])
def test_public_ceiling_matches_the_cited_combinatorial_ceiling(d, N, expected):
    assert qti.ceiling(d, N) == expected


def test_global_palindrome_matching_is_permutation_invariant_and_uses_multiplicity():
    values = np.array([10, 10.8, -10.7, -9.1], dtype=complex)
    assert palindrome_pairs(values, 0.0, tol=1.0) == 4
    assert palindrome_pairs(values[[2, 0, 3, 1]], 0.0, tol=1.0) == 4


def test_global_palindrome_matching_has_a_negative_control():
    values = np.array([10, 10.8, -10.7, -7.0], dtype=complex)
    assert palindrome_pairs(values, 0.0, tol=1.0) == 2


def test_right_hs_rayleigh_uses_conjugate_right_vectors_and_stays_real():
    diagonal = np.array([0.0, 2.0])
    right_vectors = np.array([[1.0, 1.0 + 1.0j], [1.0j, 2.0]])

    measured = right_hs_rayleigh(diagonal, right_vectors)

    expected = np.array([1.0, 4.0 / 3.0])
    np.testing.assert_allclose(measured.real, expected, atol=1e-14)
    np.testing.assert_allclose(measured.imag, 0.0, atol=1e-14)


def test_qudit_ti_producer_routes_computed_rows_through_literal_gate(capsys, monkeypatch):
    routed_rows = []
    original = qti.validate_cited_finite_rows

    def record_and_validate(rows):
        original(rows)
        routed_rows.extend(rows)

    monkeypatch.setattr(qti, "validate_cited_finite_rows", record_and_validate)
    qti.main()
    assert "TI REACHES ceiling" in capsys.readouterr().out
    assert {(d, N): ti for d, N, _cap, ti, _ceil in routed_rows} == {
        (2, 2): 16,
        (3, 2): 54,
        (3, 3): 378,
        (4, 2): 128,
        (2, 3): 64,
    }


def test_su3_producer_uses_right_hs_and_full_complex_matching(capsys):
    su3_main()
    output = capsys.readouterr().out
    assert "max |biorthogonal - right HS| = 0.025008" in output
    assert "full-complex palindrome matching = 60/81" in output
    assert "six Q=1 modes near |Im|=3.99875" in output
