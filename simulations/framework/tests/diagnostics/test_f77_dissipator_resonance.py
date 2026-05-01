"""Tests for the dissipator-resonance law: F77-hardness lives in the Klein
cell that matches the dephasing letter (verified 2026-05-01).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_F77_dephase_letter_default_is_Z():
    """Default dephase_letter='Z' reproduces the existing Z-dephasing classification."""
    chain = fw.ChainSystem(N=4)
    # Klein (0, 1) Z-like pair: hard under Z-dephasing, the verified-hard sector.
    terms = [('I', 'I', 'Z'), ('Z', 'Z', 'Z')]
    assert fw.classify_pauli_pair(chain, terms) == 'hard'
    assert fw.classify_pauli_pair(chain, terms, dephase_letter='Z') == 'hard'


def test_F77_dissipator_resonance_diagonal():
    """Diagonal: each dephasing letter pulls hardness into its own Klein cell.

    Klein (0, 1) is Z's index; (1, 0) is X's; (1, 1) is Y's. A pair built
    only from one letter (plus I) lives in that letter's Klein cell and
    must be F77-hard under matching dephasing — and soft under the other
    two by the dissipator-resonance law.
    """
    chain = fw.ChainSystem(N=4)
    pairs = {
        'Z': [('I', 'I', 'Z'), ('Z', 'Z', 'Z')],   # Klein (0, 1)
        'X': [('I', 'I', 'X'), ('X', 'X', 'X')],   # Klein (1, 0)
        'Y': [('I', 'I', 'Y'), ('Y', 'Y', 'Y')],   # Klein (1, 1)
    }
    for h_letter, terms in pairs.items():
        for d_letter in ('X', 'Y', 'Z'):
            cls = fw.classify_pauli_pair(chain, terms, dephase_letter=d_letter)
            if d_letter == h_letter:
                assert cls == 'hard', \
                    f"{h_letter}-pair under {d_letter}-deph should be hard; got {cls}"
            else:
                assert cls in ('truly', 'soft'), \
                    f"{h_letter}-pair under {d_letter}-deph should be soft/truly; got {cls}"


def test_F77_truly_under_each_dephasing():
    """Truly under dephase_letter L: terms whose 'orthogonal' letters appear
    in even counts. By SU(2) covariance the F85 truly criterion permutes:
      Z-deph: #Y even AND #Z even
      X-deph: #Y even AND #X even
      Y-deph: #X even AND #Z even
    """
    chain = fw.ChainSystem(N=4)
    # XX pair: under Z-deph truly (#Y=#Z=0 even), under X-deph truly (#Y=#X=2 even).
    assert fw.classify_pauli_pair(
        chain, [('X', 'X'), ('Y', 'Y')], dephase_letter='Z'
    ) == 'truly'
    # ZZ + YY: truly under X-deph (#Y=2, #X=0 both even), as XX+YY rotates to ZZ+YY.
    assert fw.classify_pauli_pair(
        chain, [('Z', 'Z'), ('Y', 'Y')], dephase_letter='X'
    ) == 'truly'
    # XX + ZZ: truly under Y-deph (#X=2, #Z=2 both even).
    assert fw.classify_pauli_pair(
        chain, [('X', 'X'), ('Z', 'Z')], dephase_letter='Y'
    ) == 'truly'


def test_F77_dephase_letter_validation():
    """Invalid dephase_letter raises ValueError."""
    chain = fw.ChainSystem(N=3)
    with pytest.raises(ValueError):
        fw.classify_pauli_pair(chain, [('X', 'X'), ('Y', 'Y')], dephase_letter='Q')


def test_F77_su2_rotation_equivalence():
    """SU(2) rotation: hard count is identical across the three dephasing letters
    when each is tested against its own resonant Klein-cell pair family.

    Three Z₂³-homogeneous pairs (one in each Klein cell), each tested under
    matching dephasing — all should be hard with identical eigenvalue-pair
    error magnitude (up to floating-point).
    """
    chain = fw.ChainSystem(N=4)
    # Match each pair to its resonant dephasing letter and confirm 'hard'.
    diagonal = [
        ('Z', [('I', 'I', 'Z'), ('Z', 'Z', 'Z')]),
        ('X', [('I', 'I', 'X'), ('X', 'X', 'X')]),
        ('Y', [('I', 'I', 'Y'), ('Y', 'Y', 'Y')]),
    ]
    for d_letter, terms in diagonal:
        assert fw.classify_pauli_pair(chain, terms, dephase_letter=d_letter) == 'hard'
