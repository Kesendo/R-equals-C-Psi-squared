"""Tests for PauliTerm and PauliHamiltonian raw structural classes."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_PauliTerm_klein_index():
    """Klein index = (bit_a, bit_b) parity tuple of letter sequence."""
    # Single letters
    assert fw.PauliTerm(('I',)).klein_index == (0, 0)
    assert fw.PauliTerm(('X',)).klein_index == (1, 0)
    assert fw.PauliTerm(('Y',)).klein_index == (1, 1)
    assert fw.PauliTerm(('Z',)).klein_index == (0, 1)

    # Bilinears (k=2)
    assert fw.PauliTerm(('X', 'X')).klein_index == (0, 0)  # M
    assert fw.PauliTerm(('Y', 'Y')).klein_index == (0, 0)  # M
    assert fw.PauliTerm(('Z', 'Z')).klein_index == (0, 0)  # M
    assert fw.PauliTerm(('X', 'Y')).klein_index == (0, 1)  # F_a
    assert fw.PauliTerm(('Y', 'X')).klein_index == (0, 1)  # F_a
    assert fw.PauliTerm(('X', 'Z')).klein_index == (1, 1)  # F_b
    assert fw.PauliTerm(('Z', 'X')).klein_index == (1, 1)  # F_b
    assert fw.PauliTerm(('Y', 'Z')).klein_index == (1, 0)  # C
    assert fw.PauliTerm(('Z', 'Y')).klein_index == (1, 0)  # C


def test_PauliTerm_y_parity_independent_at_k3():
    """Y-parity is determined by Klein index at k=2 (Y-par = bit_a XOR bit_b)
    but becomes an independent third Z₂ axis at k≥3."""
    # k=2: Y-par = bit_a XOR bit_b
    for letters in [('X', 'X'), ('X', 'Y'), ('X', 'Z'), ('Y', 'Z'), ('Y', 'Y')]:
        t = fw.PauliTerm(letters)
        a, b = t.klein_index
        assert t.y_parity == (a ^ b), f"k=2 Y-par should equal bit_a XOR bit_b for {letters}"

    # k=3: Y-par independent of Klein index
    # XXY: Klein = (1, 1), Y-par = 1; XOR(1,1)=0 ≠ 1
    t_xxy = fw.PauliTerm(('X', 'X', 'Y'))
    assert t_xxy.klein_index == (1, 1)
    assert t_xxy.y_parity == 1
    assert t_xxy.y_parity != (t_xxy.klein_index[0] ^ t_xxy.klein_index[1])

    # XZZ: Klein = (1, 0), Y-par = 0; XOR(1,0)=1 ≠ 0
    t_xzz = fw.PauliTerm(('X', 'Z', 'Z'))
    assert t_xzz.klein_index == (1, 0)
    assert t_xzz.y_parity == 0
    assert t_xzz.y_parity != (t_xzz.klein_index[0] ^ t_xzz.klein_index[1])


def test_PauliTerm_pi2_class():
    """Π²-class = truly / pi2_odd / pi2_even_nontruly."""
    assert fw.PauliTerm(('X', 'X')).pi2_class == 'truly'
    assert fw.PauliTerm(('Y', 'Y')).pi2_class == 'truly'
    assert fw.PauliTerm(('Z', 'Z')).pi2_class == 'truly'
    assert fw.PauliTerm(('X', 'Y')).pi2_class == 'pi2_odd'
    assert fw.PauliTerm(('Y', 'Z')).pi2_class == 'pi2_even_nontruly'


def test_PauliTerm_is_truly():
    assert fw.PauliTerm(('X', 'X')).is_truly is True
    assert fw.PauliTerm(('I', 'X')).is_truly is True   # #Y=#Z=0
    assert fw.PauliTerm(('X', 'Y')).is_truly is False
    assert fw.PauliTerm(('Y', 'Z')).is_truly is False


def test_PauliTerm_invalid_letter_raises():
    with pytest.raises(ValueError):
        fw.PauliTerm(('X', 'Q'))


def test_PauliTerm_full_z2_signature():
    """Full Z₂³ signature (bit_a, bit_b, y_parity)."""
    assert fw.PauliTerm(('X', 'X')).full_z2_signature == (0, 0, 0)
    assert fw.PauliTerm(('X', 'Y')).full_z2_signature == (0, 1, 1)
    assert fw.PauliTerm(('X', 'X', 'Y')).full_z2_signature == (1, 1, 1)
    assert fw.PauliTerm(('X', 'Z', 'Z')).full_z2_signature == (1, 0, 0)


def test_PauliHamiltonian_klein_homogeneous():
    """Klein-homogeneous Hamiltonians have all terms in the same Klein slot."""
    # Pure Π²-odd (XY+YX): both F_a → homogeneous
    H_homo = fw.PauliHamiltonian.from_letter_tuples([('X', 'Y'), ('Y', 'X')], chain_length=3)
    assert H_homo.is_klein_homogeneous is True
    assert H_homo.klein_set == {(0, 1)}

    # Truly (XX+YY): both M → homogeneous
    H_truly = fw.PauliHamiltonian.from_letter_tuples([('X', 'X'), ('Y', 'Y')], chain_length=3)
    assert H_truly.is_klein_homogeneous is True
    assert H_truly.klein_set == {(0, 0)}

    # Mixed (XY+YZ): F_a + C → inhomogeneous
    H_mixed = fw.PauliHamiltonian.from_letter_tuples([('X', 'Y'), ('Y', 'Z')], chain_length=3)
    assert H_mixed.is_klein_homogeneous is False
    assert H_mixed.klein_set == {(0, 1), (1, 0)}


def test_PauliHamiltonian_y_parity_homogeneous_k3():
    """At k=3, Y-parity homogeneity is independent from Klein homogeneity."""
    # Two k=3 terms with same Klein but different Y-parity
    # XYZ: Klein = (0, 0), Y-par = 1
    # YZX: Klein = (0, 0), Y-par = 1
    H_klein_homo_y_homo = fw.PauliHamiltonian.from_letter_tuples(
        [('X', 'Y', 'Z'), ('Y', 'Z', 'X')], chain_length=4
    )
    assert H_klein_homo_y_homo.is_klein_homogeneous is True
    assert H_klein_homo_y_homo.is_y_parity_homogeneous is True

    # XYZ (Klein 0,0, Y-par 1) + IXZ (Klein 1,1 since 1+0+0=1, 0+0+1=1; Y-par 0)
    # Klein-inhomogeneous AND Y-parity-inhomogeneous
    H_inhomo = fw.PauliHamiltonian.from_letter_tuples(
        [('X', 'Y', 'Z'), ('I', 'X', 'Z')], chain_length=4
    )
    assert H_inhomo.is_klein_homogeneous is False
    assert H_inhomo.is_y_parity_homogeneous is False

    # XXY (Klein 1,1, Y-par 1) + YYY (Klein 1,1, Y-par 1) — both same on all axes
    H_full_homo = fw.PauliHamiltonian.from_letter_tuples(
        [('X', 'X', 'Y'), ('Y', 'Y', 'Y')], chain_length=4
    )
    assert H_full_homo.is_klein_homogeneous is True
    assert H_full_homo.is_y_parity_homogeneous is True
    assert H_full_homo.is_z2_homogeneous is True


def test_PauliHamiltonian_per_term_properties():
    H = fw.PauliHamiltonian.from_letter_tuples(
        [('X', 'Y'), ('Y', 'Z'), ('X', 'X')], chain_length=3
    )
    assert H.per_term_klein_indices == [(0, 1), (1, 0), (0, 0)]
    assert H.per_term_pi2_classes == ['pi2_odd', 'pi2_even_nontruly', 'truly']
    assert H.per_term_y_parities == [1, 1, 0]
    assert H.has_truly_term is True
    assert H.is_klein_homogeneous is False


def test_PauliHamiltonian_chain_length_validation():
    """Term body count must not exceed chain length."""
    with pytest.raises(ValueError):
        fw.PauliHamiltonian.from_letter_tuples([('X', 'Y', 'Z')], chain_length=2)


def test_PauliHamiltonian_klein_homogeneity_predicts_F77_at_k2():
    """Empirical structural fact: at k=2, Klein-homogeneous Hamiltonians on
    a chain with Z-dephasing are always F77 truly or soft, never hard.

    This verifies the rule on a chain via the F77 classifier."""
    chain = fw.ChainSystem(N=3)

    # Klein-homogeneous test cases
    homogeneous = [
        [('X', 'X'), ('Y', 'Y')],   # M+M
        [('X', 'Y'), ('Y', 'X')],   # F_a + F_a
        [('X', 'Z'), ('Z', 'X')],   # F_b + F_b
        [('Y', 'Z'), ('Z', 'Y')],   # C + C
    ]
    for terms in homogeneous:
        H = fw.PauliHamiltonian.from_letter_tuples(terms, chain_length=3)
        assert H.is_klein_homogeneous is True
        f77 = fw.classify_pauli_pair(chain, terms)
        assert f77 in ('truly', 'soft'), \
            f"Klein-homogeneous {terms} should be soft or truly, got {f77}"

    # Klein-inhomogeneous: can be hard
    H_hard = fw.PauliHamiltonian.from_letter_tuples([('X', 'Y'), ('X', 'Z')], chain_length=3)
    assert H_hard.is_klein_homogeneous is False
    assert fw.classify_pauli_pair(chain, [('X', 'Y'), ('X', 'Z')]) == 'hard'
