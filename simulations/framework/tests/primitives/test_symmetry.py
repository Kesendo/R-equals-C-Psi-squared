"""Tests for symmetry primitives (chain mirror, bond mirror basis)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_chain_mirror_state_involutory():
    for N in [3, 4, 5, 6]:
        R = fw.chain_mirror_state(N)
        np.testing.assert_allclose(R @ R, np.eye(2 ** N), atol=1e-12)


def test_bond_mirror_basis_dimensions():
    """Verify the parity-tied dimensional structure."""
    for N, expected_sym, expected_asym in [(3, 1, 1), (4, 2, 1), (5, 2, 2),
                                             (6, 3, 2), (7, 3, 3), (8, 4, 3)]:
        sym, asym = fw.bond_mirror_basis(N)
        assert (len(sym), len(asym)) == (expected_sym, expected_asym), \
            f"N={N}: got ({len(sym)}, {len(asym)}), expected ({expected_sym}, {expected_asym})"
