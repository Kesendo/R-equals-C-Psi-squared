"""Tests for the three-way polarity {−1/2, 0, +1/2} coordinate decomposition.

Smoke + invariant tests for `polarity_coordinates`:
  - dict shape (expected keys present)
  - Frobenius-orthogonality invariant: ‖M‖² = ‖M_zero‖² + ‖M_plus‖² + ‖M_minus‖²
  - Heisenberg (F1 truly) gives M = 0
  - F81 cross-check: F81 M_sym norm² equals polarity M_zero norm²
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_returns_expected_keys():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_coordinates(chain, [('X', 'Y')], gamma_z=0.05)
    for key in ['M', 'M_zero', 'M_plus_half', 'M_minus_half', 'norm_sq', 'asymmetry', 'orthogonality_residual']:
        assert key in result, f"missing key: {key}"


def test_orthogonality_invariant_xy():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_coordinates(chain, [('X', 'Y')], gamma_z=0.05)
    assert result['orthogonality_residual'] < 1e-10, \
        f"orthogonality residual = {result['orthogonality_residual']}"


def test_heisenberg_M_is_zero():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_coordinates(
        chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')], gamma_z=0.05)
    assert result['norm_sq']['M'] < 1e-10


def test_f81_match_M_sym_equals_M_zero():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    f81 = fw.pi_decompose_M(chain, [('X', 'Y')], gamma_z=0.05)
    pol = fw.polarity_coordinates(chain, [('X', 'Y')], gamma_z=0.05)
    assert abs(f81['norm_sq']['M_sym'] - pol['norm_sq']['M_zero']) < 1e-12, \
        f"F81 M_sym ({f81['norm_sq']['M_sym']}) != polarity M_zero ({pol['norm_sq']['M_zero']})"
