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


@pytest.mark.parametrize("name, terms", [
    ("Heisenberg", [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]),
    ("XX_only", [('X', 'X')]),
    ("YZ_plus_ZY", [('Y', 'Z'), ('Z', 'Y')]),
    ("XY_pure", [('X', 'Y')]),
    ("XY_plus_YX", [('X', 'Y'), ('Y', 'X')]),
    ("Heisenberg_plus_XY", [('X', 'X'), ('Y', 'Y'), ('Z', 'Z'), ('X', 'Y')]),
])
def test_orthogonality_invariant_across_H_families(name, terms):
    """Frobenius-orthogonality invariant: ||M||^2 = ||M_zero||^2 + ||M_plus||^2 + ||M_minus||^2
    holds bit-exact across all bilinear H families under pure Z-dephasing."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_coordinates(chain, terms, gamma_z=0.05)
    assert result['orthogonality_residual'] < 1e-10, \
        f"{name}: orthogonality residual = {result['orthogonality_residual']}"


@pytest.mark.parametrize("name, terms", [
    ("XX_only", [('X', 'X')]),
    ("YZ_plus_ZY", [('Y', 'Z'), ('Z', 'Y')]),
    ("XY_pure", [('X', 'Y')]),
    ("XY_plus_YX", [('X', 'Y'), ('Y', 'X')]),
    ("Heisenberg", [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]),
    ("Heisenberg_plus_XY", [('X', 'X'), ('Y', 'Y'), ('Z', 'Z'), ('X', 'Y')]),
])
def test_hermitian_balance_plus_minus_half(name, terms):
    """Hermitian H + pure Z-dephasing: ||M_plus||^2 = ||M_minus||^2 (complex-conjugate
    symmetric within M_anti). Working hypothesis from spec; if the test fails for any
    family, that failure IS the discovery (update reflection doc in Task D)."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_coordinates(chain, terms, gamma_z=0.05)
    plus_norm_sq = result['norm_sq']['M_plus_half']
    minus_norm_sq = result['norm_sq']['M_minus_half']
    assert abs(plus_norm_sq - minus_norm_sq) < 1e-10, \
        f"{name}: |plus - minus| = {abs(plus_norm_sq - minus_norm_sq)}, " \
        f"plus = {plus_norm_sq}, minus = {minus_norm_sq}"


def test_t1_cooling_breaks_plus_minus_balance():
    """T1 cooling-only (gamma_t1 > 0, gamma_pump = 0) is expected to break the
    Hermitian-H balance: ||M_plus||^2 != ||M_minus||^2 measurably (F81 violation
    per F84's formula). If this test FAILS (asymmetry < 1e-3), DO NOT silently
    adjust the threshold; the failure is the discovery and Task D documents it."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_coordinates(
        chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')],
        gamma_z=0.05, gamma_t1=0.1, gamma_pump=0.0,
        strict=False  # F81 is no longer exact under T1
    )
    plus = result['norm_sq']['M_plus_half']
    minus = result['norm_sq']['M_minus_half']
    asymmetry = abs(plus - minus)
    assert asymmetry > 1e-3, \
        f"T1 cooling expected to break +/-1/2 balance; asymmetry = {asymmetry}, " \
        f"plus = {plus}, minus = {minus}"
