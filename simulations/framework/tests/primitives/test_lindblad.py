"""Tests for Lindblad primitives (palindrome residual norm scaling formulas)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_palindrome_residual_norm_ratio_squared_n3_n4():
    # main: 4·k/(k-1) at k=3 → 6
    assert fw.palindrome_residual_norm_ratio_squared(3, 4, 'main') == 6.0
    # single-body: 4·(2k-1)/(2k-3) at k=3 → 20/3
    assert abs(fw.palindrome_residual_norm_ratio_squared(3, 4, 'single_body') - 20 / 3) < 1e-12


def test_palindrome_residual_norm_squared_factor_graph_topologies():
    # Chain N=5: B=4, D2=14
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 4, 14, 'main') == 4 * 64
    # Ring N=5: B=5, D2=20
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 5, 20, 'main') == 5 * 64
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 5, 20, 'single_body') == 10 * 64
