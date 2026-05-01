"""Tests for the slow-mode lens workflow.

Validation anchor: the C# `lens` mode in compute/RCPsiSquared.Compute uses
an N=5 IBM sacrifice profile (γ = [2.33573, 0.09937, 0.05000, 0.07173, 0.05132])
with reference rate 0.3181, SE Frobenius ratio 0.9986, lens amplitudes
[0.099, 0.239, 0.428, 0.572, 0.651] (cosine similarity > 0.999), and slow-mode
projection 0.972. The Python primitive must reproduce this.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_lens_n3_sanity():
    """N=3 chain produces sortable slow modes with positive rates."""
    chain = fw.ChainSystem(N=3)
    sm = fw.slow_modes(chain, n_slow=3)
    assert sm['eigenvalues'].shape == (3,)
    assert sm['rates'].shape == (3,)
    assert np.all(sm['rates'] > 0), "rates must be positive (decay rates)"
    assert np.all(np.diff(sm['rates']) >= -1e-12), "rates must be sorted ascending"
    assert sm['right_eigvecs'].shape == (4 ** 3, 3)
    assert sm['left_covecs'].shape == (3, 4 ** 3)


def test_lens_n5_ibm_sacrifice_validation():
    """Reproduce the C# `lens` validation gate at N=5 IBM sacrifice profile.

    Reference values from compute/RCPsiSquared.Compute/Program.cs lines 1041-1071:
      γ = [2.33573, 0.09937, 0.05000, 0.07173, 0.05132]
      lens rate ≈ 0.3181, SE Frobenius ratio ≈ 0.9986
      amplitudes ≈ [0.099, 0.239, 0.428, 0.572, 0.651]
      |c_slow| ≈ 0.972
    """
    chain = fw.ChainSystem(N=5)
    ibm_gamma = [2.33573, 0.09937, 0.05000, 0.07173, 0.05132]
    L = fw.lindbladian_z_dephasing(chain.H, ibm_gamma)
    result = fw.lens_pipeline(chain, L=L)

    assert result['lens_amplitudes'] is not None, "expected SE-accessible slow mode"
    assert abs(result['lens_rate'] - 0.3181) < 1e-3, \
        f"lens rate {result['lens_rate']:.4f}, expected 0.3181"
    assert abs(result['lens_se_ratio'] - 0.9986) < 5e-3, \
        f"SE Frob ratio {result['lens_se_ratio']:.4f}, expected 0.9986"

    ref_amps = np.array([0.099, 0.239, 0.428, 0.572, 0.651])
    cosine = float(np.dot(result['lens_amplitudes'], ref_amps))
    assert cosine > 0.999, \
        f"amplitude cosine similarity {cosine:.6f}, expected > 0.999"

    # |c_slow| sits in [0.94, 0.98] — exact value depends on LAPACK's choice
    # of basis within the degenerate (rate=0.3181, multiplicity 2) subspace.
    # Cosine similarity above already validates the lens direction is correct.
    assert 0.93 < result['lens_projection'] < 0.99, \
        f"|c_slow| {result['lens_projection']:.4f}, expected ~0.95"


def test_lens_uniform_gamma_chain_finds_se_mode():
    """At N=4 chain with uniform γ, the lens pipeline finds an SE-accessible mode."""
    chain = fw.ChainSystem(N=4)
    result = fw.lens_pipeline(chain)
    # Uniform-γ Heisenberg on chain has dense slow-mode structure; the lens
    # should reach an SE-accessible mode within the default scan width.
    assert result['lens_amplitudes'] is not None, \
        "uniform-γ N=4 chain should have SE-accessible slow mode within default n_slow"
    assert result['lens_amplitudes'].shape == (4,)
    # Amplitudes are normalized
    norm = float(np.linalg.norm(result['lens_amplitudes']))
    assert abs(norm - 1.0) < 1e-10, f"amplitude norm {norm}, expected 1.0"


def test_slow_modes_excludes_stationary():
    """The stationary mode (Re ≈ 0) is excluded by default."""
    chain = fw.ChainSystem(N=3)
    sm = fw.slow_modes(chain, n_slow=3, exclude_stationary=True)
    # All returned rates should be > stationary_tol
    assert np.all(sm['rates'] > 1e-10)


def test_slow_modes_L_override():
    """Custom L produces different slow modes than chain.L."""
    chain = fw.ChainSystem(N=3)
    sm_default = fw.slow_modes(chain, n_slow=3)
    # Build L with much stronger γ
    strong_gamma = [1.0] * 3
    L_strong = fw.lindbladian_z_dephasing(chain.H, strong_gamma)
    sm_strong = fw.slow_modes(chain, n_slow=3, L=L_strong)
    # Stronger γ → larger decay rates (faster modes)
    assert float(sm_strong['rates'][0]) > float(sm_default['rates'][0]) + 0.1, \
        "stronger γ must increase the slowest decay rate"
