"""Tests for F83 closed-form anti-fraction and the predict_pi_decomposition primitive."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_F83_pi_decomposition_anti_fraction_closed_form():
    """F83: anti-fraction = 1/(2 + 4·r) where r = ‖H_even_nontruly‖²/‖H_odd‖².

    Verifies the F83 closed form across the trichotomy + mixed configurations
    at N=3 and N=4. The predicted anti-fraction matches the numerical
    M_anti²/M² ratio bit-exact.
    """
    test_cases = [
        # (label, terms, expected anti-fraction)
        ('pure odd XY+YX', [('X', 'Y'), ('Y', 'X')], 0.5),
        ('pure odd XY', [('X', 'Y')], 0.5),
        ('pure odd XZ', [('X', 'Z')], 0.5),
        ('pure even non-truly YZ+ZY', [('Y', 'Z'), ('Z', 'Y')], 0.0),
        ('pure even non-truly YZ', [('Y', 'Z')], 0.0),
        ('truly XX+YY (M=0)', [('X', 'X'), ('Y', 'Y')], 0.0),  # limiting value
        ('mixed equal XY+YZ', [('X', 'Y'), ('Y', 'Z')], 1.0 / 6),
        ('mixed equal XY+ZY', [('X', 'Y'), ('Z', 'Y')], 1.0 / 6),
        ('mixed equal XZ+YZ', [('X', 'Z'), ('Y', 'Z')], 1.0 / 6),
        ('asymmetric more odd XY+YX+YZ', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z')], 1.0 / 4),
        ('full mix XY+YX+YZ+ZY', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z'), ('Z', 'Y')], 1.0 / 6),
        ('truly+mixed XX+XY+YZ', [('X', 'X'), ('X', 'Y'), ('Y', 'Z')], 1.0 / 6),
    ]
    for N in [3, 4, 5]:
        chain = fw.ChainSystem(N=N)
        for label, terms, expected_anti in test_cases:
            predicted = fw.predict_pi_decomposition_anti_fraction(chain,terms)
            assert abs(predicted - expected_anti) < 1e-10, \
                f"N={N} {label}: predicted={predicted}, expected={expected_anti}"

            # Cross-check against numerical pi_decompose_M
            d = fw.pi_decompose_M(chain,terms, gamma_z=0.0)
            if d['norm_sq']['M'] < 1e-12:
                # truly H: M = 0, anti-fraction undefined (closed form returns 0)
                continue
            measured = d['norm_sq']['M_anti'] / d['norm_sq']['M']
            assert abs(measured - expected_anti) < 1e-9, \
                f"N={N} {label}: numerical anti-fraction={measured}, expected={expected_anti}"
            assert abs(predicted - measured) < 1e-9, \
                f"N={N} {label}: predicted closed form ≠ numerical evaluation"

    # γ-independence: predicted anti-fraction should not depend on γ_z
    chain3 = fw.ChainSystem(N=3)
    for gz in [0.0, 0.05, 0.5]:
        d = fw.pi_decompose_M(chain3,[('X', 'Y'), ('Y', 'Z')], gamma_z=gz)
        if d['norm_sq']['M'] > 1e-12:
            measured = d['norm_sq']['M_anti'] / d['norm_sq']['M']
            assert abs(measured - 1.0 / 6) < 1e-9, \
                f"γ_z={gz}: anti-fraction should be γ-independent, got {measured}"


def test_F83_predict_pi_decomposition_full_closed_form():
    """F83 forward primitive: predict_pi_decomposition returns all closed-form
    norms (‖M‖², ‖M_anti‖², ‖M_sym‖²) plus anti-fraction from H alone. The
    predictions must match the numerical pi_decompose_M output bit-exact.
    """
    test_cases = [
        ('XY+YX (pure odd)', [('X', 'Y'), ('Y', 'X')]),
        ('YZ+ZY (pure even non-truly)', [('Y', 'Z'), ('Z', 'Y')]),
        ('XY+YZ (mixed equal)', [('X', 'Y'), ('Y', 'Z')]),
        ('XY+YX+YZ (asymmetric)', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z')]),
        ('XY+YX+YZ+ZY (full mix)', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z'), ('Z', 'Y')]),
        ('XX+XY+YZ (truly + mixed)', [('X', 'X'), ('X', 'Y'), ('Y', 'Z')]),
    ]
    for N in [3, 4, 5]:
        chain = fw.ChainSystem(N=N)
        for label, terms in test_cases:
            pred = fw.predict_pi_decomposition(chain,terms)
            num = fw.pi_decompose_M(chain,terms, gamma_z=0.0)

            # All three norm-squared predictions match numerical
            assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
                f"N={N} {label} ‖M‖²: pred={pred['M_sq']}, num={num['norm_sq']['M']}"
            assert abs(pred['M_anti_sq'] - num['norm_sq']['M_anti']) < 1e-9, \
                f"N={N} {label} ‖M_anti‖²: pred={pred['M_anti_sq']}, num={num['norm_sq']['M_anti']}"
            assert abs(pred['M_sym_sq'] - num['norm_sq']['M_sym']) < 1e-9, \
                f"N={N} {label} ‖M_sym‖²: pred={pred['M_sym_sq']}, num={num['norm_sq']['M_sym']}"

            # anti_fraction matches the convenience wrapper
            wrapper = fw.predict_pi_decomposition_anti_fraction(chain,terms)
            assert abs(pred['anti_fraction'] - wrapper) < 1e-15

            # Pythagoras: M_anti_sq + M_sym_sq = M_sq
            assert abs(pred['M_anti_sq'] + pred['M_sym_sq'] - pred['M_sq']) < 1e-9

    # Special-case ratios at N=3 for canonical r values
    chain3 = fw.ChainSystem(N=3)
    pred_pure_odd = fw.predict_pi_decomposition(chain3,[('X', 'Y'), ('Y', 'X')])
    assert pred_pure_odd['r'] == 0.0
    assert abs(pred_pure_odd['anti_fraction'] - 0.5) < 1e-15

    pred_pure_even = fw.predict_pi_decomposition(chain3,[('Y', 'Z'), ('Z', 'Y')])
    assert pred_pure_even['r'] == float('inf')
    assert pred_pure_even['anti_fraction'] == 0.0
    assert pred_pure_even['M_anti_sq'] == 0.0

    pred_equal_mix = fw.predict_pi_decomposition(chain3,[('X', 'Y'), ('Y', 'Z')])
    assert abs(pred_equal_mix['r'] - 1.0) < 1e-15
    assert abs(pred_equal_mix['anti_fraction'] - 1.0/6) < 1e-15

    # All-truly: M=0, anti_fraction defaults to 0
    pred_truly = fw.predict_pi_decomposition(chain3,[('X', 'X'), ('Y', 'Y')])
    assert pred_truly['M_sq'] == 0.0
    assert pred_truly['anti_fraction'] == 0.0


def test_F83_topology_generalization():
    """F83 closed form is topology-independent: ring, star, complete K_N
    give the same closed-form prediction as numerical pi_decompose_M.

    The matrix-based F83 primitive (post-review fix) builds H_odd and
    H_even_nontruly via _build_bilinear which respects the chosen
    topology's bond graph. This test verifies that the closed form
    matches numerical L → M → pi_decompose_M at N=4 for each non-chain
    topology.
    """
    test_cases = [
        ('XY+YX (pure odd)', [('X', 'Y'), ('Y', 'X')]),
        ('YZ+ZY (pure even non-truly)', [('Y', 'Z'), ('Z', 'Y')]),
        ('XY+YZ (mixed equal)', [('X', 'Y'), ('Y', 'Z')]),
        ('XX+XY+YZ (truly + mixed)', [('X', 'X'), ('X', 'Y'), ('Y', 'Z')]),
    ]
    N = 4
    for topology in ['ring', 'star', 'complete']:
        chain = fw.ChainSystem(N=N, topology=topology)
        for label, terms in test_cases:
            pred = fw.predict_pi_decomposition(chain,terms)
            num = fw.pi_decompose_M(chain,terms, gamma_z=0.0)

            # Match closed form to numerical decomposition
            assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
                f"topology={topology} {label} ‖M‖²: pred={pred['M_sq']}, num={num['norm_sq']['M']}"
            assert abs(pred['M_anti_sq'] - num['norm_sq']['M_anti']) < 1e-9, \
                f"topology={topology} {label} ‖M_anti‖²: pred={pred['M_anti_sq']}, num={num['norm_sq']['M_anti']}"
            assert abs(pred['M_sym_sq'] - num['norm_sq']['M_sym']) < 1e-9, \
                f"topology={topology} {label} ‖M_sym‖²: pred={pred['M_sym_sq']}, num={num['norm_sq']['M_sym']}"

            # anti_fraction matches numerical (when M ≠ 0)
            if num['norm_sq']['M'] > 1e-12:
                num_anti_frac = num['norm_sq']['M_anti'] / num['norm_sq']['M']
                assert abs(pred['anti_fraction'] - num_anti_frac) < 1e-9, \
                    f"topology={topology} {label} anti-fraction"
