"""Tests for F85 k-body trichotomy classifier and predict_pi_decomposition / F81 generalization."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_F85_kbody_trichotomy_counts():
    """F85: trichotomy enumeration at k = 2, 3, 4 over {X,Y,Z}^k.

    Verifies the closed-form counts:
      Π²-odd:           (3^k − (−1)^k) / 2
      Π²-even non-truly: pure-letter triples that use {Y, Z} only
      truly: rest

    Empirical: k=2 → 3/4/2, k=3 → 7/14/6, k=4 → 21/40/20.
    """
    from itertools import product
    from collections import Counter
    from framework.core import _pauli_tuple_is_truly, _pauli_tuple_pi2_class

    expected_counts = {
        2: {'truly': 3, 'pi2_odd': 4, 'pi2_even_nontruly': 2},
        3: {'truly': 7, 'pi2_odd': 14, 'pi2_even_nontruly': 6},
        4: {'truly': 21, 'pi2_odd': 40, 'pi2_even_nontruly': 20},
    }
    for k, expected in expected_counts.items():
        counts = Counter()
        for letters in product('XYZ', repeat=k):
            counts[_pauli_tuple_pi2_class(letters)] += 1
        assert dict(counts) == expected, \
            f"k={k} trichotomy counts mismatch: got {dict(counts)}, expected {expected}"
        # Closed form for Π²-odd count: (3^k − (−1)^k) / 2
        expected_odd = (3**k - (-1)**k) // 2
        assert counts['pi2_odd'] == expected_odd

    # Backward compat: 2-body classification matches _pauli_pair_is_truly
    from framework.core import _pauli_pair_is_truly
    for a, b in product('IXYZ', repeat=2):
        pair_check = _pauli_pair_is_truly(a, b)
        tuple_check = _pauli_tuple_is_truly((a, b))
        assert pair_check == tuple_check, f"({a},{b}): pair={pair_check}, tuple={tuple_check}"


def test_F85_kbody_predict_pi_decomposition():
    """F85: predict_pi_decomposition matches numerical pi_decompose_M at k=3, 4.

    Verifies that the F83 anti-fraction closed form generalizes to k-body
    via the Π²-class trichotomy (truly / Π²-odd / Π²-even non-truly), with
    the c(k) ∈ {0, 1, 2} factor scheme.
    """
    test_cases_k3 = [
        ('XYZ (Π²-even non-truly all-3-distinct)', [('X', 'Y', 'Z')], 0.0),
        ('XXY (Π²-odd)', [('X', 'X', 'Y')], 0.5),
        ('YYY (Π²-odd, n_YZ=3)', [('Y', 'Y', 'Y')], 0.5),
        ('XYY (truly)', [('X', 'Y', 'Y')], 0.0),
        ('XYZ + XXY (mixed)', [('X', 'Y', 'Z'), ('X', 'X', 'Y')], None),  # auto-check vs numerical
    ]
    for N in [4, 5]:
        chain = fw.ChainSystem(N=N)
        for label, terms, expected_anti in test_cases_k3:
            pred = chain.predict_pi_decomposition(terms)
            num = chain.pi_decompose_M(terms, gamma_z=0.0)
            # Closed form matches numerical bit-exact
            assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
                f"N={N} {label}: M_sq pred={pred['M_sq']}, num={num['norm_sq']['M']}"
            assert abs(pred['M_anti_sq'] - num['norm_sq']['M_anti']) < 1e-9
            assert abs(pred['M_sym_sq'] - num['norm_sq']['M_sym']) < 1e-9
            if expected_anti is not None and pred['M_sq'] > 1e-12:
                assert abs(pred['anti_fraction'] - expected_anti) < 1e-10, \
                    f"N={N} {label}: anti={pred['anti_fraction']}, expected {expected_anti}"

    # 4-body verification: predict matches numerical
    chain4 = fw.ChainSystem(N=5)
    test_cases_k4 = [
        [('X', 'Y', 'Z', 'X')],         # contains all three letters → Π²-even non-truly
        [('X', 'X', 'Y', 'Y')],         # 2X 2Y arrangement → truly
        [('Y', 'Y', 'Y', 'Z')],         # 3Y 1Z arrangement → Π²-even non-truly
        [('X', 'X', 'X', 'Y')],         # 3X 1Y arrangement → Π²-odd
    ]
    for terms in test_cases_k4:
        pred = chain4.predict_pi_decomposition(terms)
        num = chain4.pi_decompose_M(terms, gamma_z=0.0)
        assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
            f"k=4 {terms}: M_sq pred {pred['M_sq']} vs num {num['norm_sq']['M']}"


def test_F85_kbody_classifier_at_k5_spot_check():
    """F85 truly criterion verified at k=5 via spot check.

    Full 243-tuple enumeration at k=5 is computationally expensive (would
    add several minutes to the suite); instead we spot-check 8 representative
    5-tuples covering all three Π²-classes plus various letter compositions.

    Together with the full enumeration at k=2, 3, 4 in
    `test_F85_kbody_trichotomy_counts`, this provides empirical evidence
    that the analytic rule "#Y even AND #Z even ⟹ truly" continues to hold
    at k=5. The closed form for the Π²-odd count (3^k − (−1)^k)/2 = 122 at
    k=5 is verifiable analytically (binomial generating function).
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_kbody_chain
    from framework.core import _pauli_tuple_pi2_class

    N = 5
    test_cases = [
        ('XXXXX', 'truly'),         # all X (#Y=0 even, #Z=0 even)
        ('YYYYY', 'pi2_odd'),       # all Y (#Y=5 odd → bit_b=1)
        ('ZZZZZ', 'pi2_odd'),       # all Z (#Z=5 odd → bit_b=1)
        ('XXXYZ', 'pi2_even_nontruly'),  # bit_b=0 but #Z=1 odd
        ('YYZZX', 'truly'),         # #Y=2 even, #Z=2 even
        ('XYZXY', 'pi2_odd'),       # bit_b=1 odd
        ('YYZZY', 'pi2_odd'),       # bit_b=1 (3Y + 2Z = 5 mod 2 = 1)
        ('XZYXZ', 'pi2_odd'),       # bit_b=1 (1Y + 2Z = 3 mod 2 = 1)
    ]
    for label, expected in test_cases:
        letters = tuple(label)
        cls_ana = _pauli_tuple_pi2_class(letters)
        assert cls_ana == expected, f'{label}: analytic={cls_ana}, expected {expected}'

        # Numerical verification
        H = _build_kbody_chain(N, [letters + (1.0,)])
        H_sq = float(np.real(np.trace(H.conj().T @ H)))
        L = lindbladian_z_dephasing(H, [0.0] * N)
        M = palindrome_residual(L, 0.0, N)
        M_sq = float(np.linalg.norm(M) ** 2)
        if M_sq < 1e-10:
            cls_num = 'truly'
        else:
            ratio = M_sq / (H_sq * 2 ** N)
            if abs(ratio - 4) < 1e-6:
                cls_num = 'pi2_odd'
            elif abs(ratio - 8) < 1e-6:
                cls_num = 'pi2_even_nontruly'
            else:
                cls_num = f'unexpected ratio {ratio}'
        assert cls_ana == cls_num, f'{label}: analytic={cls_ana}, numerical={cls_num}'


def test_F85_kbody_F81_identity_at_k3():
    """F85: F81 identity Π·M·Π⁻¹ = M − 2·L_{H_odd} generalizes verbatim to k-body.

    Tests at k=3 by computing M and Π·M·Π⁻¹ directly and verifying the
    F81 identity holds bit-exact for Π²-odd 3-body terms.
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _vec_to_pauli_basis_transform, _build_kbody_chain
    from framework.symmetry import build_pi_full

    N = 4
    chain = fw.ChainSystem(N=N)
    Pi = build_pi_full(N)
    Pi_inv = np.linalg.inv(Pi)
    T = _vec_to_pauli_basis_transform(N)
    d = 2 ** N

    # 3-body Π²-odd: H_odd contributes via L_{H_odd}, F81 holds
    H = _build_kbody_chain(N, [('X', 'X', 'Y', 1.0)])
    L = lindbladian_z_dephasing(H, [0.0]*N)
    M = palindrome_residual(L, 0.0, N)
    PiMPi = Pi @ M @ Pi_inv
    L_H_vec = -1j * (np.kron(H, np.eye(d, dtype=complex)) -
                     np.kron(np.eye(d, dtype=complex), H.T))
    L_H = (T.conj().T @ L_H_vec @ T) / d
    # F81: Π·M·Π⁻¹ = M − 2·L_H_odd (and L_H_odd = L_H since term is purely Π²-odd)
    diff = np.linalg.norm(PiMPi - (M - 2 * L_H))
    assert diff < 1e-9, f"k=3 F81 identity fails: ‖diff‖ = {diff}"
