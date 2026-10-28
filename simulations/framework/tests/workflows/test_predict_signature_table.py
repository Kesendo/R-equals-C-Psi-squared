"""Tests for predict_signature_table: the Trotter step must be the hardware's.

The Marrakesh runners (run_soft_break.py for April 26 and April 30,
run_framework_snapshots.py for the April 26 snapshots) build one
SparsePauliOp with the bond loop OUTSIDE the term loop and append one
PauliEvolutionGate per Trotter step; Qiskit's default LieTrotter synthesis
keeps that list order. The gate below compares the framework's step unitary
with the circuit Qiskit actually synthesises (`gate.definition`, not
`Operator(gate)`, which returns the exact exponential and hides the order).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw
from framework.workflows.predict_signature_table import _trotter_step_unitary


CATEGORIES = {
    'truly_unbroken':       [('X', 'X'), ('Y', 'Y')],
    'pi2_odd_pure':         [('X', 'Y'), ('Y', 'X')],
    'pi2_even_nontruly':    [('Y', 'Z'), ('Z', 'Y')],
    'mixed_anti_one_sixth': [('X', 'Y'), ('Y', 'Z')],
    'soft_broken':          [('X', 'Y'), ('Y', 'X')],
    'hard_broken':          [('X', 'X'), ('X', 'Y')],
}


def _runner_step_unitary(N, bonds, terms, J, dt):
    """The step as the hardware runner builds and Qiskit synthesises it,
    returned in the framework's site order (site 0 = leftmost kron factor)."""
    qiskit = pytest.importorskip('qiskit')
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import Operator, SparsePauliOp

    strings, coeffs = [], []
    for (i, j) in bonds:          # run_soft_break.py: bonds outside
        for (la, lb) in terms:    # terms inside
            chars = ['I'] * N
            chars[N - 1 - i] = la
            chars[N - 1 - j] = lb
            strings.append(''.join(chars))
            coeffs.append(J)
    gate = PauliEvolutionGate(SparsePauliOp(strings, coeffs=np.array(coeffs)), time=dt)
    return Operator(gate.definition).reverse_qargs().data


def _terms_outer_step_unitary(N, bonds, terms, J, dt):
    """Terms outside, bonds inside: the wrong order, kept as the negative control."""
    from scipy.linalg import expm
    U = np.eye(2 ** N, dtype=complex)
    for (P, Q) in terms:
        for (l, m) in bonds:
            ops = [fw.ur_pauli('I')] * N
            ops[l] = fw.ur_pauli(P)
            ops[m] = fw.ur_pauli(Q)
            op = ops[0]
            for o in ops[1:]:
                op = np.kron(op, o)
            U = expm(-1j * J * dt * op) @ U
    return U


# The factor sequence Qiskit 2.3.0 synthesises for the runner's XX+YY step on
# N = 3, read from `gate.definition`: one (bond, term) factor per instruction.
RECORDED_TRULY_SEQUENCE = [('rxx', [0, 1]), ('ryy', [0, 1]), ('rxx', [1, 2]), ('ryy', [1, 2])]


def _pauli_string_op(N, l, m, P, Q):
    ops = [fw.ur_pauli('I')] * N
    ops[l] = fw.ur_pauli(P)
    ops[m] = fw.ur_pauli(Q)
    op = ops[0]
    for o in ops[1:]:
        op = np.kron(op, o)
    return op


def test_trotter_step_matches_the_recorded_sequence_without_qiskit():
    """Qiskit-free pin of the order: build the step from the recorded
    sequence with exp(-iθP) = cos θ·I - i sin θ·P (P² = I) and compare.
    Error model as below; the order effect is O(dt²), 0.13 here."""
    chain = fw.ChainSystem(N=3)
    dt = 0.8 / 3
    letter = {'rxx': ('X', 'X'), 'ryy': ('Y', 'Y')}
    U = np.eye(8, dtype=complex)
    for name, (l, m) in RECORDED_TRULY_SEQUENCE:
        P = _pauli_string_op(3, l, m, *letter[name])
        U = (np.cos(chain.J * dt) * np.eye(8) - 1j * np.sin(chain.J * dt) * P) @ U
    U_tool = _trotter_step_unitary(3, chain.bonds, CATEGORIES['truly_unbroken'], chain.J, dt)
    bound = 64 * len(RECORDED_TRULY_SEQUENCE) * np.finfo(float).eps
    assert np.max(np.abs(U_tool - U)) < bound


def test_recorded_sequence_is_what_qiskit_synthesises():
    pytest.importorskip('qiskit')
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import SparsePauliOp
    gate = PauliEvolutionGate(SparsePauliOp(['IXX', 'IYY', 'XXI', 'YYI']), time=0.8 / 3)
    seq = [(ins.operation.name, [gate.definition.find_bit(q).index for q in ins.qubits])
           for ins in gate.definition.data]
    assert seq == RECORDED_TRULY_SEQUENCE


@pytest.mark.parametrize('dt', [0.01, 0.8 / 3, 1.7])
@pytest.mark.parametrize('name', sorted(CATEGORIES))
def test_trotter_step_matches_the_runner_circuit(name, dt):
    """Error model: both sides are products of n_f = len(bonds)·len(terms)
    one-string exponentials, each correct to a few ulp, so the residual is
    O(n_f · eps) entry-wise, independent of dt. The gate is 64·n_f·eps across
    three decades of dt; the measured ratio to n_f·eps stays at or below
    1.125. The order effect it must catch is O(dt²), 0.13 at dt = 0.8/3."""
    chain = fw.ChainSystem(N=3)
    terms = CATEGORIES[name]
    U_tool = _trotter_step_unitary(chain.N, chain.bonds, terms, chain.J, dt)
    U_hw = _runner_step_unitary(chain.N, chain.bonds, terms, chain.J, dt)
    n_f = len(chain.bonds) * len(terms)
    bound = 64 * n_f * np.finfo(float).eps
    residual = np.max(np.abs(U_tool - U_hw))
    assert residual < bound, f"{name}: residual {residual:.3e} vs bound {bound:.3e}"


def test_gate_can_fail_on_terms_outer():
    """Negative control: the old terms-outer step must miss the runner's
    circuit wherever the one pair the two orders swap (the second term on
    bond (0,1), the first on bond (1,2)) anticommutes (truly, mixed, hard),
    and must agree where it commutes (pi2_odd, pi2_even)."""
    chain = fw.ChainSystem(N=3)
    dt = 0.8 / 3
    for name, terms in CATEGORIES.items():
        U_old = _terms_outer_step_unitary(chain.N, chain.bonds, terms, chain.J, dt)
        U_hw = _runner_step_unitary(chain.N, chain.bonds, terms, chain.J, dt)
        residual = np.max(np.abs(U_old - U_hw))
        if name in ('truly_unbroken', 'mixed_anti_one_sixth', 'hard_broken'):
            assert residual > 1e-2, f"{name}: terms-outer should differ, residual {residual:.3e}"
        else:
            assert residual < 1e-12, f"{name}: swapped pair commutes, residual {residual:.3e}"
