"""Tests for fw.chiral_panel: when conjugation by K is a symmetry of L, and what it protects.

Under Z dephasing, conjugation by K = prod_{odd i} Z_i keeps the dissipator and sends
-i[H, .] to -i[KHK, .]: a symmetry of L exactly for K-even H. A K-even state then stays
K-even, so the K-odd Paulis keep expectation zero at every time.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
P1 = {"I": I2, "X": X, "Y": Y, "Z": Z}


def op(letters):
    out = np.array([[1.0]], dtype=complex)
    for c in letters:
        out = np.kron(out, P1[c])
    return out


def evolve(H, rho0, gamma, t):
    d = H.shape[0]; Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))          # row-major vec
    n = int(np.log2(d))
    for l in range(n):
        Zl = op("I" * l + "Z" + "I" * (n - l - 1))
        L += gamma * (np.kron(Zl, Zl.T) - np.kron(Id, Id))
    return (expm(t * L) @ rho0.reshape(-1)).reshape(d, d)


def test_k_odd_hamiltonian_is_not_a_symmetry_and_protects_nothing():
    """XY chain: KHK = -H. The panel must not claim protection; the control shows why:
    a K-odd Pauli of a K-even state leaves zero."""
    H = op("XXI") + op("YYI") + op("IXX") + op("IYY")
    psi = np.zeros(8, dtype=complex); psi[0b100] = 1; psi[0b001] = 1j; psi /= np.linalg.norm(psi)
    rho0 = np.outer(psi, psi.conj())
    out = fw.chiral_panel(H, rho0, 3)
    assert out["K_status"] == "K-odd"
    assert out["K_symmetric_L"] is False
    assert out["chiral_protected_dynamic"] == []
    rho = evolve(H, rho0, 0.5, 0.3)
    assert abs(np.trace(op("IXX") @ rho)) > 0.1                 # 0.29: K-odd, nonzero at t > 0


def test_k_even_hamiltonian_protects_exactly_the_k_odd_paulis():
    """H = X0 X2 + Z0 Z1 + Z1 commutes with K = Z_1; psi = |000> + |101> is K-even.
    Every listed Pauli is K-odd and stays at zero along the motion; a K-even one moves."""
    H = op("XIX") + op("ZZI") + op("IZI")
    psi = np.zeros(8, dtype=complex); psi[0b000] = 1; psi[0b101] = 1; psi /= np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())
    out = fw.chiral_panel(H, rho0, 3)
    assert out["K_status"] == "K-even" and out["K_symmetric_L"] is True
    prot = out["chiral_protected_dynamic"]
    K = fw.chiral_K_full(3)
    assert len(prot) == 32                                       # the K-odd strings: half of all 64
    for lab in prot:
        assert np.allclose(K @ op(lab) @ K, -op(lab))
    for t in (0.2, 0.7, 1.9):
        rho = evolve(H, rho0, 0.3, t)
        # expm in floating point: the symmetry makes these zero exactly, the error is rounding (~1e-16)
        assert max(abs(np.trace(op(lab) @ rho)) for lab in prot) < 1e-12
    assert abs(np.trace(op("XIX") @ evolve(H, rho0, 0.3, 0.7))) > 0.1   # a K-even Pauli is not protected
