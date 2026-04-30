"""ChainSystem entity: encapsulates a quantum chain at fixed (N, γ, J, topology, H_type)."""
from __future__ import annotations

import warnings

import numpy as np

from .pauli import pauli_matrix, _site_op_kron


class ChainSystem:
    """Encapsulates a quantum chain at fixed (N, γ, J, topology, H_type).

    Caches Hamiltonian and Liouvillian on first access. Pure entity:
    holds (N, γ_0, J, topology, H_type, bonds, B, degrees, D2) plus
    cached H and L. All analysis logic lives as free functions in
    `framework.diagnostics` and `framework.workflows`.

    Args:
        N: number of qubits.
        gamma_0: uniform Z-dephasing rate per site (default 0.05).
        J: bond coupling (uniform; default 1.0).
        topology: 'chain' (default), 'ring', 'star', 'complete'.
        H_type: 'heisenberg' (XX+YY+ZZ) or 'xy' (XX+YY scaled by J/2).

    F-theorem readings live as free functions in `framework.diagnostics`
    (e.g. `fw.classify_pauli_pair(chain, terms)`,
    `fw.predict_residual_norm_squared_from_terms(chain, terms, gamma_t1=...)`,
    `fw.pi_decompose_M(chain, terms, ...)`). Workflow flows that
    propagate ρ(t), run the cockpit panel, or hold cusp/probe state live
    in `framework.workflows` (e.g. `fw.cockpit_panel(chain, receiver, ...)`,
    `fw.propagate_with_hardware_noise(chain, ρ_0, t, ...)`).
    """

    _FROZEN = frozenset({'N', 'd', 'd2', 'gamma_0', 'J', 'topology', 'H_type',
                          'bonds', 'B', 'degrees', 'D2'})

    def __init__(self, N, gamma_0=0.05, J=1.0, topology='chain', H_type='heisenberg'):
        if N < 2:
            raise ValueError(f"N must be >= 2; got {N}")
        if H_type not in ('heisenberg', 'xy'):
            raise ValueError(f"H_type must be 'heisenberg' or 'xy'; got {H_type!r}")
        if N == 2:
            warnings.warn(
                "ChainSystem(N=2) is mathematically valid but structurally "
                "degenerate: F71 bond_mirror_basis has only 1 sym mode and 0 "
                "asym modes (every F71-eigenstate is capacity-suboptimal); "
                "only 1 bond exists; the drop=28 hardware anchor is "
                "unreproducible (4^2=16 Pauli strings total). Fundamental "
                "vocabulary (classify, Π, palindrome residual, Frobenius "
                "scaling) holds; structural vocabulary (F71-balance, multi-"
                "bond cockpit signatures) needs N>=3.",
                UserWarning, stacklevel=2,
            )
        # bypass __setattr__ guard during init
        object.__setattr__(self, '_initialized', False)
        self.N = N
        self.d = 2 ** N
        self.d2 = self.d * self.d
        self.gamma_0 = float(gamma_0)
        self.J = float(J)
        self.topology = topology
        self.H_type = H_type
        self._build_topology()
        self._H_cache = None
        self._L_cache = None
        object.__setattr__(self, '_initialized', True)

    def __setattr__(self, name, value):
        # Freeze structural attributes after construction. Caches (H, L) are
        # immutable in spirit because the Hamiltonian / Liouvillian depend on
        # all the frozen attrs; allowing J or gamma_0 mutation would silently
        # decouple cache from declared chain state. Make a new ChainSystem
        # instead.
        if getattr(self, '_initialized', False) and name in self._FROZEN:
            raise AttributeError(
                f"ChainSystem.{name} is immutable after construction. "
                f"Make a new ChainSystem(...) with the desired {name}."
            )
        object.__setattr__(self, name, value)

    def _build_topology(self):
        N = self.N
        if self.topology == 'chain':
            self.bonds = [(i, i + 1) for i in range(N - 1)]
        elif self.topology == 'ring':
            self.bonds = [(i, (i + 1) % N) for i in range(N)]
        elif self.topology == 'star':
            self.bonds = [(0, i) for i in range(1, N)]
        elif self.topology == 'complete':
            self.bonds = [(i, j) for i in range(N) for j in range(i + 1, N)]
        else:
            raise ValueError(f"unknown topology: {self.topology!r}")
        self.B = len(self.bonds)
        deg = [0] * N
        for i, j in self.bonds:
            deg[i] += 1
            deg[j] += 1
        self.degrees = deg
        self.D2 = sum(d * d for d in deg)

    @property
    def H(self):
        if self._H_cache is None:
            self._H_cache = self._build_H()
        return self._H_cache

    @property
    def L(self):
        if self._L_cache is None:
            self._L_cache = self._build_L()
        return self._L_cache

    def _build_H(self):
        d = self.d
        H = np.zeros((d, d), dtype=complex)
        Xm, Ym, Zm = pauli_matrix('X'), pauli_matrix('Y'), pauli_matrix('Z')
        for (i, j) in self.bonds:
            xi, xj = _site_op_kron(Xm, i, self.N), _site_op_kron(Xm, j, self.N)
            yi, yj = _site_op_kron(Ym, i, self.N), _site_op_kron(Ym, j, self.N)
            if self.H_type == 'heisenberg':
                zi, zj = _site_op_kron(Zm, i, self.N), _site_op_kron(Zm, j, self.N)
                H = H + self.J * (xi @ xj + yi @ yj + zi @ zj)
            elif self.H_type == 'xy':
                H = H + (self.J / 2.0) * (xi @ xj + yi @ yj)
        return H

    def _build_L(self):
        H = self.H
        d, d2 = self.d, self.d2
        Id = np.eye(d, dtype=complex)
        L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
        Zm = pauli_matrix('Z')
        for k in range(self.N):
            Zk = _site_op_kron(Zm, k, self.N)
            L = L + self.gamma_0 * (np.kron(Zk, Zk.conj()) - np.eye(d2))
        return L
