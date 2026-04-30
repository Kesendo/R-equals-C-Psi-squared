"""Receiver entity: state-bearing class for receiver-engineering analysis."""
from __future__ import annotations

import numpy as np

from .symmetry import f71_eigenstate_class, receiver_engineering_signature


class Receiver:
    """A quantum state vector with framework-aware F71 classification.

    Args:
        psi: complex array of length 2^N
        chain: optional ChainSystem (must match N if given)

    Properties (cached):
        N, f71_class, rho

    Methods:
        signature() → receiver_engineering_signature dict
    """

    def __init__(self, psi, chain=None, atol=1e-6):
        psi_arr = np.asarray(psi, dtype=complex)
        if psi_arr.ndim != 1:
            raise ValueError(
                f"psi must be a 1D state vector; got shape {psi_arr.shape}. "
                f"To wrap a density matrix, use Receiver.from_rho(rho, chain) instead."
            )
        self.psi = psi_arr
        self.N = int(round(np.log2(len(self.psi))))
        if 2 ** self.N != len(self.psi):
            raise ValueError(f"psi length {len(self.psi)} is not a power of 2")
        norm = float(np.linalg.norm(self.psi))
        if not np.isclose(norm, 1.0, atol=atol):
            raise ValueError(
                f"psi must be normalized (||psi||=1); got ||psi||={norm:.6g}. "
                f"Either normalize before passing or use Receiver.from_psi_unnormalized(...)."
            )
        if chain is not None and chain.N != self.N:
            raise ValueError(f"chain.N ({chain.N}) does not match psi N ({self.N})")
        self.chain = chain
        self._f71_class_cache = "unset"
        self._rho_cache = None

    @classmethod
    def from_psi_unnormalized(cls, psi, chain=None):
        """Wrap a non-normalized state vector by normalizing first."""
        psi = np.asarray(psi, dtype=complex).ravel()
        n = float(np.linalg.norm(psi))
        if n == 0.0:
            raise ValueError("psi is the zero vector; cannot normalize.")
        return cls(psi / n, chain=chain)

    @property
    def f71_class(self):
        if self._f71_class_cache == "unset":
            self._f71_class_cache = f71_eigenstate_class(self.psi)
        return self._f71_class_cache

    @property
    def rho(self):
        if self._rho_cache is None:
            rho = np.outer(self.psi, self.psi.conj())
            self._rho_cache = (rho + rho.conj().T) / 2.0
        return self._rho_cache

    def signature(self):
        return receiver_engineering_signature(self.psi)
