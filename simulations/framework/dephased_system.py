"""DephasedSystem entity: one open system as every instrument reads it.

The question it answers: what is THE system, so that two instruments (an exact block
eigensolver, the secular route, later a propagator) read the same one and can be checked
against each other where their ranges overlap?

It holds N, the Hamiltonian and the Z-watching as a real symmetric Kossakowski matrix a: the
site rates gamma_l on the diagonal, the covariances of the site fields off it
(lindbladian_correlated_z_dephasing). The Hamiltonian comes in one of two forms:
  - a dense Hermitian 2^N x 2^N matrix (any H; `from_chain` composes ChainSystem, which stays
    as it is: uniform gamma, frozen);
  - a BOND LIST of a number-conserving XXZ model, H = sum_b J_b (X X + Y Y + delta_b Z Z) over
    the bonds plus sum_l h_l Z_l (`from_bonds`): no 2^N matrix is ever formed, and a sector
    with few excitations is built directly, so the edge blocks are reachable at N = 100.

Every instrument reads the system through `sector(p)` (the popcount-p basis as a z-table and
the Hamiltonian restricted to it), `coherence_price` and `block_liouvillian`, the same way for
both forms.

The price of a coherence |i><j| is the dissipator's diagonal entry 1/2 Delta^T a Delta,
Delta = z(i) - z(j), z_l = +-1 (experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md; uniform
independent gamma gives 2 gamma * Hamming(i, j), the Absorption Theorem's sharp-bit case).

SITE ORDER: big-endian, as lindbladian_z_dephasing: site l is bit N-1-l of the basis index; a
sector lists its states in ascending basis index.
"""
from itertools import combinations
from math import comb

import numpy as np

from .lindblad import lindbladian_correlated_z_dephasing


def _check_watching(a, N):
    a = np.array(a, dtype=float)
    if a.ndim == 1:
        a = np.diag(a)
    if a.shape != (N, N):
        raise ValueError(f"a must be {N}x{N} (or a length-{N} rate vector); got {a.shape}")
    if not np.array_equal(a, a.T):
        raise ValueError("a must be real symmetric (classical site fields)")
    if np.linalg.eigvalsh(a).min() < -1e-12 * max(1.0, np.abs(a).max()):
        raise ValueError("a must be positive semidefinite for a valid generator")
    return a


def _with_covariance(gamma, covariance, N):
    a = np.diag(np.asarray(gamma, dtype=float) * np.ones(N))
    if covariance is not None:
        c = np.array(covariance, dtype=float)
        a = a + c - np.diag(np.diag(c))
    return a


class DephasedSystem:
    """N sites, a Hamiltonian (dense or bonds), Z-dephasing Kossakowski matrix a."""

    def __init__(self, H, a, label=None):
        H = np.asarray(H)
        d = H.shape[0]
        N = int(round(np.log2(d)))
        if H.shape != (d, d) or 2 ** N != d:
            raise ValueError(f"H must be 2^N x 2^N; got shape {H.shape}")
        if not np.allclose(H, H.conj().T):
            raise ValueError("H must be Hermitian")
        self.N = N
        self.H = H
        self.bonds = None
        self.a = _check_watching(a, N)
        self.label = label or f"N={N}"
        self._sectors = {}

    @classmethod
    def from_chain(cls, chain, gamma=None, covariance=None, label=None):
        """The chain's H with per-site rates gamma (default: chain.gamma_0 everywhere) and an
        optional covariance matrix whose off-diagonal entries are added to diag(gamma)."""
        g = chain.gamma_0 if gamma is None else gamma
        return cls(chain.H, _with_covariance(g, covariance, chain.N),
                   label=label or f"{chain.topology} N={chain.N} {chain.H_type}")

    @classmethod
    def from_bonds(cls, N, bonds, J=1.0, delta=1.0, h=None, gamma=0.05, covariance=None, label=None):
        """H = sum_b J_b (X X + Y Y + delta_b Z Z) + sum_l h_l Z_l on the given bonds (pairs of
        sites), with J, delta scalars or one per bond; no dense matrix is built.
        delta = 1 is Heisenberg (ChainSystem's 'heisenberg' at the same J), delta = 0 is XY."""
        self = cls.__new__(cls)
        self.N = N
        self.H = None
        nb = len(bonds)
        self.bonds = [(int(i), int(j), float(Jb), float(db)) for (i, j), Jb, db in
                      zip(bonds, np.broadcast_to(J, nb), np.broadcast_to(delta, nb))]
        self.field = np.zeros(N) if h is None else np.asarray(h, dtype=float) * np.ones(N)
        self.a = _check_watching(_with_covariance(gamma, covariance, N), N)
        self.label = label or f"bonds N={N}"
        self._sectors = {}
        return self

    @staticmethod
    def topology_bonds(N, topology):
        """The bond list of ChainSystem's topologies: chain, ring, star (hub 0), complete."""
        if topology == 'chain':
            return [(i, i + 1) for i in range(N - 1)]
        if topology == 'ring':
            return [(i, i + 1) for i in range(N - 1)] + [(N - 1, 0)]
        if topology == 'star':
            return [(0, i) for i in range(1, N)]
        if topology == 'complete':
            return [(i, j) for i in range(N) for j in range(i + 1, N)]
        raise ValueError(f"unknown topology {topology!r}")

    @property
    def gamma(self):
        """The site rates, the diagonal of a."""
        return np.diag(self.a).copy()

    @property
    def sigma(self):
        """sigma = sum of the site rates, the palindrome centre of the independent case."""
        return float(np.trace(self.a))

    def is_number_conserving(self):
        """True iff H never connects basis states of different popcount (always, for bonds)."""
        if self.H is None:
            return True
        pc = np.array([bin(i).count('1') for i in range(self.H.shape[0])])
        return not np.any(np.abs(self.H[pc[:, None] != pc[None, :]]) > 0)

    def has_xn_symmetry(self):
        """True iff [H, X^N] = 0. Then block (p, q) equals block (N-p, N-q) entry for entry after
        the complement relabelling, at any a (the price 1/2 Delta^T a Delta is even in Delta):
        bonds with no field, or a dense H that X^N maps to itself exactly."""
        if self.H is None:
            return not np.any(self.field)
        return np.array_equal(self.H, self.H[::-1, ::-1])   # X^N reverses the basis index

    def sector_size(self, p):
        return comb(self.N, p)

    def blocks(self):
        """All (p, q) with their row count C(N,p)*C(N,q)."""
        return [((p, q), comb(self.N, p) * comb(self.N, q))
                for p in range(self.N + 1) for q in range(self.N + 1)]

    def sector(self, p):
        """(z, Hp): the popcount-p states as a z-table (rows +-1, ascending basis index) and the
        Hamiltonian restricted to them, dense (small sectors; `sector_sparse` for large ones)."""
        if p in self._sectors:
            return self._sectors[p]
        if self.H is not None:
            z, idx = self._z_table_dense(p)
            out = (z, self.H[np.ix_(idx, idx)])
        else:
            z, Hs = self.sector_sparse(p)
            out = (z, Hs.toarray())
        self._sectors[p] = out
        return out

    def sector_sparse(self, p):
        """(z, Hp) as `sector`, with Hp a scipy CSR matrix; for the bond form no dense matrix of
        the sector is ever formed."""
        import scipy.sparse as sp
        key = ('sparse', p)
        if key in self._sectors:
            return self._sectors[key]
        N = self.N
        if self.H is not None:
            z, idx = self._z_table_dense(p)
            out = (z, sp.csr_matrix(self.H[np.ix_(idx, idx)]))
        else:
            # ascending basis index = descending lexicographic order of the occupied sites' bits;
            # sort by the integer index itself (Python ints, fine past 64 sites)
            occ = sorted(combinations(range(N), p), key=lambda s: sum(1 << (N - 1 - l) for l in s))
            pos = {s: k for k, s in enumerate(occ)}
            n = len(occ)
            z = np.ones((n, N), dtype=int)
            for k, s in enumerate(occ):
                z[k, list(s)] = -1
            diag = (z @ self.field).astype(float)
            rows, cols, vals = [], [], []
            for i, j, Jb, db in self.bonds:
                diag += Jb * db * z[:, i] * z[:, j]
                for k in np.flatnonzero(z[:, i] != z[:, j]):   # one of the two occupied: XX + YY hops it, amplitude 2 J
                    rows.append(pos[tuple(sorted(set(occ[k]) ^ {i, j}))])
                    cols.append(k)
                    vals.append(2.0 * Jb)
            Hs = sp.csr_matrix((vals, (rows, cols)), shape=(n, n)) + sp.diags(diag)
            out = (z, Hs.tocsr())
        self._sectors[key] = out
        return out

    def _z_table_dense(self, p):
        N = self.N
        idx = np.array([i for i in range(2 ** N) if bin(i).count('1') == p])
        bits = (idx[:, None] >> (N - 1 - np.arange(N))[None, :]) & 1
        return 1 - 2 * bits, idx

    def coherence_price(self, zi, zj):
        """Matrix R[i, j] = 1/2 Delta^T a Delta for the coherences between the rows of the z-tables
        zi (kets) and zj (bras): their decay rate under the dissipator alone."""
        zi = np.asarray(zi, dtype=float)
        zj = np.asarray(zj, dtype=float)
        qi = 0.5 * np.einsum('ik,kl,il->i', zi, self.a, zi)
        qj = 0.5 * np.einsum('jk,kl,jl->j', zj, self.a, zj)
        return qi[:, None] + qj[None, :] - zi @ self.a @ zj.T

    def block_liouvillian(self, p, q):
        """The Liouvillian restricted to the (p, q) coherences, row-major vec (ket outer, bra
        inner): -i (H_p x 1 - 1 x H_q^T) - diag(price). Needs a number-conserving H."""
        if not self.is_number_conserving():
            raise ValueError("block_liouvillian needs a number-conserving H")
        zp, Hp = self.sector(p)
        zq, Hq = self.sector(q)
        L = -1j * (np.kron(Hp, np.eye(len(Hq))) - np.kron(np.eye(len(Hp)), Hq.T))
        L -= np.diag(self.coherence_price(zp, zq).ravel())
        return L

    def block_liouvillian_sparse(self, p, q):
        """The same block as `block_liouvillian`, as a scipy CSR matrix (same vec order)."""
        import scipy.sparse as sp
        if not self.is_number_conserving():
            raise ValueError("block_liouvillian_sparse needs a number-conserving H")
        zp, Hp = self.sector_sparse(p)
        zq, Hq = self.sector_sparse(q)
        L = -1j * (sp.kron(Hp, sp.identity(Hq.shape[0]), format="csr")
                   - sp.kron(sp.identity(Hp.shape[0]), Hq.T, format="csr"))
        return (L - sp.diags(self.coherence_price(zp, zq).ravel())).tocsr()

    def dense_hamiltonian(self):
        """The 2^N matrix (built from the bonds when needed; small N only)."""
        if self.H is not None:
            return self.H
        d = 2 ** self.N
        H = np.zeros((d, d))
        for p in range(self.N + 1):
            idx = np.array([i for i in range(d) if bin(i).count('1') == p])
            H[np.ix_(idx, idx)] = self.sector(p)[1]
        return H

    def lindbladian(self):
        """The full 4^N Liouvillian (small N only), vec convention of lindblad.py."""
        return lindbladian_correlated_z_dephasing(self.dense_hamiltonian(), self.a)

    def __repr__(self):
        return f"DephasedSystem({self.label}, sigma={self.sigma:.4g})"
