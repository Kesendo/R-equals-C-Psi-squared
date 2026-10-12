"""slow_rates: how long does this system live, and where?

The first workflow of the orchestra. Several instruments read the same DephasedSystem;
each says which joint-popcount blocks it covers and gives, per block, the slowest nonzero
decay rate with its per-site weights. Where two instruments cover the same block their
readings are shown side by side (the handshake): agreement is the check, disagreement is
the size of the approximation, both reported, neither hidden.

Instruments today:
  ExactBlocks   dense eigensolver per (p, q) block, up to max_rows rows; exact.
  SparseBlocks  the same block sparse, ARPACK at the largest real part (which='LR'); exact up to
                the Krylov tolerance, reaches blocks of a few hundred thousand rows.
  SecularRoute  H eigenbasis + the coherence price 1/2 Delta^T a Delta, Bohr frequencies
                closer than `window` solved together (window 0: exact degeneracy only);
                cost 2^N-sized, so it reaches blocks the eigensolver cannot. Its first-order
                error is O(price^3 / gap^2) for well separated Bohr frequencies (F122's
                mechanism; experiments/THE_SPREAD_IS_A_RESONANCE.md for the imaginary second
                order); it reports `secular`, its slowest rate over the smallest gap between
                Bohr-frequency clusters, as a reading, not as a gate. It refuses a block whose
                largest cluster would need more than max_bytes (a degenerate H: the star's
                equivalent spokes put 143 699 of the 189 225 (2,2) coherences of N = 30 in one).

Per-site weights are basis-free: the trace of each site's disagreement projector over an
orthonormal basis of the slowest eigenspace, divided by its dimension (the mirror's
EigenspacePortfolio, same definition). Blocks (q, p) are the conjugates of (p, q) and carry
the same rates, so only p <= q is read.
"""
from dataclasses import dataclass, field

import numpy as np

_ZERO = 1e-9          # rates at or below this are the stationary kernel, not a decay
_CLUSTER = 1e-6       # relative eigenvalue distance read as one eigenspace (as the C# ledger)


@dataclass
class Reading:
    instrument: str
    block: tuple
    rate: float
    omega: float
    weights: np.ndarray        # per-site <Delta_l>, eigenspace average
    dimension: int             # eigenspace dimension at the slowest rate
    note: str = ""


def _eigenspace_weights(system, zi, zj, vecs):
    """Basis-free site weights of the span of `vecs` (columns over the (I, J) coherences)."""
    Q, R = np.linalg.qr(vecs / np.linalg.norm(vecs, axis=0))
    U, s, _ = np.linalg.svd(R)
    keep = s > 1e-6 * s.max()
    B = Q @ U[:, keep]
    P = (np.abs(B) ** 2).sum(axis=1).reshape(len(zi), len(zj))
    w = np.array([(P * (zi[:, l][:, None] != zj[:, l][None, :])).sum() for l in range(system.N)])
    return w / keep.sum(), int(keep.sum())


def _slowest_of(system, zi, zj, lam, vecs, name, block, note=""):
    rates = -lam.real
    alive = rates > _ZERO
    if not alive.any():
        return None
    lam0 = lam[alive][np.argmin(rates[alive])]
    # the most rotating member of the slowest shelf, as the mirror's Rotation voice
    shelf = alive & (np.abs(rates - (-lam0.real)) <= _CLUSTER * max(1.0, abs(lam0)))
    lam0 = lam[shelf][np.argmax(np.abs(lam[shelf].imag))]
    same = np.abs(lam - lam0) <= _CLUSTER * max(1.0, abs(lam0))
    w, dim = _eigenspace_weights(system, zi, zj, vecs[:, same])
    return Reading(name, block, float(-lam0.real), float(abs(lam0.imag)), w, dim, note)


class ExactBlocks:
    name = "exact"

    def __init__(self, max_rows=6000):
        self.max_rows = max_rows

    def covers(self, system, block, rows):
        if not system.is_number_conserving():
            return False, "H not number-conserving"
        if rows > self.max_rows:
            return False, f"{rows} rows > {self.max_rows}"
        return True, ""

    def slowest(self, system, block):
        zi, zj = system.sector(block[0])[0], system.sector(block[1])[0]
        lam, vecs = np.linalg.eig(system.block_liouvillian(*block))
        return _slowest_of(system, zi, zj, lam, vecs, self.name, block)


class SparseBlocks:
    """The exact block as a sparse matrix, its slowest modes by Arnoldi on the propagator e^{L t}.

    |e^{lambda t}| = e^{-rate t} whatever the Bohr frequency, so the largest-modulus eigenvalues of
    the propagator are exactly the slowest rates. Neither shift-invert (nearest the shift in
    modulus) nor ARPACK's which='LR' on L itself (Krylov favours the hull of the spectrum; with
    Bohr frequencies far larger than the rates the slowest mode sits inside it) targets that: 'LR'
    missed the N = 7 chain's (2, 5) block by 15 %. The propagator is applied by expm_multiply,
    never formed; the Krylov vectors are then Rayleigh-Ritz projected on L itself, which undoes
    the alias omega t mod 2 pi and returns lambda. The largest Ritz residual ||L v - lambda v|| of
    the slowest shelf is reported in the note. The kernel (rate 0) is found too, so k must exceed
    the kernel plus the slowest shelf; a shelf that fills every alive Ritz value is flagged."""
    name = "sparse"

    def __init__(self, max_rows=400_000, k=16, tol=1e-10, t=None):
        self.max_rows = max_rows
        self.k = k
        self.tol = tol
        self.t = t          # default 1 / mean site rate

    def covers(self, system, block, rows):
        if not system.is_number_conserving():
            return False, "H not number-conserving"
        if rows > self.max_rows:
            return False, f"{rows} rows > {self.max_rows}"
        if rows < self.k + 2:
            return False, f"{rows} rows: below ARPACK's reach, use exact"
        return True, ""

    def slowest(self, system, block):
        import scipy.sparse.linalg as spla
        L = system.block_liouvillian_sparse(*block)
        n, k = L.shape[0], self.k
        t = self.t if self.t is not None else 1.0 / np.mean(system.gamma)
        Lt = (L * t).tocsr()
        trace = Lt.diagonal().sum()
        E = spla.LinearOperator((n, n), dtype=complex,
                                matvec=lambda v: spla.expm_multiply(Lt, v, traceA=trace))
        ncv = min(n - 1, max(2 * k + 1, 40))
        _, krylov = spla.eigs(E, k=k, which="LM", ncv=ncv, tol=self.tol, maxiter=100_000)
        Q, _ = np.linalg.qr(krylov)
        LQ = L @ Q
        lam, y = np.linalg.eig(Q.conj().T @ LQ)
        vecs = Q @ y
        zi, zj = system.sector(block[0])[0], system.sector(block[1])[0]
        reading = _slowest_of(system, zi, zj, lam, vecs, self.name, block)
        if reading is not None:
            lam0 = -reading.rate
            shelf = np.abs(lam.real - lam0) <= _CLUSTER * max(1.0, abs(lam0))
            res = max(np.linalg.norm(LQ @ y[:, j] - lam[j] * vecs[:, j]) / np.linalg.norm(vecs[:, j])
                      for j in np.flatnonzero(shelf))
            alive = (-lam.real) > _ZERO
            reading.note = f"Ritz residual {res:.1e}"
            if shelf.sum() >= alive.sum():
                reading.note += f"; shelf fills all {alive.sum()} alive Ritz values: may be larger"
        return reading


class SecularRoute:
    name = "route"

    def __init__(self, window=0.0, tol=1e-9, max_sector=3000, max_bytes=2e9):
        self.window = window
        self.tol = tol
        self.max_sector = max_sector
        self.max_bytes = max_bytes

    def _groups(self, system, block):
        """The H eigendata of both sectors and the Bohr-frequency clusters (index groups)."""
        zi, Hp = system.sector(block[0])
        zj, Hq = system.sector(block[1])
        Ea, Ua = np.linalg.eigh(Hp)
        Eb, Ub = np.linalg.eigh(Hq)
        om = (Ea[:, None] - Eb[None, :]).ravel()
        order = np.argsort(om, kind='stable')
        cuts = np.flatnonzero(np.diff(om[order]) > max(self.window, self.tol)) + 1
        return zi, zj, Ua, Ub, om, np.split(order, cuts)

    def covers(self, system, block, rows):
        if not system.is_number_conserving():
            return False, "H not number-conserving"
        big = max(system.sector_size(block[0]), system.sector_size(block[1]))
        if big > self.max_sector:
            return False, f"sector {big} > {self.max_sector}"
        # a cluster of size c costs a rows x c matrix W and a c x c matrix M (complex)
        c = max(len(g) for g in self._groups(system, block)[-1])
        need = 16.0 * c * (rows + c)
        if need > self.max_bytes:
            return False, f"cluster {c} needs {need / 1e9:.0f} GB > {self.max_bytes / 1e9:.0f} GB"
        return True, ""

    def slowest(self, system, block):
        zi, zj, Ua, Ub, om, groups = self._groups(system, block)
        R = system.coherence_price(zi, zj)
        nb = Ub.shape[1]
        # a one-member cluster's rate is its coherence's price in the H eigenbasis; all pairs
        # (a, b) at once: P = |Ua|^2^T R |Ub|^2, one GEMM (no shared-index einsum: one core)
        P = ((np.abs(Ua) ** 2).T @ R @ (np.abs(Ub) ** 2)).ravel()
        single = np.array([g[0] for g in groups if len(g) == 1], dtype=int)
        alive = P[single] > _ZERO
        work = [g for g in groups if len(g) > 1]
        if alive.any():
            work.append(np.array([single[alive][np.argmin(P[single][alive])]]))
        best = None
        for g in work:
            X, Y = Ua[:, g // nb], Ub[:, g % nb]
            W = (X[:, None, :] * Y[None, :, :].conj()).reshape(-1, len(g))
            M = (W.conj() * R.reshape(-1, 1)).T @ W
            K = -1j * np.diag(om[g]) - M
            lam, x = np.linalg.eig(K)
            rates = -lam.real
            alive = rates > _ZERO
            if not alive.any():
                continue
            r = rates[alive].min()
            if best is None or r < best[0] - _CLUSTER * max(1.0, r):
                best = (r, g, lam, x, W)
        if best is None:
            return None
        _, g, lam, x, W = best
        centres = np.array([om[gg].mean() for gg in groups])
        gap = np.diff(np.sort(centres)).min() if len(centres) > 1 else np.inf
        reading = _slowest_of(system, zi, zj, lam, W @ x, self.name, block)
        reading.note = f"secular {reading.rate / gap:.2g}, largest cluster {max(len(gg) for gg in groups)}"
        return reading


@dataclass
class SlowRates:
    system: object
    rows: list = field(default_factory=list)     # (block, rows, {instrument: Reading or reason})

    def slowest(self):
        """The slowest reading over all blocks, preferring the exact instrument where it ran."""
        best = None
        for block, rows, by in self.rows:
            r = next((by[n] for n in ("exact", "sparse", "route") if isinstance(by.get(n), Reading)), None)
            if r is not None and (best is None or r.rate < best.rate):
                best = r
        return best

    def table(self, top=8):
        fold = "transpose + X^N" if self.system.has_xn_symmetry() else "transpose"
        lines = [f"# slow_rates: {self.system!r}  (one block per orbit: {fold})"]
        readings = []
        for block, rows, by in self.rows:
            rs = [v for v in by.values() if isinstance(v, Reading)]
            if rs:
                readings.append((min(r.rate for r in rs), block, rows, by))
        for _, block, rows, by in sorted(readings, key=lambda t: t[0])[:top]:
            parts = []
            for name, v in by.items():
                parts.append(f"{name} {v.rate:.5f}" if isinstance(v, Reading) else f"{name} -- ({v})")
            rd = [by[n] for n in ("exact", "sparse", "route") if isinstance(by.get(n), Reading)]
            for other in rd[1:]:
                parts.append(f"{other.instrument}-{rd[0].instrument} {(other.rate - rd[0].rate) / rd[0].rate:+.1e}")
            src = rd[0]
            w = " ".join(f"{100 * x:.0f}" for x in src.weights)
            lines.append(f"  {str(block):8s} {rows:6d} rows  " + "   ".join(parts)
                         + f"   |  sites % [{w}] dim {src.dimension}  {src.note}")
        return "\n".join(lines)


def slow_rates(system, instruments=None, progress=None):
    """Read one block of every orbit of `system` with every instrument that covers it.

    The orbits: transposition (p, q) -> (q, p) always (conjugate spectrum, same rates and site
    weights); and, when [H, X^N] = 0, the complement (p, q) -> (N - p, N - q), which carries the
    block over entry for entry at any a. The block read is the orbit's smallest (p, q).
    `progress`, if given, is called with one line per instrument and block (name, block, seconds)."""
    import time
    instruments = instruments or [ExactBlocks(), SparseBlocks(), SecularRoute()]
    out = SlowRates(system)
    N, xn = system.N, system.has_xn_symmetry()
    for (p, q), rows in system.blocks():
        orbit = {(p, q), (q, p)}
        if xn:
            orbit |= {(N - p, N - q), (N - q, N - p)}
        if (p, q) != min(orbit) or rows == 0:
            continue
        by = {}
        for ins in instruments:
            ok, why = ins.covers(system, (p, q), rows)
            t0 = time.time()
            by[ins.name] = (ins.slowest(system, (p, q)) or "no decaying mode") if ok else why
            if progress is not None and ok:
                progress(f"  {ins.name:6s} {str((p, q)):8s} {rows:7d} rows  {time.time() - t0:7.1f} s")
        out.rows.append(((p, q), rows, by))
    return out
