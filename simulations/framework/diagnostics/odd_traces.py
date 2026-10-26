"""The odd centred power sums of a Lindblad generator, modulo primes: the sextant's no, in Python.

The Python twin of compute/MirrorWorld/Sextant.cs OddTraceObstruction. If any similarity carries L to -L^† - c or to
-L - c, the spectrum is symmetric about the centre Tr(L)/dim, so every odd power of M = L - Tr(L)/dim has trace zero.
One odd trace that is not zero rules out every element of any kind (a no, exact); all odd traces up to the dimension
vanishing is evidence of a symmetric spectrum, not a proof (the integers could be multiples of the primes), and even a
symmetric spectrum is a pairing, not yet an element. The traces decide the complex spectrum's pairing, not whether the
rates (real parts) alone pair. The moment route the repo already holds: F117's trace face (Hardness.OddPowerSums),
F120's moment tower, F158's two-end count beside it.

Exact arithmetic: every entry of L lies in Q(i), read modulo a prime p = 1 (mod 4), where i is a number; a trace that
is not zero modulo p is not zero. The primes sit below 2^20 so that a matrix product in float64 stays exact up to
dim = 4^5 (sums below 2^50).

The setting is cube_tables.Setting, shared with the theodolite's twin (cube_letters) and cube_reflector: terms
(letters, sites, coefficient), one-site jumps (sites, 2x2 matrix, rate). Import it as
`from framework.diagnostics import odd_traces` (not exported from the package: it loads sympy).

Public API:
  PRIMES                                     the two primes, 1 (mod 4), below 2^20
  generator_mod_p(setting, p)                M = L - Tr(L)/dim over GF(p), row-major vectorization
  odd_traces_mod_p(setting, powers, p)       {k: Tr(M^k) mod p}
  odd_trace_obstruction(setting, max_power)  (k, p) for the first odd trace not zero modulo a prime, or None
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .. import cube_tables as ct

PRIMES = (1048573, 1048549)          # the two largest primes 1 (mod 4) below 2^20


def _sqrt_minus_one(p: int) -> int:
    for g in range(2, p):
        r = pow(g, (p - 1) // 4, p)
        if r * r % p == p - 1:
            return r
    raise ValueError(f"{p} has no square root of -1")


def _mod(x, p: int, i_p: int) -> int:
    """An element of Q(i) modulo p; raises when a denominator is a multiple of p."""
    x = sp.nsimplify(x)
    re, im = (sp.Rational(v) for v in x.as_real_imag())
    out = 0
    for part, unit in ((re, 1), (im, i_p)):
        if part.q % p == 0:
            raise ValueError(f"{x}: a denominator is a multiple of {p}")
        out += int(part.p) * pow(int(part.q), -1, p) * unit
    return out % p


def _letter(letter: str, p: int, i_p: int) -> np.ndarray:
    return np.array([[_mod(v, p, i_p) for v in row] for row in ct.PAULI[letter].tolist()], dtype=np.int64)


def _embed(n: int, ops: dict, p: int) -> np.ndarray:
    out = np.array([[1]], dtype=np.int64)
    for l in range(n):
        out = np.kron(out, ops.get(l, np.eye(2, dtype=np.int64))) % p
    return out


def generator_mod_p(setting: ct.Setting, p: int) -> np.ndarray:
    """M = L - Tr(L)/dim over GF(p), dim = 4^n, row-major vectorization (rho -> A rho B is A ⊗ B^T)."""
    n = setting.n
    if n > 5:
        raise ValueError("odd_traces reads up to 5 sites (dim 4^5)")
    i_p = _sqrt_minus_one(p)
    d = 2 ** n
    I = np.eye(d, dtype=np.int64)
    H = np.zeros((d, d), dtype=np.int64)
    for letters, sites, coef in ct.normalized_terms(setting.terms):
        ops = {s: _letter(c, p, i_p) for c, s in zip(letters, sites)}
        H = (H + _mod(coef, p, i_p) * _embed(n, ops, p)) % p
    L = ((p - i_p) * ((np.kron(H, I) - np.kron(I, H.T)) % p)) % p          # -i [H, rho]
    for sites, A, rate in setting.jumps:
        if len(sites) != 1:
            raise ValueError(f"jump on {sites}: one-site jumps only")
        J = _embed(n, {sites[0]: np.array([[_mod(v, p, i_p) for v in row] for row in A.tolist()], dtype=np.int64)}, p)
        Jbar = _embed(n, {sites[0]: np.array([[_mod(sp.conjugate(v), p, i_p) for v in row] for row in A.tolist()],
                                             dtype=np.int64)}, p)
        JdJ = (Jbar.T @ J) % p                                              # J^† J, J^† = conj(J)^T
        half = pow(2, -1, p)
        r = _mod(rate, p, i_p)
        L = (L + r * ((np.kron(J, Jbar) - half * ((np.kron(JdJ, I) + np.kron(I, JdJ.T)) % p)) % p)) % p
    dim = d * d
    shift = int(np.trace(L) % p) * pow(dim, -1, p) % p
    return (L - shift * np.eye(dim, dtype=np.int64)) % p


def _mul(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    return (a.astype(np.float64) @ b.astype(np.float64) % p).astype(np.int64)


def odd_traces_mod_p(setting: ct.Setting, powers, p: int) -> dict:
    """{k: Tr(M^k) mod p} for the given powers (any k >= 1)."""
    M = generator_mod_p(setting, p)
    want = sorted(set(powers))
    out, P = {}, M
    for k in range(1, want[-1] + 1):
        if k > 1:
            P = _mul(P, M, p)
        if k in want:
            out[k] = int(np.trace(P) % p)
    return out


def odd_trace_obstruction(setting: ct.Setting, max_power=None):
    """(k, p) for the first odd k with Tr(M^k) not zero modulo p, or None. By default up to the dimension when it is
    at most 256 (n <= 4) and up to 64 beyond, as the sextant reads; None below the dimension leaves the answer open."""
    dim = 4 ** setting.n
    upto = max_power if max_power is not None else (dim if dim <= 256 else 64)
    for p in PRIMES:
        traces = odd_traces_mod_p(setting, range(1, upto + 1, 2), p)
        for k in sorted(traces):
            if traces[k]:
                return k, p
    return None
