"""Sector-exact Lindblad propagator for the mediator-bridge Bell-on-vacuum runs.

The March runs (simulations/results/mediator_bridge_scale.txt, Tests 1 and 2)
start from |Phi+>_(0,1) (x) |0>^(N-2).  Both the Heisenberg bond term
J*(XX+YY+ZZ) and the Z-dephasing dissipator conserve popcount, so the whole
trajectory lives in the popcount-{0,2} sector: dimension 1 + C(N,2) instead of
2^N.  At N=11 that is 56 instead of 2048, i.e. rho is 56x56 instead of
2048x2048 -- roughly a 1300x saving on the propagation.

Conventions copied verbatim from the C# engine:
  H     = sum_bonds J * (X_a X_b + Y_a Y_b + Z_a Z_b)          Topology.cs
  drho  = -i[H,rho] + sum_k gamma_k (Z_k rho Z_k - rho)        LindbladPropagator.cs
  mask  = -2 * sum_{k: bit_k(i) != bit_k(j)} gamma_k
  qubit 0 is the MOST significant bit
  entropies in log base 2                                      DensityMatrixTools.cs
  RK4, no Hermiticity enforcement inside a step
"""

from __future__ import annotations

import itertools
import numpy as np


# ----------------------------------------------------------------------
# basis
# ----------------------------------------------------------------------

def sector_basis(n: int) -> tuple[list[int], dict[int, int]]:
    """Bitstrings of popcount 0 or 2 over n sites (qubit 0 = MSB)."""
    states = [0]
    for a, b in itertools.combinations(range(n), 2):
        states.append((1 << (n - 1 - a)) | (1 << (n - 1 - b)))
    states.sort()
    return states, {s: i for i, s in enumerate(states)}


def bit(state: int, q: int, n: int) -> int:
    return (state >> (n - 1 - q)) & 1


# ----------------------------------------------------------------------
# Hamiltonian and dephasing mask, restricted to the sector
# ----------------------------------------------------------------------

def build_h(n: int, bonds: list[tuple[int, int, float]],
            states: list[int], index: dict[int, int]) -> np.ndarray:
    """H restricted to the sector basis.

    XX+YY = 2(|01><10| + |10><01|) on the bond; ZZ is diagonal with
    z = +1 for bit 0 and z = -1 for bit 1.
    """
    m = len(states)
    h = np.zeros((m, m), dtype=complex)
    for col, s in enumerate(states):
        for (a, b, j) in bonds:
            ba, bb = bit(s, a, n), bit(s, b, n)
            za, zb = 1 - 2 * ba, 1 - 2 * bb
            h[col, col] += j * za * zb
            if ba != bb:
                flipped = s ^ (1 << (n - 1 - a)) ^ (1 << (n - 1 - b))
                h[index[flipped], col] += 2.0 * j
    return h


def build_mask(n: int, gammas, states: list[int]) -> np.ndarray:
    m = len(states)
    mask = np.zeros((m, m))
    for i, si in enumerate(states):
        for k, sk in enumerate(states):
            xor = si ^ sk
            rate = 0.0
            for q in range(n):
                if (xor >> (n - 1 - q)) & 1:
                    rate -= 2.0 * gammas[q]
            mask[i, k] = rate
    return mask


# ----------------------------------------------------------------------
# propagation
# ----------------------------------------------------------------------

class SectorPropagator:
    def __init__(self, n: int, bonds, gammas, states=None, index=None):
        if states is None:
            states, index = sector_basis(n)
        self.n = n
        self.states = states
        self.index = index
        self.h = build_h(n, bonds, states, index)
        self.mask = build_mask(n, gammas, states)

    def rhs(self, rho: np.ndarray) -> np.ndarray:
        comm = self.h @ rho - rho @ self.h
        return -1j * comm + self.mask * rho

    def propagate(self, rho0: np.ndarray, t_max: float, dt: float,
                  measure_times, callback):
        """RK4 exactly as the C# engine steps it: fixed dt, measurement fires
        when the running time first reaches a requested time."""
        rho = rho0.copy()
        times = sorted(measure_times)
        idx = 0
        if idx < len(times) and times[idx] <= dt / 2:
            callback(0.0, rho)
            idx += 1
        n_steps = int(t_max / dt) + 1
        for step in range(n_steps):
            k1 = self.rhs(rho)
            k2 = self.rhs(rho + 0.5 * dt * k1)
            k3 = self.rhs(rho + 0.5 * dt * k2)
            k4 = self.rhs(rho + dt * k3)
            rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            rho = 0.5 * (rho + rho.conj().T)          # EnforceHermiticity
            t = (step + 1) * dt
            while idx < len(times) and abs(t - times[idx]) < dt / 2:
                callback(t, rho)
                idx += 1
        return rho


    def propagate_every_step(self, rho0, t_max, dt, callback):
        """Same RK4, but the callback fires after every step (and at t=0).

        Used when the measurement grid must be as fine as the integrator, and
        when dt has to be shrunk for stability at large gamma so that a fixed
        list of measure times would no longer line up.
        """
        rho = rho0.copy()
        callback(0.0, rho)
        n_steps = int(t_max / dt) + 1
        for step in range(n_steps):
            k1 = self.rhs(rho)
            k2 = self.rhs(rho + 0.5 * dt * k1)
            k3 = self.rhs(rho + 0.5 * dt * k2)
            k4 = self.rhs(rho + dt * k3)
            rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            rho = 0.5 * (rho + rho.conj().T)
            callback((step + 1) * dt, rho)
        return rho

    def stable_dt(self, dt_wanted: float, safety: float = 1.5) -> float:
        """RK4 needs |lambda * dt| < ~2.78 for the pure decay part; the largest
        decay rate present is max|mask|.  Shrink dt when a large gamma_M would
        otherwise blow the integrator up."""
        worst = float(abs(self.mask).max())
        if worst == 0.0:
            return dt_wanted
        return min(dt_wanted, safety / worst)


def bell_on_vacuum(n: int, q0: int, q1: int, states, index) -> np.ndarray:
    psi = np.zeros(len(states), dtype=complex)
    psi[index[0]] = 1.0 / np.sqrt(2.0)
    psi[index[(1 << (n - 1 - q0)) | (1 << (n - 1 - q1))]] = 1.0 / np.sqrt(2.0)
    return np.outer(psi, psi.conj())


# ----------------------------------------------------------------------
# reduced states and mutual information
# ----------------------------------------------------------------------

def partial_trace(rho: np.ndarray, n: int, states: list[int],
                  keep: list[int]) -> np.ndarray:
    """Reduce to the sites in `keep`, on the support that is actually there."""
    keep = sorted(keep)
    k = len(keep)
    env = [q for q in range(n) if q not in keep]

    def split(s: int) -> tuple[int, int]:
        kb = 0
        for q in keep:
            kb = (kb << 1) | bit(s, q, n)
        eb = 0
        for q in env:
            eb = (eb << 1) | bit(s, q, n)
        return kb, eb

    labels = [split(s) for s in states]
    kept_vals = sorted({kb for kb, _ in labels})
    kmap = {v: i for i, v in enumerate(kept_vals)}
    red = np.zeros((len(kept_vals), len(kept_vals)), dtype=complex)
    for i, (kb_i, eb_i) in enumerate(labels):
        for j, (kb_j, eb_j) in enumerate(labels):
            if eb_i == eb_j:
                red[kmap[kb_i], kmap[kb_j]] += rho[i, j]
    return red


def entropy(rho: np.ndarray) -> float:
    ev = np.linalg.eigvalsh(rho)
    s = 0.0
    for p in ev.real:
        if p > 1e-15:
            s -= p * np.log2(p)
    return s


def mutual_information(rho: np.ndarray, n: int, states: list[int],
                       keep_a: list[int], keep_b: list[int]) -> float:
    ab = sorted(set(keep_a) | set(keep_b))
    rho_ab = partial_trace(rho, n, states, ab)
    # reduce further inside AB: rebuild an effective state list for AB
    ab_states = []
    for s in states:
        kb = 0
        for q in ab:
            kb = (kb << 1) | bit(s, q, n)
        ab_states.append(kb)
    kept_vals = sorted(set(ab_states))
    n_ab = len(ab)
    local_a = [ab.index(q) for q in sorted(keep_a)]
    local_b = [ab.index(q) for q in sorted(keep_b)]
    rho_a = partial_trace(rho_ab, n_ab, kept_vals, local_a)
    rho_b = partial_trace(rho_ab, n_ab, kept_vals, local_b)
    return entropy(rho_a) + entropy(rho_b) - entropy(rho_ab)


# ----------------------------------------------------------------------
# the mediator bridge topology (Topology.MediatorBridge, C#)
# ----------------------------------------------------------------------

def mediator_bridge(level: int, j_internal: float = 1.0,
                    j_bridge: float = 1.0, j_meta: float = 1.0):
    if level == 1:
        return [(0, 1, j_internal)]
    if level == 2:
        js = [j_internal, j_bridge, j_bridge, j_internal]
        return [(i, i + 1, js[i]) for i in range(4)]
    sub_size = bridge_size(level - 1)
    meta = sub_size
    bonds = list(mediator_bridge(level - 1, j_internal, j_bridge, j_meta))
    bonds.append((sub_size - 1, meta, j_meta))
    bonds.append((meta, meta + 1, j_meta))
    for (a, b, j) in mediator_bridge(level - 1, j_internal, j_bridge, j_meta):
        bonds.append((a + meta + 1, b + meta + 1, j))
    return bonds


def bridge_size(level: int) -> int:
    if level <= 0:
        return 1
    return 3 * (1 << (level - 1)) - 1
