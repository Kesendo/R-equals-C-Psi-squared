"""Letters for the reality-window flight: a density-matrix replay of the frozen circuits with
physical error letters, the readout forward model, the odd-Y reading vector, and the one-letter fits.

A letter is a physical error with a sign and a size, keyed by physical qubit or pair, so that arm R
(the same gates with physical 4 and 6 exchanged) and the fresh D compile see it where it physically sits.
"""
from __future__ import annotations

import math
import numpy as np

PATH = [4, 5, 6]
PAIRS = [(4, 5), (5, 6)]
LETTERS = ([("det", q) for q in PATH]          # Z detuning, rad/us, during every gate and idle gap
           + [("sxa", q) for q in PATH]        # sx/x amplitude error (fraction of the rotation angle)
           + [("sxp", q) for q in PATH]        # Z phase right after each sx/x pulse, rad
           + [("czp", p, q) for p in PAIRS for q in p]   # local Z phase on q after each cz of pair p, rad
           + [("czk", p) for p in PAIRS]       # cz conditional-phase error rzz(kappa), rad
           + [("zz", p) for p in PAIRS]        # static ZZ during all time, rad/us (rzz angle per us)
           + [("stk", (4, 5), 6), ("stk", (5, 6), 4)])   # Z phase on the spectator qubit at each cz of a pair, rad

I2 = np.eye(2, dtype=complex)
ZM = np.diag([1.0, -1.0]).astype(complex)


def on(q, M):
    ops = [I2, I2, I2]; ops[q] = M
    return np.kron(np.kron(ops[0], ops[1]), ops[2])


ZOP = [on(q, ZM) for q in range(3)]


def rz(t):
    return np.diag([np.exp(-0.5j * t), np.exp(0.5j * t)])


def rx(t):
    return np.array([[math.cos(t / 2), -1j * math.sin(t / 2)], [-1j * math.sin(t / 2), math.cos(t / 2)]])


def czm(a, b):
    d = np.ones(8, dtype=complex)
    for s in range(8):
        if (s >> (2 - a)) & 1 and (s >> (2 - b)) & 1:
            d[s] = -1
    return np.diag(d)


def rzz(a, b, t):
    return np.diag(np.exp(-0.5j * t * np.real(np.diag(ZOP[a] @ ZOP[b]))))


def decay(rho, q, tau, T1, T2):
    if tau <= 0:
        return rho
    g1 = 1 - math.exp(-tau / T1[q])
    K0 = on(q, np.array([[1, 0], [0, math.sqrt(1 - g1)]], dtype=complex))
    K1 = on(q, np.array([[0, math.sqrt(g1)], [0, 0]], dtype=complex))
    rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
    p = (1 - math.exp(-tau * max(1 / T2[q] - 0.5 / T1[q], 0.0))) / 2
    return (1 - p) * rho + p * ZOP[q] @ rho @ ZOP[q]


XP = np.array([[0, 1], [1, 0]], dtype=complex)
YP = np.array([[0, -1j], [1j, 0]])


def replay(ops, phys, T1, T2, L=None, max_step=0.05, pulse_step=0.012):
    """ops (name, logical qubits, params, duration_us) on logical 0, 1, 2, carrying the circuit's explicit
    schedule (a delay is idle time); phys[l] = physical qubit of logical l; T1/T2 indexed logically (us).
    L: {letter: value}; ideal gates when L is empty.

    Time-ordered, with finite pulses. The schedule fixes every operation's start. An sx or x is a pulse:
    H = (theta/duration)/2 (cos t X + sin t Y) on its qubit for its 36 ns, t the axis tilt (zero inside the
    library), with the continuous terms (detuning, static ZZ) acting at the same time, so a large detuning or
    ZZ turns the pulse as it would on the device. rz is a frame change at its instant; a cz acts at its start:
    the coherent terms acting on its two qubits during it are diagonal and commute with it, while decay does
    not (moving every cz from its start to its end changes an odd-Y reading by less than 1e-3). Amplitude damping and dephasing act throughout, interleaved in steps (max_step us when ZZ is on,
    pulse_step us inside pulses)."""
    from scipy.linalg import expm
    L = L or {}
    lg = {p: l for l, p in enumerate(phys)}
    det = [L.get(("det", phys[l]), 0.0) for l in range(3)]
    sxa = [L.get(("sxa", phys[l]), 0.0) for l in range(3)]
    sxp = [L.get(("sxp", phys[l]), 0.0) for l in range(3)]
    tilt = [L.get(("sxt", phys[l]), 0.0) for l in range(3)]      # outside the library: the pulse axis turned
    zz = {tuple(sorted((lg[a], lg[b]))): L.get(("zz", (a, b)), 0.0) for (a, b) in PAIRS}
    any_zz = any(zz.values())
    # the schedule: instantaneous events (time, order, kind, data) and pulses (start, end, qubit, rate)
    clock = [0.0, 0.0, 0.0]
    inst, pulses = [], []
    for i, (name, qs, prm, dur) in enumerate(ops):
        start = max(clock[q] for q in qs)
        for q in qs:
            clock[q] = start + dur
        if name == "delay":
            continue
        if name == "rz":
            inst.append((start, i, "rz", (qs[0], prm[0])))
        elif name == "cz":
            inst.append((start, i, "cz", tuple(sorted(qs))))
        elif name in ("sx", "x"):
            q = qs[0]
            theta = (math.pi / 2 if name == "sx" else math.pi) * (1 + sxa[q])
            pulses.append((start, start + dur, q, theta / dur))
            if sxp[q]:
                inst.append((start + dur, i + 0.5, "rz", (q, sxp[q])))
        else:
            raise ValueError(name)
    t_end = max(clock)
    inst.sort(key=lambda e: (e[0], e[1]))
    bounds = sorted({0.0, t_end} | {e[0] for e in inst} | {p[0] for p in pulses} | {p[1] for p in pulses})
    diag_H = np.zeros(8)
    for q in range(3):
        if det[q]:
            diag_H = diag_H + det[q] / 2 * np.real(np.diag(ZOP[q]))
    for (a, b), z in zz.items():
        if z:
            diag_H = diag_H + z / 2 * np.real(np.diag(ZOP[a] @ ZOP[b]))
    rho = np.zeros((8, 8), dtype=complex); rho[0, 0] = 1
    U_cache = {}

    def unitary(U):
        nonlocal rho
        rho = U @ rho @ U.conj().T

    def decay_all(h):
        nonlocal rho
        for q in range(3):
            rho = decay(rho, q, h, T1, T2)

    def slice_(t0, t1):
        nonlocal rho
        tau = t1 - t0
        if tau <= 0:
            return
        active = tuple(sorted((q, rate) for (s0, s1, q, rate) in pulses if s0 <= t0 + 1e-12 and t1 <= s1 + 1e-12))
        if not active:
            n = max(1, math.ceil(tau / max_step)) if any_zz else 1
            h = tau / n
            ph = np.exp(-1j * diag_H * h)
            for _ in range(n):
                rho = (ph[:, None] * rho) * ph.conj()[None, :]
                decay_all(h)
            return
        n = max(1, math.ceil(tau / pulse_step))
        h = tau / n
        key = (active, round(h, 12))
        if key not in U_cache:
            H = np.diag(diag_H).astype(complex)
            for q, rate in active:
                H = H + rate / 2 * on(q, math.cos(tilt[q]) * XP + math.sin(tilt[q]) * YP)
            U_cache[key] = expm(-1j * H * h)
        U = U_cache[key]
        for _ in range(n):
            unitary(U)
            decay_all(h)

    k = 0
    for j, t in enumerate(bounds):
        while k < len(inst) and inst[k][0] <= t + 1e-12:
            _, _, kind, data = inst[k]
            if kind == "rz":
                unitary(on(data[0], rz(data[1])))
            else:
                a, b = data
                pa, pb = sorted((phys[a], phys[b]))
                U = czm(a, b)
                kk = L.get(("czk", (pa, pb)), 0.0)
                if kk:
                    U = U @ rzz(a, b, kk)
                for l in (a, b):
                    ph_ = L.get(("czp", (pa, pb), phys[l]), 0.0)
                    if ph_:
                        U = U @ on(l, rz(ph_))
                spect = [l for l in range(3) if l not in (a, b)][0]
                st = L.get(("stk", (pa, pb), phys[spect]), 0.0)
                if st:
                    U = U @ on(spect, rz(st))
                unitary(U)
            k += 1
        if j + 1 < len(bounds):
            slice_(t, bounds[j + 1])
    return rho


# ── parities on the measured bits; bit masks (s0, s2, s1) over (q0, q2, q1) ──

MASKS = [(s0, s2, s1) for s0 in (0, 1) for s2 in (0, 1) for s1 in (0, 1) if (s0, s2, s1) != (0, 0, 0)]


def parities_rho(rho):
    out = {}
    for m in MASKS:
        O = np.eye(8)
        for s, Z in zip(m, (ZOP[0], ZOP[2], ZOP[1])):
            if s:
                O = O @ Z
        out[m] = float(np.real(np.trace(rho @ O)))
    return out


def probs_rho(rho):
    """Outcome probabilities indexed k = 4 q0 + 2 q2 + q1 (rho's basis index s has q0 = s>>2, q1 = s>>1, q2 = s)."""
    d = np.real(np.diag(rho))
    p = np.zeros(8)
    for s_ in range(8):
        q0, q1, q2 = (s_ >> 2) & 1, (s_ >> 1) & 1, s_ & 1
        p[4 * q0 + 2 * q2 + q1] += d[s_]
    return p


def parities_probs(p):
    out = {}
    for m in MASKS:
        sg = np.array([(-1) ** (m[0] * (k >> 2 & 1) + m[1] * (k >> 1 & 1) + m[2] * (k & 1)) for k in range(8)])
        out[m] = float(sg @ p)
    return out


def assignment_product(ro):
    """8x8 readout matrix A[measured, prepared] for independent bits, ro = [(p01, p10)] in order (q0, q2, q1)."""
    A = np.array([[1.0]])
    for p01, p10 in ro:
        A = np.kron(A, np.array([[1 - p01, p10], [p01, 1 - p10]]))
    return A


def counts_to_probs(counts):
    n = sum(counts.values())
    p = np.zeros(8)
    for k, c in counts.items():
        q1, q2, q0 = int(k[0]), int(k[1]), int(k[2])
        p[4 * q0 + 2 * q2 + q1] += c / n
    return p


def readout(par, ro):
    """Forward model of independent asymmetric readout. ro = [(p01, p10)] for (q0, q2, q1) as measured:
    E[s'|s] = a s + c with a = 1 - p01 - p10, c = p10 - p01."""
    a = [1 - p01 - p10 for p01, p10 in ro]
    c = [p10 - p01 for p01, p10 in ro]
    full = dict(par); full[(0, 0, 0)] = 1.0
    out = {}
    for m in MASKS:
        v = 0.0
        idx = [i for i in range(3) if m[i]]
        for sub in range(1 << len(idx)):
            keep = [idx[j] for j in range(len(idx)) if sub >> j & 1]
            coef = 1.0
            for i in idx:
                coef *= a[i] if i in keep else c[i]
            v += coef * full[tuple(1 if i in keep else 0 for i in range(3))]
        out[m] = v
    return out


def parities_counts(counts):
    """counts keys 'c2c1c0' = q1 q2 q0 -> parities over all 7 masks, and their multinomial covariance."""
    n = sum(counts.values())
    keys = list(counts)
    w = np.array([counts[k] for k in keys], dtype=float) / n
    bits = np.array([[int(k[2]), int(k[1]), int(k[0])] for k in keys])        # q0, q2, q1
    sign = {m: (-1.0) ** (bits @ np.array(m)) for m in MASKS}
    par = {m: float(sign[m] @ w) for m in MASKS}
    cov = {}
    for m1 in MASKS:
        for m2 in MASKS:
            prod = sign[m1] * sign[m2]
            cov[(m1, m2)] = (float(prod @ w) - par[m1] * par[m2]) / n
    return par, cov


def odd_y(basis, m):
    """Is parity m in basis (b0, b2) [q1 always read in Y] an odd-Y Pauli reading?"""
    b0, b2 = basis
    return ((m[0] and b0 == "Y") + (m[1] and b2 == "Y") + m[2]) % 2 == 1


def reading_index(names):
    """The ordered list of (circuit name, mask) for every odd-Y parity of every circuit."""
    out = []
    for nm in names:
        b = nm.split("_")[1]
        for m in MASKS:
            if odd_y((b[0], b[1]), m):
                out.append((nm, m))
    return out


def one_letter_fits(r, W, F):
    """r: residual vector (measured - baseline), W: inverse covariance, F: {letter: fingerprint}.
    Returns chi2_null and {letter: (amplitude, sigma_amplitude, chi2)}."""
    chi0 = float(r @ W @ r)
    out = {}
    for k, f in F.items():
        fWf = float(f @ W @ f)
        if fWf <= 0:
            out[k] = (0.0, float("inf"), chi0); continue
        a = float(f @ W @ r) / fWf
        out[k] = (a, 1 / math.sqrt(fWf), chi0 - a * a * fWf)
    return chi0, out
