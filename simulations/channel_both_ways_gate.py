"""Design gate for experiments/THE_CHANNEL_BOTH_WAYS_ON_A_PATH.md (a design record; not flown).

The object: F119 §7, readings 3 and 4 (docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md), on hardware. A payload
Phi+ on a pair, the rest in |0>, sent along a number-conserving real chain; the far pair's span entries p00,
p11 and c = <00|R|11>. Under number conservation p11 and c are the same in both readings (two excitations on
the pair force the rest empty), and the corollary makes them the same from both ends when the way back is the
transposed channel (TRANSFER). The readings differ only in p00: read with the silence it is the vacuum; read
by habit it also holds a return probability, which differs between the ends with a ZZ term and a site between
the pairs, and agrees on the XX chain through the free fermions (reading 3: READING).

Geometry: N = 5, pairs (0,1) and (3,4), site 2 between them. Bond gate exp(-i a h_b), h_b = (XX+YY)/2 +
DELTA*ZZ/2, a = J_b*dt*f, built as ONE two-qubit unitary (a UnitaryGate, synthesized by the transpiler: 3 CZ at
DELTA = 1.5, 2 CZ at DELTA = 0). Layers E = bonds (0,1),(2,3), O = bonds (1,2),(3,4). Arms (out, back, DELTA):
  main        the palindrome E/2 (O E)^(m-1) O E/2 both ways, DELTA = DELTA_ON
  main0       the same at DELTA = 0: free fermions, so the habit-read p00 is predicted to agree too
  transposed  (E O)^m out, (O E)^m back, DELTA_ON: the way back is the transposed circuit
  control     (E O)^m both ways, DELTA_ON: on this payload the first layer acts trivially, so the way back equals
              the transposed way back followed by one more O layer; a positive check of sensitivity, not a test
              specific to transposition (on this payload no depth-matched alternating way back exists but the
              transposed one)
  replica     a second copy of the main arm at REPLICA_STEP: an in-flight dispersion check
  end         zero depth: Phi+ prepared and read on the same pair, at either end: preparation-and-readout symmetry
A sufficient condition for the transfer: the read entries receive input only from the payload's (2,2) and
(0,2) entries, and on those the way back acts as the transpose of the way out.

Noise model NATIVE (every threshold): each bond gate is its synthesized native sequence (cz, sx, rz, x); idle
noise first, topping every qubit up to LAYER_US, as the as-late-as-possible schedule the runner uses places it;
then each pulse with idle noise on its own qubits for its duration (sx, x 32 ns; cz 68 ns) and its own
depolarizing (lambda_cz = 4r/3, lambda_sx = 2r; the RB error already contains the decay during the pulse, so
decay is counted twice, which widens the null). Idle: T1, T2* = T2_echo * T2STAR_FACTOR (unmeasured), gamma =
max(0, 1/T2* - 1/(2 T1))/2 (D-book), and optionally a static ZZ crosstalk on every bond. IDEALIZED placement
(rows of Part 1 and --limits, the stress row, --search): noise around whole bond gates.

Readout: the receiving pair in five settings (ZZ; XX, YY, XY, YX: Re c = (XX - YY)/4, Im c = -(XY + YX)/4),
the rest in Z. Confusion per qubit from CAL0/CAL1 pooled over the start and the end, inverted as a tensor
product, no clipping; the silence in the numerator (freqs_from_counts, here, is the estimator the runner applies to
hardware counts; the synthetic flights use a vectorized copy of the same steps). Synthetic
flights may inflate the science circuits' shot noise by a dispersion factor s (freq = p + s (sample - p)); the
thresholds are calibrated at s = DISPERSION_MARGIN.

Statistics (moduli; blind to diagonal phases that are static per shot and pass through the gates as a gauge;
the complex differences, reported for every arm, see them):
  transfer D   sqrt(dp11^2 + 2 d|c|^2) between the two directions (main, main0, transposed, control, replica)
  reading      |dp00| between the two directions by habit (main, main0), and with the silence (reported)
  spam D       sqrt(dp00^2 + dp11^2 + 2 d|c|^2) between the two zero-depth end readings (reported)
  sentinel     per direction P_silence - p11 - 1/2 (0 under number conservation; reported)
Decision:
  READING     SEEN if at every step the main arm's habit |dp00| exceeds its own null q995 AND the main0 arm's
              habit |dp00| stays within its own null q995 (reading 3: interaction separates, free fermions agree)
  CONTROL     SEEN if at every step the control's transfer D exceeds its own null q995
  TRANSFER    F = max over the main, main0 and transposed transfer D of D / its q995: HOLDS <= q99, DEPARTS >
              q999, INCONCLUSIVE between; WITHHELD if the replica check exceeds its q999
  HEADLINE    READING SEEN and CONTROL SEEN and TRANSFER HOLDS
Synthetic flights split in three: scales on the first third, F thresholds on the second, rates on the third.

Modes: --search, --certify, (default: Parts 1-3), --inputs FILE --candidate K --write-constants FILE, and
--limits: what the design can and cannot see. A number-conserving gate is its own transpose up to a relative
phase of its two sites, and on a path those phases glue into one global diagonal when each bond carries one gate,
so a static number-conserving error that stays on its bond leaves the transposed arm's moduli exactly equal; a
stochastic error, or one that reaches past its bond, needs the way back to be the exact reverse of the way out
(bonds within a layer reversed, each error before its gate), which the design's transposed arm is not; the free-chain contrast moves only at second order when the way back
stops being the transposed circuit; amplitude damping's jumps cannot feed p11 and c (F137); in the vacuum
population, and so in the silence-read fidelity, a uniform T1 is invisible by counting (the number-sector
populations close on the total excitation number when every site decays at one rate over the same intervals), a
non-uniform one moves it, the vacuum population also under unital noise that does not conserve the number;
L - L^T and L(I) are linear in F84's net-cooling vector (checked here).
"""
import argparse
import hashlib
import itertools
import json
import os

import numpy as np
from scipy.linalg import expm

N = 5
D = 1 << N
PAIR_A, PAIR_B = (0, 1), (3, 4)
J = (1.3, 2.2, 1.4, 1.4)
DELTA_ON = 1.5
JDT = 0.5
STEPS = (3, 5)
REPLICA_STEP = 5
SETTINGS = ("ZZ", "XX", "YY", "XY", "YX")
SHOTS = 16384
CAL_SHOTS_POOLED = 2 * 16384        # CAL0 and CAL1 each flown at the start and at the end
REPLICAS = 60000                    # split in three
SEED = 20261013
DISPERSION_MARGIN = 1.3
N_CZ_PREP = 1
T2STAR_FACTOR = 0.5
LAYER_US = 0.4                      # evolution-only schedule per layer at DELTA_ON, 0.396 us in the runner's certify run on the planning line (the runner checks it, +-20 %); the 2-CZ main0 layers are shorter, so their idle is overstated
NATIVE_US = {"sx": 0.032, "x": 0.032, "cz": 0.068, "rz": 0.0}
ZZ_PLANNING_KHZ = -3.9              # conditional ZZ shift measured on ibm_marrakesh (PRICE_PAIR run 4); a planning size
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results", "channel_both_ways")
DEFAULT_INPUTS = os.path.join(RESULTS, "kingston_lines_2026-10-13T0550Z.json")

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}
SWAP4 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)
ARMS = {"main": (2, 2, DELTA_ON), "main0": (2, 2, 0.0), "transposed": (1, 3, DELTA_ON), "control": (1, 1, DELTA_ON)}
AGREE = ("main", "main0", "transposed")


def op(P, q):
    return np.kron(np.kron(np.eye(1 << q), P), np.eye(1 << (N - 1 - q)))


def two(P2, q):
    return np.kron(np.kron(np.eye(1 << q), P2), np.eye(1 << (N - 2 - q)))


H_BOND = [(op(X, i) @ op(X, i + 1) + op(Y, i) @ op(Y, i + 1)) / 2 for i in range(N - 1)]
ZZ_BOND = [op(Z, i) @ op(Z, i + 1) / 2 for i in range(N - 1)]
PAULI2 = [np.kron(PAULI[a], PAULI[b]) for a in "IXYZ" for b in "IXYZ"]
PAULI1 = [PAULI[a] for a in "IXYZ"]


def gate_sha256():
    with open(os.path.abspath(__file__), "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def bond_unitary(i, f, delta):
    """The 4x4 bond unitary in chain order (site i first)."""
    a = f * JDT * J[i]
    return expm(-1j * a * ((np.kron(X, X) + np.kron(Y, Y)) / 2 + delta * np.kron(Z, Z) / 2))


# ------------------------------------------------------------------ the device
class Device:
    def __init__(self, T1, T2star, lam_cz, r_sx, lam_prep, readout, layer_us=LAYER_US):
        self.T1, self.T2star, self.lam_cz, self.r_sx = T1, T2star, lam_cz, r_sx
        self.lam_prep, self.readout, self.layer_us = dict(lam_prep), readout, layer_us
        self.t1 = self.dephasing = self.gates = self.prep = self.read = True
        self.native = True
        self.zz_khz = None            # static ZZ crosstalk per bond (conditional shift, kHz)
        self.coherent = None          # logical: f(bond) -> unitary on the chain after the bond gate
        self.native_coherent = None   # native: f(name, U) -> U'
        self.source = {}

    def without(self, **off):
        d = Device(self.T1, self.T2star, self.lam_cz, self.r_sx, self.lam_prep, self.readout, self.layer_us)
        d.source, d.coherent, d.native_coherent, d.native, d.zz_khz = (self.source, self.coherent, self.native_coherent,
                                                                     self.native, self.zz_khz)
        for k in ("t1", "dephasing", "gates", "prep", "read"):
            setattr(d, k, getattr(self, k) and not off.get(k, False))
        return d

    def idealized(self):
        d = self.without()
        d.native = False
        return d


def quiet_device():
    d = Device([1e12] * N, [1e12] * N, [0.0] * (N - 1), [0.0] * N, {PAIR_A: 0.0, PAIR_B: 0.0}, [(0, 0)] * N, 0.0)
    d.dephasing = d.t1 = d.gates = d.prep = False
    return d


def load_device(path=DEFAULT_INPUTS, candidate=0):
    with open(path, encoding="utf-8") as f:
        pull = json.load(f)
    c = pull["candidates"][candidate]
    q = c["qubits"]
    if any(v.get("sx_error") is None for v in q):
        raise ValueError("a qubit of the line has no sx_error in the properties pull; refused")
    lam_cz = [4.0 * r / 3.0 for r in c["cz_error"]]
    dev = Device([v["T1_us"] for v in q], [v["T2_us"] * T2STAR_FACTOR for v in q], lam_cz,
                 [v["sx_error"] for v in q],
                 {PAIR_A: 1 - (1 - lam_cz[0]) ** N_CZ_PREP, PAIR_B: 1 - (1 - lam_cz[3]) ** N_CZ_PREP},
                 [(v["p1_given_0"], v["p0_given_1"]) for v in q])
    dev.source = {"file": os.path.relpath(path, HERE).replace("\\", "/"), "candidate": candidate, "chain": c["chain"],
                  "properties_last_update": pull.get("properties_last_update")}
    return dev


def depolarize2(rho, i, lam):
    if lam <= 0:
        return rho
    acc = np.zeros_like(rho)
    for P in PAULI2:
        M = two(P, i)
        acc = acc + M @ rho @ M.conj().T
    return (1 - lam) * rho + lam * acc / 16


def depolarize1(rho, q, lam):
    if lam <= 0:
        return rho
    acc = np.zeros_like(rho)
    for P in PAULI1:
        M = op(P, q)
        acc = acc + M @ rho @ M.conj().T
    return (1 - lam) * rho + lam * acc / 4


def idle(rho, us, dev, qubits=None):
    if us <= 0:
        return rho
    for q in (range(N) if qubits is None else qubits):
        if dev.dephasing:
            g = max(0.0, 1 / dev.T2star[q] - 1 / (2 * dev.T1[q])) / 2
            p = (1 - np.exp(-2 * g * us)) / 2
            Zq = op(Z, q)
            rho = (1 - p) * rho + p * Zq @ rho @ Zq
        if dev.t1:
            gg = 1 - np.exp(-us / dev.T1[q])
            K0 = op(np.array([[1, 0], [0, np.sqrt(1 - gg)]], complex), q)
            K1 = op(np.array([[0, np.sqrt(gg)], [0, 0]], complex), q)
            rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
    return rho


def zz_crosstalk(rho, us, dev):
    """Static ZZ: H = sum_b 2 pi (zeta_b / 4) 1e-3 Z_b Z_b+1 (rad/us), zeta the conditional shift in kHz."""
    if dev.zz_khz is None or us <= 0:
        return rho
    Hzz = sum(2 * np.pi * (dev.zz_khz[b] / 4) * 1e-3 * op(Z, b) @ op(Z, b + 1) for b in range(N - 1))
    U = np.diag(np.exp(-1j * np.diag(Hzz) * us))
    return U @ rho @ U.conj().T


_NATIVE = {}


def native_ops(i, f, delta):
    """The bond gate's native sequence [(name, local qubits, unitary in chain order)], a basis transpile of the
    UnitaryGate; its product equals the bond unitary up to a global phase to 1e-12 (roundoff over ~15 pulses)."""
    key = (i, f, delta, J[i], JDT)
    if key not in _NATIVE:
        from qiskit import QuantumCircuit, transpile
        from qiskit.circuit.library import UnitaryGate
        from qiskit.quantum_info import Operator
        ref = bond_unitary(i, f, delta)
        qc = QuantumCircuit(2)
        qc.append(UnitaryGate(SWAP4 @ ref @ SWAP4), [0, 1])     # chain order (site i first) -> qiskit little-endian
        tq = transpile(qc, basis_gates=["cz", "sx", "rz", "x"], optimization_level=1, seed_transpiler=SEED)
        ops = []
        for inst in tq.data:
            qs = [tq.find_bit(q).index for q in inst.qubits]
            U = Operator(inst.operation).data
            if len(qs) == 2:
                U = SWAP4 @ U @ SWAP4 if qs == [0, 1] else U
                qs = [0, 1]
            ops.append((inst.operation.name, qs, U))
        W = np.eye(4, dtype=complex)
        for name, qs, U in ops:
            W = (np.kron(U, I2) if qs == [0] else np.kron(I2, U) if qs == [1] else U) @ W
        ph = np.vdot(W.reshape(-1), ref.reshape(-1))
        ph = ph / abs(ph)
        if np.abs(ph * W - ref).max() > 1e-12:
            raise RuntimeError(f"native decomposition of bond {i} deviates {np.abs(ph * W - ref).max():.1e}")
        _NATIVE[key] = ops
    return _NATIVE[key]


def native_counts(i, f, delta):
    ops = native_ops(i, f, delta)
    return sum(1 for n, _, _ in ops if n == "cz"), sum(1 for n, _, _ in ops if n in ("sx", "x"))


def layer(rho, bonds, f, delta, dev, coherent=None):
    coherent = coherent if coherent is not None else (dev.coherent if dev is not None else None)
    if dev is None:
        for i in bonds:
            U = two(bond_unitary(i, f, delta), i)
            rho = U @ rho @ U.conj().T
            if coherent is not None:
                C = coherent(i)
                rho = C @ rho @ C.conj().T
        return rho
    if dev.native:
        busy = [0.0] * N
        for i in bonds:
            for name, qs, _ in native_ops(i, f, delta):
                for q in qs:
                    busy[i + q] += NATIVE_US.get(name, 0.0)
        for q in range(N):                                   # as late as possible: the idle comes first
            rho = idle(rho, max(0.0, dev.layer_us - busy[q]), dev, [q])
        rho = zz_crosstalk(rho, dev.layer_us / 2, dev)
        for i in bonds:
            for name, qs, U in native_ops(i, f, delta):
                if dev.native_coherent is not None:
                    U = dev.native_coherent(name, U)
                M = op(U, i + qs[0]) if len(qs) == 1 else two(U, i)
                rho = M @ rho @ M.conj().T
                sites = [i + q for q in qs]
                if dev.gates:
                    if name == "cz":
                        rho = depolarize2(rho, i, dev.lam_cz[i])
                    elif name in ("sx", "x"):
                        rho = depolarize1(rho, sites[0], 2.0 * dev.r_sx[sites[0]])
                rho = idle(rho, NATIVE_US.get(name, 0.0), dev, sites)
            if coherent is not None:
                C = coherent(i)
                rho = C @ rho @ C.conj().T
        return zz_crosstalk(rho, dev.layer_us / 2, dev)
    rho = idle(rho, dev.layer_us / 2, dev)
    rho = zz_crosstalk(rho, dev.layer_us / 2, dev)
    for i in bonds:
        n_cz, n_sx = native_counts(i, f, delta)
        l2 = 1 - np.sqrt((1 - dev.lam_cz[i]) ** n_cz) if dev.gates else 0.0
        l1 = [1 - np.sqrt((1 - 2.0 * dev.r_sx[q]) ** (n_sx / 2)) if dev.gates else 0.0 for q in (i, i + 1)]
        rho = depolarize2(rho, i, l2)
        rho = depolarize1(depolarize1(rho, i, l1[0]), i + 1, l1[1])
        U = two(bond_unitary(i, f, delta), i)
        rho = U @ rho @ U.conj().T
        if coherent is not None:
            C = coherent(i)
            rho = C @ rho @ C.conj().T
        rho = depolarize1(depolarize1(rho, i, l1[0]), i + 1, l1[1])
        rho = depolarize2(rho, i, l2)
    rho = zz_crosstalk(rho, dev.layer_us / 2, dev)
    return idle(rho, dev.layer_us / 2, dev)


EVEN, ODD = (0, 2), (1, 3)


def evolve(rho, m, order, delta, dev=None, coherent=None):
    if order == 0:
        return rho
    if order in (1, 3):                              # 1: (E O)^m; 3: its transpose (O E)^m
        first, second = (EVEN, ODD) if order == 1 else (ODD, EVEN)
        for _ in range(m):
            rho = layer(rho, first, 1.0, delta, dev, coherent)
            rho = layer(rho, second, 1.0, delta, dev, coherent)
        return rho
    rho = layer(rho, EVEN, 0.5, delta, dev, coherent)     # 2: E/2 (O E)^(m-1) O E/2
    for k in range(m):
        rho = layer(rho, ODD, 1.0, delta, dev, coherent)
        rho = layer(rho, EVEN, 1.0 if k < m - 1 else 0.5, delta, dev, coherent)
    return rho


def index(bits):
    return int("".join(map(str, bits)), 2)


def bell(pair, dev=None):
    psi = np.zeros(D, complex)
    psi[0] = 1 / np.sqrt(2)
    b = [0] * N
    b[pair[0]] = b[pair[1]] = 1
    psi[index(b)] = 1 / np.sqrt(2)
    rho = np.outer(psi, psi.conj())
    if dev is not None and dev.prep:
        rho = depolarize2(rho, pair[0], dev.lam_prep[pair])
    return rho


def run_arm(arm, m, dev=None, coherent=None):
    o1, o2, delta = ARMS[arm]
    return (evolve(bell(PAIR_A, dev), m, o1, delta, dev, coherent), evolve(bell(PAIR_B, dev), m, o2, delta, dev, coherent))


# ------------------------------------------------------------------ the span entries, exact
def span_exact(rho, pair):
    """(silence, habit), each (p00, p11, c, P_silence) with c = <00|R|11>; silence unnormalized."""
    rest = [q for q in range(N) if q not in pair]

    def idx(pb, rb):
        b = [0] * N
        b[pair[0]], b[pair[1]] = pb
        for q, v in zip(rest, rb):
            b[q] = v
        return index(b)
    z = (0,) * len(rest)
    p_sil = sum(rho[idx(pb, z), idx(pb, z)].real for pb in itertools.product((0, 1), repeat=2))
    sil = (rho[idx((0, 0), z), idx((0, 0), z)].real, rho[idx((1, 1), z), idx((1, 1), z)].real,
           rho[idx((0, 0), z), idx((1, 1), z)], p_sil)
    hab = [0.0, 0.0, 0.0j, 1.0]
    for rb in itertools.product((0, 1), repeat=len(rest)):
        hab[0] += rho[idx((0, 0), rb), idx((0, 0), rb)].real
        hab[1] += rho[idx((1, 1), rb), idx((1, 1), rb)].real
        hab[2] += rho[idx((0, 0), rb), idx((1, 1), rb)]
    return np.array(sil, complex), np.array(hab, complex)


def _out(x):
    return float(x) if np.ndim(x) == 0 else x


def d_tr(a, b):
    """The transfer statistic sqrt(dp11^2 + 2 d|c|^2) on moduli (last axis = (p00, p11, c, P_silence))."""
    d = np.abs(a) - np.abs(b)
    return _out(np.sqrt(d[..., 1] ** 2 + 2 * d[..., 2] ** 2))


def d_full(a, b):
    d = np.abs(a) - np.abs(b)
    return _out(np.sqrt(d[..., 0] ** 2 + d[..., 1] ** 2 + 2 * d[..., 2] ** 2))


def d_p00(a, b):
    return _out(np.abs(np.real(a[..., 0]) - np.real(b[..., 0])))


def sentinel(s):
    """P_silence - p11 - 1/2 of one direction's silence entries."""
    return _out(np.real(s[..., 3]) - np.real(s[..., 1]) - 0.5)


# ------------------------------------------------------------------ measurement and the estimator
ROT = {"X": np.array([[1, 1], [1, -1]], complex) / np.sqrt(2),
       "Y": np.array([[1, -1j], [1, 1j]], complex) / np.sqrt(2),   # H . Sdg
       "Z": I2}


def conf(p1_given_0, p0_given_1):
    """Readout confusion, column = true bit, row = read bit."""
    return np.array([[1 - p1_given_0, p0_given_1], [p1_given_0, 1 - p0_given_1]])


def apply_confusion(v, mats):
    """v: (R, 32), chain order; mats: per qubit a (2,2) or (R,2,2) matrix. Returns (R, 32)."""
    R = v.shape[0]
    t = v.reshape((R,) + (2,) * N)
    for q in range(N):
        M = mats[q]
        if M.ndim == 2:
            t = np.moveaxis(np.tensordot(t, M, axes=([q + 1], [1])), -1, q + 1)
        else:
            t = np.moveaxis(np.einsum("r...k,rjk->r...j", np.moveaxis(t, q + 1, -1), M), -1, q + 1)
    return t.reshape(R, -1)


def true_probs(rho, pair, setting):
    """Pre-readout 32-outcome distribution, chain order (site 0 = most significant bit)."""
    U = np.eye(1, dtype=complex)
    for q in range(N):
        U = np.kron(U, ROT[setting[pair.index(q)]] if q in pair else I2)
    p = np.real(np.diag(U @ rho @ U.conj().T)).clip(0)
    return p / p.sum()


def confusion_from_cal(cal0_freq, cal1_freq):
    """Per-qubit confusion estimates from CAL0 (all |0>) and CAL1 (all |1>) frequencies, (R,32) each."""
    R = cal0_freq.shape[0]
    t0 = cal0_freq.reshape((R,) + (2,) * N)
    t1 = cal1_freq.reshape((R,) + (2,) * N)
    mats = []
    for q in range(N):
        ax = tuple(a + 1 for a in range(N) if a != q)
        p1g0 = t0.sum(axis=ax)[:, 1]
        p0g1 = t1.sum(axis=ax)[:, 0]
        M = np.zeros((R, 2, 2))
        M[:, 0, 0], M[:, 1, 0], M[:, 0, 1], M[:, 1, 1] = 1 - p1g0, p1g0, p0g1, 1 - p0g1
        mats.append(M)
    return mats


def invert(freq, mats):
    return apply_confusion(freq, [np.linalg.inv(M) for M in mats])


def _masks(pair):
    rest = [q for q in range(N) if q not in pair]
    bits = np.array([[(s >> (N - 1 - q)) & 1 for q in range(N)] for s in range(D)])
    return ~bits[:, rest].any(axis=1), bits[:, pair[0]], bits[:, pair[1]]


def span_from_freqs(freqs, pair):
    """freqs[setting] = inverted (R,32) frequencies over ALL shots. Returns (silence, habit), each (R,4)
    complex: (p00, p11, c, P_silence), c = <00|R|11>; P_silence averaged over the five settings."""
    silent, a, b = _masks(pair)
    sign = (-1.0) ** (a + b)
    out = []
    for sel in (silent, np.ones(D, bool)):
        z = freqs["ZZ"]
        p00 = z[:, sel & (a == 0) & (b == 0)].sum(axis=1)
        p11 = z[:, sel & (a == 1) & (b == 1)].sum(axis=1)
        E = {s: (freqs[s][:, sel] * sign[sel]).sum(axis=1) for s in ("XX", "YY", "XY", "YX")}
        c = (E["XX"] - E["YY"]) / 4 - 1j * (E["XY"] + E["YX"]) / 4
        p_sel = np.mean([freqs[s][:, sel].sum(axis=1) for s in SETTINGS], axis=0)
        out.append(np.stack([p00, p11, c, p_sel], axis=1).astype(complex))
    return out[0], out[1]


def counts_to_freq(counts):
    """qiskit counts (little-endian keys, clbit q = chain site q) -> chain-order frequency vector (32,)."""
    v = np.zeros(D)
    for key, n in counts.items():
        v[int(key.replace(" ", "")[::-1], 2)] += n
    return v / v.sum()


def freqs_from_counts(counts_by_key):
    """The committed estimator's first step on hardware counts: CAL0/CAL1 pooled over start and end, the per-qubit
    confusion inverted on every science key."""
    cal0 = (counts_to_freq(counts_by_key["CAL0|start"]) + counts_to_freq(counts_by_key["CAL0|end"])) / 2
    cal1 = (counts_to_freq(counts_by_key["CAL1|start"]) + counts_to_freq(counts_by_key["CAL1|end"])) / 2
    mats = confusion_from_cal(cal0[None, :], cal1[None, :])
    return {k: invert(counts_to_freq(counts_by_key[k])[None, :], mats) for k in circuit_keys()}


# ------------------------------------------------------------------ circuits: keys, probabilities, statistics
DATASETS = {"main": (2, DELTA_ON), "main0": (2, 0.0), "replica": (2, DELTA_ON), "first": (1, DELTA_ON),
            "transposed": (3, DELTA_ON), "end": (0, DELTA_ON)}


def circuit_keys():
    """Science keys in flight order: per step and setting the five-plus arms; the two zero-depth end arms first.
    First-order out is ONE dataset for the control and the transposed return."""
    keys = [f"end|A|0|{s}" for s in SETTINGS] + [f"end|B|0|{s}" for s in SETTINGS]
    for m in STEPS:
        for s in SETTINGS:
            keys += [f"main|fwd|{m}|{s}", f"main|bwd|{m}|{s}", f"main0|fwd|{m}|{s}", f"main0|bwd|{m}|{s}",
                     f"first|fwd|{m}|{s}", f"first|bwd|{m}|{s}", f"transposed|bwd|{m}|{s}"]
            if m == REPLICA_STEP:
                keys += [f"replica|fwd|{m}|{s}", f"replica|bwd|{m}|{s}"]
    return keys


def key_circuit(key):
    """(prepared pair, read pair, m, order, delta, setting) of a science key."""
    arm, direction, m, s = key.split("|")
    m = int(m)
    order, delta = DATASETS[arm]
    if arm == "end":
        prep = read = PAIR_A if direction == "A" else PAIR_B
    else:
        prep = PAIR_A if direction == "fwd" else PAIR_B
        read = PAIR_B if prep == PAIR_A else PAIR_A
    return prep, read, m, order, delta, s


def expected_cz(key):
    """The CZ count of a science circuit: preparation plus every bond gate's native count."""
    prep, read, m, order, delta, s = key_circuit(key)
    if order == 0:
        return N_CZ_PREP
    seq = []
    if order in (1, 3):
        first, second = (EVEN, ODD) if order == 1 else (ODD, EVEN)
        for _ in range(m):
            seq += [(b, 1.0) for b in first] + [(b, 1.0) for b in second]
    else:
        seq += [(b, 0.5) for b in EVEN]
        for k in range(m):
            seq += [(b, 1.0) for b in ODD] + [(b, 1.0 if k < m - 1 else 0.5) for b in EVEN]
    return N_CZ_PREP + sum(native_counts(b, f, delta)[0] for b, f in seq)


def model_probs(dev):
    """Pre-readout probabilities for every science key under the model."""
    cache, out = {}, {}
    for key in circuit_keys():
        prep, read, m, order, delta, s = key_circuit(key)
        ck = (prep, m, order, delta)
        if ck not in cache:
            cache[ck] = evolve(bell(prep, dev), m, order, delta, dev)
        out[key] = true_probs(cache[ck], read, s)
    return out


def spans_from_freqs(freqs_by_key):
    spans = {}
    for ds in ("end",):
        for direction in ("A", "B"):
            if f"end|{direction}|0|ZZ" in freqs_by_key:
                spans[("end", direction, 0)] = span_from_freqs({s: freqs_by_key[f"end|{direction}|0|{s}"] for s in SETTINGS},
                                                               PAIR_A if direction == "A" else PAIR_B)
    for m in STEPS:
        for ds, direction in (("main", "fwd"), ("main", "bwd"), ("main0", "fwd"), ("main0", "bwd"), ("first", "fwd"),
                              ("first", "bwd"), ("transposed", "bwd"), ("replica", "fwd"), ("replica", "bwd")):
            if f"{ds}|{direction}|{m}|ZZ" not in freqs_by_key:
                continue
            read = PAIR_B if direction == "fwd" else PAIR_A
            spans[(ds, direction, m)] = span_from_freqs({s: freqs_by_key[f"{ds}|{direction}|{m}|{s}"] for s in SETTINGS}, read)
    return spans


def statistics(spans, twins=None):
    st = {}
    if ("end", "A", 0) in spans:
        st[("spam", 0)] = d_full(spans[("end", "A", 0)][0], spans[("end", "B", 0)][0])
    for m in STEPS:
        mf, mb = spans[("main", "fwd", m)], spans[("main", "bwd", m)]
        zf, zb = spans[("main0", "fwd", m)], spans[("main0", "bwd", m)]
        ff, fb, tb = spans[("first", "fwd", m)], spans[("first", "bwd", m)], spans[("transposed", "bwd", m)]
        st[("main", m)] = d_tr(mf[0], mb[0])
        st[("main0", m)] = d_tr(zf[0], zb[0])
        st[("transposed", m)] = d_tr(ff[0], tb[0])
        st[("control", m)] = d_tr(ff[0], fb[0])
        st[("habit_main", m)] = d_p00(mf[1], mb[1])
        st[("habit_main0", m)] = d_p00(zf[1], zb[1])
        st[("silence_main", m)] = d_p00(mf[0], mb[0])
        st[("sentinel_fwd", m)] = sentinel(mf[0])
        st[("sentinel_bwd", m)] = sentinel(mb[0])
        if ("replica", "fwd", m) in spans:
            st[("replica", m)] = np.maximum(d_tr(mf[0], spans[("replica", "fwd", m)][0]),
                                            d_tr(mb[0], spans[("replica", "bwd", m)][0]))
        if twins is not None:
            st[("habit_main_null", m)] = d_p00(mf[1], twins[("main", "fwd", m)][1])
            st[("habit_main0_null", m)] = d_p00(zf[1], twins[("main0", "fwd", m)][1])
            st[("control_null", m)] = d_tr(ff[0], twins[("first", "fwd", m)][0])
    return st


def synthetic_flights(probs, readout, rng, R, s=1.0):
    """Counts-level flights from pre-readout probabilities: readout applied, CAL sampled and pooled, every circuit
    sampled once per flight (shared datasets shared), the science circuits' shot noise inflated by s, twin samples
    of main|fwd, main0|fwd and first|fwd for the nulls; the committed estimator applied."""
    cal_true = [conf(*readout[q]) for q in range(N)]

    def cal_freq(bit):
        p = np.zeros(D)
        p[0 if bit == 0 else D - 1] = 1
        p = apply_confusion(p[None, :], cal_true)[0]
        return rng.multinomial(CAL_SHOTS_POOLED, p / p.sum(), size=R) / CAL_SHOTS_POOLED
    mats = confusion_from_cal(cal_freq(0), cal_freq(1))
    inv = [np.linalg.inv(M) for M in mats]

    def sample(p):
        pr = apply_confusion(p[None, :], cal_true)[0]
        pr = pr / pr.sum()
        f = rng.multinomial(SHOTS, pr, size=R) / SHOTS
        return apply_confusion(pr[None, :] + s * (f - pr[None, :]), inv)
    freqs = {key: sample(p) for key, p in probs.items()}
    twin = {key: sample(p) for key, p in probs.items() if key.startswith(("main|fwd", "main0|fwd", "first|fwd"))}
    tw = {(ds, "fwd", m): span_from_freqs({st_: twin[f"{ds}|fwd|{m}|{st_}"] for st_ in SETTINGS}, PAIR_B)
          for m in STEPS for ds in ("main", "main0", "first")}
    return statistics(spans_from_freqs(freqs), tw)


def synthetic_flights_chunked(probs, readout, rng, R, s=1.0, chunk=10000):
    parts = [synthetic_flights(probs, readout, rng, min(chunk, R - k), s) for k in range(0, R, chunk)]
    return {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}


# ------------------------------------------------------------------ the decision
def scores(st, c):
    sc = c["scales"]
    F = np.max([st[(k, m)] / sc[f"{k}_{m}"] for k in AGREE for m in STEPS], axis=0)
    reading = np.all([st[("habit_main", m)] > sc[f"habit_main_null_{m}"] for m in STEPS] +
                     [st[("habit_main0", m)] <= sc[f"habit_main0_null_{m}"] for m in STEPS], axis=0)
    control = np.all([st[("control", m)] > sc[f"control_null_{m}"] for m in STEPS], axis=0)
    withheld = st[("replica", REPLICA_STEP)] > sc[f"replica_{REPLICA_STEP}"]
    return F, reading, control, withheld


def quantile_ci(x, q, z=1.96):
    """Empirical quantile with a distribution-free 95% interval from order statistics."""
    xs = np.sort(x)
    n = len(xs)
    half = z * np.sqrt(n * q * (1 - q))
    lo = int(max(0, np.floor(q * n - half)))
    hi = int(min(n - 1, np.ceil(q * n + half)))
    return float(np.quantile(xs, q)), float(xs[lo]), float(xs[hi])


def calibrate(st):
    R = len(next(iter(st.values())))
    a, b = R // 3, 2 * R // 3
    first = {k: v[:a] for k, v in st.items()}
    scales = {f"{k}_{m}": float(np.quantile(first[(k, m)], 0.995)) for k in AGREE for m in STEPS}
    for k in ("habit_main_null", "habit_main0_null", "control_null"):
        scales.update({f"{k}_{m}": float(np.quantile(first[(k, m)], 0.995)) for m in STEPS})
    scales[f"replica_{REPLICA_STEP}"] = float(np.quantile(first[("replica", REPLICA_STEP)], 0.999))
    F2 = scores({k: v[a:b] for k, v in st.items()}, {"scales": scales})[0]
    return scales, quantile_ci(F2, 0.99), quantile_ci(F2, 0.999), {k: v[b:] for k, v in st.items()}


def verdict(stats, constants):
    F, reading, control, withheld = scores({k: np.array([v]) for k, v in stats.items()}, constants)
    F, reading, control, withheld = float(F[0]), bool(reading[0]), bool(control[0]), bool(withheld[0])
    if withheld:
        transfer = "WITHHELD"
    else:
        transfer = "HOLDS" if F <= constants["F_hold"] else ("DEPARTS" if F > constants["F_depart"] else "INCONCLUSIVE")
    return {"F": F, "reading": "SEEN" if reading else "NOT SEEN", "control": "SEEN" if control else "NOT SEEN",
            "transfer": transfer, "replica_check": "EXCESS" if withheld else "OK",
            "headline": reading and control and transfer == "HOLDS"}


def rates(st, c):
    F, reading, control, withheld = scores(st, c)
    hold = (F <= c["F_hold"]) & ~withheld
    dep = (F > c["F_depart"]) & ~withheld
    return {"reading": np.mean(reading), "control": np.mean(control), "holds": np.mean(hold), "departs": np.mean(dep),
            "withheld": np.mean(withheld), "headline": np.mean(reading & control & hold)}


def side_readings(spans):
    """Reported, never gated: complex coherence differences of every arm, the zero-depth spam statistic, the
    renormalized silence-read CPsi_pur of the main arms, the sentinels and their difference."""
    def one(x):
        return np.ravel(x)[0]
    out = {}
    if ("end", "A", 0) in spans:
        out["spam_D"] = float(one(d_full(spans[("end", "A", 0)][0], spans[("end", "B", 0)][0])))
        out["spam_complex_c_diff"] = float(abs(one(spans[("end", "A", 0)][0][..., 2]) - one(spans[("end", "B", 0)][0][..., 2])))
    for m in STEPS:
        for name, (a, b) in (("main", (("main", "fwd", m), ("main", "bwd", m))), ("main0", (("main0", "fwd", m), ("main0", "bwd", m))),
                             ("transposed", (("first", "fwd", m), ("transposed", "bwd", m))),
                             ("control", (("first", "fwd", m), ("first", "bwd", m)))):
            out[f"complex_c_diff_{name}_{m}"] = float(abs(one(spans[a][0][..., 2]) - one(spans[b][0][..., 2])))
        for arm in ("main", "main0"):
            for name in ("fwd", "bwd"):
                p00, p11, c = (one(spans[(arm, name, m)][0][..., k]) for k in range(3))
                tr = (p00 + p11).real
                r = np.array([[p00.real, c], [np.conj(c), p11.real]]) / tr
                out[f"cpsi_pur_{arm}_{name}_{m}"] = float(np.trace(r @ r).real * 2 * abs(c / tr) / 3)
        out[f"sentinel_fwd_{m}"] = float(one(sentinel(spans[("main", "fwd", m)][0])))
        out[f"sentinel_bwd_{m}"] = float(one(sentinel(spans[("main", "bwd", m)][0])))
    return out


# ------------------------------------------------------------------ the printout
def exact_row(d, m, coherent=None):
    spans = {arm: [span_exact(x, p) for x, p in zip(run_arm(arm, m, d, coherent), (PAIR_B, PAIR_A))] for arm in ARMS}
    return {"main": d_tr(spans["main"][0][0], spans["main"][1][0]),
            "main0": d_tr(spans["main0"][0][0], spans["main0"][1][0]),
            "transposed": d_tr(spans["control"][0][0], spans["transposed"][1][0]),
            "control": d_tr(spans["control"][0][0], spans["control"][1][0]),
            "habit_main": d_p00(spans["main"][0][1], spans["main"][1][1]),
            "habit_main0": d_p00(spans["main0"][0][1], spans["main0"][1][1]),
            "silence_main": d_p00(spans["main"][0][0], spans["main"][1][0]),
            "c_main": abs(spans["main"][0][0][2] - spans["main"][1][0][2]),
            "c_transposed": abs(spans["control"][0][0][2] - spans["transposed"][1][0][2]),
            "sent": (sentinel(spans["main"][0][0]), sentinel(spans["main"][1][0])),
            "reach": (spans["main"][0][0][1].real, abs(spans["main"][0][0][2]), spans["main0"][0][0][1].real,
                      abs(spans["main0"][0][0][2]), spans["control"][0][0][1].real, abs(spans["control"][0][0][2]))}


def part1(dev):
    print("Part 1 (exact, no shots). Transfer D: main, main0, transposed, control | habit |dp00|: main, main0 | silence"
          " |dp00| main | sentinels fwd, bwd | reach = forward silence p11/|c| of main, main0, first order")
    print(f"   device: {dev.source}")
    for label, d in [("noiseless", None),
                     ("idealized placement, all noise", dev.idealized()),
                     ("native: dephasing only", dev.without(t1=True, gates=True, prep=True)),
                     ("native: T1 only", dev.without(dephasing=True, gates=True, prep=True)),
                     ("native: gates + preparation only", dev.without(t1=True, dephasing=True)),
                     ("native: all (readout enters in Part 2)", dev)]:
        print(f"   {label}")
        for m in STEPS:
            r = exact_row(d, m)
            print(f"      m={m}: {r['main']:.1e} {r['main0']:.1e} {r['transposed']:.1e} {r['control']:.4f} | {r['habit_main']:.4f}"
                  f" {r['habit_main0']:.1e} | {r['silence_main']:.1e} | {r['sent'][0]:+.4f}, {r['sent'][1]:+.4f} | reach"
                  f" {r['reach'][0]:.4f}/{r['reach'][1]:.3f}, {r['reach'][2]:.4f}/{r['reach'][3]:.3f}, {r['reach'][4]:.4f}/{r['reach'][5]:.3f}")
    print("   zero-depth end arms, native all: silence (p00, p11, |c|) at A and at B, and their D:")
    a = span_exact(bell(PAIR_A, dev), PAIR_A)[0]
    b = span_exact(bell(PAIR_B, dev), PAIR_B)[0]
    print(f"      A ({a[0].real:.4f}, {a[1].real:.4f}, {abs(a[2]):.4f})  B ({b[0].real:.4f}, {b[1].real:.4f}, {abs(b[2]):.4f})  D {d_full(a, b):.1e}")


def special_rows(dev):
    def lay_run(prep, m, after):
        rho = bell(prep)

        def lay(rho, bonds, f):
            for i in bonds:
                U = two(bond_unitary(i, f, DELTA_ON), i)
                rho = after(i, rho, U)
            return rho
        rho = lay(rho, EVEN, 0.5)
        for k in range(m):
            rho = lay(rho, ODD, 1.0)
            rho = lay(rho, EVEN, 1.0 if k < m - 1 else 0.5)
        return rho

    def rx_sym(i, rho, U):
        R = expm(-1j * 0.005 * op(X, i))
        W = R @ U @ R
        return W @ rho @ W.conj().T

    def deph_left(i, rho, U):
        rho = U @ rho @ U.conj().T
        Zq = op(Z, i)
        return 0.99 * rho + 0.01 * Zq @ rho @ Zq

    print("   counterexamples to a single property (noiseless otherwise; max over the steps of the main transfer D):")
    for label, after in (("Rx(0.01) on the left qubit before and after every gate (U^T = U kept, number not conserved)", rx_sym),
                         ("Z dephasing p = 0.01 on the left qubit after every gate (off the palindrome)", deph_left)):
        print(f"      {label}: {max(d_tr(span_exact(lay_run(PAIR_A, m, after), PAIR_B)[0], span_exact(lay_run(PAIR_B, m, after), PAIR_A)[0]) for m in STEPS):.1e}")
    dd = dev.idealized().without(t1=True, dephasing=True, prep=True)
    print(f"      two-qubit and one-qubit depolarizing around every gate (transpose kept, number not conserved):"
          f" {max(exact_row(dd, m)['main'] for m in STEPS):.1e}")
    xy = lambda i: expm(-1j * 0.05 * (two(np.kron(X, Y), i) - two(np.kron(Y, X), i)) / 2)
    print(f"      an XY - YX part of 0.05 rad after every gate (not its own transpose): {max(exact_row(None, m, xy)['main'] for m in STEPS):.1e}")
    print("   diagonal errors (noiseless otherwise; max over the steps): transfer D main / transposed / complex c difference main")
    rng = np.random.default_rng(11)
    zs = rng.normal(0, 0.05, N)
    zz = rng.normal(0, 0.05, N - 1)
    rows = (("a static Z profile after every gate (sd 0.05 rad)", lambda i: expm(-1j * sum(zs[q] * op(Z, q) for q in range(N)) / 2)),
            ("Rz(0.02) on the left qubit after every gate", lambda i: expm(-1j * 0.01 * op(Z, i))),
            ("a static ZZ profile after every gate (sd 0.05 rad)", lambda i: expm(-1j * sum(zz[b] * op(Z, b) @ op(Z, b + 1) for b in range(N - 1)) / 2)))
    for label, fn in rows:
        rr = [exact_row(None, m, fn) for m in STEPS]
        print(f"      {label:52s} {max(r['main'] for r in rr):.1e} / {max(r['transposed'] for r in rr):.1e} / {max(r['c_main'] for r in rr):.1e}")
    rng2 = np.random.default_rng(12)
    print(f"      a Z phase drawn anew at every gate occurrence (sd 0.05 rad): main transfer D "
          f"{max(exact_row(None, m, lambda i: expm(-1j * rng2.normal(0, 0.05) * op(Z, i) / 2))['main'] for m in STEPS):.1e}")
    d = dev.without()
    d.zz_khz = [ZZ_PLANNING_KHZ] * (N - 1)
    rr = [exact_row(d, m) for m in STEPS]
    r0 = [exact_row(dev, m) for m in STEPS]
    print(f"      static ZZ crosstalk {ZZ_PLANNING_KHZ} kHz on every bond, native device model: main transfer D"
          f" {max(r['main'] for r in rr):.1e} (without it {max(r['main'] for r in r0):.1e})")


def t1_stress(dev):
    """Idealized placement: amplitude damping's jumps do not feed p11 or c (strong uneven T1 and T2*)."""
    d = dev.idealized().without(gates=True, prep=True)
    d.T1 = [20.0 * f for f in (1.0, 1.7, 0.6, 1.3, 0.8)]
    d.T2star = [8.0 * f for f in (1.2, 0.7, 1.0, 1.4, 0.9)]
    cells = []
    for m in STEPS:
        o, b = run_arm("main", m, d)
        so, sb = span_exact(o, PAIR_B)[0], span_exact(b, PAIR_A)[0]
        cells.append(f"m={m}: {d_tr(so, sb):.1e} / {abs(so[2] - sb[2]):.1e} / P_silence {so[3].real:.4f}, {sb[3].real:.4f}"
                     f" / sentinel {sentinel(so):+.4f} / |c| {abs(so[2]):.3f}")
    print("   stress, idealized placement, T1 = 20 us x (1, 1.7, 0.6, 1.3, 0.8), T2* = 8 us x (1.2, 0.7, 1, 1.4, 0.9); main:"
          " transfer D / complex c difference / P_silence fwd, bwd / sentinel fwd / |c|:  " + " | ".join(cells))


def native_coherent():
    def cz_phase(eps):
        P = np.diag([1, 1, 1, np.exp(-1j * eps)])
        return lambda name, U: P @ U if name == "cz" else U

    def sx_over(delta):
        Rx = expm(-1j * (np.pi / 2) * (1 + delta) * X / 2)
        return lambda name, U: Rx if name == "sx" else U
    return [("CZ conditional phase 0.01 rad per CZ", cz_phase(0.01)), ("CZ conditional phase 0.02 rad per CZ", cz_phase(0.02)),
            ("SX over-rotation 0.5 %", sx_over(0.005)), ("SX over-rotation 1 %", sx_over(0.01))]


def native_table(dev):
    print("   native coherent errors (max over the steps): main / transposed transfer D, noiseless otherwise; then on top of"
          " the native device model")
    for label, fn in native_coherent():
        q = quiet_device()
        q.native_coherent = fn
        d = dev.without()
        d.native_coherent = fn
        rq = [exact_row(q, m) for m in STEPS]
        rd = [exact_row(d, m) for m in STEPS]
        print(f"      {label:38s} {max(r['main'] for r in rq):.1e} / {max(r['transposed'] for r in rq):.1e};  on the device"
              f" {max(r['main'] for r in rd):.1e} / {max(r['transposed'] for r in rd):.1e}")


def part2(dev, rng, write=None):
    print()
    print(f"Part 2 (counts level, native device model): {REPLICAS} synthetic flights split in three (scales / F thresholds /"
          f" rates), {SHOTS} shots per circuit, CAL pooled {CAL_SHOTS_POOLED} shots per type, the science circuits' shot"
          f" noise inflated by s = {DISPERSION_MARGIN} for calibration, seed {SEED}")
    probs = model_probs(dev)
    st = synthetic_flights_chunked(probs, dev.readout, rng, REPLICAS, DISPERSION_MARGIN)
    scales, f99, f999, held = calibrate(st)
    for m in STEPS:
        print(f"   m={m} (t={m * JDT:.2f}): transfer D median main {np.median(st[('main', m)]):.4f} (scale {scales[f'main_{m}']:.4f}),"
              f" main0 {np.median(st[('main0', m)]):.4f} (scale {scales[f'main0_{m}']:.4f}), transposed"
              f" {np.median(st[('transposed', m)]):.4f} (scale {scales[f'transposed_{m}']:.4f}); control"
              f" {np.median(st[('control', m)]):.4f} [q001 {np.quantile(st[('control', m)], .001):.4f}] (null {scales[f'control_null_{m}']:.4f})")
        print(f"        habit |dp00| median main {np.median(st[('habit_main', m)]):.4f} [q001 {np.quantile(st[('habit_main', m)], .001):.4f}]"
              f" (null {scales[f'habit_main_null_{m}']:.4f}); main0 {np.median(st[('habit_main0', m)]):.4f} [q999"
              f" {np.quantile(st[('habit_main0', m)], .999):.4f}] (null {scales[f'habit_main0_null_{m}']:.4f}); silence main"
              f" {np.median(st[('silence_main', m)]):.4f}")
    print(f"   zero-depth spam D median {np.median(st[('spam', 0)]):.4f}; replica threshold (q999 at m = {REPLICA_STEP})"
          f" {scales[f'replica_{REPLICA_STEP}']:.4f}")
    print(f"   F: q99 {f99[0]:.3f} [{f99[1]:.3f}, {f99[2]:.3f}], q999 {f999[0]:.3f} [{f999[1]:.3f}, {f999[2]:.3f}]"
          f" (q999 rests on about {int(0.001 * (REPLICAS // 3))} tail samples)")
    constants = {"version": "v5-draft", "gate_sha256": gate_sha256(), "device": dev.source, "J": list(J),
                 "DELTA_ON": DELTA_ON, "JDT": JDT, "steps": list(STEPS), "replica_step": REPLICA_STEP,
                 "layer_us": dev.layer_us, "t2star_factor": T2STAR_FACTOR, "dispersion_margin": DISPERSION_MARGIN,
                 "shots": SHOTS, "cal_shots_pooled": CAL_SHOTS_POOLED, "replicas": REPLICAS, "seed": SEED,
                 "scales": scales, "F_hold": f99[0], "F_hold_ci": list(f99[1:]), "F_depart": f999[0], "F_depart_ci": list(f999[1:])}
    r = rates(held, constants)
    print(f"   held-out third, s = {DISPERSION_MARGIN}: " + ", ".join(f"P({k}) {v:.4f}" for k, v in r.items()))
    if write:
        with open(write, "w", encoding="utf-8") as f:
            json.dump(constants, f, indent=1)
        print(f"   constants written: {os.path.relpath(write, HERE)}")
    return constants, probs


def part3(dev, rng, constants, probs, R=6000):
    print()
    print(f"Part 3 (power and false rates): {R} synthetic flights per alternative, native device model plus the alternative,"
          " judged with the constants above")
    alts = [("the device model, s = 1.0", dev, 1.0), ("the device model, s = 1.3", dev, 1.3),
            ("the device model, s = 1.6", dev, 1.6), ("the device model, s = 2.0", dev, 2.0)]
    d_echo = dev.without()
    d_echo.T2star = [t / T2STAR_FACTOR for t in dev.T2star]
    alts.append(("T2* = T2 (echo), s = 1.0", d_echo, 1.0))
    d_zz = dev.without()
    d_zz.zz_khz = [ZZ_PLANNING_KHZ] * (N - 1)
    alts.append((f"static ZZ crosstalk {ZZ_PLANNING_KHZ} kHz, s = 1.0", d_zz, 1.0))
    for label, fn in native_coherent():
        d = dev.without()
        d.native_coherent = fn
        alts.append((label + ", s = 1.0", d, 1.0))
    for dl in (0.02, 0.05):
        d = dev.without()
        d.lam_prep[PAIR_B] = d.lam_prep[PAIR_B] + dl
        alts.append((f"Bell preparation on (3,4), {dl} more depolarizing, s = 1.0", d, 1.0))
    for label, d, s in alts:
        pr = probs if d is dev else model_probs(d)
        r = rates(synthetic_flights(pr, d.readout, rng, R, s), constants)
        print(f"   {label:58s} " + "  ".join(f"{k} {v:.3f}" for k, v in r.items()))


# ------------------------------------------------------------------ circuits for the runner
def build_circuit(prep_pair, read_pair, m, order, delta, setting=None):
    """The logical circuit on N qubits (chain site q = qubit q = clbit q). setting=None: no measurement."""
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import UnitaryGate
    qc = QuantumCircuit(N, N if setting else 0)
    qc.h(prep_pair[0])
    qc.cx(prep_pair[0], prep_pair[1])
    qc.barrier()

    def lay(bonds, f):
        for i in bonds:
            qc.append(UnitaryGate(SWAP4 @ bond_unitary(i, f, delta) @ SWAP4, label=f"b{i}"), [i, i + 1])
        qc.barrier()
    if order in (1, 3):
        first, second = (EVEN, ODD) if order == 1 else (ODD, EVEN)
        for _ in range(m):
            lay(first, 1.0)
            lay(second, 1.0)
    elif order == 2:
        lay(EVEN, 0.5)
        for k in range(m):
            lay(ODD, 1.0)
            lay(EVEN, 1.0 if k < m - 1 else 0.5)
    if setting:
        for q, b in zip(read_pair, setting):
            if b == "X":
                qc.h(q)
            elif b == "Y":
                qc.sdg(q)
                qc.h(q)
        qc.measure(range(N), range(N))
    return qc


def key_to_circuit(key):
    prep, read, m, order, delta, s = key_circuit(key)
    return build_circuit(prep, read, m, order, delta, s)


def certify():
    """(1) logical qiskit circuits vs the model, COMPLEX density matrices; (2) probabilities over every science key;
    (3) the native sequences vs the logical model. Roundoff error model: ~1e-16 per operation over ~10^3 operations."""
    from qiskit.quantum_info import Statevector
    perm = [int(format(s, f"0{N}b")[::-1], 2) for s in range(D)]       # qiskit little-endian -> chain order
    worst = 0.0
    for prep in (PAIR_A, PAIR_B):
        for order in (1, 2, 3):
            for delta in (DELTA_ON, 0.0):
                for m in STEPS:
                    psi = Statevector(build_circuit(prep, None, m, order, delta)).data[perm]
                    worst = max(worst, np.abs(np.outer(psi, psi.conj()) - evolve(bell(prep), m, order, delta)).max())
    print(f"certify (1): max |rho_qiskit - rho_model| (complex), both preparations, three orders, both DELTA, all steps = {worst:.1e}")
    worst_p = 0.0
    for key in circuit_keys():
        qc = key_to_circuit(key)
        qc.remove_final_measurements()
        prep, read, m, order, delta, s = key_circuit(key)
        p_q = Statevector(qc).probabilities()[perm]
        worst_p = max(worst_p, np.abs(p_q - true_probs(evolve(bell(prep), m, order, delta), read, s)).max())
    print(f"certify (2): max |p_qiskit - p_model| over every science key = {worst_p:.1e}")
    nd = quiet_device()
    worst_n = max(np.abs(evolve(bell(p), m, o, dl, nd) - evolve(bell(p), m, o, dl)).max()
                  for p in (PAIR_A, PAIR_B) for o in (1, 2, 3) for dl in (DELTA_ON, 0.0) for m in STEPS)
    counts = sorted({(dl, native_counts(i, f, dl)) for i in range(N - 1) for f in (0.5, 1.0) for dl in (DELTA_ON, 0.0)})
    print(f"certify (3): the native sequences reproduce the model, max deviation {worst_n:.1e}; (DELTA, (cz, sx + x)) per bond gate {counts}")
    return worst, worst_p, worst_n


def search(dev, trials=160, seed=7):
    """Seeded random search: J in [0.5, 2.2]^4, DELTA, dt; per step the habit-read |dp00|, the control's transfer
    D and the reach (idealized placement, planning device)."""
    global J, JDT
    dev = dev.idealized()
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(trials):
        J = tuple(float(x) for x in np.round(rng.uniform(0.5, 2.2, 4), 2))
        delta = float(rng.choice([0.5, 0.75, 1.0, 1.5]))
        JDT = float(rng.choice([0.25, 0.35, 0.5]))
        for m in range(3, 10):
            o = evolve(bell(PAIR_A, dev), m, 2, delta, dev)
            b = evolve(bell(PAIR_B, dev), m, 2, delta, dev)
            co = evolve(bell(PAIR_A, dev), m, 1, delta, dev)
            cb = evolve(bell(PAIR_B, dev), m, 1, delta, dev)
            so, ho = span_exact(o, PAIR_B)
            sb, hb = span_exact(b, PAIR_A)
            sc, sd = span_exact(co, PAIR_B)[0], span_exact(cb, PAIR_A)[0]
            if sc[1].real < 0.005:
                continue
            rows.append((min(d_p00(ho, hb), d_tr(sc, sd)), m, J, delta, JDT, d_p00(ho, hb), d_tr(sc, sd), sc[1].real))
    rows.sort(key=lambda r: -r[0])
    for r in rows[:12]:
        print(f"   min {r[0]:.3f} m={r[1]} J={r[2]} DELTA={r[3]} dt={r[4]}: habit |dp00| {r[5]:.3f} control {r[6]:.3f} reach p11 {r[7]:.3f}")


# ------------------------------------------------------------------ what the design can and cannot see (--limits)
def nc_bond_gauge(G):
    """For a 4x4 two-site unitary that conserves the excitation number, the relative phase t = theta_(i+1) -
    theta_i with G^T = D G D^-1, D = diag(exp(i theta . n)); returns (t, max |G^T - D G D^-1|)."""
    t = np.angle(G[2, 1] / G[1, 2]) if abs(G[1, 2]) > 0 else 0.0
    Dl = np.diag(np.exp(1j * np.array([0.0, t, 0.0, t])))
    return t, float(np.abs(G.T - Dl @ G @ Dl.conj().T).max())


def fidelity_silence(rho, pair):
    """The silence-read Phi+ fidelity F_sil = (p00_sil + p11 + 2 Re c)/2, the operator swap's own observable
    (preparation and readout the same projector at mirrored pairs), and its shot-noise variance at SHOTS per setting
    (ZZ for the populations, XX and YY for Re c = (XX - YY)/4; the rest read in Z)."""
    rest = [q for q in range(N) if q not in pair]
    bits = [[(k >> (N - 1 - q)) & 1 for q in range(N)] for k in range(D)]
    silent = np.array([all(b[q] == 0 for q in rest) for b in bits], float)
    b0 = np.array([b[pair[0]] for b in bits])
    b1 = np.array([b[pair[1]] for b in bits])
    pz, px, py = (true_probs(rho, pair, st) for st in ("ZZ", "XX", "YY"))
    fz = 0.5 * silent * ((b0 == 0) & (b1 == 0) | (b0 == 1) & (b1 == 1))
    par = silent * (1 - 2 * ((b0 + b1) % 2))
    F = float(pz @ fz + (px @ par - py @ par) / 4)
    var = float(pz @ fz ** 2 - (pz @ fz) ** 2) + (float(px @ par ** 2 - (px @ par) ** 2)
                                                  + float(py @ par ** 2 - (py @ par) ** 2)) / 16
    return F, var / SHOTS


def limits(dev):
    rng = np.random.default_rng(SEED + 6)

    def nc_herm(eps):
        h = np.zeros((4, 4), complex)
        h[0, 0], h[3, 3] = rng.normal(), rng.normal()
        b = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        h[1:3, 1:3] = (b + b.conj().T) / 2
        return eps * h

    print("Limit 1, the tree gauge, with limit 2's free-chain rows. A two-site unitary that conserves the excitation number has |beta| = |gamma| in its")
    print("   one-excitation block, so G^T = D G D^-1 with D a relative site phase; on a path (a tree) the phases of all")
    print("   bonds glue into one global diagonal when each bond carries ONE gate. The transposed arm does; the main arm")
    print("   carries the E bonds at f = 1/2 and f = 1, two different gates. Noiseless otherwise; a fixed random")
    print("   number-conserving error after every bond gate (seeded), max over the steps:")
    for eps in (0.05, 0.15, 0.3):
        errs = [expm(-1j * nc_herm(eps)) for _ in range(N - 1)]
        gauge = max(nc_bond_gauge(errs[i] @ bond_unitary(i, 1.0, DELTA_ON))[1] for i in range(N - 1))
        rows = [exact_row(None, m, lambda i: two(errs[i], i)) for m in STEPS]
        print(f"   eps {eps:4.2f}: per-gate gauge residual {gauge:.1e}; transfer D transposed"
              f" {max(r['transposed'] for r in rows):.1e}, main {max(r['main'] for r in rows):.1e}, main0"
              f" {max(r['main0'] for r in rows):.1e}; habit |dp00| main0 {max(r['habit_main0'] for r in rows):.1e}"
              f" (complex c difference transposed {max(r['c_transposed'] for r in rows):.1e})")
    errs = [expm(-1j * nc_herm(0.05)) for _ in range(N - 1)]
    for label, bonds in (("on the O bonds only", ODD), ("on the E bonds only", EVEN)):
        rows = [exact_row(None, m, lambda i: two(errs[i], i) if i in bonds else np.eye(D)) for m in STEPS]
        print(f"   eps 0.05 {label}: transfer D main {max(r['main'] for r in rows):.1e}, transposed"
              f" {max(r['transposed'] for r in rows):.1e}")
    zz = [two(expm(-1j * 0.05 * np.kron(Z, Z)), b) for b in range(N - 1)]

    def spectator(i):
        C = np.eye(D, dtype=complex)
        for b in (i - 1, i + 1):
            if 0 <= b < N - 1:
                C = zz[b] @ C
        return C
    rows = [exact_row(None, m, spectator) for m in STEPS]
    print(f"   an error that reaches past the bond, ZZ of 0.05 rad on the neighbouring bonds after every bond gate:"
          f" transfer D main {max(r['main'] for r in rows):.1e}, transposed {max(r['transposed'] for r in rows):.1e}")
    out_seq = [0, 2, 1, 3]                       # one (E O) step, bonds in the order the gate applies them

    def seq_run(prep, m, order, err, before):
        rho = bell(prep)
        for i in order * m:
            if before:
                rho = err(rho, i)
            U = two(bond_unitary(i, 1.0, DELTA_ON), i)
            rho = U @ rho @ U.conj().T
            if not before:
                rho = err(rho, i)
        return rho
    zzn = [two(expm(-1j * 0.05 * np.kron(Z, Z)), b) for b in range(N - 1)]

    def err_spectator(rho, i):
        for b in (i - 1, i + 1):
            if 0 <= b < N - 1:
                rho = zzn[b] @ rho @ zzn[b].conj().T
        return rho

    def err_deph(rho, i):
        Zi = op(Z, i)
        return 0.99 * rho + 0.01 * Zi @ rho @ Zi
    print("   the way back that the transposition asks for is the exact reverse of the way out: the layers AND the bonds")
    print("   within each layer in reverse order, each error before its gate. The design's transposed arm reverses only")
    print("   the layers and keeps each error after its gate. Transposed transfer D, max over the steps:")
    for label, err in (("ZZ of 0.05 rad on the neighbouring bonds (static, past the bond)", err_spectator),
                       ("Z dephasing p = 0.01 on the bond's left site (stochastic, on the bond)", err_deph)):
        design = max(d_tr(span_exact(seq_run(PAIR_A, m, out_seq, err, False), PAIR_B)[0],
                          span_exact(seq_run(PAIR_B, m, [1, 3, 0, 2], err, False), PAIR_A)[0]) for m in STEPS)
        exact = max(d_tr(span_exact(seq_run(PAIR_A, m, out_seq, err, False), PAIR_B)[0],
                         span_exact(seq_run(PAIR_B, m, out_seq[::-1], err, True), PAIR_A)[0]) for m in STEPS)
        print(f"   {label}: the design's way back {design:.1e}, the exact reverse {exact:.1e}")
    print("   an XY - YX part of theta rad after every bond gate (not its own transpose), max over the steps:")
    xy = np.kron(X, Y) - np.kron(Y, X)
    for th in (0.05, 0.1, 0.2):
        C = expm(-1j * th / 2 * xy)
        rows = [exact_row(None, m, lambda i: two(C, i)) for m in STEPS]
        print(f"   theta {th:4.2f}: transfer D main {max(r['main'] for r in rows):.1e}, main0 {max(r['main0'] for r in rows):.1e},"
              f" transposed {max(r['transposed'] for r in rows):.1e}; habit |dp00| main0"
              f" {max(r['habit_main0'] for r in rows):.1e}, main at m = {STEPS[0]}, {STEPS[1]}: {rows[0]['habit_main']:.4f},"
              f" {rows[1]['habit_main']:.4f}")

    def evolve_delay(rho, m, d, tau):
        seq = [(EVEN, 0.5)]
        for k in range(m):
            seq += [(ODD, 1.0), (EVEN, 1.0 if k < m - 1 else 0.5)]
        for j, (bonds, f) in enumerate(seq):
            rho = layer(rho, bonds, f, DELTA_ON, d)
            if j < len(seq) - 1 and tau > 0:
                rho = idle(rho, tau, d)
        return rho

    def delay_row(d, m, tau):
        rf, rb = evolve_delay(bell(PAIR_A, d), m, d, tau), evolve_delay(bell(PAIR_B, d), m, d, tau)
        sf, sb = span_exact(rf, PAIR_B)[0], span_exact(rb, PAIR_A)[0]
        (ff, vf), (fb, vb) = fidelity_silence(rf, PAIR_B), fidelity_silence(rb, PAIR_A)
        return (sf[3].real - sb[3].real, d_tr(sf, sb), sf[3].real, sf[1].real, abs(sf[2]), ff - fb,
                float(np.sqrt(vf + vb)), sf[0].real - sb[0].real)

    print("Limit 4, T1 through the silence. The main arm with a delay of tau us between its layers (the sequence stays")
    print("   a palindrome). dP = P_silence forward - backward (tied only through number conservation), transfer D, and")
    print("   dF = F_sil forward - backward, the silence-read Phi+ fidelity (the operator swap's own observable), with")
    print("   its shot noise at 16384 shots per setting:")
    models = [("idealized, all", dev.idealized()), ("idealized, T1 off", dev.idealized().without(t1=True)),
              ("native, all", dev), ("native, T1 off", dev.without(t1=True))]
    for m in STEPS:
        for tau in (0.0, 5.0, 10.0, 20.0, 40.0):
            cells, reach = [], None
            for name, d in models:
                dp, dt, ps, p11, c, dF, sF, _ = delay_row(d, m, tau)
                cells.append(f"{name} {dp:+.1e} / {dt:.0e} / {dF:+.1e}")
                if name == "native, all":
                    reach = (ps, p11, c, sF)
            print(f"   m={m} tau {tau:4.1f}: " + "; ".join(cells) + f"  (native all, forward: P_silence {reach[0]:.3f},"
                  f" p11 {reach[1]:.3f}, |c| {reach[2]:.3f}; shot noise on dF {reach[3]:.1e})")
    sig = np.sqrt(2 * 0.25 / (len(SETTINGS) * SHOTS))
    print(f"   shot noise on dP (binomial at P = 1/2, pooled over the {len(SETTINGS)} settings at {SHOTS} shots, no readout"
          f" inversion): {sig:.1e}")
    dg = dev.idealized().without(t1=True, prep=True)
    print("   unital but not number-conserving (idealized placement, the gates' depolarizing and Z dephasing), dP / dF /"
          " dp00_sil at"
          " tau = 0, 10 us, m = 3: " + "; ".join("{:+.1e} / {:+.1e} / {:+.1e}".format(*(lambda r: (r[0], r[5], r[7]))(
              delay_row(dg, 3, tau))) for tau in (0.0, 10.0)))
    print("   which T1 profile sets it (idealized placement, no gate, preparation or readout error, T2* the line's mean,")
    print("   m = 3), dP at tau = 0, 5, 10 us, then dF and dp00_sil at tau = 10 us:")
    t2 = [float(np.mean(dev.T2star))] * N
    for label, t1 in [("uniform 240", [240.0] * 5), ("q4 at 137", [240.0, 240, 240, 240, 137]),
                      ("q0 at 137", [137.0, 240, 240, 240, 240]), ("q2 at 137", [240.0, 240, 137, 240, 240]),
                      ("palindromic 137, 240, 240, 240, 137", [137.0, 240, 240, 240, 137])]:
        d = Device(t1, t2, [0.0] * (N - 1), [0.0] * N, {PAIR_A: 0.0, PAIR_B: 0.0}, [(0, 0)] * N)
        d.native = d.gates = d.prep = False
        r10 = delay_row(d, 3, 10.0)
        print(f"   {label:36s} " + ", ".join(f"{delay_row(d, 3, tau)[0]:+.1e}" for tau in (0.0, 5.0)) +
              f", {r10[0]:+.1e} | dF {r10[5]:+.1e}, dp00_sil {r10[7]:+.1e}")
    saved = globals()["J"]
    globals()["J"] = (1.4, 2.2, 2.2, 1.4)
    d = Device([240.0] * N, t2, [0.0] * (N - 1), [0.0] * N, {PAIR_A: 0.0, PAIR_B: 0.0}, [(0, 0)] * N)
    d.native = d.gates = d.prep = False
    r10 = delay_row(d, 3, 10.0)
    globals()["J"] = saved
    print(f"   uniform 240 on the palindromic couplings J = (1.4, 2.2, 2.2, 1.4), tau = 10 us: dP {r10[0]:+.1e}, dF {r10[5]:+.1e}")
    zs = np.random.default_rng(11).normal(0, 0.05, N)
    for label, coh in (("a static Z profile after every gate (sd 0.05 rad, the Part 1 row)",
                        lambda i: expm(-1j * sum(zs[q] * op(Z, q) for q in range(N)) / 2)),
                       ("Rz(0.02) on the left qubit after every gate", lambda i: expm(-1j * 0.01 * op(Z, i)))):
        cells = []
        for arm in ("main", "transposed"):
            for m in STEPS:
                ro, rb = run_arm(arm, m, None, coh)
                cells.append(f"{arm} m={m} {fidelity_silence(ro, PAIR_B)[0] - fidelity_silence(rb, PAIR_A)[0]:+.1e}")
        print(f"   phases in F_sil, {label}, noiseless otherwise: dF " + ", ".join(cells))
    d = Device([240.0] * N, t2, [0.0] * (N - 1), [0.0] * N, {PAIR_A: 0.0, PAIR_B: 0.0}, [(0, 0)] * N, layer_us=3.0)
    d.native = d.gates = d.prep = False
    cells = []
    for arm in ("main", "control"):
        for m in STEPS:
            ro, rb = run_arm(arm, m, d)
            cells.append(f"{arm} m={m} {ro[0, 0].real - rb[0, 0].real:+.1e}")
    print("   the counting law: uniform T1 240 us, 3 us of idle per layer, dP_vac (the vacuum population, forward -"
          " backward): " + ", ".join(cells))
    zs = np.random.default_rng(11).normal(0, 0.05, N)
    cells = []
    for plabel, coh in (("none", None), ("Rz(0.02)", lambda i: expm(-1j * 0.01 * op(Z, i))),
                        ("static Z profile", lambda i: expm(-1j * sum(zs[q] * op(Z, q) for q in range(N)) / 2))):
        for tlabel, t1 in (("uniform 240", [240.0] * N), ("q4 at 137", [240.0, 240, 240, 240, 137])):
            d = Device(t1, t2, [0.0] * (N - 1), [0.0] * N, {PAIR_A: 0.0, PAIR_B: 0.0}, [(0, 0)] * N)
            d.native = d.gates = d.prep = False
            d.coherent = coh
            cells.append(f"{plabel}, {tlabel}: " + ", ".join(f"m={m} {delay_row(d, m, 10.0)[7]:+.2e}" for m in STEPS))
    print("   dp00_sil at tau = 10 us (idealized, no gate"
          " error), with an added phase error after every gate: " + "; ".join(cells))
    print(f"   shot noise on dp00_sil (ZZ setting alone, binomial at 1/2, both directions, no readout error):"
          f" {np.sqrt(2 * 0.25 / SHOTS):.1e}")
    base = dev.idealized().without(gates=True, prep=True)
    unif = Device([240.0] * N, list(dev.T2star), [0.0] * (N - 1), [0.0] * N, {PAIR_A: 0.0, PAIR_B: 0.0}, [(0, 0)] * N)
    unif.native = unif.gates = unif.prep = False
    cells = []
    for plabel, coh in (("none", None), ("Rz(0.02)", lambda i: expm(-1j * 0.01 * op(Z, i))),
                        ("static Z profile", lambda i: expm(-1j * sum(zs[q] * op(Z, q) for q in range(N)) / 2))):
        for tlabel, d0 in (("line T1", base), ("uniform T1 240", unif)):
            d = d0.without()
            d.coherent = coh
            cells.append(f"{plabel}, {tlabel}: " + ", ".join(f"m={m} {delay_row(d, m, 10.0)[5]:+.2e}" for m in STEPS))
    print("   phases at the signal's operating point: dF at tau = 10 us, idealized placement, the line's T2*, no gate"
          " or preparation error, with an added phase error after every gate: " + "; ".join(cells))
    print("The identity behind limits 3 and 4. For a real H with Z dephasing and thermal damping")
    print("   (cooling g_down, heating g_up per site), L - L^T = sum_l dg_l (s-_l . s+_l - s+_l . s-_l), with A . B the map")
    print("   rho -> A rho B, and L(I) = sum_l dg_l Z_l, dg_l = g_down_l - g_up_l (F84's net cooling); both vanish only at")
    print("   dg = 0. Five sites, a random integer symmetric H, random integer rates (row-stacked superoperators, exact):")
    sm = np.array([[0, 1], [0, 0]], complex)
    sp = sm.conj().T
    eye = np.eye(D)

    def sand(A, B):
        return np.kron(A, B.T)

    def dis(Jm):
        n = Jm.conj().T @ Jm
        return sand(Jm, Jm.conj().T) - 0.5 * np.kron(n, eye) - 0.5 * np.kron(eye, n.T)
    Hr = rng.integers(-3, 4, size=(D, D)).astype(float)
    Hr = Hr + Hr.T
    gd, gu, gz = (rng.integers(1, 9, N).astype(float), rng.integers(0, 5, N).astype(float),
                  rng.integers(1, 9, N).astype(float))
    L = -1j * (np.kron(Hr, eye) - np.kron(eye, Hr.T))
    for q in range(N):
        L = L + gd[q] * dis(op(sm, q)) + gu[q] * dis(op(sp, q)) + gz[q] * dis(op(Z, q))
    rhs = sum((gd[q] - gu[q]) * (sand(op(sm, q), op(sp, q)) - sand(op(sp, q), op(sm, q))) for q in range(N))
    LI = (L @ eye.reshape(-1)).reshape(D, D)
    print(f"   max |(L - L^T) - sum dg (s-.s+ - s+.s-)| = {np.abs((L - L.T) - rhs).max():.1e};"
          f" max |L(I) - sum dg Z| = {np.abs(LI - sum((gd[q] - gu[q]) * op(Z, q) for q in range(N))).max():.1e}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default=DEFAULT_INPUTS)
    ap.add_argument("--candidate", type=int, default=0)
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--certify", action="store_true")
    ap.add_argument("--limits", action="store_true", help="what the design can and cannot see")
    ap.add_argument("--no-power", action="store_true", help="skip Part 3 (for a fast day-of re-gate)")
    ap.add_argument("--write-constants", default=None)
    args = ap.parse_args()
    dev = load_device(args.inputs, args.candidate)
    print(f"channel_both_ways_gate v5: N={N}, pairs {PAIR_A}/{PAIR_B}, J={J}, DELTA_ON={DELTA_ON} and 0, dt={JDT},"
          f" steps {STEPS}, replica at m={REPLICA_STEP}, layer {LAYER_US} us, T2* = T2 x {T2STAR_FACTOR}")
    if args.search:
        search(dev)
    elif args.certify:
        certify()
    elif args.limits:
        limits(dev)
    else:
        part1(dev)
        special_rows(dev)
        t1_stress(dev)
        native_table(dev)
        rng = np.random.default_rng(SEED)
        constants, probs = part2(dev, rng, args.write_constants)
        if not args.no_power:
            part3(dev, rng, constants, probs)
