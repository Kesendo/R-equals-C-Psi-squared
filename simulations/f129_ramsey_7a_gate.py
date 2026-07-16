"""7a gate for IBM_F129_RAMSEY_FRINGE: the flown construction, simulated in Aer,
through the pinned section-5 estimator.

WHAT FLIES (per arm, compiled by the per-arm Givens compiler): GHZ cat on the
arm's differing modes (+ X layer, shared modes set to |1>) -> the 28-Givens
network (Floquet eigenbasis, per-arm column permutation) -> M Trotter steps of
the pure hopping chain -> inverse network -> un-GHZ -> seed-qubit X/Y readout.
Arms A0 (1,5,7)~(2,4,8), A1 (1,5,7)~(2,5,6), A2 (1,5,7)~(1,6,7); theta = 0.5,
M in 0..8, 16384 shots per circuit, plus 2 readout-cal circuits.

ANCHORED NOISE MODEL (Confirmations 21 + design-doc constants, v2.4 timing):
depolarizing p2 = 0.3% on CX + thermal relaxation T1 = 200 us / T2 = 70 us with
Heron-era durations (CX 150 ns, 1q 50 ns; NO idle padding, tau_step = 0.7 us
wall) + readout error 1% (corrected by the cal PUBs) + COHERENT site-dependent
NN ZZ crosstalk zeta_i ~ 3.8 kHz +-30% (fixed seed), injected as noise-free RZZ
each step in the Ramsey-shift convention H_ZZ = pi zeta Z_i Z_j (the convention
of the price_pair conditional-Ramsey measurement; rzz(2 pi zeta tau)).

ESTIMATOR (verbatim section 5): readout-corrected <X>, <Y> -> phi = atan2, V;
per-arm weighted LSQ with PREDICTED-V weights; theta-hat from A2 (pinned linear
inversion); clause (a) A0 within 3 sigma_a of the theta-hat center; clause (b)
A2 vs A0 >= 5 sigma combined; clause (c) A2 vs nominal prediction.

REQUIRED to record 7a: CONFIRMED fires with >= 3x power margin on (a) and (b);
frozen null band (parametric bootstrap under H0 through this estimator);
hierarchical bootstrap SEs vs analytic. Plus the six named checks (run with
--checks): branch decoherence, quasi-static vs Markovian T2*, mixture model,
branch-mixing zero end-to-end, readout correction efficacy, leakage bias.

Doc: experiments/IBM_F129_RAMSEY_FRINGE.md (the 7a RECORDED block). Seeds pinned.
"""

import sys
import numpy as np
from itertools import combinations

sys.path.insert(0, "simulations")
from f129_ramsey_fringe_design import zz_drift_exact  # Wick ZZ, shared pinned timing  # noqa: E402

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate, XXPlusYYGate
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (NoiseModel, depolarizing_error,
                              thermal_relaxation_error, ReadoutError)

N, THETA, MMAX, SHOTS = 8, 0.5, 8, 16384
P2, T1_US, T2S_US = 0.003, 200.0, 70.0
CX_NS, Q1_NS, TAU_STEP_US = 150.0, 50.0, 0.7   # Heron-era 2q, NO idle padding (v2.4); reconciled at 7b
ZETA_BASE_KHZ, ZETA_SPREAD = 3.8, 0.30
RO_ERR = 0.01
SEED = 20260715
ARMS = [("A0", (1, 5, 7), (2, 4, 8)),
        ("A1", (1, 5, 7), (2, 5, 6)),
        ("A2", (1, 5, 7), (1, 6, 7))]
M_DEPH = {"A0": 3.680, "A1": 3.78, "A2": 3.77}   # design-gate constants (weights)
B_QS = 0.008    # two-sided quasi-static disorder budget (checks(): -0.0075 +- 0.0037)
# B_ZZ2 is DERIVED and printed by zz2_scan(): the chiral second-order law at 2x zeta

rng = np.random.default_rng(SEED)
ZETA_KHZ = ZETA_BASE_KHZ * (1 + ZETA_SPREAD * (2 * rng.random(N - 1) - 1))
DETUNE_KHZ = None   # site-dependent single-qubit Z detuning (checks)


# ---------- single-particle layer (same as the design gate) ----------
def one_particle_step(theta):
    U = np.eye(N, dtype=complex)
    for parity in (0, 1):
        L = np.eye(N, dtype=complex)
        for i in range(parity, N - 1, 2):
            B = np.eye(N, dtype=complex)
            B[i, i] = B[i + 1, i + 1] = np.cos(theta)
            B[i, i + 1] = B[i + 1, i] = 1j * np.sin(theta)
            L = B @ L
        U = L @ U
    return U


def floquet(theta):
    w, V = np.linalg.eig(one_particle_step(theta))
    ph = np.angle(w)
    order = np.argsort(-ph)
    return ph[order], V[:, order]


def givens_decompose(A):
    A = A.copy().astype(complex)
    rots = []
    for j in range(N - 1):
        for i in range(N - 1, j, -1):
            a, b = A[i - 1, j], A[i, j]
            if abs(b) < 1e-15:
                continue
            r = np.hypot(abs(a), abs(b))
            G = np.array([[np.conj(a) / r, np.conj(b) / r],
                          [-b / r, a / r]], dtype=complex)
            A[i - 1:i + 1, :] = G @ A[i - 1:i + 1, :]
            rots.append((i - 1, G))
    return rots, np.diag(A).copy()


class ArmNetwork:
    def __init__(self, V, a, b):
        self.a, self.b = a, b
        a_priv = [k for k in a if k not in b]
        b_priv = [k for k in b if k not in a]
        shared = [k for k in a if k in b]
        rest = [k for k in range(1, N + 1) if k not in set(a) | set(b)]
        self.mode_of_qubit = a_priv + b_priv + shared + rest
        self.n_apriv, self.cat = len(a_priv), len(a_priv) + len(b_priv)
        self.shared_qubits = list(range(self.cat, self.cat + len(shared)))
        perm = [m - 1 for m in self.mode_of_qubit]
        self.rots, self.phases = givens_decompose(V[:, perm])
        self.bits_a = sum(1 << self.mode_of_qubit.index(k) for k in a)
        self.bits_b = sum(1 << self.mode_of_qubit.index(k) for k in b)


# ---------- circuit construction ----------
def nc_circuit(B):
    """Number-conserving 2-qubit gate as NATIVE-decomposable primitives:
    single-qubit phases + XXPlusYYGate (2 CX under transpile). B is the SU(2)
    block [[alpha, beta], [-conj(beta), conj(alpha)]] acting on
    (e_p, e_{p+1}) = (|01>, |10>) in little-endian; |00>, |11> untouched
    (phase-gate |11> contributions cancel since a+b+c+d = 0)."""
    alpha, beta = B[0, 0], B[0, 1]
    th = np.arctan2(np.hypot(abs(beta), 0.0), abs(alpha))
    aa, bb = np.angle(alpha), (np.angle(beta) if abs(beta) > 1e-15 else 0.0)
    # left phases (a on qubit p+1, b on qubit p), right phases (c = 0, d on p)
    b_ph = bb - np.pi / 2
    d_ph = aa - b_ph
    a_ph = -aa
    qc = QuantumCircuit(2)
    qc.p(d_ph, 0)
    qc.append(XXPlusYYGate(-2 * th, 0.0), [0, 1])
    qc.p(b_ph, 0)
    qc.p(a_ph, 1)
    return qc.to_gate(label="ncg")


def _nc_selfcheck():
    """The decomposition must equal the raw 4x4 for random SU(2) blocks."""
    from qiskit.quantum_info import Operator
    rng_c = np.random.default_rng(7)
    for _ in range(25):
        x = rng_c.normal(size=4)
        x /= np.linalg.norm(x)
        alpha = x[0] + 1j * x[1]
        beta = x[2] + 1j * x[3]
        B = np.array([[alpha, beta], [-np.conj(beta), np.conj(alpha)]])
        U = np.eye(4, dtype=complex)
        U[1, 1], U[1, 2] = B[0, 0], B[0, 1]
        U[2, 1], U[2, 2] = B[1, 0], B[1, 1]
        got = Operator(nc_circuit(B)).data
        assert np.max(np.abs(got - U)) < 1e-12, "nc decomposition broken"


_nc_selfcheck()


HOP = np.array([[np.cos(THETA), 1j * np.sin(THETA)],
                [1j * np.sin(THETA), np.cos(THETA)]])


def build_arm_circuit(net, m_steps, quad, zz=True):
    qc = QuantumCircuit(N, 1)
    # seed: GHZ over the cat block + X layers
    qc.h(0)
    for q in range(1, net.cat):
        qc.cx(q - 1, q)
    for q in range(net.n_apriv, net.cat):  # b-private bits set on the |0...> branch
        qc.x(q)                            # (seed |0> = branch b: measured slope = +dPhi)
    for q in net.shared_qubits:
        qc.x(q)
    # network (phase layer D = product of SINGLE-qubit phases e^{i phi_k n_k})
    for k in range(N):
        qc.p(np.angle(net.phases[k]), k)
    for p, G in reversed(net.rots):
        qc.append(nc_circuit(G.conj().T), [p, p + 1])
    # M Trotter steps (+ coherent site-dependent ZZ crosstalk + idle padding)
    for _ in range(m_steps):
        for parity in (0, 1):
            for i in range(parity, N - 1, 2):
                qc.append(nc_circuit(HOP), [i, i + 1])
        if zz:
            for i in range(N - 1):
                qc.rzz(2 * np.pi * ZETA_KHZ[i] * 1e3 * TAU_STEP_US * 1e-6, i, i + 1)
        if DETUNE_KHZ is not None:
            for i in range(N):
                qc.rz(2 * np.pi * DETUNE_KHZ[i] * 1e3 * TAU_STEP_US * 1e-6, i)
    # inverse network
    for p, G in net.rots:
        qc.append(nc_circuit(G), [p, p + 1])
    for k in range(N):
        qc.p(-np.angle(net.phases[k]), k)
    # un-GHZ
    for q in range(net.n_apriv, net.cat):
        qc.x(q)
    for q in range(net.cat - 1, 0, -1):
        qc.cx(q - 1, q)
    # quadrature readout on the seed qubit
    if quad == "Y":
        qc.sdg(0)
    qc.h(0)
    qc.measure(0, 0)
    return qc


# ---------- noise model ----------
def make_noise():
    nm = NoiseModel()
    t1, t2 = T1_US * 1e3, T2S_US * 1e3          # ns
    err1 = thermal_relaxation_error(t1, t2, Q1_NS)
    err2 = thermal_relaxation_error(t1, t2, CX_NS).tensor(
        thermal_relaxation_error(t1, t2, CX_NS)).compose(
        depolarizing_error(P2, 2))
    nm.add_all_qubit_quantum_error(err1, ["sx", "x"])  # rz virtual, error-free
    nm.add_all_qubit_quantum_error(err2, ["cx"])
    nm.add_all_qubit_readout_error(ReadoutError([[1 - RO_ERR, RO_ERR],
                                                 [RO_ERR, 1 - RO_ERR]]))
    return nm


# ---------- run + estimator ----------
def run_all(noise=True, zz=True, shots=SHOTS):
    ph, V = floquet(THETA)
    nets = {name: ArmNetwork(V, a, b) for name, a, b in ARMS}
    sim = AerSimulator(noise_model=make_noise() if noise else None,
                       seed_simulator=SEED)
    circuits, keys = [], []
    for name, a, b in ARMS:
        for m in range(MMAX + 1):
            for quad in ("X", "Y"):
                circuits.append(build_arm_circuit(nets[name], m, quad, zz=zz))
                keys.append((name, m, quad))
    # readout cal circuits
    cal0 = QuantumCircuit(N, 1); cal0.measure(0, 0)
    cal1 = QuantumCircuit(N, 1); cal1.x(0); cal1.measure(0, 0)
    circuits += [cal0, cal1]
    keys += [("CAL", 0, "Z"), ("CAL", 1, "Z")]
    tqc = transpile(circuits, basis_gates=["cx", "rz", "sx", "x", "delay", "rzz"],
                    optimization_level=1, seed_transpiler=SEED)
    cx_m8 = tqc[keys.index(("A0", MMAX, "X"))].count_ops().get("cx", 0)
    print(f"transpiled CX count, A0 arm at M = {MMAX}: {cx_m8} (budget model 234)")
    result = sim.run(tqc, shots=shots).result()
    counts = {k: result.get_counts(i) for i, k in enumerate(keys)}
    return counts, ph


def expval(counts):
    n0 = counts.get("0", 0)
    n1 = counts.get("1", 0)
    return (n0 - n1) / (n0 + n1)


def readout_matrix(counts):
    p0g0 = counts[("CAL", 0, "Z")].get("0", 0) / SHOTS
    p1g1 = counts[("CAL", 1, "Z")].get("1", 0) / SHOTS
    return np.array([[p0g0, 1 - p1g1], [1 - p0g0, p1g1]])


def corrected_exp(counts, key, R):
    c = counts[key]
    raw = np.array([c.get("0", 0), c.get("1", 0)], dtype=float)
    raw /= raw.sum()
    p = np.linalg.solve(R, raw)
    return p[0] - p[1]


def estimate(counts, use_correction=True):
    """Section-5 estimator: per-arm slopes with predicted-V weights."""
    R = readout_matrix(counts) if use_correction else np.eye(2)
    out = {}
    for name, a, b in ARMS:
        phis, vs = [], []
        for m in range(MMAX + 1):
            ex = corrected_exp(counts, (name, m, "X"), R)
            ey = corrected_exp(counts, (name, m, "Y"), R)
            phis.append(np.arctan2(ey, ex))
            vs.append(np.hypot(ex, ey))
        phis = np.unwrap(np.array(phis))
        ms = np.arange(MMAX + 1)
        sig = 1.0 / (np.maximum(pred_V(ms, M_DEPH[name]), 1e-3) * np.sqrt(SHOTS))
        w = 1.0 / sig ** 2
        X = np.vstack([ms, np.ones_like(ms)]).T
        W = np.diag(w)
        cov = np.linalg.inv(X.T @ W @ X)
        beta = cov @ X.T @ W @ phis
        out[name] = {"slope": beta[0], "sigma": np.sqrt(cov[0, 0]),
                     "phis": phis, "vs": np.array(vs)}
    return out


def pred_V(ms, m_deph):
    givens, bonds = N * (N - 1) // 2, N - 1
    cx = 5 + 2 * givens * 2 + ms * bonds * 2 + 5
    t_us = ms * TAU_STEP_US + 2.5   # step wall time + networks/seed (150 ns 2q)
    # 0.75 = the 7a-re-frozen 1q-gate attenuation factor (transpiled sx count
    # 336 + 42 M at 50 ns each; sim/old-model ratio 0.79-0.82, pinned conservative)
    return 0.75 * (1 - P2) ** cx * np.exp(-3 * t_us / T1_US) * np.exp(-m_deph * t_us / T2S_US)


def zz2_budget():
    """(2026-07-16: the law is now DERIVED, docs/proofs/PROOF_ZETA2_ANTI_
    PROTECTION.md + gate simulations/zeta2_anti_protection.py: second-order
    Floquet PT under the antiunitary chiral lift Theta = T.K; clean-limit
    coefficient 0.002522, which the flown estimator dresses to C2LAW below.
    Theta-mirror pairs only.)

    The one-sided clause-(a) budget from the chiral second-order law,
    derived by zz2_scan() (exact statevector): bias = C2LAW * (zeta/3.8 kHz)^2
    * (tau_step/1.2 us)^2 (1.2 us = the law's normalization point), budgeted at
    2x the transferred zeta with a 1.22 estimator-pattern uplift (the seed
    pattern gives the law value; a uniform 7.6 kHz profile gives up to 12%
    more). One-sidedness is a property of the SAME-SIGN always-on regime
    (checked by zz2_scan over all-positive draws), NOT of the chiral symmetry
    (which fixes the factor 2, not the sign); mixed-sign excursions (gate-printed
    by zz2_scan, worst ~-8e-4 at |zeta| <= 7.6 kHz) are covered by the two-sided
    B_QS. The 1.22 uplift over-covers the measured uniform-profile excess (~6%,
    gate-printed) plus rounding headroom; every slack is in the safe direction."""
    C2LAW = 0.00257   # gate-printed by zz2_scan() under the FLOWN estimator
    # weights (the law is estimator-defined: the v2.4 predicted-V weights moved
    # it from the equal-era 0.00242; zz2_scan is the authority)
    return C2LAW * (TAU_STEP_US / 1.2) ** 2 * 2.0 ** 2 * 1.22


def verdict(est, ph):
    """The pinned clauses (doc section 5, v2.4): budgeted asymmetric window."""
    lamsum = lambda t: sum(ph[k - 1] for k in t)
    drift_A0 = lamsum((1, 5, 7)) - lamsum((2, 4, 8))     # +0.0213
    slope_A2_nom = lamsum((1, 5, 7)) - lamsum((1, 6, 7))  # +0.3204
    zz_A2 = zz_drift_exact(THETA, (1, 5, 7), (1, 6, 7), N)   # Wick, current timing
    b_zz2 = zz2_budget()
    # theta-hat rule (pinned linear inversion)
    theta_hat = THETA + (est["A2"]["slope"] - slope_A2_nom) / 0.6164
    center = drift_A0 + 0.1294 * (theta_hat - THETA)
    # sigma_a: frozen projection inflated by the propagated theta-hat variance
    infl = np.sqrt(1 + (0.1294 * (est["A2"]["sigma"] / 0.6164) / est["A0"]["sigma"]) ** 2)
    sig_a = infl * est["A0"]["sigma"]
    dev = est["A0"]["slope"] - center
    dev_a = dev / sig_a
    # clause (a): center + [-(3 sig_a + B_QS), +(3 sig_a + b_zz2 + B_QS)]
    a_ok = -(3 * sig_a + B_QS) <= dev <= (3 * sig_a + b_zz2 + B_QS)
    # VIOLATED: outside center + [-(5 sig_a + B_QS), +(5 sig_a + b_zz2 + B_QS)]
    violated_window = -(5 * sig_a + B_QS) <= dev <= (5 * sig_a + b_zz2 + B_QS)
    # clause (c): A2 vs NOMINAL prediction
    c_ok = abs(est["A2"]["slope"] - (slope_A2_nom + zz_A2)) <= \
        3 * est["A2"]["sigma"] + 5 * abs(zz_A2)
    # clause (b)
    comb = np.hypot(est["A0"]["sigma"], est["A2"]["sigma"])
    sep_b = abs(est["A2"]["slope"] - est["A0"]["slope"]) / comb
    b_ok = sep_b >= 5 and np.sign(est["A2"]["slope"] - est["A0"]["slope"]) > 0
    conf = a_ok and b_ok and c_ok
    if conf:
        outcome = "CONFIRMED"
    elif not c_ok:
        outcome = "VOID (clause c)"
    elif not violated_window and c_ok:
        outcome = "VIOLATED"
    else:
        outcome = "inconclusive"
    return {"theta_hat": theta_hat, "center": center, "dev_a": dev_a,
            "sig_a": sig_a, "infl": infl, "b_zz2": b_zz2, "zz_A2": zz_A2,
            "sep_b": sep_b, "a": a_ok, "b": b_ok, "c": c_ok,
            "CONFIRMED": conf, "outcome": outcome,
            "margin_a": None, "margin_b": sep_b / 5}


def zz2_scan():
    """Exact statevector derivation of the chiral second-order ZZ law: the fitted
    A0 slope bias vs zeta scale, no shots, no noise. Prints the quadratic
    coefficient (the C2LAW constant of zz2_budget) and the budget."""
    ph, V = floquet(THETA)
    net = ArmNetwork(V, (1, 5, 7), (2, 4, 8))
    from qiskit.quantum_info import Statevector
    dPhi = sum(ph[k - 1] for k in (1, 5, 7)) - sum(ph[k - 1] for k in (2, 4, 8))
    print("zeta scale | exact fitted-slope bias | bias/scale^2 (at tau_step = 1.2 us ref)")
    ref_tau = 1.2
    for scale in (0.5, 1.0, 2.0):
        phis = []
        for m in range(MMAX + 1):
            qc = build_arm_circuit_bare(net, m, scale, ref_tau)
            sv = Statevector.from_instruction(qc)
            aA = sv.data[net.bits_a]
            aB = sv.data[net.bits_b]
            phis.append(np.angle(aA * np.conj(aB)))
        phis = np.unwrap(np.array(phis))
        ms = np.arange(MMAX + 1)
        sig = 1.0 / (np.maximum(pred_V(ms, M_DEPH["A0"]), 1e-3) * np.sqrt(SHOTS))
        w = 1.0 / sig ** 2
        X = np.vstack([ms, np.ones_like(ms)]).T
        W = np.diag(w)
        slope = (np.linalg.inv(X.T @ W @ X) @ X.T @ W @ phis)[0]
        bias = slope - dPhi
        print(f"  x{scale:4.2f}   {bias:+.6f}   {bias / scale ** 2:+.6f}")
    print(f"B_ZZ2 at the flown tau_step = {TAU_STEP_US} us, 2x zeta, 1.22 uplift: "
          f"{zz2_budget():+.5f}")

    # one-sidedness holds in the same-sign regime: minimum bias over positive draws
    def bias_for(zetas):
        phis = []
        for m in range(MMAX + 1):
            qc = build_arm_circuit_bare(net, m, 1.0, ref_tau, zetas)
            sv = Statevector.from_instruction(qc)
            phis.append(np.angle(sv.data[net.bits_a] * np.conj(sv.data[net.bits_b])))
        phis = np.unwrap(np.array(phis))
        slope = (np.linalg.inv(X.T @ W @ X) @ X.T @ W @ phis)[0]
        return slope - dPhi

    srng = np.random.default_rng(SEED + 5)
    biases = [bias_for(srng.uniform(1.0, 7.6, N - 1)) for _ in range(40)]
    print(f"same-sign check (40 positive zeta draws in [1.0, 7.6] kHz): "
          f"min bias {min(biases):+.5f}, max {max(biases):+.5f} "
          f"(one-sidedness holds iff min > 0)")
    uni = bias_for(np.full(N - 1, 7.6))
    print(f"uniform 7.6 kHz profile: bias {uni:+.5f} "
          f"(vs law x4 = {0.00257 * 4:+.5f}; the 1.22 uplift over-covers)")
    mixed = [bias_for(srng.uniform(-7.6, 7.6, N - 1)) for _ in range(40)]
    print(f"mixed-sign check (40 draws, |zeta| <= 7.6 kHz): "
          f"min bias {min(mixed):+.5f} (downward excursions fall under B_QS)")


def build_arm_circuit_bare(net, m_steps, zeta_scale, tau_ref, zetas=None):
    """Noise-free circuit variant for the zz2 scan (no measure, zz at scale or
    an explicit per-bond zeta vector in kHz)."""
    qc = QuantumCircuit(N)
    qc.h(0)
    for q in range(1, net.cat):
        qc.cx(q - 1, q)
    for q in range(net.n_apriv, net.cat):
        qc.x(q)
    for q in net.shared_qubits:
        qc.x(q)
    for k in range(N):
        qc.p(np.angle(net.phases[k]), k)
    for pp, G in reversed(net.rots):
        qc.append(nc_circuit(G.conj().T), [pp, pp + 1])
    for _ in range(m_steps):
        for parity in (0, 1):
            for i in range(parity, N - 1, 2):
                qc.append(nc_circuit(HOP), [i, i + 1])
        zv = ZETA_KHZ * zeta_scale if zetas is None else np.asarray(zetas)
        for i in range(N - 1):
            qc.rzz(2 * np.pi * zv[i] * 1e3 * tau_ref * 1e-6, i, i + 1)
    for pp, G in net.rots:
        qc.append(nc_circuit(G), [pp, pp + 1])
    for k in range(N):
        qc.p(-np.angle(net.phases[k]), k)
    # NO un-GHZ here: the scan reads the two branch amplitudes directly at
    # net.bits_a / net.bits_b after the inverse network
    return qc


def verdict_logic_selfcheck(est, ph):
    """Exercise the VIOLATED and dead-band paths on synthetic A0 slopes."""
    import copy
    print("verdict-logic self-check (synthetic A0 slopes through the pinned rule):")
    # window edges at the recorded run: upper CONFIRM edge ~ +0.0225, upper
    # VIOLATED edge ~ +0.0293; lower edges ~ -0.0182 / -0.0250
    for delta, want in ((0.0, "CONFIRMED"), (0.025, "inconclusive"),
                        (-0.022, "inconclusive"),
                        (0.05, "VIOLATED"), (-0.05, "VIOLATED")):
        e = copy.deepcopy(est)
        v0 = verdict(est, ph)
        e["A0"]["slope"] = v0["center"] + delta
        v = verdict(e, ph)
        ok = v["outcome"] == want
        print(f"  A0 = center {delta:+.3f}: outcome {v['outcome']:14s} "
              f"(expected {want}) {'OK' if ok else 'MISMATCH'}")


def bootstrap_and_nullband(counts, ph, n_boot=400):
    """Hierarchical bootstrap SEs (resample per-point binomials) and the frozen
    null band for clause (a): parametric bootstrap under H0 (A0 slope exactly at
    the theta-hat center) through THIS estimator."""
    brng = np.random.default_rng(SEED + 1)
    R = readout_matrix(counts)
    est0 = estimate(counts)
    v0 = verdict(est0, ph)

    def wls_slope(phis, name):
        phis = np.unwrap(np.array(phis))
        ms = np.arange(MMAX + 1)
        sig = 1.0 / (np.maximum(pred_V(ms, M_DEPH[name]), 1e-3) * np.sqrt(SHOTS))
        w = 1.0 / sig ** 2
        X = np.vstack([ms, np.ones_like(ms)]).T
        W = np.diag(w)
        return (np.linalg.inv(X.T @ W @ X) @ X.T @ W @ phis)[0]

    def resample_arm(name):
        phis = []
        for m in range(MMAX + 1):
            pt = []
            for quad in ("X", "Y"):
                c = counts[(name, m, quad)]
                p1 = c.get("1", 0) / SHOTS
                n1 = brng.binomial(SHOTS, p1)
                raw = np.array([SHOTS - n1, n1], dtype=float) / SHOTS
                pcorr = np.linalg.solve(R, raw)
                pt.append(pcorr[0] - pcorr[1])
            phis.append(np.arctan2(pt[1], pt[0]))
        return wls_slope(phis, name)

    print("hierarchical bootstrap vs analytic sigma_slope:")
    for n in ("A0", "A1", "A2"):
        b = np.array([resample_arm(n) for _ in range(n_boot)])
        print(f"  {n}: bootstrap SE = {b.std():.5f}, analytic = {est0[n]['sigma']:.5f}, "
              f"ratio = {b.std() / est0[n]['sigma']:.2f}")

    def null_slope():
        phis = []
        for m in range(MMAX + 1):
            phi_true = v0["center"] * m + np.pi / 4
            vmodel = pred_V(np.array([m]), M_DEPH["A0"])[0]
            pt = []
            for val in (np.cos(phi_true), np.sin(phi_true)):
                pexp = (1 + vmodel * val) / 2
                praw = R @ np.array([pexp, 1 - pexp])
                n0 = brng.binomial(SHOTS, praw[0])
                raw = np.array([n0, SHOTS - n0], dtype=float) / SHOTS
                pcorr = np.linalg.solve(R, raw)
                pt.append(pcorr[0] - pcorr[1])
            phis.append(np.arctan2(pt[1], pt[0]))
        return wls_slope(phis, "A0") - v0["center"]

    nulls = np.array([null_slope() for _ in range(n_boot)])
    lo, hi = np.quantile(nulls, [0.0013, 0.9987])
    print(f"null band (H0 slope-minus-center, 3sigma quantiles, {n_boot} draws): "
          f"[{lo:+.5f}, {hi:+.5f}] rad/step "
          f"(analytic 3 sigma_a = {3 * v0['sig_a']:.5f})")


def checks():
    """The six named 7a checks (readout, mixture, protection, quasi-static)."""
    global DETUNE_KHZ
    print()
    print("-- checks: readout correction / mixture model / protection --")
    counts, _ = run_all()
    est_c = estimate(counts, use_correction=True)
    est_n = estimate(counts, use_correction=False)
    print(f"  readout correction: A0 slope {est_c['A0']['slope']:+.5f} corrected vs "
          f"{est_n['A0']['slope']:+.5f} uncorrected "
          f"(V(8): {est_c['A0']['vs'][-1]:.3f} vs {est_n['A0']['vs'][-1]:.3f})")
    from qiskit import transpile as _tp
    ph_v, V_v = floquet(THETA)
    net_v = ArmNetwork(V_v, (1, 5, 7), (2, 4, 8))
    t8 = _tp(build_arm_circuit(net_v, MMAX, "X"),
             basis_gates=["cx", "rz", "sx", "x", "rzz"], optimization_level=1,
             seed_transpiler=SEED)
    print(f"  transpiled 1q count at M = {MMAX}: sx = {t8.count_ops().get('sx', 0)} "
          f"(the 0.75 V-model factor's provenance: sx = 336 + 42 M at 50 ns)")
    print("  mixture model, arm A0 (V_sim must sit ABOVE the conservative model, every M):")
    worst_gap = min(est_c["A0"]["vs"][m] - pred_V(np.array([m]), M_DEPH["A0"])[0]
                    for m in range(MMAX + 1))
    print(f"    min over M = 0..{MMAX} of (V_sim - V_model): {worst_gap:+.3f} "
          f"({'conservative OK' if worst_gap >= 0 else 'MODEL OPTIMISTIC'})")
    for m in (0, 4, 8):
        vs, vm = est_c["A0"]["vs"][m], pred_V(np.array([m]), M_DEPH["A0"])[0]
        print(f"    M = {m}: sim {vs:.3f}, model {vm:.3f} "
              f"({'conservative OK' if vs >= vm else 'MODEL OPTIMISTIC'})")
    DETUNE_KHZ = 2.0 * (2 * np.random.default_rng(SEED + 2).random(N) - 1)
    counts_d, _ = run_all()
    est_d = estimate(counts_d)
    DETUNE_KHZ = None
    counts_0, _ = run_all(zz=False)
    est_0 = estimate(counts_0)
    counts_nl, _ = run_all(noise=False, zz=False)
    est_nl = estimate(counts_nl, use_correction=False)
    ph_loc, _ = floquet(THETA)
    dphi0 = sum(ph_loc[k - 1] for k in (1, 5, 7)) - sum(ph_loc[k - 1] for k in (2, 4, 8))
    print(f"  noiseless estimator bias (sampling only): "
          f"{est_nl['A0']['slope'] - dphi0:+.5f} rad/step")
    print(f"  protection end-to-end: A0 slope zz+detuning {est_d['A0']['slope']:+.5f}, "
          f"zz only {est_c['A0']['slope']:+.5f}, clean {est_0['A0']['slope']:+.5f} "
          f"(zz/detuning shift {abs(est_d['A0']['slope'] - est_0['A0']['slope']) / est_0['A0']['sigma']:.2f} sigma)")
    print("-- check: quasi-static site detunings (10 draws, 2048 shots each) --")
    qrng = np.random.default_rng(SEED + 3)
    slopes = []
    for _ in range(10):
        DETUNE_KHZ = 3.0 * qrng.normal(size=N)
        c, _ = run_all(shots=2048)
        slopes.append(estimate(c)["A0"]["slope"])
    DETUNE_KHZ = None
    slopes = np.array(slopes)
    print(f"  A0 slope over quasi-static draws: mean {slopes.mean():+.5f}, "
          f"spread {slopes.std():.5f} (Markovian-run value {est_c['A0']['slope']:+.5f})")


def seeds_stability(ph, n_seeds=5):
    print()
    print("-- multi-seed stability of the verdict --")
    global SEED
    base = SEED
    for k in range(n_seeds):
        SEED = base + 100 + k
        counts, _ = run_all()
        v = verdict(estimate(counts), ph)
        print(f"  seed {SEED}: dev_a = {v['dev_a']:+.2f}, sep_b = {v['sep_b']:.1f}, "
              f"abc = {'P' if v['a'] else 'F'}{'P' if v['b'] else 'F'}"
              f"{'P' if v['c'] else 'F'} -> "
              f"{v['outcome']}")
    SEED = base


def main():
    print(f"zeta_i (kHz, seed {SEED}): {np.round(ZETA_KHZ, 2)}")
    counts, ph = run_all()
    est = estimate(counts)
    for name in ("A0", "A1", "A2"):
        e = est[name]
        print(f"{name}: slope = {e['slope']:+.5f} +- {e['sigma']:.5f}, "
              f"V(M=8) = {e['vs'][-1]:.3f} (model {pred_V(np.array([8]), M_DEPH[name])[0]:.3f})")
    v = verdict(est, ph)
    # power margin for clause (a): the impostor distance in realized sigma
    imp = 0.1092
    v["margin_a"] = abs(imp - v["center"]) / v["sig_a"] / 5
    print(f"theta_hat = {v['theta_hat']:.4f}, A0 center = {v['center']:+.5f}, "
          f"A0 deviation = {v['dev_a']:+.2f} sigma_a (sig_a = {v['sig_a']:.5f} "
          f"= {v['infl']:.3f} x sigma_slope)")
    print(f"budgets: b_zz2 = +{v['b_zz2']:.5f} (one-sided, zz2_scan law at 2x zeta), "
          f"b_qs = +-{B_QS:.3f}; zz_A2 (Wick, current timing) = {v['zz_A2']:+.5f}")
    print(f"clause (a) {'PASS' if v['a'] else 'FAIL'}, "
          f"(b) {'PASS' if v['b'] else 'FAIL'} (sep {v['sep_b']:.1f} sigma, margin {v['margin_b']:.1f}x), "
          f"(c) {'PASS' if v['c'] else 'FAIL'}")
    print(f"clause (a) impostor power margin: {v['margin_a']:.1f}x (need >= 3)")
    print(f"7a VERDICT: {v['outcome']}")
    if "--full" in sys.argv:
        verdict_logic_selfcheck(est, ph)
        zz2_scan()
        bootstrap_and_nullband(counts, ph)
        checks()
        seeds_stability(ph)


if __name__ == "__main__":
    main()
