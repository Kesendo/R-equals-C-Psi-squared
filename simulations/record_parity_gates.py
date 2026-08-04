"""The record-parity single-source module (built 2026-08-04, flight called by Tom).

SINGLE SOURCE for circuits, estimator, and verdict of the record-parity flight
(experiments/RECORD_PARITY_HARDWARE_PREDICTION.md). The sim gate (the house 7a),
the counts-level gate (7b), and the hardware runner all import THIS module; nothing
is reimplemented elsewhere. PROMOTED 2026-08-04 at the round-31 fold (prior-work
M4: the pre-reg's rehearsal numbers are this module's outputs, and the pinned
"empty round on the frozen numbers" is unexecutable against an artifact outside
the repo — F129 committed its gate scripts as part of pre-registration, the
precedent this follows); enters the repo with the next commit, before the binding
freeze.

Stage 1 of the build: the exact 4-qubit circuit physics (statevector, no sampling)
and the distribution-level estimator chain, self-tested against the law's ideal
face: rho_hat(r) = cos(r*pi/2) SIGNED, the record exactly on Y, the r = 2 flip.
Stage 2: the sampling layer (shots, pinned even/odd split at the r = 0 arms).
Stage 3 (this stage): the injection bank + the generative models, per the
pre-reg's "Sim gate specification (pinned)": a density-matrix path (16x16) for
the incoherent channels; coherent always-on ZZ per NN bond; watcher-angle offset
2*eps_cal; per-qubit RZ statics + interleave-slot drift; T1/T2 idle + per-gate
depolarizing with the lambda-map; asymmetric per-qubit confusion + CAL
finite-shot estimation -> tensor mitigation inversion (quasi-probs kept) ->
conditioning, in the pinned order; the S-correlated readout-bias injection in
BOTH modes (r-independent asymmetric additive; arm-dependent state-dependent);
the SIGNAL / NULL / three-impostor generative models. Bands, controls, and the
joint power (stage 4) build on these models; the self-tests at the bottom must
stay green through every stage (stage-1 is branch 2 of the positive control).

Qubit order (little-endian bit positions in outcome index): k=0, j=1, S=2, jp=3.
Line k - j - S - j': bonds (k,j), (j,S), (S,j').
"""

from dataclasses import dataclass, field

import numpy as np

# ---- circuit constants (pinned by the pre-reg) ----
R_ARMS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]   # 9 arms
INTERLEAVE = [1.0, 0.5, 1.5, 1.75, 0.0, 0.75, 2.0, 1.25, 0.25]  # pinned order, r=0 slot 5, r=2 slot 7
# (round 30, physics m4: slots of r=1.75 and r=0.25 SWAPPED vs the round-18
# order -- the old order put r=1.75 at MAX drift distance d=4 from the r=0
# anchor, which EASES the clause-(b) raw ordering by ~+0.012 at budget drift
# and suppresses the fallback center ~2%; now r=1.75 sits at d=1, and the
# max-distance slots d=4 carry r=1.0 and r=0.25)
Q_K, Q_J, Q_S, Q_JP = 0, 1, 2, 3
N_QUBITS = 4
DIM = 1 << N_QUBITS

# ---- exact gate layer (statevector, 16-dim) ----

def _hadamard_all(psi):
    h = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0)
    psi = psi.reshape([2] * N_QUBITS)
    for q in range(N_QUBITS):
        psi = np.tensordot(h, psi, axes=([1], [N_QUBITS - 1 - q]))
        psi = np.moveaxis(psi, 0, N_QUBITS - 1 - q)
    return psi.reshape(DIM)


def _single(psi, q, u):
    psi = psi.reshape([2] * N_QUBITS)
    psi = np.tensordot(u, psi, axes=([1], [N_QUBITS - 1 - q]))
    psi = np.moveaxis(psi, 0, N_QUBITS - 1 - q)
    return psi.reshape(DIM)


def rzz(psi, qa, qb, theta):
    """exp(-i (theta/2) Z_qa Z_qb), diagonal phases on the statevector."""
    idx = np.arange(DIM)
    za = 1 - 2 * ((idx >> qa) & 1)
    zb = 1 - 2 * ((idx >> qb) & 1)
    return psi * np.exp(-1j * (theta / 2.0) * za * zb)


def science_state(r, watcher_angle_error=0.0, write_angle_sj=np.pi / 2, write_angle_sjp=np.pi / 2):
    """The pre-measurement state of arm r (before any pre-rotation).

    H on all four; RZZ(pi/2) on (S,j) the write; RZZ(pi/2) on (S,j') the control
    write; watcher = TWO half-angle RZZ(r*pi/4) gates on (j,k) for r > 0 (the r=0
    arm omits the watcher entirely: 2 gates vs 4, the structurally privileged
    reference). watcher_angle_error = coherent per-half-gate angle offset (eps_cal
    enters twice: eps_watch = 2*eps_cal), for the injection stages.
    """
    psi = np.zeros(DIM, dtype=complex)
    psi[0] = 1.0
    psi = _hadamard_all(psi)
    psi = rzz(psi, Q_S, Q_J, write_angle_sj)
    psi = rzz(psi, Q_S, Q_JP, write_angle_sjp)
    if r > 0:
        half = r * np.pi / 4.0 + watcher_angle_error
        psi = rzz(psi, Q_J, Q_K, half)
        psi = rzz(psi, Q_J, Q_K, half)
    return psi


# ---- measurement layer ----

_H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0)
_SDG = np.array([[1, 0], [0, -1j]], dtype=complex)


def measurement_distribution(psi, variant):
    """Exact outcome distribution over the 16 basis states after the pinned
    pre-rotations: S in Z (none), k in Z (none), witnesses j and j' pre-rotated
    per variant: A = H (measures X), B = Sdg THEN H (measures Y; order
    load-bearing, the sim gate asserts Y-recovery below)."""
    for q in (Q_J, Q_JP):
        if variant == "B":
            psi = _single(psi, q, _SDG)
        psi = _single(psi, q, _H)
    return np.abs(psi) ** 2


def arm_distributions(r, **kw):
    """Both variants' exact distributions for arm r."""
    psi = science_state(r, **kw)
    return {"A": measurement_distribution(psi.copy(), "A"),
            "B": measurement_distribution(psi.copy(), "B")}


# ---- estimator chain (distribution level; the sampling layer wraps this) ----

def conditional_bloch(dists, witness):
    """b_s(w) = (<X_w>_s, <Y_w>_s) for s in {0,1}: X from variant A, Y from
    variant B, each conditioned on the S outcome of ITS OWN distribution (the
    branches are conditioned identically in every model; post-selection novelty
    owned, both branches kept)."""
    idx = np.arange(DIM)
    s_bit = (idx >> Q_S) & 1
    w_sign = 1 - 2 * ((idx >> witness) & 1)
    out = {}
    for s in (0, 1):
        comps = []
        for var in ("A", "B"):
            p = dists[var]
            sel = p[s_bit == s]
            signs = w_sign[s_bit == s]
            norm = sel.sum()
            comps.append(float((sel * signs).sum() / norm) if norm > 0 else np.nan)
        out[s] = np.array(comps)  # (X, Y)
    return out


def record_vector(dists, witness):
    """b0(w) - b1(w): the unnormalized record 2-vector of witness w."""
    b = conditional_bloch(dists, witness)
    return b[0] - b[1]


def fit_axes(dists_r0):
    """Per-witness record axes u_w = unit of b0(w) - b1(w) at the r = 0 arm.
    (The split-sample even/odd partition lives in the SAMPLING layer; at
    distribution level fit and evaluation coincide, which is exact for the
    ideal check and irrelevant to it.)"""
    axes = {}
    for name, w in (("j", Q_J), ("jp", Q_JP)):
        d = record_vector(dists_r0, w)
        n = np.linalg.norm(d)
        axes[name] = d / n if n > 0 else np.array([np.nan, np.nan])
    return axes


def signed_projection(dists, witness, axis):
    """S_hat(w; r) = 1/2 * (b0(w) - b1(w)) . u_w  (SIGNED; the folded D_hat = |S_hat|
    appears only where a magnitude is required: denominators of the inner ratio)."""
    return 0.5 * float(record_vector(dists, witness) @ axis)


def transverse_projection(dists, witness, axis):
    """T_hat(w; r) = the same construction on the perpendicular axis (guard)."""
    perp = np.array([-axis[1], axis[0]])
    return 0.5 * float(record_vector(dists, witness) @ perp)


def double_ratio(arms, axes):
    """The verdict carrier rho_hat(r) = [S(j;r)/D(j';r)] / [S(j;0)/D(j';0)] for
    every arm; arms = {r: dists}. Returns {r: rho_hat}."""
    s_j0 = signed_projection(arms[0.0], Q_J, axes["j"])
    d_jp0 = abs(signed_projection(arms[0.0], Q_JP, axes["jp"]))
    ref = s_j0 / d_jp0
    out = {}
    for r, dists in arms.items():
        s_j = signed_projection(dists, Q_J, axes["j"])
        d_jp = abs(signed_projection(dists, Q_JP, axes["jp"]))
        out[r] = (s_j / d_jp) / ref
    return out


# ---- stage 2: the sampling layer (shots, and the pinned even/odd split-sample) ----
# Pinned by the pre-reg: 16384 shots/circuit, the two r = 0 science arms at 2x = 32768;
# the split-sample partition lives ONLY at the two r = 0 arms (round-20 scope): EVEN shot
# indices in delivered memory order fit the axis, ODD indices evaluate the normalizer;
# each half is aggregated to counts BEFORE mitigation/conditioning (the inversion operates
# on distributions, so mitigate-then-split is impossible). All other arms use full shots.

SHOTS = 16384
SHOTS_R0 = 32768


def sample_shot_sequence(dist, shots, rng):
    """A delivered per-shot memory sequence (outcome indices), the sim analog of
    memory=True. The even/odd split is defined on THIS order. Only TRUE physical
    distributions are ever sampled (mitigated quasi-probs stay downstream); the
    clip removes float dust of order 1e-17 from the channel algebra, nothing
    physical."""
    d = np.clip(dist, 0.0, None)
    return rng.choice(DIM, size=shots, p=d / d.sum())


def counts_to_dist(counts):
    total = counts.sum()
    return counts / total if total > 0 else np.full(DIM, np.nan)


def sampled_arm_data(r, rng, shots=None, dist_override=None):
    """Sample both variants of arm r. Returns {"A": seq, "B": seq}."""
    dists = dist_override if dist_override is not None else arm_distributions(r)
    n = shots if shots is not None else (SHOTS_R0 if r == 0.0 else SHOTS)
    return {v: sample_shot_sequence(dists[v], n, rng) for v in ("A", "B")}


def estimator_from_shots(shot_data, mitigate=None):
    """The full estimator on sampled per-shot data. The pinned pipeline order:
    split the raw memory even/odd AT THE TWO r = 0 ARMS -> aggregate each half
    to counts -> CAL0/CAL1 inversion per half (the `mitigate` hook, a dist->dist
    map; quasi-probabilities kept) -> conditioning per half. All other arms:
    full shots -> counts -> inversion -> conditioning.

    shot_data: {r: {"A": shot sequence, "B": shot sequence}} for all 9 arms.
    Returns dict with axes, S_hat/D_hat tables, rho_hat, sigma_hat, and the
    transverse guards T_hat."""
    mit = mitigate if mitigate is not None else (lambda d: d)
    # r = 0: split even/odd per variant, aggregate each half to counts, mitigate per half
    fit_dists, eval_dists = {}, {}
    for v in ("A", "B"):
        seq = shot_data[0.0][v]
        fit_dists[v] = mit(counts_to_dist(np.bincount(seq[0::2], minlength=DIM)))
        eval_dists[v] = mit(counts_to_dist(np.bincount(seq[1::2], minlength=DIM)))
    axes = fit_axes(fit_dists)                       # axis from the EVEN half
    arm_dists = {r: {v: mit(counts_to_dist(np.bincount(shot_data[r][v], minlength=DIM)))
                     for v in ("A", "B")} for r in shot_data if r != 0.0}
    arm_dists[0.0] = eval_dists                      # normalizer from the ODD half

    s_j0 = signed_projection(arm_dists[0.0], Q_J, axes["j"])
    d_jp0 = abs(signed_projection(arm_dists[0.0], Q_JP, axes["jp"]))
    ref = s_j0 / d_jp0
    rho, t_guards = {}, {}
    for r, dists in arm_dists.items():
        s_j = signed_projection(dists, Q_J, axes["j"])
        d_jp = abs(signed_projection(dists, Q_JP, axes["jp"]))
        rho[r] = (s_j / d_jp) / ref
        t_guards[r] = {w: transverse_projection(dists, q, axes[w])
                       for w, q in (("j", Q_J), ("jp", Q_JP))}
    sigma_hat = {r: rho[r] / abs(rho[2.0]) for r in rho}
    return {"axes": axes, "rho": rho, "sigma": sigma_hat, "t": t_guards,
            "s_j0": s_j0, "d_jp0": d_jp0}


# ---- stage 3: the injection bank + generative models (density-matrix path) ----
# The incoherent channels need a density matrix; 16x16 is trivial. The statevector
# path above stays the exact reference at zero incoherent noise (asserted below).

_I2 = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULI1 = [_I2, _X, _Y, _Z]

BONDS = [(Q_K, Q_J), (Q_J, Q_S), (Q_S, Q_JP)]          # the NN bonds of the line
SLOT = {r: i + 1 for i, r in enumerate(INTERLEAVE)}    # 1-based interleave slot per arm

MODEL_KINDS = ("signal", "null", "imp1", "imp2", "imp3", "impz")
# impz (round 30, from the round-29 physics M3 find): the Z-POLARISED
# never-watching watcher -- H on k omitted, everything else identical; it
# reproduces the ENTIRE signed rho_hat curve exactly (j's coherence ROTATES by
# r*pi/2 instead of attenuating), and ONLY the |<Z_k>| guard separates it: the
# 4th negative-control kind, guard-separated (the LR cannot see it), with its
# own detection-power measurement.

# The z_k DISCRIMINATION CEILING (round 30, stats B1: impz's sole discriminator
# z_k is a TO-FREEZE band the round-28 coverage rule makes AS LARGE AS
# POSSIBLE -- the 8th protection-interaction instance, the first BETWEEN two
# consecutive amendments; the escape window m < band must be CAPPED by a pinned
# maximum, the opposite of the coverage rule). The frozen guard_z_k must
# satisfy envelope <= ceiling; an envelope above the ceiling is UNFREEZABLE
# (a design remedy, never a loosened ceiling). Rehearsal envelopes read
# 0.038-0.041; the boundary control runs impz at m = the frozen band value.
Z_K_BAND_CEILING = 0.05


def _u1(q, m):
    """16x16 embedding of a single-qubit operator m on (little-endian) qubit q."""
    out = np.array([[1.0 + 0j]])
    for qq in range(N_QUBITS - 1, -1, -1):
        out = np.kron(out, m if qq == q else _I2)
    return out


def _rho0():
    rho = np.zeros((DIM, DIM), dtype=complex)
    rho[0, 0] = 1.0
    return rho


def _apply_u(rho, u):
    return u @ rho @ u.conj().T


def _apply_1q(rho, q, m):
    """m rho m^dag for a single-qubit operator m on qubit q, via tensordot on
    the reshaped rho (no 16x16 embedding built; the seed-bank hot path)."""
    r = rho.reshape([2] * (2 * N_QUBITS))
    ax_k = N_QUBITS - 1 - q
    ax_b = 2 * N_QUBITS - 1 - q
    r = np.tensordot(m, r, axes=([1], [ax_k]))
    r = np.moveaxis(r, 0, ax_k)
    r = np.tensordot(np.conj(m), r, axes=([1], [ax_b]))
    r = np.moveaxis(r, 0, ax_b)
    return r.reshape(DIM, DIM)


def _rzz_diag(qa, qb, theta):
    idx = np.arange(DIM)
    za = 1 - 2 * ((idx >> qa) & 1)
    zb = 1 - 2 * ((idx >> qb) & 1)
    return np.exp(-1j * (theta / 2.0) * za * zb)


def _rz_diag(q, phi):
    idx = np.arange(DIM)
    z = 1 - 2 * ((idx >> q) & 1)
    return np.exp(-1j * (phi / 2.0) * z)


def _apply_diag(rho, d):
    return rho * np.outer(d, d.conj())


def _apply_kraus1(rho, q, kraus):
    out = np.zeros_like(rho)
    for k in kraus:
        out += _apply_1q(rho, q, k)
    return out


_PAIR_PAULI_CACHE = {}


def _pair_paulis(qa, qb):
    key = (qa, qb)
    if key not in _PAIR_PAULI_CACHE:
        _PAIR_PAULI_CACHE[key] = [
            _u1(qa, _PAULI1[a]) @ _u1(qb, _PAULI1[b])
            for a in range(4) for b in range(4) if not (a == 0 and b == 0)]
    return _PAIR_PAULI_CACHE[key]


def _apply_depol2(rho, qa, qb, p_total):
    """Two-qubit depolarizing: (1-p) rho + (p/15) sum_{P != II} P rho P.
    EPC (average gate error) maps to p_total = 1.25 * EPC for d = 4 (NOMINAL map;
    the pre-reg tags the EPC-to-attenuation map channel-dependent, Aer governs)."""
    if p_total == 0.0:
        return rho
    acc = np.zeros_like(rho)
    for pm in _pair_paulis(qa, qb):
        acc += pm @ rho @ pm                # 2-qubit Paulis are Hermitian
    return (1.0 - p_total) * rho + (p_total / 15.0) * acc


def _apply_idle(rho, q, t, t1, t2, lam):
    """Amplitude damping (T1) + extra pure dephasing to total T2 over idle time t,
    all incoherent rates scaled by the lambda-map (lam = 0 is the ideal end)."""
    if t <= 0.0 or lam == 0.0:
        return rho
    g = 1.0 - np.exp(-lam * t / t1)
    k0 = np.array([[1, 0], [0, np.sqrt(1.0 - g)]], dtype=complex)
    k1 = np.array([[0, np.sqrt(g)], [0, 0]], dtype=complex)
    rho = _apply_kraus1(rho, q, [k0, k1])
    rphi = max(1.0 / t2 - 1.0 / (2.0 * t1), 0.0)     # pure-dephasing rate from T2, T1
    p_phi = (1.0 - np.exp(-lam * t * rphi)) / 2.0
    if p_phi > 0.0:
        zm = _u1(q, _Z)
        rho = (1.0 - p_phi) * rho + p_phi * (zm @ rho @ zm)
    return rho


@dataclass
class Params:
    """The injection bank. All-zero defaults = the ideal device; worst_admitted()
    carries the pre-reg's pinned day-of abort ceilings as budget values."""
    # coherent injections
    eps_cal: float = 0.0                      # per half-gate watcher angle offset (rad)
    zz_static: dict = field(default_factory=dict)   # {(qa, qb): theta_static rad} always-on ZZ per bond
    rz_static: dict = field(default_factory=dict)   # {q: phi rad} static frame offsets
    rz_drift: dict = field(default_factory=dict)    # {q: rad per interleave slot} arm-time drift
    write_miss_sj: float = 0.0                # static write-angle miss on (S, j)
    write_miss_sjp: float = 0.0               # static write-angle miss on (S, j')
    # incoherent injections (lambda-scaled)
    epc_write: float = 0.0                    # per write gate (EPC)
    epc_watch: float = 0.0                    # per watcher half-gate (EPC)
    epc_write_sj: float = None                # per-edge overrides (the j-vs-j' asymmetry
    epc_write_sjp: float = None               # draw writes them; None = epc_write)
    t1: dict = field(default_factory=dict)    # {q: seconds}
    t2: dict = field(default_factory=dict)    # {q: seconds}
    idle_watcher: float = 0.0                 # extra idle on S and j' during the watcher layer (s)
    lam: float = 1.0                          # the lambda-map scale on ALL incoherent rates
    eta_slope: float = 0.0                    # the +/-s_res SHAPE axis: a linear-in-r extra
                                              # dephasing on witness j (s > 0 loads high-r arms,
                                              # s < 0 loads low-r arms; relative profile residual
                                              # ~ s*(1 - r/2) vs the r = 2 anchor)
    # readout layer
    readout: dict = field(default_factory=dict)     # {q: (e01, e10)} true asymmetric confusion
    bias_delta: float = 0.0                   # S-correlated witness assignment differential (SIGNED; the
                                              # r-independent asymmetric additive mode, 2% ceiling)
    bias_sym: float = 0.0                     # symmetric S-independent witness flip (cancels; control)
    arm_bias_delta: float = 0.0               # the ARM-DEPENDENT mode: an assignment offset whose
                                              # coefficient scales with the arm's per-branch witness
                                              # population imbalance (T1-during-readout style: the 96/4
                                              # arm feels it, the 50/50 arm does not); science-only,
                                              # structurally invisible to the tensor CAL
    cal_shots: int = SHOTS                    # CAL0/CAL1 finite shots
    impz_m: float = 1.0                       # impz's k polarisation <Z_k> = m (1 = the pure
                                              # never-watching case; the boundary control runs
                                              # at m = the frozen z_k band value, round 30)

    @classmethod
    def worst_admitted(cls):
        """The worst-admitted basis from the pinned Class-1 abort ceilings:
        watcher/write EPC 0.5%, T2* floor 70 us (T1 at 2x the floor: the
        dephasing-dominated edge), 140 ns idle per watcher half-gate, S readout
        2% both directions / witnesses+k 3%, coherent budgets eps_cal = 0.01 rad
        and eps_ZZ = 0.01 rad on the watcher bond, drift budget 0.05 rad/arm."""
        return cls(
            eps_cal=0.01,
            zz_static={(Q_K, Q_J): 0.01, (Q_J, Q_S): 0.01, (Q_S, Q_JP): 0.01},
            rz_drift={Q_J: 0.05, Q_JP: 0.05},
            write_miss_sj=0.01, write_miss_sjp=0.01,   # the V_S guard's budgeted delta = eps_cal
            epc_write=0.005, epc_watch=0.005,
            t1={q: 140e-6 for q in range(N_QUBITS)},
            t2={q: 70e-6 for q in range(N_QUBITS)},
            idle_watcher=280e-9,
            readout={Q_S: (0.015, 0.02), Q_J: (0.02, 0.03), Q_JP: (0.02, 0.03), Q_K: (0.02, 0.03)},
            bias_delta=0.02,
            arm_bias_delta=0.02,
        )


def _science_rho_once(r, p, write_sj, halves, h_qubits=None, k_flip=False):
    """One density-matrix run of the science circuit: H layer, writes (+depol),
    watcher halves (+depol, +idle on S/j'), static ZZ, RZ frames. `halves` is the
    list of signed watcher half-gate angles (empty for the r = 0 arm).
    h_qubits: the qubits receiving the initial H (default all; the impz model
    omits k's). k_flip prepares k in |1> instead of |0> (the impz partial-
    polarisation mixture's second branch; the flip sits at PREPARATION, before
    the conditional phases -- an X after the circuit would be a different
    channel)."""
    rho = _rho0()
    if k_flip:
        rho = _apply_1q(rho, Q_K, _X)
    for q in (range(N_QUBITS) if h_qubits is None else h_qubits):
        rho = _apply_1q(rho, q, _H)
    epc_sj = p.epc_write_sj if p.epc_write_sj is not None else p.epc_write
    epc_sjp = p.epc_write_sjp if p.epc_write_sjp is not None else p.epc_write
    rho = _apply_diag(rho, _rzz_diag(Q_S, Q_J, write_sj))
    rho = _apply_depol2(rho, Q_S, Q_J, 1.25 * epc_sj * p.lam)
    rho = _apply_diag(rho, _rzz_diag(Q_S, Q_JP, np.pi / 2 + p.write_miss_sjp))
    rho = _apply_depol2(rho, Q_S, Q_JP, 1.25 * epc_sjp * p.lam)
    for ang in halves:
        rho = _apply_diag(rho, _rzz_diag(Q_J, Q_K, ang))
        rho = _apply_depol2(rho, Q_J, Q_K, 1.25 * p.epc_watch * p.lam)
    if halves and p.idle_watcher > 0.0:
        for q in (Q_S, Q_JP):
            rho = _apply_idle(rho, q, p.idle_watcher,
                              p.t1.get(q, np.inf), p.t2.get(q, np.inf), p.lam)
    if p.eta_slope != 0.0:
        s = p.eta_slope
        g = 1.0 - abs(s) * (r / 2.0 if s > 0 else (2.0 - r) / 2.0)
        p_phi = (1.0 - g) / 2.0
        if p_phi > 0.0:
            rho = (1.0 - p_phi) * rho + p_phi * _apply_1q(rho, Q_J, _Z)
    for (qa, qb), th in p.zz_static.items():
        rho = _apply_diag(rho, _rzz_diag(qa, qb, th))
    slot = SLOT[r]
    for q in range(N_QUBITS):
        phi = p.rz_static.get(q, 0.0) + p.rz_drift.get(q, 0.0) * slot
        if phi != 0.0:
            rho = _apply_diag(rho, _rz_diag(q, phi))
    return rho


def science_rho(r, p, kind="signal"):
    """The pre-measurement density matrix of arm r under a generative model:
    signal = law-holds; null = the (S, j) write's conditional phase zeroed, all
    else identical; imp1 = monotone erasure (per-half-gate independent sign
    randomization -> 4-way average -> multiplier cos^2(r*pi/4)); imp2 = the
    half-period curve cos(r*pi/4) (half watcher angle); imp3 = sign-flat
    |cos(r*pi/2)| (arms r > 1 realized at the reflected angle 2 - r)."""
    write_sj = 0.0 if kind == "null" else np.pi / 2 + p.write_miss_sj
    h_qubits = (Q_J, Q_S, Q_JP) if kind == "impz" else None

    def _mix(*args):
        # impz at partial polarisation <Z_k> = m: the PREPARATION mixture of
        # k = |0> (weight (1+m)/2) and k = |1> (weight (1-m)/2) runs
        if kind == "impz" and p.impz_m < 1.0:
            w1 = (1.0 - p.impz_m) / 2.0
            return ((1.0 - w1) * _science_rho_once(*args, k_flip=False)
                    + w1 * _science_rho_once(*args, k_flip=True))
        return _science_rho_once(*args)

    if r <= 0:
        return _mix(r, p, write_sj, [], h_qubits)
    if kind == "imp2":
        half = r * np.pi / 8.0 + p.eps_cal
    elif kind == "imp3":
        r_eff = 2.0 - r if r > 1.0 else r
        half = r_eff * np.pi / 4.0 + p.eps_cal
    else:
        half = r * np.pi / 4.0 + p.eps_cal
    if kind == "imp1":
        acc = np.zeros((DIM, DIM), dtype=complex)
        for s1 in (1, -1):
            for s2 in (1, -1):
                acc += _science_rho_once(r, p, write_sj, [s1 * half, s2 * half])
        return acc / 4.0
    return _mix(r, p, write_sj, [half, half], h_qubits)


def measurement_distribution_rho(rho, variant):
    """The measured-basis distribution from a density matrix (pre-rotations per
    the pinned variants, then the Z-basis diagonal)."""
    for q in (Q_J, Q_JP):
        if variant == "B":
            rho = _apply_1q(rho, q, _SDG)
        rho = _apply_1q(rho, q, _H)
    return np.real(np.diag(rho)).copy()


# ---- readout layer: confusion, CAL, mitigation, the bias injections ----

def _confusion_matrix(e01, e10):
    """C[measured, true]: column 0 = prep |0> reads 1 w.p. e01; column 1 = prep
    |1> reads 0 w.p. e10."""
    return np.array([[1.0 - e01, e10], [e01, 1.0 - e10]])


def apply_tensor_channel(dist, mats):
    """Apply per-qubit 2x2 stochastic (or inverse) matrices to a 16-dim
    distribution; mats = {q: 2x2}."""
    v = dist.reshape([2] * N_QUBITS)
    for q, m in mats.items():
        ax = N_QUBITS - 1 - q
        v = np.tensordot(m, v, axes=([1], [ax]))
        v = np.moveaxis(v, 0, ax)
    return v.reshape(DIM)


def _cond_flip(dist, w, f0, f1):
    """Flip witness bit w with probability f0 where the MEASURED S bit is 0 and
    f1 where it is 1. f0 = f1 is the SYMMETRIC MULTIPLICATIVE mode (shrinks
    <sign_w> per branch; the mode the double ratio cancels)."""
    idx = np.arange(DIM)
    f = np.where(((idx >> Q_S) & 1) == 0, f0, f1)
    out = dist * (1.0 - f)
    np.add.at(out, idx ^ (1 << w), dist * f)
    return out


def _cond_assign_offset(dist, w, delta):
    """The r-independent ASYMMETRIC ADDITIVE mode: a fixed S-ANTISYMMETRIC
    assignment offset on witness w: the conditional probability p(w=1 | S=s)
    shifts by +delta/2 at s = 0 and -delta/2 at s = 1, i.e. the conditional
    Bloch component b_s gets a CONSTANT -/+ delta offset on the measured axis
    while the true record varies across arms (the round-9 mechanism; the
    round-27 linearization |rho_hat(1.75)| shift ~ -/+0.043 at delta = +/-2%,
    eta = 0.92, is this channel's exact arithmetic). Realized as proportional
    mass transfer between the w-sectors inside each measured-S branch.

    BASIS PINNED (round 29, physics M1: 'a percent of WHAT' was undefined and
    two doc figures sat on bases 2x apart): delta IS the S-branch assignment
    differential Delta = [p(w=1|S=0) - p(w=1|S=1)]-shift, which equals the
    induced record offset on S_hat; the 2% ceiling, this injection, and the
    in-job crosstalk pair's abort compare are all in THIS quantity. This is
    the WIDER (conservative) of the two readings the review found."""
    out = dist.copy()
    idx = np.arange(DIM)
    sbit = (idx >> Q_S) & 1
    wbit = (idx >> w) & 1
    for s, dlt in ((0, +delta / 2.0), (1, -delta / 2.0)):
        bsel = sbit == s
        p_branch = dist[bsel].sum()
        if p_branch <= 0 or dlt == 0.0:
            continue
        src = 0 if dlt > 0 else 1            # mass moves src -> flipped partner
        m_src = bsel & (wbit == src)
        p_src = dist[m_src].sum()
        move = abs(dlt) * p_branch
        if p_src <= 0:
            continue
        frac = min(move / p_src, 1.0)        # saturate rather than overdraw
        out[m_src] -= dist[m_src] * frac
        np.add.at(out, idx[m_src] ^ (1 << w), dist[m_src] * frac)
    return out


def _branch_imbalance(dist, w):
    """|p(w=1 | S=0) - p(w=1 | S=1)| of a true distribution: the per-branch
    witness population imbalance that makes state-dependent readout error
    arm-dependent (~0 at r = 1, ~1 at the fixpoints)."""
    idx = np.arange(DIM)
    sbit = (idx >> Q_S) & 1
    wbit = (idx >> w) & 1
    p1 = []
    for s in (0, 1):
        b = dist[sbit == s].sum()
        p1.append(dist[(sbit == s) & (wbit == 1)].sum() / b if b > 0 else 0.0)
    return abs(p1[0] - p1[1])


def readout_channel(dist, p, arm_imb=None):
    """The full readout layer on one true distribution: per-qubit asymmetric
    confusion, then the S-correlated bias modes on both witnesses (a property
    of the physical Z-readout: applied identically in variants A and B, per
    round-25). Mode 1 = the r-independent S-antisymmetric additive offset at
    the ceiling; mode 2 = the ARM-DEPENDENT offset, coefficient scaled by the
    arm's per-branch population imbalance (`arm_imb` = {witness: imbalance},
    computed from the arm's TRUE distribution). A witness-local arm-independent
    affine readout map is exactly invisible to the double ratio (verified in
    the stage-3 tests); the dangerous modes are exactly these two."""
    mats = {q: _confusion_matrix(*p.readout.get(q, (0.0, 0.0))) for q in range(N_QUBITS)}
    out = apply_tensor_channel(dist, mats)
    for w in (Q_J, Q_JP):
        if p.bias_sym:
            out = _cond_flip(out, w, p.bias_sym, p.bias_sym)
        if p.bias_delta:
            out = _cond_assign_offset(out, w, p.bias_delta)
        if p.arm_bias_delta and arm_imb is not None:
            out = _cond_assign_offset(out, w, p.arm_bias_delta * arm_imb[w])
    return out


def model_arm_dists(kind, r, p):
    """The post-readout outcome distributions of arm r under a generative model:
    the sampling layer's dist_override. Returns {"A": dist, "B": dist}."""
    rho = science_rho(r, p, kind)
    out = {}
    for v in ("A", "B"):
        true = measurement_distribution_rho(rho, v)
        imb = {w: _branch_imbalance(true, w) for w in (Q_J, Q_JP)}
        out[v] = readout_channel(true, p, imb)
    return out


def cal_true_dists(p):
    """CAL0/CAL1 true outcome distributions: bare preps through the TRUE
    per-qubit confusion only; the tensor CAL structurally cannot see the
    S-correlated bias or the science-only state-dependent extra (that blindness
    is the injected danger, not an omission)."""
    mats = {q: _confusion_matrix(*p.readout.get(q, (0.0, 0.0))) for q in range(N_QUBITS)}
    d0 = np.zeros(DIM); d0[0] = 1.0
    d1 = np.zeros(DIM); d1[DIM - 1] = 1.0
    return apply_tensor_channel(d0, mats), apply_tensor_channel(d1, mats)


def estimate_confusion(seq0, seq1, shots0, shots1):
    """Per-qubit confusion estimates from sampled CAL0/CAL1 shot sequences."""
    mats = {}
    for q in range(N_QUBITS):
        e01_hat = float(np.mean((seq0 >> q) & 1)) if shots0 else 0.0
        e10_hat = 1.0 - float(np.mean((seq1 >> q) & 1)) if shots1 else 0.0
        mats[q] = _confusion_matrix(e01_hat, e10_hat)
    return mats


def make_mitigator(mats_hat):
    """The CAL0/CAL1 tensor-product inversion (quasi-probabilities kept, no
    clamping) as a dist -> dist map for the estimator hook."""
    invs = {q: np.linalg.inv(m) for q, m in mats_hat.items()}
    return lambda dist: apply_tensor_channel(dist, invs)


def sampled_cal_mitigator(p, rng):
    """Sample the CAL pair at finite shots and build the mitigator (the pinned
    CAL finite-shot noise inclusion)."""
    c0, c1 = cal_true_dists(p)
    seq0 = rng.choice(DIM, size=p.cal_shots, p=c0 / c0.sum())
    seq1 = rng.choice(DIM, size=p.cal_shots, p=c1 / c1.sum())
    return make_mitigator(estimate_confusion(seq0, seq1, p.cal_shots, p.cal_shots))


def rho_signed_model(p, r):
    """The SIGNED distribution-level rho_hat(r) of the INCOHERENT-ONLY config
    (coherent injections OFF, the round-23 no-double-count guard). This is
    cos(r*pi/2) * eta(r), NOT eta itself."""
    q = Params(epc_write=p.epc_write, epc_watch=p.epc_watch, t1=dict(p.t1),
               t2=dict(p.t2), idle_watcher=p.idle_watcher, lam=p.lam,
               eta_slope=p.eta_slope)
    arms = {rr: model_arm_dists("signal", rr, q) for rr in (0.0, r)}
    axes = fit_axes(arms[0.0])
    s_j0 = signed_projection(arms[0.0], Q_J, axes["j"])
    d_jp0 = abs(signed_projection(arms[0.0], Q_JP, axes["jp"]))
    s_j = signed_projection(arms[r], Q_J, axes["j"])
    d_jp = abs(signed_projection(arms[r], Q_JP, axes["jp"]))
    return (s_j / d_jp) / (s_j0 / d_jp0)


def eta_of(p, r=2.0):
    """The incoherent-only watcher ATTENUATION at r = 2, where |cos| = 1 and
    |rho_hat(2)| = eta(2) exactly (the lambda-solve's target read). At generic
    r this is |cos * eta|, NOT eta: use eta_attenuation_profile for profiles."""
    return abs(rho_signed_model(p, r))


def eta_attenuation_profile(p_mid):
    """eta_nom(r): the per-arm watcher ATTENUATION profile of the incoherent
    config (|rho_model(r) / cos(r*pi/2)|; r = 1 interpolated from neighbors,
    where cos = 0 makes the division 0/0 and the (c)-target term cos * eta is
    0 regardless). The stage-3 build's first freeze attempt stored |cos * eta|
    here and bent the clause-(c) target into a cos*|cos| curve, collapsing
    corner power to 0.10: eta_nom is attenuation ONLY."""
    prof = {0.0: 1.0}
    for r in R_ARMS:
        c = np.cos(r * np.pi / 2.0)
        if r > 0 and abs(c) > 1e-9:
            prof[r] = abs(rho_signed_model(p_mid, r) / c)
    prof[1.0] = 0.5 * (prof[0.75] + prof[1.25])
    return prof


def solve_lambda(p, target_eta2, tol=1e-9):
    """The lambda-map: ONE loss factor scaling all incoherent rates, solved so
    the incoherent-only eta(2) hits the target (lam = 0 is the ideal end); the
    same map realizes every level point of the sweep."""
    q = Params(epc_write=p.epc_write, epc_watch=p.epc_watch, t1=dict(p.t1),
               t2=dict(p.t2), idle_watcher=p.idle_watcher)
    lo, hi = 0.0, 1.0
    q.lam = hi
    while eta_of(q) > target_eta2:            # widen until the target is bracketed
        hi *= 2.0
        q.lam = hi
        if hi > 1e4:
            raise RuntimeError("lambda solve: target below reachable eta range")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        q.lam = mid
        if eta_of(q) > target_eta2:
            lo = mid
        else:
            hi = mid
    q.lam = 0.5 * (lo + hi)
    return q.lam


# ---- stage 4a: guard statistics, extra PUB models, and the frozen verdict code ----
# The verdict machinery per the pre-reg's "Verdict rule" section: the 8-test
# pooled guard bank (5 pooled + 3 floors), clauses (a) -> (b) -> (c) with
# short-circuit, the 8-arm Gaussian LR against the three impostors, precedence
# VOID -> CONFIRMED -> DEVICE-DEVIATION -> INCONCLUSIVE. Band VALUES are inputs
# (frozen at stage 4b through THIS code); the structure is pinned here.

R_INTERIOR = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]   # clause (c)'s seven arms
R_LR = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]    # the 8 non-degenerate LR arms

IMPOSTOR_CURVES = {
    "imp1": lambda r: np.cos(r * np.pi / 4.0) ** 2,
    "imp2": lambda r: np.cos(r * np.pi / 4.0),
    "imp3": lambda r: abs(np.cos(r * np.pi / 2.0)),
}


def zleak_model_dists(kind, p):
    """The Z-leak pair: the r = 0 and r = 2 science circuits with NO
    pre-rotation on any qubit, all four measured raw in Z."""
    out = {}
    for r in (0.0, 2.0):
        rho = science_rho(r, p, kind)
        true = np.real(np.diag(rho)).copy()
        imb = {w: _branch_imbalance(true, w) for w in (Q_J, Q_JP)}
        out[r] = readout_channel(true, p, imb)
    return out


def vs_guard_model_dists(kind, p):
    """The V_S guard PUBs: S pre-rotated to X ('X') or Y ('Y'), witnesses and k
    measured raw in Z, at r in {0, 2} (4 PUBs)."""
    out = {}
    for r in (0.0, 2.0):
        rho = science_rho(r, p, kind)
        for basis, pre in (("X", (_H,)), ("Y", (_SDG, _H))):
            rr = rho
            for u in pre:
                rr = _apply_1q(rr, Q_S, u)
            true = np.real(np.diag(rr)).copy()
            imb = {w: _branch_imbalance(true, w) for w in (Q_J, Q_JP)}
            out[(r, basis)] = readout_channel(true, p, imb)
    return out


def tomograph_model_dists(kind, p):
    """The WATCHER-TOMOGRAPH PUB pair (round 30, physics BLOCKER): r in {0, 1}
    with k pre-rotated Sdg THEN H (reads Y_k) and j, S, j' measured raw in Z.
    The separating correlator is <Y_k Z_j>: +1 at r = 1 for a SUPERPOSED
    watcher (the round-10 grounding: the r = 1 watcher record is the perfect
    classical bit <Y_k Z_j> = 1), identically 0 for ANY Z-diagonal k (the impz
    class at every m -- Y of a Z-diagonal state is 0). Every other flown PUB is
    Z-diagonal on k and sees only diag(rho_k), where |+><+| and I/2 coincide:
    this pair is the ONLY place the design measures that the watcher actually
    watched. DIAGNOSTIC-ONLY, like V_S: it gates the mechanism sentence
    (watcher superposition certified), never the record verdict."""
    out = {}
    for r in (0.0, 1.0):
        rho = science_rho(r, p, kind)
        for u in (_SDG, _H):
            rho = _apply_1q(rho, Q_K, u)
        true = np.real(np.diag(rho)).copy()
        imb = {w: _branch_imbalance(true, w) for w in (Q_J, Q_JP)}
        out[r] = readout_channel(true, p, imb)
    return out


def tomograph_stat(tomo_dists):
    """<Y_k Z_j> at the r = 1 tomograph PUB (after the k pre-rotation, Y_k is
    the measured Z of k): the mechanism statistic, prediction +1 (dressed)."""
    idx = np.arange(DIM)
    yk = 1 - 2 * ((idx >> Q_K) & 1)
    zj = 1 - 2 * ((idx >> Q_J) & 1)
    d = tomo_dists[1.0]
    tot = d.sum()
    return float((d * yk * zj).sum() / tot) if tot != 0 else np.nan


def _mean_sign(dist, q):
    idx = np.arange(DIM)
    sgn = 1 - 2 * ((idx >> q) & 1)
    tot = dist.sum()
    return float((dist * sgn).sum() / tot) if tot != 0 else np.nan


def guard_stats(est, arm_dists_all, zleak_dists):
    """The 5 pooled guard statistics + the 3 floors, from the mitigated
    distributions. arm_dists_all: {r: {"A": dist, "B": dist}} (the same
    mitigated dists the estimator consumed; r = 0 is the eval half).
    zleak_dists: {0.0: dist, 2.0: dist} mitigated raw-Z distributions."""
    t_sup = {w: max(abs(est["t"][r][w]) for r in est["t"]) for w in ("j", "jp")}
    z_k = max(abs(_mean_sign(arm_dists_all[r][v], Q_K)) for r in arm_dists_all for v in ("A", "B"))
    z_s = max(abs(_mean_sign(arm_dists_all[r][v], Q_S)) for r in arm_dists_all for v in ("A", "B"))
    z_leak = max(abs(_mean_sign(zleak_dists[r], w)) for r in zleak_dists for w in (Q_J, Q_JP))
    d_jp_min = min(abs(signed_projection(arm_dists_all[r], Q_JP, est["axes"]["jp"]))
                   for r in arm_dists_all if r != 0.0)
    return {"t_sup": t_sup, "z_k": z_k, "z_s": z_s, "z_leak": z_leak,
            "s_j0": est["s_j0"], "d_jp0": est["d_jp0"], "d_jp_min": d_jp_min}


def vs_guard_stat(vs_dists):
    """|b(S | z_j)|: the j-conditional equatorial S-coherence, max over the j
    outcome, r, and pairing the X/Y bases into one 2-vector per (r, z_j)."""
    idx = np.arange(DIM)
    jbit = (idx >> Q_J) & 1
    s_sign = 1 - 2 * ((idx >> Q_S) & 1)
    worst = 0.0
    for r in (0.0, 2.0):
        for zj in (0, 1):
            comps = []
            for basis in ("X", "Y"):
                d = vs_dists[(r, basis)]
                sel = d[jbit == zj]
                sgn = s_sign[jbit == zj]
                tot = sel.sum()
                comps.append(float((sel * sgn).sum() / tot) if tot != 0 else np.nan)
            worst = max(worst, float(np.hypot(*comps)))
    return worst


@dataclass
class Bands:
    """The frozen constants the verdict code consumes. Stage 4b freezes them
    through the pinned band procedure; provisional values are for wiring tests
    only and never fly."""
    b_blind: float
    b_forgive: float
    b_curve: float
    se_contrast: float
    eta_nom: dict                      # {r: eta_nom(r)} mid-range basis targets
    lr_sigma: np.ndarray               # 8x8 empirical covariance of rho_hat over R_LR
    lr_thresholds: dict                # {impostor: threshold}
    floor_s_j0: float
    floor_d_jp0: float
    floor_d_jp_min: float
    guard_t_sup: float
    guard_z_k: float
    guard_z_s: float
    guard_z_leak: float
    b_sym: float = np.inf              # informational only
    fallback_175: tuple = None         # clause-(b) fallback band (lo, hi) if the branch engages
    use_fallback_175: bool = False
    vs_threshold: float = np.inf       # gates ONLY the MI corollary, never the verdict
    tomo_threshold: float = -np.inf    # watcher-tomograph lower cut: gates ONLY the
                                       # mechanism sentence (round 30), never the verdict


def lr_statistics(rho_vec, bands):
    """LR_i = 0.5[(x-mu_law)' Sinv (x-mu_law) - (x-mu_i)' Sinv (x-mu_i)] on the
    8-arm rho_hat vector, both hypotheses dressed with the SAME eta_nom."""
    x = np.array([rho_vec[r] for r in R_LR])
    eta = np.array([bands.eta_nom[r] for r in R_LR])
    mu_law = np.array([np.cos(r * np.pi / 2.0) for r in R_LR]) * eta
    sinv = np.linalg.inv(bands.lr_sigma)
    out = {}
    d_law = x - mu_law
    q_law = float(d_law @ sinv @ d_law)
    for name, curve in IMPOSTOR_CURVES.items():
        mu_i = np.array([curve(r) for r in R_LR]) * eta
        d_i = x - mu_i
        out[name] = 0.5 * (q_law - float(d_i @ sinv @ d_i))
    return out


def verdict(est, guards, vs_stat, bands, vs_threshold=np.inf, tomo_stat=None):
    """The frozen verdict rule. Returns (verdict, detail dict). Precedence:
    VOID -> CONFIRMED -> DEVICE-DEVIATION -> INCONCLUSIVE; clauses (a) -> (b)
    -> (c) short-circuit; the V_S guard gates ONLY the MI corollary and the
    watcher-tomograph gates ONLY the mechanism sentence (both reported, never
    verdict inputs)."""
    det = {"mi_corollary_claimable": bool(vs_stat <= vs_threshold), "vs_stat": vs_stat}
    if tomo_stat is not None:
        det["mechanism_claimable"] = bool(tomo_stat >= bands.tomo_threshold)
        det["tomo_stat"] = tomo_stat
    rho = est["rho"]
    # VOID: NaN / failed estimator
    vals = list(rho.values()) + [est["s_j0"], est["d_jp0"]]
    if any(not np.isfinite(v) for v in vals):
        det["void"] = ["nan"]
        return "VOID", det
    # VOID: the 8-test guard bank (5 pooled + 3 floors)
    trips = []
    if guards["t_sup"]["j"] > bands.guard_t_sup: trips.append("t_sup_j")
    if guards["t_sup"]["jp"] > bands.guard_t_sup: trips.append("t_sup_jp")
    if guards["z_k"] > bands.guard_z_k: trips.append("z_k")
    if guards["z_s"] > bands.guard_z_s: trips.append("z_s")
    if guards["z_leak"] > bands.guard_z_leak: trips.append("z_leak")
    if guards["s_j0"] < bands.floor_s_j0: trips.append("floor_s_j0")      # SIGNED floor > 0
    if guards["d_jp0"] < bands.floor_d_jp0: trips.append("floor_d_jp0")
    if guards["d_jp_min"] < bands.floor_d_jp_min: trips.append("floor_d_jp_min")
    if trips:
        det["void"] = trips
        return "VOID", det
    # clause (a): blindness
    a_ok = (abs(rho[1.0]) < bands.b_blind
            and min(abs(rho[0.75]), abs(rho[1.25])) > abs(rho[1.0]))
    det["a"] = a_ok
    if a_ok:
        # clause (b): forgiveness with the flip
        if bands.use_fallback_175:
            nb_ok = bands.fallback_175[0] <= abs(rho[1.75]) <= bands.fallback_175[1]
        else:
            nb_ok = abs(rho[2.0]) > abs(rho[1.75])
        b_ok = (rho[2.0] < -bands.b_forgive) and nb_ok
        det["b"] = b_ok
        if b_ok:
            # clause (c): shape SUP + contrast echo
            eta2 = bands.eta_nom[2.0]
            sup = max(abs(est["sigma"][r] - np.cos(r * np.pi / 2.0) * bands.eta_nom[r] / eta2)
                      for r in R_INTERIOR)
            contrast = (abs(rho[2.0]) - abs(rho[1.0])) / bands.se_contrast
            det["c_sup"], det["c_contrast"] = sup, contrast
            if sup < bands.b_curve and contrast > 5.0:
                return "CONFIRMED", det
    # DEVICE-DEVIATION: the frozen LR statistic
    lrs = lr_statistics(rho, bands)
    det["lr"] = lrs
    tripped = [k for k, v in lrs.items() if v > bands.lr_thresholds[k]]
    if tripped:
        det["lr_tripped"] = tripped
        return "DEVICE-DEVIATION", det
    return "INCONCLUSIVE", det


def full_pub_shot_data(kind, p, rng):
    """Sample the full science + guard PUB set of one flight under a generative
    model: 9 arms x 2 variants (r = 0 at 2x shots), the Z-leak pair, the 4 V_S
    guard PUBs, the 2 watcher-tomograph PUBs (round 30). Returns
    (shot_data, zleak_seqs, vs_seqs, tomo_seqs)."""
    shot_data = {r: sampled_arm_data(r, rng, dist_override=model_arm_dists(kind, r, p))
                 for r in R_ARMS}
    zl = {r: sample_shot_sequence(d, SHOTS, rng) for r, d in zleak_model_dists(kind, p).items()}
    vs = {key: sample_shot_sequence(d, SHOTS, rng) for key, d in vs_guard_model_dists(kind, p).items()}
    tomo = {r: sample_shot_sequence(d, SHOTS, rng) for r, d in tomograph_model_dists(kind, p).items()}
    return shot_data, zl, vs, tomo


def analyze_flight(shot_data, zleak_seqs, vs_seqs, tomo_seqs=None, mitigate=None):
    """The full analysis chain of one (real or simulated) flight: estimator with
    the pinned pipeline, mitigated guard statistics, the V_S statistic, the
    watcher-tomograph statistic. Returns (est, guards, vs_stat, tomo_stat)
    ready for verdict()."""
    mit = mitigate if mitigate is not None else (lambda d: d)
    est = estimator_from_shots(shot_data, mitigate=mitigate)
    arm_dists_all = {}
    for r in shot_data:
        if r == 0.0:
            for v in ("A", "B"):
                seq = shot_data[0.0][v]
                arm_dists_all.setdefault(0.0, {})[v] = mit(
                    counts_to_dist(np.bincount(seq[1::2], minlength=DIM)))
        else:
            arm_dists_all[r] = {v: mit(counts_to_dist(np.bincount(shot_data[r][v], minlength=DIM)))
                                for v in ("A", "B")}
    zl = {r: mit(counts_to_dist(np.bincount(s, minlength=DIM))) for r, s in zleak_seqs.items()}
    vs = {k: mit(counts_to_dist(np.bincount(s, minlength=DIM))) for k, s in vs_seqs.items()}
    guards = guard_stats(est, arm_dists_all, zl)
    tomo = np.nan
    if tomo_seqs is not None:
        td = {r: mit(counts_to_dist(np.bincount(s, minlength=DIM))) for r, s in tomo_seqs.items()}
        tomo = tomograph_stat(td)
    return est, guards, vs_guard_stat(vs), tomo


# ---- stage-1 self-test: the ideal chain against the law's face ----

def _selftest():
    arms = {r: arm_distributions(r) for r in R_ARMS}
    axes = fit_axes(arms[0.0])

    # the record lives EXACTLY on Y in the ideal circuit (zero X component):
    for name, w in (("j", Q_J), ("jp", Q_JP)):
        d = record_vector(arms[0.0], w)
        assert abs(d[0]) < 1e-12, f"ideal record of {name} must have zero X component: {d}"
        # b0 - b1 spans the full Bloch diameter at ideal contrast: |d_Y| = 2 (the 1/2 lives in S_hat)
        assert abs(abs(d[1]) - 2.0) < 1e-12, f"ideal record of {name} must be full on Y: {d}"
        s0 = signed_projection(arms[0.0], w, d / np.linalg.norm(d))
        assert abs(abs(s0) - 1.0) < 1e-12, f"ideal S_hat({name}; 0) must be unit: {s0}"

    rho = double_ratio(arms, axes)
    print("  r      rho_hat        cos(r*pi/2)    |diff|")
    worst = 0.0
    for r in R_ARMS:
        target = np.cos(r * np.pi / 2.0)
        diff = abs(rho[r] - target)
        worst = max(worst, diff)
        print(f"  {r:4.2f}  {rho[r]:+.12f}  {target:+.12f}  {diff:.2e}")
    assert worst < 1e-12, f"ideal rho_hat must equal the SIGNED law to machine precision: {worst:e}"
    assert rho[2.0] < -0.999999, f"THE FLIP: rho_hat(2) must be -1 ideal, got {rho[2.0]}"
    assert abs(rho[1.0]) < 1e-12, f"the odd fixpoint: rho_hat(1) must be 0 ideal, got {rho[1.0]}"

    # k and S sanity streams: <Z> ideal 0 on both (H then Z-diagonal gates)
    idx = np.arange(DIM)
    for q, nm in ((Q_K, "k"), (Q_S, "S")):
        z = 1 - 2 * ((idx >> q) & 1)
        for r in R_ARMS:
            for var in ("A", "B"):
                zval = float((arms[r][var] * z).sum())
                assert abs(zval) < 1e-12, f"<Z_{nm}> must be 0 ideal at r={r}/{var}: {zval}"

    # variant-B order is load-bearing: H-then-Sdg (the wrong order) must NOT recover Y
    psi = science_state(0.0)
    wrong = psi.copy()
    for q in (Q_J, Q_JP):
        wrong = _single(wrong, q, _H)
        wrong = _single(wrong, q, _SDG)
    p_wrong = np.abs(wrong) ** 2
    d_wrong = record_vector({"A": arm_distributions(0.0)["A"], "B": p_wrong}, Q_J)
    assert abs(d_wrong[1]) < 1e-12, f"the reversed pre-rotation must kill the Y record (silent-VOID trap pinned): {d_wrong}"

    print("stage-1 self-test PASS: ideal chain reproduces the SIGNED law exactly; flip at r=2;")
    print("record on Y only; k/S Z-streams zero; the reversed variant-B order is caught.")


def _selftest_sampling():
    rng = np.random.default_rng(20260804)            # deterministic; seeds recorded at freeze
    shot_data = {r: sampled_arm_data(r, rng) for r in R_ARMS}
    est = estimator_from_shots(shot_data)

    # finite-shot sanity at pinned shots on the NOISELESS circuit: the law within a loose
    # 5-sigma-ish window (SE of rho_hat ~ 1e-2 at 16384 shots); this is a smoke bound for
    # the chain's wiring, NEVER a band -- bands come from the signal model at the freeze.
    print("  r      rho_hat(sampled)  law        diff")
    for r in R_ARMS:
        d = est["rho"][r] - np.cos(r * np.pi / 2.0)
        print(f"  {r:4.2f}  {est['rho'][r]:+.6f}         {np.cos(r * np.pi / 2.0):+.6f}  {d:+.4f}")
        assert abs(d) < 0.05, f"sampled rho_hat({r}) off the law beyond smoke bound: {d}"
    assert est["rho"][2.0] < -0.9, f"THE FLIP at finite shots: {est['rho'][2.0]}"

    # determinism: the same seed reproduces the same estimator output bit for bit
    rng2 = np.random.default_rng(20260804)
    est2 = estimator_from_shots({r: sampled_arm_data(r, rng2) for r in R_ARMS})
    assert all(est2["rho"][r] == est["rho"][r] for r in R_ARMS), "seeded run must be reproducible"

    # the split is load-bearing: axis from the even half only, normalizer from the odd half
    # (fit and eval distributions differ at finite shots; identical halves would be a wiring bug)
    seqA = shot_data[0.0]["A"]
    assert not np.array_equal(np.bincount(seqA[0::2], minlength=DIM),
                              np.bincount(seqA[1::2], minlength=DIM)), "halves must differ at finite shots"
    print("stage-2 self-test PASS: sampled chain wired (split-sample axis/normalizer, determinism, flip).")


def _dist_double_ratio(kind, p):
    """Distribution-level rho_hat table under a generative model (exact, no
    sampling): the stage-3 workhorse for the exact checks."""
    arms = {r: model_arm_dists(kind, r, p) for r in R_ARMS}
    axes = fit_axes(arms[0.0])
    return double_ratio(arms, axes), arms, axes


def _selftest_stage3():
    ideal = Params()

    # (1) the density-matrix path reproduces the exact statevector reference
    worst_dm = 0.0
    for r in R_ARMS:
        sv = arm_distributions(r)
        dm = model_arm_dists("signal", r, ideal)
        for v in ("A", "B"):
            worst_dm = max(worst_dm, float(np.max(np.abs(sv[v] - dm[v]))))
    assert worst_dm < 1e-12, f"dm path must match the statevector reference: {worst_dm:e}"

    # (2) NULL: the (S, j) write zeroed kills j's record, j' stays full
    null0 = model_arm_dists("null", 0.0, ideal)
    dj = record_vector(null0, Q_J)
    djp = record_vector(null0, Q_JP)
    assert np.max(np.abs(dj)) < 1e-12, f"null must carry NO j record: {dj}"
    assert abs(abs(djp[1]) - 2.0) < 1e-12, f"null must keep the j' record full: {djp}"

    # (3) the three impostor curves, exact at distribution level
    targets = {
        "imp1": lambda r: np.cos(r * np.pi / 4.0) ** 2,
        "imp2": lambda r: np.cos(r * np.pi / 4.0),
        "imp3": lambda r: abs(np.cos(r * np.pi / 2.0)),
    }
    for kind, tgt in targets.items():
        rho, _, _ = _dist_double_ratio(kind, ideal)
        worst = max(abs(rho[r] - tgt(r)) for r in R_ARMS)
        assert worst < 1e-12, f"{kind} must realize its pinned curve exactly: {worst:e}"

    # (4a) coherent watcher offset + static ZZ on the watcher bond, exact algebra:
    # rho_hat(r>0) = cos(r*pi/2 + 2*eps + th)/cos(th); at r = 1 that is -sin(2e+th)/cos(th)
    e, th = 0.01, 0.01
    p = Params(eps_cal=e, zz_static={(Q_K, Q_J): th, (Q_J, Q_S): 0.02, (Q_S, Q_JP): 0.015})
    rho, _, _ = _dist_double_ratio("signal", p)
    for r in R_ARMS[1:]:
        expect = np.cos(r * np.pi / 2.0 + 2 * e + th) / np.cos(th)
        assert abs(rho[r] - expect) < 1e-12, f"coherent algebra at r={r}: {rho[r]} vs {expect}"
    assert abs(rho[1.0] + np.sin(2 * e + th) / np.cos(th)) < 1e-12, "the -sin(eps_w) center"

    # (4b) static write misses cancel in the double ratio (round-27: exact to 0.3 rad)
    p = Params(write_miss_sj=0.3, write_miss_sjp=-0.2)
    rho, _, _ = _dist_double_ratio("signal", p)
    worst = max(abs(rho[r] - np.cos(r * np.pi / 2.0)) for r in R_ARMS)
    assert worst < 1e-12, f"r-independent write misses must cancel exactly: {worst:e}"

    # (5a) static RZ frames are absorbed by the r = 0 axis fit; T_hat stays 0
    p = Params(rz_static={Q_K: 0.3, Q_J: 0.7, Q_S: -0.4, Q_JP: 1.1})
    rho, arms, axes = _dist_double_ratio("signal", p)
    worst = max(abs(rho[r] - np.cos(r * np.pi / 2.0)) for r in R_ARMS)
    assert worst < 1e-12, f"static frames must be absorbed by the axis fit: {worst:e}"
    tmax = max(abs(transverse_projection(arms[r], q, axes[w]))
               for r in R_ARMS for w, q in (("j", Q_J), ("jp", Q_JP)))
    assert tmax < 1e-12, f"T_hat must be 0 under statics: {tmax:e}"

    # (5b) inter-arm drift on j only: at r = 2 (slot 7, axis anchored at slot 5)
    # the tilt is 2*d -> T/S = tan(2d), rho_hat(2) = -cos(2d)
    d = 0.05
    p = Params(rz_drift={Q_J: d})
    rho, arms, axes = _dist_double_ratio("signal", p)
    s2 = signed_projection(arms[2.0], Q_J, axes["j"])
    t2 = transverse_projection(arms[2.0], Q_J, axes["j"])
    assert abs(abs(t2 / s2) - np.tan(2 * d)) < 1e-9, f"drift tan check: {t2/s2} vs {np.tan(2*d)}"
    assert abs(rho[2.0] + np.cos(2 * d)) < 1e-12, f"drifted rho_hat(2): {rho[2.0]}"
    # equal drift on BOTH witnesses: the cos factors cancel in the double ratio
    # (rho_hat blind to it) while T_hat still fires -> the guard is the detector
    p = Params(rz_drift={Q_J: d, Q_JP: d})
    rho, arms, axes = _dist_double_ratio("signal", p)
    assert abs(rho[2.0] + 1.0) < 1e-12, f"equal drift must cancel in rho_hat: {rho[2.0]}"
    t2 = transverse_projection(arms[2.0], Q_J, axes["j"])
    assert abs(t2) > 0.09, f"T_hat must fire on the drift rho_hat is blind to: {t2}"

    # (6a) exact-CAL mitigation recovers the law through worst-admitted confusion
    p = Params(readout={Q_S: (0.02, 0.02), Q_J: (0.03, 0.03), Q_JP: (0.03, 0.03), Q_K: (0.03, 0.03)})
    mats_true = {q: _confusion_matrix(*p.readout[q]) for q in range(N_QUBITS)}
    mit = make_mitigator(mats_true)
    arms = {r: {v: mit(dd) for v, dd in model_arm_dists("signal", r, p).items()} for r in R_ARMS}
    axes = fit_axes(arms[0.0])
    rho = double_ratio(arms, axes)
    worst = max(abs(rho[r] - np.cos(r * np.pi / 2.0)) for r in R_ARMS)
    assert worst < 1e-12, f"exact-CAL inversion must recover the law: {worst:e}"

    # (6b) the SYMMETRIC multiplicative readout mode cancels in the double ratio
    # even UNMITIGATED (the round-9 boundary of what the ratio can cancel)
    p = Params(bias_sym=0.02)
    rho, _, _ = _dist_double_ratio("signal", p)
    worst = max(abs(rho[r] - np.cos(r * np.pi / 2.0)) for r in R_ARMS)
    assert worst < 1e-12, f"symmetric readout mode must cancel: {worst:e}"

    # (6c) the ASYMMETRIC additive mode does NOT cancel and is ODD in the sign
    shifts = {}
    for delta in (+0.02, -0.02):
        rho, _, _ = _dist_double_ratio("signal", Params(bias_delta=delta))
        shifts[delta] = rho[1.75] - np.cos(1.75 * np.pi / 2.0)
    assert abs(shifts[0.02]) > 0.02, f"the additive mode must bias rho_hat(1.75): {shifts}"
    assert shifts[0.02] * shifts[-0.02] < 0, f"the bias must be odd in the sign: {shifts}"
    print(f"  additive S-antisymmetric assignment offset, rho_hat(1.75) shift at delta=+/-2%: "
          f"{shifts[0.02]:+.4f} / {shifts[-0.02]:+.4f} (doc: -/+0.043 at eta=0.92; ~0.039 at eta=1)")

    # (6d-i) the carrier's invariance, exact: a witness-local ARM-INDEPENDENT
    # affine readout map (asymmetric extra confusion, even UNMITIGATED) is
    # exactly invisible to the double ratio: this is what the carrier buys
    p = Params(readout={Q_J: (0.0, 0.02), Q_JP: (0.0, 0.02)})
    rho, _, _ = _dist_double_ratio("signal", p)
    worst = max(abs(rho[r] - np.cos(r * np.pi / 2.0)) for r in R_ARMS)
    assert worst < 1e-12, f"arm-independent affine readout must cancel exactly: {worst:e}"

    # (6d-ii) the ARM-DEPENDENT mode (coefficient ~ per-branch population
    # imbalance) survives exact-CAL mitigation: r-dependent, non-cancelling
    p = Params(readout={Q_S: (0.015, 0.02), Q_J: (0.02, 0.03), Q_JP: (0.02, 0.03), Q_K: (0.02, 0.03)},
               arm_bias_delta=0.02)
    mit = make_mitigator({q: _confusion_matrix(*p.readout[q]) for q in range(N_QUBITS)})
    arms = {r: {v: mit(dd) for v, dd in model_arm_dists("signal", r, p).items()} for r in R_ARMS}
    rho = double_ratio(arms, fit_axes(arms[0.0]))
    resid = abs(rho[2.0] - (-1.0))
    assert resid > 1e-3, f"the arm-dependent mode must leave a residual: {resid:e}"
    print(f"  arm-dependent mode residual on rho_hat(2) after exact-CAL mitigation: {resid:.5f}")

    # (7) the lambda-map: lam = 0 is exactly ideal; the solve hits its target
    p_inc = Params(epc_write=0.005, epc_watch=0.005,
                   t1={q: 140e-6 for q in range(N_QUBITS)},
                   t2={q: 70e-6 for q in range(N_QUBITS)}, idle_watcher=280e-9)
    p_inc.lam = 0.0
    assert abs(eta_of(p_inc) - 1.0) < 1e-12, "lam = 0 must be the ideal end"
    p_inc.lam = 1.0
    eta_worst = eta_of(p_inc)
    assert eta_worst < 1.0, f"worst-admitted incoherent config must attenuate: {eta_worst}"
    lam = solve_lambda(p_inc, 0.96)
    p_chk = Params(epc_write=0.005, epc_watch=0.005, t1=dict(p_inc.t1), t2=dict(p_inc.t2),
                   idle_watcher=280e-9, lam=lam)
    assert abs(eta_of(p_chk) - 0.96) < 1e-6, f"lambda solve must hit the target: {eta_of(p_chk)}"
    print(f"  lambda-map: eta(2) at worst-admitted lam=1 is {eta_worst:.5f}; "
          f"lam({0.96}) = {lam:.4f}")

    print("stage-3 self-test PASS: dm parity, null, 3 impostor curves exact, coherent")
    print("algebra exact, write-miss cancel, statics absorbed, drift tan + guard-only")
    print("drift, exact-CAL recovery, symmetric-mode cancel, additive-mode odd bias,")
    print("arm-dependent residual, lambda-map solve.")


def _selftest_stage3_sampled():
    """End-to-end wiring smoke: worst-admitted SIGNAL through sampling + finite
    CAL + mitigation + the pinned estimator. Loose bounds only; bands are stage
    4's job; this asserts the chain is executable and lands near the dressed law."""
    p = Params.worst_admitted()
    rng = np.random.default_rng(20260804)
    mit = sampled_cal_mitigator(p, rng)
    shot_data = {r: sampled_arm_data(r, rng, dist_override=model_arm_dists("signal", r, p))
                 for r in R_ARMS}
    est = estimator_from_shots(shot_data, mitigate=mit)
    eta2 = eta_of(p)
    print(f"  worst-admitted SIGNAL, sampled + mitigated (eta(2)_incoherent = {eta2:.4f}):")
    print("  r      rho_hat     cos*eta2ref")
    for r in R_ARMS:
        print(f"  {r:4.2f}  {est['rho'][r]:+.4f}     {np.cos(r * np.pi / 2.0) * eta2:+.4f}")
    assert est["rho"][2.0] < -0.8, f"THE FLIP must survive the full noise budget: {est['rho'][2.0]}"
    assert abs(est["rho"][1.0]) < 0.15, f"blindness smoke at r=1: {est['rho'][1.0]}"
    assert est["rho"][0.25] > 0.75, f"early arm smoke: {est['rho'][0.25]}"
    # impostor (i) through the same sampled chain: no flip, near-zero at r = 2
    shot_imp = {r: sampled_arm_data(r, rng, dist_override=model_arm_dists("imp1", r, p))
                for r in R_ARMS}
    est_i = estimator_from_shots(shot_imp, mitigate=mit)
    assert est_i["rho"][2.0] > -0.2, f"imp1 must NOT flip: {est_i['rho'][2.0]}"
    print(f"  imp1 sampled rho_hat(2) = {est_i['rho'][2.0]:+.4f} (no flip)  "
          f"signal rho_hat(2) = {est['rho'][2.0]:+.4f} (flip)")
    print("stage-3 sampled smoke PASS: full chain executable, flip present under the")
    print("full worst-admitted budget, impostor separated.")


# ---- stage 4b groundwork: the counts-level fast path for seed banks ----
# Band freezing needs 1e5 (physics) / 1e6 (guards) seeds; per-shot sequence
# sampling is ~1000x too slow for that. The even/odd halves of an iid delivered
# sequence are EXACTLY two independent multinomials of half the shots each, so
# the counts path below is the split rule's exact sampling image (the sequence
# path above stays the reference; the equivalence is distributional, asserted
# by the moment check in the stage-4b tests, and the estimator code downstream
# of the counts is THE SAME code).

def sampled_arm_counts(r, rng, dists, shots=None):
    """Multinomial counts per variant. At r = 0 returns (fit_counts,
    eval_counts) per variant (two independent multinomials of half the shots);
    other arms return one counts vector per variant."""
    n = shots if shots is not None else (SHOTS_R0 if r == 0.0 else SHOTS)
    out = {}
    for v in ("A", "B"):
        pd = np.clip(dists[v], 0.0, None)
        pd = pd / pd.sum()
        if r == 0.0:
            out[v] = (rng.multinomial(n // 2, pd).astype(float),
                      rng.multinomial(n - n // 2, pd).astype(float))
        else:
            out[v] = rng.multinomial(n, pd).astype(float)
    return out


def estimator_from_counts(counts_data, mitigate=None):
    """The identical estimator chain on counts-level samples: counts ->
    normalize -> CAL inversion -> conditioning, per the pinned order.
    counts_data: {0.0: {v: (fit_counts, eval_counts)}, r>0: {v: counts}}."""
    mit = mitigate if mitigate is not None else (lambda d: d)
    fit_dists = {v: mit(counts_to_dist(counts_data[0.0][v][0])) for v in ("A", "B")}
    eval_dists = {v: mit(counts_to_dist(counts_data[0.0][v][1])) for v in ("A", "B")}
    axes = fit_axes(fit_dists)
    arm_dists = {r: {v: mit(counts_to_dist(counts_data[r][v])) for v in ("A", "B")}
                 for r in counts_data if r != 0.0}
    arm_dists[0.0] = eval_dists
    s_j0 = signed_projection(arm_dists[0.0], Q_J, axes["j"])
    d_jp0 = abs(signed_projection(arm_dists[0.0], Q_JP, axes["jp"]))
    ref = s_j0 / d_jp0
    rho, t_guards = {}, {}
    for r, dists in arm_dists.items():
        s_j = signed_projection(dists, Q_J, axes["j"])
        d_jp = abs(signed_projection(dists, Q_JP, axes["jp"]))
        rho[r] = (s_j / d_jp) / ref
        t_guards[r] = {w: transverse_projection(dists, q, axes[w])
                       for w, q in (("j", Q_J), ("jp", Q_JP))}
    sigma_hat = {r: rho[r] / abs(rho[2.0]) for r in rho}
    return {"axes": axes, "rho": rho, "sigma": sigma_hat, "t": t_guards,
            "s_j0": s_j0, "d_jp0": d_jp0, "arm_dists": arm_dists}


def seed_bank_rho(kind, p, n_seeds, rng, mitigate=None):
    """A seed bank of rho_hat vectors over R_LR under a generative model,
    counts-level (the stage-4b workhorse). Returns (n_seeds, 8) array."""
    dists = {r: model_arm_dists(kind, r, p) for r in R_ARMS}
    bank = np.empty((n_seeds, len(R_LR)))
    for i in range(n_seeds):
        cd = {r: sampled_arm_counts(r, rng, dists[r]) for r in R_ARMS}
        est = estimator_from_counts(cd, mitigate=mitigate)
        bank[i] = [est["rho"][r] for r in R_LR]
    return bank


# ---- stage 4b: band freezing through the frozen verdict code ----
# The pinned procedure: SIGNAL model sources every fluctuation width; magnitude
# bands + floors + guard bands at the WORST-ADMITTED basis (lambda solved for
# eta(2) = eta_min); the clause-(c) target eta_nom at the MID-RANGE basis
# (lambda for (eta_min+1)/2; Aer parity on the representative line REPLACES this
# stand-in read at the runner boundary); b_curve corner-conditional over the 2-D
# level x slope grid x the readout-bias sign axis; per-statistic worst-sign
# envelopes per rounds 23-27; empirical quantiles + bootstrap SE (the round-18
# NORM branch; branch recorded); held-out halves measure power/controls.

COLS = [f"rho_{r}" for r in R_LR] + ["s_j0", "d_jp0", "d_jp_min",
                                     "t_j", "t_jp", "z_k", "z_s", "z_leak", "vs", "rsym",
                                     "tomo"]
ICOL = {c: i for i, c in enumerate(COLS)}


def _normp(d):
    v = np.clip(d, 0.0, None)
    return v / v.sum()


def draw_device(base, rng, spread=0.25):
    """Per-seed j-vs-j' asymmetry from the calibration spread (stand-in model,
    recorded for the freeze-time empty round: independent uniform +/-spread
    factors on the two write-edge EPCs, the witness/spectator T2s, and the
    witness readout error pairs)."""
    from dataclasses import replace
    f = lambda: 1.0 + rng.uniform(-spread, spread)
    t2 = {q: t / f() for q, t in base.t2.items()}
    readout = dict(base.readout)
    for w in (Q_J, Q_JP):
        e01, e10 = readout.get(w, (0.0, 0.0))
        g = f()
        readout[w] = (e01 * g, e10 * g)
    return replace(base,
                   epc_write_sj=base.epc_write * f(), epc_write_sjp=base.epc_write * f(),
                   t2=t2, readout=readout)


def flight_record(kind, p, rng):
    """One simulated flight (counts level, per-seed CAL) condensed to the COLS
    row: the 8-arm rho_hat vector + every guard statistic the verdict consumes."""
    dists = {r: model_arm_dists(kind, r, p) for r in R_ARMS}
    cd = {r: sampled_arm_counts(r, rng, dists[r]) for r in R_ARMS}
    mit = sampled_cal_mitigator(p, rng)
    est = estimator_from_counts(cd, mitigate=mit)
    zl = {r: mit(counts_to_dist(rng.multinomial(SHOTS, _normp(d)).astype(float)))
          for r, d in zleak_model_dists(kind, p).items()}
    vs_d = {k: mit(counts_to_dist(rng.multinomial(SHOTS, _normp(d)).astype(float)))
            for k, d in vs_guard_model_dists(kind, p).items()}
    guards = guard_stats(est, est["arm_dists"], zl)
    vs_stat = vs_guard_stat(vs_d)
    tomo_d = {r: mit(counts_to_dist(rng.multinomial(SHOTS, _normp(d)).astype(float)))
              for r, d in tomograph_model_dists(kind, p).items()}
    tomo = tomograph_stat(tomo_d)
    ax = est["axes"]
    dj = {r: abs(signed_projection(est["arm_dists"][r], Q_J, ax["j"])) for r in (0.0, 2.0)}
    djp = {r: abs(signed_projection(est["arm_dists"][r], Q_JP, ax["jp"])) for r in (0.0, 2.0)}
    rsym = abs(dj[2.0] / djp[2.0] - dj[0.0] / djp[0.0])
    row = [est["rho"][r] for r in R_LR]
    row += [guards["s_j0"], guards["d_jp0"], guards["d_jp_min"],
            guards["t_sup"]["j"], guards["t_sup"]["jp"], guards["z_k"], guards["z_s"],
            guards["z_leak"], vs_stat, rsym, tomo]
    return np.array(row)


def _bank_chunk(args):
    kind, p, n, seed, spread = args
    rng = np.random.default_rng(seed)
    out = np.empty((n, len(COLS)))
    for i in range(n):
        out[i] = flight_record(kind, draw_device(p, rng, spread), rng)
    return out


def run_bank(kind, p, n, seed, spread=0.25, jobs=1):
    """A seed bank of flight records. jobs > 1 fans chunks over processes
    (deterministic per-chunk seeds; the full-scale freeze runs need it)."""
    if jobs <= 1:
        return _bank_chunk((kind, p, n, seed, spread))
    from concurrent.futures import ProcessPoolExecutor
    sizes = [n // jobs + (1 if i < n % jobs else 0) for i in range(jobs)]
    args = [(kind, p, sz, seed + 1000 * i, spread) for i, sz in enumerate(sizes) if sz > 0]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        return np.vstack(list(ex.map(_bank_chunk, args)))


def _q_boot(x, q, rng, n_boot=200):
    """Empirical quantile + bootstrap SE of the quantile (the pinned NORM
    branch; the parametric ladder is recorded as not taken)."""
    v = float(np.quantile(x, q))
    n = len(x)
    bs = np.empty(n_boot)
    for i in range(n_boot):
        bs[i] = np.quantile(x[rng.integers(0, n, n)], q)
    return v, float(bs.std())


def record_verdict(row, bands):
    """THE SAME verdict() on a condensed flight record (rho_hat(0) = 1 by
    construction: the r = 0 eval half IS the reference)."""
    rho = {0.0: 1.0}
    for i, r in enumerate(R_LR):
        rho[r] = row[i]
    a2 = abs(rho[2.0])
    est = {"rho": rho, "sigma": {r: rho[r] / a2 for r in rho},
           "s_j0": row[ICOL["s_j0"]], "d_jp0": row[ICOL["d_jp0"]]}
    guards = {"t_sup": {"j": row[ICOL["t_j"]], "jp": row[ICOL["t_jp"]]},
              "z_k": row[ICOL["z_k"]], "z_s": row[ICOL["z_s"]],
              "z_leak": row[ICOL["z_leak"]], "s_j0": row[ICOL["s_j0"]],
              "d_jp0": row[ICOL["d_jp0"]], "d_jp_min": row[ICOL["d_jp_min"]]}
    return verdict(est, guards, row[ICOL["vs"]], bands, vs_threshold=bands.vs_threshold,
                   tomo_stat=row[ICOL["tomo"]])


@dataclass
class FreezeConfig:
    """Scales: the pinned freeze runs at n_phys=1e5, n_guard=1e6, n_power=1e5
    per grid point, levels=7, slopes=5; the validate mode runs the same
    machinery at reduced scale to prove wiring (its numbers NEVER fly)."""
    eta_min: float = 0.92
    signs: tuple = (-0.02, 0.0, +0.02)
    levels: int = 7
    slopes: int = 5
    s_res: float = 0.01        # stand-in until the representative-line schedule read
    spread: float = 0.25
    n_phys: int = 100_000
    n_guard: int = 1_000_000
    n_power: int = 100_000
    n_null: int = 100_000
    n_boot: int = 200
    seed: int = 20260808
    jobs: int = 1
    aer_report: str = None     # path to _record_parity_aer_report.json: when set,
                               # eta_nom comes from the Aer read (the pinned source)
                               # and s_res from its duration-slope reading; the dm
                               # stand-in is used only when this is None (validate)


def _lr_rows(bank, bands):
    """LR values per impostor for every row of a bank (vectorized)."""
    x = bank[:, :len(R_LR)]
    eta = np.array([bands.eta_nom[r] for r in R_LR])
    sinv = np.linalg.inv(bands.lr_sigma)
    mu_law = np.array([np.cos(r * np.pi / 2.0) for r in R_LR]) * eta
    d_law = x - mu_law
    q_law = np.einsum("ni,ij,nj->n", d_law, sinv, d_law)
    out = {}
    for name, curve in IMPOSTOR_CURVES.items():
        mu_i = np.array([curve(r) for r in R_LR]) * eta
        d_i = x - mu_i
        out[name] = 0.5 * (q_law - np.einsum("ni,ij,nj->n", d_i, sinv, d_i))
    return out


def _sup_shape(bank, eta_nom):
    """Clause (c)'s SUP statistic per row."""
    i2 = R_LR.index(2.0)
    a2 = np.abs(bank[:, i2])
    sup = np.zeros(len(bank))
    for r in R_INTERIOR:
        i = R_LR.index(r)
        tgt = np.cos(r * np.pi / 2.0) * eta_nom[r] / eta_nom[2.0]
        sup = np.maximum(sup, np.abs(bank[:, i] / a2 - tgt))
    return sup


def _signed_base(cfg, lam):
    p = Params.worst_admitted()
    p.lam = lam
    return p


def freeze_bands(cfg):
    """The full stage-4b procedure. Returns (Bands, report dict). Every choice
    that the pre-reg pins is executed here THROUGH the verdict code; every
    stand-in (eta_nom source, s_res, the spread model) is named in the report
    for the freeze-time empty round."""
    from dataclasses import replace
    rng = np.random.default_rng(cfg.seed + 999)
    rep = {"scale": {"n_phys": cfg.n_phys, "n_guard": cfg.n_guard, "n_power": cfg.n_power,
                     "n_null": cfg.n_null, "levels": cfg.levels, "slopes": cfg.slopes},
           "stand_ins": ["eta_nom from the dm model (Aer parity on the representative line replaces it)",
                         "s_res floor 1% (representative-line schedule read replaces it)",
                         "calibration-spread draw: uniform +/-25% on edge EPCs, T2s, witness readout"],
           "tail_branch": "empirical + bootstrap SE (parametric ladder not taken; recorded)"}

    inc = Params(epc_write=0.005, epc_watch=0.005,
                 t1={q: 140e-6 for q in range(N_QUBITS)},
                 t2={q: 70e-6 for q in range(N_QUBITS)}, idle_watcher=280e-9)
    lam_worst = solve_lambda(inc, cfg.eta_min)
    eta_mid_target = (cfg.eta_min + 1.0) / 2.0
    lam_mid = solve_lambda(inc, eta_mid_target)
    if cfg.aer_report:
        import json as _json
        aer = _json.load(open(cfg.aer_report, encoding="utf-8"))
        eta_nom = {float(k): v for k, v in aer["eta_nom"].items()}
        rep["eta_nom_source"] = f"Aer read: {cfg.aer_report} ({aer.get('pulled_utc')})"
        rep["stand_ins"] = [s for s in rep["stand_ins"] if not s.startswith("eta_nom")]
    else:
        p_mid = replace(inc, lam=lam_mid)
        eta_nom = eta_attenuation_profile(p_mid)
        rep["eta_nom_source"] = "dm-model stand-in (validate only; the freeze uses the Aer read)"
    rep["lam_worst"], rep["lam_mid"] = lam_worst, lam_mid
    rep["eta_nom"] = {str(r): eta_nom[r] for r in R_ARMS}

    # worst-basis signal banks per sign (freeze) + held-out twins (measure)
    banks, banks_ho = {}, {}
    for s in cfg.signs:
        p_s = replace(_signed_base(cfg, lam_worst), bias_delta=s, arm_bias_delta=s)
        banks[s] = run_bank("signal", p_s, cfg.n_phys, cfg.seed + int(s * 1e4),
                            cfg.spread, cfg.jobs)
        banks_ho[s] = run_bank("signal", p_s, cfg.n_phys, cfg.seed + 77777 + int(s * 1e4),
                               cfg.spread, cfg.jobs)
    guard_banks = banks if cfg.n_guard == cfg.n_phys else {
        s: run_bank("signal",
                    replace(_signed_base(cfg, lam_worst), bias_delta=s, arm_bias_delta=s),
                    cfg.n_guard, cfg.seed + 333 + int(s * 1e4), cfg.spread, cfg.jobs)
        for s in cfg.signs}

    def col(bank, name):
        return bank[:, ICOL[name]]

    i1, i2, i175 = R_LR.index(1.0), R_LR.index(2.0), R_LR.index(1.75)
    # fixed-j'-corner banks (round 29, stats M3: the per-seed asymmetry DRAW is
    # a mixture, and a device has ONE fixed j'; the floors' lower envelope must
    # include the fixed device at the draw's own j' extremes, else a bad-j'
    # valid device sees a floor frozen from the mixture: 0.03% -> ~0.26%
    # terminal false-VOID per floor at a +/-1% spread)
    corner_banks, corner_ho = {}, {}
    for s in cfg.signs:
        pc = replace(_signed_base(cfg, lam_worst), bias_delta=s, arm_bias_delta=s)
        pc.readout = dict(pc.readout)
        e01, e10 = pc.readout[Q_JP]
        pc.readout[Q_JP] = (e01 * (1 + cfg.spread), e10 * (1 + cfg.spread))
        pc.t2 = dict(pc.t2)
        pc.t2[Q_JP] /= (1 + cfg.spread)
        pc.epc_write_sjp = pc.epc_write * (1 + cfg.spread)
        corner_banks[s] = run_bank("signal", pc, max(cfg.n_phys // 4, 500),
                                   cfg.seed + 4242 + int(s * 1e4), 0.0, cfg.jobs)
        # held-out twin: the corners enter the GATES too (round 30, stats M2:
        # the round-29 fix reached the BANDS but not the gates -- power and
        # false-VOID were still measured on the DRAW MIXTURE only; general law:
        # every device-fixed nuisance enters as a CORNER, only shot noise is
        # drawn)
        corner_ho[s] = run_bank("signal", pc, max(cfg.n_phys // 4, 500),
                                cfg.seed + 424242 + int(s * 1e4), 0.0, cfg.jobs)
    # b_forgive: LOWER ENVELOPE over the sign axis (min per-sign p0.13 of
    # |rho_2|), SE to the upper edge (never-ease lives in the +SE edge; the
    # "CONFIRM-HARDENING sign" label is RETIRED, round 30 spec M5/stats M6 --
    # for a lower cut the sign axis is a coverage axis, no sign name needed)
    bf = {s: _q_boot(np.abs(banks[s][:, i2]), 0.0013, rng, cfg.n_boot) for s in cfg.signs}
    s_hard = min(bf, key=lambda s: bf[s][0])
    b_forgive = bf[s_hard][0] + bf[s_hard][1]
    rep["b_forgive_per_sign"] = {str(s): bf[s] for s in cfg.signs}
    rep["b_forgive_envelope_argmin_sign"] = s_hard
    # SE_contrast: max std over signs of (|rho_2| - |rho_1|)
    se_contrast = max(float(np.std(np.abs(banks[s][:, i2]) - np.abs(banks[s][:, i1])))
                      for s in cfg.signs)
    # LR covariance from the HELD-OUT zero-sign worst-basis bank
    sigma = np.cov(banks_ho[0.0][:, :len(R_LR)].T)

    # the 2-D level x slope grid x sign axis: freeze half / measure half
    levels = np.linspace(cfg.eta_min, 1.0, cfg.levels)
    slopes = np.linspace(-cfg.s_res, cfg.s_res, cfg.slopes)
    lam_of = {lv: (0.0 if lv >= 1.0 else solve_lambda(inc, lv)) for lv in levels}
    rep["levels"], rep["slopes"] = list(map(float, levels)), list(map(float, slopes))
    grid = {}
    for il, lv in enumerate(levels):
        for isl, sl in enumerate(slopes):
            for isg, s in enumerate(cfg.signs):
                p_g = replace(_signed_base(cfg, lam_of[lv]), eta_slope=sl,
                              bias_delta=s, arm_bias_delta=s)
                key = (round(float(lv), 4), round(float(sl), 4), s)
                grid[key] = run_bank("signal", p_g, cfg.n_power,
                                     cfg.seed + 10_000 + il * 1009 + isl * 101 + isg * 11,
                                     cfg.spread, cfg.jobs)
    # the OFF-DIAGONAL BOX CORNER joins the grid (round 30, stats M3: the
    # lambda one-knob family is a CURVE through the admitted box -- it ties
    # all incoherent rates together -- while the corner BEST COHERENCE (lam=0)
    # x WORST readout/systematics FIXED at the draw's admitted extremes is an
    # admitted device the curve-with-draw only ever visits transiently; the
    # guard-envelope arithmetic read up to ~9.6% terminal VOID there under an
    # unpinned width decomposition; as a grid member it feeds every envelope,
    # b_curve, the LR reference, and the power/false-VOID gates): slope at
    # both extremes, sign axis full, spread = 0 (device-fixed nuisance = a
    # corner, only shot noise drawn)
    for isl, sl in enumerate((-cfg.s_res, +cfg.s_res)):
        for isg, s in enumerate(cfg.signs):
            p_b = replace(_signed_base(cfg, 0.0), eta_slope=sl,
                          bias_delta=s, arm_bias_delta=s)
            p_b.readout = {q: (e01 * (1 + cfg.spread), e10 * (1 + cfg.spread))
                           for q, (e01, e10) in p_b.readout.items()}
            key = ("box", round(float(sl), 4), s)
            grid[key] = run_bank("signal", p_b, cfg.n_power,
                                 cfg.seed + 60_000 + isl * 101 + isg * 11,
                                 0.0, cfg.jobs)
    # the ETA-TOP extension past 1 (round 30, physics m3: eta = a_j / a_j' and
    # no abort enforces a_j <= a_j' -- a j'-side differential puts eta(2)
    # STRUCTURALLY above 1, excess ~0.001-0.002; the level grid's lambda map
    # cannot realize it, so the top point is realized by the j'-only
    # differential config: ideal elsewhere, j' dephasing during the watcher
    # idle only): t2(j') = 140 us at 280 ns idle -> eta(2) ~ 1.002
    for isg, s in enumerate(cfg.signs):
        p_t = Params(t2={Q_JP: 140e-6}, t1={Q_JP: 1.0}, idle_watcher=280e-9,
                     eps_cal=0.01,
                     zz_static={(Q_K, Q_J): 0.01, (Q_J, Q_S): 0.01, (Q_S, Q_JP): 0.01},
                     rz_drift={Q_J: 0.05, Q_JP: 0.05},
                     write_miss_sj=0.01, write_miss_sjp=0.01,
                     readout={Q_S: (0.015, 0.02), Q_J: (0.02, 0.03),
                              Q_JP: (0.02, 0.03), Q_K: (0.02, 0.03)},
                     bias_delta=s, arm_bias_delta=s)
        key = ("etatop", 0.0, s)
        grid[key] = run_bank("signal", p_t, cfg.n_power,
                             cfg.seed + 70_000 + isg * 11, cfg.spread, cfg.jobs)
    rep["eta_top_measured"] = eta_of(Params(t2={Q_JP: 140e-6}, t1={Q_JP: 1.0},
                                            idle_watcher=280e-9))
    # the COHERENT eps_w SIGN axis (round 31, physics MAJOR: |rho_hat(1.75)| =
    # cos(pi/8 -/+ eps_w) is FIRST order in the watcher-angle sign -- the (b)
    # raw-order gap swings ~ +/-0.0115 at budget and the NEGATIVE sign EASES
    # the CONFIRMED-bearing comparison, the same class as the drift finding;
    # rho_hat(2) is second-order flat, b_forgive untouched; the negative-eps_w
    # banks join the grid so the fallback union, the order-power report,
    # b_curve, the LR reference, and the gates all see the sign)
    for isg, s in enumerate(cfg.signs):
        p_e = replace(_signed_base(cfg, lam_worst), eps_cal=-0.01,
                      bias_delta=s, arm_bias_delta=s)
        p_e.zz_static = {k: -v for k, v in p_e.zz_static.items()}
        key = ("epsw-", 0.0, s)
        grid[key] = run_bank("signal", p_e, cfg.n_power,
                             cfg.seed + 80_000 + isg * 11, cfg.spread, cfg.jobs)

    # guard bands: p99.97 + SE, UPPER ENVELOPE over the sign axis AND the full
    # coverage range (worst-basis banks + every grid freeze half). The v29
    # single-basis pin is AMENDED here (design finding #7, caught by the
    # sign-and-level-swept power gate 2026-08-04): T-hat's drift term scales
    # with the RECORD LEVEL, so a worst-basis-only band terminally VOIDs the
    # better-than-worst valid device (corner breakdown: t_sup_j 269/300).
    # Floors stay worst-basis: for a LOWER cut the worst basis IS the envelope.
    gbands = {}
    guard_sources = [guard_banks[s] for s in cfg.signs] + \
                    [bank[: len(bank) // 2] for bank in grid.values()] + \
                    list(corner_banks.values())
    guard_source_labels = [f"worst-basis/{s}" for s in cfg.signs] + \
                          [f"grid/{k}" for k in grid] + \
                          [f"fixed-jp-corner/{s}" for s in cfg.signs]
    guard_argmax = {}
    for name in ("t_j", "t_jp", "z_k", "z_s", "z_leak", "vs"):
        vals = [_q_boot(col(b, name), 0.9997, rng, cfg.n_boot) for b in guard_sources]
        edges = [v + se for v, se in vals]
        gbands[name] = max(edges)
        # the envelope-argmax SOURCE per band (round 30, stats m11: the pinned
        # seed-cost scheme deep-seeds ONLY the envelope argmax; the report
        # names it so the binding freeze can re-run that point at n_guard if
        # it is a 1e5 grid point rather than a deep worst-basis bank)
        guard_argmax[name] = guard_source_labels[int(np.argmax(edges))]
    rep["guard_bands"] = gbands
    rep["guard_band_argmax_source"] = guard_argmax
    rep["guard_band_source"] = ("coverage-range envelope: worst-basis banks + level x slope x "
                                "sign grid freeze halves (incl. the box corner and the eta-top "
                                "point) + fixed-j'-corner banks")
    # the z_k DISCRIMINATION CEILING (round 30, stats B1): the coverage
    # envelope must sit at or below the pinned maximum -- exceedance is
    # UNFREEZABLE (design remedy, never a loosened ceiling)
    assert gbands["z_k"] <= Z_K_BAND_CEILING, \
        f"z_k envelope {gbands['z_k']:.4f} exceeds the pinned discrimination " \
        f"ceiling {Z_K_BAND_CEILING} -- unfreezable as written (round 30, stats B1)"
    # watcher-tomograph mechanism threshold (round 30, physics BLOCKER):
    # a LOWER cut, diagnostic-only; frozen from the p0.03 lower envelope
    # minus SE over the same coverage sources (a false trip only muzzles the
    # mechanism sentence, so the deep tail is cheap)
    tomo_vals = [_q_boot(col(b, "tomo"), 0.0003, rng, cfg.n_boot) for b in guard_sources]
    tomo_threshold = min(v - se for v, se in tomo_vals)
    rep["tomo_threshold"] = tomo_threshold
    # b_blind: upper cut -> the same coverage envelope (round 29, stats m1: the
    # worst-basis freeze was coverage-correct only by an unpinned cancellation;
    # the center |tan(eps_w)|*eta(1)/... rides the record level like T-hat did);
    # per-source quantile SE-shifted to its LOWER edge, envelope = max
    bb_vals = [_q_boot(np.abs(b[:, i1]), 0.9987, rng, cfg.n_boot) for b in guard_sources]
    b_blind = max(v - se for v, se in bb_vals)
    rep["b_blind"] = b_blind
    # floors: lower envelope at p0.03 - SE over sign banks, grid halves (the
    # slope axis, round 29 stats m3), and the fixed-j'-corner banks (stats M3)
    floors = {}
    for name in ("s_j0", "d_jp0", "d_jp_min"):
        fq = [_q_boot(col(b, name), 0.0003, rng, cfg.n_boot) for b in guard_sources]
        floors[name] = min(v - se for v, se in fq)
    # the round-25 pin on the s_j0 floor (>= 10 SHOT-SE below the worst valid
    # mean) is made true BY CONSTRUCTION (round 29: a p0.03 envelope alone sits
    # only ~3.4 of its own SDs below the worst source bank's mean, so the pin
    # and the envelope must be combined, not merely checked): the floor takes
    # the LOWER of the envelope and (worst valid mean - 10 shot-SE); lowering
    # is the fewer-false-void side, and the null routing keeps ~120 shot-SE
    valid_means = [float(np.mean(b[:, ICOL["s_j0"]])) for b in
                   list(banks_ho.values()) + [bk[: len(bk) // 2] for bk in grid.values()]]
    shot_se_s0 = min(float(np.std(b[:, ICOL["s_j0"]])) for b in corner_banks.values())
    floors["s_j0"] = min(floors["s_j0"], min(valid_means) - 10.0 * shot_se_s0)
    rep["floors"] = floors
    rep["s_j0_floor_construction"] = {"envelope": float(min(v - se for v, se in
                                      [_q_boot(col(b, "s_j0"), 0.0003, rng, cfg.n_boot)
                                       for b in guard_sources])),
                                      "worst_valid_mean": min(valid_means),
                                      "shot_se": shot_se_s0}

    bands = Bands(
        b_blind=b_blind, b_forgive=b_forgive, b_curve=np.inf,
        se_contrast=se_contrast, eta_nom=eta_nom, lr_sigma=sigma,
        lr_thresholds={k: np.inf for k in IMPOSTOR_CURVES},
        floor_s_j0=floors["s_j0"], floor_d_jp0=floors["d_jp0"],
        floor_d_jp_min=floors["d_jp_min"],
        guard_t_sup=max(gbands["t_j"], gbands["t_jp"]),
        guard_z_k=gbands["z_k"], guard_z_s=gbands["z_s"], guard_z_leak=gbands["z_leak"],
    )
    bands.vs_threshold = gbands["vs"]
    bands.tomo_threshold = tomo_threshold
    # the fixed-j'-corner banks join the SHAPE-band constructions too (round 30:
    # the corners are coverage-grid members for every frozen constant, not only
    # the floors/guards -- a fixed-j' differential moves the inner-ratio
    # profile the SUP and the LR read); guard_sources already carries them
    grid_and_corners = dict(grid)
    for s in cfg.signs:
        grid_and_corners[("jpcorner", 0.0, s)] = corner_banks[s]
    # b_curve: corner-conditional, sourced from the (b)-reaching subset of the
    # freeze halves, per-point p99.87 minus SE, MAX over the whole grid
    b_curve_pts = {}
    for key, bank in grid_and_corners.items():
        half = bank[: len(bank) // 2]
        reach = []
        for row in half:
            v, det = record_verdict(row, bands)     # b_curve=inf: CONFIRMED iff reaches (c)
            if v != "VOID" and det.get("a") and det.get("b"):
                reach.append(row)
        if len(reach) < 50:
            continue
        sup = _sup_shape(np.array(reach), eta_nom)
        v, se = _q_boot(sup, 0.9987, rng, cfg.n_boot)
        b_curve_pts[key] = (v, se)
    bands.b_curve = max(v - se for v, se in b_curve_pts.values())
    rep["b_curve"] = bands.b_curve
    rep["b_curve_corner"] = str(max(b_curve_pts, key=lambda k: b_curve_pts[k][0] - b_curve_pts[k][1]))

    # LR thresholds: CORNER-CONDITIONAL over the FULL level x slope x sign grid
    # (round 29, stats M1: the pooled-per-sign reference under-covered the
    # slope corner exactly like the round-17 pooled b_curve; a valid device at
    # the -1% slope corner read family DEVICE-DEVIATION at ~1% vs the 0.13%
    # pin, worsening with shots -- pooled-into-reference does NOT fix it, only
    # max-over-grid does): per-point p99.87 + SE, threshold = the grid maximum,
    # then family-tightened on the held-out halves
    lr_pts = {k: [] for k in IMPOSTOR_CURVES}
    for key, bank in grid_and_corners.items():
        vals = _lr_rows(bank[: len(bank) // 2], bands)
        for k in IMPOSTOR_CURVES:
            v, se = _q_boot(vals[k], 0.9987, rng, cfg.n_boot)
            lr_pts[k].append(v + se)
    bands.lr_thresholds = {k: max(vs) for k, vs in lr_pts.items()}
    ho_rows = np.vstack([bank[len(bank) // 2:] for bank in grid.values()])
    for _ in range(20):
        lrs = _lr_rows(ho_rows, bands)
        fam = np.zeros(len(ho_rows), dtype=bool)
        for k in IMPOSTOR_CURVES:
            fam |= lrs[k] > bands.lr_thresholds[k]
        rate = fam.mean()
        if rate <= 0.0013:
            break
        bands.lr_thresholds = {k: t * 1.2 for k, t in bands.lr_thresholds.items()}
    rep["lr_thresholds"] = dict(bands.lr_thresholds)
    rep["lr_family_exceedance"] = float(rate)

    # clause-(b) raw order-check power across the grid; fallback engagement
    order_pw = {key: float(np.mean(np.abs(bank[:, i2]) > np.abs(bank[:, i175])))
                for key, bank in grid_and_corners.items()}
    min_order = min(order_pw.values())
    rep["order_check_min_power"] = min_order
    # clause (a)'s parameter-free order check gets the same power REPORT
    # (round 31, spec m3: the structurally identical check in (b) carries a
    # power report + fallback; (a)'s margins are enormous -- ~0.9 vs ~0.03 --
    # so a report suffices, and this measures it rather than asserting it)
    i075, i125 = R_LR.index(0.75), R_LR.index(1.25)
    order_a_pw = {key: float(np.mean(
        np.minimum(np.abs(bank[:, i075]), np.abs(bank[:, i125])) > np.abs(bank[:, i1])))
        for key, bank in grid_and_corners.items()}
    rep["order_a_min_power"] = min(order_a_pw.values())
    if min_order >= 0.995:
        bands.use_fallback_175 = False
    else:
        # union-of-edges envelope over EVERY grid point (round 29, stats M4:
        # |rho(1.75)| is an unnormalized magnitude riding the level axis, and
        # "worst sign" is undefined for a two-sided band -- both edges take
        # their own envelope over the full coverage grid)
        los, his = [], []
        for key, bank in grid_and_corners.items():
            v = np.abs(bank[: len(bank) // 2, i175])
            lo, lo_se = _q_boot(v, 0.0013, rng, cfg.n_boot)
            hi, hi_se = _q_boot(v, 0.9987, rng, cfg.n_boot)
            los.append(lo - lo_se)
            his.append(hi + hi_se)
        bands.fallback_175 = (min(los), max(his))
        bands.use_fallback_175 = True
    rep["fallback_engaged"] = bands.use_fallback_175

    bands.b_sym = max(_q_boot(col(banks[s], "rsym"), 0.9987, rng, cfg.n_boot)[0]
                      for s in cfg.signs)
    return bands, rep, {"banks_ho": banks_ho, "grid": grid, "inc": inc,
                        "lam_worst": lam_worst, "corner_ho": corner_ho}


def gate_checks(bands, rep, aux, cfg):
    """The gate's measured controls on HELD-OUT material: joint power, aggregate
    false-VOID, negative controls (null + impostors, sign-swept), LR
    separation, the s_j0 floor sanity pins."""
    from dataclasses import replace
    checks = {}

    def breakdown(bank):
        """First-failure reason per row (the diagnosis the corner needs)."""
        reasons = {}
        for row in bank:
            v, det = record_verdict(row, bands)
            if v == "CONFIRMED":
                key = "CONFIRMED"
            elif v == "VOID":
                key = "VOID:" + ",".join(det["void"])
            elif not det.get("a"):
                key = "a_fail"
            elif not det.get("b"):
                key = "b_fail"
            elif det.get("c_sup", 0) >= bands.b_curve:
                key = "c_sup_fail"
            elif det.get("c_contrast", 99) <= 5:
                key = "c_contrast_fail"
            else:
                key = v
            reasons[key] = reasons.get(key, 0) + 1
        return reasons

    # joint power AND per-grid-point false-VOID (round 29, spec m5: a corner
    # false-VOID is exactly what round 28 found; the aggregate is measured at
    # EVERY coverage point, the corner maximum reported) on the measure halves
    power, fv_grid = {}, {}
    for key, bank in aux["grid"].items():
        half = bank[len(bank) // 2:]
        verds = [record_verdict(row, bands)[0] for row in half]
        power[str(key)] = sum(v == "CONFIRMED" for v in verds) / len(half)
        fv_grid[str(key)] = sum(v == "VOID" for v in verds) / len(half)
    checks["joint_power_min"] = min(power.values())
    checks["joint_power_corner"] = min(power, key=power.get)
    checks["joint_power"] = power
    checks["false_void_grid_max"] = max(fv_grid.values())
    checks["false_void_grid_corner"] = max(fv_grid, key=fv_grid.get)
    corner_key = next(k for k in aux["grid"] if str(k) == checks["joint_power_corner"])
    corner_bank = aux["grid"][corner_key]
    checks["corner_breakdown"] = breakdown(corner_bank[len(corner_bank) // 2:])
    # aggregate false-VOID on held-out worst-basis valid banks, per sign
    fv = {}
    for s, bank in aux["banks_ho"].items():
        n_void = sum(1 for row in bank if record_verdict(row, bands)[0] == "VOID")
        fv[str(s)] = n_void / len(bank)
    checks["false_void"] = fv
    checks["false_void_max"] = max(fv.values())
    # the fixed-j'-corner devices enter the GATES (round 30, stats M2: the
    # round-29 fix reached the bands but not the gates): power and false-VOID
    # measured on the held-out corner twins (device-fixed nuisance = a corner,
    # only shot noise drawn); the near-degenerate-valid certification of the
    # loosened s_j0 floor construction rides here too (round 30, stats m13:
    # certify against a VALID device near the floor, not only the null)
    cp, cfv = {}, {}
    for s, bank in aux["corner_ho"].items():
        verds = [record_verdict(row, bands)[0] for row in bank]
        cp[str(s)] = sum(v == "CONFIRMED" for v in verds) / len(bank)
        cfv[str(s)] = sum(v == "VOID" for v in verds) / len(bank)
    checks["corner_fixed_jp_power_min"] = min(cp.values())
    checks["corner_fixed_jp_false_void_max"] = max(cfv.values())
    # negative controls: null + 3 LR impostors + impz, swept over the SIGN axis
    # AND the LEVEL ends {worst-admitted, ideal} (round 30, spec M6: the
    # envelopes were widened along level/slope, so the model-independent
    # false-CONFIRMED protection sweeps the level axis too; slope's effect on a
    # non-law device is covered by the level ends, pinned); impz additionally
    # carries its own DETECTION measurement (round 30, spec M2/M3: the z_k
    # guard bears the whole discrimination for that class, so its detection
    # power is measured, never assumed)
    neg = {}
    impz_det = []
    # control CONFIGS (round 31, spec B1: the round-30 fold widened the verdict
    # bands over the fixed-j'/box/eta-top corner sources while the controls
    # swept only sign x level ends -- the <=0.1% UCB must certify the SAME band
    # set that is committed, so the corner configs join the control sweep)
    def _corner_cfgs(lam_worst):
        base_w = _signed_base(cfg, lam_worst)
        pj = replace(base_w)
        pj.readout = dict(pj.readout)
        e01, e10 = pj.readout[Q_JP]
        pj.readout[Q_JP] = (e01 * (1 + cfg.spread), e10 * (1 + cfg.spread))
        pj.t2 = dict(pj.t2)
        pj.t2[Q_JP] /= (1 + cfg.spread)
        pj.epc_write_sjp = pj.epc_write * (1 + cfg.spread)
        pb = replace(_signed_base(cfg, 0.0))
        pb.readout = {q: (a * (1 + cfg.spread), b * (1 + cfg.spread))
                      for q, (a, b) in pb.readout.items()}
        pt = replace(base_w, t2=dict(base_w.t2))
        pt.t2[Q_JP] = min(pt.t2.get(Q_JP, 70e-6), 70e-6) / (1 + cfg.spread)
        return {"jpcorner": pj, "box": pb, "etatop_null": pt}

    corner_cfgs = _corner_cfgs(aux["lam_worst"])
    for kind in ("null", "imp1", "imp2", "imp3", "impz"):
        neg[kind] = {}
        for lam_lbl, base_p in ([("worst", _signed_base(cfg, aux["lam_worst"])),
                                 ("ideal", _signed_base(cfg, 0.0))] +
                                list(corner_cfgs.items())):
            lbl_off = {"worst": 0, "ideal": 77, "jpcorner": 154,
                       "box": 231, "etatop_null": 308}[lam_lbl]
            # corner configs are device-FIXED: spread 0 (only shot noise drawn)
            spread_n = 0.0 if lam_lbl in ("jpcorner", "box") else cfg.spread
            for s in cfg.signs:
                p_n = replace(base_p, bias_delta=s, arm_bias_delta=s)
                bank = run_bank(kind, p_n, cfg.n_null,
                                cfg.seed + 555 + int(s * 1e4) + lbl_off,
                                spread_n, cfg.jobs)
                events = sum(1 for row in bank if record_verdict(row, bands)[0] == "CONFIRMED")
                try:
                    from scipy.stats import chi2
                    ucb = chi2.ppf(0.95, 2 * events + 2) / 2.0 / len(bank)
                except Exception:
                    ucb = (events + 3.0) / len(bank)
                neg[kind][f"{lam_lbl}/{s}"] = {"events": events, "n": len(bank), "ucb95": ucb}
                if kind == "impz":
                    impz_det.append(float(np.mean(bank[:, ICOL["z_k"]] > bands.guard_z_k)))
    checks["negative_controls"] = neg
    checks["neg_events_total"] = sum(v["events"] for k in neg for v in neg[k].values())
    checks["impz_zk_detection_min"] = min(impz_det)
    # impz BOUNDARY control at m = the frozen z_k band value (round 30, stats
    # B1: the escape window m < band is QUANTIFIED, not assumed away -- at the
    # boundary the z_k statistic sits ON the band and P(CONFIRMED) is the
    # honest escape rate the pre-reg states; the certificate sentence is
    # scoped to the m-polarised branch, and the watcher-tomograph PUB carries
    # the m ~ 0 discrimination as a diagnostic). REPORT, not a <= 0.1% gate.
    p_bnd = replace(_signed_base(cfg, aux["lam_worst"]), impz_m=bands.guard_z_k)
    bnd_bank = run_bank("impz", p_bnd, max(cfg.n_null // 10, 500),
                        cfg.seed + 20_260, cfg.spread, cfg.jobs)
    bnd_verds = [record_verdict(row, bands)[0] for row in bnd_bank]
    checks["impz_boundary"] = {
        "m": bands.guard_z_k,
        "p_confirmed": sum(v == "CONFIRMED" for v in bnd_verds) / len(bnd_verds),
        "p_void": sum(v == "VOID" for v in bnd_verds) / len(bnd_verds),
        "mechanism_claim_rate": float(np.mean(bnd_bank[:, ICOL["tomo"]] >=
                                              bands.tomo_threshold)),
    }
    # watcher-tomograph mechanism separation (round 30, physics BLOCKER): the
    # signal must claim the mechanism, the Z-diagonal class must never
    impz_pure = run_bank("impz", replace(_signed_base(cfg, aux["lam_worst"]), impz_m=0.0),
                         max(cfg.n_null // 10, 500), cfg.seed + 20_261,
                         cfg.spread, cfg.jobs)
    sig_tomo = np.concatenate([b[:, ICOL["tomo"]] for b in aux["banks_ho"].values()])
    checks["tomograph"] = {
        "threshold": bands.tomo_threshold,
        "signal_mean": float(np.mean(sig_tomo)),
        "signal_claim_rate": float(np.mean(sig_tomo >= bands.tomo_threshold)),
        "impz_m0_mean": float(np.mean(impz_pure[:, ICOL["tomo"]])),
        "impz_m0_claim_rate": float(np.mean(impz_pure[:, ICOL["tomo"]] >=
                                            bands.tomo_threshold)),
    }
    # LR separation >= 10 sigma per impostor (on its own zero-sign bank)
    sep = {}
    for kind in ("imp1", "imp2", "imp3"):
        p_n = _signed_base(cfg, aux["lam_worst"])
        bank = run_bank(kind, p_n, max(cfg.n_null // 10, 200), cfg.seed + 888,
                        cfg.spread, cfg.jobs)
        vals = _lr_rows(bank, bands)[kind]
        sep[kind] = float((vals.mean() - bands.lr_thresholds[kind]) / vals.std())
    checks["lr_separation_sigma"] = sep
    # guard DETECTION power (round 29, spec m6 + stats m2: the T-hat band's
    # pinned duty is two-sided -- cover the 0.05 budget AND trip beyond it;
    # probe at the 0.2 rad/arm power-cliff rate the document names)
    p_probe = replace(_signed_base(cfg, aux["lam_worst"]),
                      rz_drift={Q_J: 0.2, Q_JP: 0.2})
    probe = run_bank("signal", p_probe, max(cfg.n_null // 100, 300),
                     cfg.seed + 31337, cfg.spread, cfg.jobs)
    trips = np.mean((probe[:, ICOL["t_j"]] > bands.guard_t_sup) |
                    (probe[:, ICOL["t_jp"]] > bands.guard_t_sup))
    checks["t_probe_trip_rate_at_0p2"] = float(trips)
    # s_j0 floor sanity: >=10 SHOT-SE below the worst-sign valid mean, >=10
    # above the null. Round-25 pin: the sanity is in SHOT-SE units of S_hat(j;0)
    # (shot + CAL noise at a FIXED device), NOT the systematics-broadened bank
    # SD, which the +/-spread device draw inflates severalfold: a spread-0 bank
    # isolates the right scale.
    s_means = {s: float(np.mean(b[:, ICOL["s_j0"]])) for s, b in aux["banks_ho"].items()}
    worst_sign = min(s_means, key=s_means.get)
    worst_mean = s_means[worst_sign]
    cons = rep.get("s_j0_floor_construction")
    p_fix = replace(_signed_base(cfg, aux["lam_worst"]),
                    bias_delta=worst_sign, arm_bias_delta=worst_sign)
    fix_bank = run_bank("signal", p_fix, max(cfg.n_null // 50, 500),
                        cfg.seed + 1234, 0.0, cfg.jobs)
    shot_se = float(np.std(fix_bank[:, ICOL["s_j0"]]))
    # SE-ordering REPORT, not an assert (round 29, physics m2 -- and the
    # reviewer's "forced > 1" claim REFUTED from below at the double-ratio
    # level: the branch-variance ordering 1 vs 1-eta^2 holds for the NUMERATOR
    # records, but the double ratio weights denominator/normalizer noise by
    # the VALUE, so rho(2) ~ -1 carries three full-weight noise factors while
    # rho(1) ~ 0 carries essentially one; measured fixed-device ratio ~0.88,
    # systematics-broadened ~0.75-0.9 -- the reviewer's pinned-figure
    # consistency question stands and is answered by the gate's own frozen
    # values, never by the naive ordering)
    checks["se_ratio_rho1_rho2"] = float(np.std(fix_bank[:, R_LR.index(1.0)]) /
                                         np.std(fix_bank[:, R_LR.index(2.0)]))
    # SE-ratio PROVENANCE (round 30, physics MAJOR: the v31 "REFUTED, measured
    # 0.88" wording over-claimed -- a shot-noise-only model reads 1.20-1.26 at
    # admitted bases, and the 0.88 came from the FULL chain incl. confusion +
    # CAL + mitigation; the ratio is CONFIGURATION-DEPENDENT, not structural
    # either way; isolate by re-running the fixed device without the readout
    # layer): same config, readout/bias cleared, CAL trivial
    p_iso = replace(p_fix, readout={}, bias_delta=0.0, arm_bias_delta=0.0,
                    bias_sym=0.0)
    iso_bank = run_bank("signal", p_iso, max(cfg.n_null // 50, 500),
                        cfg.seed + 1235, 0.0, cfg.jobs)
    checks["se_ratio_no_readout_chain"] = float(np.std(iso_bank[:, R_LR.index(1.0)]) /
                                                np.std(iso_bank[:, R_LR.index(2.0)]))
    p_null = _signed_base(cfg, aux["lam_worst"])
    null_bank = run_bank("null", p_null, max(cfg.n_null // 10, 200), cfg.seed + 999,
                         cfg.spread, cfg.jobs)
    null_mean = float(np.mean(null_bank[:, ICOL["s_j0"]]))
    # the below-valid half uses THE CONSTRUCTION'S OWN (mean, shot-SE) pair
    # (round 29: check and construction briefly used different mean sets and
    # SE estimates and disagreed by a hair; one pin, one definition)
    m_ref = cons["worst_valid_mean"] if cons else worst_mean
    se_ref = cons["shot_se"] if cons else shot_se
    checks["s_j0_floor_sanity"] = {
        "floor": bands.floor_s_j0, "valid_worst_mean": m_ref,
        "shot_se_fixed_device": se_ref, "null_mean": null_mean,
        "worst_basis_mean": worst_mean, "fix_bank_shot_se": shot_se,
        "sigmas_below_valid": (m_ref - bands.floor_s_j0) / se_ref,
        "sigmas_above_null": (bands.floor_s_j0 - null_mean) / se_ref,
    }
    return checks


def run_freeze_validation(scale="validate"):
    """End-to-end machinery validation at reduced scale (its numbers never
    fly); 'full' runs the pinned scales (hours; use jobs)."""
    import json, time
    if scale == "full":
        cfg = FreezeConfig(
            jobs=max(1, (__import__("os").cpu_count() or 2) - 2),
            aer_report=(r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI"
                        r"\experiments\ibm_quantum_tomography\_record_parity_aer_report.json"))
    else:
        cfg = FreezeConfig(n_phys=3000, n_guard=3000, n_power=600, n_null=1000,
                           levels=3, slopes=3, n_boot=100)
    t0 = time.perf_counter()
    bands, rep, aux = freeze_bands(cfg)
    checks = gate_checks(bands, rep, aux, cfg)
    dt = time.perf_counter() - t0
    out = {"config": {k: (list(v) if isinstance(v, tuple) else v)
                      for k, v in cfg.__dict__.items()},
           "report": rep, "checks": {k: v for k, v in checks.items() if k != "joint_power"},
           "joint_power": checks["joint_power"], "runtime_s": dt,
           "bands": {"b_blind": bands.b_blind, "b_forgive": bands.b_forgive,
                     "b_curve": bands.b_curve, "se_contrast": bands.se_contrast,
                     "b_sym": bands.b_sym,
                     "floors": [bands.floor_s_j0, bands.floor_d_jp0, bands.floor_d_jp_min],
                     "guard_t_sup": bands.guard_t_sup,
                     "guards_z": [bands.guard_z_k, bands.guard_z_s, bands.guard_z_leak],
                     "vs_threshold": bands.vs_threshold,
                     "tomo_threshold": bands.tomo_threshold,
                     "lr_thresholds": dict(bands.lr_thresholds),
                     "fallback_engaged": bands.use_fallback_175}}
    path = "simulations/_record_parity_freeze_%s.json" % scale
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f"[{scale}] runtime {dt:.0f}s; report -> {path}")
    print(f"  b_blind={bands.b_blind:.4f}  b_forgive={bands.b_forgive:.4f}  "
          f"b_curve={bands.b_curve:.4f}  se_contrast={bands.se_contrast:.4f}")
    print(f"  floors s_j0/d_jp0/d_jp_min = {bands.floor_s_j0:.4f}/{bands.floor_d_jp0:.4f}/"
          f"{bands.floor_d_jp_min:.4f}  t_sup={bands.guard_t_sup:.4f}")
    print(f"  joint power min = {checks['joint_power_min']:.4f} at {checks['joint_power_corner']}")
    print(f"  corner breakdown: {checks['corner_breakdown']}")
    print(f"  false-VOID max = {checks['false_void_max']:.4f} (worst basis), "
          f"{checks['false_void_grid_max']:.4f} (grid corner {checks['false_void_grid_corner']})  "
          f"negative-control CONFIRMED events total = {checks['neg_events_total']}")
    print(f"  T-hat detection probe at 0.2 rad/arm: trip rate = "
          f"{checks['t_probe_trip_rate_at_0p2']:.4f}   SE(rho1)/SE(rho2) = "
          f"{checks['se_ratio_rho1_rho2']:.2f} (full chain) / "
          f"{checks['se_ratio_no_readout_chain']:.2f} (no readout layer)")
    print(f"  fixed-j' corner gates: power min = {checks['corner_fixed_jp_power_min']:.4f}, "
          f"false-VOID max = {checks['corner_fixed_jp_false_void_max']:.4f}")
    tg = checks["tomograph"]
    print(f"  tomograph: threshold {tg['threshold']:.4f}, signal claim "
          f"{tg['signal_claim_rate']:.4f} (mean {tg['signal_mean']:.4f}), impz(m=0) claim "
          f"{tg['impz_m0_claim_rate']:.4f} (mean {tg['impz_m0_mean']:.4f})")
    ib = checks["impz_boundary"]
    print(f"  impz boundary (m = z_k band {ib['m']:.4f}): P(CONFIRMED) = "
          f"{ib['p_confirmed']:.3f}, P(VOID) = {ib['p_void']:.3f}, mechanism claim = "
          f"{ib['mechanism_claim_rate']:.3f}")
    print(f"  LR separation (sigma) = " +
          ", ".join(f"{k}:{v:.0f}" for k, v in checks["lr_separation_sigma"].items()))
    print(f"  order-check min power = {rep['order_check_min_power']:.4f}  "
          f"fallback engaged = {rep['fallback_engaged']}")
    print(f"  s_j0 floor sanity: {checks['s_j0_floor_sanity']['sigmas_below_valid']:.0f} SE "
          f"below valid, {checks['s_j0_floor_sanity']['sigmas_above_null']:.0f} SE above null")
    # wiring assertions (validation-scale smoke; the FULL run enforces the pinned gates)
    assert checks["neg_events_total"] == 0, "negative controls must show zero CONFIRMED events"
    assert checks["joint_power_min"] > 0.8, f"joint power collapsed: {checks['joint_power_min']}"
    assert checks["false_void_max"] < 0.05, f"false-VOID out of range: {checks['false_void_max']}"
    assert all(v > 10 for v in checks["lr_separation_sigma"].values()), "LR separation < 10 sigma"
    fs = checks["s_j0_floor_sanity"]
    assert fs["sigmas_below_valid"] >= 10 - 1e-6 and fs["sigmas_above_null"] >= 10, \
        f"s_j0 floor sanity (shot-SE units, round-25 pin): {fs}"
    fv_limit = 0.005 if scale == "full" else 0.05     # the corner maximum GATES (round 30, spec M4)
    assert checks["false_void_grid_max"] < fv_limit, \
        f"grid false-VOID corner: {checks['false_void_grid_max']} at {checks['false_void_grid_corner']}"
    assert checks["t_probe_trip_rate_at_0p2"] >= 0.99, \
        f"T-hat detection probe at 0.2 rad/arm: {checks['t_probe_trip_rate_at_0p2']}"
    assert checks["impz_zk_detection_min"] >= 0.999, \
        f"impz z_k detection (the superposition certificate's measured duty): {checks['impz_zk_detection_min']}"
    assert checks["corner_fixed_jp_power_min"] > 0.8, \
        f"fixed-j' corner power collapsed: {checks['corner_fixed_jp_power_min']}"
    assert checks["corner_fixed_jp_false_void_max"] < fv_limit, \
        f"fixed-j' corner false-VOID: {checks['corner_fixed_jp_false_void_max']}"
    assert checks["tomograph"]["signal_claim_rate"] >= 0.99, \
        f"tomograph must certify the signal's mechanism: {checks['tomograph']}"
    assert checks["tomograph"]["impz_m0_claim_rate"] <= 0.001, \
        f"tomograph must never certify a Z-diagonal watcher: {checks['tomograph']}"
    # SE ratio is a report (see gate_checks: the naive >1 ordering is refuted
    # from below for the double ratio); no assert
    print("freeze-validation PASS (machinery; frozen numbers come from the full-scale run).")
    return bands, rep, checks


def _provisional_bands(n_bank=100, seed=20260805):
    """WIRING-ONLY bands: eta_nom = 1 everywhere (the branch-2 positive-control
    configuration), coarse guard/floor values, LR machinery from a small ideal
    signal bank. These never fly; stage 4b freezes the real ones through the
    pinned procedure."""
    rng = np.random.default_rng(seed)
    ideal = Params()
    dists = {r: model_arm_dists("signal", r, ideal) for r in R_ARMS}
    bank = []
    for _ in range(n_bank):
        sd = {r: sampled_arm_data(r, rng, dist_override=dists[r]) for r in R_ARMS}
        est = estimator_from_shots(sd)
        bank.append([est["rho"][r] for r in R_LR])
    bank = np.array(bank)
    sigma = np.cov(bank.T)
    bands = Bands(
        b_blind=0.065, b_forgive=0.5, b_curve=0.07,
        se_contrast=float(np.std(np.abs(bank[:, -1]) - np.abs(bank[:, 3]))),
        eta_nom={r: 1.0 for r in R_ARMS},
        lr_sigma=sigma, lr_thresholds={},
        floor_s_j0=0.8, floor_d_jp0=0.8, floor_d_jp_min=0.5,
        guard_t_sup=0.05, guard_z_k=0.05, guard_z_s=0.05, guard_z_leak=0.05,
    )
    # wiring LR thresholds: signal-bank empirical max with margin (in-sample; wiring only)
    lr_vals = {k: [] for k in IMPOSTOR_CURVES}
    for row in bank:
        lrs = lr_statistics({r: row[i] for i, r in enumerate(R_LR)}, bands)
        for k, v in lrs.items():
            lr_vals[k].append(v)
    bands.lr_thresholds = {k: float(np.max(v) + 3 * np.std(v)) for k, v in lr_vals.items()}
    return bands


def _selftest_stage4a():
    bands = _provisional_bands()
    ideal = Params()
    rng = np.random.default_rng(20260806)

    # branch-2 positive control: ideal signal through the ACTUAL verdict code
    # against eta = 1 targets must return CONFIRMED (clause logic + signs)
    sd, zl, vs, tomo = full_pub_shot_data("signal", ideal, rng)
    est, guards, vstat, tstat = analyze_flight(sd, zl, vs, tomo)
    v, det = verdict(est, guards, vstat, bands, tomo_stat=tstat)
    assert v == "CONFIRMED", f"branch-2 positive control must CONFIRM: {v} {det}"

    # the watcher-tomograph reads <Y_k Z_j> = +1 at r = 1 on the ideal signal
    # (round 30 physics BLOCKER: the ONLY PUB that separates a superposed
    # watcher from Z-diagonal classical noise) and 0 for the impz class at
    # every m -- exact at distribution level
    td_sig = tomograph_model_dists("signal", ideal)
    assert abs(tomograph_stat(td_sig) - 1.0) < 1e-12, \
        f"ideal tomograph must read +1 at r=1: {tomograph_stat(td_sig)}"
    for m in (1.0, 0.5, 0.0):
        td_z = tomograph_model_dists("impz", Params(impz_m=m))
        assert abs(tomograph_stat(td_z)) < 1e-12, \
            f"Z-diagonal watcher must read 0 on the tomograph (m={m}): {tomograph_stat(td_z)}"

    # each named impostor must fail CONFIRMED and trip its LR
    for kind in ("imp1", "imp2", "imp3"):
        sd, zl, vs, tomo = full_pub_shot_data(kind, ideal, rng)
        est, guards, vstat, tstat = analyze_flight(sd, zl, vs, tomo)
        v, det = verdict(est, guards, vstat, bands, tomo_stat=tstat)
        assert v == "DEVICE-DEVIATION", f"{kind} must land DEVICE-DEVIATION: {v} {det}"
        assert kind in det["lr_tripped"], f"{kind} must trip its own LR: {det['lr_tripped']}"

    # the null must be routed to VOID by the S_hat(j;0) floor (the Cauchy
    # denominator never reaches a verdict clause)
    sd, zl, vs, tomo = full_pub_shot_data("null", ideal, rng)
    est, guards, vstat, tstat = analyze_flight(sd, zl, vs, tomo)
    v, det = verdict(est, guards, vstat, bands, tomo_stat=tstat)
    assert v == "VOID" and "floor_s_j0" in det["void"], f"null must VOID on the floor: {v} {det}"

    # over-budget inter-arm drift must trip the transverse guard bank
    p_drift = Params(rz_drift={Q_J: 0.2})
    sd, zl, vs, tomo = full_pub_shot_data("signal", p_drift, rng)
    est, guards, vstat, tstat = analyze_flight(sd, zl, vs, tomo)
    v, det = verdict(est, guards, vstat, bands, tomo_stat=tstat)
    assert v == "VOID" and any(t.startswith("t_sup") for t in det["void"]), \
        f"over-budget drift must VOID on T_hat: {v} {det}"

    print("stage-4a self-test PASS: verdict machinery wired: branch-2 positive control")
    print("CONFIRMED; tomograph +1 ideal / 0 Z-diagonal (m = 1, 0.5, 0); imp1/imp2/imp3")
    print("-> DEVICE-DEVIATION via their own LR; null -> VOID on the S_hat(j;0) floor;")
    print("over-budget drift -> VOID on the transverse guard.")


def _selftest_counts_path():
    """The counts fast path against the sequence reference: same estimator code,
    exactly the split rule's sampling distribution; moment agreement within MC
    error, and the throughput that makes the 1e5/1e6 banks feasible."""
    import time
    ideal = Params()
    rng = np.random.default_rng(20260807)
    n = 400
    t0 = time.perf_counter()
    bank_c = seed_bank_rho("signal", ideal, n, rng)
    dt = time.perf_counter() - t0
    dists = {r: model_arm_dists("signal", r, ideal) for r in R_ARMS}
    bank_s = np.empty((n, len(R_LR)))
    for i in range(n):
        sd = {r: sampled_arm_data(r, rng, dist_override=dists[r]) for r in R_ARMS}
        est = estimator_from_shots(sd)
        bank_s[i] = [est["rho"][r] for r in R_LR]
    for k, r in enumerate(R_LR):
        m_c, m_s = bank_c[:, k].mean(), bank_s[:, k].mean()
        se = np.hypot(bank_c[:, k].std(), bank_s[:, k].std()) / np.sqrt(n)
        assert abs(m_c - m_s) < 5 * se, f"counts-vs-sequence moment at r={r}: {m_c} vs {m_s} (5se={5*se:.2e})"
        law = np.cos(r * np.pi / 2.0)
        assert abs(m_c - law) < 5 * bank_c[:, k].std() / np.sqrt(n) + 1e-3, \
            f"counts-path mean off the law at r={r}: {m_c} vs {law}"
    rate = n / dt
    print(f"counts-path self-test PASS: moments agree with the sequence reference; "
          f"{rate:.0f} seeds/s (1e5 bank ~ {1e5 / rate / 60:.1f} min single-core).")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in ("freeze-validate", "freeze-full"):
        run_freeze_validation("validate" if sys.argv[1] == "freeze-validate" else "full")
    else:
        _selftest()
        _selftest_sampling()
        _selftest_stage3()
        _selftest_stage3_sampled()
        _selftest_stage4a()
        _selftest_counts_path()
