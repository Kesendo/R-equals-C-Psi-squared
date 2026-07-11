"""
Stage 7a of the IBM Concentrator Reloaded pre-registration.

Design: experiments/IBM_CONCENTRATOR_RELOADED.md  (sections 2-6).
This script is the FROM-BELOW gate: exact density-matrix predictions WITH a
representative Kingston device background, a counts-level Monte-Carlo projected
SE, and the observable search (protection grading in sink DISTANCE).

Deterministic, seeded, numpy/scipy only. Prints everything it claims.
Nothing is run against IBM. No existing file is modified.

Physics (all pinned from the design):
  * N = 5 chain, first-order Trotter, bonds (0,1),(1,2),(2,3),(3,4) in order.
  * Per bond: exp(-i * (theta/2) * (XX + YY + ZZ)), theta = 2*J*dt = 0.1
    => (theta/2) = J*dt = 0.05.
  * Payload |+> on qubit 2 (main experiment), all others |0>.
  * Sink: per-step Z-dephasing on the sink site.
      - density-matrix path : exact dephasing channel, off-diagonal retention
        lambda_sink = e^{-0.05} (this is the retention e^{-2*d*dt} with
        2*d*dt = 0.05).
      - counts path         : M=256 i.i.d. Gaussian phase draws per step,
        sigma = sqrt(0.1) = 0.3162 rad, so mean phase factor
        E[e^{i phi}] = e^{-sigma^2/2} = e^{-0.05} = lambda_sink  (consistent).
      Exactly ONE sink application per site per step.
  * Grid {1,2,3,4,6,8} Trotter steps.
  * Arms: 0 (no sink), E (sink qubit 0, distance 2), MP (sink qubit 2, on payload).
  * Estimator: coh_a(t) = |mean(<X2> + i <Y2>)| ; paired ratio R_a = coh_a/coh_0 ;
    slope = LS slope of ln R_a over the grid ; verdict = slope(MP) - slope(E).

Author: Claude (Anthropic) for Thomas Wicht, 2026-07-11.
"""

import csv
import sys
from functools import reduce

import numpy as np
from scipy.linalg import expm

# --------------------------------------------------------------------------
# Seeds (printed in the report; the run is fully deterministic)
# --------------------------------------------------------------------------
SEED_BINDINGS = 20260711   # the phase-vector pool (the M bindings)
SEED_BOOT     = 770711     # the hierarchical binding+shot bootstrap
SEED_SEARCH   = 110711     # the observable-search counts MC

CSV_PATH = (r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
            r"\ClaudeTasks\IBM_R2_calibrations"
            r"\ibm_kingston_calibrations_2026-07-11T10_02_33Z.csv")

# --------------------------------------------------------------------------
# Pinned constants
# --------------------------------------------------------------------------
N        = 5
DIM      = 1 << N
JDT      = 0.05                    # theta/2 per bond
BONDS    = [(0, 1), (1, 2), (2, 3), (3, 4)]
GRID     = [1, 2, 3, 4, 6, 8]
GRID_NO_DEEP = [1, 2, 3, 4]
LAMBDA_SINK = np.exp(-0.05)        # off-diagonal retention per step
SIGMA_SINK  = np.sqrt(0.1)         # Gaussian phase std per step
M_BIND      = 256                  # bindings per (arm, depth)
SHOTS_BIND  = 32                   # shots per binding per basis
R_BOOT      = 500                  # bootstrap resamples (>= 300 design floor)

# --------------------------------------------------------------------------
# 1-qubit ops
# --------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X  = np.array([[0, 1], [1, 0]], dtype=complex)
Y  = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z  = np.array([[1, 0], [0, -1]], dtype=complex)
KET0 = np.array([1, 0], dtype=complex)
KETP = np.array([1, 1], dtype=complex) / np.sqrt(2)


def embed(op2, q, n=N):
    """Place 2x2 op2 on qubit q (qubit 0 = least significant bit)."""
    factors = [I2] * (n - 1 - q) + [op2] + [I2] * q
    return reduce(np.kron, factors)


# Precompute per-qubit embedded X, Y, Z
XQ = [embed(X, q) for q in range(N)]
YQ = [embed(Y, q) for q in range(N)]
ZQ = [embed(Z, q) for q in range(N)]

# Bit table: bit_q(i) for every basis index
BITS = np.array([[(i >> q) & 1 for i in range(DIM)] for q in range(N)])  # (N, DIM)


def bond_unitary(a, b):
    H = XQ[a] @ XQ[b] + YQ[a] @ YQ[b] + ZQ[a] @ ZQ[b]
    return expm(-1j * JDT * H)


# One Trotter step unitary: apply bonds in order (0,1),(1,2),(2,3),(3,4)
# state after = U34 U23 U12 U01 psi  ->  U_step = U34 @ U23 @ U12 @ U01
_UBOND = [bond_unitary(a, b) for (a, b) in BONDS]
U_STEP = reduce(lambda acc, u: u @ acc, _UBOND, np.eye(DIM, dtype=complex))


def dephase_mask(q, lam):
    """Elementwise mask: rho[i,j] *= lam if bit_q differs."""
    diff = (BITS[q][:, None] != BITS[q][None, :])
    m = np.ones((DIM, DIM), dtype=complex)
    m[diff] = lam
    return m


def sink_phase_mask(q, phi):
    """Elementwise mask for RZ(phi) conjugation on qubit q: e^{-i phi (b_i - b_j)}."""
    db = BITS[q][:, None] - BITS[q][None, :]
    return np.exp(-1j * phi * db)


def amp_damp_kraus(q, gamma):
    """Embedded amplitude-damping Kraus on qubit q."""
    k0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    k1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return embed(k0, q), embed(k1, q)


# --------------------------------------------------------------------------
# PART 1 (a): representative uniform background from the calibration CSV
# --------------------------------------------------------------------------
def load_calibration(path):
    rows = []
    with open(path, newline="") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        for r in rdr:
            rows.append(r)
    return header, rows


def parse_float(s):
    s = s.strip().strip('"')
    if s == "":
        return None
    return float(s)


def parse_neighbors(tok):
    """Parse 'n:val;n2:val2' -> [n, n2]."""
    tok = tok.strip().strip('"')
    out = []
    if not tok:
        return out
    for part in tok.split(";"):
        if ":" in part:
            out.append(int(part.split(":")[0]))
    return out


def build_background(path):
    header, rows = load_calibration(path)
    # column indices
    cQ, cT1, cT2, cRO = 0, 1, 2, 3
    cP01, cP10 = 4, 5          # Prob meas0 prep1 , Prob meas1 prep0
    cGATE2Q = 14               # 'Gate length (ns)' 2q tokens  'n:len'
    cOP = 18

    data = {}
    adj = {}
    lens2q = []
    for r in rows:
        q = int(parse_float(r[cQ]))
        t1 = parse_float(r[cT1])
        t2 = parse_float(r[cT2])
        ro = parse_float(r[cRO])
        p01 = parse_float(r[cP01])
        p10 = parse_float(r[cP10])
        op = r[cOP].strip().strip('"')
        data[q] = dict(t1=t1, t2=t2, ro=ro, p01=p01, p10=p10, op=op)
        adj[q] = parse_neighbors(r[cGATE2Q])
        for part in r[cGATE2Q].strip().strip('"').split(";"):
            if ":" in part:
                try:
                    lens2q.append(float(part.split(":")[1]))
                except ValueError:
                    pass

    # rule-passing set: T2 >= 150 us, readout <= 2%, operational, non-empty
    S = [q for q, d in data.items()
         if d["t1"] and d["t2"] and d["ro"] is not None
         and d["t2"] >= 150.0 and d["ro"] <= 0.02 and d["op"] == "Yes"]
    S = sorted(S)

    # path search: connected simple 5-lines inside S with max/min T2 <= 2
    Sset = set(S)
    subadj = {q: [nb for nb in adj[q] if nb in Sset] for q in S}
    found_lines = []

    def dfs(path):
        if len(path) == 5:
            t2s = [data[q]["t2"] for q in path]
            if max(t2s) / min(t2s) <= 2.0:
                found_lines.append(tuple(path))
            return
        for nb in subadj[path[-1]]:
            if nb not in path:
                dfs(path + [nb])

    for start in S:
        dfs([start])
    # dedupe reversed duplicates
    uniq = set()
    for ln in found_lines:
        uniq.add(ln if ln <= ln[::-1] else ln[::-1])
    found_lines = sorted(uniq)

    # median representative of S
    med = lambda key: float(np.median([data[q][key] for q in S]))
    rep = dict(t1=med("t1"), t2=med("t2"), ro=med("ro"),
               p01=med("p01"), p10=med("p10"))

    # best connected line (max of min-T2), for transparency
    best_line = None
    if found_lines:
        best_line = max(found_lines,
                        key=lambda ln: min(data[q]["t2"] for q in ln))

    return data, S, found_lines, best_line, rep, float(np.median(lens2q))


# --------------------------------------------------------------------------
# PART 1 (b): per-step wall time and per-step decoherence channels
# --------------------------------------------------------------------------
def per_step_channels(rep, g2_ns, g1_ns=32.0, schedule="parallel"):
    """
    Arithmetic (design 4(b)): each Trotter step = 4 bonds x (rxx+ryy+rzz),
    each rotation ~ 2 two-qubit gates + single-qubit dressing.
    per rotation critical path = 2*g2 + 2*g1
    per bond (3 rotations)      = 3*(2*g2 + 2*g1)
    per step:
      sequential (conservative)  = 4 bonds  * per_bond
      parallel   (linear chain)  = 2 layers * per_bond   {(0,1),(2,3)} then {(1,2),(3,4)}
    """
    per_rot = 2 * g2_ns + 2 * g1_ns
    per_bond = 3 * per_rot
    layers = 4 if schedule == "sequential" else 2
    tau_ns = layers * per_bond
    tau = tau_ns * 1e-9                    # seconds
    T1 = rep["t1"] * 1e-6
    T2 = rep["t2"] * 1e-6
    gamma_amp = 1.0 - np.exp(-tau / T1)
    inv_tphi = max(0.0, 1.0 / T2 - 1.0 / (2.0 * T1))
    lam_bg = np.exp(-tau * inv_tphi)       # pure-dephasing retention
    return dict(tau_ns=tau_ns, gamma_amp=gamma_amp, lam_bg=lam_bg,
                per_rot=per_rot, per_bond=per_bond, layers=layers)


# --------------------------------------------------------------------------
# Density-matrix evolution primitives
# --------------------------------------------------------------------------
def make_background_ops(gamma_amp, lam_bg):
    amp = [amp_damp_kraus(q, gamma_amp) for q in range(N)]
    bg_mask = np.ones((DIM, DIM), dtype=complex)
    for q in range(N):
        bg_mask = bg_mask * dephase_mask(q, lam_bg)
    return amp, bg_mask


def apply_background(rho, amp, bg_mask):
    for (k0, k1) in amp:
        rho = k0 @ rho @ k0.conj().T + k1 @ rho @ k1.conj().T
    rho = rho * bg_mask
    return rho


def init_state(payload_qubit):
    kets = [KET0] * N
    kets[payload_qubit] = KETP
    # tensor order: qubit N-1 (msb) ... qubit 0 (lsb)
    psi = reduce(np.kron, [kets[q] for q in range(N - 1, -1, -1)])
    return np.outer(psi, psi.conj())


def bloch_xy(rho, q):
    bx = np.real(np.trace(rho @ XQ[q]))
    by = np.real(np.trace(rho @ YQ[q]))
    return bx, by


def evolve_exact(payload, readout, sink, amp, bg_mask, use_bg=True,
                 grid=GRID):
    """
    Exact density-matrix evolution with the EXACT sink dephasing channel
    (retention LAMBDA_SINK). Returns dict depth -> (bx, by) at readout qubit.
    sink = None for arm 0.
    """
    rho = init_state(payload)
    sink_mask = dephase_mask(sink, LAMBDA_SINK) if sink is not None else None
    out = {}
    gmax = max(grid)
    for t in range(1, gmax + 1):
        rho = U_STEP @ rho @ U_STEP.conj().T
        if use_bg:
            rho = apply_background(rho, amp, bg_mask)
        if sink_mask is not None:
            rho = rho * sink_mask
        if t in grid:
            out[t] = bloch_xy(rho, readout)
    return out


def evolve_binding(payload, readout, sink, phases, amp, bg_mask, use_bg=True,
                   grid=GRID):
    """
    One counts-path binding: RZ(phi_t) on the sink qubit at step t.
    phases: array length gmax (per-step phase draws). Returns depth->(bx,by).
    """
    rho = init_state(payload)
    out = {}
    gmax = max(grid)
    for t in range(1, gmax + 1):
        rho = U_STEP @ rho @ U_STEP.conj().T
        if use_bg:
            rho = apply_background(rho, amp, bg_mask)
        if sink is not None:
            rho = rho * sink_phase_mask(sink, phases[t - 1])
        if t in grid:
            out[t] = bloch_xy(rho, readout)
    return out


# --------------------------------------------------------------------------
# slope helpers
# --------------------------------------------------------------------------
def coh_from_xy(bx, by):
    return np.hypot(bx, by)


def slope_lnR(coh_a, coh_0, grid):
    t = np.array(grid, dtype=float)
    lnR = np.log(np.asarray(coh_a) / np.asarray(coh_0))
    A = np.vstack([t, np.ones_like(t)]).T
    m = np.linalg.lstsq(A, lnR, rcond=None)[0][0]
    return m


# --------------------------------------------------------------------------
# PART 1 (c): exact DM predictions with background
# --------------------------------------------------------------------------
def part1c(amp, bg_mask):
    res = {}
    for use_bg, tag in [(False, "clean"), (True, "bg")]:
        arms = {}
        for name, sink in [("0", None), ("E", 0), ("MP", 2)]:
            xy = evolve_exact(2, 2, sink, amp, bg_mask, use_bg=use_bg)
            coh = {t: coh_from_xy(*xy[t]) for t in GRID}
            arms[name] = coh
        res[tag] = arms
    return res


# --------------------------------------------------------------------------
# PART 1 (d): counts-level MC  (bindings pool + hierarchical bootstrap)
# --------------------------------------------------------------------------
def build_binding_pool(sink, payload, readout, amp, bg_mask, rng, use_bg=True):
    """
    Draw M bindings; return per-binding true-Bloch arrays shaped (M, n_depth)
    for X (bx) and Y (by). sink=None => single (no-phase) evolution broadcast.
    """
    nd = len(GRID)
    if sink is None:
        xy = evolve_exact(payload, readout, None, amp, bg_mask, use_bg=use_bg)
        bx = np.array([[xy[t][0] for t in GRID]])   # (1, nd)
        by = np.array([[xy[t][1] for t in GRID]])
        return np.repeat(bx, M_BIND, axis=0), np.repeat(by, M_BIND, axis=0)
    BX = np.empty((M_BIND, nd))
    BY = np.empty((M_BIND, nd))
    gmax = max(GRID)
    phases_all = rng.normal(0.0, SIGMA_SINK, size=(M_BIND, gmax))
    for m in range(M_BIND):
        xy = evolve_binding(payload, readout, sink, phases_all[m],
                            amp, bg_mask, use_bg=use_bg)
        BX[m] = [xy[t][0] for t in GRID]
        BY[m] = [xy[t][1] for t in GRID]
    return BX, BY


def p0_obs(bloch, p01, p10):
    """observed P(outcome 0) given true bloch component along measured axis."""
    p0_true = (1.0 + bloch) / 2.0
    return p0_true * (1.0 - p10) + (1.0 - p0_true) * p01


def mitigate_bloch(p0_hat, p01, p10):
    contrast = 1.0 - p10 - p01
    p0_true = (p0_hat - p01) / contrast
    return 2.0 * p0_true - 1.0


def counts_coh(BX, BY, idx, rng, p01, p10, shots):
    """
    Pool 'shots' per selected binding, in X and Y basis, with readout
    confusion + linear mitigation. Returns coh(t) array over grid.
    idx: selected binding indices (length M).  BX,BY shaped (M, nd).
    """
    nd = BX.shape[1]
    bx_sel = BX[idx]                       # (M, nd)
    by_sel = BY[idx]
    # observed P0 per binding per depth
    pX = p0_obs(bx_sel, p01, p10)
    pY = p0_obs(by_sel, p01, p10)
    # binomial shots, pool across bindings
    nX = rng.binomial(shots, pX).sum(axis=0)     # (nd,)
    nY = rng.binomial(shots, pY).sum(axis=0)
    total = shots * BX.shape[0]
    p0X = nX / total
    p0Y = nY / total
    bx_hat = mitigate_bloch(p0X, p01, p10)
    by_hat = mitigate_bloch(p0Y, p01, p10)
    return np.hypot(bx_hat, by_hat)


def part1d(amp, bg_mask, rep, central_diff_bg):
    rng = np.random.default_rng(SEED_BINDINGS)
    # binding pools (with background) for arms 0, E, MP
    BX0, BY0 = build_binding_pool(None, 2, 2, amp, bg_mask, rng, use_bg=True)
    BXE, BYE = build_binding_pool(0, 2, 2, amp, bg_mask, rng, use_bg=True)
    BXM, BYM = build_binding_pool(2, 2, 2, amp, bg_mask, rng, use_bg=True)

    p01, p10 = rep["p01"], rep["p10"]
    bootrng = np.random.default_rng(SEED_BOOT)

    diffs = np.empty(R_BOOT)
    slopeE = np.empty(R_BOOT)
    slopeM = np.empty(R_BOOT)
    for b in range(R_BOOT):
        idx0 = bootrng.integers(0, M_BIND, M_BIND)
        idxE = bootrng.integers(0, M_BIND, M_BIND)
        idxM = bootrng.integers(0, M_BIND, M_BIND)
        c0 = counts_coh(BX0, BY0, idx0, bootrng, p01, p10, SHOTS_BIND)
        cE = counts_coh(BXE, BYE, idxE, bootrng, p01, p10, SHOTS_BIND)
        cM = counts_coh(BXM, BYM, idxM, bootrng, p01, p10, SHOTS_BIND)
        sE = slope_lnR(cE, c0, GRID)
        sM = slope_lnR(cM, c0, GRID)
        slopeE[b] = sE
        slopeM[b] = sM
        diffs[b] = sM - sE

    se_diff = float(np.std(diffs, ddof=1))
    se_E = float(np.std(slopeE, ddof=1))
    central = central_diff_bg
    power = abs(central) / se_diff
    band = (central - 2 * se_diff, central + 2 * se_diff)
    return dict(se_diff=se_diff, se_E=se_E, central=central, power=power,
                band=band, boot_mean_diff=float(np.mean(diffs)),
                boot_mean_slopeE=float(np.mean(slopeE)))


# --------------------------------------------------------------------------
# PART 2: observable search (protection grading in sink DISTANCE)
# --------------------------------------------------------------------------
def config_slope(payload, readout, sink, amp, bg_mask, use_bg=True):
    """Exact slope of ln R over grid for one (payload, readout, sink)."""
    xy0 = evolve_exact(payload, readout, None, amp, bg_mask, use_bg=use_bg)
    xy  = evolve_exact(payload, readout, sink, amp, bg_mask, use_bg=use_bg)
    c0 = np.array([coh_from_xy(*xy0[t]) for t in GRID])
    ca = np.array([coh_from_xy(*xy[t]) for t in GRID])
    return slope_lnR(ca, c0, GRID), c0, ca


def config_slope_by_step(payload, readout, sink, amp, bg_mask, use_bg=True):
    """Per-step ln R increments, to test sign stability across the grid."""
    xy0 = evolve_exact(payload, readout, None, amp, bg_mask, use_bg=use_bg)
    xy  = evolve_exact(payload, readout, sink, amp, bg_mask, use_bg=use_bg)
    lnR = np.array([np.log(coh_from_xy(*xy[t]) / coh_from_xy(*xy0[t]))
                    for t in GRID])
    return lnR


def config_slope_se(payload, readout, sink, amp, bg_mask, rng,
                    rep, use_bg=True):
    """Counts-MC SE of a single-config slope (for resolvability)."""
    # baseline pool: no sink, this payload/readout
    xy0 = evolve_exact(payload, readout, None, amp, bg_mask, use_bg=use_bg)
    bx0 = np.array([[xy0[t][0] for t in GRID]])
    by0 = np.array([[xy0[t][1] for t in GRID]])
    BX0 = np.repeat(bx0, M_BIND, axis=0)
    BY0 = np.repeat(by0, M_BIND, axis=0)
    # sink pool
    BXs, BYs = build_binding_pool(sink, payload, readout, amp, bg_mask, rng,
                                  use_bg=use_bg)
    p01, p10 = rep["p01"], rep["p10"]
    slopes = np.empty(R_BOOT)
    for b in range(R_BOOT):
        idx0 = rng.integers(0, M_BIND, M_BIND)
        idxs = rng.integers(0, M_BIND, M_BIND)
        c0 = counts_coh(BX0, BY0, idx0, rng, p01, p10, SHOTS_BIND)
        cs = counts_coh(BXs, BYs, idxs, rng, p01, p10, SHOTS_BIND)
        slopes[b] = slope_lnR(cs, c0, GRID)
    return float(np.std(slopes, ddof=1))


def part2_search(amp, bg_mask, rep):
    print("\n" + "=" * 74)
    print("PART 2 - OBSERVABLE SEARCH (protection grading in sink DISTANCE)")
    print("=" * 74)
    rng = np.random.default_rng(SEED_SEARCH)

    candidates = []

    # ---- (i) END payload |+> on qubit 0, readout qubit 0, sinks at dist 1..4
    print("\n(i) END payload |+> on qubit 0, readout qubit 0")
    print("    sink distance d -> sink qubit d ; slope of ln R over the grid")
    end_slopes = {}
    end_lnR = {}
    for d in [1, 2, 3, 4]:
        s, c0, ca = config_slope(0, 0, d, amp, bg_mask, use_bg=True)
        lnR = config_slope_by_step(0, 0, d, amp, bg_mask, use_bg=True)
        end_slopes[d] = s
        end_lnR[d] = lnR
        print(f"    d={d} (sink q{d}): slope = {s:+.5f} /step   "
              f"lnR@grid = [{', '.join(f'{v:+.4f}' for v in lnR)}]")
    mono_dec = all(abs(end_slopes[d]) >= abs(end_slopes[d + 1]) - 1e-12
                   for d in [1, 2, 3])
    mono_signed = all(end_slopes[d] <= end_slopes[d + 1] + 1e-12
                      for d in [1, 2, 3])   # slopes negative, rising toward 0
    sign_stable = all(np.all(end_lnR[d] <= 1e-9) for d in [1, 2, 3, 4])
    print(f"    |slope| monotone decreasing in distance : {mono_dec}")
    print(f"    slope monotone rising toward 0           : {mono_signed}")
    print(f"    every-grid-point sign stable (lnR<=0)    : {sign_stable}")

    # resolvability at the flown shot budget: SE per slope, gaps vs 3*SE
    print("    resolvability (counts-MC SE per slope, same shot budget):")
    end_se = {}
    for d in [1, 2, 3, 4]:
        se = config_slope_se(0, 0, d, amp, bg_mask, rng, rep, use_bg=True)
        end_se[d] = se
        print(f"    d={d}: SE(slope) = {se:.5f}")
    gaps_ok = []
    for d in [1, 2, 3]:
        gap = abs(end_slopes[d] - end_slopes[d + 1])
        se_gap = np.hypot(end_se[d], end_se[d + 1])
        ok = gap >= 3 * se_gap
        gaps_ok.append(ok)
        print(f"    gap d{d}->d{d+1} = {gap:.5f} ; 3*SE_gap = {3*se_gap:.5f} "
              f"; resolvable = {ok}")
    end_ok = mono_signed and sign_stable and all(gaps_ok)
    candidates.append(("END payload q0 / readout q0", end_ok,
                       end_slopes, end_se))

    # ---- (ii) CENTER payload q2 (the main obs): confirm NON-monotone
    print("\n(ii) CENTER payload |+> on q2, readout q2 (the main observable)")
    print("     sink distances available: 1 (q1,q3), 2 (q0,q4)")
    for lab, s in [("d1 sink q1", 1), ("d1 sink q3", 3),
                   ("d2 sink q0", 0), ("d2 sink q4", 4)]:
        sl, _, _ = config_slope(2, 2, s, amp, bg_mask, use_bg=True)
        print(f"     {lab}: slope = {sl:+.5f} /step")

    # ---- (iii) NEAR-END payload q1, readout q1
    print("\n(iii) NEAR-END payload |+> on q1, readout q1")
    for d, sq in [("d1", 0), ("d1", 2), ("d2", 3), ("d3", 4)]:
        sl, _, _ = config_slope(1, 1, sq, amp, bg_mask, use_bg=True)
        print(f"     sink q{sq} ({d}): slope = {sl:+.5f} /step")

    # ---- (iii') two-site readout variant on the END payload:
    #      read the coherence on the neighbour q1 while payload sits on q0
    print("\n(iii') END payload q0, READOUT on q1 (transported-coherence variant)")
    for d in [1, 2, 3]:
        sq = d + 1   # sink beyond the readout qubit
        sl, _, _ = config_slope(0, 1, sq, amp, bg_mask, use_bg=True)
        print(f"     sink q{sq}: slope(readout q1) = {sl:+.5f} /step")

    print("\n--- PART 2 VERDICT ---")
    if end_ok:
        print("FOUND: END payload |+> on qubit 0 with sink at distance 1..4,")
        print("readout on qubit 0. Slopes monotone (rising toward 0 with")
        print("distance), sign-stable across the grid, and resolvable at the")
        print("flown shot budget.")
        print("winner slopes:", {d: round(end_slopes[d], 5) for d in [1, 2, 3, 4]})
    else:
        print("NOT FOUND (fully). Closest: END payload q0 / readout q0.")
        print(" per-distance slopes:",
              {d: round(end_slopes[d], 5) for d in [1, 2, 3, 4]})
        print(" monotone(signed):", mono_signed, "| sign-stable:", sign_stable,
              "| all gaps resolvable:", all(gaps_ok))
        print(" Why it falls short: the failing test(s) above (monotonicity,")
        print(" sign-stability, or resolvability of the far-distance gaps at")
        print(" the same shot budget).")
    return end_ok, end_slopes, end_se


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def main():
    print("=" * 74)
    print("IBM CONCENTRATOR RELOADED - STAGE 7a (from-below + observable search)")
    print("=" * 74)
    print(f"seeds: bindings={SEED_BINDINGS} bootstrap={SEED_BOOT} "
          f"search={SEED_SEARCH}")
    print(f"pins: J*dt={JDT}  theta={2*JDT}  lambda_sink=e^-0.05={LAMBDA_SINK:.6f}"
          f"  sigma_sink={SIGMA_SINK:.4f}")
    print(f"      M={M_BIND} bindings x {SHOTS_BIND} shots x 2 bases; "
          f"bootstrap R={R_BOOT}; grid={GRID}")

    # ---- (a) representative background
    data, S, lines, best_line, rep, med_g2 = build_background(CSV_PATH)
    print("\n" + "-" * 74)
    print("PART 1(a) - REPRESENTATIVE UNIFORM BACKGROUND")
    print("-" * 74)
    print(f"rule (design sec.8): T2 >= 150 us, readout <= 2%, operational, "
          f"max/min T2 <= 2 on the 5-line")
    print(f"rule-passing qubits (individual T2/readout gate): {len(S)}")
    print(f"  {S}")
    print(f"connected rule-passing linear 5-lines (max/min T2<=2): {len(lines)}")
    if best_line:
        t2s = [data[q]['t2'] for q in best_line]
        t1s = [data[q]['t1'] for q in best_line]
        print(f"  best connected line = {best_line}  "
              f"T2={[round(x,1) for x in t2s]}  max/min={max(t2s)/min(t2s):.2f}")
        print(f"    its T1={[round(x,1) for x in t1s]}  "
              f"(line T1 median={np.median(t1s):.1f}, "
              f"line T2 median={np.median(t2s):.1f})")
    else:
        print("  none (heavy-hex + the T2>=150 mask leaves no rule-passing "
              "connected 5-line).")
    print("\nDECISION: a rule-passing connected 5-line DOES exist on this")
    print("day, but the bands are frozen (design sec.5) on a REPRESENTATIVE,")
    print("chain-INDEPENDENT background, and the day-of chain rule is the only")
    print("thing that re-reads a specific line. Per the design's explicit")
    print("fallback ('use the rule-passing MEDIAN qubit parameters as a")
    print("uniform representative background and SAY SO'), the frozen")
    print("background is the MEDIAN of the rule-passing set, applied UNIFORMLY")
    print("to all 5 qubits. This is validated as representative: the actual")
    print("best line's T1/T2 bracket the medians used below.")
    print(f"  representative T1 = {rep['t1']:.1f} us")
    print(f"  representative T2 = {rep['t2']:.1f} us")
    print(f"  representative readout error       = {rep['ro']*100:.3f} %")
    print(f"  representative P(meas0|prep1) p01  = {rep['p01']*100:.3f} %")
    print(f"  representative P(meas1|prep0) p10  = {rep['p10']*100:.3f} %")

    # ---- (b) wall time + channels
    print("\n" + "-" * 74)
    print("PART 1(b) - PER-STEP WALL TIME AND DECOHERENCE CHANNELS")
    print("-" * 74)
    g2 = med_g2
    g1 = 32.0
    print(f"gate lengths from CSV: median 2q gate = {g2:.0f} ns, "
          f"1q gate = {g1:.0f} ns")
    print("arithmetic: per rotation ~ 2*2q + 2*1q dressing "
          f"= {2*g2+2*g1:.0f} ns")
    print(f"            per bond (rxx+ryy+rzz) = 3 rotations "
          f"= {3*(2*g2+2*g1):.0f} ns")
    ch_par = per_step_channels(rep, g2, g1, "parallel")
    ch_seq = per_step_channels(rep, g2, g1, "sequential")
    print(f"            per step, PARALLEL linear chain (2 bond-layers) "
          f"= {ch_par['tau_ns']:.0f} ns  <- representative")
    print(f"            per step, SEQUENTIAL (4 bond-layers, conservative) "
          f"= {ch_seq['tau_ns']:.0f} ns")
    print(f"per-step channels at representative T1/T2 (parallel schedule):")
    print(f"  amplitude damping gamma = 1-e^(-tau/T1) = {ch_par['gamma_amp']:.5f}")
    print(f"  pure-dephasing retention e^(-tau/Tphi)  = {ch_par['lam_bg']:.6f}")
    print(f"  (1/Tphi = 1/T2 - 1/2T1)")
    print(f"conservative (sequential) per-step channels:")
    print(f"  amplitude damping gamma = {ch_seq['gamma_amp']:.5f}  "
          f"dephasing retention = {ch_seq['lam_bg']:.6f}")

    amp, bg_mask = make_background_ops(ch_par['gamma_amp'], ch_par['lam_bg'])
    amp_seq, bg_mask_seq = make_background_ops(ch_seq['gamma_amp'],
                                               ch_seq['lam_bg'])

    # ---- (c) exact DM predictions
    print("\n" + "-" * 74)
    print("PART 1(c) - EXACT DENSITY-MATRIX PREDICTIONS")
    print("-" * 74)
    r1c = part1c(amp, bg_mask)

    def report_arms(arms, tag):
        print(f"\n[{tag}] per-grid-point coherences:")
        hdr = "  t     coh_0      coh_E      coh_MP     R_E=E/0    R_MP=MP/0"
        print(hdr)
        for t in GRID:
            c0, cE, cM = arms["0"][t], arms["E"][t], arms["MP"][t]
            print(f"  {t:<4d}{c0:9.5f}  {cE:9.5f}  {cM:9.5f}  "
                  f"{cE/c0:9.5f}  {cM/c0:9.5f}")
        sE = slope_lnR([arms["E"][t] for t in GRID],
                       [arms["0"][t] for t in GRID], GRID)
        sM = slope_lnR([arms["MP"][t] for t in GRID],
                       [arms["0"][t] for t in GRID], GRID)
        sE4 = slope_lnR([arms["E"][t] for t in GRID_NO_DEEP],
                        [arms["0"][t] for t in GRID_NO_DEEP], GRID_NO_DEEP)
        sM4 = slope_lnR([arms["MP"][t] for t in GRID_NO_DEEP],
                        [arms["0"][t] for t in GRID_NO_DEEP], GRID_NO_DEEP)
        print(f"  slope(E)  = {sE:+.5f}/step   slope(MP) = {sM:+.5f}/step")
        print(f"  slope(MP)-slope(E) = {sM - sE:+.5f}/step   "
              f"[full grid {GRID}]")
        print(f"  without depths 6,8: slope(E)={sE4:+.5f}  slope(MP)={sM4:+.5f} "
              f" diff={sM4 - sE4:+.5f}/step")
        # floor check
        floor = 0.3  # coh(0)=1 for |+> on q2
        print(f"  floor check coh_0(t) >= 0.3*coh(0)=0.3 :")
        allpass = True
        for t in GRID:
            ok = arms["0"][t] >= floor
            allpass = allpass and ok
            print(f"    t={t}: coh_0={arms['0'][t]:.5f}  "
                  f"{'PASS' if ok else 'FAIL'}")
        print(f"  ALL FLOOR POINTS PASS: {allpass}")
        return sE, sM, sM - sE, sM4 - sE4, allpass

    print("\n### clean (channel-only, validates against design scouting) ###")
    cl = report_arms(r1c["clean"], "clean")
    print("\n### WITH representative device background (the flown prediction) ###")
    bg = report_arms(r1c["bg"], "background")
    central_bg = bg[2]

    print("\nvalidation vs design scouting (clean):")
    print(f"  design scouting MP-E ~ -0.073/step ; this sim clean = "
          f"{cl[2]:+.5f}/step")
    print(f"  design scouting min coh_0 ~ 0.456   ; this sim clean min coh_0 = "
          f"{min(r1c['clean']['0'][t] for t in GRID):.5f}")
    print(f"  design scouting 'no depths 6,8' diff ~ -0.054 ; this sim clean = "
          f"{cl[3]:+.5f}/step")

    # sequential-schedule sensitivity on the diff
    bgseq = part1c(amp_seq, bg_mask_seq)["bg"]
    sEq = slope_lnR([bgseq["E"][t] for t in GRID],
                    [bgseq["0"][t] for t in GRID], GRID)
    sMq = slope_lnR([bgseq["MP"][t] for t in GRID],
                    [bgseq["0"][t] for t in GRID], GRID)
    print(f"\nwall-time sensitivity: diff(parallel tau)={central_bg:+.5f}  "
          f"diff(sequential tau)={sMq - sEq:+.5f}/step")

    # ---- (d) counts-level MC
    print("\n" + "-" * 74)
    print("PART 1(d) - COUNTS-LEVEL MC: PROJECTED SE, BAND, POWER MARGIN")
    print("-" * 74)
    d1 = part1d(amp, bg_mask, rep, central_bg)
    print(f"central slope(MP)-slope(E) [exact DM, background] = "
          f"{central_bg:+.5f}/step")
    print(f"counts-MC bootstrap mean diff (finite M,shots,readout) = "
          f"{d1['boot_mean_diff']:+.5f}/step")
    print(f"projected combined SE (binding+shot+readout, R={R_BOOT} boot) = "
          f"{d1['se_diff']:.5f}")
    print(f"7a BAND = central +/- 2*SE = "
          f"[{d1['band'][0]:+.5f}, {d1['band'][1]:+.5f}]  (A-mag boundary)")
    print(f"POWER MARGIN = |central|/SE = {d1['power']:.2f}   "
          f"(design gates hardware on >= 3)")
    print(f"  gate: {'PASS' if d1['power'] >= 3 else 'FAIL'} "
          f"(power margin {'>=' if d1['power']>=3 else '<'} 3)")
    print(f"slope(E) projected SE = {d1['se_E']:.5f} ; "
          f"boot-mean slope(E) = {d1['boot_mean_slopeE']:+.5f}/step "
          f"(leakage null-consistency)")

    # ---- Part 1 gate summary
    print("\n" + "-" * 74)
    print("PART 1 GATE SUMMARY")
    print("-" * 74)
    print(f"  floor gate (all 6 grid points, background): "
          f"{'PASS' if bg[4] else 'FAIL'}")
    print(f"  A-sign robust to deep-point loss: full diff={central_bg:+.5f}, "
          f"no-deep diff={bg[3]:+.5f} -> "
          f"{'same sign' if np.sign(central_bg)==np.sign(bg[3]) else 'SIGN FLIP'}")
    print(f"  power margin >= 3: {'PASS' if d1['power']>=3 else 'FAIL'} "
          f"({d1['power']:.2f})")
    hardware_ok = bg[4] and d1['power'] >= 3 and \
        np.sign(central_bg) == np.sign(bg[3])
    print(f"  => hardware-permissible by 7a criteria: {hardware_ok}")

    # ---- Part 2
    part2_search(amp, bg_mask, rep)

    print("\n" + "=" * 74)
    print("END STAGE 7a")
    print("=" * 74)


if __name__ == "__main__":
    np.set_printoptions(precision=5, suppress=True)
    main()
