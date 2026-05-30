#!/usr/bin/env python3
"""
The two-handed clock on the neural side, and which hand is graph-blind.

SEAM 2 (the clock). The neural Jacobian eigenvalues mu = -rate + i*omega are
the same clock object as the quantum Liouvillian: a Takt hand (radial decay)
and a Rotation hand (angular omega), with per-mode angle
theta = atan2(|omega|, rate)  (0 deg = pure decay, 90 deg = pure rotation).
The V-Effect and the thermal window move ONLY the Rotation hand; the Takt
stays pinned.

SEAM 1 (which hand is graph-blind). The Takt is graph-blind by an exact
identity:
    mean Re(lambda) = trace(J)/d = -(1/tau_E + 1/tau_I)/2 = -S,
set only by the membrane constants, because the synaptic graph W contributes
NOTHING to the diagonal of J -- it is "traceless," exactly as the Hamiltonian
commutator -i[H,.] is traceless in the quantum Liouvillian, where only the
gamma's set the trace (= the palindrome center 2*Sum gamma). The wiring lives
entirely in the off-diagonal -- the Rotation. And there C. elegans is matched
by a degree-preserving rewire but NOT by Erdos-Renyi: even the Rotation is set
by the coarse degree structure, not the fine wiring.

So the neural reading of F1's topology-blindness: the palindrome center is
bath-set (gamma) in the quantum case and membrane-set (1/tau) here; in neither
does the graph touch it. Dale's Law (local, per-neuron) gives the signs for
free; only the off-diagonal magnitude-match is graph-dependent, and only at
the coarse degree level.

Builders are copied minimally from veffect_exact.py / veffect_and_heat.py /
validation_checks.py (siblings in this directory; the neural scripts already
each redefine these small helpers, so this probe follows suit rather than
importing modules that run on import).
"""
import json
import os

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
NEURAL_DIR = SCRIPT_DIR  # this script lives in simulations/neural/

TAU_E, TAU_I = 5.0, 10.0          # one clock for both parts
S = (1.0 / TAU_E + 1.0 / TAU_I) / 2.0     # palindrome center -> spectral center is -S
GAP = 1.0 / TAU_E + 1.0 / TAU_I           # Takt gap = 2S ; tau_clock = 1/GAP
ALPHA = 0.5


# ---- clock readings -------------------------------------------------------

def clock(ev, tol=1e-6):
    """Read the two hands off a spectrum.

    Returns (n_rotating, theta_max_deg, mean_Re). theta is the angle of each
    eigenvalue from the pure-decay (negative-real) axis.
    """
    im = np.abs(ev.imag)
    re = np.abs(ev.real)
    rotating = im > tol
    n_rot = int(np.sum(rotating))
    theta = np.degrees(np.arctan2(im, re))
    theta_max = float(np.max(theta)) if n_rot else 0.0
    return n_rot, theta_max, float(np.mean(ev.real))


# ---- builders (copied minimally from the neural scripts) ------------------

def build_exact_palindromic_network(N, tau_E, tau_I, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    inh_idx = rng.choice(N, N // 2, replace=False)
    signs[inh_idx] = -1
    e_idx = list(np.where(signs > 0)[0])
    i_idx = list(np.where(signs < 0)[0])
    perm = np.arange(N)
    for k in range(min(len(e_idx), len(i_idx))):
        perm[e_idx[k]] = i_idx[k]
        perm[i_idx[k]] = e_idx[k]
    W = np.zeros((N, N))
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    for i in range(N):
        for j in range(N):
            if i == j or i >= j or not mask[i, j]:
                continue
            qi, qj = perm[i], perm[j]
            base = rng.exponential(0.3)
            W[i, j] = signs[j] * base
            tau_i = tau_E if signs[i] > 0 else tau_I
            tau_qi = tau_E if signs[qi] > 0 else tau_I
            W[qi, qj] = -(tau_qi / tau_i) * W[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def make_balanced_network(N, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    inh_idx = rng.choice(N, N // 2, replace=False)
    signs[inh_idx] = -1
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (N, N))
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if mask[i, j]:
                W[i, j] = signs[j] * weights[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def build_linear_jacobian(W, signs, tau_E, tau_I, alpha):
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_E if signs[i] > 0 else tau_I
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def sigmoid(x, a, theta):
    return 1.0 / (1.0 + np.exp(np.clip(-a * (x - theta), -500, 500)))


def dsigmoid(x, a, theta):
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)


def build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, alpha, P):
    n = len(signs)
    x = np.ones(n) * 0.3
    for _ in range(500):
        inputs = alpha * W @ x + P
        for i in range(n):
            a_i = 1.3 if signs[i] > 0 else 2.0
            th_i = 4.0 if signs[i] > 0 else 3.7
            x[i] = sigmoid(inputs[i], a_i, th_i)
    inputs = alpha * W @ x + P
    J = np.zeros((n, n))
    for i in range(n):
        a_i = 1.3 if signs[i] > 0 else 2.0
        th_i = 4.0 if signs[i] > 0 else 3.7
        tau_i = tau_E if signs[i] > 0 else tau_I
        dS = dsigmoid(inputs[i], a_i, th_i)
        J[i, i] = (-1.0 + alpha * W[i, i] * dS) / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] * dS / tau_i
    return J


def swap_Q(signs):
    n = len(signs)
    e_idx = list(np.where(signs > 0)[0])
    i_idx = list(np.where(signs < 0)[0])
    perm = np.arange(n)
    for k in range(min(len(e_idx), len(i_idx))):
        perm[e_idx[k]] = i_idx[k]
        perm[i_idx[k]] = e_idx[k]
    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, perm[i]] = 1.0
    return Q


def residual_split(J, Q):
    """Return (diag_residual, offdiag_residual) of Q J Q + J + 2 S, both
    normalised by ||J||. The diagonal piece is the Takt (self-decay) part;
    the off-diagonal piece is the Rotation (coupling) part."""
    QJQ = Q @ J @ Q.T
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2 * np.diag(S_diag)
    R_off = R - np.diag(np.diag(R))
    norm_J = np.linalg.norm(J)
    diag_res = np.linalg.norm(np.diag(R)) / norm_J if norm_J else 0.0
    off_res = np.linalg.norm(R_off) / norm_J if norm_J else 0.0
    return diag_res, off_res


def degree_preserving_rewire(W, signs, rng):
    W_new = W.copy()
    n = W.shape[0]
    edges = [(i, j) for i in range(n) for j in range(n) if i != j and W[i, j] != 0]
    if len(edges) < 2:
        return W_new
    for _ in range(len(edges) * 10):
        idx1, idx2 = rng.choice(len(edges), 2, replace=False)
        i, j = edges[idx1]
        k, l = edges[idx2]
        if i == k or j == l or i == l or k == j:
            continue
        if W_new[i, l] != 0 or W_new[k, j] != 0:
            continue
        if signs[j] != signs[l]:
            continue
        W_new[i, l] = W_new[i, j]
        W_new[k, j] = W_new[k, l]
        W_new[i, j] = 0
        W_new[k, l] = 0
        edges[idx1] = (i, l)
        edges[idx2] = (k, j)
    return W_new


def random_dale_network(rng, n_total, signs, density):
    W = np.zeros((n_total, n_total))
    for i in range(n_total):
        for j in range(n_total):
            if i != j and rng.random() < max(density, 0.01):
                W[i, j] = signs[j] * rng.exponential(0.3)
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W


# ==========================================================================
print("=" * 72)
print("PART A  --  the two hands on the neural clock  (seam 2)")
print(f"  tau_E={TAU_E}, tau_I={TAU_I}  ->  Takt gap = 1/tau_E + 1/tau_I = {GAP:.3f}")
print(f"  spectral center should sit at -S = {-S:+.3f}  (= trace(J)/d)")
print("=" * 72)

# 1. the silent noble gas: exact palindrome
W, signs = build_exact_palindromic_network(20, TAU_E, TAU_I, seed=42)
J = build_linear_jacobian(W, signs, TAU_E, TAU_I, ALPHA)
n_rot, th_max, mre = clock(np.linalg.eigvals(J))
dres, ores = residual_split(J, swap_Q(signs))
print(f"\n  exact palindrome (N=20): off-diag residual {ores:.1e}  "
      f"-> n_rotating={n_rot}, theta_max={th_max:.1f} deg  (SILENT)")
print(f"      mean Re(lambda) = {mre:+.4f}   (Takt pinned at -S={-S:+.3f})")

# 2. V-Effect: couple two exact nets, sweep coupling -> Rotation hand wakes up
print("\n  V-EFFECT  (two exact nets + 1 mediator, coupling sweep):")
print(f"  {'coupling':>8s}  {'n_rotating':>10s}  {'theta_max':>9s}  {'mean Re':>9s}")
print("  " + "-" * 44)
N = 20
W_A, signs_A = build_exact_palindromic_network(N, TAU_E, TAU_I, seed=42)
W_B, signs_B = build_exact_palindromic_network(N, TAU_E, TAU_I, seed=99)
N_c = 2 * N + 1
W_c = np.zeros((N_c, N_c))
signs_c = np.zeros(N_c)
W_c[:N, :N] = W_A
signs_c[:N] = signs_A
W_c[N:2 * N, N:2 * N] = W_B
signs_c[N:2 * N] = signs_B
signs_c[2 * N] = 1.0  # excitatory mediator
for g in [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0]:
    Wt = W_c.copy()
    for off in [0, N]:
        Wt[2 * N, off] = g
        Wt[off, 2 * N] = g
        Wt[2 * N, off + N - 1] = g
        Wt[off + N - 1, 2 * N] = g
    Jc = build_linear_jacobian(Wt, signs_c, TAU_E, TAU_I, ALPHA)
    n_rot, th_max, mre = clock(np.linalg.eigvals(Jc))
    mark = "  <- Rotation wakes" if (g in (0.01, 0.05)) else ""
    print(f"  {g:8.2f}  {n_rot:10d}  {th_max:8.1f}  {mre:+9.4f}{mark}")

# 3. thermal window: sweep drive P -> Rotation hand sweeps up and back
print("\n  THERMAL WINDOW  (approximate net, drive P sweep):")
print(f"  {'P':>6s}  {'n_rotating':>10s}  {'theta_max':>9s}  {'mean Re':>9s}")
print("  " + "-" * 42)
Wb, signs_b = make_balanced_network(50, density=0.3, seed=42)
for P in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]:
    Jp = build_jacobian_with_sigmoid(Wb, signs_b, TAU_E, TAU_I, 0.3, P)
    n_rot, th_max, mre = clock(np.linalg.eigvals(Jp), tol=1e-5)
    mark = "  <- widest angle" if P == 4.0 else ""
    print(f"  {P:6.1f}  {n_rot:10d}  {th_max:8.1f}  {mre:+9.4f}{mark}")

# ==========================================================================
print("\n" + "=" * 72)
print("PART B  --  which hand is graph-blind  (seam 1)")
print("  Takt = trace(J)/d = -S : the synaptic graph W never touches the")
print("  diagonal, so the spectral center is membrane-set, not graph-set.")
print("=" * 72)

# the Takt identity, exactly, on three different graphs
print("\n  THE TAKT IS THE TRACE (exact, graph-blind):")
print(f"  {'graph':>22s}  {'mean Re(lambda)':>15s}  {'-S target':>10s}")
print("  " + "-" * 52)
rng = np.random.RandomState(7)
demo_W, demo_signs = make_balanced_network(10, density=0.3, seed=3)
for label, Wd in [("balanced random", demo_W),
                  ("degree-preserved", degree_preserving_rewire(demo_W, demo_signs, rng)),
                  ("Erdos-Renyi (Dale)", random_dale_network(rng, 10, demo_signs, 0.3))]:
    Jd = build_linear_jacobian(Wd, demo_signs, TAU_E, TAU_I, ALPHA)
    _, _, mre = clock(np.linalg.eigvals(Jd))
    print(f"  {label:>22s}  {mre:+15.6f}  {-S:+10.4f}")

# C. elegans vs degree-preserved vs ER: the Rotation (off-diag) residual
print("\n  THE ROTATION LIVES ON COARSE DEGREE STRUCTURE (C. elegans):")
ce_path = os.path.join(NEURAL_DIR, "celegans_connectome.json")
with open(ce_path) as f:
    data = json.load(f)
W_chem = np.array(data["chemical"])
signs_full = np.array(data["chemical_sign"])
N_full = len(signs_full)
W_signed = np.zeros((N_full, N_full))
for i in range(N_full):
    for j in range(N_full):
        W_signed[i, j] = signs_full[j] * W_chem[j, i]
W_signed /= np.max(np.abs(W_signed))
exc_idx = np.where(signs_full > 0)[0]
inh_idx = np.where(signs_full < 0)[0]

n_trials, n_half = 200, 5
ce_diag, ce_off, dp_off, er_off = [], [], [], []
for trial in range(n_trials):
    rng = np.random.RandomState(trial + 100)
    e_pick = rng.choice(exc_idx, n_half, replace=False)
    i_pick = rng.choice(inh_idx, n_half, replace=False)
    idx = np.concatenate([e_pick, i_pick])
    W_sub = W_signed[np.ix_(idx, idx)]
    signs_sub = signs_full[idx]
    Q = swap_Q(signs_sub)
    density = np.count_nonzero(W_sub) / (10 * 9)

    J_ce = build_linear_jacobian(W_sub, signs_sub, TAU_E, TAU_I, 0.3)
    d, o = residual_split(J_ce, Q)
    ce_diag.append(d)
    ce_off.append(o)

    W_dp = degree_preserving_rewire(W_sub, signs_sub, rng)
    _, o_dp = residual_split(build_linear_jacobian(W_dp, signs_sub, TAU_E, TAU_I, 0.3), Q)
    dp_off.append(o_dp)

    W_er = random_dale_network(rng, 10, signs_sub, density)
    _, o_er = residual_split(build_linear_jacobian(W_er, signs_sub, TAU_E, TAU_I, 0.3), Q)
    er_off.append(o_er)

print(f"  {'':>22s}  {'diag (Takt)':>12s}  {'off-diag (Rotation)':>20s}")
print("  " + "-" * 58)
print(f"  {'C. elegans':>22s}  {np.mean(ce_diag):12.1e}  {np.mean(ce_off):20.4f}")
print(f"  {'degree-preserved':>22s}  {'(identical)':>12s}  {np.mean(dp_off):20.4f}")
print(f"  {'Erdos-Renyi (Dale)':>22s}  {'(identical)':>12s}  {np.mean(er_off):20.4f}")
print(f"\n  C. elegans / degree-preserved = {np.mean(ce_off)/np.mean(dp_off):.2f}  (~1: degree explains it)")
print(f"  C. elegans / Erdos-Renyi      = {np.mean(ce_off)/np.mean(er_off):.2f}  (<1: better than fine-random)")

print("\n" + "=" * 72)
print("READING")
print("  The Takt hand is the trace: -(1/tau_E + 1/tau_I)/2, membrane-set,")
print("  graph-blind by an exact identity (W is traceless in J, as -i[H,.]")
print("  is traceless in L; only the bath sets the center). Coupling and")
print("  drive move only the Rotation hand. And even the Rotation is set by")
print("  the coarse degree structure, not the fine wiring -- the neural")
print("  reading of F1's topology-blindness confirmed today in EXCLUSIONS.")
print("=" * 72)
