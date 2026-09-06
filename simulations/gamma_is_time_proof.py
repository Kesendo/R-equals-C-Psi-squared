#!/usr/bin/env python3
"""
What gamma does to a two-qubit trajectory, and what tau = gamma*t does not do
============================================================================
Three finite measurements on the N=2 Heisenberg chain under Z-dephasing.

Part 1: which trajectory diagnostics change when gamma is switched on
        (recurrence of the trace distance, monotonicity of S(rho_A),
        CPsi crossings of 1/4, the late-time speed ||drho/dt||).
Part 2: the same diagnostics across (J, gamma) configurations, so the
        J-only and gamma-only corners can be read side by side.
Part 3: a falsification, run in both directions. If tau = gamma*t were the
        trajectory's own time variable, curves plotted against tau would
        collapse across gamma. At fixed J they do not, because holding J while
        sweeping gamma also moves Q = J/gamma; at fixed Q they collapse
        exactly. The generator needs two numbers and gamma supplies one.

These are finite runs at N=2, on one Hamiltonian and one channel. They
show what gamma does to these trajectories; they do not establish that
gamma IS time, and Part 3 is the measurement that blocks that reading.

Script: simulations/gamma_is_time_proof.py
Output: simulations/results/gamma_is_time_proof.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "gamma_is_time_proof.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 2; d = 4; d2 = 16


def site_op(op, k):
    return np.kron(op, I2) if k == 0 else np.kron(I2, op)


def build_H(J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * site_op(P, 0) @ site_op(P, 1)
    return H


def build_L(H, gammas):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def ptrace_A(rho):
    return np.array([[rho[0,0]+rho[1,1], rho[0,2]+rho[1,3]],
                     [rho[2,0]+rho[3,1], rho[2,2]+rho[3,3]]])


def entropy(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-15]
    return -np.sum(ev * np.log2(ev))


def trace_dist(r1, r2):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh(r1 - r2)))


def cpsi(rho):
    C = np.real(np.trace(rho @ rho))
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    return C * l1 / 3.0


def concurrence(rho):
    sy2 = np.kron(sy, sy)
    R = rho @ sy2 @ rho.conj() @ sy2
    ev = np.sort(np.sqrt(np.maximum(np.real(np.linalg.eigvals(R)), 0)))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])


def make_state(name):
    if name == 'Bell+':
        psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    elif name == '|01>':
        psi = np.kron(up, dn)
    elif name == '|+0>':
        psi = np.kron(plus, up)
    return np.outer(psi, psi.conj())



# ============================================================
log("What gamma does to a two-qubit trajectory, and what tau = gamma*t does not do")
log()

# ============================================================
# PART 1: WHAT SWITCHING GAMMA ON CHANGES
# ============================================================
log("=" * 70)
log("PART 1: WHAT SWITCHING GAMMA ON CHANGES")
log("Four trajectory diagnostics, gamma = 0 against gamma = 0.05")
log("=" * 70)
log()

H = build_H()
tlist = np.linspace(0, 50, 501)

for sname in ['|01>', 'Bell+']:
    rho0 = make_state(sname)

    for gamma in [0, 0.05]:
        L = build_L(H, [gamma, gamma])

        S_series = []
        D_series = []
        cross_down = 0
        cross_up = 0
        prev_above = cpsi(rho0) > 0.25

        for t in tlist:
            rho = evolve(L, rho0, t)
            rA = (ptrace_A(rho) + ptrace_A(rho).conj().T) / 2
            S_series.append(entropy(rA))
            D_series.append(trace_dist(rho, rho0))
            cp = cpsi(rho)
            above = cp > 0.25
            if prev_above and not above:
                cross_down += 1
            if not prev_above and above:
                cross_up += 1
            prev_above = above

        S = np.array(S_series)
        D = np.array(D_series)

        idx_1 = np.argmin(np.abs(tlist - 1))
        S_monotonic = all(S[i + 1] >= S[i] - 1e-8 for i in range(idx_1, len(S) - 1))
        D_returns = any(D[i] < 0.01 for i in range(50, len(D)))

        rho_49 = evolve(L, rho0, 49.9)
        rho_50 = evolve(L, rho0, 50.0)
        drho_norm = np.linalg.norm(rho_50 - rho_49) / 0.1

        log(f"  {sname:>6}  gamma={gamma:.2f}:")
        log(f"    S(rho_A) non-decreasing after t=1:  {'yes' if S_monotonic else 'no'}")
        log(f"    trace distance returns (< 0.01):    {'yes' if D_returns else 'no'}")
        log(f"    CPsi crossings of 1/4, down/up:     {cross_down}/{cross_up}")
        log(f"    ||drho/dt|| at t=50:                {drho_norm:.6f}")
        log()

log("  At gamma = 0 the |01> trajectory is exactly recurrent. H has spectrum {-3J, J},")
log("  a single gap 4J, so the period is pi/(2J) and the return count over the window")
log("  is set by that period, not by this grid: a coarser grid simply misses narrower")
log("  dips. At gamma = 0.05 the returns are gone and the state has all but stopped.")
log("  Recurrence is what the dephasing removes, and that is a statement about these")
log("  trajectories.")
log()
log()

# ============================================================
# PART 2: THE SAME DIAGNOSTICS ACROSS (J, gamma)
# ============================================================
log("=" * 70)
log("PART 2: THE SAME DIAGNOSTICS ACROSS (J, gamma)")
log("The J-only and the gamma-only corner side by side")
log("=" * 70)
log()

rho0 = make_state('|+0>')

log(f"  {'Config':>20}  {'S mono':>8}  {'D no return':>12}  {'CPsi 1-way':>10}  {'both':>6}")
log(f"  {'-' * 62}")

configs = [
    ('J=0.1, g=0.05', 0.1, 0.05),
    ('J=1.0, g=0.05', 1.0, 0.05),
    ('J=10,  g=0.05', 10.0, 0.05),
    ('J=0,   g=0.05', 0.0, 0.05),
    ('J=1.0, g=0', 1.0, 0.0),
]

grid_30 = np.linspace(0, 30, 301)

for name, J, gamma in configs:
    H_c = build_H(J)
    L_c = build_L(H_c, [gamma, gamma])

    S_list = []
    D_list = []
    cross_d = 0
    cross_u = 0
    prev = cpsi(rho0) > 0.25

    for t in grid_30:
        rho = evolve(L_c, rho0, t)
        rA = (ptrace_A(rho) + ptrace_A(rho).conj().T) / 2
        S_list.append(entropy(rA))
        D_list.append(trace_dist(rho, rho0))
        cp = cpsi(rho)
        above = cp > 0.25
        if prev and not above:
            cross_d += 1
        if not prev and above:
            cross_u += 1
        prev = above

    S = np.array(S_list)
    idx1 = np.argmin(np.abs(grid_30 - 1))
    S_mono = all(S[i + 1] >= S[i] - 1e-8 for i in range(idx1, len(S) - 1))
    D_noreturn = not any(D_list[i] < 0.01 for i in range(30, len(D_list)))
    one_way = cross_d > 0 and cross_u == 0
    both = S_mono and D_noreturn

    log(f"  {name:>20}  {'yes' if S_mono else 'no':>8}  "
        f"{'yes' if D_noreturn else 'no':>12}  "
        f"{'yes' if one_way else 'no':>10}  "
        f"{'yes' if both else 'no':>6}")

log()
log("  Only the J=0 row carries both diagnostics at once, and the gamma=0 row carries")
log("  neither. On this state and this grid the pair separates the pure-decay corner")
log("  from the unitary one. The last column is the conjunction of the two measured")
log("  columns and nothing beyond it; it is not a test for the presence of time.")
log()

# ============================================================
# PART 3: THE tau = gamma*t COLLAPSE, AND ITS FAILURE
# ============================================================
log("=" * 70)
log("PART 3: THE tau = gamma*t COLLAPSE, AND ITS FAILURE")
log("If tau were the trajectory's own time variable, these curves would lie on each other")
log("=" * 70)
log()

rho0 = make_state('|01>')
H = build_H()

tau_points = np.linspace(0, 1.0, 100)
OBSERVABLES = ('S(rho_A)', 'Tr(rho^2)', 'CPsi', 'Concurrence')


def collapse_table(pairs):
    """Spread of each observable across a set of (J, gamma) runs, at equal tau.

    `pairs` is a list of (J, gamma). Every run is sampled on the same tau grid,
    so t = tau/gamma differs between runs by construction; the question is
    whether the curves land on each other anyway.
    """
    curves = {name: {} for name in OBSERVABLES}
    for coupling, g in pairs:
        generator = build_L(build_H(coupling), [g, g])
        series = {name: [] for name in OBSERVABLES}
        for tau in tau_points:
            rho = evolve(generator, rho0, tau / g)
            rA = (ptrace_A(rho) + ptrace_A(rho).conj().T) / 2
            series['S(rho_A)'].append(entropy(rA))
            series['Tr(rho^2)'].append(np.real(np.trace(rho @ rho)))
            series['CPsi'].append(cpsi(rho))
            series['Concurrence'].append(concurrence(rho))
        for name in OBSERVABLES:
            curves[name][g] = np.array(series[name])

    rows = []
    for name in OBSERVABLES:
        gs = list(curves[name].keys())
        spread = max(np.max(np.abs(curves[name][gs[i]] - curves[name][gs[j]]))
                     for i in range(len(gs)) for j in range(i + 1, len(gs)))
        stacked = np.concatenate([curves[name][g] for g in gs])
        rng = float(np.max(stacked) - np.min(stacked))
        flat = rng <= 1e-12
        ratio = float('inf') if flat else spread / rng
        rows.append((name, spread, rng, ratio, (not flat) and ratio < 0.05))
    return rows


def log_table(title, rows):
    log(f"  {title}")
    log(f"  {'Observable':>20}  {'spread':>10}  {'range':>10}  {'spread/range':>13}  {'collapses?':>11}")
    log(f"  {'-' * 72}")
    for name, spread, rng, ratio, collapses in rows:
        ratio_txt = f"{'n/a (flat)':>13}" if not np.isfinite(ratio) else f"{ratio:13.4f}"
        log(f"  {name:>20}  {spread:10.6f}  {rng:10.6f}  {ratio_txt}  "
            f"{'yes' if collapses else 'no':>11}")
    log()


log("  The generator is homogeneous: L(J, gamma)*t = tau * L(J/gamma, 1) with tau = gamma*t.")
log("  So every observable is a function of tau AND of Q = J/gamma, and the test below has")
log("  to be run twice or it says nothing. Spread is the largest gap between two curves at")
log("  equal tau; range is how far the observable itself travels; collapse is called only")
log("  under 5% of the range, so an observable cannot pass by standing still.")
log()

# The negative arm: hold J and sweep gamma, which necessarily sweeps Q as well.
sweep_pairs = [(1.0, g) for g in (0.01, 0.02, 0.05, 0.10, 0.20)]
sweep_rows = collapse_table(sweep_pairs)
log_table("Fixed J = 1, gamma from 0.01 to 0.20 (so Q sweeps 100 to 5):", sweep_rows)

# The positive arm: hold Q and sweep gamma, where the homogeneity says it MUST collapse.
Q_FIXED = 20.0
control_pairs = [(Q_FIXED * g, g) for g in (0.01, 0.05, 0.20)]
control_rows = collapse_table(control_pairs)
log_table(f"Fixed Q = {Q_FIXED:g}, gamma from 0.01 to 0.20 (J moves with it):", control_rows)

collapsed_in_sweep = [name for name, _, _, _, ok in sweep_rows if ok]
missed_in_control = [name for name, _, _, _, ok in control_rows if not ok]
if collapsed_in_sweep:
    raise RuntimeError(
        "an observable collapsed while Q was moving, which the homogeneity forbids: "
        + ", ".join(collapsed_in_sweep))
if missed_in_control:
    raise RuntimeError(
        "an observable failed to collapse at fixed Q, where the homogeneity forces it: "
        + ", ".join(missed_in_control))

log("  Both arms come out as the identity says. At fixed Q the four curves land on each")
log("  other to machine precision, so tau = gamma*t IS the natural time of this generator.")
log("  At fixed J it fails, and not because tau is the wrong time variable: holding J")
log("  sweeping gamma moves Q = J/gamma, so those five runs are five different systems")
log("  compared at matched tau. What the failure measures is the second knob, not the")
log("  first. Reading the left table alone, as 'irreversible observables do not scale")
log("  with tau', mistakes a change of system for a failure of the time variable.")
log()

# ============================================================
log("=" * 70)
log("WHAT THESE THREE RUNS SHOW")
log("=" * 70)
log()
log("Part 1: switching gamma on removes recurrence and brings the trajectory to")
log("        rest. Measured on two states at N=2.")
log()
log("Part 2: on this state and this grid the two diagnostics together pick out the")
log("        pure-decay corner J=0 and reject the unitary corner gamma=0.")
log()
log("Part 3: tau = gamma*t is this generator's own time, and the two-armed test says")
log("        so: at fixed Q the curves collapse exactly, at fixed J they do not,")
log("        because holding J while sweeping gamma changes Q and therefore the")
log("        system. tau rescales the trajectory only together with Q.")
log()
log("Together: gamma sets the scale on which these trajectories stop returning, and")
log("it supplies one of the two numbers a trajectory here needs. A quantity that")
log("only sets the time variable alongside a second knob is not what 'gamma is time' claims.")
log()

log("Account: docs/GAMMA_TIME_DISTINCTION.md")
log("Results: simulations/results/gamma_is_time_proof.txt")
_outf.close()
