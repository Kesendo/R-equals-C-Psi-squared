#!/usr/bin/env python3
"""
Bridge Fingerprints: State Identification Through CΨ Crossing Signatures
=========================================================================
4-qubit system: A (qubits 0,1) = receiver in |00>, B (qubits 2,3) = sender.
Heisenberg coupling within A and B (J_internal=1.0), bridge coupling 1<->2
(J_bridge variable). Local Z-dephasing gamma=0.1 on all qubits.

Eight sender states produce distinct CΨ_A trajectories. Product states with
local coherence in the bridge qubit cross 1/4; Bell states never do.

Script:  simulations/bridge_fingerprints.py
Output:  simulations/results/bridge_fingerprints.txt
Docs:    experiments/BRIDGE_FINGERPRINTS.md
"""

import numpy as np
import os, sys, time as _time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "bridge_fingerprints.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# INFRASTRUCTURE (N=4)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)
minus = (up - dn) / np.sqrt(2)

N = 4
d = 2 ** N    # 16
d2 = d * d    # 256


def site_op(op, k, nq=N):
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H(J_internal=1.0, J_bridge=0.5):
    """Heisenberg: H_A (0-1) + H_B (2-3) + H_bridge (1-2)."""
    H = np.zeros((d, d), dtype=complex)
    bonds = [(0, 1, J_internal), (2, 3, J_internal), (1, 2, J_bridge)]
    for i, j, J in bonds:
        for P in [sx, sy, sz]:
            H += J * site_op(P, i) @ site_op(P, j)
    return H


def build_L(H, gamma):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def ptrace_keep(rho, keep, nq=N):
    """Partial trace keeping specified qubits."""
    keep = list(keep)
    trace_out = [q for q in range(nq) if q not in keep]
    dims = [2] * nq
    reshaped = rho.reshape(dims + dims)
    current_n = nq
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape(d_k, d_k)


def cpsi(rho_sub):
    """CΨ = Tr(rho^2) * L1/(d-1) where d = dim of subsystem."""
    d_sub = rho_sub.shape[0]
    purity = np.real(np.trace(rho_sub @ rho_sub))
    l1 = np.sum(np.abs(rho_sub)) - np.real(np.trace(rho_sub))
    psi = l1 / (d_sub - 1) if d_sub > 1 else 0.0
    return purity * psi


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def euler_evolve(L, rho0_vec, dt, n_steps):
    """First-order Euler integration of drho/dt = L*rho."""
    v = rho0_vec.copy()
    trajectory = np.empty((n_steps + 1, len(v)), dtype=complex)
    trajectory[0] = v
    for i in range(n_steps):
        v = v + dt * (L @ v)
        trajectory[i + 1] = v
    return trajectory


# ============================================================
# SENDER STATES
# ============================================================
bell_plus = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
bell_minus = (np.kron(up, up) - np.kron(dn, dn)) / np.sqrt(2)
psi_plus = (np.kron(up, dn) + np.kron(dn, up)) / np.sqrt(2)

sender_states = {
    "|++>":   np.kron(plus, plus),
    "|+0>":   np.kron(plus, up),
    "Bell+":  bell_plus,
    "Bell-":  bell_minus,
    "|Psi+>": psi_plus,
    "|+->":   np.kron(plus, minus),
    "|01>":   np.kron(up, dn),
    "|11>":   np.kron(dn, dn),
}

A_state = np.kron(up, up)  # |00>_A


def extract_fingerprint(cpsi_traj, times, dt):
    """Extract 8 fingerprint metrics from CΨ_A trajectory."""
    peak = np.max(cpsi_traj)
    i_peak = np.argmax(cpsi_traj)
    t_peak = times[i_peak]

    t_up = None
    t_down = None
    above_time = 0.0
    above = False

    for i in range(1, len(cpsi_traj)):
        if cpsi_traj[i] >= 0.25 and cpsi_traj[i - 1] < 0.25:
            if t_up is None:
                t_up = times[i]
            above = True
        if cpsi_traj[i] < 0.25 and cpsi_traj[i - 1] >= 0.25:
            if t_down is None:
                t_down = times[i]
            above = False
        if cpsi_traj[i] >= 0.25:
            above_time += dt

    integral = np.trapezoid(cpsi_traj, times)

    # Rise rate: slope approaching peak
    i_start = max(0, i_peak - 50)
    if i_peak > i_start:
        rise_rate = (cpsi_traj[i_peak] - cpsi_traj[i_start]) / (times[i_peak] - times[i_start])
    else:
        rise_rate = 0.0

    # Fall rate: slope after peak
    i_end = min(len(cpsi_traj) - 1, i_peak + 50)
    if i_end > i_peak:
        fall_rate = (cpsi_traj[i_end] - cpsi_traj[i_peak]) / (times[i_end] - times[i_peak])
    else:
        fall_rate = 0.0

    return {
        "peak": peak,
        "t_peak": t_peak,
        "t_up": t_up,
        "t_down": t_down,
        "above_time": above_time,
        "integral": integral,
        "rise_rate": rise_rate,
        "fall_rate": fall_rate,
        "crosses": t_up is not None,
    }


# ============================================================
# MAIN: SWEEP J_bridge VALUES
# ============================================================
log("=" * 78)
log("Bridge Fingerprints: State Identification Through CΨ Crossing Signatures")
log("=" * 78)
log()

gamma = 0.1
dt = 0.001
t_max = 5.0
n_steps = int(t_max / dt)
times = np.linspace(0, t_max, n_steps + 1)

J_bridge_values = [0.3, 0.5, 0.7, 1.0, 1.5]
J_internal = 1.0

for J_bridge in J_bridge_values:
    J_over_gamma = J_bridge / gamma
    log(f"J_bridge = {J_bridge}, J/γ = {J_over_gamma:.0f}")
    log("-" * 78)

    H = build_H(J_internal=J_internal, J_bridge=J_bridge)
    L = build_L(H, gamma)

    t0 = _time.time()

    log(f"{'B State':>8}  {'B:CΨ0':>6}  {'A:max':>6}  {'Crosses':>8}  "
        f"{'t_peak':>6}  {'t_up':>6}  {'t_down':>7}  {'above_t':>7}  {'Class':>18}")

    for name, b_state in sender_states.items():
        psi_full = np.kron(A_state, b_state)
        rho0 = ket2dm(psi_full)

        # CΨ_0 for B
        rho_B0 = ptrace_keep(rho0, [2, 3])
        cpsi_B0 = cpsi(rho_B0)

        # Euler evolution
        traj = euler_evolve(L, rho0.flatten(), dt, n_steps)

        # Extract CΨ_A trajectory
        cpsi_A = np.empty(n_steps + 1)
        for i in range(0, n_steps + 1, 10):  # sample every 10 steps for speed
            rho_t = traj[i].reshape(d, d)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_A = ptrace_keep(rho_t, [0, 1])
            cpsi_A[i] = cpsi(rho_A)
        # Interpolate skipped steps
        for i in range(n_steps + 1):
            if i % 10 != 0:
                lo = (i // 10) * 10
                hi = min(lo + 10, n_steps)
                if hi > lo:
                    frac = (i - lo) / (hi - lo)
                    cpsi_A[i] = cpsi_A[lo] + frac * (cpsi_A[hi] - cpsi_A[lo])
                else:
                    cpsi_A[i] = cpsi_A[lo]

        fp = extract_fingerprint(cpsi_A, times, dt)

        # Classify
        if fp["crosses"]:
            cls = "Local coherence"
        elif cpsi_B0 > 0.3:
            cls = "Entangled" if "Bell" in name or "Psi" in name else "Phase-sensitive"
        else:
            cls = "Classical"

        t_up_str = f"{fp['t_up']:.2f}" if fp['t_up'] else "  -"
        t_down_str = f"{fp['t_down']:.2f}" if fp['t_down'] else "   -"

        log(f"{name:>8}  {cpsi_B0:6.3f}  {fp['peak']:6.3f}  "
            f"{'YES' if fp['crosses'] else 'NEVER':>8}  "
            f"{fp['t_peak']:6.2f}  {t_up_str:>6}  {t_down_str:>7}  "
            f"{fp['above_time']:7.2f}  {cls:>18}")

    elapsed = _time.time() - t0
    log(f"  [{elapsed:.1f}s]")
    log()

# ============================================================
# ENTANGLEMENT BARRIER TEST
# ============================================================
log("Entanglement barrier ratio (product / entangled peak CΨ_A)")
log("-" * 78)

for J_bridge in [0.5, 1.0, 1.5]:
    H = build_H(J_internal=1.0, J_bridge=J_bridge)
    L = build_L(H, gamma)

    peaks = {}
    for name in ["|++>", "Bell+"]:
        psi_full = np.kron(A_state, sender_states[name])
        rho0 = ket2dm(psi_full)
        traj = euler_evolve(L, rho0.flatten(), dt, n_steps)
        max_cpsi = 0.0
        for i in range(0, n_steps + 1, 10):
            rho_t = traj[i].reshape(d, d)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_A = ptrace_keep(rho_t, [0, 1])
            c = cpsi(rho_A)
            max_cpsi = max(max_cpsi, c)
        peaks[name] = max_cpsi

    ratio = peaks["|++>"] / peaks["Bell+"] if peaks["Bell+"] > 0 else float('inf')
    log(f"  J/γ = {J_bridge/gamma:5.0f}:  |++> peak = {peaks['|++>']:.3f},  "
        f"Bell+ peak = {peaks['Bell+']:.3f},  ratio = {ratio:.1f}×")

log()
log("=" * 78)
log("Done.")
log("=" * 78)

_outf.close()
