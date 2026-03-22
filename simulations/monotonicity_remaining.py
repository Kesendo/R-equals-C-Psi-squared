#!/usr/bin/env python3
"""
Monotonicity: Remaining Gaps
==============================
Test A: General initial states — is the CΨ ENVELOPE monotonic for |01⟩ etc.?
Test B: Non-local (collective) noise — L = √γ(Z₁+Z₂)/√2
Test C: N > 2 subsystem CΨ — do all 2-qubit subsystems cross 1/4?

Script:  simulations/monotonicity_remaining.py
Output:  simulations/results/monotonicity_remaining.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "monotonicity_remaining.txt")
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


def site_op(op, k, nq):
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_liouvillian(H, collapse_ops):
    d = H.shape[0]
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for C_op in collapse_ops:
        Cd = C_op.conj().T
        CdC = Cd @ C_op
        L += np.kron(C_op, C_op.conj()) - 0.5 * (
            np.kron(CdC, Id) + np.kron(Id, CdC.T))
    return L


def evolve(L, rho0, t):
    d = rho0.shape[0]
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def partial_trace_keep(rho, keep, n_qubits):
    keep = list(keep)
    trace_out = [q for q in range(n_qubits) if q not in keep]
    dims = [2] * n_qubits
    reshaped = rho.reshape(dims + dims)
    current_n = n_qubits
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_keep = 2 ** len(keep)
    return reshaped.reshape((d_keep, d_keep))


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1) if d > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def heisenberg_H(nq, J=1.0):
    d = 2 ** nq
    H = np.zeros((d, d), dtype=complex)
    for i in range(nq - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, nq) @ site_op(P, i + 1, nq)
    return H


def envelope_monotonic(cpsi_arr, threshold=0.25):
    """Check if the LOCAL MAXIMA of CΨ above threshold are monotonically decreasing."""
    peaks = []
    for i in range(1, len(cpsi_arr) - 1):
        if cpsi_arr[i] > cpsi_arr[i-1] and cpsi_arr[i] >= cpsi_arr[i+1]:
            if cpsi_arr[i] > threshold * 0.5:  # only peaks near threshold
                peaks.append(cpsi_arr[i])
    if len(peaks) < 2:
        return True, peaks
    mono = all(peaks[i] >= peaks[i+1] for i in range(len(peaks)-1))
    return mono, peaks


# ====================================================================
# Test A: General initial states (2 qubits)
# ====================================================================

def test_A():
    log("=" * 70)
    log("TEST A: GENERAL INITIAL STATES (2 qubits)")
    log("=" * 70)
    log()

    nq = 2; d = 4
    J = 1.0; gamma = 0.05
    H = heisenberg_H(nq, J)
    c_ops = [np.sqrt(gamma) * site_op(sz, k, nq) for k in range(nq)]
    L = build_liouvillian(H, c_ops)

    t_max = 80.0; dt = 0.02
    n_steps = int(round(t_max / dt))
    t_arr = np.array([i * dt for i in range(n_steps + 1)])

    # Named states
    states = {}
    # Bell states
    for name, s in [("Bell+", (np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2)),
                    ("Bell-", (np.kron(up,up)-np.kron(dn,dn))/np.sqrt(2)),
                    ("Psi+",  (np.kron(up,dn)+np.kron(dn,up))/np.sqrt(2)),
                    ("Psi-",  (np.kron(up,dn)-np.kron(dn,up))/np.sqrt(2))]:
        states[name] = np.outer(s, s.conj())
    # Product states
    for name, s in [("|01>", np.kron(up,dn)), ("|10>", np.kron(dn,up)),
                    ("|11>", np.kron(dn,dn)), ("|+,+>", np.kron(plus,plus)),
                    ("|+,->", np.kron(plus,(up-dn)/np.sqrt(2)))]:
        states[name] = np.outer(s, s.conj())
    # Random states (10 Haar-random)
    np.random.seed(42)
    for i in range(10):
        v = np.random.randn(d) + 1j * np.random.randn(d)
        v /= np.linalg.norm(v)
        states[f"Rand{i}"] = np.outer(v, v.conj())

    log(f"  J={J}, gamma={gamma} (Z-dephasing), t_max={t_max}")
    log()
    log(f"  {'State':>10}  {'CPsi0':>7}  {'t_cross':>8}  {'#Osc':>5}  {'EnvMono':>7}  "
        f"{'#Peaks':>6}  {'Final':>8}")
    log("  " + "-" * 62)

    all_cross = True
    all_env_mono = True

    for name, rho0 in states.items():
        cpsi_arr = np.zeros(n_steps + 1)
        for i in range(n_steps + 1):
            rho = evolve(L, rho0, t_arr[i])
            cpsi_arr[i] = cpsi(rho)

        # Find crossing
        t_cross = None
        n_osc = 0  # oscillations above 1/4 (increases)
        for i in range(1, len(t_arr)):
            if cpsi_arr[i-1] > 0.25 and cpsi_arr[i] <= 0.25 and t_cross is None:
                t_cross = t_arr[i]
            if cpsi_arr[i] > cpsi_arr[i-1] and cpsi_arr[i-1] > 0.25:
                n_osc += 1

        env_mono, peaks = envelope_monotonic(cpsi_arr)

        if cpsi_arr[0] > 0.25 and t_cross is None:
            all_cross = False
        if not env_mono:
            all_env_mono = False

        tc = f"{t_cross:.2f}" if t_cross else ("N/A" if cpsi_arr[0] < 0.25 else "NEVER")
        em = "YES" if env_mono else "NO"
        log(f"  {name:>10}  {cpsi_arr[0]:7.4f}  {tc:>8}  {n_osc:5d}  {em:>7}  "
            f"{len(peaks):6d}  {cpsi_arr[-1]:8.6f}")

    log()
    log(f"  All states starting above 1/4 cross: {all_cross}")
    log(f"  All envelopes monotonic: {all_env_mono}")
    log()


# ====================================================================
# Test B: Non-local (collective) noise
# ====================================================================

def test_B():
    log("=" * 70)
    log("TEST B: NON-LOCAL (COLLECTIVE) NOISE")
    log("=" * 70)
    log()

    nq = 2; d = 4
    J = 1.0
    H = heisenberg_H(nq, J)

    bell_plus = (np.kron(up,up) + np.kron(dn,dn)) / np.sqrt(2)
    rho0 = np.outer(bell_plus, bell_plus.conj())

    t_max = 60.0; dt = 0.02
    n_steps = int(round(t_max / dt))
    t_arr = np.array([i * dt for i in range(n_steps + 1)])

    noise_configs = [
        ("Local Z",           lambda g: [np.sqrt(g)*site_op(sz,k,nq) for k in range(nq)]),
        ("Collective Z",      lambda g: [np.sqrt(g)*(site_op(sz,0,nq)+site_op(sz,1,nq))/np.sqrt(2)]),
        ("Local X",           lambda g: [np.sqrt(g)*site_op(sx,k,nq) for k in range(nq)]),
        ("Collective X",      lambda g: [np.sqrt(g)*(site_op(sx,0,nq)+site_op(sx,1,nq))/np.sqrt(2)]),
        ("Collective Z+X",    lambda g: [np.sqrt(g/2)*(site_op(sz,0,nq)+site_op(sz,1,nq))/np.sqrt(2),
                                          np.sqrt(g/2)*(site_op(sx,0,nq)+site_op(sx,1,nq))/np.sqrt(2)]),
        ("Anti-correlated Z", lambda g: [np.sqrt(g)*(site_op(sz,0,nq)-site_op(sz,1,nq))/np.sqrt(2)]),
    ]

    log(f"  Bell+ initial, J={J}")
    log(f"  {'Noise':>20}  {'gamma':>6}  {'t_cross':>8}  {'#Osc':>5}  {'Mono':>5}  {'Final':>8}")
    log("  " + "-" * 60)

    for gamma in [0.05, 0.1]:
        for name, ops_fn in noise_configs:
            c_ops = ops_fn(gamma)
            Liouv = build_liouvillian(H, c_ops)

            cpsi_arr = np.zeros(n_steps + 1)
            for i in range(n_steps + 1):
                rho = evolve(Liouv, rho0, t_arr[i])
                cpsi_arr[i] = cpsi(rho)

            t_cross = None
            n_osc = 0
            for i in range(1, len(t_arr)):
                if cpsi_arr[i-1] > 0.25 and cpsi_arr[i] <= 0.25 and t_cross is None:
                    t_cross = t_arr[i]
                if cpsi_arr[i] > cpsi_arr[i-1] and cpsi_arr[i-1] > 0.25:
                    n_osc += 1

            tc = f"{t_cross:.3f}" if t_cross else "NEVER"
            mono = "YES" if n_osc == 0 else f"NO({n_osc})"
            log(f"  {name:>20}  {gamma:6.2f}  {tc:>8}  {n_osc:5d}  {mono:>5}  "
                f"{cpsi_arr[-1]:8.6f}")
        log()

    log("  Key question: does collective noise change the monotonicity?")
    log()


# ====================================================================
# Test C: N > 2 subsystem CΨ
# ====================================================================

def test_C():
    log("=" * 70)
    log("TEST C: N > 2 SUBSYSTEM CΨ")
    log("=" * 70)
    log()

    gamma = 0.05; J = 1.0

    for nq in [3, 4, 5]:
        d = 2 ** nq
        t0 = _time.time()
        H = heisenberg_H(nq, J)
        c_ops = [np.sqrt(gamma) * site_op(sz, k, nq) for k in range(nq)]
        Liouv = build_liouvillian(H, c_ops)

        # GHZ initial state
        ghz = np.zeros(d, dtype=complex)
        ghz[0] = 1 / np.sqrt(2)     # |000...0>
        ghz[-1] = 1 / np.sqrt(2)    # |111...1>
        rho0 = np.outer(ghz, ghz.conj())

        # Also try W state for N=3,4
        if nq <= 4:
            w = np.zeros(d, dtype=complex)
            for k in range(nq):
                idx = 1 << (nq - 1 - k)  # single excitation at position k
                w[idx] = 1 / np.sqrt(nq)
            rho0_w = np.outer(w, w.conj())
            init_states = [("GHZ", rho0), ("W", rho0_w)]
        else:
            init_states = [("GHZ", rho0)]

        # All nearest-neighbor pairs
        pairs = [(i, i+1) for i in range(nq-1)]

        t_max = 40.0 if nq <= 4 else 20.0
        dt_sim = 0.05 if nq <= 4 else 0.1
        n_steps = int(round(t_max / dt_sim))
        t_arr = np.array([i * dt_sim for i in range(n_steps + 1)])

        for state_name, rho_init in init_states:
            log(f"  N={nq}, state={state_name}, J={J}, gamma={gamma}, d={d}")

            pair_results = {p: [] for p in pairs}

            for i in range(n_steps + 1):
                rho = evolve(Liouv, rho_init, t_arr[i])
                for p in pairs:
                    rho_pair = partial_trace_keep(rho, list(p), nq)
                    pair_results[p].append(cpsi(rho_pair))

            log(f"  {'Pair':>8}  {'CPsi0':>7}  {'t_cross':>8}  {'#Osc':>5}  "
                f"{'EnvMono':>7}  {'Final':>8}")
            log("  " + "-" * 52)

            all_cross = True
            for p in pairs:
                cpsi_arr = np.array(pair_results[p])
                t_cross = None
                n_osc = 0
                for idx in range(1, len(t_arr)):
                    if cpsi_arr[idx-1] > 0.25 and cpsi_arr[idx] <= 0.25 and t_cross is None:
                        t_cross = t_arr[idx]
                    if cpsi_arr[idx] > cpsi_arr[idx-1] and cpsi_arr[idx-1] > 0.25:
                        n_osc += 1

                env_mono, peaks = envelope_monotonic(cpsi_arr)
                tc = f"{t_cross:.2f}" if t_cross else ("N/A" if cpsi_arr[0] < 0.25 else "NEVER")
                em = "YES" if env_mono else "NO"

                if cpsi_arr[0] > 0.25 and t_cross is None:
                    all_cross = False

                log(f"  ({p[0]},{p[1]}){' ':>4}  {cpsi_arr[0]:7.4f}  {tc:>8}  "
                    f"{n_osc:5d}  {em:>7}  {cpsi_arr[-1]:8.6f}")

            elapsed = _time.time() - t0
            log(f"  All subsystems cross 1/4: {all_cross}  ({elapsed:.1f}s)")
            log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Monotonicity: Remaining Gaps")
    log("=" * 70)
    log()

    test_A()
    test_B()
    test_C()

    log("=" * 70)
    log("OVERALL VERDICT")
    log("=" * 70)
    log()
    log("A: General initial states — envelope monotonic?")
    log("B: Collective noise — still crosses 1/4?")
    log("C: N>2 subsystems — all pairs cross 1/4?")
    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
