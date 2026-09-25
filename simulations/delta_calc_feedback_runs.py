"""
The February 7 feedback runs in the retired tool's own reading
==============================================================
experiments/SIMULATION_EVIDENCE.md §2 reports four runs of the retired delta_calc MCP tool
(its source is kept outside the repo and is transcribed here, not imported):
simulate_dynamic_lindblad with noise_type = operator_feedback. The Claude Desktop log of the
chat's MCP holds the four calls (2026-02-07, 16:57 to 16:58 UTC) with these arguments and the
default mutual-purity bridge; the transcription below reproduces their logged <O_int> tails
digit for digit (all twelve values the log shows of each call are checked). The tool

  - builds sigma_z jumps on every site, with the feedback rate
        gamma_eff = max(0, gamma_0 * (1 - kappa * <O_int>)),  O_int = X_0 X_1,  kappa = 0.5;
  - integrates by an Euler step of dt 0.01, then takes the Hermitian part, sets negative
    eigenvalues to zero and renormalizes the trace, up to t_max = 5;
  - records every 10th step, rounded to six digits: the purity Tr(rho^2), its mutual-purity
    bridge C = (prod_k Tr(rho_k^2))^(1/n) (the geometric mean of the single-site purities),
    and Psi = l1(rho)/(d - 1) of the full state;
  - reports C_psi = C_final * Psi_final.

Both factors come from the density matrix. This script regenerates §2's four runs (CPsi_final,
purity_final), prints each run's start, largest value and end, matches the prose ranges and §4's
Psi values against the runs' turning points, and regenerates the agents' sweep of §3 and
experiments/OPERATOR_FEEDBACK.md §4: the tool's sweep routine wraps the same routine, run with its
concurrence bridge to t_max 10 (the settings under which the column regenerates), and prints
C_final * psi_approx, a fixed 0.27; read with the density matrix's own Psi (concurrence * l1/3)
the same runs stay above 1/4 at every Euler step.

It also regenerates the operator-feedback run of hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md
§3.1 and Appendix A.2: GHZ_3 on the ring at h = 0 (the tool's default), gamma_0 = 0.05, t_max 10,
read with the tool's correlation bridge C = min(1, 2 * (P - prod_k P_k)); the tool's delta is the
recorded purity minus its dephasing prediction of the initial state, diag + offdiag *
exp(-2 * n * gamma_0 * t). The feedback of that run is blind: X_0 X_1 has no matrix element
between the GHZ branches |000> and |111>, and H at h = 0 keeps the state on them.

Conventions (the tool's): H = J * sum over the listed bonds of (XX + YY + ZZ) + h * sum_k X_k
(heisenberg = open chain, heisenberg_ring = ring); Bell+ = (|00> + |11>)/sqrt(2); GHZ_3 and W_3
with qubit 0 the most significant bit.

Import-inert; prints only.
"""

from functools import reduce
import sys

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)


def site(P, k, n):
    return reduce(np.kron, [P if i == k else I2 for i in range(n)])


def hamiltonian(n, bonds, J, h):
    H = np.zeros((2**n, 2**n), dtype=complex)
    for i, j in bonds:
        for P in (X, Y, Z):
            H += J * site(P, i, n) @ site(P, j, n)
    if h != 0:
        for i in range(n):
            H += h * site(X, i, n)
    return H


def one_site(rho, n, k):
    t = rho.reshape([2] * (2 * n))
    for i in sorted((q for q in range(n) if q != k), reverse=True):
        t = np.trace(t, axis1=i, axis2=i + t.ndim // 2)
    return t.reshape(2, 2)


def purity(r):
    return float(np.real(np.trace(r @ r)))


def l1(r):
    return float(np.abs(r).sum() - np.abs(np.diag(r)).sum())


YY = np.kron(Y, Y)


def concurrence(rho):
    """The tool's concurrence bridge (two qubits)."""
    R = rho @ YY @ rho.conj() @ YY
    ev = np.sort(np.real(np.sqrt(np.maximum(0, np.abs(np.linalg.eigvals(R))))))[::-1]
    return float(max(0.0, ev[0] - ev[1] - ev[2] - ev[3]))


def bridge_value(rho, n, bridge):
    if bridge == "concurrence":
        return concurrence(rho)
    if bridge == "correlation":
        # the tool's correlation bridge: min(1, 2 * (P - prod_k P_k)), floored at 0
        product = float(np.prod([purity(one_site(rho, n, k)) for k in range(n)]))
        return float(max(0.0, min(1.0, 2 * (purity(rho) - product))))
    return float(np.prod([purity(one_site(rho, n, k)) for k in range(n)]) ** (1.0 / n))


def feedback_run(rho0, n, H, gamma0, kappa=0.5, dt=0.01, t_max=5.0, bridge="mutual_purity", every_step=False):
    """The tool's dynamic_lindblad_evolution for operator_feedback with sigma_z jumps, transcribed.
    Returns the recorded rows (t, purity, C, Psi, <O_int>), rounded as the tool rounds them; with
    every_step=True also the lowest C * Psi over every Euler step, not only the recorded ones."""
    jumps = [site(Z, k, n) for k in range(n)]
    O_int = site(X, 0, n) @ site(X, 1, n)
    d = 2**n
    rho = rho0.copy()
    n_steps = int(t_max / dt)
    rows = []
    lowest = np.inf
    for step in range(n_steps + 1):
        if every_step:
            lowest = min(lowest, bridge_value(rho, n, bridge) * l1(rho) / (d - 1))
        if step % max(1, int(0.1 / dt)) == 0:
            C = bridge_value(rho, n, bridge)
            o = float(np.clip(np.real(np.trace(rho @ O_int)), -1.0, 1.0))
            rows.append((round(step * dt, 4), round(purity(rho), 6), round(C, 6), round(l1(rho) / (d - 1), 6),
                         round(o, 6)))
        if step == n_steps:
            break
        expectation = float(np.clip(np.real(np.trace(rho @ O_int)), -1.0, 1.0))
        gamma = max(0.0, gamma0 * (1 - kappa * expectation))
        drho = -1j * (H @ rho - rho @ H)
        for L in jumps:
            LdL = L.conj().T @ L
            drho += gamma * (L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL))
        rho = rho + dt * drho
        rho = 0.5 * (rho + rho.conj().T)
        w, v = np.linalg.eigh(rho)
        rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
        trace = np.real(np.trace(rho))
        if trace > 0:
            rho /= trace
    return (rows, lowest) if every_step else rows


def pure(amplitudes):
    v = np.asarray(amplitudes, dtype=complex)
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def turning_indices(values):
    """Indices of the start, every local peak and trough, and the end of a recorded series."""
    indices = [0]
    for i in range(1, len(values) - 1):
        rising = values[i] > values[i - 1] and values[i] >= values[i + 1]
        falling = values[i] < values[i - 1] and values[i] <= values[i + 1]
        if rising or falling:
            indices.append(i)
    return indices + [len(values) - 1]


def turning_points(values):
    """Start, every local peak and trough, and end of a recorded series."""
    return [values[i] for i in turning_indices(values)]


VALUES = []
LOGGED = []
# The Claude Desktop log truncates each response and shows the last twelve <O_int> records
# (t = 3.9 to 5.0) of each call.
LOG_TAILS = {
    "Bell+, heisenberg, J 1, h 0.9": [0.963355, 0.962475, 0.961597, 0.960721, 0.959849, 0.958979,
                                      0.958112, 0.957248, 0.956387, 0.955528, 0.954672, 0.953819],
    "Bell+, J 0, h 0": [0.96102, 0.960021, 0.959023, 0.958025, 0.957026, 0.956028,
                        0.95503, 0.954032, 0.953034, 0.952037, 0.951039, 0.950042],
    "GHZ3, heisenberg_ring, J 1, h 0.9": [0.056607, 0.057887, 0.059159, 0.060421, 0.061675, 0.06292,
                                          0.064158, 0.065388, 0.06661, 0.067825, 0.069032, 0.070231],
    "W3, heisenberg_ring, J 1, h 0.9": [0.692343, 0.692967, 0.693582, 0.694188, 0.694787, 0.695379,
                                        0.695968, 0.696554, 0.697139, 0.697726, 0.698315, 0.698908],
}
CLAIMS = []


def claim(label, ok, detail):
    print(f"  {label:<44s} {detail}  {'holds' if ok else 'FAILS'}")
    CLAIMS.append((label, ok))


def check_among(label, page, candidates, digits=2):
    """A printed value that names one of a run's turning points."""
    ok = any(round(c, digits) == round(page, digits) for c in candidates)
    near = min(candidates, key=lambda c: abs(c - page))
    print(f"  {label:<44s} page {page:.{digits}f}  nearest turning point {near:.{digits + 3}f}"
          f"  {'match' if ok else 'DIFFERS'}")
    VALUES.append((label, ok))


def check(label, page, value, digits=3):
    ok = round(value, digits) == round(page, digits)
    print(f"  {label:<44s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")
    VALUES.append((label, ok))


def main():
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    bell = pure([1, 0, 0, 1])
    ghz3 = pure([1, 0, 0, 0, 0, 0, 0, 1])
    w3 = pure([0, 1, 1, 0, 1, 0, 0, 0])
    ring3 = [(0, 1), (1, 2), (2, 0)]
    runs = {
        "Bell+, heisenberg, J 1, h 0.9": (bell, 2, hamiltonian(2, [(0, 1)], 1.0, 0.9), 0.405, 0.955),
        "Bell+, J 0, h 0": (bell, 2, hamiltonian(2, [(0, 1)], 0.0, 0.0), 0.158, 0.951),
        "GHZ3, heisenberg_ring, J 1, h 0.9": (ghz3, 3, hamiltonian(3, ring3, 1.0, 0.9), 0.262, 0.875),
        "W3, heisenberg_ring, J 1, h 0.9": (w3, 3, hamiltonian(3, ring3, 1.0, 0.9), 0.413, 0.923),
    }
    series = {}
    print("SIMULATION_EVIDENCE §2, gamma_0 = 0.005, operator feedback (kappa 0.5), t_max 5")
    for name, (rho0, n, H, page_cpsi, page_purity) in runs.items():
        rows = feedback_run(rho0, n, H, 0.005)
        series[name] = rows
        _, p_final, c_final, psi_final, _ = rows[-1]
        check(f"{name}: CPsi_final", page_cpsi, c_final * psi_final)
        check(f"{name}: purity_final", page_purity, p_final)
    print("  the <O_int> records the MCP log shows of each call (t = 3.9 to 5.0):")
    for name, logged in LOG_TAILS.items():
        computed = [r[4] for r in series[name][-len(logged):]]
        hits = [a == b for a, b in zip(logged, computed)]
        print(f"    {name:<36s} {sum(hits)} of {len(hits)} match, last logged {logged[-1]:.6f}"
              f"  computed {computed[-1]:.6f}")
        LOGGED.extend((f"{name}: logged <O_int> at t = {3.9 + 0.1 * k:.1f}", ok) for k, ok in enumerate(hits))
    cpsi = {name: [r[2] * r[3] for r in rows] for name, rows in series.items()}
    idle = cpsi["Bell+, J 0, h 0"]
    claim("Bell+ idle: CPsi decays monotonically", all(b <= a for a, b in zip(idle, idle[1:])),
          f"from {idle[0]:.4f} to {idle[-1]:.4f}")
    print("  each active run's start, largest recorded value and end, and its turning points (t: CPsi):")
    for name in ("Bell+, heisenberg, J 1, h 0.9", "GHZ3, heisenberg_ring, J 1, h 0.9",
                 "W3, heisenberg_ring, J 1, h 0.9"):
        c, times = cpsi[name], [r[0] for r in series[name]]
        k = int(np.argmax(c))
        print(f"    {name.split(',')[0]:<6s} starts {c[0]:.5f}, largest {c[k]:.5f} at t = {times[k]:.1f},"
              f" ends {c[-1]:.5f}")
        print("           " + ", ".join(f"{times[i]:.1f}: {c[i]:.3f}" for i in turning_indices(c)))
    print("  the prose's ranges, against each run's turning points:")
    for name, low, high in (("Bell+, heisenberg, J 1, h 0.9", 0.17, 0.50),
                            ("GHZ3, heisenberg_ring, J 1, h 0.9", 0.07, 0.49),
                            ("W3, heisenberg_ring, J 1, h 0.9", 0.18, 0.50)):
        points = turning_points(cpsi[name])
        check_among(f"{name.split(',')[0]} CPsi range: {low:.2f}", low, points)
        check_among(f"{name.split(',')[0]} CPsi range: {high:.2f}", high, points)

    print("SIMULATION_EVIDENCE §4, Psi = l1/(d - 1) of the full state: the printed values as turning points")
    for name, printed in (("Bell+, heisenberg, J 1, h 0.9", (0.35, 0.99, 0.37, 0.81)),
                          ("W3, heisenberg_ring, J 1, h 0.9", (0.29, 0.93, 0.32, 0.90))):
        points = turning_points([r[3] for r in series[name]])
        for value in printed:
            check_among(f"{name.split(',')[0]} Psi: {value:.2f}", value, points)

    print("SIMULATION_EVIDENCE §3 and OPERATOR_FEEDBACK §4, the agents' sweep: Bell+, concurrence bridge,")
    print("t_max 10, printed as C_final * psi_approx (0.27)")
    sweep = {(0.005, 0.7): (0.909, 0.245), (0.005, 0.8): (0.912, 0.246), (0.005, 0.9): (0.914, 0.247),
             (0.005, 1.0): (0.917, 0.248), (0.006, 0.7): (0.891, 0.241), (0.006, 0.9): (0.897, 0.242)}
    lowest = 1.0
    for (g0, h), (page_c, page_cpsi) in sweep.items():
        rows, low = feedback_run(bell, 2, hamiltonian(2, [(0, 1)], 1.0, h), g0, t_max=10.0,
                                 bridge="concurrence", every_step=True)
        c_final = rows[-1][2]
        check(f"gamma_0 {g0}, h {h}: C_final", page_c, c_final)
        check(f"gamma_0 {g0}, h {h}: C_final * 0.27", page_cpsi, c_final * 0.27)
        lowest = min(lowest, low)
    claim("read as concurrence * l1/3, every run > 1/4", lowest > 0.25,
          f"lowest value over every Euler step, t <= 10: {lowest:.4f}")
    rows = feedback_run(bell, 2, hamiltonian(2, [(0, 1)], 1.0, 0.7), 0.003, t_max=10.0, bridge="concurrence")
    check("gamma_0 0.003, h 0.7: C_final * 0.27", 0.255, rows[-1][2] * 0.27)
    check("the threshold 25/27, printed 0.926", 0.926, 25 / 27)

    print("MEDIATOR_AS_QUANTUM_TRANSISTOR §3.1 and A.2: GHZ3, heisenberg_ring, J 1, h 0, gamma_0 = 0.05,")
    print("operator feedback (kappa 0.5), correlation bridge, t_max 10")
    ring_h0 = hamiltonian(3, ring3, 1.0, 0.0)
    rows = feedback_run(ghz3, 3, ring_h0, 0.05, t_max=10.0, bridge="correlation")
    unfed = feedback_run(ghz3, 3, ring_h0, 0.05, kappa=0.0, t_max=10.0, bridge="correlation")
    # The Z jumps are diagonal and H at h = 0 maps |000> and |111> to 3 times themselves, so rho stays
    # on the two branches; X_0 X_1 flips two bits and has no matrix element between them.
    O_int = site(X, 0, 3) @ site(X, 1, 3)
    branches = (0, 7)
    kept = all(ring_h0[i, j] == (3.0 if i == j else 0.0) for j in branches for i in range(8))
    unseen = all(O_int[i, j] == 0.0 for i in branches for j in branches)
    claim("<X0 X1> = 0 on the GHZ branches, which H keeps", kept and unseen,
          "exact matrix entries")
    claim("the feedback changes no recorded value", rows == unfed,
          f"{len(rows)} records at kappa 0.5 and at kappa 0, largest |<O_int>| {max(abs(r[4]) for r in rows)}")
    times, P, C = [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows]
    diag = float(np.sum(np.abs(np.diag(ghz3)) ** 2))
    off = purity(ghz3) - diag
    delta = [round(p - (diag + off * np.exp(-2 * 3 * 0.05 * s)), 6) for s, p in zip(times, P)]
    k = int(np.argmin(delta))
    check("A.2: purity_final", 0.501, P[-1])
    check("A.2: C_final", 0.752, C[-1])
    check("§3.1: C_final, printed 0.75", 0.75, C[-1], digits=2)
    check("C at its ceiling 1.0 up to t", 2.3, max(s for s, c in zip(times, C) if c == 1.0), digits=1)
    check("A.2: delta minimum", -0.125, delta[k])
    check("A.2: time of the delta minimum", 2.3, times[k], digits=1)
    print(f"  closed form, plain dephasing: P = 1/2 + exp(-0.6 t)/2 = {0.5 + 0.5 * np.exp(-6.0):.6f} at t = 10"
          f" (recorded {P[-1]:.6f}); C leaves 1.0 at P = 5/8, t = ln 4/0.6 = {np.log(4) / 0.6:.4f};")
    print(f"  delta = (exp(-0.6 t) - exp(-0.3 t))/2 has its minimum -1/8 at t = ln 2/0.3 = {np.log(2) / 0.3:.4f}"
          f" (recorded {delta[k]:.6f} at t = {times[k]:.1f})")

    matched = sum(ok for _, ok in VALUES)
    logged = sum(ok for _, ok in LOGGED)
    held = sum(ok for _, ok in CLAIMS)
    print(f"\n{matched} of {len(VALUES)} printed values match at the page's precision;"
          f" {logged} of {len(LOGGED)} logged values; {held} of {len(CLAIMS)} claims hold.")
    for label, ok in VALUES + LOGGED + CLAIMS:
        if not ok:
            print(f"  not reproduced: {label}")


if __name__ == "__main__":
    main()
