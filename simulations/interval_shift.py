"""
Interval Shift: B's Measurement Shifts A's Local Crossing Time

Key result: The shift Δt is continuous in J with NO threshold.
Any J > 0 produces a measurable interval shift.

Product states required (Bell+ has no local clock).
Bell+ with any J: rho_A = I/2, never crosses locally.

Source: Thomas Wicht + Claude (Anthropic), 2026-03-01
See: experiments/OBSERVER_GRAVITY_BRIDGE.md
"""

import numpy as np
from qutip import (
    basis, tensor, ket2dm, sigmax, sigmay, sigmaz, qeye, mesolve
)

zero = basis(2, 0)
one  = basis(2, 1)
plus = (zero + one).unit()

gamma = 0.05
t_max = 20.0
n_steps = 2000
times = np.linspace(0, t_max, n_steps)

c_ops = [np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
         np.sqrt(gamma) * tensor(qeye(2), sigmaz())]

P0_B = tensor(qeye(2), zero * zero.dag())
P1_B = tensor(qeye(2), one * one.dag())

def local_cpsi_A(rho):
    rho_A = rho.ptrace(0)
    purity = (rho_A * rho_A).tr().real
    rho_full = rho_A.full()
    l1 = abs(rho_full[0,1]) + abs(rho_full[1,0])
    return purity * l1

def find_local_crossing_A(result, times):
    for i in range(1, len(times)):
        cpsi = local_cpsi_A(result.states[i])
        cpsi_prev = local_cpsi_A(result.states[i-1])
        if cpsi < 0.25 and cpsi_prev >= 0.25:
            frac = (cpsi_prev - 0.25) / (cpsi_prev - cpsi)
            return times[i-1] + frac * (times[i] - times[i-1])
    return None

psi_pp = tensor(plus, plus)
t_B = 1.0
idx_tB = np.argmin(np.abs(times - t_B))

print("J-sweep: interval shift from B's measurement")
print(f"State: |++>, γ={gamma}, B measures Z at t={t_B}")
print(f"\n{'J':>8} {'t(nothing)':>12} {'t(B meas)':>12} {'Δt':>10} {'%':>8}")
print("-" * 52)

for J_val in [0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
    H = J_val * (tensor(sigmax(), sigmax()) +
                  tensor(sigmay(), sigmay()) +
                  tensor(sigmaz(), sigmaz()))

    r_pre = mesolve(H, ket2dm(psi_pp), times[:idx_tB+1], c_ops, [])
    rho_tB = r_pre.states[-1]

    r_n = mesolve(H, rho_tB, times, c_ops, [])
    rho_Bm = P0_B * rho_tB * P0_B.dag() + P1_B * rho_tB * P1_B.dag()
    r_m = mesolve(H, rho_Bm, times, c_ops, [])

    tc_n = find_local_crossing_A(r_n, times)
    tc_m = find_local_crossing_A(r_m, times)

    if tc_n and tc_m:
        dt = tc_m - tc_n
        print(f"{J_val:>8.3f} {tc_n:>12.4f} {tc_m:>12.4f} {dt:>10.4f} {dt/tc_n*100:>7.2f}%")
    else:
        tn = f"{tc_n:.4f}" if tc_n else "never"
        tm = f"{tc_m:.4f}" if tc_m else "never"
        print(f"{J_val:>8.3f} {tn:>12} {tm:>12}")
