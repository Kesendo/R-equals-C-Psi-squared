"""
Information Channel: CΨ crossing time as readout for weak coupling

Protocol: A and B share N product-state pairs with J > 0.
B measures ("1") or not ("0"). A reads crossing time shift.

At J=0.01, Δt = -0.218 (2.87% shift). 21 pairs for 1 bit
at 100% timing jitter. Channel capacity ~ (J/γ)².

Source: Thomas Wicht + Claude (Anthropic), 2026-03-01
See: experiments/OBSERVER_GRAVITY_BRIDGE.md §8
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
n_steps = 4000
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

def find_crossing(result, times):
    for i in range(1, len(times)):
        cpsi = local_cpsi_A(result.states[i])
        cpsi_prev = local_cpsi_A(result.states[i-1])
        if cpsi < 0.25 and cpsi_prev >= 0.25:
            frac = (cpsi_prev - 0.25) / (cpsi_prev - cpsi)
            return times[i-1] + frac * (times[i] - times[i-1])
    return None

t_B = 1.0
idx_tB = np.argmin(np.abs(times - t_B))

print(f"{'J':>8} {'t(0)':>10} {'t(1)':>10} {'Δt':>10} {'%':>8} {'N_min(σ=1)':>12}")
print("-" * 58)

for J in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
    H = J * (tensor(sigmax(), sigmax()) +
             tensor(sigmay(), sigmay()) +
             tensor(sigmaz(), sigmaz()))
    r_pre = mesolve(H, ket2dm(tensor(plus, plus)), times[:idx_tB+1], c_ops, [])
    rho_tB = r_pre.states[-1]

    r0 = mesolve(H, rho_tB, times, c_ops, [])
    rho_Bm = P0_B * rho_tB * P0_B.dag() + P1_B * rho_tB * P1_B.dag()
    r1 = mesolve(H, rho_Bm, times, c_ops, [])

    t0 = find_crossing(r0, times)
    t1 = find_crossing(r1, times)

    if t0 and t1:
        dt = t1 - t0
        N_min = (1.0 / abs(dt))**2 if abs(dt) > 0 else float('inf')
        print(f"{J:>8.3f} {t0:>10.4f} {t1:>10.4f} {dt:>10.4f} {dt/t0*100:>7.2f}% {N_min:>12.0f}")
