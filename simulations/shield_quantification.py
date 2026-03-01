"""
Shield Quantification: How coupling affects coherence lifetime

CORRECTED FINDINGS:
- Coupling does NOT protect A. It redistributes coherence into
  nonlocal pool, making A decay FASTER locally vs isolation.
- |+,0> with J=1: t_cross = 0.64 vs single |+>: t_cross = 8.58
- B's measurement cuts the return flow of the oscillation,
  causing additional acceleration within the coupled system.
- Damage is timing-dependent: max at t_B ≈ 1.0 (oscillation phase).
- |+,+> under Heisenberg: coupling has NO effect (symmetry sector).

Source: Thomas Wicht + Claude (Anthropic), 2026-03-01
See: experiments/OBSERVER_GRAVITY_BRIDGE.md §7
"""

import numpy as np
from qutip import (
    basis, tensor, ket2dm, sigmax, sigmay, sigmaz, qeye,
    mesolve, concurrence
)

zero = basis(2, 0)
one  = basis(2, 1)
plus = (zero + one).unit()

gamma = 0.05
times = np.linspace(0, 50, 5000)

def local_cpsi_A(rho):
    rho_A = rho.ptrace(0)
    purity = (rho_A * rho_A).tr().real
    rho_full = rho_A.full()
    l1 = abs(rho_full[0,1]) + abs(rho_full[1,0])
    return purity * l1

def local_cpsi_single(rho):
    purity = (rho * rho).tr().real
    rho_full = rho.full()
    l1 = abs(rho_full[0,1]) + abs(rho_full[1,0])
    return purity * l1

def find_crossing(states, times, func):
    for i in range(1, len(times)):
        c, c_p = func(states[i]), func(states[i-1])
        if c < 0.25 and c_p >= 0.25:
            f = (c_p - 0.25) / (c_p - c)
            return times[i-1] + f * (times[i] - times[i-1])
    return None

P0_B = tensor(qeye(2), zero * zero.dag())
P1_B = tensor(qeye(2), one * one.dag())
c_ops = [np.sqrt(gamma)*tensor(sigmaz(),qeye(2)),
         np.sqrt(gamma)*tensor(qeye(2),sigmaz())]

# Single qubit baseline
r_s = mesolve(0*sigmax(), ket2dm(plus), times,
              [np.sqrt(gamma)*sigmaz()], [])
t_single = find_crossing(r_s.states, times, local_cpsi_single)
print(f"Single |+>: t_cross = {t_single:.4f}")

# |+,0> with varying J
print(f"\n|+,0> coupled system:")
print(f"{'J':>6} {'t_cross':>10} {'vs single':>10}")
for J in [0, 0.05, 0.1, 0.5, 1.0, 2.0]:
    H = J*(tensor(sigmax(),sigmax())+tensor(sigmay(),sigmay())+tensor(sigmaz(),sigmaz()))
    r = mesolve(H, ket2dm(tensor(plus,zero)), times, c_ops, [])
    tc = find_crossing(r.states, times, local_cpsi_A)
    if tc: print(f"{J:>6.2f} {tc:>10.4f} {tc/t_single:>10.2f}x")

# Timing dependence
print(f"\nTiming dependence (|+,0>, J=1.0):")
H = 1.0*(tensor(sigmax(),sigmax())+tensor(sigmay(),sigmay())+tensor(sigmaz(),sigmaz()))
r_nat = mesolve(H, ket2dm(tensor(plus,zero)), times, c_ops, [])
t_nat = find_crossing(r_nat.states, times, local_cpsi_A)
print(f"Natural: {t_nat:.4f}")
print(f"{'t_B':>6} {'t_cross':>10} {'% remaining':>12}")
for t_B in [0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0]:
    idx = np.argmin(np.abs(times - t_B))
    r_pre = mesolve(H, ket2dm(tensor(plus,zero)), times[:idx+1], c_ops, [])
    rho = r_pre.states[-1]
    rho_m = P0_B*rho*P0_B.dag() + P1_B*rho*P1_B.dag()
    r_post = mesolve(H, rho_m, times, c_ops, [])
    tc = find_crossing(r_post.states, times, local_cpsi_A)
    if tc and t_nat:
        print(f"{t_B:>6.2f} {tc:>10.4f} {tc/t_nat*100:>11.1f}%")
