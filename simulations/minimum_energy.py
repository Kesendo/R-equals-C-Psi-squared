"""
Minimum Energy for CΨ Crossing

Key finding: There is no energy threshold. The crossing condition
is CΨ_max ≥ ¼, determined by the competition J/γ.

Three regimes:
  1. CΨ(0) > ¼: always crosses (decay through ¼)
  2. CΨ(0) < ¼, CΨ_max > ¼: Hamiltonian pumps, then crosses
  3. CΨ_max < ¼: never crosses (no time)

Source: Thomas Wicht + Claude (Anthropic), 2026-03-01
See: experiments/MINIMUM_CROSSING_ENERGY.md
"""

import numpy as np
from qutip import (
    basis, tensor, ket2dm, sigmax, sigmay, sigmaz, qeye,
    mesolve, concurrence, expect
)

zero = basis(2, 0)
one  = basis(2, 1)

J = 1.0
H = J * (tensor(sigmax(), sigmax()) +
         tensor(sigmay(), sigmay()) +
         tensor(sigmaz(), sigmaz()))

gamma = 0.05
c_ops = [np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
         np.sqrt(gamma) * tensor(qeye(2), sigmaz())]

t_max = 100.0
n_steps = 5000
times = np.linspace(0, t_max, n_steps)

def psi_l1(rho, d=4):
    rho_full = rho.full()
    l1 = sum(abs(rho_full[i,j]) for i in range(d) for j in range(d) if i != j)
    return l1 / (d - 1)

def cpsi_conc(rho):
    try: C = concurrence(rho)
    except: C = 0.0
    return C * psi_l1(rho)

def find_crossing_and_max(result, times):
    cpsi_max, t_max_c, t_cross = 0, 0, None
    for i in range(len(times)):
        c = cpsi_conc(result.states[i])
        if c > cpsi_max:
            cpsi_max = c
            t_max_c = times[i]
        if t_cross is None and i > 0:
            c_prev = cpsi_conc(result.states[i-1])
            if c < 0.25 and c_prev >= 0.25:
                frac = (c_prev - 0.25) / (c_prev - c)
                t_cross = times[i-1] + frac * (times[i] - times[i-1])
    return cpsi_max, t_max_c, t_cross

# Part 1: α sweep for cos(α)|00⟩ + sin(α)|11⟩
print("Family: cos(α)|00⟩ + sin(α)|11⟩  (constant ⟨H⟩ = J)")
print(f"{'α':>6} {'⟨H⟩':>6} {'CΨ(0)':>8} {'Crosses?':>10}")
for alpha_deg in [45, 40, 35, 31, 30, 25, 20]:
    alpha = np.radians(alpha_deg)
    psi = (np.cos(alpha)*tensor(zero,zero) + np.sin(alpha)*tensor(one,one)).unit()
    rho0 = ket2dm(psi)
    E = expect(H, rho0)
    result = mesolve(H, rho0, times, c_ops, [])
    _, _, tc = find_crossing_and_max(result, times)
    print(f"{alpha_deg:>6}° {E:>6.1f} {cpsi_conc(rho0):>8.4f} {'YES' if tc else 'no':>10}")

# Part 2: Product states
print("\nProduct states under Heisenberg H:")
for name, psi in [
    ("|0,1>", tensor(zero, one)),
    ("|+,0>", tensor((zero+one).unit(), zero)),
    ("|+,+>", tensor((zero+one).unit(), (zero+one).unit())),
    ("|0,0>", tensor(zero, zero)),
]:
    result = mesolve(H, ket2dm(psi), times, c_ops, [])
    cm, tm, tc = find_crossing_and_max(result, times)
    print(f"  {name:>6}: CΨ_max={cm:.4f}, crosses={'YES' if tc else 'no'}")

# Part 3: J/γ sweep for |0,1⟩
print(f"\nJ/γ sweep for |0,1⟩:")
for ratio in [0.5, 1, 2, 5, 10, 20, 50]:
    g = 0.05
    j = ratio * g
    H_t = j * (tensor(sigmax(),sigmax())+tensor(sigmay(),sigmay())+tensor(sigmaz(),sigmaz()))
    c_t = [np.sqrt(g)*tensor(sigmaz(),qeye(2)), np.sqrt(g)*tensor(qeye(2),sigmaz())]
    result = mesolve(H_t, ket2dm(tensor(zero,one)), times, c_t, [])
    cm, _, tc = find_crossing_and_max(result, times)
    print(f"  J/γ={ratio:>5.1f}: CΨ_max={cm:.4f}, crosses={'YES' if tc else 'no'}")
