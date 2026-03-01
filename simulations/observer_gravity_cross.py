"""
Observer × Gravity Crossing Matrix

Computes t_cross(Observer, Gravity) = K(Observer, State) / γ(Gravity)

Findings:
  - K is γ-invariant (CV = 0.00%)
  - K(Conc)/K(MI) = 1.2125 across all gravitational environments
  - K ratio is STATE-DEPENDENT (CV = 13.5% across initial states)
  - States with α < 30° never cross — no observer time

Source: Thomas Wicht + Claude (Anthropic), 2026-03-01
See: experiments/OBSERVER_GRAVITY_BRIDGE.md
"""

import numpy as np
from qutip import (
    basis, tensor, ket2dm, sigmax, sigmay, sigmaz, qeye,
    mesolve, entropy_vn, concurrence
)

zero = basis(2, 0)
one = basis(2, 1)
bell = (tensor(zero, zero) + tensor(one, one)).unit()

J = 1.0
H = J * (tensor(sigmax(), sigmax()) +
         tensor(sigmay(), sigmay()) +
         tensor(sigmaz(), sigmaz()))

t_max = 20.0
n_steps = 2000
times = np.linspace(0, t_max, n_steps)

gammas = {
    "Deep Space": 0.01, "Mars": 0.019, "Earth": 0.05,
    "Jupiter": 0.13, "Neutron": 0.20, "Black Hole": 0.50,
}

def c_concurrence(rho):
    try: return concurrence(rho)
    except: return 0.0

def c_mutual_info(rho):
    rho_A, rho_B = rho.ptrace(0), rho.ptrace(1)
    return (entropy_vn(rho_A, 2) + entropy_vn(rho_B, 2) - entropy_vn(rho, 2)) / 2.0

def psi_l1(rho, d=4):
    rho_full = rho.full()
    l1 = sum(abs(rho_full[i, j]) for i in range(d) for j in range(d) if i != j)
    return l1 / (d - 1)

metrics = {"Concurrence": c_concurrence, "Mutual Info": c_mutual_info}

def find_crossing(gamma_val, metric_func, H, times):
    c_ops = [np.sqrt(gamma_val) * tensor(sigmaz(), qeye(2)),
             np.sqrt(gamma_val) * tensor(qeye(2), sigmaz())]
    result = mesolve(H, ket2dm(bell), times, c_ops, [])
    for i in range(1, len(times)):
        C = metric_func(result.states[i])
        Psi = psi_l1(result.states[i])
        cpsi = C * Psi
        C_prev = metric_func(result.states[i-1])
        Psi_prev = psi_l1(result.states[i-1])
        cpsi_prev = C_prev * Psi_prev
        if cpsi < 0.25 and cpsi_prev >= 0.25:
            frac = (cpsi_prev - 0.25) / (cpsi_prev - cpsi)
            return times[i-1] + frac * (times[i] - times[i-1])
    return None

# Part 1: Crossing matrix
print("OBSERVER × GRAVITY CROSSING MATRIX")
print(f"{'':>12}", end="")
for m in metrics: print(f"  {m:>14}", end="")
print()

K_matrix = {}
for env, g in gammas.items():
    print(f"{env:>12}", end="")
    for mname, mfunc in metrics.items():
        tc = find_crossing(g, mfunc, H, times)
        if tc:
            K = g * tc
            K_matrix[(env, mname)] = K
            print(f"  {tc:>14.4f}", end="")
        else:
            print(f"  {'never':>14}", end="")
    print()

# Part 2: K invariance
print(f"\nK = γ·t_cross:")
for env in gammas:
    print(f"  {env:>12}", end="")
    for mname in metrics:
        K = K_matrix.get((env, mname))
        print(f"  {K:.6f}" if K else "  —", end="")
    print()

# Part 3: Ratio
print(f"\nK(Conc)/K(MI):")
ratios = []
for env in gammas:
    Kc = K_matrix.get((env, "Concurrence"))
    Km = K_matrix.get((env, "Mutual Info"))
    if Kc and Km:
        r = Kc / Km
        ratios.append(r)
        print(f"  {env:>12}: {r:.6f}")
print(f"  Mean: {np.mean(ratios):.6f} ± {np.std(ratios):.6f}")

# Part 4: State dependence
print(f"\nState-dependent ratio (γ=0.05):")
gamma_fixed = 0.05
c_ops_fixed = [np.sqrt(gamma_fixed) * tensor(sigmaz(), qeye(2)),
               np.sqrt(gamma_fixed) * tensor(qeye(2), sigmaz())]

for alpha_deg in [45, 40, 35, 30]:
    alpha = np.radians(alpha_deg)
    psi = (np.cos(alpha)*tensor(zero,zero) + np.sin(alpha)*tensor(one,one)).unit()
    r = mesolve(H, ket2dm(psi), times, c_ops_fixed, [])
    tc_c, tc_m = None, None
    for i in range(1, len(times)):
        Cc = c_concurrence(r.states[i])
        Cm = c_mutual_info(r.states[i])
        Psi = psi_l1(r.states[i])
        if tc_c is None and Cc*Psi < 0.25:
            Cc_p = c_concurrence(r.states[i-1])
            Psi_p = psi_l1(r.states[i-1])
            if Cc_p*Psi_p >= 0.25:
                f = (Cc_p*Psi_p - 0.25) / (Cc_p*Psi_p - Cc*Psi)
                tc_c = times[i-1] + f*(times[i]-times[i-1])
        if tc_m is None and Cm*Psi < 0.25:
            Cm_p = c_mutual_info(r.states[i-1])
            Psi_p = psi_l1(r.states[i-1])
            if Cm_p*Psi_p >= 0.25:
                f = (Cm_p*Psi_p - 0.25) / (Cm_p*Psi_p - Cm*Psi)
                tc_m = times[i-1] + f*(times[i]-times[i-1])
    if tc_c and tc_m:
        print(f"  α={alpha_deg}°: K_c/K_m = {(gamma_fixed*tc_c)/(gamma_fixed*tc_m):.6f}")
    else:
        print(f"  α={alpha_deg}°: never crosses")
