"""
K-Surface: K(Observer, State) Geometry

Computes K = gamma * t_cross for three observer metrics across
the full state family cos(alpha)|00> + sin(alpha)|11>.

Observers:
  - Concurrence (Type B): C and Psi both decay
  - Mutual Info (Type B): C and Psi both decay
  - Pure-Psi (Type A): C = 1.0, only Psi drives crossing
    (corresponds to "correlation" bridge from taxonomy)

The "correlation" bridge has C = 1.0 throughout the crossing
period, so CΨ = Ψ. This means t_cross is when Ψ alone = 0.25.

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

J = 1.0
H = J * (tensor(sigmax(), sigmax()) +
         tensor(sigmay(), sigmay()) +
         tensor(sigmaz(), sigmaz()))

gamma = 0.05
t_max = 30.0
n_steps = 3000
times = np.linspace(0, t_max, n_steps)

c_ops = [np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
         np.sqrt(gamma) * tensor(qeye(2), sigmaz())]

# === Observer metrics ===

def c_concurrence(rho):
    try: return concurrence(rho)
    except: return 0.0

def c_mutual_info(rho):
    rho_A = rho.ptrace(0)
    rho_B = rho.ptrace(1)
    SA = entropy_vn(rho_A, 2)
    SB = entropy_vn(rho_B, 2)
    SAB = entropy_vn(rho, 2)
    return max(0, (SA + SB - SAB) / 2.0)

def c_pure(rho):
    """Type A: C = 1.0 always. CΨ = Ψ alone."""
    return 1.0

def psi_l1(rho, d=4):
    rho_full = rho.full()
    l1 = sum(abs(rho_full[i, j]) for i in range(d) for j in range(d) if i != j)
    return l1 / (d - 1)

observers = {
    "Concurrence": c_concurrence,
    "Mutual Info": c_mutual_info,
    "Pure-Psi":    c_pure,
}

# === Find crossing time for CΨ = C * Ψ = 0.25 ===

def find_crossing(states_evolved, metric_func):
    for i in range(1, len(states_evolved)):
        C = metric_func(states_evolved[i])
        Psi = psi_l1(states_evolved[i])
        cpsi = C * Psi
        C_prev = metric_func(states_evolved[i-1])
        Psi_prev = psi_l1(states_evolved[i-1])
        cpsi_prev = C_prev * Psi_prev
        if cpsi < 0.25 and cpsi_prev >= 0.25:
            frac = (cpsi_prev - 0.25) / (cpsi_prev - cpsi)
            return times[i-1] + frac * (times[i] - times[i-1])
    return None

# === Sweep alpha for cos(α)|00> + sin(α)|11> ===

alphas_deg = np.arange(5, 91, 1)  # 5° to 90° in 1° steps
results = {name: [] for name in observers}
alpha_list = []

print("K-SURFACE: K(Observer, alpha)")
print(f"State family: cos(alpha)|00> + sin(alpha)|11>")
print(f"gamma = {gamma}, J = {J}")
print(f"{'alpha':>6s}", end="")
for name in observers:
    print(f"  {'K_'+name:>14s}", end="")
print(f"  {'K_c/K_m':>10s}  {'K_c/K_p':>10s}")
print("-" * 80)

for alpha_deg in alphas_deg:
    alpha = np.radians(alpha_deg)
    psi = (np.cos(alpha) * tensor(zero, zero) +
           np.sin(alpha) * tensor(one, one)).unit()
    rho0 = ket2dm(psi)

    result = mesolve(H, rho0, times, c_ops, [])

    K_vals = {}
    for name, metric in observers.items():
        tc = find_crossing(result.states, metric)
        if tc is not None:
            K = gamma * tc
            K_vals[name] = K
            results[name].append(K)
        else:
            K_vals[name] = None
            results[name].append(None)

    alpha_list.append(alpha_deg)

    # Print row
    print(f"{alpha_deg:>5.0f}°", end="")
    for name in observers:
        K = K_vals[name]
        if K is not None:
            print(f"  {K:>14.6f}", end="")
        else:
            print(f"  {'never':>14s}", end="")

    # Ratios
    Kc = K_vals.get("Concurrence")
    Km = K_vals.get("Mutual Info")
    Kp = K_vals.get("Pure-Psi")
    if Kc and Km:
        print(f"  {Kc/Km:>10.4f}", end="")
    else:
        print(f"  {'---':>10s}", end="")
    if Kc and Kp:
        print(f"  {Kc/Kp:>10.4f}", end="")
    else:
        print(f"  {'---':>10s}", end="")
    print()

# === Analysis: Ratio stability and geometry ===

print("\n" + "=" * 60)
print("ANALYSIS")
print("=" * 60)

# Compute ratios where both exist
ratios_cm = []
ratios_cp = []
ratios_mp = []
valid_alphas = []

for i, alpha_deg in enumerate(alpha_list):
    Kc = results["Concurrence"][i]
    Km = results["Mutual Info"][i]
    Kp = results["Pure-Psi"][i]

    if Kc and Km and Kp:
        ratios_cm.append(Kc / Km)
        ratios_cp.append(Kc / Kp)
        ratios_mp.append(Km / Kp)
        valid_alphas.append(alpha_deg)

print(f"\nK-ratios across {len(valid_alphas)} states "
      f"(alpha = {min(valid_alphas)}° to {max(valid_alphas)}°):")
for name, vals in [("K_Conc/K_MI", ratios_cm),
                   ("K_Conc/K_Psi", ratios_cp),
                   ("K_MI/K_Psi", ratios_mp)]:
    arr = np.array(vals)
    print(f"  {name:>14s}: mean={arr.mean():.6f} "
          f"std={arr.std():.6f} CV={100*arr.std()/arr.mean():.2f}% "
          f"min={arr.min():.6f} max={arr.max():.6f}")

# Fubini-Study distance from |00> to cos(α)|00>+sin(α)|11>
print(f"\nFubini-Study distance vs K:")
print(f"{'alpha':>6s} {'d_FS':>10s}", end="")
for name in observers:
    print(f"  {'K_'+name:>14s}", end="")
print()

for i, alpha_deg in enumerate(alpha_list):
    if alpha_deg % 10 == 0 or alpha_deg in [30, 31, 35, 45]:
        alpha = np.radians(alpha_deg)
        d_FS = alpha  # Fubini-Study = arccos(|<psi1|psi2>|), here = alpha
        print(f"{alpha_deg:>5.0f}° {d_FS:>10.4f}", end="")
        for name in observers:
            K = results[name][i]
            if K is not None:
                print(f"  {K:>14.6f}", end="")
            else:
                print(f"  {'never':>14s}", end="")
        print()

# Check if K is a function of Fubini-Study distance
print(f"\nIs K monotone in d_FS?")
for name in observers:
    vals = [(np.radians(alpha_list[i]), results[name][i])
            for i in range(len(alpha_list)) if results[name][i] is not None]
    if len(vals) < 2:
        print(f"  {name}: too few crossings")
        continue
    monotone = all(vals[i][1] >= vals[i+1][1] for i in range(len(vals)-1))
    anti_monotone = all(vals[i][1] <= vals[i+1][1] for i in range(len(vals)-1))
    if monotone:
        print(f"  {name}: MONOTONE DECREASING (K decreases with alpha)")
    elif anti_monotone:
        print(f"  {name}: MONOTONE INCREASING (K increases with alpha)")
    else:
        print(f"  {name}: NOT monotone (interesting!)")

print("\nDone.")
