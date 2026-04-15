"""
Primordial Gamma Constant: Probe 1, Effective Gamma from 3-Qubit Stack
=======================================================================

Setup: S - M - B chain, XX+YY coupling, Z-dephasing ONLY on B.
gamma_B represents the "framework constant" at the outermost layer.

Question: does the effective gamma seen by the inner observer S
decompose as gamma_eff = gamma_B * f(J_MB / gamma_B)?

Method:
  1. Build full 64x64 Liouvillian
  2. Evolve rho_S(t) = Tr_MB[exp(Lt) rho0]
  3. Extract gamma_eff from initial decay of |rho_S_{01}(t)|
  4. Sweep (J_MB, gamma_B) and test if gamma_eff / gamma_B = f(Q_MB)

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import curve_fit
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def liouvillian(H, jumps):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def partial_trace_MB(rho_8x8):
    """Trace out M and B (qubits 1,2) from a 3-qubit state."""
    r = rho_8x8.reshape(2, 2, 2, 2, 2, 2)
    return np.einsum('ibjcjc->ib', r)


def extract_gamma_eff(L_super, rho0, t_max, n_points=200):
    """Extract effective gamma from coherence decay of rho_S.

    Returns gamma_eff, fit_quality (R^2), has_rebound, rebound_amplitude.
    """
    ts = np.linspace(0, t_max, n_points + 1)
    coh = []

    for t in ts:
        v = rho0.reshape(-1, order='F')
        vt = expm(L_super * t) @ v
        rho_t = vt.reshape((8, 8), order='F')
        rho_S = partial_trace_MB(rho_t)
        coh.append(abs(rho_S[0, 1]))

    coh = np.array(coh)

    # Detect rebound
    has_rebound = False
    rebound_amp = 0.0
    min_coh = coh[0]
    for c in coh[1:]:
        if c < min_coh:
            min_coh = c
        elif c > min_coh + 0.005:
            has_rebound = True
            rebound_amp = max(rebound_amp, c - min_coh)

    # Fit initial decay: |S01|(t) = A * exp(-2*gamma_eff*t)
    # Use early portion (first 30% of time window) for the fit
    n_fit = max(10, n_points // 3)
    ts_fit = ts[:n_fit]
    coh_fit = coh[:n_fit]

    # Avoid log of zero
    valid = coh_fit > 1e-15
    if np.sum(valid) < 3:
        return 0.0, 0.0, has_rebound, rebound_amp

    try:
        def decay_model(t, A, gamma_eff):
            return A * np.exp(-2 * gamma_eff * t)

        popt, _ = curve_fit(decay_model, ts_fit[valid], coh_fit[valid],
                           p0=[coh[0], 0.01], bounds=([0, 0], [1, 100]))
        gamma_eff = popt[1]

        # R^2 for the fit region
        predicted = decay_model(ts_fit[valid], *popt)
        ss_res = np.sum((coh_fit[valid] - predicted) ** 2)
        ss_tot = np.sum((coh_fit[valid] - np.mean(coh_fit[valid])) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    except Exception:
        # Linear fit on log scale as fallback
        log_coh = np.log(coh_fit[valid])
        slope = np.polyfit(ts_fit[valid], log_coh, 1)[0]
        gamma_eff = -slope / 2
        r2 = 0.5  # mark as uncertain

    return gamma_eff, r2, has_rebound, rebound_amp


# =========================================================================
# Setup
# =========================================================================
J_SM = 1.0  # fixed inner coupling

# Sweep parameters
J_MB_values = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
gamma_B_values = [0.01, 0.03, 0.1, 0.3, 1.0]

# Initial state: S in |+>, M and B maximally mixed
rho_plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
rho_mix = 0.5 * I2
rho0 = np.kron(rho_plus, np.kron(rho_mix, rho_mix))

results_dir = Path("simulations/results/primordial_gamma")

print("=" * 80)
print("Probe 1: Effective gamma from 3-qubit S-M-B stack")
print("=" * 80)
print(f"J_SM = {J_SM} (fixed)")
print(f"J_MB sweep: {J_MB_values}")
print(f"gamma_B sweep: {gamma_B_values}")
print()

# =========================================================================
# Main sweep
# =========================================================================
all_results = []

header = f"{'J_MB':>8} {'gamma_B':>8} {'Q_MB':>8} {'gamma_eff':>10} {'g_eff/g_B':>10} {'R2':>8} {'rebound':>8} {'reb_amp':>8}"
print(header)
print("-" * len(header))

for J_MB in J_MB_values:
    for gamma_B in gamma_B_values:
        Q_MB = J_MB / gamma_B

        H = (J_SM * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2))
           + J_MB * 0.5 * (kron(I2, X, X) + kron(I2, Y, Y)))

        L_jump = np.sqrt(gamma_B) * kron(I2, I2, Z)
        L = liouvillian(H, [L_jump])

        t_max = min(200, max(20, 30 / gamma_B))
        gamma_eff, r2, has_reb, reb_amp = extract_gamma_eff(L, rho0, t_max)

        ratio = gamma_eff / gamma_B if gamma_B > 0 else 0

        result = {
            'J_MB': J_MB, 'gamma_B': gamma_B, 'Q_MB': Q_MB,
            'gamma_eff': gamma_eff, 'ratio': ratio,
            'r2': r2, 'rebound': has_reb, 'reb_amp': reb_amp
        }
        all_results.append(result)

        reb_str = f"{reb_amp:.4f}" if has_reb else "no"
        print(f"{J_MB:8.3f} {gamma_B:8.3f} {Q_MB:8.2f} {gamma_eff:10.6f} {ratio:10.6f} {r2:8.4f} {reb_str:>8}")

# =========================================================================
# Key test: does gamma_eff / gamma_B depend on Q_MB only?
# =========================================================================
print("\n" + "=" * 80)
print("Test: gamma_eff / gamma_B = f(Q_MB)?")
print("=" * 80)

# Group by Q_MB and check if ratio is constant within each group
from collections import defaultdict
q_groups = defaultdict(list)
for r in all_results:
    q_key = round(r['Q_MB'], 4)
    q_groups[q_key].append(r['ratio'])

print(f"\n{'Q_MB':>10} {'ratios':>50} {'spread':>10}")
print("-" * 75)

q_invariant = True
for q in sorted(q_groups.keys()):
    ratios = q_groups[q]
    spread = max(ratios) - min(ratios) if len(ratios) > 1 else 0
    ratios_str = ', '.join(f'{r:.6f}' for r in ratios)
    inv = "OK" if spread < 0.01 else "VARIES"
    if spread >= 0.01:
        q_invariant = False
    print(f"{q:10.4f} {ratios_str:>50} {spread:10.6f}  {inv}")

# =========================================================================
# Extract the function f(Q_MB) = gamma_eff / gamma_B
# =========================================================================
print("\n" + "=" * 80)
print("The function f(Q_MB) = gamma_eff / gamma_B")
print("=" * 80)

# Use the gamma_B = 0.1 column as reference (moderate regime)
ref_results = [r for r in all_results if abs(r['gamma_B'] - 0.1) < 0.01]
print(f"\nReference column: gamma_B = 0.1")
print(f"{'Q_MB':>8} {'f(Q_MB)':>10}")
for r in ref_results:
    print(f"{r['Q_MB']:8.2f} {r['ratio']:10.6f}")

# Try to identify closed form
# Common candidates: f(x) = x/(1+x), f(x) = x^2/(1+x^2), f(x) = 1-exp(-x), etc.
q_vals = np.array([r['Q_MB'] for r in ref_results])
f_vals = np.array([r['ratio'] for r in ref_results])

candidates = {
    'x/(1+x)': lambda x: x / (1 + x),
    'x^2/(1+x^2)': lambda x: x**2 / (1 + x**2),
    'x^2/(4+x^2)': lambda x: x**2 / (4 + x**2),
    '1-exp(-x)': lambda x: 1 - np.exp(-x),
    'x^2/(1+x)^2': lambda x: x**2 / (1 + x)**2,
    'tanh(x)/2': lambda x: np.tanh(x) / 2,
}

print(f"\nClosed-form candidates (max |residual|):")
for name, func in candidates.items():
    pred = func(q_vals)
    max_res = np.max(np.abs(f_vals - pred))
    print(f"  {name:>20}: max_res = {max_res:.6f}")

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 80)
print("VERDICT: Probe 1")
print("=" * 80)

verdict = []
# Check fit quality
r2_values = [r['r2'] for r in all_results]
min_r2 = min(r2_values)
n_rebound = sum(1 for r in all_results if r['rebound'])

if min_r2 < 0.7:
    verdict.append(f"MARKOVIAN FIT FAILS: min R^2 = {min_r2:.4f}")
    verdict.append("Effective gamma is not a sensible concept for this system.")
    verdict.append("STOP: hypothesis loses operational support at simplest test.")
else:
    verdict.append(f"Markovian fit quality: min R^2 = {min_r2:.4f} (acceptable)")
    verdict.append(f"Non-Markovian rebound detected in {n_rebound}/{len(all_results)} configs")

    if q_invariant:
        verdict.append(f"gamma_eff / gamma_B = f(Q_MB): YES, ratio is Q_MB-invariant")
        verdict.append("Refractive-index reading has operational support.")
    else:
        verdict.append(f"gamma_eff / gamma_B = f(Q_MB): NO, ratio varies at fixed Q_MB")
        verdict.append("gamma_eff depends on gamma_B and J_MB independently.")

for line in verdict:
    print(line)

# Save
with open(results_dir / 'probe1_results.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 1: Effective gamma from 3-qubit S-M-B stack\n")
    f.write("=" * 80 + "\n\n")
    f.write(header + "\n")
    for r in all_results:
        reb_str = f"{r['reb_amp']:.4f}" if r['rebound'] else "no"
        f.write(f"{r['J_MB']:8.3f} {r['gamma_B']:8.3f} {r['Q_MB']:8.2f} "
                f"{r['gamma_eff']:10.6f} {r['ratio']:10.6f} {r['r2']:8.4f} {reb_str:>8}\n")
    f.write("\nf(Q_MB) reference (gamma_B=0.1):\n")
    for r in ref_results:
        f.write(f"  Q_MB={r['Q_MB']:.2f}: f={r['ratio']:.6f}\n")
    f.write("\nVerdict:\n")
    for line in verdict:
        f.write(line + "\n")

print(f"\nResults saved to {results_dir / 'probe1_results.txt'}")
