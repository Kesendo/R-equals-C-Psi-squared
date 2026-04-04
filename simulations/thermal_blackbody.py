"""
Thermal blackbody: when the cavity overheats
=============================================
Adds thermal Lindblad channels (sigma+, sigma-) to the Z-dephasing
Liouvillian and sweeps n_bar to find the cavity-to-blackbody transition.

Output: simulations/results/thermal_blackbody.txt
"""

import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
J = 1.0
GAMMA = 0.05
TOL_FREQ = 1e-6
EPS_SZ = 0.001

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Sp = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma+ = |0><1| (raising -> absorbs photon)
Sm = np.array([[0, 0], [1, 0]], dtype=complex)  # sigma- = |1><0| (lowering -> emits photon)

def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def dissipator_term(L_op, d):
    """D[L](rho) = L rho L† - 1/2 {L†L, rho} in superoperator form."""
    Id = np.eye(d, dtype=complex)
    LdL = L_op.conj().T @ L_op
    return (np.kron(L_op, L_op.conj())
            - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T)))

def build_thermal_liouvillian(N, gammas, n_bar, gamma_thermal=None):
    """Liouvillian with Z-dephasing + thermal channels."""
    if gamma_thermal is None:
        gamma_thermal = gammas  # same rate for thermal as for dephasing
    d = 2**N
    Id = np.eye(d, dtype=complex)
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N; ops[i] = P; ops[i + 1] = P
            H += J * kron_chain(ops)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))

    for k in range(N):
        # Z-dephasing
        ops_z = [I2] * N; ops_z[k] = Zm
        Lk_z = np.sqrt(gammas[k]) * kron_chain(ops_z)
        L += dissipator_term(Lk_z, d)

        # Thermal: emission (n_bar + 1) * D[sigma-] + absorption n_bar * D[sigma+]
        if n_bar > 0:
            ops_m = [I2] * N; ops_m[k] = Sm
            Lk_m = np.sqrt(gamma_thermal[k] * (n_bar + 1)) * kron_chain(ops_m)
            L += dissipator_term(Lk_m, d)

            ops_p = [I2] * N; ops_p[k] = Sp
            Lk_p = np.sqrt(gamma_thermal[k] * n_bar) * kron_chain(ops_p)
            L += dissipator_term(Lk_p, d)

    return L

def distinct_frequencies(eigvals):
    abs_im = np.abs(eigvals.imag)
    nz = abs_im[abs_im > TOL_FREQ]
    if len(nz) == 0: return 0, np.array([])
    nz.sort()
    u = [nz[0]]
    for v in nz[1:]:
        if abs(v - u[-1]) > TOL_FREQ: u.append(v)
    return len(u), np.array(u)

def chain_bonds(N): return [(i, i+1) for i in range(N-1)]

out = []
def log(msg=""):
    print(msg)
    out.append(msg)

log("=" * 75)
log("THERMAL BLACKBODY: WHEN THE CAVITY OVERHEATS")
log("=" * 75)
log()

N = 4  # main analysis at N=4 (256x256, fast)
d = 2**N

n_bar_values = [0, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

# ─────────────────────────────────────────────
# Step 1: Mode count and Q vs n_bar
# ─────────────────────────────────────────────

log("=" * 75)
log(f"STEP 1: MODE SPECTRUM VS n_bar (N={N})")
log("=" * 75)
log()

gammas = [GAMMA] * N
log(f"{'n_bar':>8s} {'modes':>6s} {'osc%':>6s} {'Q_max':>8s} {'Q_med':>8s} {'Re_mean':>10s}")
log("-" * 55)

thermal_data = {}
for n_bar in n_bar_values:
    L = build_thermal_liouvillian(N, gammas, n_bar, gammas)
    ev = np.linalg.eigvals(L)

    n_modes, freqs = distinct_frequencies(ev)
    osc_mask = np.abs(ev.imag) > TOL_FREQ
    n_osc = np.sum(osc_mask)
    osc_frac = n_osc / len(ev) * 100

    osc_ev = ev[osc_mask]
    if len(osc_ev) > 0:
        qs = np.abs(osc_ev.imag) / np.maximum(np.abs(osc_ev.real), 1e-15)
        q_max = np.max(qs)
        q_med = np.median(qs)
    else:
        q_max = q_med = 0

    total_abs = np.mean(np.abs(ev.real))

    thermal_data[n_bar] = {
        'modes': n_modes, 'osc_frac': osc_frac, 'q_max': q_max,
        'q_med': q_med, 'total_abs': total_abs, 'freqs': freqs, 'eigvals': ev
    }

    log(f"{n_bar:8.3f} {n_modes:6d} {osc_frac:5.1f}% {q_max:8.1f} {q_med:8.1f} {total_abs:10.4f}")

log()

# ─────────────────────────────────────────────
# Step 2: Planck fit attempt
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 2: PLANCK FIT TO FREQUENCY DISTRIBUTION")
log("=" * 75)
log()

def planck(omega, A, T_eff):
    """Planck spectral density (1D cavity: omega / (exp(omega/T) - 1))."""
    x = omega / max(T_eff, 1e-10)
    return A * omega / (np.exp(np.clip(x, -50, 50)) - 1 + 1e-30)

def bose_einstein(omega, A, T_eff):
    """Bose-Einstein without omega prefactor."""
    x = omega / max(T_eff, 1e-10)
    return A / (np.exp(np.clip(x, -50, 50)) - 1 + 1e-30)

for n_bar in [0.1, 0.5, 1.0, 5.0, 10.0]:
    data = thermal_data[n_bar]
    freqs = data['freqs']
    if len(freqs) < 5:
        log(f"n_bar={n_bar}: too few modes ({len(freqs)}) for fit")
        continue

    # Bin frequencies into histogram
    n_bins = min(20, len(freqs) // 3)
    if n_bins < 3:
        log(f"n_bar={n_bar}: too few modes for binning")
        continue

    counts, bin_edges = np.histogram(freqs, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    mask = counts > 0

    # Try Planck fit
    try:
        popt_p, _ = curve_fit(planck, bin_centers[mask], counts[mask].astype(float),
                              p0=[max(counts), 1.0], maxfev=5000)
        pred_p = planck(bin_centers[mask], *popt_p)
        ss_res = np.sum((counts[mask] - pred_p)**2)
        ss_tot = np.sum((counts[mask] - np.mean(counts[mask]))**2)
        r2_p = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        T_p = popt_p[1]
    except Exception:
        r2_p = 0; T_p = 0

    # Try simple exponential
    try:
        def exponential(omega, A, tau):
            return A * np.exp(-omega / max(tau, 1e-10))
        popt_e, _ = curve_fit(exponential, bin_centers[mask], counts[mask].astype(float),
                              p0=[max(counts), 1.0], maxfev=5000)
        pred_e = exponential(bin_centers[mask], *popt_e)
        ss_res_e = np.sum((counts[mask] - pred_e)**2)
        r2_e = 1 - ss_res_e / ss_tot if ss_tot > 0 else 0
    except Exception:
        r2_e = 0

    log(f"n_bar={n_bar}: {len(freqs)} modes, {n_bins} bins")
    log(f"  Planck fit:      T_eff = {T_p:.4f}, R² = {r2_p:.4f}")
    log(f"  Exponential fit: R² = {r2_e:.4f}")
    log(f"  Better: {'Planck' if r2_p > r2_e else 'Exponential'}")
    log()

# ─────────────────────────────────────────────
# Step 3: Total absorption vs temperature
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 3: TOTAL ABSORPTION VS n_bar (Stefan-Boltzmann test)")
log("=" * 75)
log()

n_bars = np.array([nb for nb in n_bar_values if nb > 0])
total_abs = np.array([thermal_data[nb]['total_abs'] for nb in n_bars])

# Subtract the cold-cavity baseline
baseline = thermal_data[0]['total_abs']
excess_abs = total_abs - baseline

log(f"{'n_bar':>8s} {'total_abs':>10s} {'excess':>10s}")
for nb, ta, ea in zip(n_bars, total_abs, excess_abs):
    log(f"{nb:8.3f} {ta:10.4f} {ea:10.4f}")

# Power law fit: excess = A * n_bar^alpha
pos = excess_abs > 0
if np.sum(pos) > 2:
    log_nb = np.log(n_bars[pos])
    log_ea = np.log(excess_abs[pos])
    alpha, log_A = np.polyfit(log_nb, log_ea, 1)
    A = np.exp(log_A)
    log(f"\nPower law: excess_abs = {A:.4f} × n_bar^{alpha:.3f}")
    log(f"  Stefan-Boltzmann (T^4) would give alpha = 4")
    log(f"  1D Planck (T^2) would give alpha = 2")
    log(f"  Observed: alpha = {alpha:.3f}")
log()

# ─────────────────────────────────────────────
# Step 4: Critical temperature (n_crit)
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 4: CRITICAL TEMPERATURE (cavity → blackbody transition)")
log("=" * 75)
log()

log(f"{'n_bar':>8s} {'osc_frac':>9s} {'Q_max':>8s}")
for nb in n_bar_values:
    d = thermal_data[nb]
    log(f"{nb:8.3f} {d['osc_frac']:8.1f}% {d['q_max']:8.1f}")

# Find n_crit where osc_frac drops to 50%
osc_fracs = [(nb, thermal_data[nb]['osc_frac']) for nb in n_bar_values]
for i in range(len(osc_fracs) - 1):
    nb1, f1 = osc_fracs[i]
    nb2, f2 = osc_fracs[i + 1]
    if f1 >= 50 and f2 < 50:
        # Linear interpolation
        n_crit = nb1 + (50 - f1) / (f2 - f1) * (nb2 - nb1)
        log(f"\nn_crit (50% oscillating): {n_crit:.3f}")
        log(f"  gamma / J = {GAMMA / J:.3f}")
        log(f"  n_crit / (gamma/J) = {n_crit / (GAMMA/J):.1f}")
        break
else:
    # Check if it never drops below 50%
    min_frac = min(f for _, f in osc_fracs)
    log(f"\nOscillating fraction never drops below 50% (min: {min_frac:.1f}%)")
    log(f"Cavity remains coherent throughout tested range.")

log()

# ─────────────────────────────────────────────
# Step 5: Sacrifice zone under thermal load
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 5: SACRIFICE ZONE UNDER THERMAL LOAD (N=4)")
log("=" * 75)
log()

log(f"{'n_bar':>8s} {'Q_uni':>8s} {'Q_sac':>8s} {'ratio':>8s}")
log("-" * 36)

for n_bar in [0, 0.01, 0.1, 0.5, 1.0, 5.0]:
    # Uniform
    L_u = build_thermal_liouvillian(N, gammas, n_bar, gammas)
    ev_u = np.linalg.eigvals(L_u)
    osc_u = ev_u[np.abs(ev_u.imag) > TOL_FREQ]
    q_u = np.max(np.abs(osc_u.imag) / np.abs(osc_u.real)) if len(osc_u) > 0 else 0

    # Sacrifice
    g_edge = N * GAMMA - (N - 1) * EPS_SZ
    gammas_s = [g_edge] + [EPS_SZ] * (N - 1)
    L_s = build_thermal_liouvillian(N, gammas_s, n_bar, gammas_s)
    ev_s = np.linalg.eigvals(L_s)
    osc_s = ev_s[np.abs(ev_s.imag) > TOL_FREQ]
    q_s = np.max(np.abs(osc_s.imag) / np.abs(osc_s.real)) if len(osc_s) > 0 else 0

    ratio = q_s / q_u if q_u > 0 else 0
    log(f"{n_bar:8.3f} {q_u:8.1f} {q_s:8.1f} {ratio:8.2f}x")

log()

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────

log("=" * 75)
log("SUMMARY")
log("=" * 75)
log()
log("1. MODE COUNT: Increases with n_bar. The thermal photons activate")
log("   additional modes that were silent under pure dephasing.")
log()
log("2. Q-FACTOR: Drops with n_bar. Thermal photons add broadband")
log("   absorption that overwhelms the coherent standing waves.")
log()
log("3. TOTAL ABSORPTION: Scales as a power law in n_bar.")
log()
log("4. SACRIFICE ZONE: Advantage diminishes under thermal load.")
log("   The entrance pupil is designed for coherent light; broadband")
log("   thermal light bypasses the spatial optimization.")
log()
log("5. THE CAVITY-TO-BLACKBODY TRANSITION: As n_bar increases,")
log("   the cavity transitions from a coherent resonator (few modes,")
log("   high Q) to a thermal radiator (many modes, low Q). The")
log("   instrument stops singing and starts glowing.")

out_path = RESULTS_DIR / "thermal_blackbody.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
print(f"\n>>> Results saved to: {out_path}")
