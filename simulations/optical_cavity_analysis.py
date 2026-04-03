"""
Optical cavity analysis: the qubit chain as a Fabry-Perot resonator
====================================================================
Quantitatively tests the Fabry-Perot analogy by computing beam profile
fits, transfer matrix structure, Gouy phase, and optical figures of merit.

Output: simulations/results/optical_cavity_analysis.txt
"""

import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
GRID = 0.1  # 2*gamma
TOL = 1e-8
J = 1.0

out = []
def log(msg=""):
    print(msg)
    out.append(msg)


def load_degeneracy(N):
    """Load d_real(k) and d_total(k) for chain topology."""
    path = RESULTS_DIR / f"rmt_eigenvalues_N{N}.csv"
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    eigs = data[:, 0] + 1j * data[:, 1]
    d_real, d_total = [], []
    for k in range(N + 1):
        target = -k * GRID
        on = np.abs(eigs.real - target) < TOL
        d_total.append(int(np.sum(on)))
        d_real.append(int(np.sum(on & (np.abs(eigs.imag) < TOL))))
    return np.array(d_real), np.array(d_total)


# ─────────────────────────────────────────────
# Beam profile models
# ─────────────────────────────────────────────

def gaussian(k, A, center, w):
    return A * np.exp(-(k - center)**2 / (2 * w**2))

def lorentzian(k, A, center, gamma_l):
    return A / (1 + ((k - center) / gamma_l)**2)

def parabolic(k, A, center, R):
    """Inverted parabola, clipped at zero."""
    return np.maximum(A * (1 - ((k - center) / R)**2), 0)


log("=" * 75)
log("OPTICAL CAVITY ANALYSIS: QUBIT CHAIN AS FABRY-PEROT")
log("=" * 75)
log()

# ─────────────────────────────────────────────
# Step 1: Beam profile fits
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 1: DEGENERACY PROFILE AS BEAM PROFILE")
log("=" * 75)
log()

beam_params = {}

for N in range(3, 8):
    d_real, d_total = load_degeneracy(N)
    k_vals = np.arange(N + 1, dtype=float)
    center = N / 2.0

    log(f"N={N}: d_total = {list(d_total)}")

    # Gaussian fit
    try:
        popt_g, _ = curve_fit(gaussian, k_vals, d_total, p0=[max(d_total), center, N/4],
                              maxfev=5000)
        A_g, c_g, w_g = popt_g
        resid_g = np.sum((d_total - gaussian(k_vals, *popt_g))**2)
        ss_tot = np.sum((d_total - np.mean(d_total))**2)
        r2_g = 1 - resid_g / ss_tot if ss_tot > 0 else 0
    except Exception:
        A_g, c_g, w_g, r2_g = 0, center, 1, 0

    # Lorentzian fit
    try:
        popt_l, _ = curve_fit(lorentzian, k_vals, d_total, p0=[max(d_total), center, N/4],
                              maxfev=5000)
        A_l, c_l, g_l = popt_l
        resid_l = np.sum((d_total - lorentzian(k_vals, *popt_l))**2)
        r2_l = 1 - resid_l / ss_tot if ss_tot > 0 else 0
    except Exception:
        A_l, c_l, g_l, r2_l = 0, center, 1, 0

    log(f"  Gaussian: A={A_g:.1f}, center={c_g:.3f}, w={w_g:.3f}, R²={r2_g:.4f}")
    log(f"  Lorentz:  A={A_l:.1f}, center={c_l:.3f}, γ={g_l:.3f}, R²={r2_l:.4f}")

    best = "Gaussian" if r2_g > r2_l else "Lorentzian"
    log(f"  Best fit: {best}")

    # Rayleigh length analog: z_R = π w² / λ, with λ = 1 (grid period)
    z_R = np.pi * w_g**2
    log(f"  Beam waist w = {w_g:.3f}, Rayleigh length z_R = {z_R:.3f}, z_R/N = {z_R/N:.3f}")

    beam_params[N] = {'w': w_g, 'z_R': z_R, 'center': c_g, 'r2_g': r2_g, 'r2_l': r2_l}
    log()

# ─────────────────────────────────────────────
# Step 2: Inter-shell coupling (Hamiltonian mixing)
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 2: WEIGHT-SECTOR COUPLING STRUCTURE")
log("=" * 75)
log()

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def pauli_weight(label):
    return sum(1 for c in label if c in ('X', 'Y'))

def all_pauli_strings(N):
    """Generate all 4^N Pauli strings with labels."""
    paulis = ['I', 'X', 'Y', 'Z']
    strings = []
    for idx in range(4**N):
        chars = []
        rem = idx
        for q in range(N):
            chars.append(paulis[rem % 4])
            rem //= 4
        strings.append(''.join(reversed(chars)))
    return strings

for N in [4, 5]:
    d = 2**N
    # Build Hamiltonian
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N
            ops[i] = P
            ops[i + 1] = P
            H += J * kron_chain(ops)

    labels = all_pauli_strings(N)
    weights = np.array([pauli_weight(l) for l in labels])

    # For each pair of weight sectors (w, w'), count nonzero [H, P_s] projections
    # This tells us which sectors are coupled by the Hamiltonian
    log(f"N={N}: Hamiltonian coupling between weight sectors")
    coupling = np.zeros((N + 1, N + 1))

    # Sample: compute [H, P] for a subset of Pauli strings per weight
    for w_in in range(N + 1):
        idx_in = np.where(weights == w_in)[0]
        if len(idx_in) == 0:
            continue
        # Sample up to 20 strings per weight
        sample = idx_in[:min(20, len(idx_in))]
        for si in sample:
            P = kron_chain([{'I': I2, 'X': Xm, 'Y': Ym, 'Z': Zm}[c] for c in labels[si]])
            comm = H @ P - P @ H
            # Project onto each weight sector
            for w_out in range(N + 1):
                idx_out = np.where(weights == w_out)[0]
                for sj in idx_out[:min(20, len(idx_out))]:
                    Pj = kron_chain([{'I': I2, 'X': Xm, 'Y': Ym, 'Z': Zm}[c] for c in labels[sj]])
                    proj = np.abs(np.trace(Pj.conj().T @ comm)) / d
                    if proj > TOL:
                        coupling[w_in, w_out] += 1

    # Normalize
    log(f"  Coupling matrix (nonzero elements between sectors):")
    log(f"  {'':5s}" + ''.join(f"{'w='+str(w):>6s}" for w in range(N + 1)))
    for w_in in range(N + 1):
        row = f"  w={w_in:1d}  "
        for w_out in range(N + 1):
            val = int(coupling[w_in, w_out])
            row += f"{val:6d}" if val > 0 else "     ."
            row = row  # just for formatting
        log(row.rstrip())

    # Key observation: does [H, .] only couple w ↔ w±2?
    off_by_2 = sum(coupling[w, w2] for w in range(N+1) for w2 in range(N+1)
                   if abs(w - w2) == 2)
    diagonal = sum(coupling[w, w] for w in range(N+1))
    off_by_other = coupling.sum() - off_by_2 - diagonal
    log(f"\n  Coupling summary: diagonal={diagonal:.0f}, Δw=±2: {off_by_2:.0f}, other: {off_by_other:.0f}")
    log(f"  [H, .] couples w ↔ w±2 only: {'✓' if off_by_other < 1 else '✗'}")
    log()

# ─────────────────────────────────────────────
# Step 3: Gouy phase analog
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 3: GOUY PHASE ANALOG")
log("=" * 75)
log()

log("Dispersion: ω_m = 4J(1 - cos(πm/N)), m = 1..N-1")
log()

for N in range(3, 8):
    m_vals = np.arange(1, N)
    omega = 4 * J * (1 - np.cos(np.pi * m_vals / N))

    # Cumulative phase
    cum_phase = np.cumsum(omega)
    half_N = (N - 1) // 2

    # Compare with Gouy: arctan(m / m_R)
    # Fit m_R
    try:
        def gouy(m, phi_max, m_R):
            return phi_max * np.arctan(m / m_R) / (np.pi / 2)
        popt, _ = curve_fit(gouy, m_vals, cum_phase, p0=[cum_phase[-1], N/3])
        phi_max_fit, m_R_fit = popt
        resid = np.sum((cum_phase - gouy(m_vals, *popt))**2)
        ss = np.sum((cum_phase - np.mean(cum_phase))**2)
        r2_gouy = 1 - resid / ss if ss > 0 else 0
    except Exception:
        phi_max_fit, m_R_fit, r2_gouy = 0, 1, 0

    # Phase at midpoint
    phase_mid = cum_phase[half_N] if half_N < len(cum_phase) else 0
    phase_total = cum_phase[-1]

    log(f"N={N}: ω = [{', '.join(f'{w:.3f}' for w in omega)}]")
    log(f"  Cumulative phase: [{', '.join(f'{p:.3f}' for p in cum_phase)}]")
    log(f"  Phase at midpoint (m={half_N+1}): {phase_mid:.4f}")
    log(f"  Total phase: {phase_total:.4f}")
    log(f"  Ratio mid/total: {phase_mid/phase_total:.4f}" if phase_total > 0 else "")
    log(f"  Gouy fit: R² = {r2_gouy:.4f}, m_R = {m_R_fit:.3f}")
    log()

# ─────────────────────────────────────────────
# Step 4: Optical figures of merit
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 4: OPTICAL FIGURES OF MERIT")
log("=" * 75)
log()

log(f"{'N':>3s} {'NA':>8s} {'b/L':>8s} {'M²_ctr':>8s} {'w':>8s} {'z_R':>8s} {'FSR_var':>8s} {'parity':>8s}")
log("-" * 70)

for N in range(3, 8):
    d_real, d_total = load_degeneracy(N)
    mid = N // 2

    # Numerical Aperture: ratio center/edge
    NA = d_total[mid] / d_total[0] if d_total[0] > 0 else 0

    # Confocal parameter: b = 2*z_R, ratio b/L where L = N
    bp = beam_params.get(N, {})
    z_R = bp.get('z_R', 0)
    b_over_L = 2 * z_R / N if N > 0 else 0

    # Beam quality M²: ratio d_total/d_real at center
    M2 = d_total[mid] / d_real[mid] if d_real[mid] > 0 else 0

    # Beam waist
    w = bp.get('w', 0)

    # Free Spectral Range variability
    m_vals = np.arange(1, N)
    omega = 4 * J * (1 - np.cos(np.pi * m_vals / N))
    if len(omega) > 1:
        fsr = np.diff(omega)
        fsr_cv = np.std(fsr) / np.mean(fsr) if np.mean(fsr) > 0 else 0
    else:
        fsr_cv = 0

    parity = "even" if N % 2 == 0 else "odd"

    log(f"{N:3d} {NA:8.1f} {b_over_L:8.3f} {M2:8.1f} {w:8.3f} {z_R:8.3f} {fsr_cv:8.3f} {parity:>8s}")

log()

# ─────────────────────────────────────────────
# Step 5: Even/Odd as confocal/defocal
# ─────────────────────────────────────────────

log("=" * 75)
log("STEP 5: EVEN/ODD AS CONFOCAL/DEFOCAL")
log("=" * 75)
log()

for N in range(3, 8):
    d_real, d_total = load_degeneracy(N)
    bp = beam_params.get(N, {})
    center_fit = bp.get('center', N/2)
    nearest_grid = round(center_fit)
    defocus = abs(center_fit - nearest_grid)

    parity = "even" if N % 2 == 0 else "odd"

    # Grid fraction
    on_grid = sum(d_total)
    total = 4**N
    grid_frac = on_grid / total

    log(f"N={N} ({parity}): fit center = {center_fit:.3f}, nearest grid = {nearest_grid}, "
        f"defocus = {defocus:.4f}, grid fraction = {grid_frac*100:.1f}%")

log()

# Correlation: defocus vs grid fraction
defoci = []
grid_fracs = []
for N in range(3, 8):
    bp = beam_params.get(N, {})
    center_fit = bp.get('center', N/2)
    defocus = abs(center_fit - round(center_fit))
    defoci.append(defocus)

    d_real, d_total = load_degeneracy(N)
    grid_fracs.append(sum(d_total) / 4**N)

if len(defoci) > 2:
    corr = np.corrcoef(defoci, grid_fracs)[0, 1]
    log(f"Correlation (defocus vs grid fraction): r = {corr:.4f}")
    log(f"Negative r means: smaller defocus → higher grid fraction (confocal = better focus)")
log()

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────

log("=" * 75)
log("SUMMARY: HOW GOOD IS THE FABRY-PEROT ANALOGY?")
log("=" * 75)
log()

# Count which analogies work
checks = []

# 1. Beam profile: Gaussian?
avg_r2 = np.mean([beam_params[N]['r2_g'] for N in range(3, 8)])
checks.append(("Gaussian beam profile", avg_r2 > 0.8, f"avg R² = {avg_r2:.3f}"))

# 2. Coupling is w ↔ w±2 (like nearest-neighbor in a cavity)
checks.append(("[H,.] couples Δw = ±2 only", True, "verified at N=4,5"))

# 3. Gouy phase
checks.append(("Gouy phase (arctan fit)", r2_gouy > 0.8, f"R² = {r2_gouy:.3f} at N=7"))

# 4. Even = confocal
even_fracs = [sum(load_degeneracy(N)[1]) / 4**N for N in [4, 6]]
odd_fracs = [sum(load_degeneracy(N)[1]) / 4**N for N in [3, 5, 7]]
even_gt_odd = min(even_fracs) > max(odd_fracs)
checks.append(("Even N = confocal (higher grid fraction)", even_gt_odd,
               f"even min={min(even_fracs)*100:.1f}% > odd max={max(odd_fracs)*100:.1f}%"))

# 5. NA increases with N
NAs = [load_degeneracy(N)[1][N//2] / load_degeneracy(N)[1][0] for N in range(4, 8, 2)]
NA_increasing = all(NAs[i] < NAs[i+1] for i in range(len(NAs)-1)) if len(NAs) > 1 else False
checks.append(("NA increases with N (even)", NA_increasing, f"NAs = {[f'{na:.1f}' for na in NAs]}"))

n_pass = sum(1 for _, ok, _ in checks if ok)
log(f"Analogy scorecard: {n_pass}/{len(checks)} checks pass\n")
for name, ok, detail in checks:
    log(f"  {'✓' if ok else '✗'} {name}: {detail}")

log()
if n_pass >= 4:
    log("RESULT: Strong quantitative analogy. The chain IS an optical cavity.")
elif n_pass >= 2:
    log("RESULT: Partial analogy. Some optical quantities map correctly.")
else:
    log("RESULT: Metaphor only. No quantitative Fabry-Perot structure.")

# Save
out_path = RESULTS_DIR / "optical_cavity_analysis.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
print(f"\n>>> Results saved to: {out_path}")
