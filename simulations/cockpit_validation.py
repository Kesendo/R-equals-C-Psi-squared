"""
Cockpit Validation on Tomography + Shadow Data
================================================
Tests whether the cockpit instruments give consistent, useful
information on REAL hardware data.

Data sources:
  A. Q52 tomography and a separate simulator fixture with different T1,T2
  B. Shadow Q80 + Q102 (March 9); only Q80 has a simulator match
  C. Scope verdict for the readings this producer actually computes

April 2, 2026
"""
import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit
import json, os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

default_base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base = os.environ.get('RCPSI_REPO_ROOT', default_base)
results_path = os.environ.get(
    'RCPSI_COCKPIT_VALIDATION_OUTPUT',
    os.path.join(base, 'simulations', 'results', 'cockpit_validation.txt'),
)
os.makedirs(os.path.dirname(os.path.abspath(results_path)), exist_ok=True)
_lines = []

def out(s=""):
    print(s)
    _lines.append(s)

def load_json(relpath):
    with open(os.path.join(base, relpath), 'r', encoding='utf-8') as f:
        return json.load(f)


# ================================================================
# HELPERS
# ================================================================
def cpsi_from_rho(rho):
    """CPsi = Purity x Psi-norm for any d."""
    d = rho.shape[0]
    pur = float(np.real(np.trace(rho @ rho)))
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi = l1 / (d - 1) if d > 1 else 0.0
    return pur * psi, pur, psi

def theta_deg(cpsi_val):
    if cpsi_val < 0.25:
        return None
    return float(np.degrees(np.arctan(np.sqrt(4*cpsi_val - 1))))

def bures_distance(rho, sigma):
    sqrt_rho = linalg.sqrtm(rho)
    prod = sqrt_rho @ sigma @ sqrt_rho
    ev = np.real(np.linalg.eigvalsh(prod))
    fid = float(np.sum(np.sqrt(np.maximum(ev, 0))))**2
    return float(np.sqrt(max(0, 2*(1 - np.sqrt(max(0, min(1, fid)))))))

def reconstruct_rho_1q(rho01_re, rho01_im):
    """Reconstruct 1-qubit density matrix assuming populations = 1/2.
    Valid for dephasing from |+> initial state (Z-noise preserves populations)."""
    rho = np.array([[0.5, rho01_re + 1j*rho01_im],
                     [rho01_re - 1j*rho01_im, 0.5]], dtype=complex)
    return rho

def exp_decay(t, a, rate, c):
    return a * np.exp(-rate * t) + c


def q52_crossing_comparison_lines(hardware, simulator):
    """Format stored Q52 metadata and the separate simulator-fixture comparison."""
    hardware_crossing = hardware.get('crossing_us')
    simulator_t1 = float(simulator['T1_us'])
    simulator_t2 = float(simulator['T2_us'])
    simulator_prediction = float(
        simulator['analytical_prediction_generalized']['t_star_us']
    )
    stored_line = (
        f"Stored Q52 crossing metadata: {float(hardware_crossing):.1f} us "
        "(not recomputed here)"
        if hardware_crossing is not None
        else "Stored Q52 crossing metadata: not found (not recomputed here)"
    )
    simulator_line = (
        f"Separate simulator fixture uses T1={simulator_t1:.0f} us, "
        f"T2_parameter={simulator_t2:.0f} us and predicts "
        f"{simulator_prediction:.1f} us"
    )
    if hardware_crossing is None:
        difference_value = "unavailable"
    else:
        cross_fixture_difference = (
            abs(float(hardware_crossing) - simulator_prediction)
            / simulator_prediction
            * 100.0
        )
        difference_value = f"{cross_fixture_difference:.1f}%"
    difference_line = (
        f"Cross-fixture difference: {difference_value} "
        "(not the Q52 same-record generalized comparison)"
    )
    return stored_line, simulator_line, difference_line


# ================================================================
# PART A: Q52 TOMOGRAPHY AND A SEPARATE SIMULATOR FIXTURE
# ================================================================
out("=" * 70)
out("PART A: Q52 TOMOGRAPHY AND A SEPARATE SIMULATOR FIXTURE")
out("=" * 70)

hw = load_json('data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json')
sim = load_json('data/ibm_tomography_feb2026/simulator_test_20260209_125106.json')

out(f"\n  Hardware: Q52, T1={hw['T1_us']:.0f} us, T2_echo={hw['T2_us']:.0f} us, {len(hw['raw_tomography'])} points")
out(f"  Simulator fixture: T1={sim['T1_us']:.0f} us, T2_parameter={sim['T2_us']:.0f} us, {len(sim['analysis'])} points")

# Process hardware
hw_data = []
rho_prev = None
for pt in hw['raw_tomography']:
    rho = np.array(pt['density_matrix_real']) + 1j * np.array(pt['density_matrix_imag'])
    cp, pur, psi = cpsi_from_rho(rho)
    th = theta_deg(cp)
    dB = bures_distance(rho, rho_prev) if rho_prev is not None else None
    rho_prev = rho.copy()
    hw_data.append({'t': pt['delay_us'], 'cpsi': cp, 'theta': th,
                    'pur': pur, 'psi': psi, 'fid': pt['fidelity'],
                    'bures': dB, 'rho01': abs(rho[0,1])})

# Process simulator
sim_data = []
for pt in sim['analysis']:
    sim_data.append({'t': pt['delay_us'], 'cpsi': pt['cpsi_measured'],
                     'pur': pt['C_measured'], 'psi': pt['psi_measured'],
                     'rho01': pt['populations']['rho_01_abs'],
                     'cpsi_theory': pt['cpsi_theory_full']})

# Cross-fixture table at nearest sampled times
out(f"\n  CROSS-FIXTURE TABLE: Q52 hardware vs separate simulator fixture")
out("  Different T1,T2; these row differences are not a same-record agreement test.")
out(
    "  v_B=d_B/dt: backward finite-step Bures speed over "
    "[t_(i-1),t_i] (1/us); N/A at first sample"
)
out(f"\n  {'t(us)':>8} | {'CPsi_q52':>9} {'CPsi_fix':>9} {'diff':>7} | {'theta_q52(deg)':>14} {'Pur_q52':>8} {'Psi_q52':>8} {'v_B(1/us)':>11}")
out(f"  {'-'*89}")

for hpt in hw_data:
    # Find closest simulator point
    best_sim = min(sim_data, key=lambda s: abs(s['t'] - hpt['t']))
    dt_match = abs(best_sim['t'] - hpt['t'])
    if dt_match < 30:  # within 30 us
        diff = hpt['cpsi'] - best_sim['cpsi']
        dt_sample = hpt['t'] - (hw_data[hw_data.index(hpt)-1]['t'] if hw_data.index(hpt) > 0 else 0)
        vB = (
            hpt['bures'] / dt_sample
            if hpt['bures'] is not None and dt_sample > 0
            else None
        )
        theta_text = "N/A" if hpt['theta'] is None else f"{hpt['theta']:.1f}"
        vB_text = "N/A" if vB is None else f"{vB:.4f}"
        out(f"  {hpt['t']:>8.1f} | {hpt['cpsi']:>8.4f} {best_sim['cpsi']:>8.4f} {diff:>+7.4f} | "
            f"{theta_text:>14} {hpt['pur']:>7.3f} {hpt['psi']:>7.3f} {vB_text:>11}")

# Crossing metadata and cross-fixture comparison
stored_crossing_line, simulator_crossing_line, cross_fixture_line = (
    q52_crossing_comparison_lines(hw, sim)
)
out("\n  Q52 CROSSING METADATA")
out(f"    {stored_crossing_line}")
out("    This producer does not recompute a same-record generalized prediction.")
out("    No fitted Q52 prediction tuple is repeated or adjudicated here.")
out("\n  Separate simulator fixture / cross-fixture comparison (not a Q52 prediction test)")
out(f"    {simulator_crossing_line}")
out(f"    {cross_fixture_line}")

# Cross-fixture residual summary
hw_cpsi = np.array([d['cpsi'] for d in hw_data])
hw_times = np.array([d['t'] for d in hw_data])
# Interpolate simulator at hardware times
sim_times = np.array([d['t'] for d in sim_data])
sim_cpsi = np.array([d['cpsi'] for d in sim_data])
# Only compare where both have data
common_mask = hw_times <= sim_times[-1]
if np.sum(common_mask) > 3:
    sim_interp = np.interp(hw_times[common_mask], sim_times, sim_cpsi)
    residuals = hw_cpsi[common_mask] - sim_interp
    rmse = np.sqrt(np.mean(residuals**2))
    mean_cpsi = np.mean(hw_cpsi[common_mask])
    out(f"\n  Cross-fixture residual summary (different T1,T2; not an agreement metric):")
    out(f"    RMSE = {rmse:.4f}")
    out(f"    Mean |residual| = {np.mean(np.abs(residuals)):.4f}")
    out(f"    Max |residual| = {np.max(np.abs(residuals)):.4f}")
    out(f"    Relative RMSE = {rmse/mean_cpsi*100:.1f}%")


# ================================================================
# PART B: SHADOW DATA -- Q80 + Q102
# ================================================================
out(f"\n{'=' * 70}")
out("PART B: SHADOW DATA -- Q80 + Q102 (March 9, 2026)")
out("=" * 70)

shadow_hw = load_json('data/ibm_shadow_march2026/shadow_hardware_combined_20260309_181852.json')
shadow_sim = load_json('data/ibm_shadow_march2026/shadow_simulate_20260309_181709.json')

out(f"\n  Hardware: {shadow_hw['backend']}, {shadow_hw['shots']} shots")
out(f"  Delay multiples of T2*: {shadow_hw['delay_multiples']}")

for qr in shadow_hw['qubit_results']:
    qubit_id = qr['qubit']
    points = qr['points']

    out(f"\n  --- Q{qubit_id} HARDWARE ---")
    out(f"  {'t(us)':>8} {'t/T2*':>6} | {'CPsi':>7} {'theta(deg)':>10} {'Purity':>7} {'Psi':>6} | "
        f"{'|r01|':>7} {'ph(deg)':>7} | {'v_B(1/us)':>11}")
    out(f"  {'-'*84}")

    rho_prev = None
    cpsi_arr, t_arr, bures_arr = [], [], []

    for pt in points:
        t_us = pt['delay_us']
        cpsi = pt['cpsi']
        th = theta_deg(cpsi)
        purity = pt['C']
        psi = pt['psi']
        r01_abs = pt['rho01_abs']
        r01_re = pt['rho01_re']
        r01_im = pt['rho01_im']
        phase = pt['rho01_phase_deg']

        # Reconstruct density matrix for Bures
        rho = reconstruct_rho_1q(r01_re, r01_im)
        if rho_prev is not None:
            dB = bures_distance(rho, rho_prev)
            dt_us = t_us - t_prev
            vB = dB / dt_us if dt_us > 0 else 0
        else:
            dB, vB = None, None
        rho_prev = rho.copy()
        t_prev = t_us

        t_arr.append(t_us)
        cpsi_arr.append(cpsi)
        bures_arr.append(vB)

        theta_text = "N/A" if th is None else f"{th:.1f}"
        vB_text = "N/A" if vB is None else f"{vB:.4f}"
        out(f"  {t_us:>8.2f} {pt['t_over_T2star']:>6.2f} | {cpsi:>7.4f} {theta_text:>10} "
            f"{purity:>7.3f} {psi:>6.3f} | {r01_abs:>7.4f} {phase:>+7.1f} | {vB_text:>11}")

    cpsi_arr = np.array(cpsi_arr)
    t_arr = np.array(t_arr)

    # Crossing
    t_cross = None
    crossing_bracket = None
    for i in range(1, len(cpsi_arr)):
        if cpsi_arr[i-1] >= 0.25 and cpsi_arr[i] < 0.25:
            frac = (0.25 - cpsi_arr[i]) / (cpsi_arr[i-1] - cpsi_arr[i] + 1e-30)
            t_cross = t_arr[i] * (1 - frac) + t_arr[i-1] * frac
            crossing_bracket = (t_arr[i-1], t_arr[i])
            break

    if t_cross:
        out(
            f"  CΨ=1/4 linear-interpolated estimate: {t_cross:.1f} us "
            f"from sampled bracket [{crossing_bracket[0]:.2f}, "
            f"{crossing_bracket[1]:.2f}] us"
        )
    else:
        if cpsi_arr[0] > 0.25:
            out(f"  CΨ starts above 1/4 but never crosses (last CΨ={cpsi_arr[-1]:.4f})")
        else:
            out(f"  CΨ never reaches 1/4")

    # Decay fit
    valid = cpsi_arr > 0.01
    if np.sum(valid) > 3:
        try:
            popt, _ = curve_fit(exp_decay, t_arr[valid], cpsi_arr[valid],
                                p0=[cpsi_arr[0], 0.05, 0.01], maxfev=5000)
            out(f"  Finite Q{qubit_id} fit CΨ=a exp(-r t)+c: "
                f"rate r = {popt[1]:.5f} 1/us, fitted floor c = {popt[2]:.5f}")
        except Exception:
            out(f"  CΨ decay fit failed")

# --- Shadow simulator availability ---
out(f"\n  SHADOW SIMULATOR MATCH CHECK")

for qr_hw in shadow_hw['qubit_results']:
    qid = qr_hw['qubit']
    # Find matching simulator qubit
    qr_sim = None
    for qr in shadow_sim['qubit_results']:
        if qr['qubit'] == qid:
            qr_sim = qr
            break

    if qr_sim is None:
        out(f"\n  Q{qid}: No simulator match found")
        continue

    out(f"\n  Q{qid}: Hardware vs Simulation")
    out(f"  {'t/T2*':>6} | {'CPsi_hw':>8} {'CPsi_sim':>8} {'diff':>7} | "
        f"{'|r01|_hw':>8} {'|r01|_sim':>8} {'diff':>7}")
    out(f"  {'-'*65}")

    cpsi_diffs, r01_diffs = [], []
    for pt_hw, pt_sim in zip(qr_hw['points'], qr_sim['points']):
        d_cpsi = pt_hw['cpsi'] - pt_sim['cpsi']
        d_r01 = pt_hw['rho01_abs'] - pt_sim['rho01_abs']
        cpsi_diffs.append(abs(d_cpsi))
        r01_diffs.append(abs(d_r01))
        out(f"  {pt_hw['t_over_T2star']:>6.2f} | {pt_hw['cpsi']:>8.4f} "
            f"{pt_sim['cpsi']:>8.4f} {d_cpsi:>+7.4f} | {pt_hw['rho01_abs']:>8.4f} "
            f"{pt_sim['rho01_abs']:>8.4f} {d_r01:>+7.4f}")

    out(f"  Mean |CPsi diff| = {np.mean(cpsi_diffs):.4f}")
    out(f"  Mean |rho01 diff| = {np.mean(r01_diffs):.4f}")


# ================================================================
# SCOPE OF THE AVAILABLE READINGS
# ================================================================
out(f"\n{'=' * 70}")
out("SCOPE OF THE AVAILABLE READINGS")
out("=" * 70)

out(f"""
  theta(deg) is an algebraic rereading only for CΨ>=1/4 and is N/A below; Psi is an algebraic rereading of the loaded single-qubit data.
  Q52 v_B values are backward finite-step speeds over [t_(i-1),t_i] in 1/us; the first sample is N/A.
  Shadow v_B values use reconstructed populations=1/2 and the same backward interval convention; each first sample is N/A.
  Finite exponential fits are reported separately for Q80 and Q102; no cross-qubit rate consistency is inferred.
  Only Q80 has a matched shadow simulator record; Q102 explicitly does not.
  This producer computes no 5Q, MI, concurrence, curvature, or Petermann value.
  External Petermann interpretation boundary (not evidence from this run):
  The old K_P ~ 1 pure-dephasing null is refuted.
  A single-vector K_P is meaningful only for a simple isolated mode; at degeneracy use invariant-subspace or Jordan diagnostics.
""")


# ================================================================
# FINITE Q52 LATE-TIME RECORD
# ================================================================
out(f"{'=' * 70}")
out("FINITE Q52 LATE-TIME RECORD")
out("=" * 70)

hw_r01_arr = np.array([d['rho01'] for d in hw_data])
late_mask = hw_times > 300
if np.sum(late_mask) > 3:
    late_r01 = hw_r01_arr[late_mask]
    out(f"\n  Late-time coherence samples (t > 300 us):")
    out(f"    Measured mean |rho01| = {np.mean(late_r01):.5f}")
    out(
        "    These finite samples establish neither a nonzero asymptote nor "
        "a Q52 mechanism."
    )


# ================================================================
# VERDICT
# ================================================================
out(f"\n{'=' * 70}")
out("VERDICT: WHAT DOES THIS DATASET TEST?")
out("=" * 70)

out(f"""
  ANSWER: INCOMPLETELY TESTED

  - Q52 supplies a qualitative crossing record, but it is not a precision prediction match.
  - This producer does not recompute a same-record Q52 prediction.
  - The separate simulator fixture and cross-fixture residuals are not a Q52 prediction test.
  - Q80 is the only shadow record with a matched simulator; Q102 has none here.

  THE CURRENT DATASET LEAVES THE COCKPIT INCOMPLETELY TESTED.
""")

out("=" * 70)
out("ANALYSIS COMPLETE")
out("=" * 70)

# Save
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_lines) + '\n')
print(f"\nResults saved to {results_path}")
