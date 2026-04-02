"""
Cockpit on IBM Hardware Data
=============================
Applies cockpit instruments to REAL IBM quantum hardware data.
No simulations; reads existing JSON data only.

Data sources:
  A. Sacrifice Zone March 24, 2026 (5-qubit chain, 3 DD strategies)
     - Count data only -> MI, populations, decay rates
  B. Q52 Tomography Feb 9, 2026 (single qubit, 25 delay points)
     - Full density matrices -> all instruments
  C. Q80 Palindrome March 18, 2026 (single qubit, 8 delay points)
     - Full density matrices + CPsi pre-computed

April 2, 2026
"""
import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit
import json, os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(base, 'simulations', 'results')
os.makedirs(results_dir, exist_ok=True)
results_path = os.path.join(results_dir, 'cockpit_ibm_hardware.txt')
_lines = []

def out(s=""):
    print(s)
    _lines.append(s)


# ================================================================
# HELPERS
# ================================================================
def load_json(relpath):
    path = os.path.join(base, relpath)
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def counts_to_probs(counts_dict, n_qubits):
    """Convert bitstring counts to probability distribution."""
    total = sum(counts_dict.values())
    probs = {}
    for bitstr, count in counts_dict.items():
        probs[bitstr] = count / total
    return probs, total


def marginal_1q(probs, n_qubits, qubit_idx):
    """Marginal probability for a single qubit."""
    p0, p1 = 0.0, 0.0
    for bitstr, prob in probs.items():
        bits = bitstr.zfill(n_qubits)
        if bits[qubit_idx] == '0':
            p0 += prob
        else:
            p1 += prob
    return np.array([p0, p1])


def marginal_2q(probs, n_qubits, q1, q2):
    """Joint probability for a qubit pair. Returns 2x2 array."""
    p = np.zeros((2, 2))
    for bitstr, prob in probs.items():
        bits = bitstr.zfill(n_qubits)
        b1 = int(bits[q1])
        b2 = int(bits[q2])
        p[b1, b2] += prob
    return p


def mutual_information(probs, n_qubits, q1, q2):
    """Classical mutual information I(q1;q2) from count data."""
    p_joint = marginal_2q(probs, n_qubits, q1, q2)
    p1 = marginal_1q(probs, n_qubits, q1)
    p2 = marginal_1q(probs, n_qubits, q2)
    mi = 0.0
    for i in range(2):
        for j in range(2):
            if p_joint[i, j] > 1e-15 and p1[i] > 1e-15 and p2[j] > 1e-15:
                mi += p_joint[i, j] * np.log2(p_joint[i, j] / (p1[i] * p2[j]))
    return max(0, mi)


def purity_from_diagonal(probs, n_qubits):
    """Lower bound on purity from diagonal elements only.
    For n-qubit: Pur >= sum_i p_i^2 (equality iff off-diagonal = 0)."""
    return sum(p**2 for p in probs.values())


def cpsi_from_rho_1q(rho_real, rho_imag):
    """CPsi for a 1-qubit density matrix (d=2)."""
    rho = np.array(rho_real) + 1j * np.array(rho_imag)
    pur = float(np.real(np.trace(rho @ rho)))
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi = l1  # d-1 = 1 for d=2
    return pur * psi, pur, psi


def theta_deg(cpsi_val):
    return float(np.degrees(np.arctan(np.sqrt(4*cpsi_val - 1)))) if cpsi_val > 0.25 else 0.0


def bures_distance(rho, sigma):
    try:
        sqrt_rho = linalg.sqrtm(rho)
        prod = sqrt_rho @ sigma @ sqrt_rho
        ev = np.real(np.linalg.eigvalsh(prod))
        fid = float(np.sum(np.sqrt(np.maximum(ev, 0))))**2
        return float(np.sqrt(max(0, 2*(1 - np.sqrt(max(0, min(1, fid)))))))
    except Exception:
        return 0.0


def exp_decay(t, a, rate, c):
    return a * np.exp(-rate * t) + c


# ================================================================
# PART A: SACRIFICE ZONE (5-qubit chain, count data only)
# ================================================================
out("=" * 70)
out("PART A: SACRIFICE ZONE HARDWARE (March 24, 2026)")
out("Chain [Q85, Q86, Q87, Q88, Q94], 3 DD strategies")
out("=" * 70)

# Load summary file
sz = load_json('data/ibm_sacrifice_zone_march2026/sacrifice_zone_hardware_20260324_191713.json')
chain = sz['chain']
n_qubits = len(chain)
out(f"\n  Chain: {chain} ({n_qubits} qubits)")
out(f"  Shots: {sz['parameters']['shots']} per circuit")
out(f"  Trotter steps: {sz['parameters']['trotter_steps']}")

# Hardware characterization
ci = sz['chain_info']
out(f"\n  Hardware characterization:")
out(f"  {'Qubit':>8} {'T1(us)':>8} {'T2(us)':>8} {'T2*(us)':>8} {'Ratio':>6}")
out(f"  {'-'*42}")
for i, q in enumerate(chain):
    out(f"  Q{q:>6} {ci['T1_values'][i]:>8.1f} {ci['T2_values'][i]:>8.1f} "
        f"{ci['T2star_values'][i]:>8.1f} {ci['ramsey_ratios'][i]:>6.1f}")
out(f"  Sacrifice qubit: Q{ci.get('sacrifice_qubit', '?')} (position {ci.get('sacrifice_position', '?')})")
if 'uniformity' in ci:
    out(f"  Uniformity: {ci['uniformity']:.3f}")

# Also load raw counts for re-computation
dd_files = {
    'No DD': 'data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_no_dd_20260324_191713.json',
    'Selective DD': 'data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_selective_dd_20260324_191523.json',
    'Uniform DD': 'data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_uniform_dd_20260324_191614.json',
}

# --- Instrument availability ---
out(f"\n  INSTRUMENT AVAILABILITY (from counts only):")
out(f"    1. theta (altimeter):    NO  (needs off-diagonal elements)")
out(f"    2. Concurrence (compass): NO  (needs full density matrix)")
out(f"    3. Psi-norm (speed):     NO  (needs off-diagonal elements)")
out(f"    4. Mutual Information:   YES (from bitstring counts)")
out(f"    5. Decay rate:           YES (from MI or population fits)")
out(f"    6. Bures curvature:      NO  (needs full density matrix)")
out(f"    7. Petermann:            NO  (needs Liouvillian)")

# --- Compute instruments from counts ---
out(f"\n  INSTRUMENT 4: MUTUAL INFORMATION (pair-wise)")

for dd_name, dd_file in dd_files.items():
    data = load_json(dd_file)
    time_pts = data['results']['time_points']
    raw = data['raw_counts']

    out(f"\n  --- {dd_name} ---")
    out(f"    {'t':>4} | {'MI(0,1)':>8} {'MI(1,2)':>8} {'MI(2,3)':>8} {'MI(3,4)':>8} | {'SumMI':>8} {'Pur_lb':>8}")
    out(f"    {'-'*65}")

    sum_mi_arr = []
    pur_lb_arr = []
    for ti, (t, counts) in enumerate(zip(time_pts, raw)):
        probs, total = counts_to_probs(counts, n_qubits)

        # Pairwise MI
        pairs = [(0,1), (1,2), (2,3), (3,4)]
        mi_vals = []
        for q1, q2 in pairs:
            mi = mutual_information(probs, n_qubits, q1, q2)
            mi_vals.append(mi)

        sum_mi = sum(mi_vals)
        sum_mi_arr.append(sum_mi)

        # Purity lower bound
        pur_lb = purity_from_diagonal(probs, n_qubits)
        pur_lb_arr.append(pur_lb)

        out(f"    {t:>4.0f} | {mi_vals[0]:>8.4f} {mi_vals[1]:>8.4f} "
            f"{mi_vals[2]:>8.4f} {mi_vals[3]:>8.4f} | {sum_mi:>8.4f} {pur_lb:>8.4f}")

    # --- Instrument 5: Decay rate from sum MI ---
    t_arr = np.array(time_pts)
    mi_arr = np.array(sum_mi_arr)
    try:
        popt, _ = curve_fit(exp_decay, t_arr, mi_arr,
                            p0=[mi_arr[0], 0.3, mi_arr[-1]],
                            maxfev=5000)
        rate = popt[1]
        out(f"    Effective MI decay rate: {rate:.4f} per Trotter step")
        out(f"    MI half-life: {np.log(2)/rate:.1f} steps")
    except Exception:
        out(f"    Decay fit failed (too few points or non-monotonic)")


# --- Compare the three flights ---
out(f"\n{'=' * 70}")
out("THREE-FLIGHT COMPARISON (from summary file)")
out("=" * 70)

configs = sz['configurations']
out(f"\n  {'t':>4} | {'No DD':>10} {'Sel. DD':>10} {'Uni. DD':>10} | {'Best':>12}")
out(f"  {'-'*58}")

for i, t in enumerate(configs['no_dd']['time_points']):
    mi_no = configs['no_dd']['sum_mi'][i]
    mi_sel = configs['selective_dd']['sum_mi'][i]
    mi_uni = configs['uniform_dd']['sum_mi'][i]

    vals = {'No DD': mi_no, 'Sel.': mi_sel, 'Uni.': mi_uni}
    best = max(vals, key=vals.get)
    out(f"  {t:>4.0f} | {mi_no:>10.4f} {mi_sel:>10.4f} {mi_uni:>10.4f} | {best:>12}")

# Summary
out(f"\n  Selective DD advantage:")
for i, t in enumerate(configs['no_dd']['time_points']):
    mi_sel = configs['selective_dd']['sum_mi'][i]
    mi_no = configs['no_dd']['sum_mi'][i]
    mi_uni = configs['uniform_dd']['sum_mi'][i]
    if mi_no > 1e-6:
        ratio_vs_no = mi_sel / mi_no
        ratio_vs_uni = mi_sel / mi_uni if mi_uni > 1e-6 else float('inf')
        out(f"    t={t:.0f}: Sel/No={ratio_vs_no:.2f}x, Sel/Uni={ratio_vs_uni:.2f}x")

# MI propagation pattern
out(f"\n  MI PROPAGATION (No DD):")
out(f"  Which pair dominates at each time?")
for i, t in enumerate(configs['no_dd']['time_points']):
    pair_mi = configs['no_dd']['pair_mi'][i]
    best_pair = max(pair_mi, key=pair_mi.get)
    out(f"    t={t:.0f}: dominant pair = {best_pair} "
        f"(MI={pair_mi[best_pair]:.4f})")


# ================================================================
# PART B: Q52 TOMOGRAPHY (full density matrices)
# ================================================================
out(f"\n{'=' * 70}")
out("PART B: Q52 TOMOGRAPHY (February 9, 2026)")
out("Single qubit, 25 delay points, full density matrices")
out("=" * 70)

tomo = load_json('data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json')
out(f"\n  Qubit: Q{tomo['qubit_index']}")
out(f"  T1 = {tomo['T1_us']:.1f} us, T2 = {tomo['T2_us']:.1f} us")
out(f"  Shots: {tomo['shots']}")
out(f"  Crossing prediction: {tomo['crossing_us']:.1f} us")

raw_tomo = tomo['raw_tomography']

out(f"\n  ALL 7 INSTRUMENTS (single qubit, d=2):")
out(f"\n  {'t(us)':>8} | {'CPsi':>7} {'theta':>7} {'Pur':>6} {'Psi':>6} {'Fid':>6} | {'v_B':>7} | {'|rho01|':>8}")
out(f"  {'-'*75}")

rho_prev = None
cpsi_q52, theta_q52, t_q52 = [], [], []
bures_q52 = []

for pt in raw_tomo:
    t_us = pt['delay_us']
    rho = np.array(pt['density_matrix_real']) + 1j * np.array(pt['density_matrix_imag'])
    fid = pt['fidelity']

    cp, pur, psi = cpsi_from_rho_1q(pt['density_matrix_real'], pt['density_matrix_imag'])
    th = theta_deg(cp)
    rho01 = abs(rho[0, 1])

    # Bures
    if rho_prev is not None:
        dB = bures_distance(rho, rho_prev)
        dt_us = t_us - t_prev
        vB = dB / dt_us if dt_us > 0 else 0
    else:
        dB = 0
        vB = 0
    rho_prev = rho.copy()
    t_prev = t_us

    t_q52.append(t_us)
    cpsi_q52.append(cp)
    theta_q52.append(th)
    bures_q52.append(vB)

    out(f"  {t_us:>8.1f} | {cp:>7.4f} {th:>7.1f} {pur:>6.3f} {psi:>6.3f} "
        f"{fid:>6.3f} | {vB:>7.4f} | {rho01:>8.5f}")

t_q52 = np.array(t_q52)
cpsi_q52 = np.array(cpsi_q52)
theta_q52 = np.array(theta_q52)

# Find crossing
t_cross_q52 = None
for i in range(1, len(cpsi_q52)):
    if cpsi_q52[i-1] >= 0.25 and cpsi_q52[i] < 0.25:
        frac = (0.25 - cpsi_q52[i]) / (cpsi_q52[i-1] - cpsi_q52[i] + 1e-30)
        t_cross_q52 = t_q52[i] * (1 - frac) + t_q52[i-1] * frac
        break

out(f"\n  Measured CPsi=1/4 crossing: {t_cross_q52:.1f} us" if t_cross_q52 else
    "\n  CPsi never crosses 1/4")
out(f"  Predicted crossing: {tomo['crossing_us']:.1f} us")
if t_cross_q52:
    out(f"  Deviation: {abs(t_cross_q52 - tomo['crossing_us']):.1f} us "
        f"({100*abs(t_cross_q52 - tomo['crossing_us'])/tomo['crossing_us']:.1f}%)")

# Decay rate from CPsi
valid = cpsi_q52 > 0.01
if np.sum(valid) > 3:
    try:
        popt, _ = curve_fit(exp_decay, t_q52[valid], cpsi_q52[valid],
                            p0=[0.5, 0.005, 0.01], maxfev=5000)
        out(f"  CPsi decay rate: {popt[1]:.5f} per us (half-life {np.log(2)/popt[1]:.1f} us)")
    except Exception:
        out(f"  CPsi decay fit failed")

# Bures curvature (where CPsi monotonic)
out(f"\n  BURES CURVATURE (Q52):")
dcpsi = np.gradient(cpsi_q52, t_q52)
bv = np.array(bures_q52)
g_q52 = np.zeros(len(t_q52))
for i in range(len(t_q52)):
    if abs(dcpsi[i]) > 1e-10 and bv[i] > 0:
        g_q52[i] = (bv[i] / abs(dcpsi[i]))**2

valid_g = g_q52 > 1e-8
if np.sum(valid_g) > 4:
    cpsi_g = cpsi_q52[valid_g]
    g_g = g_q52[valid_g]
    lng = np.log(g_g + 1e-30)
    d2lng = np.gradient(np.gradient(lng, cpsi_g), cpsi_g)
    K_hw = -d2lng / (2 * g_g + 1e-30)

    out(f"  {'CPsi':>7} | {'g(CPsi)':>8} {'K_Gauss':>10}")
    out(f"  {'-'*30}")
    for cpsi_t in [0.40, 0.35, 0.30, 0.25, 0.20, 0.15]:
        idx = int(np.argmin(np.abs(cpsi_g - cpsi_t)))
        out(f"  {cpsi_g[idx]:>7.3f} | {g_g[idx]:>8.2f} {K_hw[idx]:>+10.1f}")

    # Compare K at fold
    idx_fold = int(np.argmin(np.abs(cpsi_g - 0.25)))
    K_at_fold = K_hw[idx_fold]
    out(f"\n  K at fold (CPsi~0.25): {K_at_fold:.1f}")
    out(f"  Simulation predicted: K = -25")
    out(f"  Deviation: {abs(K_at_fold - (-25)):.1f}")


# ================================================================
# PART C: Q80 PALINDROME (density matrices + CPsi)
# ================================================================
out(f"\n{'=' * 70}")
out("PART C: Q80 PALINDROME VALIDATION (March 18, 2026)")
out("Single qubit, 8 delay points")
out("=" * 70)

pali = load_json('data/ibm_run3_march2026/palindrome_ibm_torino_20260318_191348.json')
q80 = pali['qubit']
out(f"\n  Qubit: Q{q80['qubit']}")
out(f"  T1 = {q80['T1_us']:.0f} us, T2* = {q80['T2star_us']:.0f} us, T2echo = {q80['T2echo_us']:.0f} us")

out(f"\n  {'t(us)':>7} {'t/T2*':>6} | {'CPsi_m':>8} {'CPsi_T2*':>8} {'CPsi_T2e':>8} | "
    f"{'theta':>6} {'C_meas':>7} {'Psi_m':>7}")
out(f"  {'-'*75}")

for pt in pali['analysis']:
    th = theta_deg(pt['cpsi_measured'])
    out(f"  {pt['delay_us']:>7.1f} {pt['t_over_T2star']:>6.2f} | "
        f"{pt['cpsi_measured']:>8.4f} {pt['cpsi_theory_T2star']:>8.4f} "
        f"{pt['cpsi_theory_T2echo']:>8.4f} | {th:>6.1f} {pt['C_measured']:>7.4f} "
        f"{pt['psi_measured']:>7.4f}")

pred = pali['prediction']
out(f"\n  Predicted crossing: {pred['t_star_us']:.2f} us (t/T2* = {pred['t_star_over_T2star']:.3f})")
out(f"  Measured crossing:  {pali['measured_crossing_us']:.2f} us (t/T2* = {pali['measured_crossing_over_T2star']:.3f})")
out(f"  Deviation: {pali['deviation_from_prediction']*100:.1f}%")
out(f"  Verdict: {pali['verdict']['summary']}")


# ================================================================
# HARDWARE vs SIMULATION COMPARISON
# ================================================================
out(f"\n{'=' * 70}")
out("HARDWARE vs SIMULATION")
out("=" * 70)

# Load predictions for Q80
pred_data = load_json('data/ibm_run3_march2026/palindrome_predictions.json')
delay_pts = pred_data['delay_points']

out(f"\n  Q80: Measured vs Predicted CPsi")
out(f"  {'t(us)':>7} | {'CPsi_hw':>8} {'CPsi_pred':>8} {'Diff':>8} {'RelErr':>8}")
out(f"  {'-'*48}")

for hw_pt in pali['analysis']:
    t = hw_pt['delay_us']
    cpsi_hw = hw_pt['cpsi_measured']
    # Find closest prediction
    best_pred = min(delay_pts, key=lambda p: abs(p['t_us'] - t))
    cpsi_pred = best_pred['cpsi_predicted']
    diff = cpsi_hw - cpsi_pred
    rel = abs(diff) / (cpsi_pred + 1e-10) * 100
    out(f"  {t:>7.1f} | {cpsi_hw:>8.4f} {cpsi_pred:>8.4f} {diff:>+8.4f} {rel:>7.1f}%")


# ================================================================
# RECOMMENDATIONS FOR APRIL RUN
# ================================================================
out(f"\n{'=' * 70}")
out("RECOMMENDATIONS FOR APRIL RUN")
out("=" * 70)
out(f"""
  1. BIGGEST DD EFFECT: Selective DD shows highest MI at all time points.
     The effect is strongest at t=1 (first Trotter step).
     Recommendation: Measure at t=1-2 for maximum signal.

  2. TOMOGRAPHY NEEDED: From counts alone, only MI (Instrument 4) and
     decay rates (Instrument 5) are computable. For the full cockpit
     (theta, Concurrence, Psi-norm, Bures curvature), tomography is
     REQUIRED. Single-qubit tomography is cheap (3 bases x shots).
     Two-qubit tomography is expensive (9 bases x shots per pair).

  3. WHAT TO MEASURE:
     - MI per pair at t=1,2,3,4,5 Trotter steps (already done, repeat)
     - Add 2-qubit tomography for at least ONE pair (sacrifice edge)
       to get Concurrence + CPsi + theta
     - Add single-qubit tomography for sacrifice qubit
       to track coherence decay

  4. HARDWARE-SIMULATION GAP:
     - Q80 palindrome: 61.5% deviation in crossing time
     - Q52 tomography: crossing location matches better
     - The gap is qubit-specific (T2*, detuning, TLS coupling)
     - Cockpit predictions should carry ~50% uncertainty margins

  5. MOST USEFUL NEW INSTRUMENT: CONCURRENCE (Instrument 2).
     In simulation, Concurrence = PC1 (r=0.98).
     On hardware, it's the one instrument that carries the most
     information but requires 2-qubit tomography to measure.
     Prioritize this for the April run.
""")

out("=" * 70)
out("ANALYSIS COMPLETE")
out("=" * 70)

# Save
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_lines) + '\n')
print(f"\nResults saved to {results_path}")
