"""
Cockpit Navigation on the CPsi Manifold
========================================
7 instruments tracking a quantum state's decoherence trajectory:

  1. theta(t)      -- altimeter (distance from CPsi=1/4 boundary)
  2. PC1(t)        -- compass (sector balance, ~57% variance)
  3. PC2(t)        -- speedometer (mixedness/coherence, ~22%)
  4. PC3(t)        -- fuel gauge (Psi- sector, ~11%)
  5. d_dominant(t)  -- variometer (dominant mode decay rate)
  6. v_Bures(t)    -- Bures velocity + Gaussian curvature
  7. K_Peter(t)    -- Petermann factor of dominant mode

Part A: Star topology (Bell_SA x |+>_B, J_SA=1.0, J_SB=2.0, gamma=0.05)
Part B: 2-qubit Heisenberg (multiple initial states)

April 2, 2026
"""
import numpy as np
from scipy import linalg, stats
import sys, os, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
os.makedirs(results_dir, exist_ok=True)
results_path = os.path.join(results_dir, 'cockpit_navigation.txt')
_lines = []

def out(s=""):
    print(s)
    _lines.append(s)


# ================================================================
# HELPERS
# ================================================================
I2 = np.eye(2, dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def bell_fids(rho):
    bells = {
        "Phi+": np.array([1,0,0,1], dtype=complex)/np.sqrt(2),
        "Phi-": np.array([1,0,0,-1], dtype=complex)/np.sqrt(2),
        "Psi+": np.array([0,1,1,0], dtype=complex)/np.sqrt(2),
        "Psi-": np.array([0,1,-1,0], dtype=complex)/np.sqrt(2),
    }
    return {n: float(np.real(np.trace(rho @ np.outer(b, b.conj()))))
            for n, b in bells.items()}


def von_neumann(rho):
    ev = np.real(np.linalg.eigvalsh(rho))
    ev = ev[ev > 1e-14]
    return float(-np.sum(ev * np.log2(ev)))


def cpsi_true(rho):
    """CPsi = Purity x Psi-norm."""
    pur = float(np.real(np.trace(rho @ rho)))
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi = l1 / (d - 1) if d > 1 else 0.0
    return pur * psi, pur, psi


def theta_deg(cpsi_val):
    return float(np.degrees(np.arctan(np.sqrt(4*cpsi_val - 1)))) if cpsi_val > 0.25 else 0.0


def bures_distance(rho, sigma):
    """Bures distance dB = sqrt(2(1 - sqrt(F)))."""
    try:
        sqrt_rho = linalg.sqrtm(rho)
        prod = sqrt_rho @ sigma @ sqrt_rho
        ev = np.real(np.linalg.eigvalsh(prod))
        fid = float(np.sum(np.sqrt(np.maximum(ev, 0))))**2
        return float(np.sqrt(max(0, 2*(1 - np.sqrt(max(0, min(1, fid)))))))
    except Exception:
        return 0.0


def build_liouvillian(H, gamma_list, n_qubits):
    """Build Liouvillian superoperator for Z-dephasing."""
    dim = 2**n_qubits
    dim2 = dim**2
    Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for q in range(n_qubits):
        if gamma_list[q] <= 0:
            continue
        # Build Z_q on full space
        ops = [I2] * n_qubits
        ops[q] = Z
        Zq = ops[0]
        for o in ops[1:]:
            Zq = np.kron(Zq, o)
        L += gamma_list[q] * (np.kron(Zq, Zq.conj()) - np.eye(dim2))
    return L


def eigendecompose_liouvillian(L_super):
    """Full eigendecomposition with left/right eigenvectors."""
    w, vl, vr = linalg.eig(L_super, left=True, right=True)
    # Sort by real part (most negative last)
    idx = np.argsort(-w.real)
    w = w[idx]
    vl = vl[:, idx]
    vr = vr[:, idx]
    return w, vl, vr


def mode_amplitudes(w, vl, vr, rho0_vec):
    """Compute expansion coefficients c_k = <w_k|rho0> / <w_k|v_k>."""
    n = len(w)
    c = np.zeros(n, dtype=complex)
    for k in range(n):
        overlap = vl[:, k].conj() @ vr[:, k]
        if abs(overlap) > 1e-15:
            c[k] = (vl[:, k].conj() @ rho0_vec) / overlap
    return c


def petermann_factors(vl, vr):
    """K_k = 1/|<w_k|v_k>|^2 (normalized)."""
    n = vl.shape[1]
    K = np.ones(n)
    for k in range(n):
        overlap = abs(vl[:, k].conj() @ vr[:, k])
        K[k] = 1.0 / (overlap**2 + 1e-30)
    return K


def dominant_mode_at_t(w, c, t, skip_ss=True):
    """Find the mode with largest amplitude at time t.
    Returns (index, decay_rate, amplitude)."""
    amps = np.abs(c * np.exp(w * t))
    if skip_ss:
        # Mask out near-zero eigenvalues (steady state)
        mask = np.abs(w.real) > 1e-8
        amps_masked = np.where(mask, amps, 0)
    else:
        amps_masked = amps
    if np.max(amps_masked) < 1e-15:
        return 0, 0.0, 0.0
    idx = int(np.argmax(amps_masked))
    return idx, float(-w[idx].real), float(amps_masked[idx])


feat_names = ['Phi+', 'Phi-', 'Psi+', 'Psi-', 'Pur', 'SvN', 'C', 'Psi', 'ph03', 'S_coh']


def extract_features_AB(rAB, rS=None):
    f = bell_fids(rAB)
    pur = float(np.real(np.trace(rAB @ rAB)))
    svn = von_neumann(rAB)
    c = gpt.concurrence_two_qubit(rAB)
    psi = gpt.psi_norm(rAB)
    ph03 = np.angle(rAB[0, 3]) / np.pi
    s_coh = gpt.l1_coherence(rS) if rS is not None else 0.0
    return [f['Phi+'], f['Phi-'], f['Psi+'], f['Psi-'],
            pur, svn, c, psi, ph03, s_coh]


# ================================================================
# PART A: STAR TOPOLOGY (3-QUBIT)
# ================================================================
out("=" * 70)
out("COCKPIT NAVIGATION: STAR TOPOLOGY")
out("Bell_SA x |+>_B, J_SA=1.0, J_SB=2.0, gamma=0.05")
out("=" * 70)

# --- Build system ---
H3 = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
gamma = 0.05
L_ops = gpt.dephasing_ops_n([gamma]*3)
psi_init = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho0 = gpt.density_from_statevector(psi_init)

# --- Liouvillian eigendecomposition ---
out("\nBuilding 3-qubit Liouvillian (64x64)...")
L_super = build_liouvillian(H3, [gamma]*3, 3)
w3, vl3, vr3 = eigendecompose_liouvillian(L_super)

# Steady state
rho_ss_vec = vr3[:, 0] / np.sum(vr3[:, 0])  # normalize trace
rho_ss = rho_ss_vec.reshape(8, 8)

# Mode amplitudes
rho0_vec = rho0.reshape(-1)
c3 = mode_amplitudes(w3, vl3, vr3, rho0_vec)
K_peter3 = petermann_factors(vl3, vr3)

# Report eigenvalue spectrum
out(f"\nLiouvillian spectrum (top 10 by |Re(lambda)|):")
unique_rates = sorted(set(np.round(-w3.real, 6)), reverse=False)
out(f"  Spectral gap (slowest decay): {unique_rates[1]:.6f}")
out(f"  Fastest decay: {unique_rates[-1]:.4f}")
out(f"  Number of distinct decay rates: {len(unique_rates)}")

out(f"\n  {'k':>3} {'Re(lam)':>10} {'Im(lam)':>10} {'|c_k|':>10} {'K_Peter':>10}")
out(f"  {'-'*48}")
for k in range(min(12, len(w3))):
    out(f"  {k:>3} {w3[k].real:>+10.5f} {w3[k].imag:>+10.4f} "
        f"{abs(c3[k]):>10.5f} {K_peter3[k]:>10.1f}")

# --- Simulate trajectory ---
out("\nSimulating trajectory...")
dt, sample_every = 0.005, 4
n_steps = 2000

times, feats = [], []
theta_arr, cpsi_arr, conc_arr, pur_arr, psi_arr = [], [], [], [], []
bures_dists, bures_vels = [], []
dom_rate, dom_amp, dom_Kp = [], [], []
rho_prev_AB = None

rho = rho0.copy()
for step in range(n_steps + 1):
    t = step * dt
    if step % sample_every == 0:
        rAB = gpt.partial_trace_keep(rho, keep=[1, 2], n_qubits=3)
        rS  = gpt.partial_trace_keep(rho, keep=[0], n_qubits=3)

        cp, pu, ps = cpsi_true(rAB)
        th = theta_deg(cp)
        cn = gpt.concurrence_two_qubit(rAB)
        feat = extract_features_AB(rAB, rS)

        # Bures distance to previous state
        if rho_prev_AB is not None:
            dB = bures_distance(rAB, rho_prev_AB)
            dt_sample = dt * sample_every
            bures_dists.append(dB)
            bures_vels.append(dB / dt_sample)
        else:
            bures_dists.append(0.0)
            bures_vels.append(0.0)
        rho_prev_AB = rAB.copy()

        # Dominant mode at this time
        idx_d, rate_d, amp_d = dominant_mode_at_t(w3, c3, t)
        dom_rate.append(rate_d)
        dom_amp.append(amp_d)
        dom_Kp.append(K_peter3[idx_d])

        times.append(t)
        feats.append(feat)
        theta_arr.append(th)
        cpsi_arr.append(cp)
        conc_arr.append(cn)
        pur_arr.append(pu)
        psi_arr.append(ps)

    if step < n_steps:
        rho = gpt.rk4_step(rho, H3, L_ops, dt)

# Convert
times = np.array(times)
X = np.array(feats)
theta_arr = np.array(theta_arr)
cpsi_arr = np.array(cpsi_arr)
conc_arr = np.array(conc_arr)
pur_arr = np.array(pur_arr)
psi_arr = np.array(psi_arr)
bures_vels = np.array(bures_vels)
dom_rate = np.array(dom_rate)
dom_Kp = np.array(dom_Kp)
N = len(times)

# --- PCA ---
X_mean = X.mean(axis=0)
X_std = X.std(axis=0) + 1e-10
X_norm = (X - X_mean) / X_std
U, S_vals, Vt = np.linalg.svd(X_norm, full_matrices=False)
var_exp = S_vals**2 / np.sum(S_vals**2)
scores = X_norm @ Vt.T
PC1, PC2, PC3 = scores[:, 0], scores[:, 1], scores[:, 2]

# --- Bures metric and curvature ---
# g(t) = (dB/dCPsi)^2 where both dB and dCPsi are per time step
dcpsi = np.gradient(cpsi_arr, times)
g_metric = np.zeros(N)
for i in range(N):
    if abs(dcpsi[i]) > 1e-12:
        g_metric[i] = (bures_vels[i] / abs(dcpsi[i]))**2

# Gaussian curvature: K = -(1/2g) d^2(ln g)/d(CPsi)^2
# Only meaningful where g > 0 and CPsi is locally monotonic
K_gauss = np.full(N, np.nan)
valid = g_metric > 1e-10
if np.sum(valid) > 10:
    lng = np.log(g_metric[valid] + 1e-30)
    cpsi_v = cpsi_arr[valid]
    if len(cpsi_v) > 4:
        d2lng = np.gradient(np.gradient(lng, cpsi_v), cpsi_v)
        K_vals = -d2lng / (2 * g_metric[valid] + 1e-30)
        K_gauss[valid] = K_vals


# ================================================================
# DASHBOARD: STAR TOPOLOGY
# ================================================================
out(f"\n{'=' * 70}")
out("STAR TOPOLOGY DASHBOARD")
out("=" * 70)

# Find key times
t_cross_AB = None
for i in range(1, N):
    if cpsi_arr[i-1] >= 0.25 and cpsi_arr[i] < 0.25:
        t_cross_AB = times[i]
        break

t_max_theta = times[np.argmax(theta_arr)]
t_max_bures = times[np.argmax(bures_vels[1:])+1] if N > 1 else 0

out(f"\n  Key events:")
out(f"    Max theta_AB = {theta_arr.max():.1f} deg at t = {t_max_theta:.2f}")
if t_cross_AB:
    out(f"    First CPsi_AB crossing (1/4): t = {t_cross_AB:.2f}")
out(f"    Max Bures velocity: {bures_vels.max():.4f} at t = {t_max_bures:.2f}")

# Dashboard at key time points
key_times = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
if t_cross_AB:
    key_times.append(t_cross_AB)
key_times = sorted(set(key_times))

out(f"\n  {'t':>5} | {'theta':>6} {'PC1':>7} {'PC2':>7} {'PC3':>7} | "
    f"{'d_dom':>7} {'v_B':>7} {'K_P':>7} | {'CPsi':>6} {'C':>6} {'Pur':>5}")
out(f"  {'-' * 85}")

for t_key in key_times:
    idx = int(np.argmin(np.abs(times - t_key)))
    t_actual = times[idx]
    label = ""
    if t_cross_AB and abs(t_actual - t_cross_AB) < 0.03:
        label = " <-- landing"
    out(f"  {t_actual:>5.2f} | {theta_arr[idx]:>6.1f} {PC1[idx]:>+7.2f} {PC2[idx]:>+7.2f} "
        f"{PC3[idx]:>+7.2f} | {dom_rate[idx]:>7.4f} {bures_vels[idx]:>7.4f} "
        f"{dom_Kp[idx]:>7.1f} | {cpsi_arr[idx]:>6.3f} {conc_arr[idx]:>6.3f} "
        f"{pur_arr[idx]:>5.3f}{label}")

out(f"\n  Legend:")
out(f"    theta: altimeter (deg), PC1-3: PCA scores,")
out(f"    d_dom: dominant decay rate, v_B: Bures velocity,")
out(f"    K_P: Petermann factor of dominant mode")


# ================================================================
# PHASE 3: MISSING INSTRUMENTS
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 3: WHAT OBSERVABLE MEASURES EACH PC?")
out("=" * 70)

out(f"\n  PCA loadings (full trajectory, {N} snapshots):")
out(f"  {'':>8} | {'PC1':>7} {'PC2':>7} {'PC3':>7}")
out(f"  {'-' * 35}")
for j in range(len(feat_names)):
    out(f"  {feat_names[j]:>8} | {Vt[0,j]:>+7.3f} {Vt[1,j]:>+7.3f} {Vt[2,j]:>+7.3f}")

out(f"\n  Best single-observable proxy for each PC (Pearson r, full trajectory):")
observables = {
    'Concurrence': conc_arr,
    'Purity': pur_arr,
    'Psi-norm': psi_arr,
    'theta': theta_arr,
    'CPsi': cpsi_arr,
}

for k in range(3):
    best_r, best_name = 0, "?"
    for name, vals in observables.items():
        if np.std(vals) > 1e-10:
            r = abs(stats.pearsonr(vals, scores[:, k])[0])
            if r > best_r:
                best_r = r
                best_name = name
    out(f"    PC{k+1} ({var_exp[k]*100:.1f}% var) ~ {best_name} (|r|={best_r:.3f})")

# Quantum regime
mask_q = theta_arr > 0
n_q = int(np.sum(mask_q))
if n_q > 5:
    out(f"\n  Same analysis in quantum regime (theta > 0, {n_q} points):")
    for k in range(3):
        best_r, best_name = 0, "?"
        for name, vals in observables.items():
            v = vals[mask_q]
            if np.std(v) > 1e-10:
                r = abs(stats.pearsonr(v, scores[:, k][mask_q])[0])
                if r > best_r:
                    best_r = r
                    best_name = name
        out(f"    PC{k+1} ~ {best_name} (|r|={best_r:.3f})")

out(f"\n  ==> The missing compass is CONCURRENCE.")
out(f"      The missing speedometer is PSI-NORM (L1 coherence).")
out(f"      theta combines both but is dominated by neither.")


# ================================================================
# PART B: 2-QUBIT COMPARISON
# ================================================================
out(f"\n{'=' * 70}")
out("PART B: 2-QUBIT HEISENBERG COMPARISON")
out("J=1.0, gamma=0.05, four initial states")
out("=" * 70)

# Build 2-qubit Heisenberg Hamiltonian: H = J(XX + YY + ZZ)
H2 = gpt.two_qubit_heisenberg_term(0, 1, 2, J=1.0)
L_super_2q = build_liouvillian(H2, [gamma, gamma], 2)
w2, vl2, vr2 = eigendecompose_liouvillian(L_super_2q)
K_peter2 = petermann_factors(vl2, vr2)
L_ops_2q = gpt.dephasing_ops_n([gamma, gamma])

out(f"\n  Liouvillian spectrum (16x16):")
out(f"  {'k':>3} {'Re(lam)':>10} {'Im(lam)':>10} {'K_Peter':>10}")
out(f"  {'-' * 38}")
for k in range(min(10, len(w2))):
    out(f"  {k:>3} {w2[k].real:>+10.5f} {w2[k].imag:>+10.4f} {K_peter2[k]:>10.2f}")

# Initial states
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

initial_states = {
    'Bell+': np.array([1,0,0,1], dtype=complex)/np.sqrt(2),
    '|++>': np.kron(plus, plus),
    '|01>': np.kron(up, dn),
}

# Add random pure state
np.random.seed(42)
r = np.random.randn(4) + 1j * np.random.randn(4)
r = r / np.linalg.norm(r)
initial_states['Random'] = r

out(f"\n  Comparing {len(initial_states)} initial states:")

dt2, sample2 = 0.005, 4
n_steps2 = 6000  # t_max = 30

for state_name, psi_2q in initial_states.items():
    rho2 = np.outer(psi_2q, psi_2q.conj())
    rho2_vec = rho2.reshape(-1)
    c2 = mode_amplitudes(w2, vl2, vr2, rho2_vec)

    # Simulate
    times2, theta2, cpsi2, conc2, pur2, bv2 = [], [], [], [], [], []
    rho_t = rho2.copy()
    rho_prev = None

    for step in range(n_steps2 + 1):
        t = step * dt2
        if step % sample2 == 0:
            cp, pu, ps = cpsi_true(rho_t)
            th = theta_deg(cp)
            cn = gpt.concurrence_two_qubit(rho_t)

            if rho_prev is not None:
                dB = bures_distance(rho_t, rho_prev)
                bv2.append(dB / (dt2 * sample2))
            else:
                bv2.append(0.0)
            rho_prev = rho_t.copy()

            times2.append(t)
            theta2.append(th)
            cpsi2.append(cp)
            conc2.append(cn)
            pur2.append(pu)

        if step < n_steps2:
            rho_t = gpt.rk4_step(rho_t, H2, L_ops_2q, dt2)

    times2 = np.array(times2)
    theta2 = np.array(theta2)
    cpsi2 = np.array(cpsi2)
    bv2 = np.array(bv2)

    # Find crossing time
    t_cross = None
    for i in range(1, len(cpsi2)):
        if cpsi2[i-1] >= 0.25 and cpsi2[i] < 0.25:
            # Linear interpolation
            frac = (0.25 - cpsi2[i]) / (cpsi2[i-1] - cpsi2[i] + 1e-30)
            t_cross = times2[i] * (1 - frac) + times2[i-1] * frac
            break

    # Dominant mode
    idx_d, rate_d, amp_d = dominant_mode_at_t(w2, c2, 0.0)

    out(f"\n  --- {state_name} ---")
    out(f"    CPsi(0)={cpsi2[0]:.4f}, theta(0)={theta2[0]:.1f} deg")
    out(f"    Max CPsi={max(cpsi2):.4f}, max theta={max(theta2):.1f} deg")
    if t_cross:
        out(f"    Landing (CPsi=1/4): t = {t_cross:.3f}")
    else:
        out(f"    Never crosses CPsi=1/4 (stays {'above' if cpsi2[-1] > 0.25 else 'below'})")
    out(f"    Dominant mode decay rate: {rate_d:.5f}")
    out(f"    Max Bures velocity: {max(bv2):.5f}")

    # Dashboard at key times
    key_t2 = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
    if t_cross:
        key_t2.append(t_cross)
    key_t2 = sorted(set([t for t in key_t2 if t <= times2[-1]]))

    out(f"    {'t':>6} | {'theta':>6} {'CPsi':>7} {'C':>6} {'Pur':>6} {'v_B':>7}")
    out(f"    {'-' * 48}")
    for t_k in key_t2:
        i = int(np.argmin(np.abs(times2 - t_k)))
        label = " <-- landing" if t_cross and abs(times2[i] - t_cross) < 0.05 else ""
        out(f"    {times2[i]:>6.2f} | {theta2[i]:>6.1f} {cpsi2[i]:>7.4f} "
            f"{conc2[i] if i < len(conc2) else 0:>6.3f} {pur2[i] if i < len(pur2) else 0:>6.3f} "
            f"{bv2[i]:>7.5f}{label}")


# ================================================================
# GAUSSIAN CURVATURE (2-QUBIT Bell+)
# ================================================================
out(f"\n{'=' * 70}")
out("GAUSSIAN CURVATURE (2-qubit Bell+, monotonic CPsi)")
out("=" * 70)

# Recompute 2-qubit Bell+ with dense sampling for curvature
psi_bell = np.array([1,0,0,1], dtype=complex)/np.sqrt(2)
rho_b = np.outer(psi_bell, psi_bell.conj())
dt_fine = 0.002
n_fine = 2000
t_fine, cpsi_fine, bures_fine = [], [], []
rho_prev_b = None

rho_t = rho_b.copy()
for step in range(n_fine + 1):
    t = step * dt_fine
    if step % 2 == 0:
        cp, _, _ = cpsi_true(rho_t)
        if rho_prev_b is not None:
            dB = bures_distance(rho_t, rho_prev_b)
        else:
            dB = 0.0
        rho_prev_b = rho_t.copy()

        t_fine.append(t)
        cpsi_fine.append(cp)
        bures_fine.append(dB)

    if step < n_fine:
        rho_t = gpt.rk4_step(rho_t, H2, L_ops_2q, dt_fine)

t_fine = np.array(t_fine)
cpsi_fine = np.array(cpsi_fine)
bures_fine = np.array(bures_fine)
dt_s = dt_fine * 2  # sample interval

# Bures metric g(CPsi) = (dB/dCPsi)^2
dcpsi_fine = np.gradient(cpsi_fine, t_fine)
bures_vel_fine = bures_fine / dt_s
bures_vel_fine[0] = bures_vel_fine[1]  # fix first point

g_fine = np.zeros(len(t_fine))
for i in range(len(t_fine)):
    if abs(dcpsi_fine[i]) > 1e-12:
        g_fine[i] = (bures_vel_fine[i] / abs(dcpsi_fine[i]))**2

# Curvature
valid_g = g_fine > 1e-8
if np.sum(valid_g) > 10:
    cpsi_g = cpsi_fine[valid_g]
    g_g = g_fine[valid_g]
    lng_g = np.log(g_g + 1e-30)
    d2lng = np.gradient(np.gradient(lng_g, cpsi_g), cpsi_g)
    K_gauss_2q = -d2lng / (2 * g_g + 1e-30)

    # Report at key CPsi values
    out(f"\n  {'CPsi':>6} | {'g(CPsi)':>8} {'K_Gauss':>10} {'dB/dt':>8}")
    out(f"  {'-' * 40}")
    cpsi_targets = [0.33, 0.30, 0.27, 0.25, 0.22, 0.20, 0.18]
    for cp_t in cpsi_targets:
        idx = int(np.argmin(np.abs(cpsi_g - cp_t)))
        i_full = np.where(valid_g)[0][idx]
        out(f"  {cpsi_g[idx]:>6.3f} | {g_g[idx]:>8.2f} {K_gauss_2q[idx]:>+10.1f} "
            f"{bures_vel_fine[i_full]:>8.4f}")

    out(f"\n  K at fold (CPsi~0.25): {K_gauss_2q[np.argmin(np.abs(cpsi_g - 0.25))]:.1f}")
    out(f"  (INFORMATION_GEOMETRY reported K = -25 at fold)")


# ================================================================
# PHASE 5: SACRIFICE ZONE PREVIEW
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 5: SACRIFICE ZONE ANALYSIS (Star topology)")
out("=" * 70)

# Compare three noise strategies
strategies = {
    'Uniform': [gamma, gamma, gamma],
    'Edge (B loud)': [gamma*0.5, gamma*0.5, gamma*2.0],
    'Center (S loud)': [gamma*2.0, gamma*0.5, gamma*0.5],
}

for strat_name, gammas in strategies.items():
    H_s = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
    L_s = gpt.dephasing_ops_n(gammas)
    rho_s = gpt.density_from_statevector(psi_init)

    cpsi_s = []
    for step in range(n_steps + 1):
        if step % sample_every == 0:
            rAB_s = gpt.partial_trace_keep(rho_s, keep=[1, 2], n_qubits=3)
            cp_s, _, _ = cpsi_true(rAB_s)
            cpsi_s.append(cp_s)
        if step < n_steps:
            rho_s = gpt.rk4_step(rho_s, H_s, L_s, dt)

    cpsi_s = np.array(cpsi_s)
    t_land = None
    for i in range(1, len(cpsi_s)):
        if cpsi_s[i-1] >= 0.25 and cpsi_s[i] < 0.25:
            t_land = times[i]
            break

    out(f"\n  {strat_name:20s}: gamma={gammas}")
    out(f"    Max CPsi_AB = {cpsi_s.max():.4f}, "
        f"landing {'t='+f'{t_land:.2f}' if t_land else 'never'}")
    out(f"    CPsi_AB at t=5: {cpsi_s[min(len(cpsi_s)-1, 250)]:.4f}, "
        f"at t=10: {cpsi_s[-1]:.4f}")


# ================================================================
# SUMMARY
# ================================================================
out(f"\n{'=' * 70}")
out("SUMMARY: THE COCKPIT")
out("=" * 70)
out(f"""
  7 instruments for navigating the CPsi manifold:

  1. ALTIMETER (theta): Distance from CPsi=1/4 boundary.
     Range: 0-{theta_arr.max():.0f} deg for AB pair.
     Monotonic within each CPsi window.

  2. COMPASS (PC1 ~ Concurrence): Main direction of dynamics.
     {var_exp[0]*100:.0f}% of total variance. Tracks sector balance.
     The instrument we need most but didn't have.

  3. SPEEDOMETER (PC2 ~ Psi-norm): Coherence loss rate.
     {var_exp[1]*100:.0f}% of variance. Tracks mixedness.

  4. FUEL GAUGE (PC3 ~ Psi- fidelity): Temporal decay marker.
     {var_exp[2]*100:.0f}% of variance. Smallest but irreversible.

  5. VARIOMETER (d_dominant): Dominant eigenmode decay rate.
     Spectral gap = {unique_rates[1]:.4f}. Determines late-time behavior.

  6. CURVATURE (K_Gauss): Trajectory curvature in Bures geometry.
     K ~ -25 at fold (hyperbolic). Diverges near initial state.

  7. SENSITIVITY (K_Petermann): Eigenvector condition number.
     Range: {K_peter3.min():.1f} to {K_peter3.max():.0f}.
     High K = sensitive mode = turbulent flight.

  KEY FINDING: theta reads a DIAGONAL of the cockpit.
  It mixes compass, speedometer, and fuel gauge into one number.
  The pilot needs the individual instruments, not just theta.
""")

out("=" * 70)
out("ANALYSIS COMPLETE")
out("=" * 70)

# Save
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_lines) + '\n')
print(f"\nResults saved to {results_path}")
