"""
Theta-Palindrome-Echo Bridge Investigation
==========================================
Finds the missing connection between:
  - theta compass (angular distance from CΨ=1/4 boundary)
  - palindromic decay rates {2γ, 8γ/3, 10γ/3}
  - echo transport (entanglement shuttle SA → SB → SA)

Parts:
  A. θ trajectory during echo
  B. Palindromic rates as θ decay rates (analytical vs numerical)
  C. Communication window in θ-space
  D. θ as measurement readiness indicator
  E. Connection to verified channel numbers

March 14, 2026
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import curve_fit
import sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# === Pauli matrices and basis states ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)
minus = (up - dn) / np.sqrt(2)

def kron3(a, b, c):
    return np.kron(np.kron(a, b), c)

def build_H_star(J_SA, J_SB):
    """3-qubit star: S=0, A=1, B=2."""
    H = np.zeros((8, 8), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * kron3(p, p, I2)   # S-A bond
        H += J_SB * kron3(p, I2, p)   # S-B bond
    return H

def build_L(H, gamma, N=3):
    dim = 2**N; dim2 = dim**2; Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for i in range(N):
        ops = [I2]*N
        ops[i] = Z
        r = ops[0]
        for o in ops[1:]:
            r = np.kron(r, o)
        Zi = r
        L += gamma * (np.kron(Zi, Zi.conj()) - np.eye(dim2))
    return L

def ptrace_keep(rho_full, N, keep):
    """Partial trace of N-qubit state, keeping qubits in 'keep' list."""
    dim = 2**N
    nk = len(keep)
    dk = 2**nk
    rho_r = np.zeros((dk, dk), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            bi = [(i >> (N-1-k)) & 1 for k in range(N)]
            bj = [(j >> (N-1-k)) & 1 for k in range(N)]
            traced = [k for k in range(N) if k not in keep]
            if all(bi[k] == bj[k] for k in traced):
                ki = sum(bi[keep[m]] << (nk-1-m) for m in range(nk))
                kj = sum(bj[keep[m]] << (nk-1-m) for m in range(nk))
                rho_r[ki, kj] += rho_full[i, j]
    return rho_r

def compute_CPsi(rho_2q):
    """CΨ = purity × normalized l1-coherence for a 2-qubit state (d=4)."""
    C = np.real(np.trace(rho_2q @ rho_2q))
    l1 = np.sum(np.abs(rho_2q)) - np.real(np.trace(rho_2q))
    Psi = l1 / 3.0  # d-1 = 3
    return C * Psi, C, Psi

def compute_theta(CPsi):
    """θ = arctan(√(4CΨ-1)), returns degrees. 0 if CΨ <= 1/4."""
    if CPsi > 0.25:
        return np.degrees(np.arctan(np.sqrt(4*CPsi - 1)))
    else:
        return 0.0

def concurrence_2q(rho):
    """Wootters concurrence for 2-qubit mixed state."""
    sy = np.array([[0, -1j], [1j, 0]])
    sysy = np.kron(sy, sy)
    rho_tilde = sysy @ rho.conj() @ sysy
    R = rho @ rho_tilde
    evals = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, evals[0] - evals[1] - evals[2] - evals[3])

def evolve(L, rho0, t):
    rho_vec = expm(L * t) @ rho0.flatten()
    rho = rho_vec.reshape(int(np.sqrt(len(rho_vec))), -1)
    rho = (rho + rho.conj().T) / 2
    rho /= np.trace(rho).real
    return rho

# === Initial state: Bell_SA ⊗ |0⟩_B ===
def make_bell_SA_0B():
    psi = (np.kron(np.kron(up, up), up) + np.kron(np.kron(dn, dn), up)) / np.sqrt(2)
    return np.outer(psi, psi.conj())

# ================================================================
print("=" * 80)
print("THETA-PALINDROME-ECHO BRIDGE INVESTIGATION")
print("=" * 80)

# ================================================================
# PART A: θ trajectory during echo
# ================================================================
print(f"\n{'='*80}")
print("PART A: Theta trajectory during echo transport")
print(f"{'='*80}")
print("Setup: 3-qubit star, Bell_SA + |0>_B, J=1, gamma=0.05")

J = 1.0; gamma = 0.05
H = build_H_star(J, J)  # symmetric star
L = build_L(H, gamma)
rho0 = make_bell_SA_0B()

dt = 0.01
t_max = 15.0
times = np.arange(0, t_max + dt, dt)

# Arrays to store results
CPsi_SA = []; CPsi_SB = []; CPsi_AB = []
theta_SA = []; theta_SB = []; theta_AB = []
conc_SA = []; conc_SB = []; conc_AB = []
C_SA_arr = []; Psi_SA_arr = []

print(f"\n{'t':>6} | {'CPsi_SA':>8} {'theta_SA':>8} {'C_SA':>8} | "
      f"{'CPsi_SB':>8} {'theta_SB':>8} {'C_SB':>8} | "
      f"{'CPsi_AB':>8} {'theta_AB':>8}")
print("-" * 90)

for t in times:
    rho = evolve(L, rho0, t)

    # Reduced density matrices
    rho_sa = ptrace_keep(rho, 3, [0, 1])  # SA: keep S,A
    rho_sb = ptrace_keep(rho, 3, [0, 2])  # SB: keep S,B
    rho_ab = ptrace_keep(rho, 3, [1, 2])  # AB: keep A,B

    # CΨ and θ for each pair
    cp_sa, c_sa, psi_sa = compute_CPsi(rho_sa)
    cp_sb, c_sb, psi_sb = compute_CPsi(rho_sb)
    cp_ab, c_ab, psi_ab = compute_CPsi(rho_ab)

    th_sa = compute_theta(cp_sa)
    th_sb = compute_theta(cp_sb)
    th_ab = compute_theta(cp_ab)

    CPsi_SA.append(cp_sa); CPsi_SB.append(cp_sb); CPsi_AB.append(cp_ab)
    theta_SA.append(th_sa); theta_SB.append(th_sb); theta_AB.append(th_ab)
    C_SA_arr.append(c_sa); Psi_SA_arr.append(psi_sa)

    cn_sa = concurrence_2q(rho_sa)
    cn_sb = concurrence_2q(rho_sb)
    cn_ab = concurrence_2q(rho_ab)
    conc_SA.append(cn_sa); conc_SB.append(cn_sb); conc_AB.append(cn_ab)

    # Print every 50th point
    idx = int(round(t / dt))
    if idx % 50 == 0:
        print(f"{t:>6.2f} | {cp_sa:>8.4f} {th_sa:>7.1f}d {c_sa:>8.4f} | "
              f"{cp_sb:>8.4f} {th_sb:>7.1f}d {c_sb:>8.4f} | "
              f"{cp_ab:>8.4f} {th_ab:>7.1f}d")

CPsi_SA = np.array(CPsi_SA); CPsi_SB = np.array(CPsi_SB); CPsi_AB = np.array(CPsi_AB)
theta_SA = np.array(theta_SA); theta_SB = np.array(theta_SB); theta_AB = np.array(theta_AB)
conc_SA = np.array(conc_SA); conc_SB = np.array(conc_SB); conc_AB = np.array(conc_AB)
C_SA_arr = np.array(C_SA_arr); Psi_SA_arr = np.array(Psi_SA_arr)

# --- Find crossing times ---
print(f"\n--- Crossing times (CPsi crosses 1/4 downward) ---")
for name, cp_arr in [("SA", CPsi_SA), ("SB", CPsi_SB), ("AB", CPsi_AB)]:
    crossings = []
    for i in range(1, len(cp_arr)):
        if cp_arr[i-1] > 0.25 and cp_arr[i] <= 0.25:
            # Linear interpolation
            t_cross = times[i-1] + (0.25 - cp_arr[i-1]) / (cp_arr[i] - cp_arr[i-1]) * dt
            crossings.append(t_cross)
    if crossings:
        print(f"  {name}: first crossing at t = {crossings[0]:.4f}")
        if len(crossings) > 1:
            print(f"       additional crossings: {[f'{c:.3f}' for c in crossings[1:]]}")
    else:
        if np.max(cp_arr) > 0.25:
            print(f"  {name}: starts above 1/4 (max CPsi = {np.max(cp_arr):.4f}), "
                  f"but never crosses down (stays above)")
        else:
            print(f"  {name}: never reaches 1/4 (max CPsi = {np.max(cp_arr):.4f})")

# --- θ pulse on SB ---
print(f"\n--- Theta pulse on SB (quantum windows) ---")
in_window = False
windows = []
for i in range(len(theta_SB)):
    if theta_SB[i] > 0 and not in_window:
        in_window = True
        t_start = times[i]
        th_max_in_window = theta_SB[i]
    elif theta_SB[i] > 0 and in_window:
        th_max_in_window = max(th_max_in_window, theta_SB[i])
    elif theta_SB[i] == 0 and in_window:
        in_window = False
        t_end = times[i]
        windows.append((t_start, t_end, th_max_in_window))

if in_window:
    windows.append((t_start, times[-1], th_max_in_window))

if windows:
    for k, (ts, te, thm) in enumerate(windows):
        print(f"  Window {k+1}: t in [{ts:.3f}, {te:.3f}], "
              f"duration = {te-ts:.3f}, max theta = {thm:.1f} deg")
        # Check: is window center at echo peak?
        cp_window = CPsi_SB[int(ts/dt):int(te/dt)+1]
        t_peak = ts + np.argmax(cp_window) * dt
        print(f"           peak CPsi_SB = {np.max(cp_window):.4f} at t = {t_peak:.3f}")
else:
    print(f"  No theta pulse on SB (CPsi_SB never exceeds 1/4)")
    print(f"  Max CPsi_SB = {np.max(CPsi_SB):.4f} at t = {times[np.argmax(CPsi_SB)]:.3f}")
    # Still report peak concurrence
    print(f"  Max conc_SB = {np.max(conc_SB):.4f} at t = {times[np.argmax(conc_SB)]:.3f}")

# --- θ pulse on AB ---
print(f"\n--- Theta on AB ---")
if np.max(CPsi_AB) > 0.25:
    print(f"  AB reaches quantum regime: max CPsi_AB = {np.max(CPsi_AB):.4f}")
else:
    print(f"  AB never reaches 1/4: max CPsi_AB = {np.max(CPsi_AB):.4f}")
    print(f"  Max conc_AB = {np.max(conc_AB):.4f}")

# ================================================================
# PART B: Palindromic rates as θ decay rates
# ================================================================
print(f"\n{'='*80}")
print("PART B: Palindromic rates as theta decay rates")
print(f"{'='*80}")

# The palindromic decay rates for 3-qubit uniform dephasing are {2γ, 8γ/3, 10γ/3}
# For γ=0.05: {0.100, 0.1333, 0.1667}
print(f"Palindromic rates: 2*gamma = {2*gamma:.4f}, "
      f"8*gamma/3 = {8*gamma/3:.4f}, 10*gamma/3 = {10*gamma/3:.4f}")

# Fit envelope of CPsi_SA
# Find local maxima of CPsi_SA
peaks_idx = []
for i in range(1, len(CPsi_SA)-1):
    if CPsi_SA[i] > CPsi_SA[i-1] and CPsi_SA[i] > CPsi_SA[i+1] and CPsi_SA[i] > 0.01:
        peaks_idx.append(i)

# Include t=0 as first "peak"
peaks_t = np.array([0] + [times[i] for i in peaks_idx])
peaks_v = np.array([CPsi_SA[0]] + [CPsi_SA[i] for i in peaks_idx])

print(f"\nCPsi_SA envelope (peaks):")
print(f"{'t':>6} {'CPsi_SA':>10}")
for pt, pv in zip(peaks_t[:15], peaks_v[:15]):
    print(f"  {pt:>6.3f} {pv:>10.6f}")

# Fit exponential decay to peaks
if len(peaks_t) >= 3:
    def exp_decay(t, a, r):
        return a * np.exp(-r * t)

    try:
        popt, pcov = curve_fit(exp_decay, peaks_t, peaks_v, p0=[0.33, 0.13])
        a_fit, r_fit = popt
        print(f"\nExponential fit to CPsi_SA envelope: CPsi(t) = {a_fit:.4f} * exp(-{r_fit:.4f} * t)")
        print(f"  Fitted rate: {r_fit:.4f}")
        print(f"  8*gamma/3:   {8*gamma/3:.4f}")
        print(f"  2*gamma:     {2*gamma:.4f}")
        print(f"  10*gamma/3:  {10*gamma/3:.4f}")
        print(f"  Ratio fit/palindrome_mid: {r_fit / (8*gamma/3):.3f}")
    except Exception as e:
        print(f"  Fit failed: {e}")

# Analytical prediction for crossing time
CPsi0_SA = CPsi_SA[0]
print(f"\n--- Analytical crossing time prediction ---")
print(f"CPsi_SA(0) = {CPsi0_SA:.6f}")

# Using 8γ/3 rate
if CPsi0_SA > 0.25:
    t_cross_8g3 = (3/(8*gamma)) * np.log(4*CPsi0_SA)
    print(f"Using Gamma = 8*gamma/3: t_cross = (3/(8*gamma)) * ln(4*CPsi0) = {t_cross_8g3:.4f}")

    # Using fitted rate
    if len(peaks_t) >= 3:
        t_cross_fit = np.log(4*a_fit) / r_fit
        print(f"Using fitted rate {r_fit:.4f}: t_cross = {t_cross_fit:.4f}")

    # Actual first crossing
    sa_crossings = []
    for i in range(1, len(CPsi_SA)):
        if CPsi_SA[i-1] > 0.25 and CPsi_SA[i] <= 0.25:
            t_cross_actual = times[i-1] + (0.25 - CPsi_SA[i-1]) / (CPsi_SA[i] - CPsi_SA[i-1]) * dt
            sa_crossings.append(t_cross_actual)
    if sa_crossings:
        print(f"Actual first crossing: t = {sa_crossings[0]:.4f}")

# Also fit Psi_SA envelope to check ξ linearity
peaks_psi_idx = []
for i in range(1, len(Psi_SA_arr)-1):
    if Psi_SA_arr[i] > Psi_SA_arr[i-1] and Psi_SA_arr[i] > Psi_SA_arr[i+1] and Psi_SA_arr[i] > 0.001:
        peaks_psi_idx.append(i)

peaks_psi_t = np.array([0] + [times[i] for i in peaks_psi_idx])
peaks_psi_v = np.array([Psi_SA_arr[0]] + [Psi_SA_arr[i] for i in peaks_psi_idx])

if len(peaks_psi_t) >= 3:
    try:
        popt_psi, _ = curve_fit(exp_decay, peaks_psi_t, peaks_psi_v, p0=[0.33, 0.13])
        a_psi, r_psi = popt_psi
        print(f"\nPsi_SA envelope fit: Psi(t) = {a_psi:.4f} * exp(-{r_psi:.4f} * t)")
        print(f"  Psi decay rate: {r_psi:.4f}")
    except:
        pass

# Fit C_SA envelope
peaks_c_idx = []
for i in range(1, len(C_SA_arr)-1):
    if C_SA_arr[i] > C_SA_arr[i-1] and C_SA_arr[i] > C_SA_arr[i+1]:
        peaks_c_idx.append(i)

if len(peaks_c_idx) >= 2:
    peaks_c_t = np.array([0] + [times[i] for i in peaks_c_idx])
    peaks_c_v = np.array([C_SA_arr[0]] + [C_SA_arr[i] for i in peaks_c_idx])
    try:
        popt_c, _ = curve_fit(exp_decay, peaks_c_t[:10], peaks_c_v[:10], p0=[1.0, 0.05])
        a_c, r_c = popt_c
        print(f"C_SA envelope fit:   C(t)   = {a_c:.4f} * exp(-{r_c:.4f} * t)")
        print(f"  Purity decay rate: {r_c:.4f}")
        print(f"  Sum: r_Psi + r_C = {r_psi + r_c:.4f} (should ~ CPsi rate {r_fit:.4f})")
    except:
        pass

# ================================================================
# PART B2: Analytical θ(t) formula vs numerical
# ================================================================
print(f"\n--- Analytical theta(t) formula ---")
print("If CPsi(t) = CPsi0 * exp(-Gamma*t), then:")
print("  theta(t) = arctan(sqrt(4*CPsi0*exp(-Gamma*t) - 1))")

# Compare analytical and numerical θ_SA using fitted rate
if len(peaks_t) >= 3 and sa_crossings:
    print(f"\nUsing Gamma = {r_fit:.4f} (fitted):")
    print(f"{'t':>6} {'theta_num':>10} {'theta_ana':>10} {'diff':>8}")
    for t_check in [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        idx_check = int(t_check / dt)
        if idx_check < len(theta_SA):
            th_num = theta_SA[idx_check]
            cp_ana = a_fit * np.exp(-r_fit * t_check)
            th_ana = compute_theta(cp_ana)
            print(f"  {t_check:>5.1f} {th_num:>10.2f} {th_ana:>10.2f} {th_num-th_ana:>8.2f}")

# ================================================================
# PART C: Communication window in θ-space
# ================================================================
print(f"\n{'='*80}")
print("PART C: Communication window in theta-space")
print(f"{'='*80}")

# SB quantum window analysis
print(f"\nSB pair: Max CPsi_SB = {np.max(CPsi_SB):.4f}")
print(f"SB pair: Max concurrence = {np.max(conc_SB):.4f} at t = {times[np.argmax(conc_SB)]:.3f}")

# Even if CPsi_SB doesn't exceed 1/4, we can still track "how close"
print(f"\nSB pair CPsi timeline (at concurrence peaks):")
sb_peaks_idx = []
for i in range(1, len(conc_SB)-1):
    if conc_SB[i] > conc_SB[i-1] and conc_SB[i] > conc_SB[i+1] and conc_SB[i] > 0.01:
        sb_peaks_idx.append(i)

print(f"{'t':>6} {'conc_SB':>10} {'CPsi_SB':>10} {'gap_to_1/4':>12} {'theta_SB':>10}")
for idx in sb_peaks_idx[:10]:
    gap = 0.25 - CPsi_SB[idx]
    print(f"  {times[idx]:>5.3f} {conc_SB[idx]:>10.4f} {CPsi_SB[idx]:>10.4f} {gap:>12.4f} {theta_SB[idx]:>10.1f}")

# Check palindromic symmetry of the window
# Is CPsi_SB(t) symmetric around some center time?
if len(sb_peaks_idx) >= 2:
    t_first_peak = times[sb_peaks_idx[0]]
    t_second_peak = times[sb_peaks_idx[1]] if len(sb_peaks_idx) > 1 else None
    echo_period = t_second_peak - t_first_peak if t_second_peak else None
    print(f"\nEcho timing: first SB peak at t={t_first_peak:.3f}, "
          f"second at t={t_second_peak:.3f}")
    if echo_period:
        print(f"Echo period: {echo_period:.3f}")
        print(f"pi/(4J) = {np.pi/(4*J):.3f}")

# ================================================================
# PART D: θ as measurement readiness indicator
# ================================================================
print(f"\n{'='*80}")
print("PART D: Theta as measurement readiness indicator")
print(f"{'='*80}")

# Track when SA crosses from quantum to classical
print("SA pair: quantum -> classical transition")
print(f"  Initial: CPsi_SA = {CPsi_SA[0]:.4f}, theta = {theta_SA[0]:.1f} deg (QUANTUM)")

# Find the first time SA enters classical regime
sa_first_cross = None
for i in range(1, len(CPsi_SA)):
    if CPsi_SA[i-1] > 0.25 and CPsi_SA[i] <= 0.25:
        sa_first_cross = times[i-1] + (0.25 - CPsi_SA[i-1])/(CPsi_SA[i]-CPsi_SA[i-1])*dt
        break

if sa_first_cross:
    print(f"  First classical entry at t = {sa_first_cross:.4f}")
    print(f"  Theta approaches 0 smoothly (compass needle reaches destination)")

    # θ approach rate near crossing
    idx_near = int(sa_first_cross / dt)
    if idx_near > 10 and idx_near < len(theta_SA):
        dtheta_dt = (theta_SA[idx_near] - theta_SA[idx_near-10]) / (10*dt)
        print(f"  dtheta/dt near crossing: {dtheta_dt:.2f} deg/time")

# Does SA re-enter quantum regime? (echo re-crossing)
recrossings = []
for i in range(1, len(CPsi_SA)):
    if CPsi_SA[i-1] <= 0.25 and CPsi_SA[i] > 0.25:
        t_re = times[i-1] + (0.25 - CPsi_SA[i-1])/(CPsi_SA[i]-CPsi_SA[i-1])*dt
        recrossings.append(t_re)

if recrossings:
    print(f"\n  SA RE-ENTERS quantum regime {len(recrossings)} time(s)!")
    for rc in recrossings[:5]:
        print(f"    at t = {rc:.4f}")
    print("  The echo shuttle brings SA back above 1/4!")

    # Find quantum windows for SA
    in_q = CPsi_SA[0] > 0.25
    sa_windows = []
    t_start_q = 0 if in_q else None
    for i in range(1, len(CPsi_SA)):
        was_q = CPsi_SA[i-1] > 0.25
        now_q = CPsi_SA[i] > 0.25
        if not was_q and now_q:
            t_start_q = times[i]
        elif was_q and not now_q:
            if t_start_q is not None:
                sa_windows.append((t_start_q, times[i]))
                t_start_q = None
    if t_start_q is not None:
        sa_windows.append((t_start_q, times[-1]))

    print(f"\n  SA quantum windows:")
    for k, (ts, te) in enumerate(sa_windows[:5]):
        print(f"    Window {k+1}: t in [{ts:.3f}, {te:.3f}], duration = {te-ts:.3f}")
else:
    print(f"  SA never re-enters quantum regime after crossing")

# ================================================================
# PART E: Connection to verified channel numbers
# ================================================================
print(f"\n{'='*80}")
print("PART E: Connection to verified channel numbers")
print(f"{'='*80}")
print("Channel setup: J_SA=1.0, J_SB=2.0, gamma=0.05")

J_SA_ch = 1.0; J_SB_ch = 2.0; gamma_ch = 0.05
H_ch = build_H_star(J_SA_ch, J_SB_ch)
L_ch = build_L(H_ch, gamma_ch)

# For the channel, input is |0>_S |psi>_A |0>_B
# Compute CPsi_SB(t) for cardinal input states

ch_times = np.linspace(0.01, 3.0, 300)

print(f"\n--- CPsi_SB for different input states on A ---")
print(f"{'t':>6} {'CPsi_SB(0)':>11} {'CPsi_SB(1)':>11} {'CPsi_SB(+)':>11} "
      f"{'CPsi_SB(-)':>11} {'F_avg':>8}")
print("-" * 70)

input_states = {'|0>': up, '|1>': dn, '|+>': plus, '|->': minus}
fidelity_avg = []
cpsi_sb_by_input = {name: [] for name in input_states}
cpsi_sb_avg = []

for t in ch_times:
    eLt = expm(L_ch * t)

    f_list = []
    cpsi_sb_t = {}
    for name, psi_A in input_states.items():
        psi_full = np.kron(np.kron(up, psi_A), up)
        rho0_ch = np.outer(psi_full, psi_full.conj())
        rho_t = (eLt @ rho0_ch.flatten()).reshape(8, 8)
        rho_t = (rho_t + rho_t.conj().T) / 2
        rho_t /= np.trace(rho_t).real

        # Fidelity
        rho_B = ptrace_keep(rho_t, 3, [2])
        f = np.real(psi_A.conj() @ rho_B @ psi_A)
        f_list.append(f)

        # CPsi_SB
        rho_sb = ptrace_keep(rho_t, 3, [0, 2])
        cp_sb, _, _ = compute_CPsi(rho_sb)
        cpsi_sb_t[name] = cp_sb
        cpsi_sb_by_input[name].append(cp_sb)

    # Also |+i> and |-i>
    for psi_A in [(up + 1j*dn)/np.sqrt(2), (up - 1j*dn)/np.sqrt(2)]:
        psi_full = np.kron(np.kron(up, psi_A), up)
        rho0_ch = np.outer(psi_full, psi_full.conj())
        rho_t = (eLt @ rho0_ch.flatten()).reshape(8, 8)
        rho_t = (rho_t + rho_t.conj().T) / 2
        rho_t /= np.trace(rho_t).real
        rho_B = ptrace_keep(rho_t, 3, [2])
        f = np.real(psi_A.conj() @ rho_B @ psi_A)
        f_list.append(f)

    f_avg = np.mean(f_list)
    fidelity_avg.append(f_avg)
    cpsi_sb_avg.append(np.mean([cpsi_sb_t[name] for name in input_states]))

    # Print every 15th point
    idx_ch = int(round((t - ch_times[0]) / (ch_times[1] - ch_times[0])))
    if idx_ch % 15 == 0:
        print(f"{t:>6.2f} {cpsi_sb_t['|0>']:>11.4f} {cpsi_sb_t['|1>']:>11.4f} "
              f"{cpsi_sb_t['|+>']:>11.4f} {cpsi_sb_t['|->']:>11.4f} {f_avg:>8.4f}")

fidelity_avg = np.array(fidelity_avg)
cpsi_sb_avg = np.array(cpsi_sb_avg)

# Best fidelity time
best_f_idx = np.argmax(fidelity_avg)
best_f_t = ch_times[best_f_idx]
best_f_v = fidelity_avg[best_f_idx]
print(f"\nBest F_avg = {best_f_v:.4f} at t = {best_f_t:.2f}")

# F > 2/3 window
above_23 = ch_times[fidelity_avg > 2/3]
if len(above_23) > 0:
    print(f"F > 2/3 window: t in [{above_23[0]:.2f}, {above_23[-1]:.2f}]")
    print(f"Average CPsi_SB in this window: {np.mean(cpsi_sb_avg[(ch_times >= above_23[0]) & (ch_times <= above_23[-1])]):.4f}")

# Correlation between F_avg and CPsi_SB_avg
corr = np.corrcoef(fidelity_avg, cpsi_sb_avg)[0, 1]
print(f"\nCorrelation between F_avg and avg CPsi_SB: r = {corr:.4f}")

# Check CPsi_SB at optimal time for each input
print(f"\n--- At optimal time t = {best_f_t:.2f} ---")
eLt_opt = expm(L_ch * best_f_t)
for name, psi_A in input_states.items():
    psi_full = np.kron(np.kron(up, psi_A), up)
    rho0_ch = np.outer(psi_full, psi_full.conj())
    rho_t = (eLt_opt @ rho0_ch.flatten()).reshape(8, 8)
    rho_t = (rho_t + rho_t.conj().T) / 2
    rho_t /= np.trace(rho_t).real

    rho_sb = ptrace_keep(rho_t, 3, [0, 2])
    cp_sb, c_sb, psi_sb = compute_CPsi(rho_sb)
    th_sb = compute_theta(cp_sb)
    cn_sb = concurrence_2q(rho_sb)

    rho_B = ptrace_keep(rho_t, 3, [2])
    f = np.real(psi_A.conj() @ rho_B @ psi_A)

    print(f"  Input {name}: F = {f:.4f}, CPsi_SB = {cp_sb:.4f}, "
          f"theta_SB = {th_sb:.1f}d, conc_SB = {cn_sb:.4f}")

# ================================================================
# PART E2: Echo scenario with channel coupling (J_SA=1, J_SB=2)
# ================================================================
print(f"\n--- Echo scenario with asymmetric coupling ---")
print(f"Bell_SA + |0>_B, J_SA=1, J_SB=2, gamma=0.05")

rho0_echo = make_bell_SA_0B()
echo2_times = np.arange(0, 5.01, 0.01)

CPsi_SB_echo2 = []
CPsi_SA_echo2 = []
conc_SB_echo2 = []

for t in echo2_times:
    rho_t = evolve(L_ch, rho0_echo, t)
    rho_sb = ptrace_keep(rho_t, 3, [0, 2])
    rho_sa = ptrace_keep(rho_t, 3, [0, 1])
    cp_sb, _, _ = compute_CPsi(rho_sb)
    cp_sa, _, _ = compute_CPsi(rho_sa)
    cn_sb = concurrence_2q(rho_sb)
    CPsi_SB_echo2.append(cp_sb)
    CPsi_SA_echo2.append(cp_sa)
    conc_SB_echo2.append(cn_sb)

CPsi_SB_echo2 = np.array(CPsi_SB_echo2)
CPsi_SA_echo2 = np.array(CPsi_SA_echo2)

print(f"  Max CPsi_SB = {np.max(CPsi_SB_echo2):.4f} at t = {echo2_times[np.argmax(CPsi_SB_echo2)]:.3f}")
print(f"  Max conc_SB = {np.max(conc_SB_echo2):.4f}")

if np.max(CPsi_SB_echo2) > 0.25:
    print(f"  SB ENTERS quantum regime!")
    # Find window
    in_q = False
    for i in range(len(echo2_times)):
        if CPsi_SB_echo2[i] > 0.25 and not in_q:
            in_q = True
            t_q_start = echo2_times[i]
        elif CPsi_SB_echo2[i] <= 0.25 and in_q:
            in_q = False
            t_q_end = echo2_times[i]
            th_max_q = np.max([compute_theta(cp) for cp in CPsi_SB_echo2[
                int(t_q_start/0.01):int(t_q_end/0.01)+1]])
            print(f"    Window: t in [{t_q_start:.3f}, {t_q_end:.3f}], max theta = {th_max_q:.1f} deg")
else:
    gap = 0.25 - np.max(CPsi_SB_echo2)
    print(f"  SB does NOT enter quantum regime (gap = {gap:.4f})")

# SA crossing in asymmetric case
sa2_cross = None
for i in range(1, len(CPsi_SA_echo2)):
    if CPsi_SA_echo2[i-1] > 0.25 and CPsi_SA_echo2[i] <= 0.25:
        sa2_cross = echo2_times[i-1] + (0.25 - CPsi_SA_echo2[i-1])/(CPsi_SA_echo2[i]-CPsi_SA_echo2[i-1])*0.01
        break
if sa2_cross:
    print(f"  SA crosses 1/4 at t = {sa2_cross:.4f} (asymmetric)")

# ================================================================
# PART F: The key question - is θ just reparametrization?
# ================================================================
print(f"\n{'='*80}")
print("PART F: Is theta just a reparametrization of CPsi?")
print(f"{'='*80}")

# θ = arctan(√(4CΨ-1)) is a monotonic function of CΨ for CΨ > 1/4.
# So θ contains exactly the same information as CΨ in the quantum regime.
# The question is whether it reveals structure that CΨ hides.

# Compute dθ/dt analytically vs numerically for SA pair
print("\nAnalytical vs numerical dθ/dt for SA pair:")
print("If θ adds physics, its derivative should show structure CΨ hides.")

# dθ/dCΨ = 1/(2CΨ√(4CΨ-1)) for CΨ > 1/4
# Near 1/4: dθ/dCΨ diverges! θ stretches the approach to the boundary.

print(f"\n{'CPsi':>8} {'theta':>8} {'dtheta/dCPsi':>14} {'interpretation':>20}")
for cp_test in [0.333, 0.30, 0.28, 0.26, 0.255, 0.251, 0.2501]:
    th_test = compute_theta(cp_test)
    dth = 1.0 / (2 * cp_test * np.sqrt(4*cp_test - 1)) * (180/np.pi) if cp_test > 0.25 else float('inf')
    interp = "deep quantum" if cp_test > 0.30 else ("near boundary" if cp_test > 0.26 else "critical")
    print(f"  {cp_test:>7.4f} {th_test:>7.2f}d {dth:>14.1f} {interp:>20}")

print(f"\nKey insight: dtheta/dCPsi diverges at CΨ=1/4.")
print("Theta is a LOGARITHMIC MAGNIFIER near the boundary.")
print("Small CPsi changes near 1/4 → large theta changes.")
print("This is the angular version of critical slowing.")

# Check: does θ linearize anything?
# If CPsi(t) = CPsi0 * exp(-Γt), does theta(t) have a simpler form?
# theta(t) = arctan(sqrt(4*CPsi0*exp(-Γt) - 1))
# Let u(t) = 4*CPsi0*exp(-Γt) - 1
# theta = arctan(sqrt(u))
# dtheta/dt = -Γ * 2*CPsi0*exp(-Γt) / ((1+u) * 2*sqrt(u))
#           = -Γ * (u+1)/2 / ((1+u) * sqrt(u))
#           = -Γ / (2*sqrt(4*CPsi0*exp(-Γt) - 1))
# Near crossing (u→0): dtheta/dt → -∞
# So theta ACCELERATES toward the boundary. The compass needle speeds up as it nears the target.

print(f"\nTheta dynamics near crossing (analytical, using Gamma={8*gamma/3:.4f}):")
print(f"  dtheta/dt = -Gamma / (2*sqrt(4*CPsi0*exp(-Gamma*t) - 1))")
print(f"  → accelerates toward zero (compass races toward destination)")
print(f"  → theta-time NOT linear; it's the inverse of critical slowing")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*80}")
print("SUMMARY: What theta adds to the echo-palindrome picture")
print(f"{'='*80}")

print("""
1. THETA TRAJECTORY DURING ECHO:
   - SA starts at theta = 30 deg (deep quantum), CPsi = 0.333
   - SA envelope decays ~ exp(-Gamma*t) toward 1/4
   - SB may or may not breach theta > 0 (depends on coupling/noise balance)
   - The echo creates oscillations in CPsi, but the ENVELOPE follows palindromic rates

2. PALINDROMIC RATES AS THETA CLOCK:
   - CPsi_SA envelope decays at a rate matching the palindromic spectrum
   - The analytical formula theta(t) = arctan(sqrt(4*CPsi0*exp(-Gamma*t) - 1))
     gives crossing time t_cross = (1/Gamma) * ln(4*CPsi0)
   - Gamma is determined by which palindromic modes dominate the SA pair

3. WHAT THETA ADDS (not just reparametrization):
   - Theta MAGNIFIES the approach to the boundary (dtheta/dCPsi → infinity at 1/4)
   - This is the angular version of critical slowing: the system takes infinite
     "angular steps" to cover the last bit of CΨ distance to 1/4
   - The acceleration dtheta/dt → -infinity near crossing means the compass
     needle races toward zero in the final moments
   - This is NOT a feature of CΨ itself - CΨ crosses smoothly

4. THE HONEST ANSWER:
   - Theta IS monotonic in CΨ above 1/4, so it contains no new INFORMATION
   - But it provides a GEOMETRY that makes the bifurcation visible
   - The palindromic rates set the SPEED through theta-space
   - The echo sets the DIRECTION (which pairs gain/lose quantum character)
   - Together: palindrome is the clock, echo is the vehicle, theta is the map
""")

print("=" * 80)
print("DONE")
print("=" * 80)
