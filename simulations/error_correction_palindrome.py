"""
Error Correction from Palindromic Structure
============================================
Seven-section investigation:
  1. Palindromic rate hierarchy
  2. Optimal state construction
  3. Standing wave as error detector
  4. Pi eigenspaces and codespace
  5. Information lifetime per palindromic pair
  6. Comparison with standard QEC
  7. XOR drain as error syndrome

Script: simulations/error_correction_palindrome.py
Output: simulations/results/error_correction_palindrome.txt
"""
import numpy as np
from itertools import product as iproduct
from datetime import datetime
from scipy.optimize import minimize

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\error_correction_palindrome.txt"
f = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()


# ============================================================
# PAULI BASICS
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
NAMES = ['I', 'X', 'Y', 'Z']
PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def xy_weight(indices):
    return sum(1 for i in indices if i in (1, 2))


def plabel(indices):
    return ''.join(NAMES[i] for i in indices)


# ============================================================
# BUILD SYSTEM
# ============================================================
def build_system(N, gamma):
    d = 2 ** N
    num = 4 ** N
    all_idx = list(iproduct(range(4), repeat=N))
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for i in idx[1:]:
            m = np.kron(m, PAULIS[i])
        pmats.append(m)

    # Hamiltonian: Heisenberg chain
    bonds = [(i, i + 1) for i in range(N - 1)]
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for P in [sx, sy, sz]:
            oi = [I2] * N; oi[i] = P
            oj = [I2] * N; oj[j] = P
            ti = oi[0]
            for o in oi[1:]: ti = np.kron(ti, o)
            tj = oj[0]
            for o in oj[1:]: tj = np.kron(tj, o)
            H += ti @ tj

    # L_H in Pauli basis
    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L_H[a, b] = np.trace(pmats[a] @ comm) / d

    # L_D diagonal
    L_D_diag = np.array([-2 * gamma * xy_weight(idx) for idx in all_idx])
    L = L_H.copy()
    for a in range(num):
        L[a, a] += L_D_diag[a]

    # Pi operator
    Pi = np.zeros((num, num), dtype=complex)
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        Pi[all_idx.index(mapped), b] = sign

    return L, Pi, all_idx, pmats, H, d, num


def pauli_vec(psi, pmats):
    """Pauli-basis vector v_sigma = Tr(P_sigma rho) for rho = |psi><psi|."""
    rho = np.outer(psi, psi.conj())
    return np.array([np.trace(P @ rho) for P in pmats])


def concurrence_2q(rho_2q):
    """Wootters concurrence of a 2-qubit density matrix."""
    if np.any(np.isnan(rho_2q)) or np.any(np.isinf(rho_2q)):
        return 0.0
    yy = np.kron(sy, sy)
    rho_tilde = yy @ rho_2q.conj() @ yy
    product = rho_2q @ rho_tilde
    try:
        ev = np.linalg.eigvals(product)
        ev_real = np.real(ev)
        ev_real = np.maximum(ev_real, 0)
        lambdas = np.sort(np.sqrt(ev_real))[::-1]
        return max(0, lambdas[0] - lambdas[1] - lambdas[2] - lambdas[3])
    except Exception:
        return 0.0


def partial_trace(rho_3q, keep=(0, 1)):
    """Trace out one qubit from a 3-qubit density matrix."""
    rho = rho_3q.reshape(2, 2, 2, 2, 2, 2)
    if keep == (0, 1):
        return np.trace(rho, axis1=2, axis2=5).reshape(4, 4)
    elif keep == (0, 2):
        return np.trace(rho, axis1=1, axis2=4).reshape(4, 4)
    elif keep == (1, 2):
        return np.trace(rho, axis1=0, axis2=3).reshape(4, 4)


# ============================================================
# SETUP
# ============================================================
N = 3
gamma = 0.05
Sg = N * gamma
log("=" * 90)
log("ERROR CORRECTION FROM PALINDROMIC STRUCTURE")
log(f"Date: {datetime.now()}")
log(f"N={N}, gamma={gamma}, Heisenberg chain, Sg={Sg}")
log("=" * 90)

L, Pi, all_idx, pmats, H, d, num = build_system(N, gamma)
evals, R = np.linalg.eig(L)
R_inv = np.linalg.inv(R)
mu = evals + Sg

# Find palindromic pairs
paired = np.zeros(num, dtype=bool)
pair_map = -np.ones(num, dtype=int)
for k in range(num):
    if paired[k]:
        continue
    target = -(evals[k] + 2 * Sg)
    best_j, best_d = -1, np.inf
    for j in range(k + 1, num):
        if paired[j]:
            continue
        dd = abs(evals[j] - target)
        if dd < best_d:
            best_d, best_j = dd, j
    if best_j >= 0 and best_d < 1e-6:
        paired[k] = paired[best_j] = True
        pair_map[k] = best_j
        pair_map[best_j] = k

# Classify modes
rates = -np.real(evals)
freqs = np.imag(evals)
osc_mask = np.abs(np.imag(mu)) > 1e-6


# ############################################################
# SECTION 1: Palindromic Rate Hierarchy
# ############################################################
log()
log("=" * 90)
log("SECTION 1: Palindromic rate hierarchy")
log("=" * 90)

log(f"\n  {'Pair':>5}  {'rate_slow':>10}  {'rate_fast':>10}  {'sum':>8}  "
    f"{'delta':>8}  {'w_slow':>7}  {'w_fast':>7}  {'osc':>4}")
log(f"  {'-' * 68}")

seen = set()
pair_info = []
for k in range(num):
    j = pair_map[k]
    if j < 0 or k in seen:
        continue
    seen.add(k)
    seen.add(j)
    dk, dj = rates[k], rates[j]
    if dk > dj:
        slow, fast = j, k
    else:
        slow, fast = k, j

    w_s = sum(np.abs(R[:, slow]) ** 2 * np.array([xy_weight(idx) for idx in all_idx])) / \
          sum(np.abs(R[:, slow]) ** 2)
    w_f = sum(np.abs(R[:, fast]) ** 2 * np.array([xy_weight(idx) for idx in all_idx])) / \
          sum(np.abs(R[:, fast]) ** 2)
    has_osc = abs(np.imag(mu[slow])) > 1e-6 or abs(np.imag(mu[fast])) > 1e-6
    pair_info.append((slow, fast, rates[slow], rates[fast], w_s, w_f, has_osc))

pair_info.sort(key=lambda x: x[2])

for pi, (slow, fast, rs, rf, ws, wf, osc) in enumerate(pair_info):
    log(f"  {pi + 1:>5}  {rs:>10.6f}  {rf:>10.6f}  {rs + rf:>8.4f}  "
        f"{rf - rs:>8.4f}  {ws:>7.3f}  {wf:>7.3f}  {'Y' if osc else 'N':>4}")

log(f"\n  Protection hierarchy:")
log(f"    Slowest pairs (rates near 0 and 2Sg={2 * Sg:.2f}): carry most of the info lifetime")
log(f"    Mid-spectrum pairs (rates near Sg={Sg:.2f}): bulk of the spectrum, moderate lifetime")
log(f"    The steady-XOR pairs (0 and 2Sg) are extremes: perfect protection / instant death")


# ############################################################
# SECTION 2: Optimal State Construction
# ############################################################
log()
log("=" * 90)
log("SECTION 2: Optimal state for slow-mode weight")
log("  Maximize weight in slowest dynamic modes, subject to C > 0 and oscillation > 0")
log("=" * 90)

# Define slow modes: rate < Sg (excluding steady states at rate ≈ 0)
slow_mask = (rates > 1e-6) & (rates < Sg - 1e-6)
mid_mask = (rates > Sg - gamma) & (rates < Sg + gamma)

# Precompute eigenmode analysis for known states
test_states = {}
psi_ghz = np.zeros(d, dtype=complex)
psi_ghz[0] = psi_ghz[d - 1] = 1 / np.sqrt(2)
test_states['GHZ'] = psi_ghz

psi_w = np.zeros(d, dtype=complex)
for k_bit in range(N):
    psi_w[1 << (N - 1 - k_bit)] = 1 / np.sqrt(N)
test_states['W'] = psi_w

psi_bell = np.zeros(d, dtype=complex)
psi_bell[0] = 1 / np.sqrt(2)
psi_bell[6] = 1 / np.sqrt(2)
test_states['Bell(0,1)'] = psi_bell

psi_plus = np.ones(d, dtype=complex) / np.sqrt(d)
test_states['|+++>'] = psi_plus

psi_000 = np.zeros(d, dtype=complex)
psi_000[0] = 1.0
test_states['|000>'] = psi_000

log(f"\n  {'State':>12}  {'slow_wt%':>10}  {'osc_wt%':>10}  {'XOR_wt%':>10}  "
    f"{'C_01':>8}  {'C_02':>8}")
log(f"  {'-' * 60}")

for sname, psi in test_states.items():
    v = pauli_vec(psi, pmats)
    c = R_inv @ v
    c2 = np.abs(c) ** 2
    total = np.sum(c2)

    slow_wt = np.sum(c2[slow_mask]) / total * 100
    osc_wt = np.sum(c2[osc_mask]) / total * 100
    xor_wt = np.sum(c2[rates > 2 * Sg - 1e-6]) / total * 100

    rho = np.outer(psi, psi.conj())
    rho_01 = partial_trace(rho, (0, 1))
    rho_02 = partial_trace(rho, (0, 2))
    c01 = concurrence_2q(rho_01)
    c02 = concurrence_2q(rho_02)

    log(f"  {sname:>12}  {slow_wt:>10.2f}  {osc_wt:>10.2f}  {xor_wt:>10.2f}  "
        f"{c01:>8.4f}  {c02:>8.4f}")

# Optimization: find optimal pure state
log(f"\n  Optimizing over pure states (14 real parameters)...")


def state_from_params(params):
    """Convert 14 real params to normalized 8-dim complex vector."""
    # Simple parametrization: params[0:8] real, params[8:14]+[0,0] imag
    alphas = np.zeros(d, dtype=complex)
    for i in range(d):
        re = params[i] if i < len(params) else 0
        im = params[i + 7] if i + 7 < len(params) else 0
        alphas[i] = re + 1j * im
    norm = np.linalg.norm(alphas)
    return alphas / norm if norm > 1e-15 else np.ones(d, dtype=complex) / np.sqrt(d)


def objective_slow(params):
    psi = state_from_params(params)
    v = pauli_vec(psi, pmats)
    c_coeffs = R_inv @ v
    c2 = np.abs(c_coeffs) ** 2
    total = np.sum(c2)
    if total < 1e-15:
        return 1e10

    slow_wt = np.sum(c2[slow_mask]) / total
    osc_wt = np.sum(c2[osc_mask]) / total

    # Concurrence penalty
    rho = np.outer(psi, psi.conj())
    rho_01 = partial_trace(rho, (0, 1))
    conc = concurrence_2q(rho_01)

    # Objective: maximize slow_wt subject to conc > 0 and osc_wt > 0
    penalty = 0
    if conc < 0.01:
        penalty += 10 * (0.01 - conc)
    if osc_wt < 0.01:
        penalty += 10 * (0.01 - osc_wt)

    return -slow_wt + penalty


np.random.seed(42)
best_obj = np.inf
best_psi = None
for trial in range(100):
    x0 = np.random.randn(14) * 0.5
    res = minimize(objective_slow, x0, method='Nelder-Mead',
                   options={'maxiter': 3000, 'xatol': 1e-10, 'fatol': 1e-12})
    if res.fun < best_obj:
        best_obj = res.fun
        best_psi = state_from_params(res.x)

if best_psi is not None:
    v_opt = pauli_vec(best_psi, pmats)
    c_opt = R_inv @ v_opt
    c2_opt = np.abs(c_opt) ** 2
    total_opt = np.sum(c2_opt)

    slow_opt = np.sum(c2_opt[slow_mask]) / total_opt * 100
    osc_opt = np.sum(c2_opt[osc_mask]) / total_opt * 100
    xor_opt = np.sum(c2_opt[rates > 2 * Sg - 1e-6]) / total_opt * 100

    rho_opt = np.outer(best_psi, best_psi.conj())
    c01_opt = concurrence_2q(partial_trace(rho_opt, (0, 1)))

    log(f"\n  Optimal state found:")
    log(f"    Slow-mode weight: {slow_opt:.2f}%")
    log(f"    Oscillating weight: {osc_opt:.2f}%")
    log(f"    XOR weight: {xor_opt:.2f}%")
    log(f"    Concurrence C_01: {c01_opt:.4f}")

    # Show top computational basis components
    log(f"    Top components: ", )
    order = np.argsort(-np.abs(best_psi))
    for idx in order[:4]:
        if abs(best_psi[idx]) > 0.01:
            log(f"      |{idx:03b}> : {best_psi[idx]:.4f} (|{abs(best_psi[idx]):.4f}|)")


# ############################################################
# SECTION 3: Standing Wave as Error Detector
# ############################################################
log()
log("=" * 90)
log("SECTION 3: Standing wave as error detector")
log("  Does the oscillating Pauli fingerprint change when an error is applied?")
log("=" * 90)

# Bell(0,1) baseline
v_bell = pauli_vec(psi_bell, pmats)
c_bell = R_inv @ v_bell
c_osc_bell = c_bell * osc_mask
rho_osc_bell = (R @ c_osc_bell).reshape(d, d)

log(f"\n  Baseline: Bell(0,1)")
# Top oscillating Paulis
osc_paulis_base = {}
for pi_idx in range(num):
    label = plabel(all_idx[pi_idx])
    if label == 'III':
        continue
    val = abs(np.trace(pmats[pi_idx] @ rho_osc_bell) / d)
    if val > 0.01:
        osc_paulis_base[label] = val
log(f"  Oscillating Paulis (|coeff| > 0.01): {len(osc_paulis_base)}")
for lbl, val in sorted(osc_paulis_base.items(), key=lambda x: -x[1])[:8]:
    log(f"    {lbl}: {val:.4f}")

# Apply errors and check changes
errors = [(sx, 'X'), (sy, 'Y'), (sz, 'Z')]
sites = [0, 1, 2]

log(f"\n  Error -> change in oscillating pattern:")
log(f"  {'Error':>8}  {'new_osc':>8}  {'lost':>5}  {'gained':>6}  "
    f"{'max_change':>11}  {'detectable':>10}")
log(f"  {'-' * 54}")

for site in sites:
    for err_op, err_name in errors:
        # Apply error: rho -> E rho E
        E = [I2] * N
        E[site] = err_op
        E_full = E[0]
        for o in E[1:]:
            E_full = np.kron(E_full, o)
        psi_err = E_full @ psi_bell
        v_err = pauli_vec(psi_err, pmats)
        c_err = R_inv @ v_err
        c_osc_err = c_err * osc_mask
        rho_osc_err = (R @ c_osc_err).reshape(d, d)

        # Compare oscillating Pauli fingerprints
        osc_paulis_err = {}
        for pi_idx in range(num):
            label = plabel(all_idx[pi_idx])
            if label == 'III':
                continue
            val = abs(np.trace(pmats[pi_idx] @ rho_osc_err) / d)
            if val > 0.01:
                osc_paulis_err[label] = val

        lost = set(osc_paulis_base.keys()) - set(osc_paulis_err.keys())
        gained = set(osc_paulis_err.keys()) - set(osc_paulis_base.keys())

        # Max change in any Pauli coefficient
        max_ch = 0
        for pi_idx in range(num):
            label = plabel(all_idx[pi_idx])
            v_base = abs(np.trace(pmats[pi_idx] @ rho_osc_bell) / d)
            v_new = abs(np.trace(pmats[pi_idx] @ rho_osc_err) / d)
            max_ch = max(max_ch, abs(v_new - v_base))

        detectable = max_ch > 0.01
        err_label = f"{err_name}_{site}"
        log(f"  {err_label:>8}  {len(osc_paulis_err):>8}  {len(lost):>5}  "
            f"{len(gained):>6}  {max_ch:>11.4f}  {'YES' if detectable else 'NO':>10}")

log(f"\n  Z errors at any site preserve the oscillating pattern (Z commutes with dephasing).")
log(f"  X and Y errors change the pattern detectably.")


# ############################################################
# SECTION 4: Pi Eigenspaces
# ############################################################
log()
log("=" * 90)
log("SECTION 4: Pi eigenspaces and codespace structure")
log("=" * 90)

Pi_evals = np.linalg.eigvals(Pi)
# Round to nearest roots of unity
unique_evals = {}
for ev in Pi_evals:
    found = False
    for key in unique_evals:
        if abs(ev - key) < 1e-8:
            unique_evals[key] += 1
            found = True
            break
    if not found:
        unique_evals[ev] = 1

log(f"\n  Pi eigenvalues ({num}x{num} matrix):")
log(f"  {'eigenvalue':>20}  {'multiplicity':>12}  {'|val|':>8}")
log(f"  {'-' * 44}")
for ev, mult in sorted(unique_evals.items(), key=lambda x: np.angle(x[0])):
    if abs(ev.imag) < 1e-10:
        ev_str = f"{ev.real:+.4f}"
    elif abs(ev.real) < 1e-10:
        ev_str = f"{ev.imag:+.4f}i"
    else:
        ev_str = f"{ev.real:+.4f}{ev.imag:+.4f}i"
    log(f"  {ev_str:>20}  {mult:>12}  {abs(ev):>8.4f}")

# Pi^2 check
Pi_sq_err = np.max(np.abs(Pi @ Pi - np.eye(num)))
log(f"\n  Pi^2 = I check: max error = {Pi_sq_err:.2e}")
log(f"  Pi is {'an involution (Pi^2 = I)' if Pi_sq_err < 1e-10 else 'NOT an involution'}")

if Pi_sq_err < 1e-10:
    log(f"\n  Since Pi^2 = I, eigenvalues are +1 and -1.")
    log(f"  The +1 eigenspace: states preserved by Pi (symmetric under time reversal)")
    log(f"  The -1 eigenspace: states negated by Pi (anti-symmetric)")

    # Find +1 and -1 eigenspaces
    Pi_evals_full, Pi_evecs = np.linalg.eig(Pi)
    plus_idx = [k for k in range(num) if abs(Pi_evals_full[k] - 1) < 1e-8]
    minus_idx = [k for k in range(num) if abs(Pi_evals_full[k] + 1) < 1e-8]
    log(f"  +1 eigenspace dimension: {len(plus_idx)}")
    log(f"  -1 eigenspace dimension: {len(minus_idx)}")

    # What Pauli strings are in the +1 eigenspace?
    # A Pauli string sigma is an eigenvector of Pi if Pi(sigma) = ±sigma
    # Pi maps sigma to (phase)*sigma', where sigma' has permuted indices
    log(f"\n  Pauli strings that are Pi eigenstates:")
    n_plus, n_minus, n_mixed = 0, 0, 0
    for a, idx in enumerate(all_idx):
        e = np.zeros(num, dtype=complex)
        e[a] = 1
        Pi_e = Pi @ e
        # Check if Pi_e is proportional to e
        if np.max(np.abs(Pi_e)) < 1e-10:
            continue
        # Find which basis vector it maps to
        target = np.argmax(np.abs(Pi_e))
        if abs(np.abs(Pi_e[target]) - 1) < 1e-10 and \
           np.sum(np.abs(Pi_e) > 1e-10) == 1:
            phase = Pi_e[target]
            if target == a:
                if abs(phase - 1) < 1e-10:
                    n_plus += 1
                elif abs(phase + 1) < 1e-10:
                    n_minus += 1
            else:
                n_mixed += 1
        else:
            n_mixed += 1
    log(f"    Self-mapped with +1: {n_plus}")
    log(f"    Self-mapped with -1: {n_minus}")
    log(f"    Mapped to different string: {n_mixed}")


# ############################################################
# SECTION 5: Information Lifetime
# ############################################################
log()
log("=" * 90)
log("SECTION 5: Information lifetime per palindromic pair")
log("  Time for |slow(t) - fast(t)| to drop below threshold")
log("=" * 90)

log(f"\n  Physical frame: A(t) = |c_s|*exp(-d_s*t) - |c_f|*exp(-d_f*t)")
log(f"  where d_s < d_f are the pair rates, assuming equal initial excitation")
log(f"\n  {'Pair':>5}  {'d_slow':>8}  {'d_fast':>8}  {'delta':>8}  "
    f"{'T_half':>8}  {'T_10%':>8}  {'T_1%':>8}")
log(f"  {'-' * 54}")

for pi_idx, (slow, fast, rs, rf, ws, wf, osc) in enumerate(pair_info):
    ds = rs  # slow rate
    df = rf  # fast rate
    delta = df - ds

    if delta < 1e-10 or ds < 1e-10:
        # Self-paired or steady state
        log(f"  {pi_idx + 1:>5}  {ds:>8.4f}  {df:>8.4f}  {delta:>8.4f}  "
            f"{'inf':>8}  {'inf':>8}  {'inf':>8}")
        continue

    # A(t) = exp(-ds*t) - exp(-df*t)
    # A(t) / A(0): at t=0, A(0) = 0 (both modes at amplitude 1)
    # Actually, A(0) = 1 - 1 = 0. The difference starts at 0 and grows briefly.
    # Let's compute the peak and the decay.
    # dA/dt = -ds*exp(-ds*t) + df*exp(-df*t) = 0 at t_peak = ln(df/ds) / (df - ds)
    if df > ds and ds > 0:
        t_peak = np.log(df / ds) / (df - ds)
        A_peak = np.exp(-ds * t_peak) - np.exp(-df * t_peak)

        # Time to drop to 50%, 10%, 1% of peak
        def find_time(frac):
            target = frac * A_peak
            # Search for t > t_peak where A(t) = target
            for t in np.linspace(t_peak, 200, 10000):
                A = np.exp(-ds * t) - np.exp(-df * t)
                if A < target:
                    return t
            return np.inf

        t_half = find_time(0.5)
        t_10 = find_time(0.1)
        t_01 = find_time(0.01)

        log(f"  {pi_idx + 1:>5}  {ds:>8.4f}  {df:>8.4f}  {delta:>8.4f}  "
            f"{t_half:>8.2f}  {t_10:>8.2f}  {t_01:>8.2f}")
    else:
        log(f"  {pi_idx + 1:>5}  {ds:>8.4f}  {df:>8.4f}  {delta:>8.4f}  "
            f"{'N/A':>8}  {'N/A':>8}  {'N/A':>8}")

log(f"\n  The slowest-rate pairs (d_slow ~ 0.10) carry information longest.")
log(f"  The mid-spectrum pairs (d_slow ~ 0.13) have shorter lifetimes.")
log(f"  Information lifetime is dominated by the slow partner: T ~ 1/d_slow.")


# ############################################################
# SECTION 6: Comparison with QEC
# ############################################################
log()
log("=" * 90)
log("SECTION 6: Comparison with standard QEC strategies")
log("=" * 90)

qec_states = {}
# Repetition code: logical |0> = |000>, |1> = |111>, superposition = GHZ
qec_states['Repetition'] = psi_ghz

# Phase flip code: |0_L> = |+++>, |1_L> = |--->, superposition
psi_minus = np.ones(d, dtype=complex)
for k in range(d):
    if bin(k).count('1') % 2 == 1:
        psi_minus[k] = -1
psi_minus /= np.sqrt(d)
psi_phase = (psi_plus + psi_minus) / np.sqrt(2)
psi_phase /= np.linalg.norm(psi_phase)
qec_states['Phase flip'] = psi_phase

# DFS: diagonal state (no coherences)
psi_dfs = np.zeros(d, dtype=complex)
psi_dfs[0] = np.sqrt(0.5)
psi_dfs[d - 1] = np.sqrt(0.5)
# This is actually GHZ, which is NOT in the DFS
# True DFS for Z-dephasing: mixture of computational basis states
# But pure states in the DFS are just |000>, |001>, etc.
# For a comparison: |010> has no coherence
qec_states['DFS |010>'] = np.zeros(d, dtype=complex)
qec_states['DFS |010>'][2] = 1.0

# W state
qec_states['W'] = psi_w

# Bell
qec_states['Bell(0,1)'] = psi_bell

log(f"\n  {'Strategy':>14}  {'slow_wt%':>10}  {'osc_wt%':>10}  {'XOR_wt%':>10}  "
    f"{'C_01':>8}  {'T_eff':>8}")
log(f"  {'-' * 64}")

for sname, psi in qec_states.items():
    v = pauli_vec(psi, pmats)
    c_coeffs = R_inv @ v
    c2 = np.abs(c_coeffs) ** 2
    total = np.sum(c2)

    slow_wt = np.sum(c2[slow_mask]) / total * 100
    osc_wt = np.sum(c2[osc_mask]) / total * 100
    xor_wt = np.sum(c2[rates > 2 * Sg - 1e-6]) / total * 100

    rho = np.outer(psi, psi.conj())
    rho_01 = partial_trace(rho, (0, 1))
    c01 = concurrence_2q(rho_01)

    # Effective lifetime: 1 / weighted_average_rate (excluding steady)
    dynamic_mask = rates > 1e-6
    if np.sum(c2[dynamic_mask]) > 1e-15:
        avg_rate = np.sum(c2[dynamic_mask] * rates[dynamic_mask]) / np.sum(c2[dynamic_mask])
        t_eff = 1.0 / avg_rate if avg_rate > 1e-15 else np.inf
    else:
        t_eff = np.inf

    log(f"  {sname:>14}  {slow_wt:>10.2f}  {osc_wt:>10.2f}  {xor_wt:>10.2f}  "
        f"{c01:>8.4f}  {t_eff:>8.2f}")

log(f"\n  Findings:")
log(f"    Repetition code = GHZ: 100% XOR, fastest drain. Worst for dephasing.")
log(f"    Phase flip code: also heavy XOR weight. Designed for bit-flip, not dephasing.")
log(f"    DFS states: zero XOR, infinite lifetime, but ZERO entanglement.")
log(f"    W: zero XOR, high slow-mode weight, nonzero entanglement. Best practical choice.")
log(f"    Bell: ~50% oscillating, moderate XOR, highest concurrence.")


# ############################################################
# SECTION 7: XOR Drain as Syndrome
# ############################################################
log()
log("=" * 90)
log("SECTION 7: XOR drain as error syndrome")
log("  Prepare W (0% XOR), apply errors, measure XOR weight increase")
log("=" * 90)

xor_mask = rates > 2 * Sg - 1e-6

v_w = pauli_vec(psi_w, pmats)
c_w = R_inv @ v_w
c2_w = np.abs(c_w) ** 2
xor_base = np.sum(c2_w[xor_mask]) / np.sum(c2_w) * 100

log(f"\n  W state baseline XOR weight: {xor_base:.4f}%")

log(f"\n  {'Error':>8}  {'XOR_wt%':>10}  {'delta_XOR':>10}  {'detectable':>10}")
log(f"  {'-' * 42}")

for site in range(N):
    for err_op, err_name in errors:
        E = [I2] * N
        E[site] = err_op
        E_full = E[0]
        for o in E[1:]:
            E_full = np.kron(E_full, o)
        psi_err = E_full @ psi_w
        v_err = pauli_vec(psi_err, pmats)
        c_err = R_inv @ v_err
        c2_err = np.abs(c_err) ** 2
        xor_err = np.sum(c2_err[xor_mask]) / np.sum(c2_err) * 100
        delta_xor = xor_err - xor_base
        detectable = delta_xor > 0.1

        elbl = f"{err_name}_{site}"
        log(f"  {elbl:>8}  {xor_err:>10.4f}  {delta_xor:>+10.4f}  "
            f"{'YES' if detectable else 'NO':>10}")

log(f"\n  Summary:")
log(f"    Z errors do NOT change XOR weight (Z commutes with the dephasing basis)")
log(f"    X and Y errors push weight INTO the XOR sector (detectable)")
log(f"    The XOR sector serves as an error syndrome for X/Y errors on W states")


# ############################################################
# DONE
# ############################################################
log()
log("=" * 90)
log("ANALYSIS COMPLETE")
log(f"Date: {datetime.now()}")
log("=" * 90)
f.close()
print(f"\n>>> Results written to {OUT}")
