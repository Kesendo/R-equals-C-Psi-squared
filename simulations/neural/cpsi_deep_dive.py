#!/usr/bin/env python3
"""
CΨ deep dive: two leads from the first exploration.

Lead 1: CΨ_5 (boundary amplitude) crosses 1/4 near alpha=0.25.
         Is this crossing point parameter-independent?

Lead 2: Quantum CΨ = purity * coherence (a PRODUCT of two factors).
         None of the five candidates are products. Try:
         CΨ_6 = (1 - palindromic residual) * (E_freq / E_total)
         = "palindromic quality" * "oscillation fraction"
         This is zero when either factor vanishes.

Lead 3: The quantum V-Effect lives in d^2 space (correlations).
         Wilson-Cowan activities are d-dimensional.
         The CORRELATION matrix C = <x x^T> is d^2-dimensional.
         Build the correlation Jacobian L_C = J x I + I x J^T
         and look for palindromic structure there.
"""
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def build_linear_jacobian(W, signs, tau_E, tau_I, alpha):
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_E if signs[i] > 0 else tau_I
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def make_balanced_network(N, n_exc, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    inh_idx = rng.choice(N, N - n_exc, replace=False)
    signs[inh_idx] = -1
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (N, N))
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if mask[i, j]:
                W[i, j] = signs[j] * weights[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def cpsi_5_boundary(ev, evec, signs):
    """Boundary amplitude: oscillation energy at E/I interface."""
    e_idx = np.where(signs > 0)[0]
    i_idx = np.where(signs < 0)[0]
    boundary_energy = 0
    total_energy = 0
    for k in range(len(ev)):
        osc = abs(ev[k].imag)
        vec = evec[:, k]
        amp_E = np.sum(np.abs(vec[e_idx])**2)
        amp_I = np.sum(np.abs(vec[i_idx])**2)
        total_amp = amp_E + amp_I
        if total_amp > 0:
            boundary = 4 * amp_E * amp_I / (total_amp**2)
            boundary_energy += osc * boundary
            total_energy += osc
    return boundary_energy / total_energy if total_energy > 0 else 0


def find_crossing(alphas, values, target=0.25):
    """Find alpha where values crosses target (linear interpolation)."""
    for i in range(len(values) - 1):
        if (values[i] < target and values[i+1] >= target) or \
           (values[i] > target and values[i+1] <= target):
            frac = (target - values[i]) / (values[i+1] - values[i])
            return alphas[i] + frac * (alphas[i+1] - alphas[i])
    return None


# ================================================================
# LEAD 1: Is CΨ_5 = 1/4 crossing point parameter-independent?
# ================================================================
print("=" * 65)
print("LEAD 1: CΨ_5 crossing point stability")
print("Does CΨ_5 = 1/4 always occur at the same alpha?")
print("=" * 65)

alphas = np.linspace(0.01, 2.0, 100)

print(f"\n{'Config':>30s}  {'alpha at CΨ_5=1/4':>18s}")
print("-" * 55)

for N, n_exc, density, tau_E, tau_I, seed in [
    (50, 25, 0.3, 5.0, 10.0, 42),
    (50, 25, 0.3, 5.0, 15.0, 42),    # tau_ratio 3.0
    (50, 25, 0.3, 5.0, 20.0, 42),    # tau_ratio 4.0
    (50, 25, 0.3, 10.0, 20.0, 42),   # different tau_E
    (50, 25, 0.1, 5.0, 10.0, 42),    # lower density
    (50, 25, 0.5, 5.0, 10.0, 42),    # higher density
    (20, 10, 0.3, 5.0, 10.0, 42),    # smaller N
    (100, 50, 0.3, 5.0, 10.0, 42),   # larger N
    (50, 25, 0.3, 5.0, 10.0, 99),    # different network
    (50, 25, 0.3, 5.0, 10.0, 7),     # different network
]:
    W, signs = make_balanced_network(N, n_exc, density, seed)
    values = []
    for a in alphas:
        J = build_linear_jacobian(W, signs, tau_E, tau_I, a)
        ev = np.linalg.eigvals(J)
        _, evec = np.linalg.eig(J)
        values.append(cpsi_5_boundary(ev, evec, signs))

    crossing = find_crossing(alphas, values, 0.25)
    label = f"N={N} E={n_exc} d={density} t={tau_I/tau_E:.0f}x s={seed}"
    if crossing is not None:
        print(f"  {label:>30s}  alpha = {crossing:.4f}")
    else:
        print(f"  {label:>30s}  no crossing")


# ================================================================
# LEAD 2: Product candidate CΨ_6
# ================================================================
print("\n" + "=" * 65)
print("LEAD 2: Product candidate CΨ_6 = quality * oscillation")
print("Analog of quantum CΨ = purity * coherence")
print("=" * 65)

def palindrome_quality(J, signs):
    """1 - normalized palindrome residual. Range [0, 1]."""
    n = len(signs)
    e_local = list(np.where(signs > 0)[0])
    i_local = list(np.where(signs < 0)[0])
    n_pairs = min(len(e_local), len(i_local))
    perm = np.arange(n)
    for k in range(n_pairs):
        perm[e_local[k]] = i_local[k]
        perm[i_local[k]] = e_local[k]
    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, perm[i]] = 1.0

    QJQ = Q @ J @ Q.T
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2 * np.diag(S_diag)
    R_off = R - np.diag(np.diag(R))
    norm_J = np.linalg.norm(J)
    residual = np.linalg.norm(R_off) / norm_J if norm_J > 0 else 1.0
    return max(0, 1.0 - residual)


def cpsi_6_product(J, signs):
    """CΨ_6 = palindrome_quality * (E_freq / E_total).
    Product of two factors, like quantum CΨ = purity * coherence."""
    ev = np.linalg.eigvals(J)
    E_freq = np.sum(np.abs(np.imag(ev)))
    E_decay = np.sum(np.abs(np.real(ev)))
    osc_fraction = E_freq / (E_freq + E_decay) if (E_freq + E_decay) > 0 else 0

    quality = palindrome_quality(J, signs)
    return quality * osc_fraction


N, n_exc = 50, 25
tau_E, tau_I = 5.0, 10.0
W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)

print(f"\n{'alpha':>6s}  {'quality':>8s}  {'osc_frac':>8s}  {'CΨ_6':>8s}  {'cross?':>6s}")
print("-" * 50)

for alpha in [0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]:
    J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)
    ev = np.linalg.eigvals(J)
    E_freq = np.sum(np.abs(np.imag(ev)))
    E_decay = np.sum(np.abs(np.real(ev)))
    osc = E_freq / (E_freq + E_decay) if (E_freq + E_decay) > 0 else 0
    qual = palindrome_quality(J, signs)
    c6 = qual * osc

    cross = " <-- 1/4" if abs(c6 - 0.25) < 0.02 else ""
    print(f"  {alpha:4.2f}  {qual:8.4f}  {osc:8.4f}  {c6:8.4f}{cross}")


# ================================================================
# LEAD 3: Correlation space (d^2 dimensional)
# ================================================================
print("\n" + "=" * 65)
print("LEAD 3: Palindrome in correlation space (N^2 dimensional)")
print("The quantum palindrome lives in d^2, not d.")
print("=" * 65)

# Use smaller N for d^2 feasibility
N_small = 10
n_exc_small = 5
W_s, signs_s = make_balanced_network(N_small, n_exc_small, density=0.3, seed=42)

print(f"\nN={N_small}, correlation space dimension = {N_small**2}")

for alpha in [0.1, 0.3, 0.5, 1.0, 2.0]:
    J_s = build_linear_jacobian(W_s, signs_s, tau_E, tau_I, alpha)

    # Correlation Liouvillian: L_C = J x I + I x J^T
    # This governs dC/dt = J*C + C*J^T (without noise driving term)
    eye_n = np.eye(N_small)
    L_C = np.kron(J_s, eye_n) + np.kron(eye_n, J_s.T)

    ev_C = np.linalg.eigvals(L_C)

    # Check palindromic pairing in correlation space
    rates_C = sorted(-ev.real for ev in ev_C if abs(ev.real) > 1e-10)
    if len(rates_C) >= 2:
        center = (min(rates_C) + max(rates_C)) / 2.0
        tol = 0.03 * (max(rates_C) - min(rates_C))
        paired = 0
        used = [False] * len(rates_C)
        for i in range(len(rates_C)):
            if used[i]:
                continue
            target = 2 * center - rates_C[i]
            for j in range(len(rates_C)):
                if i != j and not used[j] and abs(rates_C[j] - target) < tol:
                    used[i] = True
                    used[j] = True
                    paired += 2
                    break
        pct = paired / len(rates_C) * 100

        # Energy ratio in correlation space
        E_freq_C = np.sum(np.abs(np.imag(ev_C)))
        E_decay_C = np.sum(np.abs(np.real(ev_C)))
        osc_C = E_freq_C / (E_freq_C + E_decay_C) if (E_freq_C + E_decay_C) > 0 else 0

        print(f"  alpha={alpha:.1f}:  {paired}/{len(rates_C)} paired ({pct:.0f}%)  "
              f"osc_frac={osc_C:.4f}  "
              f"n_freq={len(set(np.round(np.abs(np.imag(ev_C)), 4)))}")
    else:
        print(f"  alpha={alpha:.1f}:  insufficient rates")

# Build Q in correlation space: Q_C = Q x Q
print(f"\n  Algebraic test in correlation space:")
for alpha in [0.3, 1.0]:
    J_s = build_linear_jacobian(W_s, signs_s, tau_E, tau_I, alpha)
    eye_n = np.eye(N_small)
    L_C = np.kron(J_s, eye_n) + np.kron(eye_n, J_s.T)

    # Build Q for activity space
    e_local = list(np.where(signs_s > 0)[0])
    i_local = list(np.where(signs_s < 0)[0])
    perm = np.arange(N_small)
    for k in range(min(len(e_local), len(i_local))):
        perm[e_local[k]] = i_local[k]
        perm[i_local[k]] = e_local[k]
    Q = np.zeros((N_small, N_small))
    for i in range(N_small):
        Q[i, perm[i]] = 1.0

    # Q in correlation space: Q_C = Q x Q
    Q_C = np.kron(Q, Q)

    # Palindrome test in correlation space
    QJQ_C = Q_C @ L_C @ Q_C.T
    S_diag_C = -(np.diag(QJQ_C) + np.diag(L_C)) / 2.0
    R_C = QJQ_C + L_C + 2 * np.diag(S_diag_C)
    R_C_off = R_C - np.diag(np.diag(R_C))
    norm_LC = np.linalg.norm(L_C)
    residual_C = np.linalg.norm(R_C_off) / norm_LC if norm_LC > 0 else 0

    # Compare with activity space
    residual_act = 1 - palindrome_quality(J_s, signs_s)

    print(f"  alpha={alpha:.1f}:  activity residual = {residual_act:.4f}  "
          f"correlation residual = {residual_C:.4f}  "
          f"ratio = {residual_C/residual_act:.3f}" if residual_act > 0
          else f"  alpha={alpha:.1f}:  activity residual = 0")
