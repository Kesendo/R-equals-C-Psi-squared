"""
Bures metric and degeneracy structure correlation
===================================================
Propagates density matrices under Lindblad dynamics for N=2..5,
computes the Bures speed at each time step, decomposes into
eigenmode weight sectors, and correlates speed with degeneracy.

Output: simulations/results/bures_degeneracy_correlation.txt
"""

import numpy as np
from scipy.linalg import sqrtm, expm
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
GAMMA = 0.05
J = 1.0
GRID = 2 * GAMMA  # 0.1
TOL = 1e-10

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_liouvillian(N):
    """Build full Liouvillian superoperator as d²×d² matrix."""
    d = 2**N
    Id = np.eye(d, dtype=complex)

    # Hamiltonian: Heisenberg chain
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N
            ops[i] = P
            ops[i + 1] = P
            H += J * kron_chain(ops)

    # Liouvillian: L = -i(H⊗I - I⊗H^T) + dissipator
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))

    # Z-dephasing dissipator
    for k in range(N):
        ops = [I2] * N
        ops[k] = Zm
        Lk = np.sqrt(GAMMA) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T))

    return L


def initial_state_plus(N):
    """ρ = |+⟩⟨+|^⊗N as a d×d density matrix."""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    rho_1 = np.outer(plus, plus.conj())
    rho = rho_1
    for _ in range(N - 1):
        rho = np.kron(rho, rho_1)
    return rho


def bures_fidelity(rho, sigma):
    """F(ρ,σ) = (Tr √(√ρ σ √ρ))²."""
    # Regularize near-zero eigenvalues
    rho_half = sqrtm(rho)
    M = rho_half @ sigma @ rho_half
    # Eigenvalues of M (should be real non-negative)
    eigvals = np.linalg.eigvalsh((M + M.conj().T) / 2)  # force Hermitian
    eigvals = np.maximum(eigvals, 0)
    return np.sum(np.sqrt(eigvals))**2


def bures_distance(rho, sigma):
    """d_B = arccos(√F), clamped for numerical safety."""
    F = bures_fidelity(rho, sigma)
    F = np.clip(F.real, 0, 1)
    return np.arccos(np.sqrt(F))


def assign_grid(re_val, N):
    """Assign a real-part value to grid position k=0..N."""
    k = round(-re_val / GRID)
    if 0 <= k <= N:
        return k
    return -1


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

out = []
def log(msg=""):
    print(msg)
    out.append(msg)


log("=" * 75)
log("BURES METRIC AND DEGENERACY CORRELATION")
log("=" * 75)
log(f"Parameters: J={J}, gamma={GAMMA}, Chain topology, initial state |+⟩^N")
log()

for N in range(2, 6):
    d = 2**N
    d2 = d * d

    log("=" * 75)
    log(f"N={N}: d={d}")
    log("=" * 75)

    # Build Liouvillian and eigendecompose
    L = build_liouvillian(N)
    eigvals, eigvecs_right = np.linalg.eig(L)

    # Also need left eigenvectors for projection: L^T u = λ* u
    _, eigvecs_left = np.linalg.eig(L.T)
    # Normalize: ⟨u_i | v_j⟩ = δ_ij
    overlap = eigvecs_left.conj().T @ eigvecs_right
    for j in range(d2):
        eigvecs_left[:, j] /= overlap[j, j]

    # Initial state
    rho0 = initial_state_plus(N)
    rho0_vec = rho0.ravel()

    # Project onto eigenmodes: c_j = ⟨u_j | ρ0⟩
    coeffs = eigvecs_left.conj().T @ rho0_vec

    # Assign eigenvalues to grid positions
    grid_assign = np.array([assign_grid(eigvals[j].real, N) for j in range(d2)])

    # Time grid
    t_max = 10.0 / GAMMA  # = 200
    n_steps = 400
    times = np.linspace(0, t_max, n_steps + 1)
    dt = times[1] - times[0]

    # Propagate and compute Bures speed
    log(f"  Propagating t=0..{t_max:.0f} ({n_steps} steps, dt={dt:.2f})...")

    bures_speeds = []
    weight_profiles = []
    dominant_weights = []
    rho_prev = rho0.copy()

    for ti in range(n_steps + 1):
        t = times[ti]

        # ρ(t) via eigendecomposition
        exp_lambda = np.exp(eigvals * t)
        rho_vec = eigvecs_right @ (coeffs * exp_lambda)
        rho = rho_vec.reshape(d, d)

        # Force Hermitian and trace-1
        rho = (rho + rho.conj().T) / 2
        rho /= np.trace(rho).real

        # Weight profile: contribution from each grid position
        w_profile = np.zeros(N + 1)
        for k in range(N + 1):
            mask = grid_assign == k
            w_profile[k] = np.sum(np.abs(coeffs[mask] * exp_lambda[mask])**2)
        if w_profile.sum() > 0:
            w_profile /= w_profile.sum()
        weight_profiles.append(w_profile)

        # Dominant weight
        dom_k = np.argmax(w_profile)
        dominant_weights.append(dom_k)

        # Bures speed (from second step onward)
        if ti > 0:
            d_B = bures_distance(rho_prev, rho)
            v_B = d_B / dt
            bures_speeds.append(v_B)

        rho_prev = rho.copy()

    bures_speeds = np.array(bures_speeds)
    times_mid = (times[:-1] + times[1:]) / 2  # midpoints

    # ─── Bures speed profile ───
    log(f"\n  Bures speed profile (sampled):")
    log(f"  {'t':>8s} {'v_B':>10s} {'dom_k':>6s} {'w_profile':>30s}")
    sample_idx = np.linspace(0, len(bures_speeds) - 1, min(20, len(bures_speeds)), dtype=int)
    for si in sample_idx:
        t = times_mid[si]
        v = bures_speeds[si]
        dk = dominant_weights[si + 1]
        wp = weight_profiles[si + 1]
        wp_str = ' '.join(f'{w:.3f}' for w in wp)
        log(f"  {t:8.1f} {v:10.6f} {dk:6d} [{wp_str}]")

    # ─── Peak speed and its weight sector ───
    peak_idx = np.argmax(bures_speeds)
    peak_t = times_mid[peak_idx]
    peak_v = bures_speeds[peak_idx]
    peak_wp = weight_profiles[peak_idx + 1]
    peak_dk = dominant_weights[peak_idx + 1]

    log(f"\n  Peak Bures speed: v_B = {peak_v:.6f} at t = {peak_t:.1f}")
    log(f"  Dominant weight at peak: k = {peak_dk}")
    log(f"  Weight profile at peak: [{' '.join(f'{w:.3f}' for w in peak_wp)}]")

    # ─── Correlation: average Bures speed per dominant weight ───
    log(f"\n  Average Bures speed by dominant weight sector:")
    for k in range(N + 1):
        mask = np.array(dominant_weights[1:]) == k
        if mask.any():
            avg_v = np.mean(bures_speeds[mask])
            max_v = np.max(bures_speeds[mask])
            log(f"    k={k}: avg_v={avg_v:.6f}, max_v={max_v:.6f}, "
                f"time range=[{times_mid[mask][0]:.1f}, {times_mid[mask][-1]:.1f}]")

    # ─── d_total vs peak speed per weight sector ───
    # Load d_total from eigenvalue data
    eig_path = RESULTS_DIR / f"rmt_eigenvalues_N{N}.csv"
    if eig_path.exists():
        eig_data = np.loadtxt(eig_path, delimiter="\t", skiprows=1)
        eig_complex = eig_data[:, 0] + 1j * eig_data[:, 1]

        d_total = []
        d_real_seq = []
        for k in range(N + 1):
            target = -k * GRID
            on_grid = np.abs(eig_complex.real - target) < 1e-8
            d_total.append(int(np.sum(on_grid)))
            real_on_grid = on_grid & (np.abs(eig_complex.imag) < 1e-8)
            d_real_seq.append(int(np.sum(real_on_grid)))

        log(f"\n  Degeneracy vs Bures speed by weight sector:")
        log(f"  {'k':>3s} {'d_real':>7s} {'d_total':>8s} {'max_v_B':>10s}")
        for k in range(N + 1):
            mask = np.array(dominant_weights[1:]) == k
            max_v = np.max(bures_speeds[mask]) if mask.any() else 0
            log(f"  {k:3d} {d_real_seq[k]:7d} {d_total[k]:8d} {max_v:10.6f}")

    log()

# ─────────────────────────────────────────────
# Cross-N summary
# ─────────────────────────────────────────────

log("=" * 75)
log("CROSS-N CORRELATION SUMMARY")
log("=" * 75)
log()
log("Does higher degeneracy correlate with faster Bures speed?")
log("If the palindrome shapes geometry, d(k) should predict v_B.")
log()

# For each N, compute Pearson correlation between d_total(k) and max_v_B(k)
for N in range(2, 6):
    d = 2**N
    d2 = d * d

    # Rebuild quickly for correlation
    L_mat = build_liouvillian(N)
    eigvals, R = np.linalg.eig(L_mat)
    _, Lf = np.linalg.eig(L_mat.T)
    ov = Lf.conj().T @ R
    for j in range(d2):
        Lf[:, j] /= ov[j, j]

    rho0 = initial_state_plus(N)
    c = Lf.conj().T @ rho0.ravel()
    ga = np.array([assign_grid(eigvals[j].real, N) for j in range(d2)])

    eig_path = RESULTS_DIR / f"rmt_eigenvalues_N{N}.csv"
    if not eig_path.exists():
        continue
    ed = np.loadtxt(eig_path, delimiter="\t", skiprows=1)
    ec = ed[:, 0] + 1j * ed[:, 1]

    d_total = []
    for k in range(N + 1):
        on = np.abs(ec.real + k * GRID) < 1e-8
        d_total.append(int(np.sum(on)))

    # Compute max Bures speed when each weight dominates
    t_max = 10.0 / GAMMA
    times_c = np.linspace(0, t_max, 401)
    dt_c = times_c[1] - times_c[0]

    max_v_per_k = np.zeros(N + 1)
    rho_prev = initial_state_plus(N)

    for ti in range(1, len(times_c)):
        t = times_c[ti]
        exp_l = np.exp(eigvals * t)
        rv = R @ (c * exp_l)
        rho = rv.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        rho /= np.trace(rho).real

        wp = np.zeros(N + 1)
        for k in range(N + 1):
            m = ga == k
            wp[k] = np.sum(np.abs(c[m] * exp_l[m])**2)
        if wp.sum() > 0:
            wp /= wp.sum()
        dom = np.argmax(wp)

        d_B = bures_distance(rho_prev, rho)
        v_B = d_B / dt_c
        if v_B > max_v_per_k[dom]:
            max_v_per_k[dom] = v_B

        rho_prev = rho.copy()

    # Pearson correlation
    d_arr = np.array(d_total, dtype=float)
    v_arr = max_v_per_k
    if np.std(d_arr) > 0 and np.std(v_arr) > 0:
        corr = np.corrcoef(d_arr, v_arr)[0, 1]
    else:
        corr = 0

    log(f"N={N}: d_total = {d_total}")
    log(f"      max_v_B = [{', '.join(f'{v:.6f}' for v in max_v_per_k)}]")
    log(f"      Pearson correlation: r = {corr:.4f}")

    # Palindrome check for Bures speed
    v_palindrome = all(abs(max_v_per_k[k] - max_v_per_k[N - k]) < 0.001 * max(max_v_per_k)
                       for k in range(N + 1))
    log(f"      Bures speed palindromic: {'✓' if v_palindrome else '✗'}")
    log()

# Save
out_path = RESULTS_DIR / "bures_degeneracy_correlation.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
print(f"\n>>> Results saved to: {out_path}")
