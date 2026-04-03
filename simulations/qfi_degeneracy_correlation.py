"""
QFI-based degeneracy-geometry correlation (fixes Bures sqrtm issues)
=====================================================================
Replaces the sqrtm-based Bures speed with the Quantum Fisher Information
metric, which only needs eigendecomposition of rho (stable even for
near-singular matrices). Uses logarithmic time grid and analytical drho/dt.

Output: simulations/results/qfi_degeneracy_correlation.txt
"""

import numpy as np
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
GAMMA = 0.05
J = 1.0
GRID = 2 * GAMMA
EPS = 1e-14

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
    d = 2**N
    Id = np.eye(d, dtype=complex)
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N
            ops[i] = P
            ops[i + 1] = P
            H += J * kron_chain(ops)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        ops = [I2] * N
        ops[k] = Zm
        Lk = np.sqrt(GAMMA) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T))
    return L


# ─────────────────────────────────────────────
# Initial states
# ─────────────────────────────────────────────

def state_plus(N):
    """ρ = |+⟩⟨+|^⊗N"""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    rho = np.outer(plus, plus.conj())
    for _ in range(N - 1):
        rho = np.kron(rho, np.outer(plus, plus.conj()))
    return rho


def state_ghz(N):
    """ρ = |GHZ⟩⟨GHZ|, GHZ = (|00..0⟩ + |11..1⟩)/√2"""
    d = 2**N
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[d - 1] = 1 / np.sqrt(2)
    return np.outer(psi, psi.conj())


def state_w(N):
    """ρ = |W⟩⟨W|, W = (|00..01⟩ + |00..10⟩ + ... + |10..00⟩)/√N"""
    d = 2**N
    psi = np.zeros(d, dtype=complex)
    for k in range(N):
        psi[1 << k] = 1 / np.sqrt(N)
    return np.outer(psi, psi.conj())


def state_random(N, seed=42):
    """Random pure state."""
    rng = np.random.RandomState(seed)
    d = 2**N
    psi = rng.randn(d) + 1j * rng.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


# ─────────────────────────────────────────────
# QFI speed
# ─────────────────────────────────────────────

def qfi_speed(rho, drho_dt):
    """Quantum Fisher Information speed: √(Σ_{ij} (p_i-p_j)²/(p_i+p_j) |⟨i|ρ̇|j⟩|²)."""
    eigvals, eigvecs = np.linalg.eigh(rho)
    d = len(eigvals)

    # drho in eigenbasis
    drho_eig = eigvecs.conj().T @ drho_dt @ eigvecs

    speed_sq = 0.0
    for i in range(d):
        for j in range(d):
            denom = eigvals[i] + eigvals[j]
            if denom > EPS:
                speed_sq += 2.0 * (eigvals[i] - eigvals[j])**2 / denom * abs(drho_eig[i, j])**2
    return np.sqrt(max(speed_sq, 0.0))


# ─────────────────────────────────────────────
# Time grid
# ─────────────────────────────────────────────

def make_time_grid():
    t_early = np.logspace(-2, 0, 60)       # 0.01 to 1.0
    t_mid = np.linspace(1.0, 20.0, 60)[1:]  # 1.0 to 20
    t_late = np.linspace(20.0, 200.0, 30)[1:]  # 20 to 200
    return np.concatenate([t_early, t_mid, t_late])


def assign_grid(re_val, N):
    k = round(-re_val / GRID)
    return k if 0 <= k <= N else -1


# ─────────────────────────────────────────────
# Core computation
# ─────────────────────────────────────────────

def compute_qfi_trajectory(N, rho0, L_mat, eigvals, R, L_left, coeffs, grid_assign, times):
    """Compute QFI speed and weight profile along trajectory."""
    d = 2**N

    speeds = []
    weight_profiles = []
    dominant_ks = []

    for t in times:
        exp_l = np.exp(eigvals * t)

        # ρ(t) from eigendecomposition
        rho_vec = R @ (coeffs * exp_l)
        rho = rho_vec.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        tr = np.trace(rho).real
        if tr > EPS:
            rho /= tr

        # dρ/dt = L(ρ) = Σ λ_j c_j exp(λ_j t) v_j (analytical)
        drho_vec = R @ (coeffs * eigvals * exp_l)
        drho = drho_vec.reshape(d, d)
        drho = (drho + drho.conj().T) / 2
        if tr > EPS:
            drho /= tr

        # QFI speed
        v = qfi_speed(rho, drho)
        speeds.append(v)

        # Weight profile
        wp = np.zeros(N + 1)
        for k in range(N + 1):
            mask = grid_assign == k
            wp[k] = np.sum(np.abs(coeffs[mask] * exp_l[mask])**2)
        s = wp.sum()
        if s > 0:
            wp /= s
        weight_profiles.append(wp)
        dominant_ks.append(int(np.argmax(wp)))

    return np.array(speeds), np.array(weight_profiles), np.array(dominant_ks)


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

out = []
def log(msg=""):
    print(msg)
    out.append(msg)


log("=" * 75)
log("QFI-BASED DEGENERACY-GEOMETRY CORRELATION")
log("=" * 75)
log(f"Parameters: J={J}, gamma={GAMMA}, Chain topology")
log()

times = make_time_grid()
log(f"Time grid: {len(times)} points, t ∈ [{times[0]:.3f}, {times[-1]:.0f}]")
log()

# Load d_total for each N
d_totals = {}
for N in range(2, 6):
    path = RESULTS_DIR / f"rmt_eigenvalues_N{N}.csv"
    if path.exists():
        ed = np.loadtxt(path, delimiter="\t", skiprows=1)
        ec = ed[:, 0] + 1j * ed[:, 1]
        dt = []
        dr = []
        for k in range(N + 1):
            on = np.abs(ec.real + k * GRID) < 1e-8
            dt.append(int(np.sum(on)))
            dr.append(int(np.sum(on & (np.abs(ec.imag) < 1e-8))))
        d_totals[N] = (dt, dr)

# ─── N=2..5 with |+⟩ initial state ───

correlation_results = {}

for N in range(2, 6):
    d = 2**N
    d2 = d * d

    log("=" * 75)
    log(f"N={N}: d={d}, initial state |+⟩^{N}")
    log("=" * 75)

    L_mat = build_liouvillian(N)
    eigvals, R = np.linalg.eig(L_mat)
    _, Lf = np.linalg.eig(L_mat.T)
    ov = Lf.conj().T @ R
    for j in range(d2):
        Lf[:, j] /= ov[j, j]

    rho0 = state_plus(N)
    coeffs = Lf.conj().T @ rho0.ravel()
    grid_assign = np.array([assign_grid(eigvals[j].real, N) for j in range(d2)])

    speeds, wps, dom_ks = compute_qfi_trajectory(N, rho0, L_mat, eigvals, R, Lf, coeffs, grid_assign, times)

    # Sampled output
    log(f"\n  QFI speed profile (sampled):")
    log(f"  {'t':>8s} {'v_QFI':>12s} {'dom_k':>6s}")
    sample_idx = np.linspace(0, len(times) - 1, 20, dtype=int)
    for si in sample_idx:
        log(f"  {times[si]:8.3f} {speeds[si]:12.6f} {dom_ks[si]:6d}")

    # Peak speed
    peak_i = np.argmax(speeds)
    log(f"\n  Peak QFI speed: {speeds[peak_i]:.6f} at t = {times[peak_i]:.3f}")
    log(f"  Weight profile at peak: [{' '.join(f'{w:.3f}' for w in wps[peak_i])}]")

    # Max speed per dominant weight sector
    max_v_per_k = np.zeros(N + 1)
    for k in range(N + 1):
        mask = dom_ks == k
        if mask.any():
            max_v_per_k[k] = np.max(speeds[mask])

    log(f"\n  Max QFI speed by dominant weight sector:")
    dt, dr = d_totals.get(N, ([], []))
    log(f"  {'k':>3s} {'d_real':>7s} {'d_total':>8s} {'max_v_QFI':>12s}")
    for k in range(N + 1):
        d_r = dr[k] if k < len(dr) else 0
        d_t = dt[k] if k < len(dt) else 0
        log(f"  {k:3d} {d_r:7d} {d_t:8d} {max_v_per_k[k]:12.6f}")

    # Pearson correlation
    if dt:
        d_arr = np.array(dt, dtype=float)
        v_arr = max_v_per_k
        if np.std(d_arr) > 0 and np.std(v_arr) > 0:
            corr = np.corrcoef(d_arr, v_arr)[0, 1]
        else:
            corr = 0.0
        log(f"\n  Pearson r(d_total, max_v_QFI): {corr:.4f}")
        correlation_results[N] = corr

    log()

# ─── N=4 with multiple initial states ───

log("=" * 75)
log("N=4: INITIAL STATE COMPARISON")
log("=" * 75)

N = 4
d = 2**N
d2 = d * d
L_mat = build_liouvillian(N)
eigvals, R = np.linalg.eig(L_mat)
_, Lf = np.linalg.eig(L_mat.T)
ov = Lf.conj().T @ R
for j in range(d2):
    Lf[:, j] /= ov[j, j]
grid_assign = np.array([assign_grid(eigvals[j].real, N) for j in range(d2)])

states = [
    ("|+⟩^4", state_plus(N)),
    ("GHZ", state_ghz(N)),
    ("W", state_w(N)),
    ("Random", state_random(N)),
]

dt_n4, dr_n4 = d_totals.get(N, ([], []))

log(f"\n  d_total(N=4) = {dt_n4}\n")

for name, rho0 in states:
    coeffs = Lf.conj().T @ rho0.ravel()
    speeds, wps, dom_ks = compute_qfi_trajectory(N, rho0, L_mat, eigvals, R, Lf, coeffs, grid_assign, times)

    max_v_per_k = np.zeros(N + 1)
    for k in range(N + 1):
        mask = dom_ks == k
        if mask.any():
            max_v_per_k[k] = np.max(speeds[mask])

    peak_i = np.argmax(speeds)

    if dt_n4:
        d_arr = np.array(dt_n4, dtype=float)
        v_arr = max_v_per_k
        corr = np.corrcoef(d_arr, v_arr)[0, 1] if np.std(d_arr) > 0 and np.std(v_arr) > 0 else 0.0
    else:
        corr = 0.0

    log(f"  {name:10s}: peak v_QFI = {speeds[peak_i]:.4f} at t={times[peak_i]:.3f}  "
        f"max_v = [{', '.join(f'{v:.4f}' for v in max_v_per_k)}]  r={corr:.4f}")

log()

# ─── Summary ───

log("=" * 75)
log("SUMMARY")
log("=" * 75)
log()
log("Pearson correlation d_total(k) vs max_v_QFI(k), initial state |+⟩^N:")
for N in sorted(correlation_results):
    log(f"  N={N}: r = {correlation_results[N]:.4f}")
log()

all_strong = all(abs(r) > 0.7 for r in correlation_results.values())
all_positive = all(r > 0 for r in correlation_results.values())
log(f"All correlations > 0.7: {'✓' if all_strong else '✗'}")
log(f"All correlations positive: {'✓' if all_positive else '✗'}")

if all_strong and all_positive:
    log("RESULT: Degeneracy shapes geometry. Higher d(k) → faster QFI speed.")
elif all_positive:
    log("RESULT: Positive trend but weak at some N. Partially confirmed.")
else:
    log("RESULT: No consistent correlation. Degeneracy does not directly shape geometry.")

# Save
out_path = RESULTS_DIR / "qfi_degeneracy_correlation.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
print(f"\n>>> Results saved to: {out_path}")
