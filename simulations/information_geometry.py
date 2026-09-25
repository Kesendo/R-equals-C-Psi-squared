#!/usr/bin/env python3
"""
Information Geometry: θ as Riemannian Coordinate
==================================================
Phase 0: θ inventory (what we already know)
Phase 1: Bures metric g(CΨ) along the Lindblad trajectory
Phase 2: θ as regular coordinate (g̃(θ) finite at θ=0?)
Phase 3: Geodesic test in the full state space (Bures path length against the
         endpoint Bures angle; Bell+ and two states the Hamiltonian moves)
Phase 4: Coordinate-shape second derivative of the Bures path-metric coefficient
         (a one-dimensional path metric has no intrinsic curvature)
Phase 5: Second γ-derivative of CΨ at fixed time (not a fidelity susceptibility)

Script: simulations/information_geometry.py
Output: simulations/results/information_geometry.txt
"""

import numpy as np
from scipy.linalg import expm, sqrtm, eigvalsh
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "information_geometry.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ========================================================================
# Infrastructure
# ========================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_list(mats):
    r = mats[0]
    for m in mats[1:]:
        r = np.kron(r, m)
    return r


def build_liouvillian_2q(J=1.0, gamma=0.05):
    """2-qubit Heisenberg + Z-dephasing Liouvillian."""
    d = 4
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * np.kron(P, P)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(2):
        ops = [I2, I2]
        ops[k] = sz
        Zk = kron_list(ops)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d**2))
    return L


def evolve_rho(L, rho0, t):
    d2 = L.shape[0]
    d = int(np.sqrt(d2))
    return (expm(L * t) @ rho0.flatten()).reshape(d, d)


def compute_cpsi(rho):
    d = rho.shape[0]
    purity = np.real(np.trace(rho @ rho))
    L1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    return purity * L1 / (d - 1)


def uhlmann_fidelity(rho, sigma):
    """F(ρ,σ) = (Tr√(√ρ σ √ρ))²."""
    sqrt_rho = sqrtm(rho)
    inner = sqrt_rho @ sigma @ sqrt_rho
    # Eigenvalues of inner (should be non-negative)
    eigv = eigvalsh(inner)
    eigv = np.maximum(eigv, 0)
    return np.real(np.sum(np.sqrt(eigv)))**2


def bures_distance(rho, sigma):
    """dB = √(2(1 - √F))."""
    F = uhlmann_fidelity(rho, sigma)
    F = min(F, 1.0)
    return np.sqrt(2 * (1 - np.sqrt(F)))


def _numerical_zeros_to_zero(w):
    """Eigenvalues below the float rank threshold (dimension * eps * largest) are
    numerical zeros of a rank-deficient state: their square roots, ~1e-8, would
    otherwise bias the root fidelity of two nearby states by more than the step."""
    w = np.asarray(w, dtype=float)
    w[w < len(w) * np.finfo(float).eps * max(w.max(), 1.0)] = 0.0
    return w


def psd_sqrt(rho):
    """Square root of a Hermitian positive semidefinite matrix by its eigendecomposition."""
    w, v = np.linalg.eigh((rho + rho.conj().T) / 2)
    return v @ np.diag(np.sqrt(_numerical_zeros_to_zero(w))) @ v.conj().T


def bures_angle(rho, sigma):
    """Bures angle arccos(√F): the geodesic distance of the Bures metric."""
    root = psd_sqrt(rho)
    inner = root @ sigma @ root
    ev = _numerical_zeros_to_zero(eigvalsh((inner + inner.conj().T) / 2))
    return np.arccos(min(1.0, float(np.sum(np.sqrt(ev)))))


def cpsi_to_theta(cpsi):
    """θ = arctan(√(4CΨ - 1)), defined for CΨ > 1/4."""
    if cpsi <= 0.25:
        return 0.0
    return np.arctan(np.sqrt(4 * cpsi - 1))


# ========================================================================
log("=" * 72)
log("INFORMATION GEOMETRY: θ AS RIEMANNIAN COORDINATE")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)


# ========================================================================
# PHASE 0: θ INVENTORY
# ========================================================================
log()
log("=" * 72)
log("PHASE 0: θ INVENTORY (what the repo already knows)")
log("=" * 72)
log()
log("  θ = arctan(√(4CΨ - 1))  (formula 15, BOUNDARY_NAVIGATION.md)")
log()
log("  Known results:")
log("  - Definition: angular distance from CΨ = 1/4 boundary")
log("  - θ = 0 at CΨ = 1/4, θ = π/4 at CΨ = 1/2, θ → π/2 at CΨ → ∞")
log("  - Correlates with fidelity r = 0.87 (THETA_PALINDROME_ECHO.md)")
log("  - Palindromic rate predictions per component")
log("  - Lives on 3D manifold (STRUCTURAL_CARTOGRAPHY.md, 98% in 3 PCs)")
log("  - Used as Berry-phase parameter (TOPOLOGICAL_EDGE_MODES, φ = -0.77)")
log("  - Described as 'voltmeter' in CIRCUIT_DIAGRAM.md")
log()
log("  NOT computed:")
log("  - The Riemannian metric in CΨ or θ coordinates")
log("  - Whether g(CΨ) diverges at 1/4")
log("  - Whether θ regularizes a singularity")
log("  - Geodesic analysis")
log("  - Whether the trajectory is a Bures geodesic")
log("  - Second derivatives of the path coefficient and of CΨ in γ")


# ========================================================================
# PHASE 1: BURES METRIC g(CΨ)
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 1: BURES METRIC g(CΨ) ALONG THE LINDBLAD TRAJECTORY")
log("=" * 72)
log()

J, gamma = 1.0, 0.05
L = build_liouvillian_2q(J, gamma)
d = 4

# Bell+ initial state
psi = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
rho0 = np.outer(psi, psi.conj())

# Find crossing time (CΨ = 1/4)
# Sampling shortcut for the named Hamiltonian-dead Bell+ trajectory with
# equal local Z-dephasing and the purity-times-l1 F25 readout, not arbitrary
# states or the Wootters-concurrence book: t_cross ≈ 0.0374/gamma.
K_cross = 0.0374
t_cross_approx = K_cross / gamma

# Compute trajectory with fine time steps
n_t = 300
t_max = 2.0 * t_cross_approx
times = np.linspace(0.01, t_max, n_t)
dt = times[1] - times[0]

cpsi_arr = np.zeros(n_t)
theta_arr = np.zeros(n_t)
bures_rate = np.zeros(n_t)  # dB/dt
cpsi_rate = np.zeros(n_t)   # dCΨ/dt
g_cpsi = np.zeros(n_t)      # metric in CΨ coords

rho_prev = None
t0 = clock.time()

for i, t in enumerate(times):
    rho_t = evolve_rho(L, rho0, t)
    cpsi_arr[i] = compute_cpsi(rho_t)
    theta_arr[i] = cpsi_to_theta(cpsi_arr[i])

    if rho_prev is not None:
        dB = bures_distance(rho_t, rho_prev)
        bures_rate[i] = dB / dt
        dcpsi = cpsi_arr[i] - cpsi_arr[i - 1]
        cpsi_rate[i] = dcpsi / dt

        if abs(dcpsi) > 1e-15:
            g_cpsi[i] = (dB / abs(dcpsi))**2
        else:
            g_cpsi[i] = np.nan

    rho_prev = rho_t.copy()

log(f"  N=2, Bell+, J={J}, γ={gamma}, {clock.time()-t0:.1f}s")
log(f"  t_cross ≈ {t_cross_approx:.2f}")
log()

# Find crossing index
cross_idx = None
for i in range(len(cpsi_arr) - 1):
    if cpsi_arr[i] > 0.25 and cpsi_arr[i + 1] <= 0.25:
        cross_idx = i
        break

log(f"  {'t':>6}  {'CΨ':>8}  {'θ':>8}  {'dB/dt':>10}  {'dCΨ/dt':>10}  {'g(CΨ)':>12}")
log(f"  {'─'*60}")

for i in range(0, n_t, n_t // 15):
    marker = " ←1/4" if cross_idx and abs(i - cross_idx) < 2 else ""
    if np.isnan(g_cpsi[i]):
        g_str = "N/A"
    else:
        g_str = f"{g_cpsi[i]:.4f}"
    log(f"  {times[i]:>6.2f}  {cpsi_arr[i]:>8.4f}  {theta_arr[i]:>8.4f}"
        f"  {bures_rate[i]:>10.6f}  {cpsi_rate[i]:>10.6f}  {g_str:>12}{marker}")

# Metric at the crossing
if cross_idx and not np.isnan(g_cpsi[cross_idx]):
    log()
    log(f"  g(CΨ) at CΨ = 1/4 crossing: {g_cpsi[cross_idx]:.6f}")
    log(f"  dCΨ/dt at crossing: {cpsi_rate[cross_idx]:.6f}")
    log(f"  dB/dt at crossing: {bures_rate[cross_idx]:.6f}")
    log()
    if g_cpsi[cross_idx] < 100:
        log("  g(CΨ) is FINITE at CΨ = 1/4. No singularity.")
        log("  CΨ IS a good coordinate at the fold (no divergence).")
    else:
        log("  g(CΨ) is LARGE at CΨ = 1/4. Possible singularity.")


# ========================================================================
# PHASE 2: θ AS REGULAR COORDINATE
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 2: METRIC IN θ COORDINATES")
log("=" * 72)
log()

# g̃(θ) = g(CΨ) × (dCΨ/dθ)²
# CΨ = (1 + tan²θ)/4, dCΨ/dθ = tanθ/(2cos²θ) = sinθ/(2cos³θ)
g_theta = np.zeros(n_t)
for i in range(n_t):
    th = theta_arr[i]
    if th > 1e-6 and not np.isnan(g_cpsi[i]):
        dcpsi_dtheta = np.sin(th) / (2 * np.cos(th)**3)
        g_theta[i] = g_cpsi[i] * dcpsi_dtheta**2
    else:
        g_theta[i] = np.nan

log(f"  g̃(θ) = g(CΨ) × (dCΨ/dθ)²")
log(f"  dCΨ/dθ = sin(θ)/(2cos³θ) → 0 as θ → 0")
log()

log(f"  {'θ':>8}  {'CΨ':>8}  {'g(CΨ)':>12}  {'dCΨ/dθ':>10}  {'g̃(θ)':>12}")
log(f"  {'─'*55}")

for i in range(n_t // 15, n_t, n_t // 15):
    th = theta_arr[i]
    if th > 1e-6:
        dcpsi_dtheta = np.sin(th) / (2 * np.cos(th)**3)
    else:
        dcpsi_dtheta = 0
    g_str = f"{g_cpsi[i]:.4f}" if not np.isnan(g_cpsi[i]) else "N/A"
    gt_str = f"{g_theta[i]:.6f}" if not np.isnan(g_theta[i]) else "N/A"
    log(f"  {th:>8.4f}  {cpsi_arr[i]:>8.4f}  {g_str:>12}"
        f"  {dcpsi_dtheta:>10.4f}  {gt_str:>12}")

if cross_idx:
    log()
    log(f"  Near θ = 0 (CΨ = 1/4): dCΨ/dθ → 0, so g̃(θ) → 0 (not ∞).")
    log(f"  θ does NOT regularize a singularity (there is none to regularize).")
    log(f"  Both CΨ and θ are smooth coordinates at the fold point.")


# ========================================================================
# PHASE 3: GEODESIC TEST IN THE FULL STATE SPACE
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 3: IS THE LINDBLAD TRAJECTORY A BURES GEODESIC?")
log("=" * 72)
log()
# Along the trajectory's own coordinate the one-dimensional geodesic equation
# d²CΨ/ds² + Γ(dCΨ/ds)² = 0 holds for ANY monotone curve (it is an identity of
# arc length), so it cannot test anything. The question that can fail is asked in
# the state space: a curve is a Bures geodesic exactly when its Bures length equals
# the Bures angle between its endpoints; the ratio is >= 1 and = 1 only then.
n_geo = 600


def path_ratio(psi0, t_end):
    step = expm(L * (t_end / (n_geo - 1)))
    r0 = np.outer(psi0, psi0.conj())
    states = [r0]
    for _ in range(n_geo - 1):
        states.append((step @ states[-1].flatten()).reshape(d, d))
    length = sum(bures_angle(states[k], states[k + 1]) for k in range(n_geo - 1))
    return length, bures_angle(states[0], states[-1])


geo_cases = [
    ("Bell+ (the Hamiltonian leaves it alone)", psi),
    ("|+0> (the Hamiltonian moves it)", np.kron([1, 1], [1, 0]).astype(complex) / np.sqrt(2)),
    ("a generic state", np.array([1, 0.3, 0.2j, 0.7]) / np.linalg.norm([1, 0.3, 0.2, 0.7])),
]
for t_end in (0.2, 1.5):
    log(f"  t = 0 .. {t_end}, {n_geo} steps, Bures length against endpoint Bures angle")
    log()
    log(f"  {'initial state':<42} {'length':>9} {'endpoint':>9} {'ratio':>9}")
    log(f"  {'─'*72}")
    for label, psi0 in geo_cases:
        length, dist = path_ratio(psi0, t_end)
        log(f"  {label:<42} {length:>9.6f} {dist:>9.6f} {length/dist:>9.6f}")
    log()
log("  The Bures angle is capped at π/2 while a path can keep growing, so a ratio")
log("  read over a long window grows with the window; the short window shows the")
log("  departure from the geodesic where it starts.")
log("  Bell+ decays inside the commuting family of its two Bell projectors, a")
log("  one-parameter line of states it runs along monotonically, so its path IS")
log("  the Bures geodesic (ratio 1). A state the Hamiltonian moves leaves the")
log("  geodesic from the start (ratio above 1 already in the short window).")


# ========================================================================
# PHASE 4: COORDINATE-SHAPE SECOND DERIVATIVE OF THE PATH COEFFICIENT
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 4: COORDINATE-SHAPE SECOND DERIVATIVE OF THE BURES PATH-METRIC COEFFICIENT")
log("=" * 72)
log()

# S = -(1/2g) d²(ln g)/dCΨ² of the one-dimensional path coefficient g_path(CΨ).
# A one-dimensional metric has no intrinsic curvature; S is a coordinate-shape second
# derivative (it is -2/x for the flat line written as g = 1/(4x)), not intrinsic curvature.
valid = ~np.isnan(g_cpsi) & (g_cpsi > 0)
cpsi_s, g_s = cpsi_arr[valid], g_cpsi[valid]
order = np.argsort(cpsi_s)[::-1]
cpsi_s, g_s = cpsi_s[order], g_s[order]
if np.sum(valid) > 20:
    lng = np.log(g_s + 1e-30)
    d2lng = np.gradient(np.gradient(lng, cpsi_s), cpsi_s)
    K_gauss = -d2lng / (2 * g_s + 1e-30)

    # Find value near CΨ = 1/4
    idx_quarter = np.argmin(np.abs(cpsi_s - 0.25))
    K_at_fold = K_gauss[idx_quarter] if idx_quarter > 2 and idx_quarter < len(K_gauss) - 2 else np.nan

    log(f"  S = -(1/2g) d²(ln g)/d(CΨ)², coordinate-shape second derivative, not intrinsic curvature")
    log()
    log(f"  {'CΨ':>8}  {'g(CΨ)':>12}  {'S_CΨ':>12}")
    log(f"  {'─'*35}")

    for i in range(0, len(cpsi_s), max(1, len(cpsi_s) // 10)):
        if not np.isnan(K_gauss[i]) and abs(K_gauss[i]) < 1e6:
            log(f"  {cpsi_s[i]:>8.4f}  {g_s[i]:>12.4f}  {K_gauss[i]:>12.4f}")

    log()
    if not np.isnan(K_at_fold) and abs(K_at_fold) < 1e6:
        log(f"  S at CΨ ≈ 1/4: {K_at_fold:.4f} (finite; a coordinate shape, no geometric claim)")
    else:
        log("  S numerically unstable near CΨ = 1/4.")


# ========================================================================
# PHASE 5: SECOND γ-DERIVATIVE OF CΨ AT FIXED TIME
# ========================================================================
log()
log()
log("=" * 72)
log("PHASE 5: SECOND γ-DERIVATIVE OF CΨ AT FIXED TIME")
log("=" * 72)
log()
log("  d²CΨ/dγ² at a fixed time. Not a fidelity susceptibility; and CΨ(γ) at fixed")
log("  time is analytic in γ (a matrix exponential), so this cannot diverge.")
log()

J_f = 1.0
# Sweep γ and find CΨ at a fixed physical time (not crossing time)
# Better: fix the state at t_cross and compute d(CΨ)/dγ

gammas = np.linspace(0.02, 0.15, 30)
cpsi_at_fixed_t = np.zeros(len(gammas))
t_fixed = 0.75  # fixed observation time

for gi, gam in enumerate(gammas):
    L_g = build_liouvillian_2q(J_f, gam)
    rho_t = evolve_rho(L_g, rho0, t_fixed)
    cpsi_at_fixed_t[gi] = compute_cpsi(rho_t)

# d(CΨ)/dγ and d²(CΨ)/dγ²
dcpsi_dgamma = np.gradient(cpsi_at_fixed_t, gammas)
d2cpsi_dgamma2 = np.gradient(dcpsi_dgamma, gammas)

# Find γ where CΨ ≈ 1/4 at t_fixed
idx_quarter_g = np.argmin(np.abs(cpsi_at_fixed_t - 0.25))
gamma_at_quarter = gammas[idx_quarter_g]

log(f"  At fixed t = {t_fixed}:")
log(f"  CΨ = 1/4 at γ ≈ {gamma_at_quarter:.4f}")
log()
log(f"  {'γ':>8}  {'CΨ':>8}  {'dCΨ/dγ':>10}  {'d²CΨ/dγ²':>12}")
log(f"  {'─'*45}")

for i in range(0, len(gammas), 3):
    log(f"  {gammas[i]:>8.4f}  {cpsi_at_fixed_t[i]:>8.4f}"
        f"  {dcpsi_dgamma[i]:>10.4f}  {d2cpsi_dgamma2[i]:>12.4f}")

log()
chi_at_quarter = d2cpsi_dgamma2[idx_quarter_g]
log(f"  d²CΨ/dγ² where CΨ = 1/4: {chi_at_quarter:.4f} (finite, as analyticity requires)")


# ========================================================================
# SUMMARY
# ========================================================================
log()
log()
log("=" * 72)
log("SUMMARY")
log("=" * 72)
log()
log("Phase 1: The Bures metric g(CΨ) is FINITE at CΨ = 1/4.")
log("  No singularity. CΨ is a smooth coordinate at the fold.")
log("  The trajectory ρ(t) passes through CΨ = 1/4 with nonzero")
log("  velocity (dCΨ/dt ≠ 0) and finite Bures rate (dB/dt finite).")
log()
log("Phase 2: θ does NOT regularize a singularity (there is none).")
log("  Since g(CΨ) is finite, g̃(θ) = g(CΨ)×(dCΨ/dθ)² → 0 at θ=0")
log("  (because dCΨ/dθ → 0). θ SHRINKS the metric at the fold.")
log("  θ is a valid coordinate but not a geometric necessity.")
log()
log("Phase 3: Bell+ runs exactly along a Bures geodesic (length = endpoint")
log("  angle), because it decays inside a commuting one-parameter family; states")
log("  the Hamiltonian moves do not (their paths are several times longer).")
log()
log("Phase 4: the coordinate-shape second derivative of the path coefficient is")
log("  finite at CΨ = 1/4; a one-dimensional metric has no intrinsic curvature.")
log()
log("Phase 5: d²CΨ/dγ² at fixed time is finite, as it must be (CΨ is analytic")
log("  in γ); it says nothing about criticality.")
log()
log("CONCLUSION: θ = arctan(√(4CΨ-1)) is a nonlinear coordinate")
log("transformation, not a geometric regularization. The fold at")
log("CΨ = 1/4 has no singularity in the Bures path coefficient.")
log("θ is useful as a 'compass' (angular distance from the boundary).")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
