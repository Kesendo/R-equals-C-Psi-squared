"""F86 bare-doubled-PTF closed-form derivation of the universal HWHM_left/x_peak
floor (= BareDoubledPtfHwhmRatio = 0.671535 in C2HwhmRatio.cs).

Setting: γ₀ = 1, t = t_peak = 1/4.

Closed forms derived analytically (verified bit-exact against the L_2 Duhamel
numerical evaluator):

POST-EP (x > 1, ξ = √(x²−1)):
    K_b(ξ) = e^(−2) · x · [(ξ² + 2)·cos(ξ) − 2] / ξ⁴       with x = √(ξ² + 1)

PRE-EP (x < 1, μ = √(1−x²)):
    K_b(μ) = e^(−2) · x · [(2 − μ²)·cosh(μ) − 2] / μ⁴      with x = √(1 − μ²)

EP limit (x = 1, ξ = μ = 0):  K_b(EP) = −5·e^(−2)/12

The two regimes connect via analytic continuation ξ = iμ (so cos(ξ) = cosh(μ),
ξ² = −μ², etc.). The post-EP form is the canonical one for x_peak; the pre-EP
form is needed for the HWHM_left search because x_half (≈ 0.722) lies in the
pre-EP regime.

Result:
- x_peak = √(ξ_peak² + 1) where ξ_peak ≈ 1.956122 solves an implicit equation
  from dK/dξ = 0 (transcendental, no clean algebraic closed form found).
- HWHM_left/x_peak = 0.671535355... (matches C2HwhmRatio.BareDoubledPtfHwhmRatio
  to all 6 surfaced decimal places).

So the bare floor is now derived analytically as an implicit-equation pair plus
a (transcendental) numerical evaluation. The CONSTANT 0.671535 itself is
transcendental and has no obvious closed form among standard π/e/√n candidates
tested.
"""

from __future__ import annotations

import math
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import sympy as sp
from scipy.optimize import brentq

# --- Step 1: derive K_b post-EP via sympy, factor cleanly ---

xi = sp.symbols("xi", positive=True, real=True)
x_of_xi = sp.sqrt(xi ** 2 + 1)
lam_plus = -4 + 2 * sp.I * xi
lam_minus = -4 - 2 * sp.I * xi
v_plus = sp.Matrix([1, (xi + sp.I) / x_of_xi])
v_minus = sp.Matrix([1, (-xi + sp.I) / x_of_xi])
R = sp.Matrix.hstack(v_plus, v_minus)
Rinv = sp.simplify(R.inv())

t_peak = sp.Rational(1, 4)
phi = xi / 2
exp_lam_plus_t = sp.exp(-1) * sp.exp(sp.I * phi)
exp_lam_minus_t = sp.exp(-1) * sp.exp(-sp.I * phi)
I_pp = t_peak * exp_lam_plus_t
I_mm = t_peak * exp_lam_minus_t
I_pm = (exp_lam_minus_t - exp_lam_plus_t) / (lam_minus - lam_plus)
I_mp = (exp_lam_plus_t - exp_lam_minus_t) / (lam_plus - lam_minus)
I_mat = sp.Matrix([[I_pp, I_pm], [I_mp, I_mm]])

V_x = sp.Matrix([[0, 2 * sp.I], [2 * sp.I, 0]])
M_eig = sp.simplify(Rinv * V_x * R)
rho_0 = sp.Matrix([1, 0])
c = sp.simplify(Rinv * rho_0)
drho_dx_eig = sp.Matrix([sum(M_eig[i, j] * I_mat[i, j] * c[j] for j in range(2)) for i in range(2)])
drho_dx = sp.simplify(R * drho_dx_eig)
rho_t = sp.simplify(R * sp.Matrix([exp_lam_plus_t * c[0], exp_lam_minus_t * c[1]]))
inner = sum(sp.conjugate(rho_t[i]) * drho_dx[i] for i in range(2))
K_b_post = sp.simplify(2 * sp.re(sp.simplify(inner)))

print("=" * 70)
print("K_b(ξ) post-EP closed form (sympy result):")
print(K_b_post)

# Factor: pull out (ξ²+1) from numerator
K_b_post_factored = sp.simplify(K_b_post * xi ** 4 * x_of_xi / sp.exp(-2))
print(f"\nK_b · ξ⁴ · x / e^(−2) = {sp.expand(K_b_post_factored)}")
# Expected: (ξ²+2)·cos(ξ) − 2 expanded = ξ²·cos(ξ) + 2·cos(ξ) − 2


# --- Step 2: stable numerical K_b in both regimes ---

def k_b_post_ep(xi_val: float) -> float:
    """Post-EP closed form, stable for ξ ≥ 0.05."""
    if xi_val < 1e-8:
        return -5.0 * math.exp(-2) / 12.0  # EP limit
    if xi_val < 0.05:
        # Taylor expansion near 0 to avoid cancellation
        # (ξ²+2)cos(ξ) − 2 = −5ξ⁴/12 + 11ξ⁶/360 + O(ξ⁸)
        num = -5 * xi_val ** 4 / 12 + 11 * xi_val ** 6 / 360
    else:
        num = (xi_val ** 2 + 2) * math.cos(xi_val) - 2
    x_val = math.sqrt(xi_val ** 2 + 1)
    return math.exp(-2) * x_val * num / xi_val ** 4


def k_b_pre_ep(mu_val: float) -> float:
    """Pre-EP closed form, mu = √(1−x²), x = √(1−μ²), valid for 0 ≤ μ < 1."""
    if mu_val < 1e-8:
        return -5.0 * math.exp(-2) / 12.0
    if mu_val < 0.05:
        num = -5 * mu_val ** 4 / 12 - 11 * mu_val ** 6 / 360  # sign flips for sinh terms
    else:
        num = (2 - mu_val ** 2) * math.cosh(mu_val) - 2
    x_val = math.sqrt(1 - mu_val ** 2)
    return math.exp(-2) * x_val * num / mu_val ** 4


def k_b(x_val: float) -> float:
    """K_b at dimensionless x ∈ (0, ∞), via post-EP or pre-EP closed form."""
    if x_val > 1.0:
        xi_val = math.sqrt(x_val ** 2 - 1)
        return k_b_post_ep(xi_val)
    elif x_val < 1.0:
        mu_val = math.sqrt(1 - x_val ** 2)
        return k_b_pre_ep(mu_val)
    else:
        return -5.0 * math.exp(-2) / 12.0


# --- Step 3: verify against Python script's brute K_b ---

print("\n--- Cross-check vs original brute K_b at sample x values ---")
sample_x = [0.5, 0.722, 1.0, 1.5, 2.196910, 3.0, 5.0]
for xv in sample_x:
    k_closed = k_b(xv)
    # Brute eval: directly diagonalize L(x) and integrate Duhamel
    L_num = np.array([[-2.0, 2j * xv], [2j * xv, -6.0]], dtype=complex)
    Vb_num = np.array([[0.0, 2j], [2j, 0.0]], dtype=complex)
    evals, R_num = np.linalg.eig(L_num)
    Rinv_num = np.linalg.inv(R_num)
    rho0_num = np.array([1.0, 0.0], dtype=complex)
    t_num = 0.25
    expL = R_num @ np.diag(np.exp(evals * t_num)) @ Rinv_num
    rho_t_num = expL @ rho0_num
    # Duhamel
    n = 2
    I_brute = np.zeros((n, n), dtype=complex)
    expLam = np.exp(evals * t_num)
    for i in range(n):
        for j in range(n):
            diff = evals[j] - evals[i]
            if abs(diff) < 1e-10:
                I_brute[i, j] = t_num * expLam[i]
            else:
                I_brute[i, j] = (expLam[j] - expLam[i]) / diff
    Vb_eig = Rinv_num @ Vb_num @ R_num
    c0 = Rinv_num @ rho0_num
    drho_eig = np.zeros(n, dtype=complex)
    for i in range(n):
        for j in range(n):
            drho_eig[i] += Vb_eig[i, j] * I_brute[i, j] * c0[j]
    drho_brute = R_num @ drho_eig
    k_brute = 2.0 * np.real(np.vdot(rho_t_num, drho_brute))
    print(f"  x = {xv:.6f}: closed = {k_closed:+.6e}, brute = {k_brute:+.6e}, diff = {abs(k_closed - k_brute):.2e}")


# --- Step 4: find x_peak via dK_b/dξ = 0 (post-EP) ---

print("\n--- x_peak via dK_b/dξ = 0 (implicit equation, transcendental) ---")
dKdxi_post = sp.diff(K_b_post, xi)
dKdxi_func = sp.lambdify(xi, dKdxi_post, modules=["numpy"])

xi_peak = brentq(lambda xv: float(dKdxi_func(xv)), 1.5, 2.5, xtol=1e-14)
x_peak = math.sqrt(xi_peak ** 2 + 1)
K_at_peak = k_b(x_peak)
print(f"  ξ_peak                  = {xi_peak:.12f}")
print(f"  x_peak = √(ξ_peak²+1)   = {x_peak:.12f}")
print(f"  |K_b(x_peak)|           = {abs(K_at_peak):.12f}")
print(f"  C# constant             = 2.196910")


# --- Step 5: find x_half via |K_b(x_half)| = |K_b(x_peak)|/2 ---

print("\n--- x_half via |K_b(x_half)| = |K_b(x_peak)|/2 in (0, x_peak) ---")
half = abs(K_at_peak) / 2.0

# Search downward from x_peak: |K_b| is 0.0851 at peak, then decreases to half
# somewhere in (0, x_peak). Per original Python script x_half ≈ 0.722 (pre-EP).
def f_minus_half(xv: float) -> float:
    return abs(k_b(xv)) - half

# Sample to find the crossing closest to x_peak going leftward
x_grid = np.linspace(0.05, x_peak - 1e-4, 2000)
crossings = []
for k in range(len(x_grid) - 1):
    a, b = x_grid[k], x_grid[k + 1]
    if f_minus_half(a) * f_minus_half(b) < 0:
        crossings.append((a, b))

print(f"  found {len(crossings)} half-max crossings in (0, x_peak)")
for (a, b) in crossings:
    root = brentq(f_minus_half, a, b, xtol=1e-14)
    print(f"    crossing near x = {root:.6f}  (regime: {'pre-EP' if root < 1 else 'post-EP'})")

# HWHM_left is from the rightmost crossing (closest to x_peak)
if crossings:
    a, b = crossings[-1]
    x_half = brentq(f_minus_half, a, b, xtol=1e-14)
    hwhm_left = x_peak - x_half
    ratio = hwhm_left / x_peak
    print(f"\n  x_half (closest to peak) = {x_half:.12f}")
    print(f"  HWHM_left                = {hwhm_left:.12f}")
    print(f"  HWHM_left / x_peak       = {ratio:.12f}")
    print(f"  C# constant              = 0.671535")
    print(f"  Residual                 = {abs(ratio - 0.671535):.2e}")
else:
    print("  no crossings found")


# --- Step 6: try algebraic identification of the ratio ---

print("\n--- Algebraic candidates for ratio = 0.671535355... ---")
target = ratio
candidates = [
    ("(π − 1)/π",            float((sp.pi - 1) / sp.pi)),
    ("1 − 1/(π + 1/π)",       float(1 - 1 / (sp.pi + 1 / sp.pi))),
    ("e⁻¹ + (1 − e⁻¹)/2",     float(sp.exp(-1) + (1 - sp.exp(-1)) / 2)),
    ("ξ_peak/x_peak (= sin θ)", xi_peak / x_peak),
    ("2/√(2π/e)",            float(2 / sp.sqrt(2 * sp.pi / sp.E))),
    ("√(1 − e⁻²)",           float(sp.sqrt(1 - sp.exp(-2)))),
    ("erf(1)",               float(sp.erf(1))),
]
for name, val in candidates:
    print(f"  {name:30s} = {val:.10f},  diff = {val - target:+.2e}")
