"""F99: depth-3 dyadic anchor (α = 1/8) via a non-uniform Dicke
superposition at γ = √3/2 = cos(30°).

The F86b derivation parametrises
  α_total(γ) = (1 − γ²)/2   with γ = ⟨ψ|X⊗N|ψ⟩
on the X⊗N-eigenbasis decomposition of the Dicke superposition
(|D_n⟩ + |D_{n+1}⟩)/√2. The clean uniform-Dicke superposition gives only
γ ∈ {0, 1/2, 1} (Generic / KIntermediate / Mirror cases).

Tonight we generalise to NON-UNIFORM Dicke:
  ψ = (|D_n⟩ + c·|D_{n+1}⟩)/√(1+c²)
Then ⟨ψ|X⊗N|ψ⟩ = c²/(1+c²) for the canonical case n = N/2 − 1 at even N.
Rearranging: c² = γ/(1−γ). With γ = cos(θ) and the half-angle identity
1 − cos(θ) = 2sin²(θ/2):
  c² = cos(θ)/(2·sin²(θ/2))

F99 selects the five angles {0°, 30°, 45°, 60°, 90°}, which yield
{0, 1/8, 1/4, 3/8, 1/2}. This script evaluates the named states in floating
point at N = 4, 6, 8 and prints numerical residuals against the closed form.

The optional period-table display is a comparison with curated input labels
from `period_2_at_framework_anchors.py`. Equal rational values do not map an
atomic shell to this Dicke state, a Π² axis, a carbon Hamiltonian, a bath, or
an observable. A material claim would require specified degrees of freedom,
Hamiltonian/coupling, bath/dissipator, preparation, measurement, and producer
chain.

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/depth_3_anchor_derivation.py
"""
from __future__ import annotations

import sys
import numpy as np
from math import comb, sqrt
from fractions import Fraction

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ---------------------- Pauli + Π² split ----------------------

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, SX, SY, SZ]
BIT_B = [0, 0, 1, 1]  # bit_b parity for Y, Z = 1 (Z-dephasing convention)


def kron_n(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def dicke_state(N: int, k: int) -> np.ndarray:
    """Symmetric Dicke state |D_k^N⟩ = normalised sum of all binary strings
    of popcount k."""
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    norm = 1.0 / sqrt(comb(N, k))
    for b in range(d):
        if bin(b).count("1") == k:
            psi[b] = norm
    return psi


def dicke_nonuniform_superposition(N: int, n: int, c: float) -> np.ndarray:
    """ψ = (|D_n⟩ + c·|D_{n+1}⟩) / √(1+c²)."""
    psi = dicke_state(N, n) + c * dicke_state(N, n + 1)
    psi /= sqrt(1 + c * c)
    return psi


def x_global_overlap(psi: np.ndarray, N: int) -> float:
    """γ = ⟨ψ|X⊗N|ψ⟩."""
    XN = kron_n([SX] * N)
    return float(np.real(psi.conj() @ XN @ psi))


def pi2_odd_fraction(rho: np.ndarray, N: int) -> float:
    """Π²-odd Frobenius² fraction of rho in the Z-dephasing convention
    (bit_b parity = number of Y or Z letters mod 2)."""
    d = 2 ** N
    inv_d = 1.0 / d
    total = 0.0
    odd = 0.0
    for k in range(4 ** N):
        kk = k
        idxs = []
        for _ in range(N):
            idxs.append(kk & 3)
            kk >>= 2
        parity = sum(BIT_B[i] for i in idxs) & 1
        sigma = kron_n([PAULIS[i] for i in idxs])
        coeff = np.trace(sigma @ rho) * inv_d
        contrib = abs(coeff) ** 2 * d
        total += contrib
        if parity == 1:
            odd += contrib
    return odd / total if total > 0 else 0.0


# ---------------------- F86b formula + new depth-3 case ----------------------

def f86b_alpha_from_gamma(gamma: float) -> float:
    """F86b closed form: α_total = (1 − γ²)/2."""
    return (1 - gamma ** 2) / 2


# Anchor catalogue: canonical trigonometric angles for the F86b
# non-uniform Dicke parametrisation. Five anchors corresponding to the
# canonical {0°, 30°, 45°, 60°, 90°} trig angles produce α ∈ {0, 1/8, 1/4, 3/8, 1/2}.
TRIG_ANCHORS = [
    # (θ in degrees, γ = cos(θ), c², α = (1−γ²)/2, anchor name)
    (0,   1.0,          float("inf"),       0.0,      "Mirror (uniform Dicke at N odd, n=(N-1)/2)"),
    (30,  sqrt(3) / 2,  2 * sqrt(3) + 3,    0.125,    "DEPTH-3 / Pi2DyadicLadder a_3 (NEW tonight)"),
    (45,  sqrt(2) / 2,  1 + sqrt(2),        0.25,     "QuarterAsBilinearMaxval (non-uniform γ=√2/2)"),
    (60,  0.5,          1.0,                3.0 / 8.0, "KIntermediate (uniform Dicke c=1; today morning)"),
    (90,  0.0,          0.0,                0.5,      "Generic / HalfAsStructuralFixedPoint"),
]


# ---------------------- Verification loop ----------------------

def verify_all_canonical_anchors():
    """Numerically check the five F99 angles at N = 4, 6, 8."""
    print("=" * 92)
    print("  Verifying all FIVE canonical trigonometric anchors via non-uniform Dicke")
    print("=" * 92)
    print()

    for theta_deg, gamma_pred, c_sq, alpha_pred, name in TRIG_ANCHORS:
        if c_sq == float("inf"):
            print(f"  θ = {theta_deg}° ({name}): c → ∞ (Mirror = single Dicke state |D_{{N/2}}⟩), skipping non-uniform verification")
            continue
        c = sqrt(c_sq)
        gamma_pred_frac = Fraction(gamma_pred).limit_denominator(1000) if gamma_pred in (0, 0.5, 1) else None
        print(f"  θ = {theta_deg}°  ({name}):")
        print(f"    Predicted: γ = {gamma_pred:.10f}, α = {alpha_pred:.10f} (= {Fraction(alpha_pred).limit_denominator(100)}), c² = {c_sq:.6f}")
        print(f"    {'N':>3} {'n':>2} {'γ obs':>14} {'Δγ':>11} {'α obs':>14} {'Δα':>11}")
        for N in [4, 6, 8]:
            n = N // 2 - 1
            psi = dicke_nonuniform_superposition(N, n, c)
            gamma_obs = x_global_overlap(psi, N)
            rho = np.outer(psi, psi.conj())
            alpha_obs = pi2_odd_fraction(rho, N)
            gamma_dev = abs(gamma_obs - gamma_pred)
            alpha_dev = abs(alpha_obs - alpha_pred)
            g_ok = "✓" if gamma_dev < 1e-12 else "✗"
            a_ok = "✓" if alpha_dev < 1e-12 else "✗"
            print(f"    {N:>3} {n:>2} {gamma_obs:>14.10f} {gamma_dev:>10.2e}{g_ok} {alpha_obs:>14.10f} {alpha_dev:>10.2e}{a_ok}")
        print()


def display_trig_anchor_pattern():
    """Show the 30°-60°-90° standard-trigonometry pattern."""
    print("=" * 84)
    print("  F86b anchors at canonical trigonometric angles (30°-60°-90° triangle)")
    print("=" * 84)
    print()
    print(f"  {'θ':>3} {'γ = cos(θ)':<12} {'c² (Dicke)':<14} {'α = (1−γ²)/2':<14} {'Anchor'}")
    print(f"  {'-'*3} {'-'*12} {'-'*14} {'-'*14} {'-'*40}")
    for theta, gamma, c_sq, alpha, name in TRIG_ANCHORS:
        c_str = "∞" if c_sq == float("inf") else f"{c_sq:.6f}"
        alpha_frac = Fraction(alpha).limit_denominator(100)
        gamma_str = f"{gamma:.6f}" if gamma not in (0, 1, 0.5) else str(Fraction(gamma).limit_denominator(10))
        if abs(gamma - sqrt(3) / 2) < 1e-12:
            gamma_str = "√3/2"
        print(f"  {theta:>3}° {gamma_str:<12} {c_str:<14} {str(alpha_frac):<14} {name}")
    print()
    print("  The F86b α = (1 − γ²)/2 formula is parameterised by γ = ⟨ψ|X⊗N|ψ⟩.")
    print("  For γ ∈ {0, 1/2, √3/2, 1} = {cos(90°), cos(60°), cos(30°), cos(0°)} —")
    print("  the four canonical trigonometric angles of the 30°-60°-90° standard")
    print("  triangle — the formula produces the four polarity-anchor α values")
    print("  {1/2, 3/8, 1/8, 0} on the Pi2 dyadic ladder.")
    print()
    print("  Uniform Dicke (c = 1) realises the 60° case (KIntermediate, α = 3/8).")
    print("  Non-uniform Dicke c² = 2√3 + 3 realises the 30° case (depth-3, α = 1/8).")
    print("  The c² formula c² = cos(θ)/(2·sin²(θ/2)) parameterises the X⊗N-angle θ.")
    print()


def display_periodic_table_comparison():
    print("=" * 84)
    print("  Curated period-table fraction comparison — scope limited")
    print("=" * 84)
    print()
    print("  F99 derives α = 1/8 from the named non-uniform Dicke construction:")
    print("  c² = 2√3 + 3, γ = √3/2, α = (1 − γ²)/2 = 1/8.")
    print("  The calculation is independent of period-table labels.")
    print()
    print("  The companion script groups curated H--Ar occupation/slot labels as")
    print("  exact fractions. The following labels share the displayed fractions:")
    print()
    print("    α = 1/8  at θ = 30°  ↔ curated labels Li, Na, F, Cl")
    print("    α = 1/4  at θ = 45°  ↔ curated labels Be, Mg")
    print("    α = 3/8  at θ = 60°  ↔ curated labels B, Al, N, P")
    print("    α = 1/2  at θ = 90°  ↔ curated labels H, C, Si")
    print()
    print("  This reports equal rational labels only. It is not a material bridge, a")
    print("  spectral or dynamical prediction, or a producer/measurement result.")
    print("  Any such translation needs a specified degree of freedom, Hamiltonian,")
    print("  bath, preparation, measurement, and producer chain.")
    print()


def main():
    print()
    print("=" * 84)
    print("  F99: depth-3 anchor derivation via non-uniform Dicke")
    print("=" * 84)
    print()
    print("  Named construction: c² = 2√3 + 3, γ = √3/2 = cos(30°), α = 1/8.")
    print("  The finite-N output is a floating-point residual check of the closed form.")
    print()
    print()
    verify_all_canonical_anchors()
    display_trig_anchor_pattern()
    display_periodic_table_comparison()


if __name__ == "__main__":
    main()
