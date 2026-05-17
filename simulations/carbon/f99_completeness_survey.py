"""Are the five F99 canonical anchors {0, 1/8, 1/4, 3/8, 1/2} ALL of them?

Tom: "ich kann mir fast nicht vorstellen das es alle sind."

Two questions:

  (1) For the F86b α = (1−γ²)/2 = sin²(θ)/2 formula with γ = cos(θ):
      which θ (rational multiples of π) give α ∈ ℚ (a dyadic fraction)?

  (2) For OTHER state classes (non-Dicke superpositions giving different
      α(γ) functions): can the framework reach depth-4 dyadic anchors
      {1/16, 3/16, ..., 15/16} via different algebraic mechanisms?

ANSWER TO (1) — Niven's theorem (1956):
  For θ a rational multiple of π, cos(θ) is rational only for θ ∈ {0°,
  60°, 90°, 120°, 180°} (modulo reflections). Applied to 2θ: sin²(θ) =
  (1 − cos(2θ))/2 is rational iff 2θ ∈ {0°, 60°, 90°, 120°, 180°}, i.e.
  θ ∈ {0°, 30°, 45°, 60°, 90°}. THE FIVE F99 ANCHORS ARE NIVEN-COMPLETE.

  No other rational-multiple-of-π angle produces a rational α via the F86b
  formula. The 30°-60°-90° + 45°-45°-90° standard triangles ARE the
  Niven-rational set for 2θ.

ANSWER TO (2) — alternative α(γ) formulas:
  The F86b α(γ) = (1−γ²)/2 came from the SPECIFIC X⊗N-eigenbasis
  decomposition of the symmetric two-Dicke superposition. Other quantum
  state classes give DIFFERENT α(γ) functions, potentially with different
  rationality structure. We survey W-state, GHZ-state, and Bell-pair
  classes here.

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/f99_completeness_survey.py
"""
from __future__ import annotations

import sys
import numpy as np
from fractions import Fraction
from math import comb, cos, sin, pi, sqrt

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ---------------------- Pauli + Π² split (same as depth_3_anchor_derivation.py) ----------------------

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, SX, SY, SZ]
BIT_B = [0, 0, 1, 1]


def kron_n(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def dicke_state(N, k):
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    norm = 1.0 / sqrt(comb(N, k))
    for b in range(d):
        if bin(b).count("1") == k:
            psi[b] = norm
    return psi


def x_global_overlap(psi, N):
    XN = kron_n([SX] * N)
    return float(np.real(psi.conj() @ XN @ psi))


def pi2_odd_fraction(rho, N):
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


# ---------------------- Part (1): Niven survey on F86b formula ----------------------

def niven_completeness_survey():
    """Survey θ = k·15° for k = 0..6 (covers all elementary constructible angles
    in [0°, 90°] at 15° resolution). For each, compute sin²(θ)/2 and check if it's
    rational (= dyadic in the framework's algebra)."""
    print("=" * 86)
    print("  PART (1): Niven survey — which θ give rational α = sin²(θ)/2?")
    print("=" * 86)
    print()
    print("  Survey at 15° resolution covers all elementary constructible angles in [0°, 90°]:")
    print("  {0°, 15°, 30°, 45°, 60°, 75°, 90°}")
    print()
    print(f"  {'θ':>4} {'sin²(θ)':<28} {'α = sin²(θ)/2':<28} {'Niven-rational?':<14}")
    print(f"  {'-'*4} {'-'*28} {'-'*28} {'-'*14}")

    for deg in range(0, 91, 15):
        theta = deg * pi / 180
        sin_sq = sin(theta) ** 2
        alpha = sin_sq / 2
        # Try to identify as a rational with small denominator
        f_sinsq = Fraction(sin_sq).limit_denominator(1000)
        f_alpha = Fraction(alpha).limit_denominator(1000)
        # "Niven-rational" = matches a clean rational at very tight tolerance
        sin_sq_rational = abs(float(f_sinsq) - sin_sq) < 1e-12
        # Convert to nice display
        if sin_sq_rational:
            sin_str = str(f_sinsq)
            alpha_str = str(f_alpha)
            marker = "✓ rational"
        else:
            # Identify algebraic form
            sin_str = f"{sin_sq:.10f}"
            alpha_str = f"{alpha:.10f}"
            if abs(sin_sq - (2 - sqrt(3))/4) < 1e-12:
                sin_str = "(2-√3)/4"
            elif abs(sin_sq - (2 + sqrt(3))/4) < 1e-12:
                sin_str = "(2+√3)/4"
            marker = "✗ irrational"
        print(f"  {deg:>3}° {sin_str:<28} {alpha_str:<28} {marker:<14}")

    print()
    print("  Result: exactly {0°, 30°, 45°, 60°, 90°} produce rational α.")
    print("  This matches Niven's theorem (1956): cos(θ) for θ ∈ ℚ·π is rational only")
    print("  for θ ∈ {0°, 60°, 90°, 120°, 180°} ∪ reflections. Apply to 2θ: sin²(θ) =")
    print("  (1 − cos(2θ))/2 is rational ⟺ cos(2θ) rational ⟺ θ ∈ {0°, 30°, 45°, 60°, 90°}.")
    print()
    print("  THE FIVE F99 ANCHORS ARE NIVEN-COMPLETE for the α = (1−γ²)/2 formula.")
    print("  Tom's intuition 'kann mir fast nicht vorstellen das es alle sind' is correct")
    print("  in the spirit but wrong in the letter: they ARE all of them — for THIS formula.")
    print("  To reach further anchors (depth-4 dyadic 1/16, 3/16, ...), the framework")
    print("  needs a STRUCTURALLY DIFFERENT α(γ) mechanism, not just more angles.")
    print()


def survey_finer_resolution():
    """Sanity check at finer resolution: 7.5° steps. Should still hit only the same five."""
    print("=" * 86)
    print("  Finer-resolution sanity check at 7.5° (constructible 'half' angles included)")
    print("=" * 86)
    print()
    print(f"  {'θ':>5} {'sin²(θ)':<28} {'α':<22} {'Status'}")
    print(f"  {'-'*5} {'-'*28} {'-'*22} {'-'*30}")
    rational_hits = []
    for half_deg in range(0, 181):  # 0, 0.5, 1, ..., 90 degrees in 0.5° steps
        deg = half_deg / 2
        if deg > 90:
            continue
        theta = deg * pi / 180
        sin_sq = sin(theta) ** 2
        alpha = sin_sq / 2
        f_sinsq = Fraction(sin_sq).limit_denominator(10000)
        if abs(float(f_sinsq) - sin_sq) < 1e-12 and f_sinsq.denominator <= 16:
            rational_hits.append((deg, sin_sq, alpha, f_sinsq))
    for deg, sin_sq, alpha, f_sinsq in rational_hits:
        f_alpha = Fraction(alpha).limit_denominator(100)
        print(f"  {deg:>5}° {str(f_sinsq):<28} {str(f_alpha):<22} ✓ rational hit")
    print()
    print(f"  Total rational hits in [0°, 90°] at 0.5° resolution: {len(rational_hits)}")
    print(f"  Expected per Niven (constructible angles): 5 ({{0°, 30°, 45°, 60°, 90°}})")
    print(f"  Match: {'✓' if len(rational_hits) == 5 else '✗'}")
    print()


# ---------------------- Part (2): alternative α(γ) formulas ----------------------

def w_state(N):
    """W-state |W_N⟩ = (1/√N) Σ_k |0...010...0⟩ (one excitation, equal superposition)."""
    return dicke_state(N, 1)


def ghz_state(N):
    """GHZ-state |GHZ_N⟩ = (|0...0⟩ + |1...1⟩)/√2."""
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / sqrt(2)
    psi[d - 1] = 1 / sqrt(2)
    return psi


def bell_pair_like(N):
    """For even N: pairs of qubits each in Bell-+: ψ = (|00⟩+|11⟩)/√2 ⊗ ... ⊗ same."""
    if N % 2 != 0:
        raise ValueError("Bell-pair requires even N.")
    bell = np.array([1, 0, 0, 1], dtype=complex) / sqrt(2)
    psi = bell
    for _ in range(N // 2 - 1):
        psi = np.kron(psi, bell)
    return psi


def alternative_alpha_survey():
    """Compute γ and α for W-state, GHZ-state, Bell-pair-product for N = 2, 4, 6.
    Compare to F86b's α = (1 − γ²)/2 prediction. If the state class follows the
    F86b formula, γ → α is determined. If a state gives α ≠ (1−γ²)/2 at the
    measured γ, that state class has a DIFFERENT α(γ) function — a new algebraic
    mechanism. Such mechanisms could open up new dyadic anchors at depth-4+."""
    print("=" * 86)
    print("  PART (2): alternative state-class α(γ) survey")
    print("=" * 86)
    print()
    print("  For non-Dicke state classes, compute γ = ⟨ψ|X⊗N|ψ⟩ + α = Π²-odd Frobenius²")
    print("  fraction. Compare to F86b prediction α_F86b = (1 − γ²)/2. Deviation reveals")
    print("  a structurally different α(γ) function.")
    print()
    print(f"  {'State':<22} {'N':>3} {'γ':>14} {'α (observed)':>14} {'α F86b pred':>14} {'Δ':>10}")
    print(f"  {'-'*22} {'-'*3} {'-'*14} {'-'*14} {'-'*14} {'-'*10}")

    for N in [2, 4, 6]:
        # W-state
        psi = w_state(N)
        rho = np.outer(psi, psi.conj())
        gamma_w = x_global_overlap(psi, N)
        alpha_w = pi2_odd_fraction(rho, N)
        alpha_pred_w = (1 - gamma_w ** 2) / 2
        print(f"  {'W-state':<22} {N:>3} {gamma_w:>14.10f} {alpha_w:>14.10f} {alpha_pred_w:>14.10f} {abs(alpha_w - alpha_pred_w):>10.2e}")

        # GHZ-state
        psi = ghz_state(N)
        rho = np.outer(psi, psi.conj())
        gamma_g = x_global_overlap(psi, N)
        alpha_g = pi2_odd_fraction(rho, N)
        alpha_pred_g = (1 - gamma_g ** 2) / 2
        print(f"  {'GHZ-state':<22} {N:>3} {gamma_g:>14.10f} {alpha_g:>14.10f} {alpha_pred_g:>14.10f} {abs(alpha_g - alpha_pred_g):>10.2e}")

        # Bell-pair product (N even)
        if N % 2 == 0:
            psi = bell_pair_like(N)
            rho = np.outer(psi, psi.conj())
            gamma_b = x_global_overlap(psi, N)
            alpha_b = pi2_odd_fraction(rho, N)
            alpha_pred_b = (1 - gamma_b ** 2) / 2
            print(f"  {'Bell-pair product':<22} {N:>3} {gamma_b:>14.10f} {alpha_b:>14.10f} {alpha_pred_b:>14.10f} {abs(alpha_b - alpha_pred_b):>10.2e}")
        print()

    print("  Reading:")
    print("    - W-state: deviates from F86b prediction iff γ_W ≠ 0 (when N=2 gives Bell-like)")
    print("    - GHZ-state: γ = 1 exactly (X⊗N is GHZ stabiliser), so F86b predicts α = 0")
    print("      regardless of N. Empirical α should match.")
    print("    - Bell-pair product: γ = 1 (each pair is X⊗X eigenstate), F86b predicts α = 0.")
    print()
    print("  If all states match α_F86b = (1−γ²)/2: these classes ARE F86b-derivable,")
    print("  just at different γ values than uniform Dicke. The same algebraic family.")
    print("  If any state has α ≠ (1−γ²)/2: that class has a different α(γ) function and")
    print("  is a candidate for opening depth-4 anchors via a new mechanism.")
    print()


# ---------------------- Part (3): outlook ----------------------

def outlook():
    print("=" * 86)
    print("  OUTLOOK: where to look for depth-4 dyadic anchors")
    print("=" * 86)
    print()
    print("  F99 is provably complete for the F86b α(γ) = (1−γ²)/2 mechanism via Niven.")
    print("  To reach depth-4 dyadic anchors {1/16, 3/16, 5/16, ...}, the framework needs")
    print("  a structurally different α(γ) function whose rational set lands at depth-4.")
    print()
    print("  Candidate mechanisms to explore (future sessions):")
    print()
    print("  (i)  Higher-degree polynomial in γ. E.g., α = γ²(1−γ²)/something, or")
    print("       α = (1 − γ⁴)/k. These have richer rationality structure: γ⁴ rational")
    print("       is implied by γ² rational, but the inverse map γ → α is different.")
    print()
    print("  (ii) Multi-amplitude state class. Three-Dicke superposition")
    print("       ψ = (a|D_n⟩ + b|D_{n+1}⟩ + c|D_{n+2}⟩)/N has two free amplitude ratios,")
    print("       giving a 2D parameter space. Possibly some 2D rationality structure")
    print("       lands at new dyadic anchors.")
    print()
    print("  (iii) Non-pure states (mixed-state α). For ρ a mixture rather than pure")
    print("       state, the F86b derivation doesn't directly apply. A mixed-state")
    print("       analog could parametrise α via thermal weight.")
    print()
    print("  (iv) Different Lindblad classes. F86b was derived for chain XY + Z-deph")
    print("       'truly' class. Other Lindblad classes (T1 amplitude damping, depolarising)")
    print("       give different α(γ) on the same state class. F82-F84 already cover the")
    print("       T1 correction; explicit depth-4 anchors might fall out of those.")
    print()
    print("  Tom's intuition stands at a meta-level: there are likely MORE anchors")
    print("  reachable by the framework, just not via THIS algebraic mechanism (Niven-")
    print("  closes it). The question shifts from 'are there more F99-like angles?' to")
    print("  'are there other α(γ) formulas with their own complete anchor sets?'")
    print()


def main():
    print()
    print("=" * 86)
    print("  F99 completeness survey — is the five-anchor set really all of them?")
    print("=" * 86)
    print()
    niven_completeness_survey()
    survey_finer_resolution()
    alternative_alpha_survey()
    outlook()


if __name__ == "__main__":
    main()
