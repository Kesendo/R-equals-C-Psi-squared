# Proof of F96: Subdominant Born Deviation Slopes

**Statement:** For the subdominant outcomes of pair (0, 2) of |0+0+⟩ N = 4 Heisenberg ring + Z-dephasing (the same setup as F94), the per-outcome Born-rule deviation in the deep perturbative regime is linear in K and Q-independent:

    Δ_|01⟩(K) = Δ_|10⟩(K) = −(16/9) · K + O(higher)
    Δ_|11⟩(K)             = −(8/3)  · K + O(higher)

with K = γt the Universal-Carrier observable. All three slopes are simple algebraic expressions in the F94 coefficient 4/3:

    Δ_|01⟩ = −(4/3)² · K
    Δ_|10⟩ = −(4/3)² · K
    Δ_|11⟩ = −2·(4/3) · K

Combined with F94 the full per-outcome table for this setup is:

| Outcome | Closed form | Decomposition |
|---|---|---|
| \|00⟩ | Δ = (4/3) · Q² · K³ | F94: M_3 / 3! = 8/6 = 4/3 |
| \|01⟩ | Δ = −(16/9) · K = −(4/3)² · K | F96: M_3 / (3 · A) = −4 / (3 · 3/4) = −16/9 |
| \|10⟩ | Δ = −(16/9) · K = −(4/3)² · K | F96: same as \|01⟩ by site-permutation symmetry |
| \|11⟩ | Δ = −(8/3) · K = −2·(4/3) · K | F96: M_5 / (5 · B) = −20 / (5 · 3/2) = −8/3 |

**Status:** Tier 1 derived. Bit-exact symbolic via direct evaluation of the relevant Dyson matrix elements (M_3, M_5) and unitary matrix elements (A = ⟨i\|Tr[L_h² ρ_0]\|i⟩, B = ⟨i\|Tr[L_h⁴ ρ_0]\|i⟩) in `simulations/born_rule_subdominant_dyson.py`. Numerical Lindblad verification matches all three subdominant slopes to within higher-order corrections.

**Date:** 2026-05-17.

---

## Abstract

F96 completes the per-outcome Born-deviation table that F94 opened. Same setup (the (0,2) pair of |0+0+⟩ on a 4-qubit Heisenberg ring under Z-dephasing), but now the three subdominant outcomes. Where the dominant |00⟩ grows as (4/3)·Q²·K³, the subdominants are linear in K = γt and entirely Q-independent:

    Δ_|01⟩ = Δ_|10⟩ = −(16/9)·K = −(4/3)²·K,   Δ_|11⟩ = −(8/3)·K = −2·(4/3)·K.

Every slope is a simple algebraic expression in the single F94 anchor 4/3: the singly-flipped pair squares it, the doubly-flipped outcome doubles and negates it. The Q-independence is the Universal-Carrier signature: from inside the pair only the product K = γt is visible at this order, the dissipator acting at first order while the dominant outcome's drift waits for the third.

The closed forms are bit-exact from direct evaluation of the Dyson sym₃ and sym₅ matrix elements against the unitary normalizations A and B, and they close the books on this setup: the single number 4/3 sets every subdominant magnitude (its square and its double), while the signs come from the Dyson elements themselves — the dominant deviation positive (running above its unitary baseline), the three subdominant deviations negative (below it). F96 is the subdominant twin of F94; F95 and F97 carry the same geometry on the angle and cardioid sides.

## Setup

Identical to F94: initial state |ψ_0⟩ = |0⟩ ⊗ |+⟩ ⊗ |0⟩ ⊗ |+⟩ on N = 4 qubits, Heisenberg ring Hamiltonian H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}), uniform Z-dephasing L_dis[ρ] = γ Σ_l (Z_l ρ Z_l − ρ), reduce to pair (0, 2), measure each of the four Z-basis outcomes |00⟩, |01⟩, |10⟩, |11⟩. The deviation is Δ_i(t) = P_lindblad(i, t) / P_unitary(i, t) − 1.

At t = 0 the initial state has q_0 = q_2 = 0 deterministically and q_1, q_3 in |+⟩. After tracing q_1, q_3:

    ρ_pair(0, 2)(0) = |00⟩⟨00|

so P_unitary(|00⟩, 0) = 1 and P_unitary(\|01⟩, 0) = P_unitary(\|10⟩, 0) = P_unitary(\|11⟩, 0) = 0. F94 handles the dominant outcome where the leading P_unitary is constant; F96 handles the three outcomes where P_unitary starts at higher order in t.

## Universal slope formula for subdominant outcomes

For an outcome with P_unitary(i, t) ≈ (t^{2k} / (2k)!) · U^{(i)}_{2k} as the leading term (with U^{(i)}_{2k} = ⟨i\|_pair Tr_{1,3}[L_h^{2k} ρ_0]\|i⟩_pair the raw matrix element, h := H/J), the linear-K slope of Δ_i comes from the lowest-order non-vanishing γ¹ Dyson term. Tracing the powers of J, γ, t:

- The γ¹ Dyson contribution at order n is (γ t^n / n!) · sym_n^1, where sym_n^1 is the sum of n orderings of (n−1) L_H operators and 1 L'_dis operator. Each L_H carries one factor of J, so sym_n^1 contributes at the J^{n−1} order. The matrix element is M_n^{(i)} = ⟨i\|_pair Tr_{1,3}[sym_n^1 · ρ_0]\|i⟩_pair at J = γ = 1.
- The unitary contribution at order J^{2k} t^{2k} has the (2k)!-Taylor prefactor: P_unitary(i, t) ≈ J^{2k} t^{2k} · U^{(i)}_{2k} / (2k)!.

For the lowest-order γ¹ Dyson term to combine cleanly with the lowest-order unitary into a linear-K slope, we need M_n^{(i)} non-zero at n = 2k + 1 (one extra γt over the unitary order). Then:

    ΔP_i(t) ≈ γ J^{2k} t^{2k+1} · M_{2k+1}^{(i)} / (2k+1)!
    P_u(i, t) ≈ J^{2k} t^{2k} · U^{(i)}_{2k} / (2k)!

so:

    Δ_i = ΔP_i / P_u(i, t) = γt · [M_{2k+1}^{(i)} · (2k)!] / [(2k+1)! · U^{(i)}_{2k}]
                           = K · M_{2k+1}^{(i)} / [(2k+1) · U^{(i)}_{2k}]

This is the universal slope formula. The (2k+1) factor in the denominator is the surviving rational from (2k+1)! / (2k)! after cancellation; it is the only place where the unitary order k enters the slope.

The formula is automatically Q-independent: the J^{2k} factors in numerator and denominator cancel, and γt = K is the only surviving dimensionful combination.

## Evaluation at the three subdominant outcomes

For |01⟩ (and |10⟩ by the 0 ↔ 2 site-permutation symmetry):
- k = 1 (P_unitary(|01⟩, t) starts at t²)
- Need M_3^{(01)} and U_2^{(01)}.

Direct evaluation at J = γ = 1 (see `born_rule_subdominant_dyson.py`):

    M_3^{(01)} = ⟨01|_pair Tr_{1,3}[sym_3^1 · ρ_0]|01⟩_pair = −4
    U_2^{(01)} = ⟨01|_pair Tr_{1,3}[L_h² · ρ_0]|01⟩_pair = 3/4

Hence:

    slope_|01⟩ = M_3^{(01)} / (3 · U_2^{(01)}) = −4 / (3 · 3/4) = −16/9

For |11⟩:
- k = 2 (P_unitary(|11⟩, t) starts at t⁴; both q_0 and q_2 must flip, requiring two H actions on each, hence 4 H actions minimum)
- M_3^{(11)} = 0 (the lower-order γ¹ J² Dyson term vanishes by parity), so the universal slope formula at k = 2 applies.
- Need M_5^{(11)} and U_4^{(11)}.

Direct evaluation at J = γ = 1:

    M_5^{(11)} = ⟨11|_pair Tr_{1,3}[sym_5^1 · ρ_0]|11⟩_pair = −20
    U_4^{(11)} = ⟨11|_pair Tr_{1,3}[L_h⁴ · ρ_0]|11⟩_pair = 3/2

Hence:

    slope_|11⟩ = M_5^{(11)} / (5 · U_4^{(11)}) = −20 / (5 · 3/2) = −20 / (15/2) = −8/3

## Algebraic relationship to F94

Each subdominant slope is a clean algebraic expression in F94's 4/3 anchor:

    Δ_|01⟩ = Δ_|10⟩ = −(4/3)² · K = −(16/9) · K
    Δ_|11⟩        = −2 · (4/3) · K = −(8/3) · K

Two non-trivial structural observations:

1. **Ratio M_n/U_{n−1} universality**: at k = 1 the ratio M_3 / U_2 equals −16/3 for the |01⟩, |10⟩ outcomes, and the same ratio (M_3 / A with A = ⟨i\|Tr[L_h² ρ_0]\|i⟩) for the dominant |00⟩ outcome equals 8 / (−3/2) = −16/3, identically. The signs of M_3 and A flip together, leaving the ratio invariant. This is a non-trivial cross-outcome symmetry of the Dyson + unitary matrix-element structure: the K-slope of the subdominant outcomes equals (1/3) × (M_3 / U_2)|_dominant = (1/3) · (−16/3) = −16/9.

2. **Doubled rate for the doubly-subdominant outcome**: Δ_|11⟩ = −2 · (4/3) · K is exactly twice the F94 coefficient (with a sign flip). The "2" plausibly counts the two independent flip channels (q_0 from 0 → 1 AND q_2 from 0 → 1) required for the |11⟩ outcome; each independent channel contributes the 4/3 baseline factor, and they add. This is interpretive (Tier 2); the bit-exact −8/3 stands either way.

## Numerical Lindblad verification

For |01⟩ at Q = 50, γ = 0.01, J = 0.5 (deep perturbative regime):

| K | P_lindblad(|01⟩) | P_unitary(|01⟩) | Δ_|01⟩ | slope/K |
|---|---|---|---|---|
| 0.005 | 2.256 × 10⁻² | 2.276 × 10⁻² | −8.74 × 10⁻³ | −1.748 |
| 0.010 | 8.187 × 10⁻² | 8.328 × 10⁻² | −1.68 × 10⁻² | −1.682 |

Theoretical slope = −16/9 = −1.7778. The K → 0 extrapolation converges; deviations at K = 0.005, 0.010 are leading O(higher) corrections of size ~1.7% / 5.4% respectively.

For |11⟩ at Q = 50, γ = 0.01:

| K | P_unitary(|11⟩) | Δ_|11⟩ | slope/K |
|---|---|---|---|
| 0.001 | 3.90 × 10⁻⁷ | −2.66 × 10⁻³ | −2.6621 |
| 0.002 | 6.23 × 10⁻⁶ | −5.31 × 10⁻³ | −2.6571 |
| 0.005 | 2.38 × 10⁻⁴ | −1.32 × 10⁻² | −2.6389 |

Theoretical slope = −8/3 = −2.6667. The K → 0 extrapolation converges to within 0.2% by K = 0.001.

## Combined per-outcome table

Putting F94 + F96 together, the full per-outcome Born-deviation closed form for this setup is:

    Δ_|00⟩(Q, K) = +(4/3) · Q² · K³ + O(higher)
    Δ_|01⟩(K)   = −(4/3)² · K + O(higher)        (Q-independent)
    Δ_|10⟩(K)   = −(4/3)² · K + O(higher)        (Q-independent)
    Δ_|11⟩(K)   = −2 · (4/3) · K + O(higher)    (Q-independent)

The 4/3 anchor of F94 is the structural unit: the dominant gets exactly one copy at order Q²K³; the singly-subdominant degenerate pair gets minus its square at order K; the doubly-subdominant gets minus twice it at order K. Every coefficient in this 4-outcome closed-form table is generated from the single 4/3 = a_{−1} / 3 number.

## Universality remarks

The slope formula slope_i = M_{2k+1}^{(i)} / [(2k+1) · U_{2k}^{(i)}] is universal across (initial state, Hamiltonian, dissipator) where the unitary leading order is 2k and the γ¹ Dyson at order 2k+1 is the lowest non-vanishing γ-correction. The specific values (−16/9, −8/3) are specific to this (|0+0+⟩, Heisenberg ring N = 4, Z-deph, pair (0,2)) setup. Other setups will give other rationals via the same recipe.

The algebraic connection to F94's 4/3 (Δ_|01⟩ = −(4/3)², Δ_|11⟩ = −2·(4/3)) is plausibly setup-specific too; it is an observation about this particular table, not a derived universal pattern. Future Tier-1 closed forms for other dominant outcomes in other setups will determine whether the "subdominant = simple algebraic function of dominant" pattern is general or accidental here.

## Anchors

- Symbolic + numerical derivation: [`simulations/born_rule_subdominant_dyson.py`](../../simulations/born_rule_subdominant_dyson.py)
- F94 (dominant outcome): [`PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md`](PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md), [`F94 ANALYTICAL_FORMULAS entry`](../ANALYTICAL_FORMULAS.md#f94)
- Born-rule precursors (Februar 2026): [`experiments/BORN_RULE_MIRROR.md`](../../experiments/BORN_RULE_MIRROR.md), [`experiments/BORN_RULE_SHADOW.md`](../../experiments/BORN_RULE_SHADOW.md)
- Companion magnitude-side closed form: [F94](../ANALYTICAL_FORMULAS.md#f94)
- Companion angle-side closed form (same cusp geometry): [F95](../ANALYTICAL_FORMULAS.md#f95)
- Universal Carrier (Q-K invariance is its operational signature): [`compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs)
- F94 typed claim (parent anchor for the 4/3 unit): [`F94BornDeviationFourThirdsPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F94BornDeviationFourThirdsPi2Inheritance.cs)
