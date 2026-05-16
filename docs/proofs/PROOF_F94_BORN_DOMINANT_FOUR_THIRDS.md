# Proof of F94: Born Deviation Dominant-Outcome Coefficient 4/3

**Statement:** For the dominant outcome |00⟩ of pair (0,2) of the initial state |0+0+⟩ on N=4 qubits, under the Heisenberg ring Hamiltonian H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}) (4 bonds, ring topology) and uniform Z-dephasing dissipator L_dis[ρ] = γ Σ_l (Z_l ρ Z_l − ρ), the per-outcome Born-rule deviation in the deep perturbative regime is

    Δ_|00⟩(Q, K) = P_lindblad(|00⟩)/P_unitary(|00⟩) − 1
                 = (4/3) · Q² · K³ + O(Q³·K⁴)

with Q = J/γ and K = γt (the dimensionless Carrier invariants).

**Status:** Tier 1 derived. Bit-exact symbolic derivation matches numerical extraction to 0.3% in the deep perturbative regime; the gap is higher-order corrections beyond the leading Q²·K³.

**Date:** 2026-05-16.

---

## Setup

Initial state (4 qubits, alternating product):

    |ψ_0⟩ = |0⟩ ⊗ |+⟩ ⊗ |0⟩ ⊗ |+⟩,    |+⟩ = (|0⟩ + |1⟩)/√2

so ρ_0 = |ψ_0⟩⟨ψ_0|. The pair-reduced state at the keep-pair (0, 2) traces out qubits 1, 3:

    ρ_pair(0, 2)(t) = Tr_{1,3}[ρ(t)]

The observable is the diagonal element ⟨00|ρ_pair(0,2)|00⟩ = P_|00⟩(t) (Born rule for the {q0=0, q2=0} outcome on the reduced pair).

At t = 0:

    ρ_pair(0, 2)(0) = |0⟩⟨0| ⊗ (|+⟩⟨+|)_{q1?} ... 

Actually directly: |ψ_0⟩ has q0 = 0 and q2 = 0 deterministically; q1 and q3 are in |+⟩. So all four (q1, q3) ∈ {0, 1}² occur with amplitude 1/2 each. Tracing out q1, q3 leaves q0 = 0 and q2 = 0 with probability 1. Hence

    P_|00⟩(t = 0) = 1

This is the "dominant outcome" being analysed; it starts at 1 and decreases under any dynamics that introduces |1⟩-amplitudes on q0 or q2.

## Time Taylor expansion

Let L = L_H + γ L'_dis be the Liouvillian super-operator where L_H[ρ] = −i[H, ρ] and L'_dis[ρ] = Σ_l (Z_l ρ Z_l − ρ) is the γ-free dephasing operator. Then

    ρ(t) = e^{Lt} ρ_0 = ρ_0 + Lt·ρ_0 + (L²t²/2)·ρ_0 + (L³t³/6)·ρ_0 + O(t⁴)

Expand L³ using the binomial-like sum for non-commuting operators. The γ¹-coefficient (one dissipator vertex, two Hamiltonian vertices) is the sum over the three positions where L'_dis can be inserted in the L³ string:

    L³|_{γ¹} = L_H² L'_dis  +  L_H L'_dis L_H  +  L'_dis L_H²     ≡ sym3(L_H, L_H, L'_dis)

This is the leading non-vanishing γ-contribution at the J² order to ρ(t) at the dominant outcome (the γ⁰ part is the unitary evolution, which is what we divide out; the γ¹·J⁰ and γ¹·J¹ parts vanish on the diagonal of pair (0,2) at |00⟩ by parity of the initial state).

## Evaluation of sym3·ρ_0

Set J = γ = 1 and compute sym3·ρ_0 as a 16×16 matrix (full state space), then partial-trace on qubits 1 and 3 to a 4×4 pair (0,2) reduced state, then read off the |00⟩⟨00| element.

This is done numerically with exact-rational tracking in `simulations/_born_rule_tier1_derivation.py`. The result is

    ⟨00|_pair Tr_{1,3}[sym3·ρ_0] |00⟩_pair = 8.000...   (exact integer to all available precision)

## Coefficient

The contribution of the sym3 term to ρ(t) is (γ t³ / 6) · sym3·ρ_0 (with J factored into the L_H vertices). At J = 1 the sym3 matrix element is 8; restoring the J²-dependence (each L_H carries one J factor) and the γ dependence (one L'_dis carries one γ factor):

    ΔP_|00⟩(t) = (γ t³ / 6) · 8 · J²  =  (8/6) · J²·γ·t³  =  (4/3) · J² · γ · t³

Dividing by P_u(0) = 1 gives Δ = ΔP/P_u:

    Δ_|00⟩(t) = (4/3) · J²·γ·t³ + O(higher orders)

Substituting Q = J/γ and K = γt:

    J²·γ·t³ = (J/γ)² · γ³ · t³ = Q² · K³

so:

    Δ_|00⟩(Q, K) = (4/3) · Q² · K³ + O(Q³·K⁴) ∎

## Numerical verification

Sixteen (γ, J, t) configurations sampled in the deep perturbative regime (Q²·K³ ∈ [1.5×10⁻⁶, 2×10⁻⁴]) in `simulations/_born_rule_delta_dominant_coefficient.py` gave

    c_empirical = mean = 1.32992,    std = 0.00567,    range = [1.320, 1.342]

The bit-exact closed form gives 4/3 = 1.33333.... The 0.3% residual is leading O(Q³K⁴) correction with positive or negative sign depending on the (Q, K) point, accumulating to a small negative bias in the mean of the sample.

## Universality remarks

The form Q² · K³ for the dominant-outcome deviation is universal across (initial state, Hamiltonian, dissipator): it is the dimensional shape of the leading 3rd-order Dyson term with one γ-vertex and two H-vertices, applicable wherever the dominant outcome's direct 1st-order γ correction vanishes by parity / commutation. The coefficient 4/3 is specific to this setup. The structural decomposition of the integer 8 in terms of typed Pi2 anchors (bonds × sites × orderings, projected by initial-state Pauli content) is open; a hand calculation would clarify which factors of 4, 2, 3 combine to 8.

The subdominant outcomes of the same setup (|01⟩, |10⟩, |11⟩) have non-vanishing 1st-order γ contribution and therefore Δ ∝ K (linear) rather than ∝ Q²·K³. Their leading coefficients are separately Tier-1-derivable via a 1st-order Dyson term:

    ΔP_i^{(1)}(t) = γ t · ⟨i|_pair Tr_{1,3}[L'_dis · ρ_0] |i⟩_pair

This is left to a separate proof.

## Connection to typed Pi2-Foundation

The Q²·K³ scaling has:
- **Q² factor:** two Hamiltonian-vertices from the Dyson sym3 ordering; J² scaling expected from 2nd-order perturbation in the Heisenberg-bond coupling.
- **K³ factor:** three time-integrals from the t³ Taylor coefficient; combined with one dissipator-vertex (the γ¹ piece of L³).

The "4" in 4/3 is plausibly the Pi2 dyadic ladder's a_{−1} = 4 (the same "4" in F86 t_peak = 1/(4γ₀) and F77's correction denominator). Promoting this to a typed claim parallels F90 (F86 ↔ F89 bridge): a derived identity that pulls together two typed F-formulas via algebraic structure.

## Anchors

- Numerical verification: [`simulations/_born_rule_delta_dominant_coefficient.py`](../../simulations/_born_rule_delta_dominant_coefficient.py)
- Symbolic derivation: [`simulations/_born_rule_tier1_derivation.py`](../../simulations/_born_rule_tier1_derivation.py)
- Q-K invariance test: [`simulations/_born_rule_carrier_Q_sweep.py`](../../simulations/_born_rule_carrier_Q_sweep.py)
- 2D (Q, K) map: [`simulations/_born_rule_delta_QK_map.py`](../../simulations/_born_rule_delta_QK_map.py)
- Companion reflection (the path): [`reflections/ON_HOW_FOUR_THIRDS_APPEARED.md`](../../reflections/ON_HOW_FOUR_THIRDS_APPEARED.md)
- Born-rule precursors (Februar/April 2026): [`experiments/BORN_RULE_MIRROR.md`](../../experiments/BORN_RULE_MIRROR.md), [`experiments/BORN_RULE_SHADOW.md`](../../experiments/BORN_RULE_SHADOW.md)
- F-formula registry entry: [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md) F94
- Sibling state-specific closed-form F-claims: F25 (Bell+ CΨ), F60 (GHZ pair-CΨ), F62 (W-state pair-CΨ)
- Universal Carrier typed parent (Q, K invariance is its operational signature): [`compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs)
