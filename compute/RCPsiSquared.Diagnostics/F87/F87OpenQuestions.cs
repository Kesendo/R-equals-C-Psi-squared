using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Open theoretical items for F87 and the F-chain it indexes: directions where the
/// trichotomy or its downstream F-formulas (F80-F85) need extension or proof.</summary>
public static class F87OpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "F80 cluster sign-walk at k≥3",
            "F80 derives Spec(M) = ±2i·Spec(H_non-truly) for chain Π²-odd 2-body bilinears. " +
            "Verified bit-exact at k=3 (N=4,5,6) and k=4 (N=5,6) for 17 Π²-odd cases, but the " +
            "closed-form Bloch sign-walk enumeration at k≥3 is not done.",
            "Generalise the F80 OBC sine-mode derivation to k≥3 monomials; check whether the " +
            "single-particle Spec(H_non-truly) at higher k still arises from a tractable matrix.",
            "docs/proofs/PROOF_F85_KBODY_GENERALIZATION.md k≥3 cluster section"),
        new OpenQuestion(
            "F81/F82/F84 two-qubit dissipators",
            "F81's identity extends to Π·M·Π⁻¹ = M − 2L_{H_odd} − 2D_odd for any dissipator, two-site " +
            "ones included, with the Π²-odd part D_odd read off each jump as the cross terms " +
            "between Pauli strings of opposite bit_b (F82 and PROOF_F82, 'Other dissipators'). " +
            "Jumps of one parity have none, by that rule: ZZ-dephasing, the swap jump " +
            "(II + XX + YY + ZZ)/2, and Z dephasing with correlated site fields, jumps built " +
            "linearly from the Z_k (no odd part for any coefficients; for a real Kossakowski matrix " +
            "a over the Z_k the F1 residual of the dissipator against 2·tr(a) has " +
            "‖M_D‖² = 8·4^N·Σ_{k<l} a_kl², experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md, typed as CorrelatedZDephasingPriceClaim). The " +
            "pair-decay jump σ⁻⊗σ⁻ and its siblings σ⁺⊗σ⁺, σ⁻⊗σ⁺, σ⁺⊗σ⁻, each on one pair of " +
            "sites, have ‖D_odd‖² = (3/2)·γ²·4^(N−2) (derived for σ⁻⊗σ⁻ in PROOF_F84 by PROOF_F82's " +
            "'Other dissipators' rule, the value simulations/2qubit_dissipator_exploration.py measures for " +
            "σ⁻⊗σ⁻, σ⁺⊗σ⁺ and σ⁻⊗σ⁺, σ⁺⊗σ⁻ following by the site swap). " +
            "What stays open is a closed " +
            "form where such jumps overlap, their odd entries landing on the same transitions so " +
            "that they do not add (PROOF_F84), and for other parity-mixing two-site jumps such as " +
            "collective emission σ⁻_k + σ⁻_l.",
            "Count the odd entries of overlapping pair jumps and of collective emission by the " +
            "cross-term rule and derive ‖D_odd‖² in closed form, as F82 does for single-site T1.",
            "docs/ANALYTICAL_FORMULAS.md F82 'Other dissipators' + " +
            "docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md 'Other dissipators' + " +
            "docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md (two-qubit dissipators)"),
        new OpenQuestion(
            "F83 topology generalization at higher body",
            "F83 anti-fraction = 1/(2+4r) is verified bit-exact for chain/ring/star/K_4 at " +
            "N=4 (k=2-body). Higher-body coefficients beyond n_YZ ∈ {0, 1, 2} are empirical-only.",
            "Sample mixed-Π²-class Hamiltonians at k=3, 4 across multiple topologies; check " +
            "whether the 1/(2+4r) form still holds with the F85-style Π²-class redefinition.",
            "docs/proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md line 151"),
        new OpenQuestion(
            "EQ-030 per-backend Δ(soft−truly) amplification predictor",
            "Marrakesh measured Δ = −0.72, Kingston −0.92, Fez −0.81. The ratification across " +
            "three backends is itself open: what predicts the per-backend amplification factor? " +
            "Likely depends on the Trotter discretisation parameter δt and the backend-specific " +
            "‖H‖_op profile (the original interpretation that T1 amplifies Δ was REFUTED; T1 " +
            "monotonically attenuates).",
            "Sweep δt vs Δ predictions on the simulator; check correlation with reported " +
            "backend two-qubit error rates and median pulse durations.",
            "review/EMERGING_QUESTIONS.md EQ-030 + lebensader_skeleton_trace_decoupling Confirmation"),
        new OpenQuestion(
            "4-way Π²-class enumeration as a typed C# enum",
            "The Python framework distinguishes pi2_odd_pure / pi2_even_nontruly / mixed / truly " +
            "as a 4-way refinement of F87's 3-way trichotomy (Truly / Soft / Hard). C# currently " +
            "exposes only TrichotomyClass (3-way) and PauliPairBondTerm.Pi2Parity (binary).",
            "Add Pi2Class4Way enum to RCPsiSquared.Core.Pauli and a classifier that dispatches " +
            "across a term list; surface as typed witnesses parallel to F87CanonicalWitness.",
            "Gap identified during F87KB cleanup 2026-05-03"),
    };
}
