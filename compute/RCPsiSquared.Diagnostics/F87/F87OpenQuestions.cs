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
            "F81 (Π·M·Π⁻¹ = M − 2L_{H_odd}) extends to T1 amplitude damping (F82) and thermal " +
            "amplitude damping (F84), but only for single-qubit Lindblad operators. Two-qubit " +
            "dissipators (correlated decay, ZZ-cross-channel, swap channels) remain analytically " +
            "open.",
            "Apply the F81 commutator decomposition to two-qubit Lindblad operators; check " +
            "whether the same anti-commutator identity holds with adjusted coefficients.",
            "docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md + PROOF_F82_T1_DISSIPATOR_CORRECTION.md"),
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
