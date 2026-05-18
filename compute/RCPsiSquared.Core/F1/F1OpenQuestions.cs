using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>Open theoretical items for the F1 family: directions where the palindrome
/// identity breaks, partially survives, or generalises. Anchored at the F1 "Breaks for"
/// clause and the OPERATOR_RIGIDITY scaling experiment.
///
/// <para>The earlier "T1 amplitude damping: full closed form" item was closed on
/// 2026-05-18 by <c>docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md</c> and now lives
/// as the Tier-1-derived <see cref="F1T1ResidualClosedForm"/> claim on
/// <see cref="F1KnowledgeBase"/>.</para>
///
/// <para>The earlier "depolarizing noise: residual scaling" item was closed on
/// 2026-05-18 by <c>docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md</c> and now lives
/// as the Tier-1-derived <see cref="F1DepolResidualClosedForm"/> claim on
/// <see cref="F1KnowledgeBase"/>.</para>
///
/// <para>The earlier "non-uniform γ_i: site-dependent dephasing" item was closed on
/// 2026-05-18 by <c>docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md</c> as a NEGATIVE result:
/// the H-block residual scaling factor F(N, G) is γ-independent because the per-site
/// Z-dephasing dissipator kernel is proportional to the identity, so the F1 σ-shift
/// `2Σγ·I` cancels the per-Pauli-string diagonal of `Π·L_D·Π⁻¹ + L_D` exactly for
/// arbitrary γ patterns. The conjectured Σγ_l² scaling does not occur; no formula
/// change to <see cref="PalindromeResidualScalingClaim"/> is required.</para></summary>
public static class F1OpenQuestions
{
    private const string Anchor = "docs/ANALYTICAL_FORMULAS.md F1 \"Breaks for\" clause";

    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "general topology beyond chain/ring/star/K_N",
            "PalindromeResidualScaling is bit-exact verified on chain, ring, star, K_N at N = 4, 5. " +
            "Open: extend the verification (and prove generality) for arbitrary connected graphs, " +
            "including disconnected components and weighted edges.",
            "Sample random graphs at N = 5..7; check the F(N, G) prediction against full L numerical residual; " +
            "if confirmed, the proof is straightforward by graph induction.",
            "experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md"),
    };
}
