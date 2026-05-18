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
/// <see cref="F1KnowledgeBase"/>. The original item asked about graph-parameter (B, D2)
/// dependence; the proof shows depol is per-site only, so the residual scales purely
/// with (Σγ², (Σγ)²) without any B / D2 dependence.</para></summary>
public static class F1OpenQuestions
{
    private const string Anchor = "docs/ANALYTICAL_FORMULAS.md F1 \"Breaks for\" clause";

    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "non-uniform γ_i: site-dependent dephasing",
            "F1 holds for site-dependent γ_i (the identity is per-site additive in the Klein parities). " +
            "Open: characterise the residual scaling when γ_i has site-dependent structure; does the " +
            "main / single-body class scaling factor become Σ_i γ_i² instead of (Σγ)²?",
            "Replace 2Σγ·I → 2Σ_i γ_i·(per-site projector); derive the analogous scaling lemma.",
            Anchor),
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
