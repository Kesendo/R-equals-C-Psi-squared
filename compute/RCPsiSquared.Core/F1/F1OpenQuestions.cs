using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>Open theoretical items for the F1 family: directions where the palindrome
/// identity breaks, partially survives, or generalises. Anchored at the F1 "Breaks for"
/// clause and the OPERATOR_RIGIDITY scaling experiment.</summary>
public static class F1OpenQuestions
{
    private const string Anchor = "docs/ANALYTICAL_FORMULAS.md F1 \"Breaks for\" clause";

    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "depolarizing noise: residual scaling",
            "F1 breaks under depolarizing noise with residual error (2/3)Σγ, linear in γ and N. " +
            "Open: derive the closed-form ‖M_depol‖² scaling and its dependence on graph parameters " +
            "(B, D2), analogous to the main / single-body scaling for non-truly Hamiltonians.",
            "Repeat the bond-sum + spectator-variance + disjoint-supports argument with the depolarizing " +
            "Lindblad operators in place of Z-dephasing.",
            Anchor),
        new OpenQuestion(
            "T1 amplitude damping: full closed form",
            "Memory entry project_palindrome_frobenius_scaling: ‖M‖² = 2^(N+2)·n_YZ·‖H‖²_F + " +
            "4^(N−1)·[3·Σγ_T1² + 4·(Σγ_T1)²]. T1 part is H-independent, γ_Z-independent, M-orthogonal; " +
            "verified N = 3..6. Open: derive this closed form analytically (currently empirical).",
            "Track the T1 dissipator's Π-conjugation directly; the H-independence and γ_Z-independence " +
            "suggest a clean per-site decomposition.",
            "memory: project_palindrome_frobenius_scaling (full closed form for ‖M‖²); " +
            "docs/ANALYTICAL_FORMULAS.md F82 + F84 give the F81-anti-component closed form ‖D_{T1, odd}‖_F = √(Σγ²)·2^(N−1) only"),
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
