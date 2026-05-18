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
        // CLOSED 2026-05-18 by docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md.
        // ‖M(T1)‖²_F = 4^(N−1) · [3·Σγ²_T1 + 4·(Σγ_T1)²] derived from the per-site
        // M_l = Π·D_T1·Π⁻¹ + D_T1 kernel: ‖M_l‖²_F = 7, tr(M_l) = −4. Multi-site
        // assembly via tr(M_l† M_l′) = |tr(M_l)|² · 4^(N−2) for l ≠ l′ gives
        // (7 − 4)·Σγ² + 4·(Σγ)². Verified N = 2..5 in simulations/_f1_t1_residual_verify.py.
        // Kept as a registered entry for now so the proof reference stays discoverable
        // via the typed knowledge layer; once promoted to a Tier1Derived claim, remove this.
        new OpenQuestion(
            "T1 amplitude damping: full closed form (CLOSED 2026-05-18)",
            "‖M(T1)‖²_F = 4^(N−1)·[3·Σγ²_T1 + 4·(Σγ_T1)²] in the framework's orthonormal Pauli " +
            "basis. T1 part is H-independent, γ_Z-independent, Frobenius-orthogonal to the H and " +
            "Z-dephasing blocks. The (3, 4) pair derives from the single-site kernel " +
            "‖M_l‖²_F = 7 and |tr(M_l)|² = 16, combined via per-site / identity-elsewhere assembly.",
            "Closed; see docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md (Steps 1-5) for the " +
            "derivation and simulations/_f1_t1_residual_verify.py for bit-exact verification N=2..5.",
            "docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md; " +
            "memory: project_palindrome_frobenius_scaling (T1 extension section, derivation 2026-05-18)"),
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
