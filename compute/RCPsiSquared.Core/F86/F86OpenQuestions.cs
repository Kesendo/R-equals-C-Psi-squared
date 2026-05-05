using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>The substantive open theoretical items remaining for F86: the "what's missing
/// for full Tier 1 promotion" list from <c>docs/proofs/PROOF_F86_QPEAK.md</c>.</summary>
public static class F86OpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "Item 1' (c=2): closed-form HWHM_left/Q_peak constant",
            "Empirical anchor reproduced within tolerance 0.005 (typical residual ≤ 0.001) " +
            "(C2HwhmRatio Tier1Candidate, c=2 N=5..8). " +
            "Directional Endpoint > Interior split derived empirically (gap ≈ 0.022). " +
            "Closed-form constant for HWHM_left/Q_peak ratios per bond class NOT pinned this " +
            "session. Three next directions ranked: (a) first-order perturbation in cross-block " +
            "(most promising — leverages B2's Endpoint < Interior cross-block Frobenius split, " +
            "ε ~ ‖V_b cross‖_F/σ_0 ~ O(0.1)); (b) projector-overlap lift of A3's |u_0⟩, |v_0⟩; " +
            "(c) symbolic char-poly factorisation at Q_EP.",
            "Direction (a): perturb the K-resonance K_b(Q, t) around the Statement-1 2×2 EP form. " +
            "The cross-block contribution at first order in ε should produce the directional split.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 1 (c=2); F86OpenQuestions.cs (this file)"),
        new OpenQuestion(
            "Item 4': extend 4-mode construction to c≥3",
            "Each adjacent-channel pair (HD=2k−1, HD=2k+1) for k ∈ {1, …, c−1} contributes its own (|c_{2k−1}⟩, |c_{2k+1}⟩, |u_0^{(k)}⟩, |v_0^{(k)}⟩) quartet → full effective L is 4·(c−1)-dimensional. Verify slowest pair k=1 still dominates K_b response.",
            "Iterate the 4-mode construction over k; concatenate the orthonormal subspaces; project full block-L; numerical verification at c=3, c=4 first.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 4'"),
        new OpenQuestion(
            "Item 5: derive σ_0 → 2√2 asymptote (c=2)",
            "The trajectory σ_0(N) = {2.7651, 2.8023, 2.8284, 2.8393} at N=5..8 converges to 2√2 = 2.8284. The 2√2 has the look of an XY-chain matrix element √(2/(N+1))·sin(πk·b/(N+1)) but the closed form is not yet derived.",
            "OBC sine-mode algebra applied to the inter-HD-channel coupling matrix; Bogoliubov / JW transform of the SE chain may make this transparent.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 5"),
    };
}
