using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>The substantive open theoretical items remaining for F86: the "what's missing
/// for full Tier 1 promotion" list from <c>docs/proofs/PROOF_F86_QPEAK.md</c>.</summary>
public static class F86OpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "Item 1': derive 4×4 effective L_eff(Q, b)",
            "Compute cross-coupling ⟨c_α | M_H_per_bond[b] | u_0/v_0⟩ as analytical expressions in (N, n, b); diagonalise; identify which eigenvalue pair gives the Q_peak observed in K_CC_pr (NOT the SVD top pair); derive f_class(x) and HWHM_left/Q_peak as closed forms.",
            "Multi-page algebra; XY single-particle structure of the (n, n+1) block; OBC sine-mode matrix elements.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 1'"),
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
