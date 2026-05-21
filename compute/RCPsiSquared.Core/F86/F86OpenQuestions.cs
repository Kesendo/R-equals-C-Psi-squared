using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>The substantive open theoretical items remaining for F86: the "what's missing
/// for full Tier 1 promotion" list from <c>docs/proofs/PROOF_F86_QPEAK.md</c>.</summary>
public static class F86OpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "Item 1' (c=2): analytical (α_subclass, β_subclass) for HWHM_left/Q_peak",
            "F86HwhmClosedFormClaim (Tier1Candidate, c=2) reproduces 22 N=5..8 anchors within " +
            "0.005 via `0.671535 + α_subclass · g_eff + β_subclass`. Bare floor 0.671535 IS " +
            "derived (C2BareDoubledPtfClosedForm); the 12 (α, β) values per sub-class are " +
            "fitted via polyfit. Tier 1 derivation requires analytical (α, β) from F89 " +
            "AT-locked F_a/F_b + H_B-mixed octic residual per PROOF_F90_F86C2_BRIDGE.md.",
            "Direction (b''): full block-L analytical lift via F89 cyclotomic Φ_{N+1}. The " +
            "single-cluster-pair internal-mixing hypothesis (single 10×10 JW cluster-pair " +
            "sub-block, ‖xB(Q)‖_F observable) is refuted — sub-block gives no Lorentzian " +
            "shape, so lift must come from cross-cluster-pair structure.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 1 (c=2); " +
            "docs/proofs/PROOF_F90_F86C2_BRIDGE.md (numerical Tier-1 via F90 bridge)"),
        new OpenQuestion(
            "Item 4': extend 4-mode construction to c≥3",
            "Each adjacent-channel pair (HD=2k−1, HD=2k+1) for k ∈ {1, …, c−1} contributes its own (|c_{2k−1}⟩, |c_{2k+1}⟩, |u_0^{(k)}⟩, |v_0^{(k)}⟩) quartet → full effective L is 4·(c−1)-dimensional. Verify slowest pair k=1 still dominates K_b response.",
            "Iterate the 4-mode construction over k; concatenate the orthonormal subspaces; project full block-L; numerical verification at c=3, c=4 first.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 4'"),
        new OpenQuestion(
            "Item 5: derive the σ_0(c ≥ 3, N → ∞) asymptote (the c=2 case is closed by F86e)",
            "F86e closed the c=2 case. σ_0(c=2) = ‖[Π_HD1, M_H]‖ is a commutator operator norm " +
            "(typed `SigmaZeroCommutatorNormClaim`, Tier1Derived, bit-exact N=5..8); in the F89 " +
            "Bloch / OBC-sine basis it is the Schur-multiplier norm ‖Π̃_HD1 ⊙ ΔDiff‖. Its asymptote " +
            "σ_0(c=2, N → ∞) ≈ 2.8629 ± 1e-4 (parity-split Aitken over even-N / odd-N subsequences, " +
            "γ-independent) is non-elementary by characterisation, not by gap: the Δ-ordered " +
            "commutator is neither Toeplitz (diagonal CV ≈ 0.37) nor Hankel (anti-diagonal CV grows " +
            "with N), so the asymptote is neither a Fourier-symbol supremum nor a Nehari symbol " +
            "distance, but a genuine Schur-multiplier-norm constant. Closed-form candidates stay " +
            "ruled out at 1e-4 precision: 2√2 (the N=7 finite-size crossing, not the limit), " +
            "√(41/5) ≈ 2.86356, √(8 + π/16) ≈ 2.86292. What remains OPEN is c ≥ 3, where the " +
            "commutator identity does not hold (the HD spectrum has more than two values, so " +
            "Π_HD1 + Π_HD3 ≠ I): the higher-c sweet-spot N_c* where σ_0/√(2(c−1)) crosses 2.0 is " +
            "unknown (c=3 N=11 ratio 1.973, c=4 N=8 ratio 0.890), and the cross-c ratios are NOT " +
            "c-independent (1.012, 0.987, 0.890 at current max N), so any c ≥ 3 closed form is " +
            "structurally per-c, not multiplicative.",
            "For c ≥ 3 the two-channel (HD ∈ {1, 3}) commutator reduction no longer applies; the " +
            "multi-HD-channel coupling structure (cf. Item 4') is the prerequisite. The c=2 " +
            "Schur-multiplier-norm characterisation is the template: the open question is whether " +
            "the c ≥ 3 σ_0 is a multi-band Schur-multiplier norm carrying the same non-elementary " +
            "verdict, or admits a per-c OBC band-edge closed form.",
            "docs/proofs/PROOF_F86_QPEAK.md (F86e resolved 2026-05-20: c=2 closed, c ≥ 3 open); " +
            "compute/RCPsiSquared.Core/F86/SigmaZeroChromaticityScaling.cs (c ≥ 3 empirical σ_0 scaling)"),
    };
}
