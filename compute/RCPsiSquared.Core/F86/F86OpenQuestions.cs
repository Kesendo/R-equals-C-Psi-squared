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
            "Item 5: derive the true σ_0(c, N → ∞) asymptote (earlier 2√(2(c−1)) refuted as crossing)",
            "σ_0(c=2) grows monotonically past 2√2: σ_0(N=7..18) = 2.8284, 2.8393, 2.8483, 2.8525, " +
            "2.8561, 2.8579, 2.8596, 2.8604, 2.8612, 2.8616, 2.8620, 2.8623 (the earlier 2√(2(c−1)) " +
            "conjecture was a trajectory crossing at N=7, not an asymptote). Parity-split Aitken on " +
            "even-N and odd-N subsequences converges to σ_0(c=2, N → ∞) ≈ 2.8628 ± 1e-4 — sharpening " +
            "the earlier ~2.85..2.89 band by ~50×. Verified γ-independent (bit-exact across γ ∈ {0.01, " +
            "0.5, 5.0}). Surviving structure: monotone growth within each c; σ_0·√(3/8) bridge to " +
            "g_eff_E (Δ ≤ 0.01 for N ≥ 6); newly noted even/odd-N parity in the convergence rate " +
            "(successive Δ's halve in pairs). Closed-form candidates ruled out at 1e-4 precision: " +
            "2√2 (REJECTED as asymptote, retained as N=7 sweet-spot crossing), √(41/5) ≈ 2.86356, " +
            "√(8 + π/16) ≈ 2.86292.",
            "OBC sine-mode algebra applied to inter-HD-channel coupling; Bogoliubov / JW free-fermion " +
            "transform may make finite-size corrections transparent. The true asymptote likely involves " +
            "an OBC band-edge factor. Polynomial-in-1/N fits are unstable (leading coefficient drifts " +
            "2.85 → 2.91 across degrees 1–5), so simple algebraic-correction ansätze are inadequate " +
            "and the band-edge integral is the natural next step. Higher-c sweet-spot N_c* (where " +
            "σ_0/√(2(c−1)) crosses 2.0) remains unknown: c=3 N=11 ratio is 1.973 (still below 2.0); " +
            "c=4 N=8 ratio is 0.890 (far below). Cross-c ratios are NOT c-independent (1.012, 0.987, " +
            "0.890 at current max N), so the closed form is structurally per-c, not multiplicative.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 5 (open); " +
            "compute/RCPsiSquared.Core.Tests/F86/SigmaZeroAsymptoteReconTests.cs (numerical recon, this session)"),
    };
}
