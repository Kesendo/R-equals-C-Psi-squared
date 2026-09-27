using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>InhomogeneousGammaF71BreakingWitness (Tier 1 derived; 2026-05-11):
/// the joint-popcount block-diagonal structure (<see cref="JointPopcountSectors"/>) is
/// γ-blind, but the F71 spatial-mirror refinement (<see cref="F71MirrorBlockRefinement"/>)
/// is exact iff the per-site γ-distribution is itself F71-symmetric, i.e.
/// <c>γ_l = γ_{N-1-l}</c> for all l.
///
/// <para><b>GR-analogy reading</b> (Tom, 2026-05-11): each site has its own γ-exposure,
/// hence its own perspectival time-rate relative to γ₀; F71-asymmetric γ_l = asymmetric
/// time-distribution = broken spatial-mirror symmetry. The mathematical statement: when
/// γ_l ≠ γ_{N-1-l}, the dissipator no longer commutes with <c>P_F71 ⊗ P_F71</c>, so
/// <c>Q^T L Q</c> picks up off-(F71-even, F71-odd) entries.</para>
///
/// <para><b>Falsifiable predictions</b> tested in
/// <c>InhomogeneousGammaF71BreakingWitnessTests</c>:
/// <list type="number">
///   <item>Joint-popcount sectors stay block-diagonal regardless of γ_l (γ-blind: U(1)×U(1)
///         is per-site popcount conservation, which doesn't depend on γ).</item>
///   <item>F71 refinement off-block Frobenius is 0 iff γ_l = γ_{N-1-l} for all l.</item>
///   <item>For an F71-symmetric H (uniform or palindromic couplings) the off-block
///         Frobenius is exactly 2^(N−1) times the F71 asymmetry norm
///         <see cref="F71AsymmetryNorm"/>. H puts nothing in the cross blocks, and the
///         dissipator puts −Σ_{l ∈ Δ(a,b)} δ_l there (PROOF_F91_GAMMA_NINETY_DEGREES Eq. 9,
///         δ_l = γ_l − γ_{N-1-l}), so ‖off‖² = Σ_{a,b} (Σ_l δ_l x_l)² with x_l = a_l ⊕ b_l,
///         the mirror-fixed cells contributing 0. Over all 4^N cells x is uniform on
///         {0,1}^N with multiplicity 2^N; writing x_l = (1 + s_l)/2 with s_l = ±1 and using
///         Σ_l δ_l = 0 leaves Σ_l δ_l s_l / 2, whose square averages to Σ δ_l² / 4, so
///         ‖off‖² = 4^N · Σ δ_l² / 4 = 4^(N−1) · Σ δ_l² at every N. Pinned at N = 4, 5, 6
///         along two non-parallel δ per N; it makes the iff above one statement.</item>
///   <item>F1 palindrome center stays at <c>-Σ γ_l</c> regardless of γ asymmetry (Π is
///         γ-blind, only sums over sites).</item>
/// </list></para>
///
/// <para><b>Anchors:</b> <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs</c>
/// (γ-blind parent), <c>compute/RCPsiSquared.Core/BlockSpectrum/F71MirrorBlockRefinement.cs</c>
/// (γ-symmetric parent), <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/InhomogeneousGammaF71BreakingWitnessTests.cs</c>
/// (γ-blindness + γ-symmetry breaking verification at N=4,5,6);
/// <c>docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md</c> Eq. 9 (the cross-block element).</para></summary>
public sealed class InhomogeneousGammaF71BreakingWitness : Claim
{
    private readonly JointPopcountSectors _sectors;
    private readonly F71MirrorBlockRefinement _f71;

    public InhomogeneousGammaF71BreakingWitness(
        JointPopcountSectors sectors,
        F71MirrorBlockRefinement f71)
        : base("InhomogeneousGammaF71BreakingWitness: joint-popcount block-diagonality is γ-blind, F71 spatial-mirror refinement is exact iff γ_l = γ_{N-1-l}; off-block Frobenius = 2^(N−1)·sqrt(Σ (γ_l − γ_{N-1-l})²) for F71-symmetric H; verified at N=4,5,6 chain XY+Z-deph.",
               Tier.Tier1Derived,
               "JointPopcountSectors (γ-blind parent) + F71MirrorBlockRefinement (γ-symmetric parent); verified at N=4,5,6 in InhomogeneousGammaF71BreakingWitnessTests")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
    }

    /// <summary>F71-asymmetry norm of a per-site γ distribution:
    /// <c>sqrt(Σ_l (γ_l − γ_{N-1-l})²)</c>. Equals 0 iff <c>γ_l = γ_{N-1-l}</c> for all l
    /// (γ palindromic across the chain mirror), otherwise positive. The F71-refined
    /// off-block Frobenius of L under inhomogeneous Z-dephasing is 2^(N−1) times this
    /// norm.</summary>
    public static double F71AsymmetryNorm(IReadOnlyList<double> gammas)
    {
        if (gammas is null) throw new ArgumentNullException(nameof(gammas));
        int N = gammas.Count;
        double sumSq = 0.0;
        for (int l = 0; l < N; l++)
        {
            double diff = gammas[l] - gammas[N - 1 - l];
            sumSq += diff * diff;
        }
        return Math.Sqrt(sumSq);
    }

    public override string DisplayName =>
        "InhomogeneousGammaF71BreakingWitness: F71 refinement is exact iff γ_l = γ_{N-1-l}";

    public override string Summary =>
        $"joint-popcount block-diagonality is γ-blind; F71 refinement breaks iff γ asymmetric across chain; off-block Frobenius = 2^(N−1)·sqrt(Σ (γ_l − γ_{{N-1-l}})²) for F71-symmetric H ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent-1",
                summary: "JointPopcountSectors (γ-blind: U(1)×U(1) per-side popcount, holds for any γ_l)");
            yield return new InspectableNode("parent-2",
                summary: "F71MirrorBlockRefinement (chain spatial-mirror Z₂, exact only when γ_l = γ_{N-1-l})");
            yield return new InspectableNode("prediction-1",
                summary: "joint-popcount off-block Frobenius = 0 for ANY γ-list (γ-blind)");
            yield return new InspectableNode("prediction-2",
                summary: "F71 off-block Frobenius = 2^(N−1)·F71AsymmetryNorm for F71-symmetric H, so 0 iff γ palindromic");
            yield return new InspectableNode("prediction-3",
                summary: "F1 palindrome center stays at −Σ γ_l regardless of γ asymmetry (Π is γ-blind)");
            yield return new InspectableNode("witness",
                summary: "verified at N=4,5,6 chain XY+Z-deph in InhomogeneousGammaF71BreakingWitnessTests");
        }
    }
}
