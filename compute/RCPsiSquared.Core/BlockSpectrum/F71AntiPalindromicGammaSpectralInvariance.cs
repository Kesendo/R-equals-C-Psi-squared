using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>F71AntiPalindromicGammaSpectralInvariance (Tier 1 candidate; 2026-05-11):
/// for chain XY + Z-dephasing Liouvillian on N qubits, the <b>F71-refined diagonal-block
/// spectrum</b> (the multiset of eigenvalues of the (F71-even, F71-odd) diagonal sub-blocks
/// of <c>Q^T L Q</c>) is invariant under any γ-distribution satisfying
/// <c>γ_l + γ_{N-1-l} = 2·γ_avg</c> for all l ∈ {0..N−1},
/// where <c>γ_avg = (1/N)·Σ γ_l</c>. We call such γ-distributions <i>F71-anti-palindromic
/// around their mean</i>.
///
/// <para><b>Important caveat — this is NOT full-L spectral invariance.</b> Under
/// anti-palindromic γ, the F71-rotated Liouvillian <c>Q^T L Q</c> is generally NOT
/// block-diagonal (the (F71-even ↔ F71-odd) cross blocks carry the F71-asymmetry, see
/// <see cref="InhomogeneousGammaF71BreakingWitness"/>). Diagonalising only the diagonal
/// sub-blocks (which is what <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>
/// does) gives the same multiset across all anti-palindromic γ-profiles, but this multiset
/// is NOT the true full-L spectrum when γ breaks F71 symmetry. The full L spectrum
/// (computed via <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> or direct
/// diagonalisation of <c>L</c>) differs across anti-palindromic γ-profiles by the
/// perturbation induced by the cross blocks.</para>
///
/// <para><b>What anti-palindromy controls.</b> The diagonal blocks of <c>Q^T L Q</c> in
/// the F71-refined basis depend on γ only through the <i>per-pair sums</i>
/// <c>γ_l + γ_{N-1-l}</c>. When all pair sums equal <c>2·γ_avg</c>, the diagonal blocks
/// reduce to the same form as under uniform γ = γ_avg. The off-diagonal F71-even ↔ F71-odd
/// coupling depends on per-pair <i>differences</i> <c>γ_l − γ_{N-1-l}</c> (see
/// <see cref="InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm"/>), independent of
/// the pair-sum invariance.</para>
///
/// <para>This is a <b>strict relaxation</b> of the trivial uniform-γ case (γ_l = γ_avg)
/// and a <b>distinct condition</b> from F71 symmetry (γ_l = γ_{N-1-l}, captured by
/// <see cref="F71MirrorBlockRefinement"/>): anti-palindromic and F71-palindromic
/// γ-profiles overlap only at uniform γ. Both are <b>strictly stronger</b> than F1 alone
/// (which only requires Σγ_l invariant). For odd N, the middle site
/// <c>l = (N-1)/2</c> must equal γ_avg.</para>
///
/// <para><b>Empirical witness at N=4, 5, 6.</b> The F71-refined diagonal-sub-block
/// eigenvalue multiset is bit-identical (binned at 1e-7) across:
/// <list type="bullet">
///   <item>uniform γ = 0.45,</item>
///   <item>monotonic linear anti-palindromic γ (e.g. [0.2, 0.3, 0.4, 0.5, 0.6, 0.7] at N=6),</item>
///   <item>non-monotonic anti-palindromic γ (e.g. [0.3, 0.5, 0.4, 0.5, 0.4, 0.6] at N=6).</item>
/// </list>
/// The <i>full</i> L spectrum (joint-popcount block-eig, no F71 truncation) DIFFERS
/// across these profiles, confirming that the invariance is a property of the F71-projected
/// dynamics only. Permuted non-anti-palindromic γ
/// ([0.7, 0.2, 0.5, 0.3, 0.6, 0.4]) and concentrated γ ([0.1, 0.1, 0.1, 0.1, 0.1, 2.2])
/// produce distinct F71-refined diagonal spectra.</para>
///
/// <para>Tier outcome <see cref="Tier.Tier1Candidate"/> rather than Derived: the
/// closed-form proof that the F71-refined diagonal blocks reduce to a γ_avg-only form when
/// γ is anti-palindromic is plausible but not yet written. The empirical witness is strong
/// (bit-exact across multiple anti-palindromic profiles at N=4, 5, 6).</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs</c>
/// (γ-blind parent), <c>compute/RCPsiSquared.Core/BlockSpectrum/F71MirrorBlockRefinement.cs</c>
/// (γ-symmetric parent — anti-palindromy is a distinct condition that overlaps only at uniform γ),
/// <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/F71AntiPalindromicGammaSpectralInvarianceTests.cs</c>
/// (empirical witness at N=4, 5, 6 + explicit verification that full-L spectrum DIFFERS).</para></summary>
public sealed class F71AntiPalindromicGammaSpectralInvariance : Claim
{
    private readonly JointPopcountSectors _sectors;
    private readonly F71MirrorBlockRefinement _f71;

    public F71AntiPalindromicGammaSpectralInvariance(
        JointPopcountSectors sectors,
        F71MirrorBlockRefinement f71)
        : base("F71AntiPalindromicGammaSpectralInvariance: chain XY+Z-deph F71-refined diagonal-block eigenvalue multiset is invariant under any γ-distribution satisfying γ_l + γ_{N-1-l} = 2·γ_avg; full L spectrum differs (cross blocks carry the asymmetry); diagonal blocks depend only on per-pair sums = 2·γ_avg under anti-palindromy; bit-exact verified at N=4,5,6.",
               Tier.Tier1Candidate,
               "JointPopcountSectors (γ-blind parent) + F71MirrorBlockRefinement (γ-symmetric parent; anti-palindromy distinct, both reduce to uniform-γ at intersection); algebraic proof that F71-rotated diagonal blocks depend only on per-pair sums is plausible but not yet written; empirical bit-exact witness at N=4,5,6 in F71AntiPalindromicGammaSpectralInvarianceTests")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
    }

    /// <summary>F71-anti-palindromic deviation norm of a per-site γ distribution:
    /// <c>sqrt(Σ_l ((γ_l + γ_{N-1-l}) − 2·γ_avg)²)</c> where
    /// <c>γ_avg = (1/N)·Σ_l γ_l</c>. Equals 0 iff every F71-pair sum equals
    /// <c>2·γ_avg</c> (γ-distribution F71-anti-palindromic around its mean),
    /// otherwise positive. The empirical scaling proxy for spectral deviation of the
    /// F71-refined diagonal sub-blocks from the uniform-γ_avg reference.</summary>
    public static double AntiPalindromicDeviation(IReadOnlyList<double> gammas)
    {
        if (gammas is null) throw new ArgumentNullException(nameof(gammas));
        int N = gammas.Count;
        if (N == 0) return 0.0;
        double avg = 0.0;
        for (int l = 0; l < N; l++) avg += gammas[l];
        avg /= N;
        double sumSq = 0.0;
        for (int l = 0; l < N; l++)
        {
            double pair = gammas[l] + gammas[N - 1 - l];
            double diff = pair - 2.0 * avg;
            sumSq += diff * diff;
        }
        return Math.Sqrt(sumSq);
    }

    /// <summary>True iff γ is F71-anti-palindromic around its mean within
    /// <paramref name="tolerance"/>. For odd N this implicitly requires
    /// <c>γ_{(N-1)/2} = γ_avg</c> (middle site self-paired).</summary>
    public static bool IsAntiPalindromic(IReadOnlyList<double> gammas, double tolerance = 1e-10)
        => AntiPalindromicDeviation(gammas) < tolerance;

    public override string DisplayName =>
        "F71AntiPalindromicGammaSpectralInvariance: F71-refined diagonal-block spectrum invariant under γ_l + γ_{N-1-l} = 2·γ_avg";

    public override string Summary =>
        $"chain XY+Z-deph F71-refined diagonal-block spectrum invariant whenever γ is anti-palindromic around its mean; full L spectrum differs (cross blocks carry the asymmetry); bit-exact at N=4,5,6 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent-1",
                summary: "JointPopcountSectors (γ-blind: U(1)×U(1) per-side popcount conservation, holds for any γ_l)");
            yield return new InspectableNode("parent-2",
                summary: "F71MirrorBlockRefinement (γ-symmetric: exact iff γ_l = γ_{N-1-l}; distinct from anti-palindromy, overlap only at uniform γ)");
            yield return new InspectableNode("condition",
                summary: "γ_l + γ_{N-1-l} = 2·γ_avg for all l (distinct from F71 symmetry; strictly stronger than F1 alone)");
            yield return new InspectableNode("scope",
                summary: "applies to F71-refined DIAGONAL sub-block spectrum only; full L spectrum differs because (F71-even ↔ F71-odd) cross blocks are nonzero under F71-asymmetric γ");
            yield return new InspectableNode("prediction-1",
                summary: "anti-palindromic γ → F71-refined diagonal sub-block spectrum = uniform-γ_avg sub-block spectrum (bit-identical)");
            yield return new InspectableNode("prediction-2",
                summary: "non-anti-palindromic γ with same Σγ → F71-refined diagonal sub-block spectrum DIFFERS (F1-only protection insufficient)");
            yield return new InspectableNode("prediction-3",
                summary: "odd N: middle site γ_{(N-1)/2} must equal γ_avg (forced by self-pairing)");
            yield return new InspectableNode("prediction-4",
                summary: "full L spectrum DIFFERS across anti-palindromic γ-profiles (cross-block perturbation lifts the diagonal-block degeneracy)");
            yield return new InspectableNode("witness",
                summary: "bit-exact F71-refined diagonal-block spectrum equality across anti-palindromic γ-profiles at N=4,5,6 chain XY+Z-deph");
            yield return new InspectableNode("open",
                summary: "closed-form proof that F71-rotated diagonal blocks depend on γ only through per-pair sums (γ_l + γ_{N-1-l}) not yet written");
        }
    }
}
