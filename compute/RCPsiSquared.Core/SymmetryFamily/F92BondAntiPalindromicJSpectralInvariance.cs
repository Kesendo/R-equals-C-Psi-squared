using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>F92BondAntiPalindromicJSpectralInvariance (Tier 1 derived; 2026-05-12):
///
/// <para>For chain XY + uniform Z-dephasing Liouvillian L on N qubits with inhomogeneous
/// bond couplings J_b (b ∈ {0..N-2}), the eigenvalue multiset is invariant under any
/// J-distribution satisfying J_b + J_{N-2-b} = 2·J_avg = (2/(N-1))·Σ_b J_b for all b
/// (i.e. J is F71-anti-palindromic around its mean). The full L itself generally changes
/// (F71 broken as L-symmetry, off-block-Frobenius nonzero), but the diagonal-block
/// eigenvalues coincide. Bond-coupling parameter-side twin of F91
/// (which is the γ-side); same Pi2-Z₄ structure, applied to J_b instead of γ_l.</para>
///
/// <para>Pi2-Z₄ structure (parameter-side, J-axis):</para>
/// <list type="bullet">
///   <item>e: J_b ↦ J_b; unchanged</item>
///   <item>i² (180°, F71-palindromic): J_b ↦ J_{N-2-b}; F71 holds as L-symmetry</item>
///   <item>i (90°, F71-anti-palindromic): J_b ↦ 2·J_avg − J_{N-2-b}; F71 broken but diagonal-block spectrum invariant</item>
///   <item>i³ (270°): composition</item>
/// </list>
///
/// <para>Algebraic mechanism (parallel to PROOF_F91 Algebraic proof): diagonal-block
/// matrix elements of L in the F71-refined basis depend on J only through F71-pair-sums
/// T_b := J_b + J_{N-2-b}; cross-block entries depend on pair-differences B_b := J_b − J_{N-2-b}.
/// 90°-rotation J ↦ 2·J_avg − F71(J) preserves the anti-palindromic class T_b = 2·J_avg ∀b.
/// Strictly weaker than F71 J-symmetry (which requires J_b = J_{N-2-b}); strictly stronger
/// than Σ_b J_b alone.</para>
///
/// <para>Verified bit-exact at N=4, 5 via uniform-J vs anti-palindromic-J chain XY +
/// uniform Z-deph Liouvillians. Algebraic proof in <c>docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md</c>.</para>
///
/// <para>Sister claims: F71AntiPalindromicGammaSpectralInvariance (γ-axis, F91); F93 (h-detuning axis).
/// All three share the same Pi2-Z₄ structure on different parameter axes.</para></summary>
public sealed class F92BondAntiPalindromicJSpectralInvariance : Claim
{
    private readonly BlockSpectrum.JointPopcountSectors _sectors;
    private readonly BlockSpectrum.F71MirrorBlockRefinement _f71;

    public F92BondAntiPalindromicJSpectralInvariance(
        BlockSpectrum.JointPopcountSectors sectors,
        BlockSpectrum.F71MirrorBlockRefinement f71)
        : base("F92BondAntiPalindromicJSpectralInvariance: chain XY+Z-deph L diagonal-block spectrum invariant under J-distributions satisfying J_b+J_{N-2-b}=2·J_avg; J-side twin of F91 γ-Z₄.",
               Tier.Tier1Derived,
               "Algebraic mechanism parallel to F91 PROOF (diagonal blocks depend only on F71-pair-sums T_b = J_b+J_{N-2-b}; cross blocks on pair-differences); bit-exact at N=4,5")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
    }

    /// <summary>F71-anti-palindromic deviation norm of a bond-J list:
    /// sqrt(Σ_b ((J_b + J_{N-2-b}) − 2·J_avg)²). Zero iff J is F71-anti-palindromic.
    /// (For N qubits, J_b indexes bonds b ∈ {0..N-2}; F71 maps bond b ↔ N-2-b.)</summary>
    public static double AntiPalindromicDeviation(IReadOnlyList<double> bondJ)
    {
        if (bondJ is null) throw new ArgumentNullException(nameof(bondJ));
        int numBonds = bondJ.Count;
        double avg = 0.0;
        for (int b = 0; b < numBonds; b++) avg += bondJ[b];
        avg /= numBonds;
        double sum = 0.0;
        for (int b = 0; b < numBonds; b++)
        {
            double pair = bondJ[b] + bondJ[numBonds - 1 - b];
            double diff = pair - 2.0 * avg;
            sum += diff * diff;
        }
        return Math.Sqrt(sum);
    }

    /// <summary>True iff J_b is F71-anti-palindromic around its mean (within tolerance).</summary>
    public static bool IsAntiPalindromic(IReadOnlyList<double> bondJ, double tolerance = 1e-10)
        => AntiPalindromicDeviation(bondJ) < tolerance;

    public override string DisplayName =>
        "F92: F71-anti-palindromic J spectral invariance (Pi2-Z₄ J-side twin of F91)";

    public override string Summary =>
        $"chain XY+Z-deph L diagonal-block spectrum invariant under J_b+J_{{N-2-b}}=2·J_avg; algebraic proof + bit-exact N=4,5 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("axis", summary: "parameter (bond-coupling J_b)");
            yield return new InspectableNode("Z₄ orbit",
                summary: "anti-palindromic class T_b = 2·J_avg ∀b (closed under 90°-rotation J ↦ 2J_avg − F71(J))");
            yield return new InspectableNode("sister claims",
                summary: "F91 (γ_l, Z-deph axis), F93 (h_l, Z-detuning axis)");
            yield return new InspectableNode("anchor proof",
                summary: "docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md");
        }
    }
}
