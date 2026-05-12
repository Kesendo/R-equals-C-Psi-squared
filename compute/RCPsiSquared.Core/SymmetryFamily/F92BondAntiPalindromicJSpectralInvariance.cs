using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>F92BondAntiPalindromicJSpectralInvariance (Tier 1 derived; 2026-05-12):
///
/// <para>For chain XY + uniform Z-dephasing Liouvillian L on N qubits with inhomogeneous
/// bond couplings J_b (b ∈ {0..N-2}), the F71-refined diagonal-block eigenvalue multiset
/// is invariant under any J-distribution satisfying J_b + J_{N-2-b} = 2·J_avg (i.e. J is
/// F71-anti-palindromic around its mean). The full L generally changes (F71 broken as
/// L-symmetry, off-block-Frobenius nonzero), but diagonal-block eigenvalues coincide;
/// the breaking lives in eigenvectors only.</para>
///
/// <para>Pi2-Z₄ structure on the J-axis: identical to F91 (γ-axis) and F93 (h-axis).
/// Diagonal blocks depend on J only through F71-pair-sums T_b = J_b + J_{N-2-b};
/// cross blocks on pair-differences. See <c>docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md</c>
/// for the full Z₄ table and algebraic derivation.</para>
///
/// <para>Algebraic proof + empirical witness (bit-exact at N=4, 5):
/// <c>docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md</c>.</para></summary>
/// <seealso cref="BlockSpectrum.F71AntiPalindromicGammaSpectralInvariance"/>
/// <seealso cref="F93DetuningAntiPalindromicSpectralInvariance"/>
public sealed class F92BondAntiPalindromicJSpectralInvariance : Claim
{
    private readonly BlockSpectrum.JointPopcountSectors _sectors;
    private readonly BlockSpectrum.F71MirrorBlockRefinement _f71;
    private readonly SymmetryFamilyInventory _inventory;

    public F92BondAntiPalindromicJSpectralInvariance(
        BlockSpectrum.JointPopcountSectors sectors,
        BlockSpectrum.F71MirrorBlockRefinement f71,
        SymmetryFamilyInventory inventory)
        : base("F92BondAntiPalindromicJSpectralInvariance: chain XY+Z-deph L diagonal-block spectrum invariant under J-distributions satisfying J_b+J_{N-2-b}=2·J_avg; J-side twin of F91 γ-Z₄.",
               Tier.Tier1Derived,
               "Algebraic mechanism parallel to F91 PROOF (diagonal blocks depend only on F71-pair-sums T_b = J_b+J_{N-2-b}; cross blocks on pair-differences); bit-exact at N=4,5")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
        _inventory = inventory ?? throw new ArgumentNullException(nameof(inventory));
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
