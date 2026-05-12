using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>F93DetuningAntiPalindromicSpectralInvariance (Tier 1 derived; 2026-05-12):
///
/// <para>For chain XY + uniform Z-dephasing Liouvillian L on N qubits with per-site
/// Z-detuning h_l Z_l added to the Hamiltonian (longitudinal field; preserves joint-popcount),
/// the F71-refined diagonal-block eigenvalue multiset is invariant under any h-distribution
/// satisfying h_l + h_{N-1-l} = 2·h_avg = (2/N)·Σ_l h_l for all l (i.e. h is F71-anti-palindromic
/// around its mean). The full L itself generally changes (F71 broken as L-symmetry), but
/// diagonal-block eigenvalues coincide; breaking lives in eigenvectors.</para>
///
/// <para>Note: only longitudinal h_l Z_l is in scope. Transverse h_l X_l or h_l Y_l would
/// flip popcount and break joint-popcount conservation, taking us out of the BlockSpectrum
/// framework entirely.</para>
///
/// <para>Pi2-Z₄ structure (parameter-side, h-axis): identical to F91 (γ) and F92 (J).
/// h-side detuning twin in the family.</para>
///
/// <para>Algebraic proof + empirical witness (bit-exact at N=4, 5):
/// <c>docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md</c>.</para></summary>
/// <seealso cref="BlockSpectrum.F71AntiPalindromicGammaSpectralInvariance"/>
/// <seealso cref="F92BondAntiPalindromicJSpectralInvariance"/>
public sealed class F93DetuningAntiPalindromicSpectralInvariance : Claim
{
    private readonly BlockSpectrum.JointPopcountSectors _sectors;
    private readonly BlockSpectrum.F71MirrorBlockRefinement _f71;
    private readonly SymmetryFamilyInventory _inventory;

    public F93DetuningAntiPalindromicSpectralInvariance(
        BlockSpectrum.JointPopcountSectors sectors,
        BlockSpectrum.F71MirrorBlockRefinement f71,
        SymmetryFamilyInventory inventory)
        : base("F93DetuningAntiPalindromicSpectralInvariance: chain XY+Z-deph+h_l Z_l L diagonal-block spectrum invariant under h-distributions satisfying h_l+h_{N-1-l}=2·h_avg; h-side twin of F91 γ-Z₄ + F92 J-Z₄.",
               Tier.Tier1Derived,
               "Algebraic mechanism parallel to F91/F92 (diagonal blocks depend only on F71-pair-sums; cross blocks on pair-differences); bit-exact at N=4,5")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
        _inventory = inventory ?? throw new ArgumentNullException(nameof(inventory));
    }

    /// <summary>F71-anti-palindromic deviation norm of a per-site h_l list:
    /// sqrt(Σ_l ((h_l + h_{N-1-l}) − 2·h_avg)²). Zero iff h is F71-anti-palindromic.</summary>
    public static double AntiPalindromicDeviation(IReadOnlyList<double> hPerSite)
    {
        if (hPerSite is null) throw new ArgumentNullException(nameof(hPerSite));
        int N = hPerSite.Count;
        double avg = 0.0;
        for (int l = 0; l < N; l++) avg += hPerSite[l];
        avg /= N;
        double sum = 0.0;
        for (int l = 0; l < N; l++)
        {
            double pair = hPerSite[l] + hPerSite[N - 1 - l];
            double diff = pair - 2.0 * avg;
            sum += diff * diff;
        }
        return Math.Sqrt(sum);
    }

    /// <summary>True iff h_l is F71-anti-palindromic around its mean (within tolerance).</summary>
    public static bool IsAntiPalindromic(IReadOnlyList<double> hPerSite, double tolerance = 1e-10)
        => AntiPalindromicDeviation(hPerSite) < tolerance;

    public override string DisplayName =>
        "F93: F71-anti-palindromic h spectral invariance (Pi2-Z₄ h-side twin of F91/F92)";

    public override string Summary =>
        $"chain XY+Z-deph+h_l Z_l L diagonal-block spectrum invariant under h_l+h_{{N-1-l}}=2·h_avg; algebraic proof + bit-exact N=4,5 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("axis", summary: "parameter (longitudinal Z-detuning h_l)");
            yield return new InspectableNode("scope",
                summary: "longitudinal h_l Z_l only; transverse h_l X_l / h_l Y_l breaks joint-popcount and is out of scope");
            yield return new InspectableNode("Z₄ orbit",
                summary: "anti-palindromic class h_l + h_{N-1-l} = 2·h_avg ∀l (closed under 90°-rotation)");
            yield return new InspectableNode("sister claims",
                summary: "F91 (γ_l, Z-deph), F92 (J_b, bond-coupling)");
            yield return new InspectableNode("anchor proof",
                summary: "docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md");
        }
    }
}
