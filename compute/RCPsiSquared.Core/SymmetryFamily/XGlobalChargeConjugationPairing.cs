using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>XGlobalChargeConjugationPairing (Tier 1 derived; 2026-05-12): the global
/// X-string operator X⊗N = ⊗_l X_l flips every bit on each Hilbert side, mapping a
/// computational basis state |a⟩ to |~a⟩ = |2^N - 1 - a⟩. On joint-popcount labels this
/// pairs sector (p_c, p_r) with sector (N - p_c, N - p_r). Under chain XY+Z-deph L,
/// X⊗N commutes with L, so paired sectors share spectra: full eig over one sector
/// gives the spectrum of its X⊗N-pair for free.
///
/// <para>Algorithmic gain: halves the number of distinct eigendecompositions needed at
/// any N. Block sizes are unchanged; this is sector-pairing, not sector-splitting.</para></summary>
public sealed class XGlobalChargeConjugationPairing : Claim
{
    private readonly SymmetryFamilyInventory _inventory;

    public XGlobalChargeConjugationPairing(SymmetryFamilyInventory inventory)
        : base("XGlobalChargeConjugationPairing: X⊗N pairs sector (p_c, p_r) with (N-p_c, N-p_r); paired sectors share spectra (chain XY+Z-deph commutes with X⊗N).",
               Tier.Tier1Derived,
               "X⊗N · σ_α · X⊗N = (-1)^{n_Y+n_Z}·σ_α; chain XY commutes; Z-deph commutes (since (-Z)·(-Z) = +Z·Z); joint-popcount labels reflect bit-count, so X⊗N flips both popcount labels")
    {
        _inventory = inventory ?? throw new ArgumentNullException(nameof(inventory));
    }

    /// <summary>X⊗N image of joint-popcount sector (p_c, p_r): (N - p_c, N - p_r).</summary>
    public static (int PairCol, int PairRow) PairSector(int N, int pCol, int pRow)
    {
        if (N < 0 || pCol < 0 || pCol > N || pRow < 0 || pRow > N)
            throw new ArgumentOutOfRangeException(
                $"N={N}, pCol={pCol}, pRow={pRow}: require 0 ≤ pCol, pRow ≤ N");
        return (N - pCol, N - pRow);
    }

    /// <summary>True if sector (p_c, p_r) is X⊗N-self-paired (paired with itself).
    /// At even N: (p_c, p_r) = (N/2, N/2) is the only self-paired sector. At odd N: never.</summary>
    public static bool IsSelfPaired(int N, int pCol, int pRow)
    {
        var (pairC, pairR) = PairSector(N, pCol, pRow);
        return pairC == pCol && pairR == pRow;
    }

    /// <summary>Number of distinct spectral classes after X⊗N pairing at given N.
    /// Self-paired sectors count once; non-self-paired sectors count once per pair.</summary>
    public static int DistinctSpectralClasses(int N)
    {
        int total = (N + 1) * (N + 1);
        int selfPaired = (N % 2 == 0) ? 1 : 0;
        return (total + selfPaired) / 2;
    }

    public override string DisplayName =>
        "XGlobalChargeConjugationPairing: X⊗N pairs (p_c, p_r) ↔ (N-p_c, N-p_r); halves number of eig-calls";

    public override string Summary =>
        $"X⊗N sector-pairing under chain XY+Z-deph; (N+1)² sectors collapse to ≈ (N+1)²/2 distinct spectral classes ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("pair formula",
                summary: "(p_c, p_r) ↔ (N - p_c, N - p_r)");
            yield return new InspectableNode("N=8 distinct classes",
                summary: $"{DistinctSpectralClasses(8)} (vs 81 unpaired sectors)");
            yield return new InspectableNode("N=10 distinct classes",
                summary: $"{DistinctSpectralClasses(10)} (vs 121 unpaired sectors)");
        }
    }
}
