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
/// any N. Block sizes are unchanged; this is sector-pairing, not sector-splitting.</para>
///
/// <para>Empirical anchor (April 2026, retroactively explained): <c>experiments/DEGENERACY_HUNT.md</c>
/// observed at N=5 that the 14 degenerate eigenmodes at Re(λ)=−0.400 spread across joint-popcount
/// sectors (1,1)↔(4,4), (2,2)↔(3,3), (1,3)↔(4,2). These pairings are X⊗N images by the rule above;
/// the April observation is empirical confirmation of this primitive one month before it was typed.</para></summary>
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

    /// <summary>Partition a list of joint-popcount sectors into "primary" sectors (compute eig)
    /// and "follower" sectors (copy spectrum from X⊗N pair). Self-paired sectors are always
    /// primary; non-self-paired pairs use lex-smaller (PCol, PRow) as primary.
    ///
    /// <para>Returns indices into the input sectors list. The dictionary maps follower-index
    /// → primary-index. Used by LiouvillianBlockSpectrum and F71MirrorBlockRefinement to
    /// avoid duplicate eig computations across X⊗N-paired sectors.</para>
    ///
    /// <para>Optional <paramref name="sectorSize"/> projection: when supplied, primary sectors
    /// are sorted descending by size so the most expensive eig starts first under
    /// <see cref="Parallel.ForEach"/>; this overlaps the largest block's wall-time with smaller
    /// blocks running concurrently.</para></summary>
    public static (List<int> Primaries, Dictionary<int, int> FollowerToPrimary) PartitionByXNPairing<TSector>(
        int N,
        IReadOnlyList<TSector> sectors,
        Func<TSector, (int PCol, int PRow)> getPair,
        Func<TSector, int>? sectorSize = null)
    {
        if (sectors is null) throw new ArgumentNullException(nameof(sectors));
        if (getPair is null) throw new ArgumentNullException(nameof(getPair));

        var sectorIndexByPair = new Dictionary<(int, int), int>(sectors.Count);
        for (int i = 0; i < sectors.Count; i++)
            sectorIndexByPair[getPair(sectors[i])] = i;

        var primaries = new List<int>(sectors.Count);
        var followerToPrimary = new Dictionary<int, int>();
        for (int i = 0; i < sectors.Count; i++)
        {
            var (pCol, pRow) = getPair(sectors[i]);
            if (IsSelfPaired(N, pCol, pRow))
            {
                primaries.Add(i);
                continue;
            }
            var (pairCol, pairRow) = PairSector(N, pCol, pRow);
            bool isPrimary = pCol < pairCol || (pCol == pairCol && pRow < pairRow);
            if (isPrimary) primaries.Add(i);
            else followerToPrimary[i] = sectorIndexByPair[(pairCol, pairRow)];
        }

        if (sectorSize is not null)
            primaries.Sort((a, b) => sectorSize(sectors[b]).CompareTo(sectorSize(sectors[a])));

        return (primaries, followerToPrimary);
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
