using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>One F71 orbit's Q_peak trajectory across an N-range, Tier2Verified.
/// <see cref="IsEscaping"/> classifies the trajectory as escaping when Q_peak grows
/// monotonically and the latest witness sits at the grid edge.</summary>
public sealed class OrbitKTrend : Claim
{
    public int OrbitIndex { get; }
    public IReadOnlyList<OrbitKTrendEntry> Trend { get; }

    public OrbitKTrend(int orbitIndex, IReadOnlyList<OrbitKTrendEntry> trend)
        : base($"F71-orbit #{orbitIndex} K-resonance trend across N range",
               Tier.Tier2Verified,
               F86Anchors.Statement2Plus3)
    {
        OrbitIndex = orbitIndex;
        Trend = trend;
    }

    public static OrbitKTrend BuildTrend(int orbitIndex, IEnumerable<CoherenceBlock> blocks)
    {
        if (blocks is null) throw new ArgumentNullException(nameof(blocks));
        var blockList = blocks.ToList();
        if (blockList.Count == 0)
            throw new ArgumentException("blocks list must be non-empty", nameof(blocks));

        var trend = new List<OrbitKTrendEntry>(blockList.Count);
        foreach (var block in blockList)
        {
            var table = PerF71OrbitKTable.Build(block);
            if (orbitIndex < 0 || orbitIndex >= table.OrbitWitnesses.Count)
                throw new ArgumentOutOfRangeException(nameof(orbitIndex),
                    $"orbitIndex {orbitIndex} out of range for N={block.N} ({table.OrbitWitnesses.Count} orbits)");
            trend.Add(new OrbitKTrendEntry(block.N, table.OrbitWitnesses[orbitIndex]));
        }
        return new OrbitKTrend(orbitIndex, trend);
    }

    public bool IsEscaping(IReadOnlyList<double> qGrid)
    {
        if (qGrid is null) throw new ArgumentNullException(nameof(qGrid));
        for (int i = 1; i < Trend.Count; i++)
        {
            if (Trend[i].Witness.QPeak < Trend[i - 1].Witness.QPeak)
                return false;
        }
        return Trend[^1].Witness.IsEscaped(qGrid);
    }

    public override string DisplayName =>
        $"OrbitKTrend #{OrbitIndex} (N={Trend[0].N}..{Trend[^1].N})";

    public override string Summary =>
        $"orbit #{OrbitIndex} Q_peak across N={Trend[0].N}..{Trend[^1].N}: " +
        $"{Trend[0].Witness.QPeak:F4} → {Trend[^1].Witness.QPeak:F4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("orbit index", summary: OrbitIndex.ToString());
            yield return new InspectableNode("N range",
                summary: $"{Trend[0].N}..{Trend[^1].N}");
            foreach (var entry in Trend)
                yield return new InspectableNode($"N={entry.N}", summary: entry.Witness.Summary);
        }
    }
}

/// <summary>One entry in an <see cref="OrbitKTrend"/>: the chain length and its
/// orbit witness.</summary>
public sealed record OrbitKTrendEntry(int N, OrbitKWitness Witness);
