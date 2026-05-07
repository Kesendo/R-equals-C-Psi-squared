using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F71-orbit Q_peak trend across an N-range. Tier2Verified.
/// <see cref="BuildTrend"/> iterates the given blocks, builds per-N
/// <see cref="PerF71OrbitKTable"/>, picks the orbit at <paramref name="orbitIndex"/>
/// from each, returns the (N, OrbitKWitness) sequence. <see cref="IsEscaping"/>
/// reads true when the trend's Q_peak is monotonically non-decreasing AND the
/// latest N's witness is grid-escaped — together signal that the orbit's resonance
/// peak has migrated past the scanned Q range.</summary>
public sealed class OrbitKTrend : Claim
{
    public int OrbitIndex { get; }
    public IReadOnlyList<(int N, OrbitKWitness Witness)> Trend { get; }

    public OrbitKTrend(int orbitIndex, IReadOnlyList<(int N, OrbitKWitness Witness)> trend)
        : base($"F71-orbit #{orbitIndex} K-resonance trend across N range",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Statement 2 + Statement 3")
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

        var trend = new List<(int N, OrbitKWitness)>(blockList.Count);
        foreach (var block in blockList)
        {
            var table = PerF71OrbitKTable.Build(block);
            if (orbitIndex < 0 || orbitIndex >= table.OrbitWitnesses.Count)
                throw new ArgumentOutOfRangeException(nameof(orbitIndex),
                    $"orbitIndex {orbitIndex} out of range for N={block.N} ({table.OrbitWitnesses.Count} orbits)");
            trend.Add((block.N, table.OrbitWitnesses[orbitIndex]));
        }
        return new OrbitKTrend(orbitIndex, trend);
    }

    public bool IsEscaping(IReadOnlyList<double> qGrid)
    {
        if (qGrid is null) throw new ArgumentNullException(nameof(qGrid));
        // Monotone non-decreasing check (vacuous for 1-element trend) ...
        for (int i = 1; i < Trend.Count; i++)
        {
            if (Trend[i].Witness.QPeak < Trend[i - 1].Witness.QPeak)
                return false;
        }
        // ... AND the latest witness must be grid-escaped. BuildTrend guarantees the
        // trend is non-empty.
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
            foreach (var (n, witness) in Trend)
                yield return new InspectableNode($"N={n}", summary: witness.Summary);
        }
    }
}
