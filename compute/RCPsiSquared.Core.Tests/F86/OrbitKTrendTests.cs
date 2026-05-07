using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

public class OrbitKTrendTests
{
    [Fact]
    public void BuildTrend_Endpoint_QPeakInRange_AcrossN5To8()
    {
        var blocks = Enumerable.Range(5, 4)
            .Select(N => new CoherenceBlock(N, n: 1, gammaZero: 0.05))
            .ToList();
        var trend = OrbitKTrend.BuildTrend(orbitIndex: 0, blocks);

        Assert.Equal(4, trend.Trend.Count);
        Assert.All(trend.Trend, e =>
            Assert.True(e.Witness.QPeak > 2.49 && e.Witness.QPeak < 2.55,
                $"N={e.N}: Q_peak={e.Witness.QPeak:F4} outside expected Endpoint range [2.49, 2.55]"));
    }

    [Fact]
    public void BuildTrend_FlankingOrbit_IsEscaping_AtN9_WithDefaultGrid()
    {
        var blocks = new[] { 5, 6, 7, 8, 9 }
            .Select(N => new CoherenceBlock(N, n: 1, gammaZero: 0.05))
            .ToList();
        var trend = OrbitKTrend.BuildTrend(orbitIndex: 1, blocks);

        var defaultGrid = ResonanceScan.DefaultQGrid();
        Assert.True(trend.IsEscaping(defaultGrid),
            $"orbit 1 trend should be escaping at N≥9; last Q_peak = {trend.Trend[^1].Witness.QPeak:F4}");
    }

    [Fact]
    public void BuildTrend_CenterOrbit_QPeakStable_AcrossN5To10()
    {
        var blocks = new[] { 5, 6, 7, 8, 9, 10 }
            .Select(N => new CoherenceBlock(N, n: 1, gammaZero: 0.05))
            .ToList();
        // Center orbit is the largest-Index orbit at each N. Different N has different
        // count; verify each block's center witness Q_peak in [1.4, 1.7].
        foreach (var block in blocks)
        {
            var table = PerF71OrbitKTable.Build(block);
            var center = table.OrbitWitnesses[^1];
            Assert.True(center.QPeak > 1.4 && center.QPeak < 1.7,
                $"N={block.N}: center orbit Q_peak={center.QPeak:F4} outside expected mid-chain range [1.4, 1.7]");
        }
    }

    [Fact]
    public void BuildTrend_ThrowsIfBlocksEmpty()
    {
        Assert.Throws<ArgumentException>(() =>
            OrbitKTrend.BuildTrend(orbitIndex: 0, Enumerable.Empty<CoherenceBlock>()));
    }

    [Fact]
    public void BuildTrend_ThrowsIfOrbitIndexInvalid()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        // N=5 has 2 orbits (indices 0, 1); index 5 is out of range.
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            OrbitKTrend.BuildTrend(orbitIndex: 5, new[] { block }));
    }
}
