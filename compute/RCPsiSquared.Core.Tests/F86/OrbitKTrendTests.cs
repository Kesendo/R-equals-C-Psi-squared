using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

public class OrbitKTrendTests
{
    /// <summary>The endpoint orbit's envelope Q_peak along N = 5..8, the values and gate of
    /// <see cref="PerF71OrbitKTableTests.Build_OrbitZero_QPeakAndHwhmRatio_MatchPinnedAnchor"/>.
    /// The trend is not monotone: 2.517 at N = 5, its maximum 2.581 at N = 6, then falling.</summary>
    [Fact]
    public void BuildTrend_Endpoint_QPeak_AcrossN5To8()
    {
        var blocks = Enumerable.Range(5, 4)
            .Select(N => new CoherenceBlock(N, n: 1, gammaZero: 0.05))
            .ToList();
        var trend = OrbitKTrend.BuildTrend(orbitIndex: 0, blocks);
        double[] expected = { 2.517071370613098, 2.580604998577294, 2.560073427805408, 2.543127121122935 };

        Assert.Equal(4, trend.Trend.Count);
        for (int i = 0; i < expected.Length; i++)
            Assert.True(Math.Abs(trend.Trend[i].Witness.QPeak - expected[i]) <= 1e-10,
                $"N={trend.Trend[i].N}: Q_peak={trend.Trend[i].Witness.QPeak:R}, expected {expected[i]:R}");
    }

    [Fact]
    public void BuildTrend_FlankingOrbit_CanonicalPeakRises_ButDoesNotEscape_N5To9()
    {
        // The table reads the canonical (first) maximum: orbit 1 rises 1.483 → 1.637 → 1.699 →
        // 1.727 → 1.749 over N=5..9 and stays inside the default grid. Its far maximum, which
        // does leave the grid from N=7 on, is the global reading of PerF71OrbitObservation.
        var blocks = new[] { 5, 6, 7, 8, 9 }
            .Select(N => new CoherenceBlock(N, n: 1, gammaZero: 0.05))
            .ToList();
        var trend = OrbitKTrend.BuildTrend(orbitIndex: 1, blocks);

        for (int i = 1; i < trend.Trend.Count; i++)
            Assert.True(trend.Trend[i].Witness.QPeak > trend.Trend[i - 1].Witness.QPeak,
                $"N={trend.Trend[i].N}: {trend.Trend[i].Witness.QPeak:R} not above {trend.Trend[i - 1].Witness.QPeak:R}");
        Assert.False(trend.IsEscaping(ResonanceScan.DefaultQGrid()),
            $"last Q_peak = {trend.Trend[^1].Witness.QPeak:F4}");

        // N=5..8 against the frozen orbit table (the scan's resolution, 1e-10).
        var frozen = new PerF71OrbitObservation().Witnesses.Where(w => w.Chromaticity == 2).ToDictionary(w => w.N);
        foreach (var entry in trend.Trend.Where(e => e.N <= 8))
            Assert.Equal(frozen[entry.N].QPeakPerOrbit[1], entry.Witness.QPeak, 1e-10);
    }

    [Fact]
    public void BuildTrend_CenterOrbit_QPeakStable_AcrossN5To10()
    {
        // Center orbit index = N/2 - 1 (orbit count is N/2, center is the largest index).
        // Single-block trend per N exercises BuildTrend without a redundant Build call.
        int[] ns = { 5, 6, 7, 8, 9, 10 };
        var defaultGrid = ResonanceScan.DefaultQGrid();
        foreach (int n in ns)
        {
            var block = new CoherenceBlock(n, n: 1, gammaZero: 0.05);
            int centerIndex = n / 2 - 1;
            var trend = OrbitKTrend.BuildTrend(centerIndex, new[] { block });

            Assert.Single(trend.Trend);
            var center = trend.Trend[0].Witness;
            Assert.True(center.QPeak > 1.4 && center.QPeak < 1.7,
                $"N={n}: center orbit Q_peak={center.QPeak:F4} outside expected mid-chain range [1.4, 1.7]");
            Assert.False(trend.IsEscaping(defaultGrid),
                $"N={n}: center orbit trend should NOT be escaping; Q_peak={center.QPeak:F4}, grid_max={defaultGrid[^1]}");
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
