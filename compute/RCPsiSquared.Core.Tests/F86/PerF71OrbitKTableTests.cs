using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

public class PerF71OrbitKTableTests
{
    [Theory]
    [InlineData(5, 2)]
    [InlineData(6, 3)]
    [InlineData(7, 3)]
    [InlineData(8, 4)]
    [InlineData(9, 4)]
    [InlineData(10, 5)]
    public void Build_OrbitCount_MatchesF71BondOrbitDecomposition(int N, int expectedOrbits)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        Assert.Equal(expectedOrbits, table.OrbitWitnesses.Count);
        Assert.Equal(expectedOrbits, new F71BondOrbitDecomposition(N).Orbits.Count);
    }

    /// <summary>Pinned anchor for orbit 0 (Endpoint) Q_peak / HWHM/Q* across N=5..10
    /// at γ₀=0.05. Sources: PolarityInheritanceLink._polarityWitnesses for N=5..8
    /// (live-pipeline pinned 2026-05-06), 2026-05-07 c2hwhm CLI runs for N=9, N=10.
    /// Tolerance 0.005 matches existing C2HwhmRatioTests.</summary>
    [Theory]
    [InlineData(5, 2.5008, 0.7700)]
    [InlineData(6, 2.5470, 0.7738)]
    [InlineData(7, 2.5299, 0.7738)]
    [InlineData(8, 2.5145, 0.7734)]
    [InlineData(9, 2.5082, 0.7733)]
    [InlineData(10, 2.5039, 0.7733)]
    public void Build_OrbitZero_QPeakAndHwhmRatio_MatchPinnedAnchor(
        int N, double expectedQPeak, double expectedHwhmRatio)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        var w = table.OrbitWitnesses[0];

        Assert.False(w.Orbit.IsSelfPaired);
        Assert.Equal(0, w.Orbit.BondA);
        Assert.Equal(N - 2, w.Orbit.BondB);
        Assert.True(Math.Abs(w.QPeak - expectedQPeak) <= 0.005,
            $"orbit 0 Q_peak at N={N}: got {w.QPeak:F4}, expected {expectedQPeak:F4} ± 0.005");
        Assert.True(Math.Abs(w.HwhmLeftOverQPeak - expectedHwhmRatio) <= 0.005,
            $"orbit 0 HWHM/Q* at N={N}: got {w.HwhmLeftOverQPeak:F4}, expected {expectedHwhmRatio:F4} ± 0.005");
    }

    /// <summary>Pinned anchor for the chain-center (largest-Index) orbit at N=9, N=10
    /// from the 2026-05-07 c2hwhm CLI runs. N=9: orbit 3 = mirror pair {3,4}; N=10:
    /// orbit 4 = self-paired {4}.</summary>
    [Theory]
    [InlineData(9, 1.5655, 0.7503, false)]
    [InlineData(10, 1.5829, 0.7518, true)]
    public void Build_CenterOrbit_QPeakAndHwhmRatio_MatchPinnedAnchor(
        int N, double expectedQPeak, double expectedHwhmRatio, bool expectSelfPaired)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        var center = table.OrbitWitnesses[^1];

        Assert.Equal(expectSelfPaired, center.Orbit.IsSelfPaired);
        Assert.True(Math.Abs(center.QPeak - expectedQPeak) <= 0.005,
            $"center orbit Q_peak at N={N}: got {center.QPeak:F4}, expected {expectedQPeak:F4}");
        Assert.True(Math.Abs(center.HwhmLeftOverQPeak - expectedHwhmRatio) <= 0.005,
            $"center orbit HWHM/Q* at N={N}: got {center.HwhmLeftOverQPeak:F4}, expected {expectedHwhmRatio:F4}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(9)]
    public void Build_AllWitnessValues_AreFiniteAndPositive(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        Assert.NotEmpty(table.OrbitWitnesses);
        Assert.All(table.OrbitWitnesses,
            w => Assert.True(double.IsFinite(w.QPeak) && w.QPeak > 0));
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => PerF71OrbitKTable.Build(block));
    }

    [Fact]
    public void IsEscaped_FlagsFlankingOrbit_AtN9_WithDefaultGrid()
    {
        var block = new CoherenceBlock(N: 9, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        // Orbit 1 at N=9 is the first-flanking pair {1, 6} which escapes the default
        // [0.20, 4.00] grid (Q_peak hits grid edge 4.0).
        var orbit1 = table.OrbitWitnesses[1];
        var defaultGrid = ResonanceScan.DefaultQGrid();
        Assert.True(orbit1.IsEscaped(defaultGrid),
            $"orbit 1 at N=9 should be escaped; got Q_peak={orbit1.QPeak:F4}, grid_max={defaultGrid[^1]}");
    }
}
