using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Symmetry;
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
        // Expected formula: (N-1+1)/2 orbits for an N-qubit chain with N-1 bonds —
        // pinned via F71BondOrbitDecomposition's Tier1Derived contract. We verify
        // PerF71OrbitKTable mirrors that count.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        Assert.Equal(expectedOrbits, table.OrbitWitnesses.Count);
    }

    /// <summary>Orbit 0 (the endpoint pair) on the F86 envelope max_t |K_b| across N = 5..10 at
    /// γ₀ = 0.05. The values come from the second route of
    /// <see cref="Build_N6_Orbits_QPeak_MatchTheEnvelopeSecondRoute"/> (the F89 eigendecomposition
    /// operator, an independent t-maximiser, the same Q grid and parabola), averaged over the two
    /// mirror bonds as the table does. On that route the two mirror bonds differ by at most 1.9e-11
    /// in Q_peak through N = 10 (eigensolver rounding, growing with the block); the gate 1e-10 sits
    /// above that and far below any grid or window effect (a 21-point t grid on [3, 8] reads 0.016 to
    /// 0.034 lower in Q_peak and 0.006 to 0.008 lower in the ratio).</summary>
    [Theory]
    [InlineData(5, 2.517071370613098, 0.7765595000530079)]
    [InlineData(6, 2.580604998577294, 0.7815203313098789)]
    [InlineData(7, 2.560073427805408, 0.7813218454967338)]
    [InlineData(8, 2.543127121122935, 0.7807879816314429)]
    [InlineData(9, 2.536319040954015, 0.7807174593467199)]
    [InlineData(10, 2.531819984071979, 0.780704262070788)]
    public void Build_OrbitZero_QPeakAndHwhmRatio_MatchPinnedAnchor(
        int N, double expectedQPeak, double expectedHwhmRatio)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        var w = table.OrbitWitnesses[0];

        Assert.False(w.Orbit.IsSelfPaired);
        Assert.Equal(0, w.Orbit.BondA);
        Assert.Equal(N - 2, w.Orbit.BondB);
        Assert.True(Math.Abs(w.QPeak - expectedQPeak) <= 1e-10,
            $"orbit 0 Q_peak at N={N}: got {w.QPeak:R}, expected {expectedQPeak:R}");
        Assert.True(Math.Abs(w.HwhmLeftOverQPeak - expectedHwhmRatio) <= 1e-10,
            $"orbit 0 HWHM/Q* at N={N}: got {w.HwhmLeftOverQPeak:R}, expected {expectedHwhmRatio:R}");
    }

    /// <summary>The c=2 N=6 orbits on the F86 envelope max_t |K_b|: endpoint {0,4}, the flanking
    /// pair {1,3} above the self-paired central bond {2}. Second route:
    /// simulations/f86_envelope_second_route.py (`6 2 4.0`), the F89 eigendecomposition operator
    /// maximised over t by an independent maximiser (a 1600-point scan on [0.3, 40] and a bounded
    /// Brent step, not the golden-section seed bracket used here), then the same [0.2, 4] grid of 153
    /// points and three-point parabola. Averaged over each mirror pair it reads 2.580604998577294,
    /// 1.6366804166414261, 1.428310837271295. The two routes compute the same
    /// quantity; what separates them is the eigensolvers' rounding of K divided by its curvature over
    /// the grid step and the t-maximisers' tolerance, which enters K only at second order. The gate
    /// 1e-11 sits a decade above that; reading K on a fixed t grid instead of the envelope moves the
    /// value by 1e-3 or more (a 21-point t grid on [3, 8] reads 2.5470, 1.6157, 1.4355).</summary>
    [Fact]
    public void Build_N6_Orbits_QPeak_MatchTheEnvelopeSecondRoute()
    {
        var block = new CoherenceBlock(6, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        var end = table.OrbitWitnesses[0];
        var flank = table.OrbitWitnesses[1];
        var centre = table.OrbitWitnesses[2];
        Assert.Equal(1, flank.Orbit.BondA);
        Assert.Equal(3, flank.Orbit.BondB);
        Assert.True(centre.Orbit.IsSelfPaired);
        Assert.Equal(2, centre.Orbit.BondA);
        Assert.True(Math.Abs(end.QPeak - 2.580604998577294) <= 1e-11, $"endpoint Q_peak {end.QPeak:R}");
        Assert.True(Math.Abs(flank.QPeak - 1.6366804166414261) <= 1e-11, $"flanking Q_peak {flank.QPeak:R}");
        Assert.True(Math.Abs(centre.QPeak - 1.428310837271295) <= 1e-11, $"central Q_peak {centre.QPeak:R}");
    }

    /// <summary>The chain-centre (largest-index) orbit on the F86 envelope at N = 9 (the mirror
    /// pair {3, 4}) and N = 10 (the self-paired bond {4}), from the same second route and gate as
    /// <see cref="Build_OrbitZero_QPeakAndHwhmRatio_MatchPinnedAnchor"/>.</summary>
    [Theory]
    [InlineData(9, 1.5880454575113698, 0.7599841733222457, false)]
    [InlineData(10, 1.6075349142696183, 0.7616106042680529, true)]
    public void Build_CenterOrbit_QPeakAndHwhmRatio_MatchPinnedAnchor(
        int N, double expectedQPeak, double expectedHwhmRatio, bool expectSelfPaired)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        var center = table.OrbitWitnesses[^1];

        Assert.Equal(expectSelfPaired, center.Orbit.IsSelfPaired);
        Assert.True(Math.Abs(center.QPeak - expectedQPeak) <= 1e-10,
            $"center orbit Q_peak at N={N}: got {center.QPeak:R}, expected {expectedQPeak:R}");
        Assert.True(Math.Abs(center.HwhmLeftOverQPeak - expectedHwhmRatio) <= 1e-10,
            $"center orbit HWHM/Q* at N={N}: got {center.HwhmLeftOverQPeak:R}, expected {expectedHwhmRatio:R}");
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
        Assert.All(table.OrbitWitnesses, w =>
        {
            Assert.True(double.IsFinite(w.QPeak) && w.QPeak > 0,
                $"QPeak invalid: {w.QPeak}");
            Assert.True(double.IsFinite(w.HwhmLeft) && w.HwhmLeft > 0,
                $"HwhmLeft invalid: {w.HwhmLeft}");
            Assert.True(double.IsFinite(w.HwhmLeftOverQPeak) && w.HwhmLeftOverQPeak > 0,
                $"HwhmLeftOverQPeak invalid: {w.HwhmLeftOverQPeak}");
            Assert.True(double.IsFinite(w.KMax) && w.KMax > 0,
                $"KMax invalid: {w.KMax}");
        });
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => PerF71OrbitKTable.Build(block));
    }

    [Fact]
    public void IsEscaped_FlankingOrbit_AtN9_CanonicalPeakStaysInsideTheDefaultGrid()
    {
        var block = new CoherenceBlock(N: 9, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        // Orbit 1 at N=9 is the first-flanking pair {1, 6}. Its global maximum lies past the
        // default [0.20, 4.00] grid, but the table reads the canonical first maximum, which
        // lies inside it.
        var orbit1 = table.OrbitWitnesses[1];
        var defaultGrid = ResonanceScan.DefaultQGrid();
        Assert.False(orbit1.IsEscaped(defaultGrid),
            $"orbit 1 at N=9: Q_peak={orbit1.QPeak:F4}, grid_max={defaultGrid[^1]}");
    }
}
