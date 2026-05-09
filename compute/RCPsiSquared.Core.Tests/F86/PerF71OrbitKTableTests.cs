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

    /// <summary>The chain-center F71 orbit (deepest interior under spatial mirror)
    /// sits at HWHM_left/Q_peak ≈ 1 − 1/4 = 3/4 within tol 0.005 across N=9, 10. The
    /// 1/4 = a_3 = QuarterAsBilinearMaxvalClaim on Pi2DyadicLadder is the mirror
    /// partner of 3/4 (= 1 − a_3) via inversion symmetry a_n · a_{2−n} = 1.
    ///
    /// <para>Two perspectives, same geometric fact:</para>
    /// <list type="bullet">
    ///   <item>Inside view (Q_peak going outward): HWHM_left/Q_peak ≈ 3/4.</item>
    ///   <item>Outside view (Q=0 going inward): half-max-left at Q ≈ Q_peak/4.</item>
    /// </list>
    ///
    /// <para><b>Orbit-escape signature:</b> while the center orbit stays anchored at
    /// 3/4 within ≤ 0.2% across N=9, 10, the first-flanking orbit at N=9 escapes
    /// the default scan grid (see <see cref="IsEscaped_FlagsFlankingOrbit_AtN9_WithDefaultGrid"/>).
    /// The class-mean Interior HwhmLeftOverQPeak drift in
    /// <c>PolarityInheritanceLink._polarityWitnesses</c> across N=5..8 (max ±0.0045
    /// from 3/4) is the visible signature of orbit-mixing: at N≤8 all Interior orbits
    /// are still bounded but flanking orbits are already drifting toward escape; at
    /// N≥9 flanking orbits hit the grid and the class-mean breaks down. The center
    /// orbit's 3/4 ≈ 1 − a_3 plateau is what remains stable.</para>
    ///
    /// <para>Tier 1 candidate. Mirror-partner reading is exact in the structural
    /// sense (HwhmRatio_center = 1 − a_3); empirical drift ≤ 0.2% at N=9, 10 leaves
    /// the closed-form Tier 1 promotion open in F86b Direction (a''-d'') in
    /// <c>C2HwhmRatio.PendingDerivationNote</c>.</para></summary>
    [Theory]
    [InlineData(9)]
    [InlineData(10)]
    public void CenterOrbit_HwhmRatio_AnchorsToQuarterMirrorPartner(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKTable.Build(block);
        var center = table.OrbitWitnesses[^1];

        var ladder = new Pi2DyadicLadderClaim();
        double quarter = ladder.Term(3);          // a_3 = 1/4
        double mirrorPartner = 1.0 - quarter;     // 3/4
        const double tol = 0.005;

        // Inside view: HWHM_left/Q_peak ≈ 1 − a_3 = 3/4
        Assert.True(Math.Abs(center.HwhmLeftOverQPeak - mirrorPartner) <= tol,
            $"Inside view (center orbit) at N={N}: HwhmLeftOverQPeak {center.HwhmLeftOverQPeak:F4} should sit within {tol} of 1 − a_3 = {mirrorPartner}");

        // Outside view: 1 − HWHM_left/Q_peak ≈ a_3 = 1/4
        double outsideView = 1.0 - center.HwhmLeftOverQPeak;
        Assert.True(Math.Abs(outsideView - quarter) <= tol,
            $"Outside view at N={N}: 1 − HwhmLeftOverQPeak {outsideView:F4} should sit within {tol} of a_3 = {quarter}");

        // Distance flip identity: same magnitude, opposite sign
        double distInside = center.HwhmLeftOverQPeak - mirrorPartner;
        double distOutside = outsideView - quarter;
        Assert.Equal(distInside, -distOutside, precision: 14);
    }
}
