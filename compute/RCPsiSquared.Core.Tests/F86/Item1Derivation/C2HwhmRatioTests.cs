using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), Stage D2: HWHM_left/Q_peak ratio empirical anchor +
/// bond-class signature pin. The empirical-anchor test is the load-bearing one — the closed
/// form (analytical pass within the 90-min time-box) is registered via the Tier on the
/// returned Claim.</summary>
public class C2HwhmRatioTests
{
    /// <summary>Class-averaged HWHM_left/Q_peak on the F86 envelope max_t |K_b| (bond-class curves
    /// averaged first, then the peak reading), c = 2, γ₀ = 0.05. The values come from an independent
    /// route: K_b from the F89 eigendecomposition operator, maximised over t by a dense scan plus
    /// bounded Brent, then the same Q grid and parabola. On that route the mirror bonds of a class
    /// agree to about 1e-11 through N = 8; the gate 1e-10 sits above that. The endpoint class holds
    /// near 0.781 from N = 6 on; the interior class rises with N.</summary>
    [Theory]
    [InlineData(5, BondClass.Interior, 0.7523433331321849)]
    [InlineData(5, BondClass.Endpoint, 0.7765595000526813)]
    [InlineData(6, BondClass.Interior, 0.7569384805030636)]
    [InlineData(6, BondClass.Endpoint, 0.7815203313098792)]
    [InlineData(7, BondClass.Interior, 0.7591357778625846)]
    [InlineData(7, BondClass.Endpoint, 0.7813218454965299)]
    [InlineData(8, BondClass.Interior, 0.7601773594897803)]
    [InlineData(8, BondClass.Endpoint, 0.7807879816311578)]
    public void HwhmRatio_C2_MatchesEnvelopeAnchor(int N, BondClass cls, double expected)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        double r = ratio.HwhmLeftOverQPeakMean(cls);
        Assert.True(Math.Abs(r - expected) <= 1e-10,
            $"HWHM_left/Q_peak at N={N}, {cls}: got {r:R}, expected {expected:R}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Witnesses_AllBondsTagged(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        Assert.Equal(block.NumBonds, ratio.Witnesses.Count);
        Assert.Equal(BondClass.Endpoint, ratio.Witnesses[0].BondClass);
        Assert.Equal(BondClass.Endpoint, ratio.Witnesses[block.NumBonds - 1].BondClass);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void HwhmRatio_EndpointGreaterThanInterior_AtCEquals2(int N)
    {
        // On the envelope at c=2 N=5..8 the endpoint class exceeds the interior class by 0.021 to 0.025.
        // This is the bond-class signature.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        double interior = ratio.HwhmLeftOverQPeakMean(BondClass.Interior);
        double endpoint = ratio.HwhmLeftOverQPeakMean(BondClass.Endpoint);
        Assert.True(endpoint > interior,
            $"Endpoint ({endpoint:F4}) should exceed Interior ({interior:F4}) at c=2 N={N}");
    }

    [Fact]
    public void Tier_IsValidOutcome()
    {
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        Assert.True(ratio.Tier == Tier.Tier1Derived
                    || ratio.Tier == Tier.Tier1Candidate
                    || ratio.Tier == Tier.Tier2Verified);
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2HwhmRatio.Build(block));
    }

    [Fact]
    public void Build_ThrowsGridEdgeEscape_WhenTheGridEndsBeforeTheFirstMaximum()
    {
        // Strict mode raises GridEdgeEscapeException when a bond's Q_peak sits within one dQ
        // of the grid's upper edge, so the caller cannot consume a grid-snap reading. The
        // per-bond witness is the canonical peak (first local maximum); a grid that ends at
        // Q = 1.2 lies below every bond's first maximum at N=5 (1.48 and 2.52), so each bond
        // rises to the edge and the first one in the scan, bond 0, is flagged.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var shortGrid = ResonanceScan.LinearQGrid(0.20, 1.20, 41);
        var ex = Assert.Throws<GridEdgeEscapeException>(
            () => C2HwhmRatio.Build(block, qGrid: shortGrid, throwOnGridEdgeSnap: true));
        Assert.Equal(0, ex.Bond);
        Assert.Equal(1.2, ex.GridUpper, precision: 12);
        Assert.True(ex.QPeak >= ex.GridUpper - ex.DQ,
            $"Q_peak {ex.QPeak} should be within one dQ {ex.DQ} of grid upper {ex.GridUpper}");
    }

    [Fact]
    public void Build_AtN9_DefaultGrid_CanonicalPeaksStayInside()
    {
        // The flanking pair b=1 / b=6 at N=9 has a far maximum past Q = 4 (the global
        // reading snaps it to the edge), but its canonical first maximum lies inside the
        // default grid, so strict mode does not throw.
        var block = new CoherenceBlock(N: 9, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block, qGrid: null, throwOnGridEdgeSnap: true);
        Assert.All(ratio.Witnesses, w => Assert.True(w.QPeak < 4.0 - 0.025, $"b={w.Bond}: {w.QPeak}"));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void ComputePerBond_ReturnsFiniteAndPositiveValues(int N)
    {
        // Per-bond pipeline must yield finite Q_peak, K_max, HWHM_left for every bond.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        for (int b = 0; b < block.NumBonds; b++)
        {
            var (qPeak, kMax, hwhmLeft) = ratio.ComputePerBond(b);
            Assert.True(double.IsFinite(qPeak) && qPeak > 0, $"b={b}: Q_peak={qPeak}");
            Assert.True(double.IsFinite(kMax) && kMax > 0, $"b={b}: K_max={kMax}");
            Assert.True(double.IsFinite(hwhmLeft) && hwhmLeft > 0, $"b={b}: HWHM_left={hwhmLeft}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void Witnesses_HaveConsistentRatio(int N)
    {
        // Each witness's HwhmLeftOverQPeak must equal HwhmLeft / QPeak.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        foreach (var w in ratio.Witnesses)
        {
            double recomputed = w.HwhmLeft / w.QPeak;
            Assert.Equal(recomputed, w.HwhmLeftOverQPeak, 12);
        }
    }

    [Fact]
    public void BareDoubledPtfConstants_AreUniversal()
    {
        // Tier 1 derived analytical constants from the 2026-05-06 doubled-PTF Ansatz.
        // Don't depend on N, n, or any block parameter. Pinned to 6-digit derivation precision.
        Assert.InRange(C2HwhmRatio.BareDoubledPtfXPeak, 2.196909, 2.196911);
        Assert.InRange(C2HwhmRatio.BareDoubledPtfHwhmRatio, 0.671534, 0.671536);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EmpiricalHwhmRatios_AboveBareDoubledPtfFloor_ProbeBlockContribution(int N)
    {
        // The doubled-PTF Ansatz gives the SVD-block-only floor; on the envelope, the observable
        // the class values are read on, it is C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak
        // = 0.7265 (0.671535 is the fixed-τ = 1 slice). Both bond classes sit above it, by
        // +0.026 to +0.036 (Interior) and +0.050 to +0.055 (Endpoint) over N = 5..12.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        var interior = ratio.HwhmLeftOverQPeakMean(BondClass.Interior);
        var endpoint = ratio.HwhmLeftOverQPeakMean(BondClass.Endpoint);
        double floor = C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak;
        Assert.True(interior > floor, $"Interior {interior:F4} must exceed the bare envelope floor {floor:F4} at N={N}");
        Assert.True(endpoint > floor, $"Endpoint {endpoint:F4} must exceed the bare envelope floor {floor:F4} at N={N}");
    }

    [Fact]
    public void HwhmLeftOverQPeakMean_EqualsPerBondAverage_AtN5_ByChainMirror()
    {
        // The class-mean is computed via average-curves-first-then-peak-find (canonical
        // Python pipeline). For N=5 this is reasonably close to the per-bond ratio average
        // because the only Interior bonds are the symmetric pair (1, 2), but the two
        // approaches diverge at higher N (N=7, 8 Interior includes off-symmetric bonds).
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);

        double endpointMean = ratio.HwhmLeftOverQPeakMean(BondClass.Endpoint);
        double endpointFromWitnessAverage =
            (ratio.Witnesses[0].HwhmLeftOverQPeak + ratio.Witnesses[3].HwhmLeftOverQPeak) / 2.0;
        // For N=5 the chain mirror symmetry between bonds 0 and 3 makes their K curves
        // identical (Endpoint pair), so the curve-mean = per-bond-mean.
        Assert.Equal(endpointFromWitnessAverage, endpointMean, 4);

        double interiorMean = ratio.HwhmLeftOverQPeakMean(BondClass.Interior);
        double interiorFromWitnessAverage =
            (ratio.Witnesses[1].HwhmLeftOverQPeak + ratio.Witnesses[2].HwhmLeftOverQPeak) / 2.0;
        // Same chain-mirror argument for the Interior pair (1, 2) at N=5.
        Assert.Equal(interiorFromWitnessAverage, interiorMean, 4);
    }
}
