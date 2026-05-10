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
    /// <summary>Empirical anchor table (PROOF_F86_QPEAK Statement 2 + section "HWHM_left/Q_peak
    /// across all tested cases (two bond classes)" at γ₀ = 0.05). Tolerance 0.005 is generous
    /// against the 0.025 Q-grid spacing — the empirical pipeline lands within 1e-3 of the
    /// canonical Python pipeline values.</summary>
    [Theory]
    [InlineData(5, BondClass.Interior, 0.7455, 0.005)]
    [InlineData(5, BondClass.Endpoint, 0.7700, 0.005)]
    [InlineData(6, BondClass.Interior, 0.7529, 0.005)]
    [InlineData(6, BondClass.Endpoint, 0.7738, 0.005)]
    [InlineData(7, BondClass.Interior, 0.7507, 0.005)]
    [InlineData(7, BondClass.Endpoint, 0.7738, 0.005)]
    [InlineData(8, BondClass.Interior, 0.7531, 0.005)]
    [InlineData(8, BondClass.Endpoint, 0.7734, 0.005)]
    public void HwhmRatio_C2_MatchesEmpiricalAnchor(int N, BondClass cls, double expected, double tol)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        double r = ratio.HwhmLeftOverQPeakMean(cls);
        Assert.True(Math.Abs(r - expected) <= tol,
            $"HWHM_left/Q_peak at N={N}, {cls}: got {r:F4}, expected {expected:F4} +/- {tol}");
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
        // Empirical fact: at c=2 N=5..8, Endpoint HWHM/Q* > Interior HWHM/Q* by ~0.022.
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
    public void Build_ThrowsGridEdgeEscape_AtN9_WithDefaultGrid_AndExplicitStrictMode()
    {
        // Tom 2026-05-10: "Lass uns auch eine Exception einbauen, wenn der Peak exakt
        // DefaultQGrid Grenze ist, das ist ja Fehlerquelle pur."
        //
        // At N=9 c=2 with the default Q-grid [0.20, 4.00], the first-flanking interior
        // orbit (b=1 ↔ b=N−2−1=6) has its true Q_peak above 4.0 — confirmed empirically
        // 2026-05-10 via extended-grid scan at N=11 where Q_peak=8.79. The default-grid
        // peak finder snaps it to Q=4.0 which is a grid artefact, not physical.
        //
        // With throwOnGridEdgeSnap: true the snap raises GridEdgeEscapeException so the
        // caller cannot silently consume meaningless Q_peak/HWHM data. Default behavior
        // (false) preserves backwards compatibility for tests that explicitly study escape
        // (see PerF71OrbitKTableTests.IsEscaped_FlagsFlankingOrbit_AtN9_WithDefaultGrid).
        var block = new CoherenceBlock(N: 9, n: 1, gammaZero: 0.05);
        var ex = Assert.Throws<GridEdgeEscapeException>(
            () => C2HwhmRatio.Build(block, qGrid: null, throwOnGridEdgeSnap: true));
        // Sanity: the escaped bond is one of the flanking-1 pair (b=1 or b=6 at N=9
        // numBonds=8; F71 mirror pairs b=1 ↔ b=6). The peak finder picks whichever
        // shows up first in the scan; bit-identical mirror pair both at grid edge.
        Assert.True(ex.Bond == 1 || ex.Bond == 6,
            $"escaped bond should be flanking-1 pair (b=1 or b=6) at N=9; got b={ex.Bond}");
        Assert.Equal(4.0, ex.GridUpper, precision: 12);
        Assert.True(ex.QPeak >= ex.GridUpper - ex.DQ,
            $"Q_peak {ex.QPeak} should be within one dQ {ex.DQ} of grid upper {ex.GridUpper}");
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
        // The 2026-05-06 Direction (b) doubled-PTF Ansatz gives 0.671535 as the SVD-block-only
        // FLOOR. Empirical Interior + Endpoint HWHM ratios sit above this floor; the gap
        // (~0.08-0.10) is the probe-block 2-level sub-resonance contribution. Both bond classes
        // must be above the floor for the structural explanation to hold.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var ratio = C2HwhmRatio.Build(block);
        var interior = ratio.HwhmLeftOverQPeakMean(BondClass.Interior);
        var endpoint = ratio.HwhmLeftOverQPeakMean(BondClass.Endpoint);
        Assert.True(interior > C2HwhmRatio.BareDoubledPtfHwhmRatio,
            $"Interior {interior:F4} must exceed bare-doubled-PTF floor {C2HwhmRatio.BareDoubledPtfHwhmRatio:F4} at N={N}");
        Assert.True(endpoint > C2HwhmRatio.BareDoubledPtfHwhmRatio,
            $"Endpoint {endpoint:F4} must exceed bare-doubled-PTF floor {C2HwhmRatio.BareDoubledPtfHwhmRatio:F4} at N={N}");
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
