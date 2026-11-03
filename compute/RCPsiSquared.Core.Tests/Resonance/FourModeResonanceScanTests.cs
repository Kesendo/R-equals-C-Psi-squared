using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.Resonance;

/// <summary>F86 documentary tests for the 4-mode minimal effective model.
///
/// <para>The structural tests verify that the 4-mode scan produces a well-formed K-curve.
/// The numerical fingerprint tests record the empirical answer to F86's open question
/// (PROOF_F86_QPEAK "What's missing for full Tier 1" item 1):</para>
/// <list type="bullet">
///   <item>4-mode <b>Interior</b> at odd N: Q_peak 1.90 (N = 5) and 3.08 (N = 7) against the full
///         block's 1.48 and 1.60, HWHM_left/Q_peak 0.740 and 0.762. At even N (6, 8) the interior
///         curve is still rising at Q = 4, the grid edge.</item>
///   <item>4-mode <b>Endpoint</b>: still rising at Q = 4 for N = 5, 7, 8; at N = 6 it peaks inside
///         the grid at Q = 3.78, against the full block's 2.58.</item>
///   <item>4-mode <b>K_max</b> is 3–10× smaller than full block-L — the 4-mode subspace
///         carries only a fraction of the observable signal.</item>
/// </list>
/// <para>The values are pinned against an independent route: the same four-mode construction in
/// Python (the R-even top pair at even N, see <see cref="RCPsiSquared.Core.Decomposition.FourModeBasis"/>),
/// maximised over t by a dense scan on [0.3, 40] plus bounded Brent:
/// simulations/f86_envelope_second_route.py --fourmode N.</para>
/// <para>Conclusion: the 4-mode minimal effective is a structural witness (probe lives in
/// channel-uniform 2D, EP partners in SVD-top 2D, with cross-coupling carrying bond-position
/// information) but NOT sufficient to reproduce the universal shape numerically. More modes
/// are needed for Tier 1 promotion of the K_class(Q) = f_class(Q/Q_EP) statement.</para>
/// </summary>
public class FourModeResonanceScanTests
{
    [Fact]
    public void KCurve_DefaultGrids_ProducesShapeCompatibleWithFullScan()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var curve = new FourModeResonanceScan(block).ComputeKCurve();

        Assert.Equal(block.NumBonds, curve.NumBonds);
        Assert.Equal(ResonanceScan.DefaultQGrid().Length, curve.QGrid.Count);
        Assert.Equal(0.20, curve.QGrid[0], 12);
        Assert.Equal(4.00, curve.QGrid[^1], 12);
    }

    [Theory]
    [InlineData(5, 1.9005775355788428, 0.739558136499575)]
    [InlineData(7, 3.078068625836869, 0.7615384762563832)]
    public void Interior_FourMode_OddN_QPeakAndHwhmRatio(int N, double qPeak, double hwhmRatio)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var interior = new FourModeResonanceScan(block).ComputeKCurve().Peak(BondClass.Interior);
        Assert.True(Math.Abs(interior.QPeak - qPeak) <= 1e-10, $"N={N} Q_peak {interior.QPeak:R}");
        Assert.NotNull(interior.HwhmLeftOverQPeak);
        Assert.True(Math.Abs(interior.HwhmLeftOverQPeak!.Value - hwhmRatio) <= 1e-10,
            $"N={N} HWHM_left/Q_peak {interior.HwhmLeftOverQPeak.Value:R}");
    }

    [Theory]
    [InlineData(6, BondClass.Interior)]
    [InlineData(8, BondClass.Interior)]
    [InlineData(5, BondClass.Endpoint)]
    [InlineData(7, BondClass.Endpoint)]
    [InlineData(8, BondClass.Endpoint)]
    public void FourMode_QPeakBeyondDefaultGrid(int N, BondClass cls)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var peak = new FourModeResonanceScan(block).ComputeKCurve().Peak(cls);
        Assert.True(peak.QPeak >= 3.99, $"c=2 N={N} {cls} 4-mode peaks at Q={peak.QPeak:R} inside the grid");
    }

    [Fact]
    public void Endpoint_FourMode_N6_PeaksInsideTheGrid()
    {
        var block = new CoherenceBlock(6, n: 1, gammaZero: 0.05);
        var endpoint = new FourModeResonanceScan(block).ComputeKCurve().Peak(BondClass.Endpoint);
        Assert.True(Math.Abs(endpoint.QPeak - 3.7812635958294947) <= 1e-10, $"Q_peak {endpoint.QPeak:R}");
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(8, 1)]
    public void Interior_FourMode_KMax_FractionOfFullBlockL(int N, int n)
    {
        // Documented finding: 4-mode K_max is a fraction of the full block-L K_max —
        // the 4-mode subspace carries only part of the K_CC_pr signal. We assert the
        // fraction is small (< 0.5) to lock in the "more modes needed" conclusion.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var fullCurve = new ResonanceScan(block).ComputeKCurve();
        var fourCurve = new FourModeResonanceScan(block).ComputeKCurve();

        double kFull = fullCurve.Peak(BondClass.Interior).KMax;
        double kFour = fourCurve.Peak(BondClass.Interior).KMax;
        double ratio = kFour / kFull;
        Assert.True(ratio > 0 && ratio < 0.5,
            $"c=2 N={N} Interior K_max ratio (4-mode / full) = {ratio:F3}; expected (0, 0.5).");
    }
}
