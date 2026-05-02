using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.Resonance;

public class ResonanceScanTests
{
    [Fact]
    public void DefaultQGrid_Has153Points_From020To400_Step0025()
    {
        var grid = ResonanceScan.DefaultQGrid();
        Assert.Equal(153, grid.Length);
        Assert.Equal(0.20, grid[0], 12);
        Assert.Equal(4.00, grid[^1], 12);
        Assert.Equal(0.025, grid[1] - grid[0], 12);
    }

    [Fact]
    public void DefaultTGrid_Spans_06_To_16_TPeak()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var scan = new ResonanceScan(block);
        var t = scan.DefaultTGrid();
        // t_peak = 1/(4·0.05) = 5
        Assert.Equal(21, t.Length);
        Assert.Equal(0.6 * 5.0, t[0], 12);
        Assert.Equal(1.6 * 5.0, t[^1], 12);
    }

    [Fact]
    public void KCurve_C2_N5_Interior_MatchesPythonStepF()
    {
        // step_f c=2 N=5 Interior: Q_peak = 1.4821, HWHM-/Q* = 0.7455, |K|max = 0.14681
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var scan = new ResonanceScan(block);
        var curve = scan.ComputeKCurve();

        var peak = curve.Peak(BondClass.Interior);
        Assert.InRange(peak.QPeak, 1.481, 1.484);
        Assert.NotNull(peak.HwhmLeftOverQPeak);
        Assert.InRange(peak.HwhmLeftOverQPeak!.Value, 0.744, 0.747);
        Assert.InRange(peak.KMax, 0.146, 0.148);
    }

    [Fact]
    public void KCurve_C2_N5_Endpoint_MatchesPythonStepF()
    {
        // step_f c=2 N=5 Endpoint: Q_peak = 2.5008, HWHM-/Q* = 0.7700, |K|max = 0.27673
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var scan = new ResonanceScan(block);
        var curve = scan.ComputeKCurve();

        var peak = curve.Peak(BondClass.Endpoint);
        Assert.InRange(peak.QPeak, 2.499, 2.503);
        Assert.NotNull(peak.HwhmLeftOverQPeak);
        Assert.InRange(peak.HwhmLeftOverQPeak!.Value, 0.769, 0.771);
        Assert.InRange(peak.KMax, 0.275, 0.278);
    }
}
