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
    public void DefaultTGrid_Spans_02_To_24_TPeak()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var scan = new ResonanceScan(block);
        var t = scan.DefaultTGrid();
        // t_peak = 1/(4·0.05) = 5
        Assert.Equal(45, t.Length);
        Assert.Equal(0.2 * 5.0, t[0], 12);
        Assert.Equal(2.4 * 5.0, t[^1], 12);
    }

    /// <summary>The c = 2, N = 5 bond-class peaks of the envelope max_t |K_b(Q, t)|, pinned
    /// against an independent route, <c>simulations/f86_envelope_second_route.py</c>: K_b from the
    /// F89 eigendecomposition of <c>simulations/f89_to_f86_kbond_via_eigendecomp.py</c> (its J being
    /// half of F86's), the maximum over t by a dense scan on [0.3, 40] plus bounded Brent rather than this scan's
    /// seed grid plus golden section, then the same parabolic peak reading on the default Q grid.
    /// The two mirror bonds of each class agree on that route to about 3e-12, which sets the gate.</summary>
    [Theory]
    [InlineData(BondClass.Interior, 1.4831232668733452, 0.7523433331321849, 0.14681468195933925)]
    [InlineData(BondClass.Endpoint, 2.517071370613079, 0.7765595000526813, 0.2767357997650519)]
    public void KCurve_C2_N5_ClassPeaks_MatchTheEnvelopeSecondRoute(
        BondClass cls, double qPeak, double hwhmRatio, double kMax)
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var curve = new ResonanceScan(block).ComputeKCurve();

        var peak = curve.Peak(cls);
        Assert.True(Math.Abs(peak.QPeak - qPeak) <= 1e-11, $"Q_peak {peak.QPeak:R}");
        Assert.NotNull(peak.HwhmLeftOverQPeak);
        Assert.True(Math.Abs(peak.HwhmLeftOverQPeak!.Value - hwhmRatio) <= 1e-11,
            $"HWHM_left/Q_peak {peak.HwhmLeftOverQPeak.Value:R}");
        Assert.True(Math.Abs(peak.KMax - kMax) <= 1e-11, $"KMax {peak.KMax:R}");
    }

    /// <summary>Q = J/γ₀ and the t window scales with t_peak = 1/(4γ₀), so at c=2 N=5 the class
    /// Q_peak and HWHM_left/Q_peak are γ₀-invariant: the generator at γ₀ is γ₀ times the generator
    /// at γ₀ = 1 and the same Q, and the seed grid in units of t_peak is identical. What remains is
    /// the eigensolver's round-off on differently scaled matrices, which the parabolic vertex turns
    /// into ~10⁻¹¹ (the scan's reproducibility); the gate 1e-10 holds that.</summary>
    [Theory]
    [InlineData(0.025)]
    [InlineData(0.10)]
    public void ClassPeaks_AreGammaZeroInvariant(double gammaZero)
    {
        var reference = new ResonanceScan(new CoherenceBlock(5, n: 1, gammaZero: 0.05)).ComputeKCurve();
        var scaled = new ResonanceScan(new CoherenceBlock(5, n: 1, gammaZero: gammaZero)).ComputeKCurve();
        foreach (var cls in new[] { BondClass.Interior, BondClass.Endpoint })
        {
            var a = reference.Peak(cls);
            var b = scaled.Peak(cls);
            Assert.True(Math.Abs(a.QPeak - b.QPeak) <= 1e-10, $"{cls} Q_peak {a.QPeak:R} vs {b.QPeak:R}");
            Assert.True(Math.Abs(a.HwhmLeftOverQPeak!.Value - b.HwhmLeftOverQPeak!.Value) <= 1e-10,
                $"{cls} ratio {a.HwhmLeftOverQPeak:R} vs {b.HwhmLeftOverQPeak:R}");
        }
    }
}
