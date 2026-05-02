using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Tests.Numerics;

public class ParabolicPeakFinderTests
{
    [Fact]
    public void GaussianPeak_ParabolicRefinement_HitsTrueApex()
    {
        // Synthetic peak at Q* = 1.7 with width σ. Sampled on dQ = 0.025 grid.
        double qStar = 1.7;
        double sigma = 0.5;
        var qs = new List<double>();
        var ks = new List<double>();
        for (double q = 0.2; q <= 4.0 + 1e-9; q += 0.025)
        {
            qs.Add(q);
            ks.Add(Math.Exp(-Math.Pow((q - qStar) / sigma, 2)));
        }
        var info = ParabolicPeakFinder.Find(qs, ks);
        Assert.InRange(info.QPeak, qStar - 0.01, qStar + 0.01);
        Assert.InRange(info.KMax, 0.99, 1.001);
    }

    [Fact]
    public void Hwhm_OfStandardGaussian_MatchesAnalyticalFormula()
    {
        // For y = exp(-(x/σ)²), HWHM = σ·√(ln 2). At σ=0.5: HWHM = 0.5·√(ln 2) ≈ 0.4163.
        var qs = new List<double>();
        var ks = new List<double>();
        double sigma = 0.5;
        for (double q = 0.2; q <= 4.0 + 1e-9; q += 0.005)
        {
            qs.Add(q);
            ks.Add(Math.Exp(-Math.Pow((q - 1.7) / sigma, 2)));
        }
        var info = ParabolicPeakFinder.Find(qs, ks);
        double expected = sigma * Math.Sqrt(Math.Log(2.0));
        Assert.NotNull(info.HwhmLeft);
        Assert.NotNull(info.HwhmRight);
        Assert.InRange(info.HwhmLeft!.Value, expected - 0.005, expected + 0.005);
        Assert.InRange(info.HwhmRight!.Value, expected - 0.005, expected + 0.005);
    }
}
