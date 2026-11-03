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

    [Fact]
    public void FirstLocalMaximum_TakesTheLowPeak_FindTakesTheHigherFarOne()
    {
        // Two parabolic bumps on an integer grid: a lower one with vertex at x = 2 (height 3)
        // and a higher one at x = 8 (height 5). Parabolas are fitted exactly by the 3-point fit,
        // so both vertices come back exactly (small integers, no rounding anywhere in the fit).
        double[] qs = Enumerable.Range(0, 11).Select(i => (double)i).ToArray();
        double[] ks = qs.Select(x => Math.Max(3.0 - (x - 2) * (x - 2), 5.0 - (x - 8) * (x - 8))).ToArray();
        var canonical = ParabolicPeakFinder.FindFirstLocalMaximum(qs, ks);
        var global = ParabolicPeakFinder.Find(qs, ks);
        Assert.Equal(2.0, canonical.QPeak);
        Assert.Equal(3.0, canonical.KMax);
        Assert.Equal(8.0, global.QPeak);
        Assert.Equal(5.0, global.KMax);
    }

    [Fact]
    public void FirstLocalMaximum_WithoutAnInteriorMaximum_IsTheGlobalOne()
    {
        double[] qs = { 0, 1, 2, 3, 4 };
        double[] ks = { 0, 1, 2, 3, 4 };
        Assert.Equal(ParabolicPeakFinder.Find(qs, ks), ParabolicPeakFinder.FindFirstLocalMaximum(qs, ks));
    }
}
