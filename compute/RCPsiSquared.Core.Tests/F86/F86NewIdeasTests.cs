using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>Tests for the two new F86 findings discovered 2026-05-03 in the OOP-driven
/// exploration: σ_0 chromaticity scaling and F71 mirror invariance of per-bond Q_peak.
/// </summary>
public class F86NewIdeasTests
{
    [Fact]
    public void SigmaZeroChromaticityScaling_AsymptoteFormula_IsTwoSqrtTwoCMinusOne()
    {
        Assert.Equal(2.0 * Math.Sqrt(2.0), SigmaZeroChromaticityScaling.Asymptote(c: 2), 12);
        Assert.Equal(4.0, SigmaZeroChromaticityScaling.Asymptote(c: 3), 12);
        Assert.Equal(2.0 * Math.Sqrt(6.0), SigmaZeroChromaticityScaling.Asymptote(c: 4), 12);
        // Q_EP asymptote
        Assert.Equal(1.0 / Math.Sqrt(2.0), SigmaZeroChromaticityScaling.QEpAsymptote(c: 2), 12);
        Assert.Equal(0.5, SigmaZeroChromaticityScaling.QEpAsymptote(c: 3), 12);
        Assert.Equal(1.0 / Math.Sqrt(6.0), SigmaZeroChromaticityScaling.QEpAsymptote(c: 4), 12);
    }

    [Fact]
    public void SigmaZeroChromaticityScaling_C2N7_HitsAsymptoteWithinTenMicro()
    {
        // Sweet-spot observation: c=2 N=7 σ_0 = 2√2 to within 10⁻⁵.
        var w = new SigmaZeroScalingWitness(chromaticity: 2, n: 7, gammaZero: 0.05);
        Assert.InRange(w.NormalisedRatio, 2.0 - 1e-5, 2.0 + 1e-5);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void SigmaZeroChromaticityScaling_NormalisedRatio_ConvergesMonotonicallyTo2(int c)
    {
        // Within each chromaticity, σ_0/√(2(c−1)) should grow monotonically with N
        // toward the asymptote 2.0.
        int[] Ns = c switch
        {
            2 => new[] { 5, 6, 7, 8 },
            3 => new[] { 5, 6, 7, 8 },
            4 => new[] { 7, 8 },
            _ => throw new ArgumentException(),
        };

        var ratios = Ns
            .Select(N => new SigmaZeroScalingWitness(c, N, 0.05).NormalisedRatio)
            .ToArray();

        for (int i = 0; i < ratios.Length - 1; i++)
            Assert.True(ratios[i] < ratios[i + 1],
                $"σ_0/√(2(c−1)) should grow monotonically with N at c={c}: " +
                $"N={Ns[i]}→{ratios[i]:F4}, N={Ns[i + 1]}→{ratios[i + 1]:F4}");

        // All values approach 2.0 from below at the tested envelope (0.65 < ratio < 2.01).
        Assert.All(ratios, r => Assert.InRange(r, 1.65, 2.01));
    }

    [Theory]
    [InlineData(2, 5)]
    [InlineData(2, 6)]
    [InlineData(3, 5)]
    [InlineData(3, 6)]
    public void F71MirrorInvariance_PerBondQPeak_BitExactSymmetricUnderBondMirror(int c, int N)
    {
        // For every bond pair (b, N−2−b), Q_peak(b) = Q_peak(N−2−b) bit-exactly. This is
        // the F71 spatial-mirror algebraic identity applied to F86's per-bond observable.
        int n = c - 1;
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        // Wide grid so Endpoint Q_peak (~2.5) lands inside.
        var qGrid = ResonanceScan.LinearQGrid(0.20, 6.00, 60);
        var curve = new ResonanceScan(block).ComputeKCurve(qGrid);

        double maxDev = F71MirrorInvariance.MaxMirrorDeviation(curve);
        Assert.True(maxDev < 1e-10,
            $"F71 mirror invariance violated at c={c} N={N}: max |Q_peak(b) − Q_peak(N−2−b)| = {maxDev:E3}");
    }

    [Fact]
    public void F86KnowledgeBase_HasNewSigma0ScalingAndF71MirrorClaims()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);

        Assert.NotNull(kb.Sigma0Scaling);
        Assert.NotNull(kb.F71Mirror);
        Assert.Equal(Tier.Tier1Candidate, kb.Sigma0Scaling.Tier);
        Assert.Equal(Tier.Tier1Derived, kb.F71Mirror.Tier);
    }
}
