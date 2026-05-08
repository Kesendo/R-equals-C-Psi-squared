using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Spectrum;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Spectrum;

public class W1DispersionTests
{
    private readonly ITestOutputHelper _out;

    public W1DispersionTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var d = new W1Dispersion(N: 5, J: 1.0, gammaZero: 0.05);
        Assert.Equal(Tier.Tier1Derived, d.Tier);
    }

    [Fact]
    public void Anchor_References_AnalyticalSpectrum_AndD10()
    {
        var d = new W1Dispersion(N: 5, J: 1.0, gammaZero: 0.05);
        Assert.Contains("ANALYTICAL_SPECTRUM.md", d.Anchor);
        Assert.Contains("D10_W1_DISPERSION.md", d.Anchor);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(5)]
    [InlineData(8)]
    public void NumModesIsNMinusOne(int N)
    {
        var d = new W1Dispersion(N, J: 1.0, gammaZero: 0.05);
        Assert.Equal(N - 1, d.Frequencies.Count);
        Assert.Equal(N - 1, d.QFactors.Count);
    }

    [Fact]
    public void NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new W1Dispersion(N: 1, J: 1.0, gammaZero: 0.05));
    }

    [Fact]
    public void NonPositiveGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new W1Dispersion(N: 5, J: 1.0, gammaZero: 0));
    }

    /// <summary>Bit-exact reproduction of the 15 pinned frequencies tabulated in
    /// experiments/ANALYTICAL_SPECTRUM.md (J = 1). The pinned table was verified to machine
    /// precision against numerical Liouvillian eigenvalues; we reproduce it here as a
    /// drift-check anchor for the closed-form ω_k = 4J·(1 − cos(πk/N)).</summary>
    [Theory]
    [InlineData(2, 1, 4.000000)]
    [InlineData(3, 1, 2.000000)]
    [InlineData(3, 2, 6.000000)]
    [InlineData(4, 1, 1.171573)]
    [InlineData(4, 2, 4.000000)]
    [InlineData(4, 3, 6.828427)]
    [InlineData(5, 1, 0.763932)]
    [InlineData(5, 2, 2.763932)]
    [InlineData(5, 3, 5.236068)]
    [InlineData(5, 4, 7.236068)]
    [InlineData(6, 1, 0.535898)]
    [InlineData(6, 2, 2.000000)]
    [InlineData(6, 3, 4.000000)]
    [InlineData(6, 4, 6.000000)]
    [InlineData(6, 5, 7.464102)]
    public void Frequency_MatchesPinnedAnalyticalSpectrumTable(int N, int k, double expected)
    {
        var d = new W1Dispersion(N, J: 1.0, gammaZero: 0.05);
        double actual = d.Frequencies[k - 1];
        Assert.Equal(expected, actual, precision: 6);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(5)]
    [InlineData(8)]
    public void UniformDecayRate_EqualsTwoGamma(int N)
    {
        double γ = 0.05;
        var d = new W1Dispersion(N, J: 1.0, gammaZero: γ);
        Assert.Equal(2.0 * γ, d.UniformDecayRate, precision: 12);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(10)]
    public void MeanQ_EqualsTwoJOverGamma(int N)
    {
        double J = 1.0;
        double γ = 0.05;
        var d = new W1Dispersion(N, J, γ);
        // Σ_{k=1}^{N-1} cos(πk/N) = 0 (telescoping identity), so Σ ω_k = (N-1)·4J → mean Q = 2J/γ.
        double meanQEmpirical = d.QFactors.Average();
        Assert.Equal(2.0 * J / γ, d.MeanQ, precision: 12);
        Assert.Equal(d.MeanQ, meanQEmpirical, precision: 10);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(10)]
    public void QSpread_EqualsCotSquaredOfPiOverTwoN(int N)
    {
        var d = new W1Dispersion(N, J: 1.0, gammaZero: 0.05);
        double cot = Math.Cos(Math.PI / (2.0 * N)) / Math.Sin(Math.PI / (2.0 * N));
        Assert.Equal(cot * cot, d.QSpread, precision: 10);
        Assert.Equal(d.MaxQ / d.MinQ, d.QSpread, precision: 10);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(8)]
    public void Bandwidth_EqualsEightJCosPiOverN(int N)
    {
        double J = 1.0;
        var d = new W1Dispersion(N, J, gammaZero: 0.05);
        double expected = 8.0 * J * Math.Cos(Math.PI / N);
        Assert.Equal(expected, d.Bandwidth, precision: 12);
        Assert.Equal(d.Frequencies[N - 2] - d.Frequencies[0], d.Bandwidth, precision: 10);
    }

    [Fact]
    public void Reconnaissance_EmitsSpectrumAcrossN2To6()
    {
        _out.WriteLine("  N | k | ω_k        | Q_k        | MaxQ      | MinQ      | MeanQ    | QSpread");
        _out.WriteLine("  --|---|------------|------------|-----------|-----------|----------|---------");
        foreach (int N in new[] { 2, 3, 4, 5, 6 })
        {
            var d = new W1Dispersion(N, J: 1.0, gammaZero: 0.05);
            for (int k = 1; k <= N - 1; k++)
            {
                _out.WriteLine($"  {N} | {k} | {d.Frequencies[k - 1],10:F6} | " +
                               $"{d.QFactors[k - 1],10:F6} | {d.MaxQ,9:F4} | {d.MinQ,9:F4} | " +
                               $"{d.MeanQ,8:F4} | {d.QSpread,7:F4}");
            }
        }
    }
}
