using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.Resonance;

public class EpAlgebraTests
{
    [Theory]
    [InlineData(0.025, 10.0)]
    [InlineData(0.05, 5.0)]
    [InlineData(0.10, 2.5)]
    [InlineData(1.0, 0.25)]
    public void TPeak_IsInverseFourGamma(double gamma, double expected)
    {
        Assert.Equal(expected, EpAlgebra.TPeak(gamma), 12);
    }

    [Fact]
    public void TPeak_RejectsNonPositiveGamma()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => EpAlgebra.TPeak(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => EpAlgebra.TPeak(-0.05));
    }

    [Theory]
    [InlineData(1.0, 2.0)]
    [InlineData(2.0, 1.0)]
    [InlineData(2.8284271247461903, 0.7071067811865476)]   // 2/(2√2) = 1/√2
    public void QEp_IsTwoOverGEff(double gEff, double expected)
    {
        Assert.Equal(expected, EpAlgebra.QEp(gEff), 12);
    }

    [Fact]
    public void SlowestPair_AtJZero_GivesPureRates()
    {
        // λ_±(k=1) at J=0: -4γ₀ ± 2γ₀ = -2γ₀ (HD=1) and -6γ₀ (HD=3). With γ₀=0.05: -0.1, -0.3.
        var (lp, lm) = EpAlgebra.SlowestPairEigenvalues(gammaZero: 0.05, j: 0.0, gEff: 1.0);
        Assert.Equal(-0.1, lp, 12);
        Assert.Equal(-0.3, lm, 12);
    }

    [Fact]
    public void SlowestPair_AtEp_Coalesces()
    {
        // EP: J·g_eff = 2γ₀ → eigenvalues coalesce at -4γ₀.
        double gamma = 0.05;
        double gEff = 1.0;
        double jEp = 2.0 * gamma / gEff;
        var (lp, lm) = EpAlgebra.SlowestPairEigenvalues(gamma, jEp, gEff);
        Assert.Equal(-4.0 * gamma, lp, 12);
        Assert.Equal(-4.0 * gamma, lm, 12);
    }
}
