using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F26CPsiPauliChannelsPi2InheritanceTests
{
    private static F26CPsiPauliChannelsPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, new QuarterAsBilinearMaxvalClaim());
        return new F26CPsiPauliChannelsPi2Inheritance(ladder, f25);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void DecayRateCoefficient_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().DecayRateCoefficient, precision: 14);
    }

    [Fact]
    public void NormalizationDenominator_IsExactlyTwelve()
    {
        Assert.Equal(12.0, BuildClaim().NormalizationDenominator, precision: 14);
    }

    [Theory]
    // F26 rate coefficients: α = 4(γ_y+γ_z), β = 4(γ_x+γ_z), δ = 4(γ_x+γ_y)
    [InlineData(0.05, 0, 0, 0, 0.20, 0.20)]   // pure X: α=0, β=δ=4·0.05
    [InlineData(0, 0.05, 0, 0.20, 0, 0.20)]   // pure Y
    [InlineData(0, 0, 0.05, 0.20, 0.20, 0)]   // pure Z (= F25)
    [InlineData(0.05, 0.05, 0.05, 0.40, 0.40, 0.40)]  // all equal
    public void RateCoefficients_MatchClosedForm(double gx, double gy, double gz, double expectedAlpha, double expectedBeta, double expectedDelta)
    {
        var (a, b, d) = BuildClaim().RateCoefficients(gx, gy, gz);
        Assert.Equal(expectedAlpha, a, precision: 12);
        Assert.Equal(expectedBeta, b, precision: 12);
        Assert.Equal(expectedDelta, d, precision: 12);
    }

    [Fact]
    public void CPsiAtTime_AtTZero_IsOneThird()
    {
        // CΨ(0) = 1·(1+1+1+1)/12 = 4/12 = 1/3
        Assert.Equal(1.0 / 3.0, BuildClaim().CPsiAtTime(0.05, 0.05, 0.05, 0.0), precision: 12);
    }

    [Theory]
    [InlineData(0.05, 0.5)]
    [InlineData(0.10, 1.0)]
    [InlineData(0.01, 5.0)]
    public void RecoversF25AtSingleChannelLimit_AcrossGammaT(double gammaZ, double t)
    {
        Assert.True(BuildClaim().RecoversF25AtSingleChannelLimit(gammaZ, t));
    }

    [Fact]
    public void BellPlusInitialIsOneThird_HoldsExactly()
    {
        Assert.True(BuildClaim().BellPlusInitialIsOneThird());
    }

    [Fact]
    public void CPsiAtTime_DecreasesMonotonically()
    {
        var f = BuildClaim();
        double cpsi0 = f.CPsiAtTime(0.05, 0.05, 0.05, 0.0);
        double cpsi1 = f.CPsiAtTime(0.05, 0.05, 0.05, 0.5);
        double cpsi2 = f.CPsiAtTime(0.05, 0.05, 0.05, 5.0);

        Assert.True(cpsi1 < cpsi0);
        Assert.True(cpsi2 < cpsi1);
    }

    [Fact]
    public void DepolarizingCPsi_AtTZero_IsOneThird()
    {
        Assert.Equal(1.0 / 3.0, BuildClaim().DepolarizingCPsi(0.05, 0.0), precision: 12);
    }

    [Fact]
    public void DepolarizingCPsi_LargeT_GoesToZero()
    {
        Assert.Equal(0.0, BuildClaim().DepolarizingCPsi(0.05, 1000.0), precision: 6);
    }

    [Fact]
    public void RateCoefficients_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().RateCoefficients(-0.05, 0.05, 0.05));
    }

    [Fact]
    public void CPsiAtTime_NegativeT_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().CPsiAtTime(0.05, 0.05, 0.05, -0.1));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, new QuarterAsBilinearMaxvalClaim());
        Assert.Throws<ArgumentNullException>(() =>
            new F26CPsiPauliChannelsPi2Inheritance(null!, f25));
    }

    [Fact]
    public void Constructor_NullF25_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F26CPsiPauliChannelsPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
