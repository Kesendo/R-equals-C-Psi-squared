using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F74ChromaticityPi2InheritanceTests
{
    private static F74ChromaticityPi2Inheritance BuildClaim() =>
        new F74ChromaticityPi2Inheritance(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void RateCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().RateCoefficient, precision: 14);
    }

    [Theory]
    // c(n, N) = min(n, N-1-n) + 1
    // N=5: c(0)=1, c(1)=2, c(2)=3, c(3)=2, c(4)=1
    // N=6: c(0)=1, c(1)=2, c(2)=3, c(3)=3, c(4)=2, c(5)=1
    // N=3: c(0)=1, c(1)=2, c(2)=1
    [InlineData(0, 5, 1)]
    [InlineData(1, 5, 2)]
    [InlineData(2, 5, 3)]
    [InlineData(3, 5, 2)]
    [InlineData(4, 5, 1)]
    [InlineData(0, 6, 1)]
    [InlineData(2, 6, 3)]
    [InlineData(3, 6, 3)]
    [InlineData(0, 3, 1)]
    [InlineData(1, 3, 2)]
    [InlineData(2, 3, 1)]
    public void Chromaticity_MatchesClosedForm(int n, int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().Chromaticity(n, N));
    }

    [Theory]
    [InlineData(2, 5, new[] { 1, 3, 5 })]
    [InlineData(0, 5, new[] { 1 })]
    [InlineData(1, 6, new[] { 1, 3 })]
    public void HammingDistanceValues_MatchesOddLadder(int n, int N, int[] expected)
    {
        var hds = BuildClaim().HammingDistanceValues(n, N);
        Assert.Equal(expected, hds);
    }

    [Theory]
    [InlineData(0.05, 1, 0.1)]    // 2·0.05·1 = 0.1
    [InlineData(0.05, 3, 0.3)]    // 2·0.05·3 = 0.3
    [InlineData(1.0, 5, 10.0)]    // 2·1·5 = 10
    public void PureRate_Equals2GammaTimesHd(double gammaZero, int hd, double expected)
    {
        Assert.Equal(expected, BuildClaim().PureRate(gammaZero, hd), precision: 14);
    }

    [Theory]
    [InlineData(0, 5, true)]
    [InlineData(4, 5, true)]   // N-1
    [InlineData(2, 5, false)]
    public void IsMonochromatic_TrueAtEnds(int n, int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsMonochromatic(n, N));
    }

    [Theory]
    [InlineData(2, 1)]    // N=2 even: max = N/2 = 1
    [InlineData(3, 2)]    // N=3 odd: max = (N+1)/2 = 2
    [InlineData(5, 3)]
    [InlineData(6, 3)]
    [InlineData(7, 4)]
    public void MaxChromaticity_EqualsCeilN_2(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().MaxChromaticity(N));
    }

    [Theory]
    [InlineData(2, 5, true)]    // unique center for odd N
    [InlineData(0, 5, false)]
    [InlineData(2, 6, true)]    // pair centers for even N
    [InlineData(3, 6, true)]
    public void IsCenterBlock_TrueAtMaxChromaticity(int n, int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsCenterBlock(n, N));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(10)]
    public void MirrorPalindromicityHolds(int N)
    {
        Assert.True(BuildClaim().MirrorPalindromicityHolds(N));
    }

    [Fact]
    public void Chromaticity_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Chromaticity(0, N: 1));
    }

    [Fact]
    public void Chromaticity_nOutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Chromaticity(n: -1, N: 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Chromaticity(n: 5, N: 3));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F74ChromaticityPi2Inheritance(null!));
    }
}
