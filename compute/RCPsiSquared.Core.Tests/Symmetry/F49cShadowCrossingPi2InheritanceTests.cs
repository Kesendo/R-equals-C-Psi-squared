using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F49cShadowCrossingPi2InheritanceTests
{
    private static F49cShadowCrossingPi2Inheritance BuildClaim() =>
        new F49cShadowCrossingPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(2, 4.0)]    // 4^1
    [InlineData(3, 16.0)]   // 4^2
    [InlineData(4, 64.0)]   // 4^3
    [InlineData(5, 256.0)]  // 4^4
    [InlineData(6, 1024.0)] // 4^5
    public void FourPowerNMinus1Factor_MatchesLadder(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().FourPowerNMinus1Factor(N), precision: 12);
    }

    [Theory]
    [InlineData(2, -1)]
    [InlineData(3, -3)]
    [InlineData(4, -5)]
    [InlineData(7, -11)]
    public void LadderIndexFor_Is3Minus2N(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().LadderIndexFor(N));
    }

    [Theory]
    [InlineData(2, 1)]
    [InlineData(3, 2)]
    [InlineData(6, 5)]
    public void OperatorSpaceQubitCountFor_IsNMinus1(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().OperatorSpaceQubitCountFor(N));
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void MirrorPinnedFourPowerNMinus1_AgreesWithLadder(int N)
    {
        var f = BuildClaim();
        Assert.Equal(f.FourPowerNMinus1Factor(N), f.MirrorPinnedFourPowerNMinus1(N), precision: 12);
    }

    [Theory]
    [InlineData(2, 1)]
    [InlineData(3, 2)]
    [InlineData(4, 3)]
    [InlineData(7, 6)]
    public void SpectatorCount_IsNMinus1(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().SpectatorCount(N));
    }

    [Theory]
    [InlineData(2, 1.0 / 2.0)]
    [InlineData(3, 2.0 / 3.0)]
    [InlineData(4, 3.0 / 4.0)]
    [InlineData(6, 5.0 / 6.0)]
    public void VarianceRatio_IsNMinus1OverN(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().VarianceRatio(N), precision: 12);
    }

    [Theory]
    // R²(N) = (N-1) / (N · 4^(N-1))
    [InlineData(2, 1.0 / 8.0)]                     // 1 / (2 · 4) = 1/8
    [InlineData(3, 2.0 / 48.0)]                    // 2 / (3 · 16) = 2/48 = 1/24
    [InlineData(4, 3.0 / 256.0)]                   // 3 / (4 · 64) = 3/256
    [InlineData(5, 4.0 / 1280.0)]                  // 4 / (5 · 256) = 4/1280 = 1/320
    public void RSquared_MatchesClosedForm(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().RSquared(N), precision: 12);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void R_IsSqrtOfRSquared(int N)
    {
        var f = BuildClaim();
        Assert.Equal(Math.Sqrt(f.RSquared(N)), f.R(N), precision: 14);
    }

    [Fact]
    public void F49ShadowBalancedRSquared_AtN2_IsZero()
    {
        // F49 has (N-2) numerator: at N=2 the spectator count is 0,
        // exact Pythagorean decomposition (R = 0).
        Assert.Equal(0.0, BuildClaim().F49ShadowBalancedRSquared(2), precision: 14);
    }

    [Theory]
    // F49: R²(N) = (N-2) / (N · 4^(N-1))
    [InlineData(3, 1.0 / 48.0)]   // 1 / 48
    [InlineData(4, 2.0 / 256.0)]  // 2 / 256 = 1/128
    [InlineData(5, 3.0 / 1280.0)] // 3 / 1280
    public void F49ShadowBalancedRSquared_MatchesF49ClosedForm(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().F49ShadowBalancedRSquared(N), precision: 12);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ShadowCrossingMinusBalancedGap_EqualsOneOverNTimesFourPowerNMinus1(int N)
    {
        var f = BuildClaim();
        // R²_F49c − R²_F49 = ((N-1) - (N-2)) / (N · 4^(N-1)) = 1 / (N · 4^(N-1))
        double expected = 1.0 / (N * f.FourPowerNMinus1Factor(N));
        Assert.Equal(expected, f.ShadowCrossingMinusBalancedGap(N), precision: 14);
        Assert.Equal(f.RSquared(N) - f.F49ShadowBalancedRSquared(N), f.ShadowCrossingMinusBalancedGap(N), precision: 14);
    }

    [Fact]
    public void FourPowerNMinus1Factor_ThrowsForNLessThanTwo()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().FourPowerNMinus1Factor(1));
    }

    [Fact]
    public void RSquared_ThrowsForNLessThanTwo()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().RSquared(1));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F49cShadowCrossingPi2Inheritance(
                ladder: null!,
                mirror: new Pi2OperatorSpaceMirrorClaim()));
    }

    [Fact]
    public void Constructor_NullMirror_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F49cShadowCrossingPi2Inheritance(
                ladder: new Pi2DyadicLadderClaim(),
                mirror: null!));
    }
}
