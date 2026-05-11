using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89PathKAtLockMechanismClaimTests
{
    private static F89PathKAtLockMechanismClaim BuildClaim() =>
        new F89PathKAtLockMechanismClaim(new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim()));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(2, 1)]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 2)]
    [InlineData(6, 3)]
    [InlineData(7, 3)]
    [InlineData(8, 4)]
    public void FaCount_IsFloorNBlockHalf(int nBlock, int expected)
    {
        Assert.Equal(expected, F89PathKAtLockMechanismClaim.FaCount(nBlock));
    }

    [Theory]
    [InlineData(4, 2)]
    [InlineData(5, 2)]
    [InlineData(6, 3)]
    [InlineData(7, 3)]
    public void SeAntiBlochOrbit_HasFloorNBlockHalfElements(int nBlock, int expected)
    {
        var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(nBlock);
        Assert.Equal(expected, orbit.Count);
        Assert.All(orbit, n => Assert.True(n % 2 == 0, $"SE-anti orbit element {n} should be even"));
        Assert.All(orbit, n => Assert.True(n >= 2 && n <= nBlock, $"SE-anti orbit element {n} out of range"));
    }

    [Theory]
    [InlineData(4, 2, 1.2360679774997898)]
    [InlineData(4, 4, -3.2360679774997894)]
    [InlineData(5, 2, 2.0)]
    [InlineData(5, 4, -2.0)]
    [InlineData(7, 4, 0.0)]
    public void BlochEigenvalueY_MatchesAnalytical(int nBlock, int n, double expected)
    {
        double y = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, n);
        Assert.Equal(expected, y, precision: 12);
    }

    [Fact]
    public void BlochEigenvalueY_PathSixZeroMode_IsExactlyZero()
    {
        Assert.Equal(0.0, F89PathKAtLockMechanismClaim.BlochEigenvalueY(7, 4), precision: 14);
    }

    [Fact]
    public void FaCount_NBlockBelowTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKAtLockMechanismClaim.FaCount(1));
    }

    [Fact]
    public void SeAntiBlochOrbit_NBlockBelowTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(1));
    }
}
