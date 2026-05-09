using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F4StationaryModeCountPi2InheritanceTests
{
    private static F4StationaryModeCountPi2Inheritance BuildClaim() =>
        new F4StationaryModeCountPi2Inheritance(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void IrrepDimensionCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().IrrepDimensionCoefficient, precision: 14);
    }

    [Theory]
    // SpinMultiplicity verified against Clebsch-Gordan tables:
    // N=2: J=0 (twoJ=0) m=1, J=1 (twoJ=2) m=1
    // N=3: J=1/2 (twoJ=1) m=2, J=3/2 (twoJ=3) m=1
    // N=4: J=0 (twoJ=0) m=2, J=1 (twoJ=2) m=3, J=2 (twoJ=4) m=1
    // N=5: J=1/2 (twoJ=1) m=5, J=3/2 (twoJ=3) m=4, J=5/2 (twoJ=5) m=1
    [InlineData(2, 0, 1)]
    [InlineData(2, 2, 1)]
    [InlineData(3, 1, 2)]
    [InlineData(3, 3, 1)]
    [InlineData(4, 0, 2)]
    [InlineData(4, 2, 3)]
    [InlineData(4, 4, 1)]
    [InlineData(5, 1, 5)]
    [InlineData(5, 3, 4)]
    [InlineData(5, 5, 1)]
    public void SpinMultiplicity_MatchesClebschGordan(int N, int twoJ, int expected)
    {
        Assert.Equal(expected, BuildClaim().SpinMultiplicity(N, twoJ));
    }

    [Theory]
    [InlineData(0, 1)]
    [InlineData(1, 2)]
    [InlineData(2, 3)]
    [InlineData(5, 6)]
    public void IrrepDimension_EqualsTwoJPlusOne(int twoJ, int expected)
    {
        Assert.Equal(expected, BuildClaim().IrrepDimension(twoJ));
    }

    [Theory]
    // Verified: manual Clebsch-Gordan computation
    // N=6: J=0 (m=5, c=5), J=1 (m=9, c=81), J=2 (m=5, c=125), J=3 (m=1, c=49) → 260
    [InlineData(2, 10)]
    [InlineData(3, 24)]
    [InlineData(4, 54)]
    [InlineData(5, 120)]
    [InlineData(6, 260)]
    public void StationaryModeCount_MatchesVerifiedTable(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().StationaryModeCount(N));
    }

    [Theory]
    // Schur-Weyl identity: Σ_J m(J)·(2J+1) = 2^N
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(10)]
    public void SchurWeylDimensionIdentityHolds(int N)
    {
        Assert.True(BuildClaim().SchurWeylDimensionIdentityHolds(N));
    }

    [Fact]
    public void SpinDecomposition_N4_HasThreeJValues()
    {
        // N=4: twoJ ∈ {0, 2, 4} → three J-multiplets.
        var decomp = BuildClaim().SpinDecomposition(4);
        Assert.Equal(3, decomp.Count);
        // Sorted descending in twoJ
        Assert.Equal(4, decomp[0].TwoJ);   // J=2: m=1, (2J+1)²=25, c=25
        Assert.Equal(2, decomp[1].TwoJ);   // J=1: m=3, (2J+1)²=9, c=27
        Assert.Equal(0, decomp[2].TwoJ);   // J=0: m=2, (2J+1)²=1, c=2
        Assert.Equal(25 + 27 + 2, decomp.Sum(t => t.Contribution));
    }

    [Fact]
    public void SpinMultiplicity_WrongParity_Throws()
    {
        // N=4 even → twoJ must be even; twoJ=1 (odd) is wrong parity.
        Assert.Throws<ArgumentException>(() => BuildClaim().SpinMultiplicity(N: 4, twoJ: 1));
        // N=3 odd → twoJ must be odd; twoJ=2 is wrong parity.
        Assert.Throws<ArgumentException>(() => BuildClaim().SpinMultiplicity(N: 3, twoJ: 2));
    }

    [Fact]
    public void SpinMultiplicity_TwoJOutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SpinMultiplicity(N: 3, twoJ: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SpinMultiplicity(N: 3, twoJ: 4));
    }

    [Fact]
    public void StationaryModeCount_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().StationaryModeCount(N: 0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F4StationaryModeCountPi2Inheritance(null!));
    }
}
