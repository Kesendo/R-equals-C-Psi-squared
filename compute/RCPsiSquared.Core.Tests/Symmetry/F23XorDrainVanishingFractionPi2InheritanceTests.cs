using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F23XorDrainVanishingFractionPi2InheritanceTests
{
    private static F23XorDrainVanishingFractionPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        return new F23XorDrainVanishingFractionPi2Inheritance(ladder, mirror);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void BaseFactor_IsExactlyFour()
    {
        // 4 = a_{-1} on dyadic ladder.
        Assert.Equal(4.0, BuildClaim().BaseFactor, precision: 14);
    }

    [Theory]
    [InlineData(1, 4.0)]
    [InlineData(2, 16.0)]
    [InlineData(3, 64.0)]
    [InlineData(4, 256.0)]
    [InlineData(5, 1024.0)]
    [InlineData(6, 4096.0)]
    [InlineData(8, 65536.0)]
    [InlineData(10, 1048576.0)]
    public void OperatorSpaceDim_Equals4ToTheN(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().OperatorSpaceDim(N), precision: 8);
    }

    [Theory]
    [InlineData(1, 2)]
    [InlineData(3, 4)]
    [InlineData(5, 6)]
    [InlineData(8, 9)]
    [InlineData(20, 21)]
    public void XorDrainCount_EqualsNPlusOne(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().XorDrainCount(N));
    }

    [Theory]
    // F23 closed form: fraction(XOR) = (N+1)/4^N
    // N=3: 4/64 = 0.0625 = 6.25%
    // N=5: 6/1024 ≈ 0.005859 = 0.586%
    // N=8: 9/65536 ≈ 0.0001373 = 0.0137%
    [InlineData(3, 0.0625)]
    [InlineData(5, 6.0 / 1024.0)]
    [InlineData(8, 9.0 / 65536.0)]
    public void XorDrainFraction_MatchesClosedForm(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().XorDrainFraction(N), precision: 12);
    }

    [Fact]
    public void XorDrainFraction_DecreasesExponentiallyInN()
    {
        // Each successive N should reduce the fraction.
        var f = BuildClaim();
        double prev = f.XorDrainFraction(1);
        for (int N = 2; N <= 20; N++)
        {
            double curr = f.XorDrainFraction(N);
            Assert.True(curr < prev, $"fraction at N={N} ({curr}) should be < fraction at N-1 ({prev})");
            prev = curr;
        }
    }

    [Fact]
    public void XorDrainFraction_AtN20_IsApproximately10ToMinus11()
    {
        // ANALYTICAL_FORMULAS line 191 says "N=20: ~10^-11".
        // Exact: 21/4^20 = 21/1099511627776 ≈ 1.91e-11.
        double frac = BuildClaim().XorDrainFraction(20);
        Assert.True(frac > 1e-12 && frac < 1e-10,
            $"N=20 fraction should be in [1e-12, 1e-10] range; got {frac}");
    }

    [Fact]
    public void IsMacroscopicallyNegligible_FalseAtSmallN()
    {
        var f = BuildClaim();
        Assert.False(f.IsMacroscopicallyNegligible(N: 3));    // 6.25% >> 1e-6
        Assert.False(f.IsMacroscopicallyNegligible(N: 5));    // 0.59% >> 1e-6
        Assert.False(f.IsMacroscopicallyNegligible(N: 8));    // 0.014% >> 1e-6
    }

    [Fact]
    public void IsMacroscopicallyNegligible_TrueAtLargeN()
    {
        // Default threshold = 1e-6. At N=12: 13/4^12 ≈ 7.75e-7 (below).
        var f = BuildClaim();
        Assert.True(f.IsMacroscopicallyNegligible(N: 12));    // ~7.75e-7
        Assert.True(f.IsMacroscopicallyNegligible(N: 15));    // ~10^-8
        Assert.True(f.IsMacroscopicallyNegligible(N: 20));    // ~10^-11
    }

    [Fact]
    public void IsMacroscopicallyNegligible_NonPositiveThreshold_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().IsMacroscopicallyNegligible(N: 5, threshold: 0.0));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(3)]
    [InlineData(6)]
    public void MatchesMirrorTable_HoldsForN1To6(int N)
    {
        // F23's OperatorSpaceDim should match Pi2OperatorSpaceMirror's pinned table.
        Assert.True(BuildClaim().MatchesMirrorTable(N));
    }

    [Fact]
    public void MatchesMirrorTable_OutOfRange_Throws()
    {
        // Mirror table covers N=1..6 only.
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MatchesMirrorTable(N: 7));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MatchesMirrorTable(N: 0));
    }

    [Fact]
    public void OperatorSpaceDim_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().OperatorSpaceDim(N: 0));
    }

    [Fact]
    public void XorDrainCount_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().XorDrainCount(N: 0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F23XorDrainVanishingFractionPi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
    }

    [Fact]
    public void Constructor_NullMirror_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F23XorDrainVanishingFractionPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
