using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F85KBodyFChainPi2InheritanceTests
{
    private static F85KBodyFChainPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var f49 = new F49Pi2Inheritance(ladder, mirror);
        var f83 = new F83AntiFractionPi2Inheritance(ladder);
        return new F85KBodyFChainPi2Inheritance(ladder, f49, f83);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void FactorPi2Odd_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().FactorPi2Odd, precision: 14);
    }

    [Fact]
    public void FactorPi2EvenNonTruly_IsExactlyEight()
    {
        Assert.Equal(8.0, BuildClaim().FactorPi2EvenNonTruly, precision: 14);
    }

    [Fact]
    public void OddCountDenominator_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().OddCountDenominator, precision: 14);
    }

    [Theory]
    [InlineData("truly", 0)]
    [InlineData("pi2_odd", 1)]
    [InlineData("pi2_odd_nontruly", 1)]
    [InlineData("pi2_even_nontruly", 2)]
    public void FrobeniusFactorPerClass_MatchesC(string pi2Class, int expected)
    {
        Assert.Equal(expected, BuildClaim().FrobeniusFactorPerClass(pi2Class));
    }

    [Theory]
    // F85 verified Π²-odd counts (3^k − (−1)^k) / 2:
    // k=1: (3 − (−1))/2 = 2
    // k=2: (9 − 1)/2 = 4
    // k=3: (27 + 1)/2 = 14
    // k=4: (81 − 1)/2 = 40
    [InlineData(1, 2)]
    [InlineData(2, 4)]
    [InlineData(3, 14)]
    [InlineData(4, 40)]
    public void Pi2OddCount_MatchesF85VerifiedTable(int k, long expected)
    {
        Assert.Equal(expected, BuildClaim().Pi2OddCount(k));
    }

    [Theory]
    // F85 trichotomy classifier:
    // XX: truly (#Y=0 even, #Z=0 even)
    // XY: Π²-odd (#Y=1, #Z=0; bit_b=1)
    // YY: truly (#Y=2 even, #Z=0 even)
    // YZ: Π²-even non-truly (#Y=1, #Z=1; bit_b=0; not truly because Y odd)
    // ZZ: truly (#Y=0, #Z=2)
    // YYY: Π²-odd (#Y=3 odd, but bit_b = (3+0) mod 2 = 1)
    [InlineData("XX", "truly")]
    [InlineData("XY", "pi2_odd")]
    [InlineData("YY", "truly")]
    [InlineData("YZ", "pi2_even_nontruly")]
    [InlineData("ZZ", "truly")]
    [InlineData("YYY", "pi2_odd")]
    [InlineData("XYZ", "pi2_even_nontruly")]   // #Y=1, #Z=1, bit_b=0, not truly
    [InlineData("XYY", "truly")]                // #Y=2 even, #Z=0 even → truly
    public void Pi2ClassFromLetters_ClassifiesCorrectly(string letters, string expected)
    {
        // Recompute expected for some edge cases
        var letterArray = letters.ToCharArray();
        Assert.Equal(expected, BuildClaim().Pi2ClassFromLetters(letterArray));
    }

    [Theory]
    [InlineData("truly", 3, 0.0)]                // c=0 → coefficient 0
    [InlineData("pi2_odd", 3, 32.0)]             // 4·1·8 = 32
    [InlineData("pi2_even_nontruly", 3, 64.0)]   // 4·2·8 = 64
    [InlineData("pi2_odd", 4, 64.0)]             // 4·1·16
    public void TotalCoefficientForClass_MatchesClosedForm(string pi2Class, int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().TotalCoefficientForClass(pi2Class, N), precision: 12);
    }

    [Fact]
    public void MatchesF83Coefficients_HoldsExactly()
    {
        Assert.True(BuildClaim().MatchesF83Coefficients());
    }

    [Fact]
    public void F49IsBaseCase_HoldsTrue()
    {
        Assert.True(BuildClaim().F49IsBaseCase());
    }

    [Fact]
    public void Pi2OddCount_KLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Pi2OddCount(0));
    }

    [Fact]
    public void FrobeniusFactorPerClass_UnknownClass_Throws()
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().FrobeniusFactorPerClass("unknown"));
    }

    [Fact]
    public void Pi2ClassFromLetters_InvalidLetter_Throws()
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().Pi2ClassFromLetters("AB".ToCharArray()));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var f49 = new F49Pi2Inheritance(ladder, mirror);
        var f83 = new F83AntiFractionPi2Inheritance(ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F85KBodyFChainPi2Inheritance(null!, f49, f83));
    }
}
