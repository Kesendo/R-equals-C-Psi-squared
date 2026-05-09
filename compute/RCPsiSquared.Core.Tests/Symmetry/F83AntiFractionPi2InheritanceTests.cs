using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F83AntiFractionPi2InheritanceTests
{
    private static F83AntiFractionPi2Inheritance BuildClaim() =>
        new F83AntiFractionPi2Inheritance(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void MaximumAntiFraction_IsExactlyOneHalf_FromBilinearApex()
    {
        // F83's maximum anti-fraction = 1/2 at r=0 (pure Π²-odd).
        // = a_2 on dyadic ladder = BilinearApex apex value.
        Assert.Equal(0.5, BuildClaim().MaximumAntiFraction, precision: 14);
    }

    [Fact]
    public void DenominatorConstantCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().DenominatorConstantCoefficient, precision: 14);
    }

    [Fact]
    public void DenominatorLinearCoefficient_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().DenominatorLinearCoefficient, precision: 14);
    }

    [Fact]
    public void MNormCoefficientForOdd_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().MNormCoefficientForOdd, precision: 14);
    }

    [Fact]
    public void MNormCoefficientForEvenNontruly_IsExactlyEight()
    {
        // 8 = a_{-2} on the dyadic ladder; encodes n_YZ=2 for Π²-even non-truly
        // via F49's 2^(N+2)·n_YZ mechanism.
        Assert.Equal(8.0, BuildClaim().MNormCoefficientForEvenNontruly, precision: 14);
    }

    [Theory]
    // 2^N = a_{1-N} on the dyadic ladder.
    // N=1: a_0 = 2
    // N=2: a_{-1} = 4
    // N=3: a_{-2} = 8
    [InlineData(1, 2.0)]
    [InlineData(2, 4.0)]
    [InlineData(3, 8.0)]
    [InlineData(4, 16.0)]
    [InlineData(5, 32.0)]
    public void HilbertSpaceDimension_IsTwoPowerN(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().HilbertSpaceDimension(N), precision: 12);
    }

    [Theory]
    // F83 verified table from ANALYTICAL_FORMULAS:
    // r = 0       → 1/2 (50/50 split, F81 Step 8)
    // r = 1/2     → 1/4 (asymmetric)
    // r = 1       → 1/6 (equal Frobenius mix)
    // r → ∞       → 0
    [InlineData(0.0, 1.0 / 2.0)]
    [InlineData(0.5, 1.0 / 4.0)]
    [InlineData(1.0, 1.0 / 6.0)]
    [InlineData(2.0, 1.0 / 10.0)]
    [InlineData(10.0, 1.0 / 42.0)]
    [InlineData(100.0, 1.0 / 402.0)]
    public void AntiFraction_MatchesClosedForm(double r, double expected)
    {
        Assert.Equal(expected, BuildClaim().AntiFraction(r), precision: 12);
    }

    [Fact]
    public void IsAtBilinearApex_TrueOnlyAtRZero()
    {
        var f = BuildClaim();
        Assert.True(f.IsAtBilinearApex(0.0));
        Assert.False(f.IsAtBilinearApex(0.1));
        Assert.False(f.IsAtBilinearApex(0.5));
    }

    [Fact]
    public void RAtQuarterCrossover_IsExactlyOneHalf()
    {
        Assert.Equal(0.5, BuildClaim().RAtQuarterCrossover, precision: 14);
    }

    [Fact]
    public void QuarterCrossoverHolds_BitExactly()
    {
        // At r = 1/2, anti-fraction = 1/4 = a_3 = QuarterAsBilinearMaxval.
        // Cross-anchor between BilinearApex (max) and QuarterAsBilinearMaxval.
        Assert.True(BuildClaim().QuarterCrossoverHolds());
    }

    [Theory]
    // From F83 verified instances table at N=3, J=1, γ_z=0:
    // |   H              | ‖H_odd‖² | ‖H_even‖² | predicted ‖M‖² |
    // |   XY+YX (pure)   |    32    |     0     |      1024      |
    // |   YZ+ZY (pure)   |     0    |    32     |      2048      |
    // |   XY+YZ (mixed)  |    16    |    16     |      1536      |
    // |   XY+YX+YZ       |    32    |    16     |      2048      |
    // (At N=3: 2^N = 8; coefficients are 4, 8)
    [InlineData(32.0, 0.0, 3, 1024.0)]    // 4·32·8 = 1024
    [InlineData(0.0, 32.0, 3, 2048.0)]    // 8·32·8 = 2048
    [InlineData(16.0, 16.0, 3, 1536.0)]   // 4·16·8 + 8·16·8 = 512 + 1024 = 1536
    [InlineData(32.0, 16.0, 3, 2048.0)]   // 4·32·8 + 8·16·8 = 1024 + 1024 = 2048
    public void MNormSquared_MatchesF83VerifiedTable(double hOddSq, double hEvenSq, int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().MNormSquared(hOddSq, hEvenSq, N), precision: 8);
    }

    [Theory]
    [InlineData(32.0, 3, 512.0)]    // 2·32·8 = 512
    [InlineData(16.0, 3, 256.0)]    // 2·16·8 = 256
    public void MAntiNormSquared_MatchesF83Mechanism(double hOddSq, int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().MAntiNormSquared(hOddSq, N), precision: 8);
    }

    [Fact]
    public void AntiFraction_NegativeR_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AntiFraction(-0.1));
    }

    [Fact]
    public void HilbertSpaceDimension_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().HilbertSpaceDimension(0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F83AntiFractionPi2Inheritance(null!));
    }
}
