using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F64CavityModeExposurePi2InheritanceTests
{
    private static F64CavityModeExposurePi2Inheritance BuildClaim() =>
        new F64CavityModeExposurePi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new QuarterAsBilinearMaxvalClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void EigenvalueConventionCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().EigenvalueConventionCoefficient, precision: 14);
    }

    [Fact]
    public void UniformJSpecialValue_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, BuildClaim().UniformJSpecialValue, precision: 14);
    }

    [Theory]
    [InlineData(0.05, 0.5, 0.025)]
    [InlineData(1.0, 0.25, 0.25)]
    [InlineData(0.1, 0.0, 0.0)]    // zero amplitude → no decoherence
    public void LorentzianHalfWidth_EqualsGammaTimesABSquared(double gammaB, double aBSquared, double expected)
    {
        Assert.Equal(expected, BuildClaim().LorentzianHalfWidth(gammaB, aBSquared), precision: 14);
    }

    [Theory]
    [InlineData(0.05, 0.5, 0.05)]   // α = 2·0.05·0.5 = 0.05
    [InlineData(1.0, 0.25, 0.5)]
    [InlineData(0.1, 1.0, 0.2)]     // |a_B|²=1 (full localisation) → α = 2γ
    public void LiouvillianDecayConstant_EqualsTwoGammaTimesABSquared(double gammaB, double aBSquared, double expected)
    {
        Assert.Equal(expected, BuildClaim().LiouvillianDecayConstant(gammaB, aBSquared), precision: 14);
    }

    [Theory]
    // N=3 g(r): zero mode r²/(r²+1) for r < 1/√2, bonding 1/(2(r²+1)) for r ≥ 1/√2.
    // r=0.5: 0.25/1.25 = 0.2 (zero mode, r=0.5 < 1/√2 ≈ 0.707)
    // r=1: 1/(2·2) = 0.25 (bonding mode)
    // r=2: 1/(2·5) = 0.1 (bonding mode)
    [InlineData(0.5, 0.2)]
    [InlineData(1.0, 0.25)]
    [InlineData(2.0, 0.1)]
    [InlineData(0.0, 0.0)]
    public void N3GValue_MatchesClosedForm(double r, double expected)
    {
        Assert.Equal(expected, BuildClaim().N3GValue(r), precision: 12);
    }

    [Fact]
    public void N3CrossoverRatio_IsOneOverSqrt2()
    {
        Assert.Equal(1.0 / Math.Sqrt(2.0), F64CavityModeExposurePi2Inheritance.N3CrossoverRatio, precision: 14);
    }

    [Fact]
    public void N3GValueAtCrossover_IsOneThird()
    {
        Assert.Equal(1.0 / 3.0, F64CavityModeExposurePi2Inheritance.N3GValueAtCrossover, precision: 14);
    }

    [Fact]
    public void N3GValueAtUniformJ_IsExactlyOneQuarter()
    {
        // F64 says g(r=1) = 1/4. Should match QuarterAsBilinearMaxval.
        Assert.Equal(0.25, BuildClaim().N3GValueAtUniformJ, precision: 14);
        Assert.Equal(0.25, BuildClaim().N3GValue(1.0), precision: 14);
    }

    [Fact]
    public void UniformJSpecialValueMatchesQuarter_IsTrue()
    {
        // Drift check: F64's g(1) = QuarterAsBilinearMaxval anchor.
        Assert.True(BuildClaim().UniformJSpecialValueMatchesQuarter());
    }

    [Fact]
    public void N3CrossoverContinuityHolds_BothBranchesGiveOneThird()
    {
        Assert.True(BuildClaim().N3CrossoverContinuityHolds());
    }

    [Fact]
    public void LorentzianHalfWidth_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().LorentzianHalfWidth(-0.05, 0.5));
    }

    [Fact]
    public void LorentzianHalfWidth_OutOfRangeAmplitude_Throws()
    {
        // |a_B|² ∈ [0, 1] (probability).
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().LorentzianHalfWidth(0.05, -0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().LorentzianHalfWidth(0.05, 1.5));
    }

    [Fact]
    public void N3GValue_NegativeR_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().N3GValue(-1.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F64CavityModeExposurePi2Inheritance(null!, new QuarterAsBilinearMaxvalClaim()));
    }

    [Fact]
    public void Constructor_NullQuarter_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F64CavityModeExposurePi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
