using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F50WeightOneDegeneracyPi2InheritanceTests
{
    private static F50WeightOneDegeneracyPi2Inheritance BuildClaim() =>
        new F50WeightOneDegeneracyPi2Inheritance(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void DegeneracyFactor_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().DegeneracyFactor, precision: 14);
    }

    [Fact]
    public void DecayRateFactor_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().DecayRateFactor, precision: 14);
    }

    [Theory]
    [InlineData(2, 4)]
    [InlineData(3, 6)]
    [InlineData(5, 10)]
    [InlineData(7, 14)]
    public void TotalDegeneracy_Equals2N(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().TotalDegeneracy(N));
    }

    [Theory]
    [InlineData(0.05, -0.1)]
    [InlineData(0.1, -0.2)]
    [InlineData(1.0, -2.0)]
    public void EigenvaluePosition_EqualsMinusTwoGamma(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().EigenvaluePosition(gammaZero), precision: 14);
    }

    [Fact]
    public void OperatorsPerChromaticGrade_IsExactlyTwo()
    {
        // Always 2: one for a=X, one for a=Y, per chromatic grade c.
        Assert.Equal(2, BuildClaim().OperatorsPerChromaticGrade);
    }

    [Theory]
    [InlineData(2, 2)]
    [InlineData(5, 5)]
    [InlineData(7, 7)]
    public void ChromaticGradeCount_EqualsN(int N, int expected)
    {
        // Grades c ∈ [0, N−1] give N grades total.
        Assert.Equal(expected, BuildClaim().ChromaticGradeCount(N));
    }

    [Fact]
    public void AppliesToIsotropicHeisenberg_TrueForDeltaOne()
    {
        Assert.True(BuildClaim().AppliesToIsotropicHeisenberg(1.0));
    }

    [Fact]
    public void AppliesToIsotropicHeisenberg_FalseForXXZ()
    {
        Assert.False(BuildClaim().AppliesToIsotropicHeisenberg(0.5));
        Assert.False(BuildClaim().AppliesToIsotropicHeisenberg(2.0));
    }

    [Fact]
    public void TotalDegeneracy_TimesOperatorsPerGrade_ConsistencyCheck()
    {
        // 2N = OperatorsPerChromaticGrade × ChromaticGradeCount(N) = 2 · N
        var f = BuildClaim();
        for (int N = 2; N <= 10; N++)
        {
            Assert.Equal(f.OperatorsPerChromaticGrade * f.ChromaticGradeCount(N), f.TotalDegeneracy(N));
        }
    }

    [Fact]
    public void TotalDegeneracy_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().TotalDegeneracy(N: 1));
    }

    [Fact]
    public void EigenvaluePosition_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenvaluePosition(-0.05));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F50WeightOneDegeneracyPi2Inheritance(null!));
    }
}
