using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F33ExactN3DecayRatesPi2InheritanceTests
{
    private static F33ExactN3DecayRatesPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f50 = new F50WeightOneDegeneracyPi2Inheritance(ladder);
        return new F33ExactN3DecayRatesPi2Inheritance(ladder, f50);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void WeightOneRateCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().WeightOneRateCoefficient, precision: 14);
    }

    [Theory]
    [InlineData(0.05, 0.1)]
    [InlineData(0.5, 1.0)]
    [InlineData(1.0, 2.0)]
    public void Rate1_Equals2Gamma(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().Rate1(gammaZero), precision: 14);
    }

    [Theory]
    [InlineData(0.05, 8.0 * 0.05 / 3.0)]
    [InlineData(1.0, 8.0 / 3.0)]
    [InlineData(3.0, 8.0)]   // 8·3/3 = 8
    public void Rate2_Equals8GammaOverThree(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().Rate2(gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(0.05, 10.0 * 0.05 / 3.0)]
    [InlineData(1.0, 10.0 / 3.0)]
    [InlineData(3.0, 10.0)]   // 10·3/3 = 10
    public void Rate3_Equals10GammaOverThree(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().Rate3(gammaZero), precision: 12);
    }

    [Fact]
    public void Rate2_DivByRate3_IsExactlyFourFifths()
    {
        // 8γ/3 ÷ 10γ/3 = 8/10 = 4/5 = 0.8 exactly.
        var f = BuildClaim();
        double ratio = f.Rate2(1.0) / f.Rate3(1.0);
        Assert.Equal(0.8, ratio, precision: 14);
    }

    [Fact]
    public void RatesAreRationalRatios_HoldsAtAnyGamma()
    {
        var f = BuildClaim();
        Assert.True(f.RatesAreRationalRatios(0.05));
        Assert.True(f.RatesAreRationalRatios(1.0));
        Assert.True(f.RatesAreRationalRatios(0.001));
    }

    [Fact]
    public void XorBoundaryRate_Equals4Gamma()
    {
        // 2(N-1)γ at N=3 = 4γ.
        Assert.Equal(0.2, BuildClaim().XorBoundaryRate(0.05), precision: 14);
        Assert.Equal(4.0, BuildClaim().XorBoundaryRate(1.0), precision: 14);
    }

    [Theory]
    // Absorption Theorem: ⟨n_XY⟩ = rate / (2γ)
    // rate_1 = 2γ at γ=0.05 → 0.1; ⟨n_XY⟩ = 0.1/0.1 = 1
    // rate_2 = 8γ/3 at γ=0.05 → 0.4/3 ≈ 0.1333; ⟨n_XY⟩ = (0.4/3)/0.1 = 4/3
    // rate_3 = 10γ/3 at γ=0.05 → 0.5/3 ≈ 0.1667; ⟨n_XY⟩ = (0.5/3)/0.1 = 5/3
    [InlineData(0.1, 0.05, 1.0)]
    [InlineData(0.4 / 3.0, 0.05, 4.0 / 3.0)]
    [InlineData(0.5 / 3.0, 0.05, 5.0 / 3.0)]
    public void NXyExpectationFromRate_MatchesAbsorptionTheorem(double rate, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().NXyExpectationFromRate(rate, gammaZero), precision: 12);
    }

    [Fact]
    public void AbsorptionTheoremTable_HasThreeEntries()
    {
        var table = BuildClaim().AbsorptionTheoremTable;
        Assert.Equal(3, table.Count);
        Assert.Equal((1, 1.0), table[0]);
        Assert.Equal((2, 4.0 / 3.0), table[1]);
        Assert.Equal((3, 5.0 / 3.0), table[2]);
    }

    [Fact]
    public void Rate1_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Rate1(-0.05));
    }

    [Fact]
    public void NXyExpectationFromRate_NonPositiveGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().NXyExpectationFromRate(0.1, gammaZero: 0.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f50 = new F50WeightOneDegeneracyPi2Inheritance(ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F33ExactN3DecayRatesPi2Inheritance(null!, f50));
    }

    [Fact]
    public void Constructor_NullF50_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F33ExactN3DecayRatesPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
