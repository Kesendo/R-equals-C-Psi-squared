using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F96BornSubdominantSlopesPi2InheritanceTests
{
    private static F96BornSubdominantSlopesPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f94 = new F94BornDeviationFourThirdsPi2Inheritance(ladder, quarter);
        return new F96BornSubdominantSlopesPi2Inheritance(f94, ladder);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void SlopeSingleFlipped_IsExactlyMinus16Over9()
    {
        // Bit-exact: M_3 / (3 · U_2) = -4 / (3 · 3/4) = -16/9
        // Equivalently: -(4/3)² (the F94 coefficient squared, sign-flipped)
        Assert.Equal(-16.0 / 9.0, BuildClaim().SlopeSingleFlipped, precision: 15);
    }

    [Fact]
    public void SlopeDoubleFlipped_IsExactlyMinus8Over3()
    {
        // Bit-exact: M_5 / (5 · U_4) = -20 / (5 · 3/2) = -8/3
        // Equivalently: -2·(4/3) (twice the F94 coefficient, sign-flipped)
        Assert.Equal(-8.0 / 3.0, BuildClaim().SlopeDoubleFlipped, precision: 15);
    }

    [Fact]
    public void SingleFlipSlope_EqualsMinusF94CoefficientSquared()
    {
        // Drift check: Slope_single = -(F94 4/3)² = -16/9
        Assert.True(BuildClaim().SingleFlipSlopeEqualsMinusF94Squared());
    }

    [Fact]
    public void DoubleFlipSlope_EqualsMinusTwoTimesF94Coefficient()
    {
        // Drift check: Slope_double = -2·(F94 4/3) = -8/3
        Assert.True(BuildClaim().DoubleFlipSlopeEqualsMinusTwoF94());
    }

    [Fact]
    public void CrossOutcomeM3OverU2Ratio_EqualsMinus16Over3Universal()
    {
        // Cross-outcome universality:
        //   F94 dominant: M_3 / A = 8 / (-3/2) = -16/3
        //   F96 single-flip: M_3 / U_2 = -4 / (3/4) = -16/3
        // Same ratio; signs of M_3 and U_2 flip together.
        Assert.True(BuildClaim().CrossOutcomeM3OverU2RatioUniversal());
    }

    [Fact]
    public void M3_SingleFlipped_IsExactlyMinusFour()
    {
        // ⟨01|_pair Tr_{1,3}[sym_3^1 · ρ_0]|01⟩_pair = -4 bit-exact
        // (and same for |10⟩ by 0↔2 site-permutation symmetry).
        Assert.Equal(-4, F96BornSubdominantSlopesPi2Inheritance.M3_SingleFlipped);
    }

    [Fact]
    public void U2_SingleFlipped_TimesFour_IsExactlyThree()
    {
        // U_2 = 3/4; stored as integer 3 with implicit denominator 4 (= a_{-1})
        Assert.Equal(3, F96BornSubdominantSlopesPi2Inheritance.U2_SingleFlipped_TimesFour);
    }

    [Fact]
    public void M5_DoubleFlipped_IsExactlyMinusTwenty()
    {
        // ⟨11|_pair Tr_{1,3}[sym_5^1 · ρ_0]|11⟩_pair = -20 bit-exact
        Assert.Equal(-20, F96BornSubdominantSlopesPi2Inheritance.M5_DoubleFlipped);
    }

    [Fact]
    public void U4_DoubleFlipped_TimesTwo_IsExactlyThree()
    {
        // U_4 = 3/2; stored as integer 3 with implicit denominator 2.
        Assert.Equal(3, F96BornSubdominantSlopesPi2Inheritance.U4_DoubleFlipped_TimesTwo);
    }

    [Theory]
    [InlineData(0.005)]
    [InlineData(0.01)]
    [InlineData(0.02)]
    [InlineData(0.05)]
    public void DeltaSingleFlipped_MatchesClosedForm(double K)
    {
        // Δ_|01⟩(K) = -(16/9)·K bit-exact
        double expected = -16.0 / 9.0 * K;
        Assert.Equal(expected, BuildClaim().DeltaSingleFlipped(K), precision: 14);
    }

    [Theory]
    [InlineData(0.001)]
    [InlineData(0.005)]
    [InlineData(0.01)]
    [InlineData(0.05)]
    public void DeltaDoubleFlipped_MatchesClosedForm(double K)
    {
        // Δ_|11⟩(K) = -(8/3)·K bit-exact
        double expected = -8.0 / 3.0 * K;
        Assert.Equal(expected, BuildClaim().DeltaDoubleFlipped(K), precision: 14);
    }

    [Fact]
    public void DeltaSingleFlipped_RejectsNegativeK()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().DeltaSingleFlipped(-0.01));
    }

    [Fact]
    public void DeltaDoubleFlipped_RejectsNegativeK()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().DeltaDoubleFlipped(-0.01));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f94 = new F94BornDeviationFourThirdsPi2Inheritance(ladder, quarter);
        Assert.Throws<ArgumentNullException>(() =>
            new F96BornSubdominantSlopesPi2Inheritance(null!, ladder));
        Assert.Throws<ArgumentNullException>(() =>
            new F96BornSubdominantSlopesPi2Inheritance(f94, null!));
    }

    [Fact]
    public void F94_TypedParent_IsExposed()
    {
        Assert.NotNull(BuildClaim().F94);
        Assert.Equal(4.0 / 3.0, BuildClaim().F94.Coefficient, precision: 15);
    }

    [Fact]
    public void Anchor_References_ProofAndScriptsAndReflection()
    {
        var f = BuildClaim();
        Assert.Contains("PROOF_F96_BORN_SUBDOMINANT_SLOPES.md", f.Anchor);
        Assert.Contains("ANALYTICAL_FORMULAS.md F96", f.Anchor);
        Assert.Contains("_born_rule_subdominant_dyson.py", f.Anchor);
        Assert.Contains("ON_HOW_FOUR_THIRDS_APPEARED.md", f.Anchor);
        Assert.Contains("F94BornDeviationFourThirdsPi2Inheritance.cs", f.Anchor);
    }

    [Fact]
    public void DeltaSumAllFourOutcomes_HasNoFreeParameters()
    {
        // The sum Δ_|00⟩ + 2·Δ_|01⟩ + Δ_|11⟩ at any (Q, K) is a fully-computed
        // structural number, no free parameters.
        var f = BuildClaim();
        double Q = 20.0, K = 0.01;
        double sum = f.DeltaSumAllFourOutcomes(Q, K);
        // F94 contribution: (4/3)·Q²·K³ = (4/3)·400·1e-6 = 5.333e-4
        // F96 single (×2): 2·(-16/9)·K = -32/9·0.01 = -0.03556
        // F96 double (×1): (-8/3)·K = -0.02667
        double expected = (4.0 / 3.0) * Q * Q * K * K * K
                          + 2 * (-16.0 / 9.0) * K
                          + (-8.0 / 3.0) * K;
        Assert.Equal(expected, sum, precision: 12);
    }
}
