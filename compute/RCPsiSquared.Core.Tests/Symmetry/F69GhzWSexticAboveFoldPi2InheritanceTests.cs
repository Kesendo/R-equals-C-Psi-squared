using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F69GhzWSexticAboveFoldPi2InheritanceTests
{
    private static F69GhzWSexticAboveFoldPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f60 = new F60GhzBornBelowFoldPi2Inheritance(
            ladder,
            new PolarityLayerOriginClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new ArgmaxMaxvalPairClaim());
        var f62 = new F62WStateBornBelowFoldPi2Inheritance(
            ladder,
            new QuarterAsBilinearMaxvalClaim(),
            F62WStateBornBelowFoldPi2InheritanceTests.BuildF61(ladder));
        return new F69GhzWSexticAboveFoldPi2Inheritance(ladder, f60, f62, new QuarterAsBilinearMaxvalClaim());
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void FoldPosition_IsExactlyOneQuarter()
    {
        // 1/4 = a_3 on the dyadic ladder. Same anchor as F60/F62.
        Assert.Equal(0.25, BuildClaim().FoldPosition, precision: 14);
    }

    [Fact]
    public void OptimumPairCpsi_ExceedsOneQuarter()
    {
        // The structural F69 claim: the optimum sits ABOVE 1/4.
        Assert.True(F69GhzWSexticAboveFoldPi2Inheritance.OptimumPairCpsi > 0.25);
    }

    [Fact]
    public void LiftRatioAboveFold_IsApproximately1281()
    {
        // 0.320411541127025 / 0.25 ≈ 1.281646
        var f = BuildClaim();
        Assert.Equal(1.281646, f.LiftRatioAboveFold, precision: 5);
    }

    [Fact]
    public void OptimumExceedsFold_IsTrue()
    {
        Assert.True(BuildClaim().OptimumExceedsFold);
    }

    [Fact]
    public void GhzBaselineBelowFold_TrueAtN3()
    {
        // GHZ_3: CΨ(0) = 1/(2^3 − 1) = 1/7 ≈ 0.143 < 0.25.
        Assert.True(BuildClaim().GhzBaselineBelowFold(N: 3));
    }

    [Fact]
    public void WBaselineBelowFold_TrueAtN3()
    {
        // W_3: CΨ(0) = 10/81 ≈ 0.1235 < 0.25.
        Assert.True(BuildClaim().WBaselineBelowFold(N: 3));
    }

    [Fact]
    public void LiftExistsAtN_TrueOnlyForN3()
    {
        var f = BuildClaim();
        Assert.True(f.LiftExistsAtN(3));
        Assert.False(f.LiftExistsAtN(4));
        Assert.False(f.LiftExistsAtN(5));
        Assert.False(f.LiftExistsAtN(8));
    }

    [Fact]
    public void SexticCoefficients_MatchClosedForm()
    {
        // P(x) = 2900x⁶ − 8060x⁵ + 4211x⁴ + 3832x³ − 2428x² − 512x + 300.
        var coeffs = BuildClaim().SexticCoefficients;
        Assert.Equal(7, coeffs.Count);
        Assert.Equal(2900, coeffs[0]);
        Assert.Equal(-8060, coeffs[1]);
        Assert.Equal(4211, coeffs[2]);
        Assert.Equal(3832, coeffs[3]);
        Assert.Equal(-2428, coeffs[4]);
        Assert.Equal(-512, coeffs[5]);
        Assert.Equal(300, coeffs[6]);
    }

    [Fact]
    public void SexticOptimum_IsRoot()
    {
        // Verify P(α²_opt) ≈ 0 to confirm the optimum value is consistent with the sextic.
        // P(x) = 2900x⁶ − 8060x⁵ + 4211x⁴ + 3832x³ − 2428x² − 512x + 300
        double x = F69GhzWSexticAboveFoldPi2Inheritance.AlphaSquaredOptimum;
        double p = 2900 * Math.Pow(x, 6)
                 - 8060 * Math.Pow(x, 5)
                 + 4211 * Math.Pow(x, 4)
                 + 3832 * Math.Pow(x, 3)
                 - 2428 * Math.Pow(x, 2)
                 - 512 * x
                 + 300;
        // The constant α²_opt is given to 15-digit precision; allow O(10⁻⁹) residual
        // due to the rapid coefficient growth (sextic with coefficients up to 2900 amplifies
        // the input precision by ~3 orders).
        Assert.True(Math.Abs(p) < 1e-9, $"sextic P(α²_opt) should be ≈ 0 but got {p}");
    }

    [Fact]
    public void AlphaOptimum_IsSquareRootOfAlphaSquaredOptimum()
    {
        Assert.Equal(F69GhzWSexticAboveFoldPi2Inheritance.AlphaOptimum,
                     Math.Sqrt(F69GhzWSexticAboveFoldPi2Inheritance.AlphaSquaredOptimum),
                     precision: 12);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f60 = new F60GhzBornBelowFoldPi2Inheritance(
            ladder,
            new PolarityLayerOriginClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new ArgmaxMaxvalPairClaim());
        var f62 = new F62WStateBornBelowFoldPi2Inheritance(
            ladder,
            new QuarterAsBilinearMaxvalClaim(),
            F62WStateBornBelowFoldPi2InheritanceTests.BuildF61(ladder));
        var quarter = new QuarterAsBilinearMaxvalClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F69GhzWSexticAboveFoldPi2Inheritance(null!, f60, f62, quarter));
        Assert.Throws<ArgumentNullException>(() =>
            new F69GhzWSexticAboveFoldPi2Inheritance(ladder, null!, f62, quarter));
        Assert.Throws<ArgumentNullException>(() =>
            new F69GhzWSexticAboveFoldPi2Inheritance(ladder, f60, null!, quarter));
        Assert.Throws<ArgumentNullException>(() =>
            new F69GhzWSexticAboveFoldPi2Inheritance(ladder, f60, f62, null!));
    }

    [Fact]
    public void TypedParents_QuarterIsExposed()
    {
        Assert.NotNull(BuildClaim().Quarter);
    }
}
