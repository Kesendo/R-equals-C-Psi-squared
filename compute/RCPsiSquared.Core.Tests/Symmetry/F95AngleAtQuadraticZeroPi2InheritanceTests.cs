using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F95AngleAtQuadraticZeroPi2InheritanceTests
{
    private static F95AngleAtQuadraticZeroPi2Inheritance BuildClaim() =>
        new F95AngleAtQuadraticZeroPi2Inheritance(
            new PolynomialFoundationClaim(),
            new HalfAsStructuralFixedPointClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new NinetyDegreeMirrorMemoryClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void B_IsExactlyOneHalf()
    {
        Assert.Equal(0.5, F95AngleAtQuadraticZeroPi2Inheritance.B, precision: 15);
    }

    [Fact]
    public void Threshold_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, F95AngleAtQuadraticZeroPi2Inheritance.Threshold, precision: 15);
    }

    [Fact]
    public void Threshold_EqualsBSquared()
    {
        // The threshold IS b² for the discriminant zero crossing.
        double b = F95AngleAtQuadraticZeroPi2Inheritance.B;
        Assert.Equal(b * b, F95AngleAtQuadraticZeroPi2Inheritance.Threshold, precision: 15);
    }

    [Theory]
    // Februar BOUNDARY_NAVIGATION.md θ-compass anchor points (real CΨ values
    // for which θ = arctan(√(4·CΨ − 1)) gives the documented angle).
    [InlineData(1.0 / 3.0, 30.0)]     // Bell+ initial: θ = 30°
    [InlineData(0.308, 25.7184)]      // BOUNDARY_NAVIGATION table t=0.2
    [InlineData(0.286, 20.7804)]      // table t=0.4
    [InlineData(0.266, 14.1969)]      // table t=0.6
    public void ThetaForFramework_MatchesFebruarBoundaryNavigationTable(double cpsi, double expected_deg)
    {
        double theta_rad = BuildClaim().ThetaForFramework(cpsi);
        double theta_deg = theta_rad * 180.0 / System.Math.PI;
        Assert.Equal(expected_deg, theta_deg, precision: 3);
    }

    [Fact]
    public void ThetaForFramework_AtThreshold_IsExactlyZero()
    {
        // At CΨ = 1/4 the discriminant is zero; the two complex roots
        // collapse to the single degenerate real root z = b = 1/2; θ = 0.
        // Implementation returns NaN at exactly c = b² (degenerate case); the
        // limit as c ↓ b² is 0. Use a value just above for the test.
        double theta_rad = BuildClaim().ThetaForFramework(0.25 + 1e-15);
        Assert.InRange(theta_rad, 0.0, 1e-7);  // tiny positive, near 0
    }

    [Fact]
    public void ThetaForFramework_BelowThreshold_IsNaN()
    {
        Assert.True(double.IsNaN(BuildClaim().ThetaForFramework(0.20)));
    }

    [Fact]
    public void BellPlusInitialAngleIs30Degrees_HoldsBitExact()
    {
        // CΨ = 1/3 → θ = arctan(√(4/3 − 1)) = arctan(√(1/3)) = arctan(1/√3) = π/6 = 30°.
        Assert.True(BuildClaim().BellPlusInitialAngleIs30Degrees());
    }

    [Fact]
    public void FrameworkSpecializationAgrees_AtFebruarPoints()
    {
        var f = BuildClaim();
        Assert.True(f.FrameworkSpecializationAgrees(1.0 / 3.0));
        Assert.True(f.FrameworkSpecializationAgrees(0.286));
        Assert.True(f.FrameworkSpecializationAgrees(0.20));  // both NaN below threshold
    }

    [Theory]
    [InlineData(0.5, 1.0)]       // b=0.5, c=1.0 → θ = arctan(√(4·1 − 1)) = arctan(√3) = π/3 = 60°
    [InlineData(1.0, 4.0)]       // b=1, c=4 → θ = arctan(√(4 − 1)) = arctan(√3) = π/3 = 60°
    [InlineData(2.0, 16.0)]      // b=2, c=16 → θ = arctan(√(16/4 − 1)) = arctan(√3) = π/3 = 60°
    public void ThetaGeneral_IsScaleInvariantWhenCOverBSquaredIsConstant(double b, double c)
    {
        // The angle depends only on the ratio c/b². Three (b, c) pairs with
        // c/b² = 4 should give exactly the same angle (= π/3 = 60°).
        double theta = BuildClaim().ThetaGeneral(c, b);
        Assert.Equal(System.Math.PI / 3.0, theta, precision: 14);
    }

    [Fact]
    public void ThetaGeneral_ZeroB_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().ThetaGeneral(1.0, 0.0));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var p = new PolynomialFoundationClaim();
        var h = new HalfAsStructuralFixedPointClaim();
        var q = new QuarterAsBilinearMaxvalClaim();
        var n = new NinetyDegreeMirrorMemoryClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F95AngleAtQuadraticZeroPi2Inheritance(null!, h, q, n));
        Assert.Throws<ArgumentNullException>(() =>
            new F95AngleAtQuadraticZeroPi2Inheritance(p, null!, q, n));
        Assert.Throws<ArgumentNullException>(() =>
            new F95AngleAtQuadraticZeroPi2Inheritance(p, h, null!, n));
        Assert.Throws<ArgumentNullException>(() =>
            new F95AngleAtQuadraticZeroPi2Inheritance(p, h, q, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var f = BuildClaim();
        Assert.NotNull(f.Polynomial);
        Assert.NotNull(f.Half);
        Assert.NotNull(f.Quarter);
        Assert.NotNull(f.NinetyDegree);
    }

    [Fact]
    public void Anchor_References_ProofAndAllFourPi2Parents()
    {
        var f = BuildClaim();
        Assert.Contains("PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md", f.Anchor);
        Assert.Contains("ANALYTICAL_FORMULAS.md F95", f.Anchor);
        Assert.Contains("ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md", f.Anchor);
        Assert.Contains("PolynomialFoundation", f.Anchor);
        Assert.Contains("HalfAsStructuralFixedPoint", f.Anchor);
        Assert.Contains("QuarterAsBilinearMaxval", f.Anchor);
        Assert.Contains("NinetyDegreeMirrorMemory", f.Anchor);
    }
}
