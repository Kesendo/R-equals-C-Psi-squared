using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F95AngleAtQuadraticZeroPi2InheritanceTests
{
    private static F95AngleAtQuadraticZeroPi2Inheritance BuildClaim() => new();

    [Fact]
    public void PositiveBFormulaHasCorrectBranches()
    {
        var f = BuildClaim();
        Assert.True(double.IsNaN(f.ThetaGeneral(0.9, 1.0)));
        Assert.Equal(0.0, f.ThetaGeneral(1.0, 1.0));
        Assert.Equal(Math.PI / 6.0, f.ThetaGeneral(1.0 / 3.0, 0.5), precision: 14);
        Assert.Equal(Math.Atan(0.5), f.ThetaGeneral(5.0, 2.0), precision: 14);
    }

    [Theory]
    [InlineData(0.5, 1.0)]    // c/b² = 4: θ = arctan(√3) = 60°
    [InlineData(1.0, 4.0)]
    [InlineData(2.0, 16.0)]
    public void FormulaDependsOnlyOnCOverBSquared(double b, double c)
    {
        Assert.Equal(Math.PI / 3.0, BuildClaim().ThetaGeneral(c, b), precision: 14);
    }

    // experiments/BOUNDARY_NAVIGATION.md's trajectory table prints CΨ to three decimals and θ to
    // one; θ steepens toward ¼, so only the rows away from the boundary are reproducible from the
    // printed CΨ within the 0.1° of the printed θ (the t = 0.7 row, CΨ = 0.254, is not).
    [Theory]
    [InlineData(1.0 / 3.0, 30.0)]    // t = 0.0
    [InlineData(0.308, 25.7)]        // t = 0.2
    [InlineData(0.285, 20.5)]        // t = 0.4
    [InlineData(0.264, 13.4)]        // t = 0.6
    public void BoundaryNavigationTableRowsAreReproducedFromTheirPrintedCPsi(double c, double degrees)
    {
        double theta = BuildClaim().ThetaForFramework(c) * 180.0 / Math.PI;
        Assert.True(Math.Abs(theta - degrees) <= 0.1, $"c={c}: {theta} vs {degrees}");
    }

    [Fact]
    public void BoundaryNavigationRealRootRowHasNoAngle()
    {
        Assert.True(double.IsNaN(BuildClaim().ThetaForFramework(0.245)));   // t = 0.8, real roots
    }

    [Fact]
    public void FrameworkBranchIsUndefinedBelowTheQuarter()
    {
        var f = BuildClaim();
        Assert.True(double.IsNaN(f.ThetaForFramework(0.2)));
        Assert.True(f.FrameworkSpecializationAgrees(0.2));
    }

    [Fact]
    public void FrameworkSpecializationIsOneQuarterPoint()
    {
        var f = BuildClaim();
        Assert.Equal(0.5, F95AngleAtQuadraticZeroPi2Inheritance.B);
        Assert.Equal(0.25, F95AngleAtQuadraticZeroPi2Inheritance.Threshold);
        Assert.Equal(0.0, f.ThetaForFramework(0.25));
        Assert.True(f.FrameworkSpecializationAgrees(0.286));
        Assert.True(f.BellPlusInitialAngleIs30Degrees());
    }

    [Theory]
    [InlineData(1.0, 0.0)]
    [InlineData(1.0, -1.0)]
    [InlineData(1.0, double.PositiveInfinity)]
    [InlineData(1.0, double.NegativeInfinity)]
    [InlineData(1.0, double.NaN)]
    [InlineData(double.NaN, 1.0)]
    public void NonFiniteCOrNonPositiveBIsRejected(double c, double b)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().ThetaGeneral(c, b));
    }

    [Fact]
    public void TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void ClaimIsParentlessAndHasNoPolarityAxis()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F95AngleAtQuadraticZeroPi2Inheritance)));
        var parents = typeof(F95AngleAtQuadraticZeroPi2Inheritance)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
        Assert.NotNull(typeof(F95AngleAtQuadraticZeroPi2Inheritance).GetConstructor(Type.EmptyTypes));
    }
}
