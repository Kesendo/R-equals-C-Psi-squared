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

    [Fact]
    public void FormulaIsScaleInvariantForPositiveScaling()
    {
        var f = BuildClaim();
        double theta = f.ThetaGeneral(5.0, 2.0);
        Assert.Equal(theta, f.ThetaGeneral(45.0, 6.0), precision: 14);
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
    [InlineData(double.NaN, 1.0)]
    public void NonFiniteCOrNonPositiveBIsRejected(double c, double b)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().ThetaGeneral(c, b));
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
