using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F97CardioidHalfFixedPointPi2InheritanceTests
{
    private static F97CardioidHalfFixedPointPi2Inheritance BuildClaim() => new();

    [Theory]
    [InlineData(0.0)]
    [InlineData(1.5707963267948966)]
    [InlineData(3.141592653589793)]
    [InlineData(6.283185307178586)]
    public void SelectedPeriodOneRootIsMarginalAndSatisfiesQuadratic(double phi)
    {
        var f = BuildClaim();
        Assert.True(f.SelectedMultiplierIsMarginal(phi));
        Assert.True(f.AlgebraicIdentityHolds(phi));
        Assert.Equal(0.5, f.FixedPointMagnitude(phi), precision: 14);
    }

    [Fact]
    public void ParametrizationPinsCuspAndTail()
    {
        var f = BuildClaim();
        Assert.True(f.CuspAgreesWithF95Threshold());
        Assert.True(f.TailAtMinusThreeQuarters());
        var atQuarterTurn = f.CardioidC(Math.PI / 2);
        Assert.Equal(0.25, atQuarterTurn.Real, precision: 14);
        Assert.Equal(0.5, atQuarterTurn.Imaginary, precision: 14);
    }

    [Fact]
    public void OtherRootIsNotGenerallyMarginal()
    {
        var f = BuildClaim();
        Complex other = f.OtherFixedPoint(Math.PI / 2);
        Assert.Equal(Math.Sqrt(5.0), F97CardioidHalfFixedPointPi2Inheritance.MultiplierMagnitude(other), precision: 14);
        Assert.NotEqual(1.0, F97CardioidHalfFixedPointPi2Inheritance.MultiplierMagnitude(other));
    }

    [Fact]
    public void OriginHasTwoDistinctFixedPoints()
    {
        Assert.True(BuildClaim().OriginRootsAreDistinct());
    }

    [Fact]
    public void ArgumentUsesWrappedPrincipalBranch()
    {
        var f = BuildClaim();
        Assert.Equal(Math.PI / 2, f.FixedPointArgument(Math.PI / 2), precision: 14);
        Assert.Equal(-Math.PI / 2, f.FixedPointArgument(3 * Math.PI / 2), precision: 14);
        Assert.Equal(-1e-6, f.FixedPointArgument(2 * Math.PI - 1e-6), precision: 10);
    }

    [Fact]
    public void ClaimIsParentlessAndHasNoPolarityAxis()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F97CardioidHalfFixedPointPi2Inheritance)));
        var parents = typeof(F97CardioidHalfFixedPointPi2Inheritance)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
        Assert.NotNull(typeof(F97CardioidHalfFixedPointPi2Inheritance).GetConstructor(Type.EmptyTypes));
    }
}
