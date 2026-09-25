using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F96BornSubdominantSlopesPi2InheritanceTests
{
    private static F96BornSubdominantSlopesPi2Inheritance BuildClaim() => new();

    [Fact]
    public void SlopesComeFromOwnExactMAndUElements()
    {
        var f = BuildClaim();
        Assert.Equal(-4, F96BornSubdominantSlopesPi2Inheritance.M3_SingleFlipped);
        Assert.Equal(3, F96BornSubdominantSlopesPi2Inheritance.U2_SingleFlipped_TimesFour);
        Assert.Equal(-20, F96BornSubdominantSlopesPi2Inheritance.M5_DoubleFlipped);
        Assert.Equal(3, F96BornSubdominantSlopesPi2Inheritance.U4_DoubleFlipped_TimesTwo);
        Assert.Equal(-16.0 / 9.0, f.SlopeSingleFlipped, precision: 15);
        Assert.Equal(-8.0 / 3.0, f.SlopeDoubleFlipped, precision: 15);
    }

    [Fact]
    public void LocalFourThirdsComparisonsAreNotDependencies()
    {
        var f = BuildClaim();
        Assert.True(f.SingleFlipSlopeEqualsMinusFourThirdsSquared());
        Assert.True(f.DoubleFlipSlopeEqualsMinusTwoTimesFourThirds());
    }

    [Fact]
    public void AbsoluteThirdOrderDiagonalIsTheTraceControl()
    {
        Assert.Equal(new[] { 8, -4, -4, 0 }, F96BornSubdominantSlopesPi2Inheritance.AbsoluteThirdOrderDiagonal());
        Assert.True(F96BornSubdominantSlopesPi2Inheritance.AbsoluteThirdOrderDiagonalSumsToZero());
    }

    [Fact]
    public void ClaimIsParentlessAndHasNoPolarityAxis()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F96BornSubdominantSlopesPi2Inheritance)));
        var parents = typeof(F96BornSubdominantSlopesPi2Inheritance)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
        Assert.NotNull(typeof(F96BornSubdominantSlopesPi2Inheritance).GetConstructor(Type.EmptyTypes));
    }

    [Theory]
    [InlineData(-1.0)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void DeltaHelpersRejectInvalidK(double k)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().DeltaSingleFlipped(k));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().DeltaDoubleFlipped(k));
    }
}
