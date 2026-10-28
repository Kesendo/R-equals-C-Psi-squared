using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F96BornSubdominantSlopesTests
{
    private static F96BornSubdominantSlopes BuildClaim() => new();

    [Fact]
    public void SlopesComeFromOwnExactMAndUElements()
    {
        var f = BuildClaim();
        Assert.Equal(-4, F96BornSubdominantSlopes.M3_SingleFlipped);
        Assert.Equal(3, F96BornSubdominantSlopes.U2_SingleFlipped_TimesFour);
        Assert.Equal(-20, F96BornSubdominantSlopes.M5_DoubleFlipped);
        Assert.Equal(3, F96BornSubdominantSlopes.U4_DoubleFlipped_TimesTwo);
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
        Assert.Equal(new[] { 8, -4, -4, 0 }, F96BornSubdominantSlopes.AbsoluteThirdOrderDiagonal());
        Assert.True(F96BornSubdominantSlopes.AbsoluteThirdOrderDiagonalSumsToZero());
    }

    [Fact]
    public void ClaimIsParentlessAndHasNoPolarityAxis()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F96BornSubdominantSlopes)));
        var parents = typeof(F96BornSubdominantSlopes)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
        Assert.NotNull(typeof(F96BornSubdominantSlopes).GetConstructor(Type.EmptyTypes));
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
