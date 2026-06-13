using System;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class SymmetricGammaSliceTests
{
    [Fact]
    public void Profile_N5_ReproducesTheThreeAnchorsExactly()
    {
        Assert.Equal(new[] { 1.0, 1.0, 1.0, 1.0, 1.0 }, SymmetricGammaSlice.Profile(5, 1.0, 1.0));
        Assert.Equal(new[] { 0.25, 0.75, 3.0, 0.75, 0.25 }, SymmetricGammaSlice.Profile(5, 0.25, 3.0));
        Assert.Equal(new[] { 0.25, 1.5, 1.5, 1.5, 0.25 }, SymmetricGammaSlice.Profile(5, 0.25, 1.5));
    }

    [Theory]
    [InlineData(5, 0.3, 2.0)]
    [InlineData(6, 0.5, 1.0)]
    public void Profile_SumsToN(int n, double wEdge, double wCenter)
    {
        double sum = 0;
        foreach (double w in SymmetricGammaSlice.Profile(n, wEdge, wCenter)) sum += w;
        Assert.Equal(n, sum, 10);
    }

    [Fact]
    public void Profile_EvenN_PutsCenterWeightOnBothMiddleSites()
    {
        // N=6: [e, b, c, c, b, e], b = 3 - e - c = 1.5
        Assert.Equal(new[] { 0.5, 1.5, 1.0, 1.0, 1.5, 0.5 }, SymmetricGammaSlice.Profile(6, 0.5, 1.0));
    }

    [Fact]
    public void Profile_N4_Throws_NoBulkReservoir()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SymmetricGammaSlice.Profile(4, 0.5, 1.0));
    }

    [Fact]
    public void Profile_ThrowsWhenBulkNonPositive()
    {
        // 2*1.5 + 2.5 = 5.5 > 5  ->  w_bulk = (5 - 5.5)/2 < 0
        Assert.Throws<ArgumentException>(() => SymmetricGammaSlice.Profile(5, 1.5, 2.5));
    }

    [Fact]
    public void IsAdmissible_MatchesTheBulkPositivityConstraint()
    {
        Assert.True(SymmetricGammaSlice.IsAdmissible(5, 0.25, 1.5));   // bulk 1.5 > 0
        Assert.False(SymmetricGammaSlice.IsAdmissible(5, 1.5, 2.5));   // bulk < 0
        Assert.False(SymmetricGammaSlice.IsAdmissible(4, 0.5, 1.0));   // no bulk reservoir
    }
}
