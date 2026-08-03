using MirrorWorld;
using Xunit;

namespace MirrorWorld.Tests;

// From below: every pin here is an integer fact, and the two faces of the typed
// NivenRationalityRootClaim must fall out of the SAME routine rather than being restated.
public class CyclotomyTests
{
    [Theory]
    [InlineData(1, 1)] [InlineData(2, 1)] [InlineData(3, 2)] [InlineData(4, 2)]
    [InlineData(6, 2)] [InlineData(8, 4)] [InlineData(12, 4)] [InlineData(14, 6)]
    public void TotientIsTheEulerTotient(int n, int expected)
        => Assert.Equal(expected, Cyclotomy.Totient(n));

    [Theory]
    [InlineData(1, 4, 4)] [InlineData(2, 4, 2)] [InlineData(3, 4, 4)]
    [InlineData(2, 6, 3)] [InlineData(3, 6, 2)] [InlineData(4, 8, 2)]
    public void OrderReducesTheTurnFraction(int j, int n, int expected)
        => Assert.Equal(expected, Cyclotomy.Order(j, n));

    // Niven, as the object states it: exactly five orders give a rational 2*cos.
    [Fact]
    public void RationalOrdersAreExactlyOneTwoThreeFourSix()
    {
        for (int d = 1; d <= 60; d++)
            Assert.Equal(Cyclotomy.Degree(d) == 1, Cyclotomy.IsRational(d));
        Assert.Equal(new[] { 1, 2, 3, 4, 6 }, Cyclotomy.RationalOrders);
    }

    // The IM face of the typed claim, recomputed here from the turn fraction alone: the band edge is
    // k = 1 on the path, so order 2(N+1) and degree phi(2(N+1))/2. Rational only to N=2, quadratic
    // through N=5 (sqrt2, golden, sqrt3), first cubic at N=6.
    [Theory]
    [InlineData(1, 1)] [InlineData(2, 1)] [InlineData(3, 2)] [InlineData(4, 2)]
    [InlineData(5, 2)] [InlineData(6, 3)] [InlineData(7, 4)]
    public void BandEdgeDegreeMatchesTheTypedClaim(int n, int expectedDegree)
        => Assert.Equal(expectedDegree, Cyclotomy.Degree(Cyclotomy.Order(1, 2 * (n + 1))));

    [Fact]
    public void BandEdgeIsRationalExactlyThroughNTwo()
    {
        for (int n = 1; n <= 30; n++)
            Assert.Equal(n <= 2, Cyclotomy.IsRational(Cyclotomy.Order(1, 2 * (n + 1))));
    }

    // The RE face: the rates are all rational exactly when N+1 is in {1,2,3,4,6}, i.e. N in
    // {1,2,3,5} over this range -- N=3 the last before the gap, N=5 the island. Same routine,
    // doubled angle.
    [Theory]
    [InlineData(1, true)] [InlineData(2, true)] [InlineData(3, true)]
    [InlineData(4, false)] [InlineData(5, true)] [InlineData(6, false)]
    [InlineData(7, false)] [InlineData(11, false)]
    public void RatesAllRationalExactlyWhenNPlusOneIsANivenOrder(int n, bool expected)
        => Assert.Equal(expected, Cyclotomy.AllRational(Cyclotomy.PathRateOrders(n)));

    [Fact]
    public void TheRateIslandAtNFiveIsRealAndIsolated()
    {
        Assert.True(Cyclotomy.AllRational(Cyclotomy.PathRateOrders(5)));
        Assert.False(Cyclotomy.AllRational(Cyclotomy.PathRateOrders(4)));
        Assert.False(Cyclotomy.AllRational(Cyclotomy.PathRateOrders(6)));
        Assert.Equal(new[] { 6, 3, 2, 3, 6 }, Cyclotomy.PathRateOrders(5));
    }

    // The ring is a DIFFERENT own output, and it is cleaner: k/N never doubles, so an even N keeps
    // the order-2 fraction and the whole spectrum of a Niven-order ring is rational.
    [Theory]
    [InlineData(3, true)] [InlineData(4, true)] [InlineData(6, true)]
    [InlineData(5, false)] [InlineData(7, false)] [InlineData(8, false)]
    public void RingIsRationalExactlyAtNivenOrderRings(int n, bool expected)
        => Assert.Equal(expected, Cyclotomy.AllRational(Cyclotomy.RingOrders(n)));

    // The point the object exists to make: path and ring are different OWN outputs, so the same N
    // can be rational on one and not the other. A boundary condition is not a physics change, and
    // any reading that calls one of these 'rocks' and the other 'breaking' is reading our choice.
    [Fact]
    public void PathAndRingDisagreeAtTheSameN()
    {
        Assert.True(Cyclotomy.AllRational(Cyclotomy.RingOrders(4)));
        Assert.False(Cyclotomy.AllRational(Cyclotomy.PathRateOrders(4)));
    }

    // The two-bucket ontology, exercised rather than described: the turn fractions are ours, and
    // there is nothing above us to inherit the arithmetic from.
    [Fact]
    public void OwnHoldsTheTurnFractionsAndInheritedIsEmpty()
    {
        var c = new Cyclotomy();
        Assert.Equal(2, c.Own.Count);
        Assert.All(c.Own, o => Assert.Contains("turn fractions", o));
        Assert.Empty(c.Inherited);
    }
}
