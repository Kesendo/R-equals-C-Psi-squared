using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of F168 under a twisted coupling (<c>experiments/THE_CAT_UNDER_A_TURNED_FIELD.md</c>,
/// section "Under a twisted coupling"): twists and own field angles at once, one bit unwatched, the odd
/// moments read exactly in the Pauli basis. Every comparison is exact (rational arithmetic).</summary>
public class CatTwistedCouplingWitnessTests
{
    [Theory]
    [InlineData(3, false)]
    [InlineData(4, false)]
    [InlineData(5, false)]
    [InlineData(3, true)]
    [InlineData(4, true)]
    [InlineData(5, true)]
    public void OrderSeven_IsNinetySixTimesTheWatchedDivergenceSquared(int n, bool ring)
    {
        var w = new CatTwistedCouplingWitness(n, ring);
        Assert.True(w.LawHolds, $"m_7 = {w.MatchedMoments[7]} against {w.LawValue}; {w.Summary}");
    }

    [Theory]
    [InlineData(3, false)]
    [InlineData(4, false)]
    [InlineData(5, false)]
    [InlineData(4, true)]
    [InlineData(5, true)]
    public void OtherPhase_SameLawFromFourOn_DifferentAtThree(int n, bool ring)
    {
        var w = new CatTwistedCouplingWitness(n, ring);
        Assert.True(w.OffPhaseReadingHolds, $"off-phase m_7 = {w.OffPhaseMoments[7]}, law {w.LawValue}");
    }

    [Theory]
    [InlineData(3, false)]
    [InlineData(4, false)]
    [InlineData(3, true)]
    [InlineData(4, true)]
    public void FollowingTheTwist_IsFoneSixtyEightWithTheSumOfAngles(int n, bool ring)
    {
        var w = new CatTwistedCouplingWitness(n, ring);
        Assert.True(w.FollowingVanishes, w.Summary);
        Assert.True(w.FollowingLawHolds, $"m_{2 * n + 1} = {w.FollowingOffMoments[2 * n + 1]} against {w.FollowingLawValue}");
    }

    [Theory]
    [InlineData(4, false)]
    [InlineData(4, true)]
    public void Control_BondByBondSquaresDiffer(int n, bool ring) =>
        Assert.True(new CatTwistedCouplingWitness(n, ring).ControlDiffers);

    [Theory]
    [InlineData(3, false)]
    [InlineData(4, false)]
    [InlineData(3, true)]
    [InlineData(4, true)]
    public void ConstantCurrent_EndSourcesOnTheChain_InvisibleUntilNineOnTheRing(int n, bool ring)
    {
        var w = new CatTwistedCouplingWitness(n, ring);
        Assert.True(w.ConstantCurrentReadingHolds,
            $"m_7 = {w.ConstantCurrentMoments[7]} against {w.EndSourceValue}, m_9 = {w.ConstantCurrentMoments[9]}");
    }

    [Fact]
    public void ClaimHelpers_ReadByHand()
    {
        // δ = 0, φ_a = 0, φ_b = π/2: β = −π/2, sin β = −1
        var zero = (BigRational.One, BigRational.Zero);
        var quarter = (BigRational.Zero, BigRational.One);
        Assert.Equal(-BigRational.One, CatTwistedCouplingClaim.SinResidualTwist(zero, zero, quarter));
        // a chain 0 → 1 → 2 with currents 2 and 2: sources at the ends, none in the middle
        var div = CatTwistedCouplingClaim.Divergence(3, new[] { (0, 1, new BigRational(2)), (1, 2, new BigRational(2)) });
        Assert.Equal(new[] { new BigRational(2), BigRational.Zero, new BigRational(-2) }, div);
        // 96·(1·4 + 0 + 1·4) = 768
        Assert.Equal(new BigRational(768),
            CatTwistedCouplingClaim.OrderSevenMoment(new[] { BigRational.One, BigRational.One, BigRational.One }, div));
        Assert.Throws<System.ArgumentException>(() =>
            CatTwistedCouplingClaim.SinResidualTwist((new BigRational(1, 2), new BigRational(1, 2)), zero, zero));
    }

    [Fact]
    public void Guards()
    {
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CatTwistedCouplingWitness(2));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CatTwistedCouplingWitness(6));
    }
}
