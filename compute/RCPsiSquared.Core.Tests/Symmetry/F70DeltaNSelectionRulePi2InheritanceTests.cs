using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F70DeltaNSelectionRulePi2InheritanceTests
{
    private static F70DeltaNSelectionRulePi2Inheritance BuildClaim() =>
        new F70DeltaNSelectionRulePi2Inheritance(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void SingleSiteMaxDeltaN_IsExactlyOne_FromSelfMirrorPivot()
    {
        // a_1 = 1 = self-mirror pivot, F77 anchor.
        Assert.Equal(1.0, BuildClaim().SingleSiteMaxDeltaN, precision: 14);
    }

    [Fact]
    public void PairMaxDeltaN_IsExactlyTwo_FromPolynomialRoot()
    {
        // a_0 = 2 = polynomial root d, F1/F66 anchor.
        Assert.Equal(2.0, BuildClaim().PairMaxDeltaN, precision: 14);
    }

    [Theory]
    [InlineData(1, 1)]
    [InlineData(2, 2)]
    [InlineData(3, 3)]
    [InlineData(5, 5)]
    public void PartialTraceMaxDeltaN_EqualsK(int k, int expected)
    {
        Assert.Equal(expected, BuildClaim().PartialTraceMaxDeltaN(k));
    }

    [Theory]
    [InlineData(1, 0, true)]    // single-site sees |ΔN| = 0
    [InlineData(1, 1, true)]    // single-site sees |ΔN| = 1
    [InlineData(1, 2, false)]   // single-site does NOT see |ΔN| = 2
    [InlineData(2, 1, true)]    // pair sees |ΔN| = 1
    [InlineData(2, 2, true)]    // pair sees |ΔN| = 2
    [InlineData(2, 3, false)]   // pair does NOT see |ΔN| = 3
    [InlineData(3, 3, true)]    // triple sees |ΔN| = 3
    [InlineData(3, 4, false)]   // triple does NOT see |ΔN| = 4
    public void IsVisibleToKLocalTrace_FollowsKBoundary(int kLocal, int deltaN, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsVisibleToKLocalTrace(kLocal, deltaN));
    }

    [Theory]
    [InlineData(1, true)]     // a_1 = 1
    [InlineData(2, true)]     // a_0 = 2
    [InlineData(3, false)]    // not on ladder
    [InlineData(5, false)]
    public void MaxDeltaNIsLadderAnchor_HoldsForKOneAndTwo(int k, bool expected)
    {
        Assert.Equal(expected, BuildClaim().MaxDeltaNIsLadderAnchor(k));
    }

    [Theory]
    [InlineData(1, 1)]
    [InlineData(2, 0)]
    [InlineData(3, null)]
    [InlineData(7, null)]
    public void LadderIndexForKLocalThreshold_LandsOnAnchorOrNull(int k, int? expected)
    {
        Assert.Equal(expected, BuildClaim().LadderIndexForKLocalThreshold(k));
    }

    [Fact]
    public void SingleSiteThresholdMatchesSelfMirrorPivot_HoldsExactly()
    {
        Assert.True(BuildClaim().SingleSiteThresholdMatchesSelfMirrorPivot());
    }

    [Fact]
    public void PairThresholdMatchesPolynomialRoot_HoldsExactly()
    {
        Assert.True(BuildClaim().PairThresholdMatchesPolynomialRoot());
    }

    [Fact]
    public void PartialTraceMaxDeltaN_KLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PartialTraceMaxDeltaN(0));
    }

    [Fact]
    public void IsVisibleToKLocalTrace_NegativeDeltaN_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().IsVisibleToKLocalTrace(1, -1));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F70DeltaNSelectionRulePi2Inheritance(null!));
    }
}
