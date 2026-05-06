using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.ObjectManager;

public class TierStrengthTests
{
    [Fact]
    public void Of_Tier1Derived_IsHighest()
    {
        Assert.Equal(5, TierStrength.Of(Tier.Tier1Derived));
    }

    [Fact]
    public void Of_Retracted_IsLowest()
    {
        Assert.Equal(0, TierStrength.Of(Tier.Retracted));
    }

    [Fact]
    public void Of_VerifiedStrongerThanEmpirical()
    {
        Assert.True(TierStrength.Of(Tier.Tier2Verified) > TierStrength.Of(Tier.Tier2Empirical));
    }

    [Fact]
    public void IsAtLeastAsStrong_Tier1DerivedParent_Tier1DerivedChild_True()
    {
        Assert.True(TierStrength.IsAtLeastAsStrong(Tier.Tier1Derived, Tier.Tier1Derived));
    }

    [Fact]
    public void IsAtLeastAsStrong_Tier1CandidateParent_Tier1DerivedChild_False()
    {
        Assert.False(TierStrength.IsAtLeastAsStrong(Tier.Tier1Candidate, Tier.Tier1Derived));
    }
}
