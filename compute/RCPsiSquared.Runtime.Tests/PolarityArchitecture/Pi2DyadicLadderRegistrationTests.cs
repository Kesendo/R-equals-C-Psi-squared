using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2DyadicLadderRegistrationTests
{
    [Fact]
    public void RegisterPi2DyadicLadder_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .Build();

        Assert.True(registry.Contains<Pi2DyadicLadderClaim>());
    }

    [Fact]
    public void RegisterPi2DyadicLadder_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<Pi2DyadicLadderClaim>().Tier);
    }

    [Fact]
    public void RegisterPi2DyadicLadder_AncestorsContainsAllThreeKnownAnchors()
    {
        // The ladder declares one edge per known anchor. AncestorsOf the ladder must
        // include all three: QubitDimensionalAnchorClaim (n=0), HalfAsStructuralFixedPointClaim
        // (n=2), QuarterAsBilinearMaxvalClaim (n=3).
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .Build();

        var ancestors = registry.AncestorsOf<Pi2DyadicLadderClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
    }

    [Fact]
    public void RegisterPi2DyadicLadder_LadderTermsMatchAnchorValues()
    {
        // Cross-verification through the registry: the ladder's KnownAnchors[i].Value
        // matches Term(KnownAnchors[i].N) for each entry. Drift between the ladder and
        // its anchors is a structural inconsistency that the registry now catches.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .Build();

        var ladder = registry.Get<Pi2DyadicLadderClaim>();
        foreach (var a in ladder.KnownAnchors)
            Assert.Equal(ladder.Term(a.N), a.Value, precision: 12);
    }

    [Fact]
    public void RegisterPi2DyadicLadder_WithoutPi2Family_Throws()
    {
        // Missing all three Pi2 anchor Claims surfaces as MissingParent at Build().
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2DyadicLadder()
                .Build());
    }
}
