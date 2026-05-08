using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F88StaticDyadicAnchorRegistrationTests
{
    [Fact]
    public void RegisterF88StaticDyadicAnchor_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .Build();

        Assert.True(registry.Contains<F88StaticDyadicAnchor>());
    }

    [Fact]
    public void RegisterF88StaticDyadicAnchor_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<F88StaticDyadicAnchor>().Tier);
    }

    [Fact]
    public void RegisterF88StaticDyadicAnchor_AncestorsContainsBothParents()
    {
        // The inheritance reading: F88StaticDyadicAnchor descends from BOTH the Pi2 ladder
        // (the lattice positions) AND the F88 closed form (the binomial structure that
        // realises 1/(2N) at dyadic N). Both parent edges must surface in the registry.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .Build();

        var ancestors = registry.AncestorsOf<F88StaticDyadicAnchor>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(PopcountCoherenceClaim), ancestors);
    }

    [Fact]
    public void RegisterF88StaticDyadicAnchor_DyadicWitnessesMatchLadderTerms()
    {
        // The structural identity through the registry: each pinned witness's StaticFraction
        // equals Pi2DyadicLadder.Term(k+2). If either side drifts, this test catches it.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .Build();

        var anchor = registry.Get<F88StaticDyadicAnchor>();
        var ladder = registry.Get<Pi2DyadicLadderClaim>();

        foreach (var w in anchor.Witnesses)
            Assert.Equal(ladder.Term(w.LadderIndex), w.StaticFraction, precision: 12);
    }

    [Fact]
    public void RegisterF88StaticDyadicAnchor_WithoutParents_Throws()
    {
        // Without the Pi2 ladder registered, the inheritance edge cannot be drawn.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterF88PopcountCoherence()
                // Missing: RegisterPi2DyadicLadder
                .RegisterF88StaticDyadicAnchor()
                .Build());
    }
}
