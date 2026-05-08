using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2OperatorSpaceMirrorRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor();

    [Fact]
    public void RegisterPi2OperatorSpaceMirror_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2OperatorSpaceMirror()
            .Build();

        Assert.True(registry.Contains<Pi2OperatorSpaceMirrorClaim>());
    }

    [Fact]
    public void RegisterPi2OperatorSpaceMirror_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2OperatorSpaceMirror()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<Pi2OperatorSpaceMirrorClaim>().Tier);
    }

    [Fact]
    public void RegisterPi2OperatorSpaceMirror_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2OperatorSpaceMirror()
            .Build();

        var ancestors = registry.AncestorsOf<Pi2OperatorSpaceMirrorClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(F88StaticDyadicAnchor), ancestors);
    }

    [Fact]
    public void RegisterPi2OperatorSpaceMirror_PairsCrossVerifyWithLadderTerms()
    {
        // The Schicht-1 wiring exposes the inversion identity through both parents:
        // the registered Pi2OperatorSpaceMirrorClaim's pinned values must match
        // Pi2DyadicLadderClaim.Term at the same indices when both come through the
        // registry. Drift between them is caught here.
        var registry = BuildBaseRegistry()
            .RegisterPi2OperatorSpaceMirror()
            .Build();

        var mirror = registry.Get<Pi2OperatorSpaceMirrorClaim>();
        var ladder = registry.Get<Pi2DyadicLadderClaim>();

        foreach (var p in mirror.Pairs)
        {
            Assert.Equal(ladder.Term(p.UpperIndex), p.OperatorSpace, precision: 12);
            Assert.Equal(ladder.Term(p.LowerIndex), p.MirrorMass, precision: 12);
            Assert.Equal(1.0, p.OperatorSpace * p.MirrorMass, precision: 12);
        }
    }

    [Fact]
    public void RegisterPi2OperatorSpaceMirror_N1Pair_LowerSideEqualsQuarterAsBilinearMaxval()
    {
        // Specific verification: the N=1 mirror-pair's lower side (1/4) IS the value
        // pinned by QuarterAsBilinearMaxvalClaim. The two Tier1Derived facts coincide
        // at this point on the ladder.
        var registry = BuildBaseRegistry()
            .RegisterPi2OperatorSpaceMirror()
            .Build();

        var mirror = registry.Get<Pi2OperatorSpaceMirrorClaim>();
        var n1 = mirror.PairAt(1)!;
        Assert.Equal(0.25, n1.MirrorMass, precision: 12);
        Assert.Equal(4.0, n1.OperatorSpace, precision: 12);
    }

    [Fact]
    public void RegisterPi2OperatorSpaceMirror_WithoutF88Anchor_Throws()
    {
        // Without F88StaticDyadicAnchor registered, the inheritance edge to N≥2 anchors
        // cannot be drawn.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88PopcountCoherence()
                // Missing: RegisterF88StaticDyadicAnchor
                .RegisterPi2OperatorSpaceMirror()
                .Build());
    }
}
