using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F89F88aKleinPpAnchorRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterHalfIntegerMirror(N: 5)
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .RegisterF89TopologyOrbitClosure();

    [Fact]
    public void RegisterF89F88aKleinPpAnchor_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F88aKleinPpAnchor()
            .Build();

        Assert.True(registry.Contains<F89F88aKleinPpAnchor>());
    }

    [Fact]
    public void RegisterF89F88aKleinPpAnchor_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F88aKleinPpAnchor()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<F89F88aKleinPpAnchor>().Tier);
    }

    [Fact]
    public void RegisterF89F88aKleinPpAnchor_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F88aKleinPpAnchor()
            .Build();

        var ancestors = registry.AncestorsOf<F89F88aKleinPpAnchor>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
        Assert.Contains(typeof(F89TopologyOrbitClosure), ancestors);
    }

    [Fact]
    public void RegisterF89F88aKleinPpAnchor_WithoutF89TopologyOrbitClosure_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterAbsorptionTheoremClaim()
                .RegisterHalfIntegerMirror(N: 5)
                .RegisterF70DeltaNSelectionRulePi2Inheritance()
                .RegisterF72BlockDiagonalPurityPi2Inheritance()
                .RegisterF73SpatialSumPurityClosurePi2Inheritance()
                .RegisterF71MirrorSymmetryPi2Inheritance()
                // Missing: RegisterF89TopologyOrbitClosure
                .RegisterF89F88aKleinPpAnchor()
                .Build());
    }
}
