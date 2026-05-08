using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class QubitNecessityPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterQubitNecessityPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterQubitNecessityPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<QubitNecessityPi2Inheritance>());
    }

    [Fact]
    public void RegisterQubitNecessityPi2Inheritance_AncestorsContainAllFourPi2Anchors()
    {
        // The QubitNecessity per-site basis split inherits from all four typed Pi2
        // anchors: PolynomialFoundation, HalfAsStructuralFixedPoint, Pi2DyadicLadder,
        // Pi2OperatorSpaceMirror.
        var registry = BuildBaseRegistry()
            .RegisterQubitNecessityPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<QubitNecessityPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Fact]
    public void RegisterQubitNecessityPi2Inheritance_BijectionHolds()
    {
        // Cross-registry verification: the bijection equation 4 = 2 · 2 holds through
        // the registered claim, with cross-check against operator-space mirror.
        var registry = BuildBaseRegistry()
            .RegisterQubitNecessityPi2Inheritance()
            .Build();

        var q = registry.Get<QubitNecessityPi2Inheritance>();
        Assert.True(q.BijectionHolds);
        Assert.Equal(4.0, q.TotalPauliOpsPerSite, precision: 14);
        Assert.Equal(2.0, q.ImmuneOpsPerSite, precision: 14);
        Assert.Equal(2.0, q.DecayingOpsPerSite, precision: 14);
        Assert.Equal(0.5, q.BalancedFraction, precision: 14);
    }

    [Fact]
    public void RegisterQubitNecessityPi2Inheritance_WithoutOperatorSpaceMirror_Throws()
    {
        // The mirror anchor is required for the d² = 4 cross-check.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing F88* + Pi2OperatorSpaceMirror
                .RegisterQubitNecessityPi2Inheritance()
                .Build());
    }
}
