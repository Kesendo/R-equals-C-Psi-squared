using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F81Pi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF81Pi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF81Pi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F81Pi2Inheritance>());
    }

    [Fact]
    public void RegisterF81Pi2Inheritance_AncestorsContainAllThreeParents()
    {
        // Including the Mirror Space connection: Pi2OperatorSpaceMirror appears as
        // an ancestor, not just the dyadic ladder + half-fixed-point.
        var registry = BuildBaseRegistry()
            .RegisterF81Pi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F81Pi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Fact]
    public void RegisterF81Pi2Inheritance_TwoFactorAndHalfFactor_AreMirrorPartners()
    {
        // Cross-registry verification: F81's "2" and "1/2" multiply to 1 exactly.
        var registry = BuildBaseRegistry()
            .RegisterF81Pi2Inheritance()
            .Build();

        var f = registry.Get<F81Pi2Inheritance>();
        Assert.Equal(2.0, f.TwoFactor, precision: 14);
        Assert.Equal(0.5, f.HalfFactor, precision: 14);
        Assert.Equal(1.0, f.TwoTimesHalf, precision: 14);
    }

    [Fact]
    public void RegisterF81Pi2Inheritance_OperatorSpaceDimensionAtN3_Is64()
    {
        // The Mirror Space connection: F81's M lives in d² = 4^N. At N=3 → 64.
        var registry = BuildBaseRegistry()
            .RegisterF81Pi2Inheritance()
            .Build();

        var f = registry.Get<F81Pi2Inheritance>();
        Assert.Equal(64.0, f.OperatorSpaceDimension(3), precision: 12);
    }

    [Fact]
    public void RegisterF81Pi2Inheritance_WithoutMirror_Throws()
    {
        // Mirror Space connection is required; without it the inheritance edge fails.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: F88* + Pi2OperatorSpaceMirror
                .RegisterF81Pi2Inheritance()
                .Build());
    }
}
