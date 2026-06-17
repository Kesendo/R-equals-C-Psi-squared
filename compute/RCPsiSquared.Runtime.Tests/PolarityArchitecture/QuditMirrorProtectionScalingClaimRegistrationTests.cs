using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class QuditMirrorProtectionScalingClaimRegistrationTests
{
    /// <summary>Minimal registry: the claim plus its parent QuditProductMirrorCap and the transitive
    /// Pi2-Foundation chain (copied from QuditProductMirrorCapRegistrationTests.BuildMinimalRegistry,
    /// with the scaling corollary appended).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterQubitNecessityPi2Inheritance()
            .RegisterQuditPartialPalindromeCeiling()
            .RegisterQuditProductMirrorCap()
            .RegisterQuditMirrorProtectionScalingClaim()
            .Build();

    [Fact]
    public void Register_AddsClaim()
    {
        Assert.True(BuildMinimalRegistry().Contains<QuditMirrorProtectionScalingClaim>());
    }

    [Fact]
    public void Register_TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived,
            BuildMinimalRegistry().Get<QuditMirrorProtectionScalingClaim>().Tier);
    }

    [Fact]
    public void Ancestors_ContainTheCapParent()
    {
        var ancestors = BuildMinimalRegistry().AncestorsOf<QuditMirrorProtectionScalingClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(QuditProductMirrorCap), ancestors);
        // the cap's own parents stay in the chain
        Assert.Contains(typeof(QuditPartialPalindromeCeiling), ancestors);
    }

    [Fact]
    public void BuildDefault_ContainsClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<QuditMirrorProtectionScalingClaim>());
        Assert.Equal(Tier.Tier1Derived,
            registry.Get<QuditMirrorProtectionScalingClaim>().Tier);
    }
}
