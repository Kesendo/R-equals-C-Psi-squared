using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class QuditProductMirrorCapRegistrationTests
{
    /// <summary>Minimal registry: the claim plus its parents QuditPartialPalindromeCeiling and
    /// QubitNecessityPi2Inheritance and their transitive Pi2-Foundation chain (copied verbatim
    /// from QuditPartialPalindromeCeilingRegistrationTests.BuildMinimalRegistry).</summary>
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
            .Build();

    [Fact]
    public void Register_AddsClaim()
    {
        Assert.True(BuildMinimalRegistry().Contains<QuditProductMirrorCap>());
    }

    [Fact]
    public void Register_TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildMinimalRegistry().Get<QuditProductMirrorCap>().Tier);
    }

    [Fact]
    public void Ancestors_ContainBothTypedParents()
    {
        var ancestors = BuildMinimalRegistry().AncestorsOf<QuditProductMirrorCap>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(QuditPartialPalindromeCeiling), ancestors);
        Assert.Contains(typeof(QubitNecessityPi2Inheritance), ancestors);
    }

    [Fact]
    public void Battery_AllPass_InRegisteredClaim()
    {
        var claim = BuildMinimalRegistry().Get<QuditProductMirrorCap>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void BuildDefault_ContainsClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<QuditProductMirrorCap>());
        var claim = registry.Get<QuditProductMirrorCap>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
