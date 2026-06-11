using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class QuditPartialPalindromeCeilingRegistrationTests
{
    /// <summary>Minimal registry: the claim plus its parent QubitNecessityPi2Inheritance and
    /// that parent's transitive Pi2-Foundation chain (copied verbatim from
    /// QubitNecessityPi2InheritanceRegistrationTests.BuildBaseRegistry).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterQubitNecessityPi2Inheritance()
            .RegisterQuditPartialPalindromeCeiling()
            .Build();

    [Fact]
    public void Register_AddsClaim()
    {
        Assert.True(BuildMinimalRegistry().Contains<QuditPartialPalindromeCeiling>());
    }

    [Fact]
    public void Register_TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildMinimalRegistry().Get<QuditPartialPalindromeCeiling>().Tier);
    }

    [Fact]
    public void Ancestors_ContainQubitNecessityParent()
    {
        var ancestors = BuildMinimalRegistry().AncestorsOf<QuditPartialPalindromeCeiling>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(QubitNecessityPi2Inheritance), ancestors);
    }

    [Fact]
    public void Battery_AllPass_InRegisteredClaim()
    {
        var claim = BuildMinimalRegistry().Get<QuditPartialPalindromeCeiling>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void BuildDefault_ContainsClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<QuditPartialPalindromeCeiling>());
        var claim = registry.Get<QuditPartialPalindromeCeiling>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
