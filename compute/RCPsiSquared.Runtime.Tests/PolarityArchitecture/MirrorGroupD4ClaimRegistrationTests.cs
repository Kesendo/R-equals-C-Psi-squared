using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class MirrorGroupD4ClaimRegistrationTests
{
    /// <summary>Minimal standalone registry: the claim plus its three typed parents
    /// (and KleinFourCellClaim, the cube's own Klein2 parent).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            .Register<KleinFourCellClaim>(_ => new KleinFourCellClaim())
            .RegisterKleinEightCellClaim()
            .RegisterPi2KleinV4DephaseSwapGroup()
            .RegisterCommutatorDConjugationSign()
            .RegisterMirrorGroupD4Claim()
            .Build();

    [Fact]
    public void RegisterMirrorGroupD4Claim_AddsClaim()
    {
        var registry = BuildMinimalRegistry();
        Assert.True(registry.Contains<MirrorGroupD4Claim>());
    }

    [Fact]
    public void RegisterMirrorGroupD4Claim_TierIsTier1Derived()
    {
        var registry = BuildMinimalRegistry();
        Assert.Equal(Tier.Tier1Derived, registry.Get<MirrorGroupD4Claim>().Tier);
    }

    [Fact]
    public void Ancestors_ContainAllThreeTypedParents()
    {
        var registry = BuildMinimalRegistry();
        var ancestors = registry.AncestorsOf<MirrorGroupD4Claim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(KleinEightCellClaim), ancestors);
        Assert.Contains(typeof(CommutatorDConjugationSign), ancestors);
        Assert.Contains(typeof(Pi2KleinV4DephaseSwapGroup), ancestors);
        // Transitive: the cube's own Klein2 parent.
        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
    }

    [Fact]
    public void Battery_AllPass_InRegisteredClaim()
    {
        var registry = BuildMinimalRegistry();
        var claim = registry.Get<MirrorGroupD4Claim>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void BuildDefault_ContainsMirrorGroupD4Claim()
    {
        // The claim is wired into the production registry via KnowledgeRegistryFactory
        // (registered directly after its later parent KleinEightCellClaim, 2026-06-10).
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<MirrorGroupD4Claim>());
        var claim = registry.Get<MirrorGroupD4Claim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
